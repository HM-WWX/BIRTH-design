rm(list=ls())
# 加载必要的包
library(doParallel)
library(foreach)
library(psrwe)
library(MASS)
library(boot)
library(coda)
library(rjags)
library(RBesT)
library(abind)
library(R2jags)
library(tictoc)
library(rstan)
library(lattice)
library(tidyverse)
library(dplyr)
library(do)
library(patchwork) 
library(stringr)
library(ggplot2)
library(ggrepel)
library(scales)
library(invgamma)
library(doSNOW)
library(progress)
library(parallel)
library(Hmisc)
library(writexl)
library(openxlsx)


MAP_only <- function(trt.effect_f,mu.cur_f,mu.history_f){
# 设置并行集群，使用3个核
num_cores <- 40
cl <- makeSOCKcluster(num_cores)
registerDoSNOW(cl)

#init_values <- list(
#  .RNG.name = "base::Mersenne-Twister",  # 指定随机数生成器
#  .RNG.seed = 1234                      # 指定种子
#)

ssize <- 1000


set.seed(1234)  # 主种子
task_seeds <- sample(1:1e6, ssize)


pb <- txtProgressBar(max=ssize, style=3, char = "*",)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

mu.control.all <- c()
hist.numbers <- c(150,150,150,150,150)
con.numbers <- c(20,20,20,20,20)
trt.numbers <- c(40,40,40,40,40)
################################################################################
trt.effect <- trt.effect_f
mu.cur <- mu.cur_f 
mu.history <- mu.history_f
################################################################################
n.subtrial <- length(mu.cur)
Y.cur.control <- list()
Y.cur.treat <- list()
Hist.Prior.theta <- c()
Hist.Prior.sd.theta <- c()
All.trials_con <- c()
All.trials_trt <- c()
Prior.theta.Hist <- list()
Prior.theta.sd.Hist <- list()
subtrial_trt_p <- c()
subtrial_con_p <- c()
subtrial_trt_p_95CI <- c()
subtrial_con_p_95CI <- c()
subtrial_trt <- list()
subtrial_con <- list()
subtrial_trt_95CI <- list()
subtrial_con_95CI <- list()
power_R <- list()
power_counter <- c(0,0,0,0,0)
ESS_con_CPP <- c()
ESS <- list()

 results <- foreach(n = 1:ssize,.combine= "c", .packages = c('psrwe', "do","dplyr",'MASS', 'boot',
                            "lattice","tidyverse","patchwork","R2jags",
                            "scales","invgamma","stringr","rstan","tictoc",
                            'coda', 'rjags', 'RBesT', 'abind',"doSNOW","Hmisc"), .options.snow = opts) %dopar% {
  init_values <- list(
  .RNG.name = "base::Mersenne-Twister",  # 指定随机数生成器
  .RNG.seed =  task_seeds[n]                     # 指定种子
  )
  power_counter <- c(0,0,0,0,0)
  DAW_MAP.prior <- list()
  
  tryCatch({   
  
     for (i in 1: n.subtrial){
    setwd("//home//wwx")
    
    set.seed(task_seeds[n]+i)
    #######################
    ###                 ###
    ###  generate data  ###
    ###                 ###
    #######################

      
      X_convert = function(X, binary_col) {
        
        for (i in 1:length(binary_col)) {
          current_col = X[, i]
          X[which(current_col > 10), i] = 1
          X[which(current_col <= 10), i] = 0
        }
        return(t(X))
      }
      
      
      gen_Y_C = function(X, beta, var_eps, study.effect = 0) {
        
        N = dim(X)[1]
        Y = matrix(0, N)
        for (i in 1:N) {
          xb = t(beta) %*% X[i, ]
          logit_P = xb  + study.effect+ rnorm(1,0,0.25)
          P = boot::inv.logit(logit_P)
          Y[i, ] = rbinom(1,1,P)
        }
        
        return(Y)
      }
      
      
      gen_Y_T = function(X, beta, var_eps, trt.effect) {
        
        N = dim(X)[1]
        Y = matrix(0, N)
        for (i in 1:N) {
          xb = t(beta) %*% X[i, ]
          logit_P = xb + trt.effect+ rnorm(1,0,0.25)
          P = boot::inv.logit(logit_P) 
          Y[i, ] = rbinom(1,1,P)
        }
        
        return(Y)
      }
      
      DAW_weight <- function(X, N_borrow) {
        DAW_weight = c()
        
        weight = X / (1 - X)
        total_weight = sum(weight)
        
        for (j in 1:length(X)) {
          DAW_weight[j] = weight[j] * N_borrow / total_weight
        }
        
        
        return(DAW_weight)
      }
      
      n.hist.samples <- hist.numbers[i] # number of patients in each trial
      # Current 
      n.con.samples <- con.numbers[i] #number of patients in control
      n.trt.samples <- trt.numbers[i]
      
      target.borrow.ESS <- con.numbers[i]
      #Generate covariates X
      p <- 9
      rho <- 0.1 #correlations of the covariates
      beta_0 <- rep(1,p)
      
      ##### PSMAP simulation scenario #####
      mu.hist <- mu.history[i] 
      mu.cur.con <- mu.cur[i] # p=9, muc=0.5
      
      S0 = diag(x = rep(0.25^2, p))
      S0[which(S0 == 0)] = rho * (0.25^2)
      
      hist.list <- list()
      X.hist <- X.hist.init <- list()
      X.hist.init = MASS::mvrnorm(n.hist.samples, mu = rep(mu.hist, p), Sigma = S0)
      X.hist = t(X_convert(X.hist.init, c(1, 2)))
      hist.list = gen_Y_C(X.hist, beta = beta_0, study.effect = 0)
      
      X.hist.combine = abind(X.hist, along=1)
      Y.hist.combine = abind(hist.list, along=1)
      hist.ybar = sapply(hist.list, FUN = sum)
      
      set.seed(task_seeds[n]+i+1)
      ### Current 
      X.cur.con.init = MASS::mvrnorm(n.con.samples, mu = rep(mu.cur.con, p), Sigma = S0)
      X.cur.con = t(X_convert(X.cur.con.init, c(1, 2)))
      
      set.seed(task_seeds[n]+i)
      ### Treat
      X.cur.trt.init = MASS::mvrnorm(n.trt.samples, mu = rep(mu.cur.con, p), Sigma = S0)
      X.cur.trt = t(X_convert(X.cur.trt.init, c(1, 2)))
      
      set.seed(task_seeds[n]+i+1)
      #generate outcomes for control
      y.control = gen_Y_C(X.cur.con, beta = beta_0, study.effect = 0)
      ybar.control = sum(y.control)
      Y.cur.con=y.control
      Y.cur.control[[i]] <-Y.cur.con
      
      set.seed(task_seeds[n]+i)
      #generate outcomes for treat
      y.treat = gen_Y_T(X.cur.trt, beta = beta_0, trt.effect = trt.effect[i])
      ybar.treat = sum(y.treat)
      Y.cur.trt = y.treat
      Y.cur.treat[[i]] <-Y.cur.trt     
      
      S = 1
      Y.combine = rbind(Y.hist.combine, Y.cur.con) 
      dta = rbind(X.hist.combine, X.cur.con)
      covs = NULL
      for (v in 1:p) {
        if(v<10){
          covs[v] = paste("V", v, sep = "")
        }else{
          covs[1:9] = paste("V0", 1:9, sep = "")
          covs[10:v] = paste("V", 10:v, sep = "")
        }
      }
      colnames(dta) = covs
      group = c(rep(0, dim(X.hist.combine)[1]), rep(1, dim(X.cur.con)[1])) #historical data =0; current control=1
      dta = cbind(group, dta)
      dta = as.data.frame(dta)
      
      ## estimate PS by logistic regression
      ana.ps_fun <-function(ps_method){
        ana.ps <-psrwe_est(
          dta,
          v_grp = "group",
          v_covs = covs,
          nstrata = S,
          ps_method = ps_method
        )
        return(ana.ps)
      }
      ana.ps_log <- ana.ps_fun("logistic")
      ana.ovl  <- summary(ana.ps_log, metric = "ovl")
      if(v<10){
        ana.ps_log$data$V1 <- factor(ana.ps_log$data$V1)
        ana.ps_log$data$V2 <- factor(ana.ps_log$data$V2)
      }else{
        ana.ps_log$data$V01 <- factor(ana.ps_log$data$V01)
        ana.ps_log$data$V02 <- factor(ana.ps_log$data$V02)
      }
      ###########################
      ###                     ###
      ###  DAW-MAP vs PS-MAP  ###
      ###                     ###
      ###########################
      ###########################################################################PS-MAP###########################################################
      rS <- ana.ovl$Summary$Distance[1:S] 
      overlap_coefficients =rS
      PS_data <- ana.ps_log[[1]]
      PS_data <- cbind(PS_data,Y.combine)
      PS_data <- na.omit(PS_data)
      DAW.hist <- PS_data[PS_data$group == 0,]
      DAW.cur <- PS_data[PS_data$group == 1,]
      #target ESS
      target.ESS <- con.numbers[i]
      #DAW_weight 
      DAW.hist$DAW_weight <- DAW_weight(DAW.hist$`_ps_`,N_borrow = target.ESS)
      Ybar.hist <- round(sum(DAW.hist$DAW_weight[DAW.hist$Y.combine == 1]))
      Ybar.hist.nresp <- round(sum(DAW.hist$DAW_weight[DAW.hist$Y.combine == 0]))
      n.cur.stratum <- 1
      n.hist.stratum <- Ybar.hist+Ybar.hist.nresp
      y.control <- DAW.cur$Y.combine
      DAW_MAP.fit = function(tau.init, target.ESS, std_heterogeneity,
                             n.cur.stratum, Ybar.hist, n.hist.stratum, 
                             overlap_coefficients, niter,lim, data.indiv) {
        # Obtain prior samples when using the PS-MAP prior with the specified target ESS
        n.strata = 1
        WS1 = 1
        ess.res = c(0)
        stop = FALSE
        tau.scale = tau.init
        low = lim[1]
        high = lim[2]
        tau.balance <- c()
        iter_count <- 0
        while(stop==FALSE){
          theta.pred = list()
          if((high-low) <= 0.0001){
            stop = TRUE
          }
          #estimate stratum-specific MAP prior
          for(i in 1:n.strata){
            if (overlap_coefficients[i]==0){
              tau.balance[i] <- 0.1
            }else {
              tau.balance[i] <- overlap_coefficients[i]
            }
            dataTemp = list(
              "resp" = Ybar.hist[i],
              "n" = n.hist.stratum[i],
              "tau.scale"=tau.scale,
              "tau.balance"=tau.balance[i]
            )
            dataTemp$std_heterogeneity <- std_heterogeneity
            
            model <-
              jags.model(
                file = "//home//wwx//NEW-WWX.txt",
                data = dataTemp,
                n.chains = 1,
                n.adapt = 5000,
                quiet = TRUE,
                inits = init_values
              )
            update(model, n.iter = niter, progress.bar = 'none',seed=888) # burn in
            MAP_model <-
              coda.samples(model,
                           variable.names = "p.pred.h",
                           thin = 10,
                           n.iter = 10000, 
                           progress.bar = 'none',
                           na.rm = TRUE)
            theta.pred[[i]] = c(MAP_model[[1]][,1])
          }
          theta.pred = do.call(cbind, theta.pred)
          theta = WS1 %*% t(theta.pred)
          mix.res <- automixfit(theta[1,], Nc = 1:5, thresh = 0, type = "beta")
          posterior.indiv <- postmix(mix.res, n = length(data.indiv), r = sum(data.indiv))
          mu.con_fat <- rmix(posterior.indiv, n=10000)
          mu.con_p <- mu.con_fat[seq(1, length(mu.con_fat), by=10)]
          #calculate ESS
          ess.res = c(ess.res, ess(mix.res, method = "elir"))
          if(abs(ess.res[length(ess.res)] - target.ESS) <= 1 ){
            stop = TRUE
          }else if(ess.res[length(ess.res)]-target.ESS > 1 && iter_count <= 50) {
            low = tau.scale
            tau.scale = low+0.03*(high-low)
            iter_count <- iter_count+1
          }else if(ess.res[length(ess.res)]-target.ESS < -1 && iter_count <= 50){
            tau.scale = tau.scale-0.03*(high-low)
            low=lim[1]
            iter_count <- iter_count+1
          }else if (iter_count > 50){
            stop = TRUE
            
          }
          #print(ess.res)
          #print(tau.scale)
          #print(low)
          #print(high)
        }
        results = list(mix.res, ess.res[length(ess.res)],theta,mu.con_p)
        return(results)
      }
      DAW_MAP <- DAW_MAP.fit(tau.init = 0.8, target.ESS = target.ESS, std_heterogeneity = 0.4,
                             n.cur.stratum=n.cur.stratum, Ybar.hist=Ybar.hist, n.hist.stratum=n.hist.stratum, 
                             overlap_coefficients = overlap_coefficients, niter = 15000,lim=c(0.01,1.26),data.indiv = y.control)  
    
    All.p.pred <- as.numeric(DAW_MAP[[3]])
    All.theta <- logit(All.p.pred)
    DAW_MAP.prior[[i]] <- All.theta
    Hist.Prior.theta[i] <- mean(All.theta)
    Hist.Prior.sd.theta[i] <- var(All.theta)
  }
  
##############################Treat group CPP###################################
  resp_trt <- sapply(Y.cur.treat, FUN = sum)
  kMod <- length(resp_trt)
  kMod <- as.numeric(kMod)
  resp_con <- sapply(Y.cur.control, FUN = sum)  
###########################Treat group posterior################################    
  Prior.theta_trt <- rep(0, kMod)
  Prior.sd.theta_trt <- rep(0.0001, kMod)
    dataTemp_p_trt <-list("resp_p"=resp_trt,"n_p"=trt.numbers,"Prior.theta_p"=Prior.theta_trt,
                          "Prior.sd.theta_p"=Prior.sd.theta_trt,"kMod"=kMod) 
    model_posterior_trt <-
      jags.model(
        file = "//home//wwx//Posterior.txt",
        data = dataTemp_p_trt,
        n.chains = 1,
        n.adapt = 8000,
        quiet = TRUE,
        inits=init_values
      )
    update(model_posterior_trt, n.iter = 18000, progress.bar = 'none') # burn in
    MAP_model_p_trt<-
      coda.samples(model_posterior_trt,
                   variable.names = "theta_p",
                   thin = 10,
                   n.iter = 10000, 
                   progress.bar = 'none',
                   na.rm = TRUE)
#############################Control group posterior############################ 
    Prior.theta_con <- Hist.Prior.theta
    Prior.sd.theta_con <- 1/Hist.Prior.sd.theta
################################################################################
    dataTemp_p_con <-list("resp_p"=resp_con,"n_p"=con.numbers,"Prior.theta_p"=Prior.theta_con,
                          "Prior.sd.theta_p"=Prior.sd.theta_con,"kMod"=kMod) 
    model_posterior_con <-
      jags.model(
        file = "//home//wwx//Posterior.txt",
        data = dataTemp_p_con,
        n.chains = 1,
        n.adapt = 8000,
        quiet = TRUE,
        inits=init_values
      )
    update(model_posterior_con, n.iter = 18000, progress.bar = 'none') # burn in
    MAP_model_p_con<-
      coda.samples(model_posterior_con,
                   variable.names = "theta_p",
                   thin = 10,
                   n.iter = 10000, 
                   progress.bar = 'none',
                   na.rm = TRUE)
    for (j in 1:kMod){
    subtrial_con_p[j]<-mean(plogis(MAP_model_p_con[[1]][,j]))
    subtrial_con_p_95CI[j] <- smean.cl.normal(plogis(MAP_model_p_con[[1]][,j]))[3]-smean.cl.normal(plogis(MAP_model_p_con[[1]][,j]))[2]
    No_borrowing_con_p <- rbeta(1000,0.5+resp_con[j],0.5+con.numbers[j]-resp_con[j])
    ESS_con_CPP[j] <- con.numbers[j]*var(No_borrowing_con_p)/var(plogis(MAP_model_p_con[[1]][,j])) 
    
    subtrial_trt_p[j]<-mean(plogis(MAP_model_p_trt[[1]][,j]))
    subtrial_trt_p_95CI[j] <-  smean.cl.normal(plogis(MAP_model_p_trt[[1]][,j]))[3]-smean.cl.normal(plogis(MAP_model_p_trt[[1]][,j]))[2]
    
    trt.response <- as.matrix(MAP_model_p_trt[[1]][,j])
    con.response <- as.matrix(MAP_model_p_con[[1]][,j])
    trt.difference <- trt.response - con.response
    power <- sum(trt.difference >= 0)
   if (power >= 950){
      power_counter[j] <- power_counter[j]+1
    }else{
      power_counter[j] <- power_counter[j]
  }
    
}

  list(subtrial_trt = subtrial_trt_p, subtrial_con = subtrial_con_p, power_R = power_counter, 
       subtrial_trt_95CI=subtrial_trt_p_95CI, subtrial_con_95CI=subtrial_con_p_95CI,ESS=ESS_con_CPP)

  }, error = function(e) {
    # 捕获错误并继续下一个种子
    message("Error in iteration ", n, ": ", e$message)
    message("Skipping to next seed...")
    return(NULL)
  })
}
 
   
  subtrial_trt_all <- matrix(NA,ncol=n.subtrial,nrow=1)
  subtrial_con_all <- matrix(NA,ncol=n.subtrial,nrow=1)
  Power_all <- matrix(NA,ncol=n.subtrial,nrow=1)
  subtrial_trt_95CI_all <- matrix(NA,ncol=n.subtrial,nrow=1)
  subtrial_con_95CI_all <- matrix(NA,ncol=n.subtrial,nrow=1)
  ESS.ALL <- matrix(NA,ncol=n.subtrial,nrow=1)
  
  trial_numbers <- length(results)/6
  
 for (i in 0:(trial_numbers-1)){
 
   subtrial_trt_all <- rbind(subtrial_trt_all,results[[1+6*i]])  
   subtrial_con_all <- rbind(subtrial_con_all,results[[2+6*i]])
   Power_all <-rbind(Power_all,results[[3+6*i]])
   subtrial_trt_95CI_all <- rbind(subtrial_trt_95CI_all,results[[4+6*i]])
   subtrial_con_95CI_all <- rbind(subtrial_con_95CI_all,results[[5+6*i]])
   ESS.ALL <- rbind(ESS.ALL,results[[6+6*i]])
 }


subtrial_trt_all <- subtrial_trt_all[-1,]
subtrial_con_all <- subtrial_con_all[-1,]
Power_all <- Power_all[-1,]
subtrial_trt_95CI_all <- subtrial_trt_95CI_all[-1,]
subtrial_con_95CI_all <- subtrial_con_95CI_all[-1,]
ESS.ALL <- ESS.ALL[-1,]
ESS.ALL <- ESS.ALL[!apply(ESS.ALL, 1, function(x) any(is.infinite(x))), ]

aa1 <- mean(subtrial_con_all[,1])
bb1 <- mean(subtrial_con_all[,2])
cc1 <- mean(subtrial_con_all[,3])
dd1 <- mean(subtrial_con_all[,4])
ee1 <- mean(subtrial_con_all[,5])

subtrial_con_mean <- c(aa1,bb1,cc1,dd1,ee1)

aa2 <- mean(subtrial_trt_all[,1])
bb2 <- mean(subtrial_trt_all[,2])
cc2 <- mean(subtrial_trt_all[,3])
dd2 <- mean(subtrial_trt_all[,4])
ee2 <- mean(subtrial_trt_all[,5])

subtrial_trt_mean <- c(aa2,bb2,cc2,dd2,ee2)

aa3 <- mean((subtrial_con_all[,1]-plogis(mu.cur[1]*7))**2)
bb3 <- mean((subtrial_con_all[,2]-plogis(mu.cur[2]*7))**2)
cc3 <- mean((subtrial_con_all[,3]-plogis(mu.cur[3]*7))**2)
dd3 <- mean((subtrial_con_all[,4]-plogis(mu.cur[4]*7))**2)
ee3 <- mean((subtrial_con_all[,5]-plogis(mu.cur[5]*7))**2)

subtrial_con_mse <- c(aa3,bb3,cc3,dd3,ee3)

aa4 <- mean(subtrial_trt_95CI_all[,1])
bb4 <- mean(subtrial_trt_95CI_all[,2])
cc4 <- mean(subtrial_trt_95CI_all[,3])
dd4 <- mean(subtrial_trt_95CI_all[,4])
ee4 <- mean(subtrial_trt_95CI_all[,5])

subtrial_trt_95CI <- c(aa4,bb4,cc4,dd4,ee4)

aa5 <- mean(subtrial_con_95CI_all[,1])
bb5 <- mean(subtrial_con_95CI_all[,2])
cc5 <- mean(subtrial_con_95CI_all[,3])
dd5 <- mean(subtrial_con_95CI_all[,4])
ee5 <- mean(subtrial_con_95CI_all[,5])

subtrial_con_95CI <- c(aa5,bb5,cc5,dd5,ee5)

aa6 <- mean(Power_all[,1])
bb6 <- mean(Power_all[,2])
cc6 <- mean(Power_all[,3])
dd6 <- mean(Power_all[,4])
ee6 <- mean(Power_all[,5])

subtrial_Power <- c(aa6,bb6,cc6,dd6,ee6) 

aa7 <- mean(ESS.ALL[,1])
bb7 <- mean(ESS.ALL[,2])
cc7 <- mean(ESS.ALL[,3])
dd7 <- mean(ESS.ALL[,4])
ee7 <- mean(ESS.ALL[,5])

subtrial_ESS <- c(aa7,bb7,cc7,dd7,ee7) 

results_final = list(subtrial_con_mean, subtrial_trt_mean,subtrial_con_mse,subtrial_trt_95CI,
               subtrial_con_95CI,subtrial_Power,subtrial_ESS)

return(results_final)

#MAP ONLY
close(pb)
stopImplicitCluster()
stopCluster(cl)

}

MAP_only_results <- createWorkbook()
MAP_only_results_1 <- createWorkbook()
MAP_only_results_2 <- createWorkbook()
MAP_only_results_3 <- createWorkbook()
#Scenario 1
Scenario_1 <- MAP_only(trt.effect <- c(1.386294,1.386294,1.386294,1.386294,1.386294),
                           mu.cur <- c(0,0,0,0,0), 
                           mu.history <- c(0,0,0,0,0))

Scenario_1_results <- do.call(rbind, lapply(Scenario_1, function(x) as.data.frame(t(x))))

#Scenario 2
Scenario_2 <- MAP_only(trt.effect <- c(0.9808293,0.9808293,0.9808293,0.9808293,0.9808293),
                           mu.cur <- c(0.05792359,0.05792359,0.05792359,0.05792359,0.05792359),
                           mu.history <- c(0.05792359,0.05792359,0.05792359,0.05792359,0.05792359))

Scenario_2_results <- do.call(rbind, lapply(Scenario_2, function(x) as.data.frame(t(x))))

#Scenario 3
Scenario_3 <- MAP_only(trt.effect <- c(0.5389965,0.5389965,0.5389965,0.5389965,0.5389965),
                           mu.cur <- c(0.1210426,0.1210426,0.1210426,0.1210426,0.1210426), 
                           mu.history <- c(0.1210426,0.1210426,0.1210426,0.1210426,0.1210426))

Scenario_3_results <- do.call(rbind, lapply(Scenario_3, function(x) as.data.frame(t(x))))

#Scenario 4
Scenario_4 <- MAP_only(trt.effect <- c(0,0,0,0,0),
                           mu.cur <- c(0.1980421,0.1980421,0.1980421,0.1980421,0.1980421), 
                           mu.history <- c(0.1980421,0.1980421,0.1980421,0.1980421,0.1980421))

Scenario_4_results <- do.call(rbind, lapply(Scenario_4, function(x) as.data.frame(t(x))))

#############################################################################################
addWorksheet(MAP_only_results_1, "Sheet1")
writeData(MAP_only_results_1, "Sheet1", x = NULL)
writeData(MAP_only_results_1, sheet = "Sheet1", c(Scenario_1_results,Scenario_2_results,
                                                  Scenario_3_results,Scenario_4_results))
saveWorkbook(MAP_only_results_1, "MAP_only_1_output.xlsx", overwrite = TRUE)
##############################################################################################

#Scenario 5
Scenario_5 <- MAP_only(trt.effect <- c(0.9808293,0.8979416,0.8472979,0.8197099,0.8109302),
                           mu.cur <- c(0.05792359,0.02866724,0,0.02866724,-0.05792359), 
                           mu.history <- c(0.05792359,0.02866724,0,0.02866724,-0.05792359))

Scenario_5_results <- do.call(rbind, lapply(Scenario_5, function(x) as.data.frame(t(x))))

#Scenario 6
Scenario_6 <- MAP_only(trt.effect <- c(0.9808293,0.8472979,0.8109302,0.8472979,0.9808293),
                           mu.cur <- c(0.05792359,0,-0.05792359,-0.1210426,-0.1980421),
                           mu.history <- c(0.05792359,0,-0.05792359,-0.1210426,-0.1980421))

Scenario_6_results <- do.call(rbind, lapply(Scenario_6, function(x) as.data.frame(t(x))))

#Scenario 7
Scenario_7 <- MAP_only(trt.effect <- c(0,0,0,0,0),
                           mu.cur <- c(0.1980421,0.1569446,0.1210426,0.08843417,0.05792359), 
                           mu.history <- c(0.1980421,0.1569446,0.1210426,0.08843417,0.05792359))

Scenario_7_results <- do.call(rbind, lapply(Scenario_7, function(x) as.data.frame(t(x))))

#Scenario 8
Scenario_8 <- MAP_only(trt.effect <- c(0,0,0,0,0),
                           mu.cur <- c(0.1980421,0.1210426,0.05792359,0,-0.05792359), 
                           mu.history <- c(0.1980421,0.1210426,0.05792359,0,-0.05792359))

Scenario_8_results <- do.call(rbind, lapply(Scenario_8, function(x) as.data.frame(t(x))))

#############################################################################################
addWorksheet(MAP_only_results_2, "Sheet1")
writeData(MAP_only_results_2, "Sheet1", x = NULL)
writeData(MAP_only_results_2, sheet = "Sheet1", c(Scenario_5_results,Scenario_6_results,
                                                      Scenario_7_results,Scenario_8_results))
saveWorkbook(MAP_only_results_2, "MAP_only_2_output.xlsx", overwrite = TRUE)
##############################################################################################

#Scenario 9
Scenario_9 <- MAP_only(trt.effect <- c(0.9808293,0.9808293,0.9808293,0.9808293,0),
                           mu.cur <- c(0.05792359,0.05792359,0.05792359,0.05792359,0.1980421), 
                           mu.history <- c(0.05792359,0.05792359,0.05792359,0.05792359,0.1980421))

Scenario_9_results <- do.call(rbind, lapply(Scenario_9, function(x) as.data.frame(t(x))))

#Scenario 10
Scenario_10 <- MAP_only(trt.effect <- c(0.9808293,0.9808293,0.9808293,0,0),
                            mu.cur <- c(0.05792359,0.05792359,0.05792359,0.1980421,0.1980421), 
                            mu.history <- c(0.05792359,0.05792359,0.05792359,0.19804219,0.1980421))

Scenario_10_results <- do.call(rbind, lapply(Scenario_10, function(x) as.data.frame(t(x))))

#Scenario 11
Scenario_11 <- MAP_only(trt.effect <- c(0.9808293,0,0.8109302,0,0.9808293),
                            mu.cur <- c(0.05792359,0.1210426,-0.05792359,0,-0.1980421), 
                            mu.history <- c(0.05792359,0.1210426,-0.05792359,0,-0.1980421))

Scenario_11_results <- do.call(rbind, lapply(Scenario_11, function(x) as.data.frame(t(x))))

#Scenario 12
Scenario_12 <- MAP_only(trt.effect <- c(0.9808293,0,0.8472979,0,0.8109302),
                            mu.cur <- c(0.05792359,0.1569446,0,0.08843417,-0.05792359), 
                            mu.history <- c(0.05792359,0.1569446,0,0.08843417,-0.05792359))

Scenario_12_results <- do.call(rbind, lapply(Scenario_12, function(x) as.data.frame(t(x))))

#############################################################################################
addWorksheet(MAP_only_results_3, "Sheet1")
writeData(MAP_only_results_3, "Sheet1", x = NULL)
writeData(MAP_only_results_3, sheet = "Sheet1", c(Scenario_9_results,Scenario_10_results,
                                                      Scenario_11_results,Scenario_12_results))
saveWorkbook(MAP_only_results_3, "MAP_only_3_output.xlsx", overwrite = TRUE)
##############################################################################################

#Scenario 13
Scenario_13 <- MAP_only(trt.effect <- c(0,0,0,0,0),
                            mu.cur <- c(0.1980421,0.1980421,0.1980421,0.1980421,0.1980421), 
                            mu.history <- c(0.2478002,0.2478002,0.2478002,0.2478002,0.2478002))

Scenario_13_results <- do.call(rbind, lapply(Scenario_13, function(x) as.data.frame(t(x))))

#Scenario 14
Scenario_14 <- MAP_only(trt.effect <- c(0,0,0,0,0),
                            mu.cur <- c(0.1980421,0.1980421,0.1980421,0.1980421,0.1980421), 
                            mu.history <- c(0.1569446,0.1569446,0.1569446,0.1569446,0.1569446))

Scenario_14_results <- do.call(rbind, lapply(Scenario_14, function(x) as.data.frame(t(x))))

#Scenario 15
Scenario_15 <- MAP_only(trt.effect <- c(0,0,0,0,0),
                            mu.cur <- c(0.1980421,0.1210426,0.05792359,0,-0.05792359),
                            mu.history <- c(0.2478002,0.1569446,0.08843417,0.02866724,-0.02866724))

Scenario_15_results <- do.call(rbind, lapply(Scenario_15, function(x) as.data.frame(t(x))))

#Scenario 16
Scenario_16 <- MAP_only(trt.effect <- c(0,0,0,0,0),
                            mu.cur <- c(0.1980421,0.1210426,0.05792359,0,-0.05792359), 
                            mu.history <- c(0.1569446,0.08843417,0.02866724,-0.02866724,-0.08843417))

Scenario_16_results <- do.call(rbind, lapply(Scenario_16, function(x) as.data.frame(t(x))))
addWorksheet(MAP_only_results, "Sheet1")
writeData(MAP_only_results, "Sheet1", x = NULL)
writeData(MAP_only_results, sheet = "Sheet1", c(Scenario_1_results,Scenario_2_results,Scenario_3_results,
                                                    Scenario_4_results,Scenario_5_results,Scenario_6_results,
                                                    Scenario_7_results,Scenario_8_results,Scenario_9_results,
                                                    Scenario_10_results,Scenario_11_results,Scenario_12_results,
                                                    Scenario_13_results,Scenario_14_results,Scenario_15_results,
                                                    Scenario_16_results))
saveWorkbook(MAP_only_results, "MAP_only_output.xlsx", overwrite = TRUE)




