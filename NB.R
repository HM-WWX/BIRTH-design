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

No_borrowing <- function(trt.effect_f,mu.cur_f,mu.history_f){

# 设置并行集群，使用3个核
num_cores <- 60
cl <- makeSOCKcluster(num_cores)
registerDoSNOW(cl)

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
subtrial_trt_p <- c()
subtrial_con_p <- c()
subtrial_trt_p_95CI <- c()
subtrial_con_p_95CI <- c()
subtrial_trt_95CI <- list()
subtrial_con_95CI <- list()
power_R <- list()
subtrial_trt <- list()
subtrial_con <- list()
power_counter <- c(0,0,0,0,0)

results <- foreach(n = 1:ssize,.combine= "c", .packages = c('psrwe', "do","dplyr",'MASS', 'boot',
                            "lattice","tidyverse","patchwork","R2jags",
                            "scales","invgamma","stringr","rstan","tictoc",
                            'coda', 'rjags', 'RBesT', 'abind',"doSNOW","Hmisc"), .options.snow = opts) %dopar% {
  init_values <- list(
  .RNG.name = "base::Mersenne-Twister",  # 指定随机数生成器
  .RNG.seed =  task_seeds[n]                     # 指定种子
  )
  power_counter <- c(0,0,0,0,0)
  
  tryCatch({   
  
   for (i in 1: n.subtrial){
    setwd("//home//wwx")
    
    set.seed(task_seeds[n])
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
          logit_P = xb + study.effect + rnorm(1,0,0.25)
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
          logit_P = xb + trt.effect + rnorm(1,0,0.25)
          P = boot::inv.logit(logit_P) 
          Y[i, ] = rbinom(1,1,P)
        }
        
        return(Y)
      }
      # Current 
      n.con.samples <- con.numbers[i] #number of patients in control
      n.trt.samples <- trt.numbers[i]
      #Generate covariates X
      p <- 9
      rho <- 0.1 #correlations of the covariates
      beta_0 <- rep(1,p)
      
      ##### PSMAP simulation scenario #####
      mu.hist <- mu.history[i] 
      mu.cur.con <- mu.cur[i] # p=9, muc=0.5
      
      S0 = diag(x = rep(0.25^2, p))
      S0[which(S0 == 0)] = rho * (0.25^2)
      
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

      ###########################
      ###                     ###
      ###  DAW-MAP vs PS-MAP  ###
      ###                     ###
      ###########################
      ##########################################################################   
      
  }
  
  resp_trt <- sapply(Y.cur.treat, FUN = sum)
  resp_con <- sapply(Y.cur.control, FUN = sum)
  kMod <- length(resp_trt)
  kMod <- as.numeric(kMod)
  for (j in 1:kMod){
    MAP_model_CPP_trt <- rbeta(1000,0.5+resp_trt[j],0.5+trt.numbers[j]-resp_trt[j]) 
    MAP_model_CPP_con <- rbeta(1000,0.5+resp_con[j],0.5+con.numbers[j]-resp_con[j])
    
    subtrial_trt_p[j]<-mean(MAP_model_CPP_trt)
    subtrial_trt_p_95CI[j] <- smean.cl.normal(MAP_model_CPP_trt)[3]-smean.cl.normal(MAP_model_CPP_trt)[2]
    
    subtrial_con_p[j]<-mean(MAP_model_CPP_con)
    subtrial_con_p_95CI[j] <- smean.cl.normal(MAP_model_CPP_con)[3]-smean.cl.normal(MAP_model_CPP_con)[2]
    
    
    trt.response <- as.matrix(MAP_model_CPP_trt)
    con.response <- as.matrix(MAP_model_CPP_con)
    trt.difference <- trt.response - con.response
    power <- sum(trt.difference > 0)
    
    if (power >= 950){
      power_counter[j] <- power_counter[j]+1
    }else{
      power_counter[j] <- power_counter[j]
      
    }
  
  }
  
  list(subtrial_trt = subtrial_trt_p, subtrial_con = subtrial_con_p, power_R = power_counter, 
       subtrial_trt_95CI=subtrial_trt_p_95CI, subtrial_con_95CI=subtrial_con_p_95CI)
  
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
 
 trial_numbers <- length(results)/5
 
 for (i in 0:(trial_numbers-1)){
   
   subtrial_trt_all <- rbind(subtrial_trt_all,results[[1+5*i]])  
   subtrial_con_all <- rbind(subtrial_con_all,results[[2+5*i]])
   Power_all <-rbind(Power_all,results[[3+5*i]])
   subtrial_trt_95CI_all <- rbind(subtrial_trt_95CI_all,results[[4+5*i]])
   subtrial_con_95CI_all <- rbind(subtrial_con_95CI_all,results[[5+5*i]])
 }
 
 
 subtrial_trt_all <- subtrial_trt_all[-1,]
 subtrial_con_all <- subtrial_con_all[-1,]
 Power_all <- Power_all[-1,]
 subtrial_trt_95CI_all <- subtrial_trt_95CI_all[-1,]
 subtrial_con_95CI_all <- subtrial_con_95CI_all[-1,]
 
 
aa1 <- mean(subtrial_con_all[,1])-0
bb1 <- mean(subtrial_con_all[,2])-0
cc1 <- mean(subtrial_con_all[,3])-0
dd1 <- mean(subtrial_con_all[,4])-0
ee1 <- mean(subtrial_con_all[,5])-0
 
subtrial_con_mean <- c(aa1,bb1,cc1,dd1,ee1)

aa2 <- mean(subtrial_trt_all[,1])-0
bb2 <- mean(subtrial_trt_all[,2])-0
cc2 <- mean(subtrial_trt_all[,3])-0
dd2 <- mean(subtrial_trt_all[,4])-0
ee2 <- mean(subtrial_trt_all[,5])-0

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

results_final = list(subtrial_con_mean, subtrial_trt_mean,subtrial_con_mse,subtrial_trt_95CI,subtrial_con_95CI,subtrial_Power)
return(results_final)
 
 #NO BORROWING
 close(pb)
 stopImplicitCluster()
 stopCluster(cl)
 
}

No_borrowing_results <- createWorkbook()
No_borrowing_results_1 <- createWorkbook()
No_borrowing_results_2 <- createWorkbook()
No_borrowing_results_3 <- createWorkbook()
#Scenario 1
Scenario_1 <- No_borrowing(trt.effect <- c(1.386294,1.386294,1.386294,1.386294,1.386294),
                           mu.cur <- c(0,0,0,0,0), 
                           mu.history <- c(0,0,0,0,0))

Scenario_1_results <- do.call(rbind, lapply(Scenario_1, function(x) as.data.frame(t(x))))

#Scenario 2
Scenario_2 <- No_borrowing(trt.effect <- c(0.9808293,0.9808293,0.9808293,0.9808293,0.9808293),
                           mu.cur <- c(0.05792359,0.05792359,0.05792359,0.05792359,0.05792359),
                           mu.history <- c(0.05792359,0.05792359,0.05792359,0.05792359,0.05792359))

Scenario_2_results <- do.call(rbind, lapply(Scenario_2, function(x) as.data.frame(t(x))))

#Scenario 3
Scenario_3 <- No_borrowing(trt.effect <- c(0.5389965,0.5389965,0.5389965,0.5389965,0.5389965),
                           mu.cur <- c(0.1210426,0.1210426,0.1210426,0.1210426,0.1210426), 
                           mu.history <- c(0.1210426,0.1210426,0.1210426,0.1210426,0.1210426))

Scenario_3_results <- do.call(rbind, lapply(Scenario_3, function(x) as.data.frame(t(x))))

#Scenario 4
Scenario_4 <- No_borrowing(trt.effect <- c(0,0,0,0,0),
                           mu.cur <- c(0.1980421,0.1980421,0.1980421,0.1980421,0.1980421), 
                           mu.history <- c(0.1980421,0.1980421,0.1980421,0.1980421,0.1980421))

Scenario_4_results <- do.call(rbind, lapply(Scenario_4, function(x) as.data.frame(t(x))))

#############################################################################################
addWorksheet(No_borrowing_results_1, "Sheet1")
writeData(No_borrowing_results_1, "Sheet1", x = NULL)
writeData(No_borrowing_results_1, sheet = "Sheet1", c(Scenario_1_results,Scenario_2_results,
                                                      Scenario_3_results,Scenario_4_results))
saveWorkbook(No_borrowing_results_1, "No_borrowing_1_output.xlsx", overwrite = TRUE)
##############################################################################################

#Scenario 5
Scenario_5 <- No_borrowing(trt.effect <- c(0.9808293,0.8979416,0.8472979,0.8197099,0.8109302),
                           mu.cur <- c(0.05792359,0.02866724,0,0.02866724,-0.05792359), 
                           mu.history <- c(0.05792359,0.02866724,0,0.02866724,-0.05792359))

Scenario_5_results <- do.call(rbind, lapply(Scenario_5, function(x) as.data.frame(t(x))))

#Scenario 6
Scenario_6 <- No_borrowing(trt.effect <- c(0.9808293,0.8472979,0.8109302,0.8472979,0.9808293),
                           mu.cur <- c(0.05792359,0,-0.05792359,-0.1210426,-0.1980421),
                           mu.history <- c(0.05792359,0,-0.05792359,-0.1210426,-0.1980421))

Scenario_6_results <- do.call(rbind, lapply(Scenario_6, function(x) as.data.frame(t(x))))

#Scenario 7
Scenario_7 <- No_borrowing(trt.effect <- c(0,0,0,0,0),
                           mu.cur <- c(0.1980421,0.1569446,0.1210426,0.08843417,0.05792359), 
                           mu.history <- c(0.1980421,0.1569446,0.1210426,0.08843417,0.05792359))

Scenario_7_results <- do.call(rbind, lapply(Scenario_7, function(x) as.data.frame(t(x))))

#Scenario 8
Scenario_8 <- No_borrowing(trt.effect <- c(0,0,0,0,0),
                           mu.cur <- c(0.1980421,0.1210426,0.05792359,0,-0.05792359), 
                           mu.history <- c(0.1980421,0.1210426,0.05792359,0,-0.05792359))

Scenario_8_results <- do.call(rbind, lapply(Scenario_8, function(x) as.data.frame(t(x))))

#############################################################################################
addWorksheet(No_borrowing_results_2, "Sheet1")
writeData(No_borrowing_results_2, "Sheet1", x = NULL)
writeData(No_borrowing_results_2, sheet = "Sheet1", c(Scenario_5_results,Scenario_6_results,
                                                      Scenario_7_results,Scenario_8_results))
saveWorkbook(No_borrowing_results_2, "No_borrowing_2_output.xlsx", overwrite = TRUE)
##############################################################################################

#Scenario 9
Scenario_9 <- No_borrowing(trt.effect <- c(0.9808293,0.9808293,0.9808293,0.9808293,0),
                           mu.cur <- c(0.05792359,0.05792359,0.05792359,0.05792359,0.1980421), 
                           mu.history <- c(0.05792359,0.05792359,0.05792359,0.05792359,0.1980421))

Scenario_9_results <- do.call(rbind, lapply(Scenario_9, function(x) as.data.frame(t(x))))

#Scenario 10
Scenario_10 <- No_borrowing(trt.effect <- c(0.9808293,0.9808293,0.9808293,0,0),
                            mu.cur <- c(0.05792359,0.05792359,0.05792359,0.1980421,0.1980421), 
                            mu.history <- c(0.05792359,0.05792359,0.05792359,0.19804219,0.1980421))

Scenario_10_results <- do.call(rbind, lapply(Scenario_10, function(x) as.data.frame(t(x))))

#Scenario 11
Scenario_11 <- No_borrowing(trt.effect <- c(0.9808293,0,0.8109302,0,0.9808293),
                            mu.cur <- c(0.05792359,0.1210426,-0.05792359,0,-0.1980421), 
                            mu.history <- c(0.05792359,0.1210426,-0.05792359,0,-0.1980421))

Scenario_11_results <- do.call(rbind, lapply(Scenario_11, function(x) as.data.frame(t(x))))

#Scenario 12
Scenario_12 <- No_borrowing(trt.effect <- c(0.9808293,0,0.8472979,0,0.8109302),
                            mu.cur <- c(0.05792359,0.1569446,0,0.08843417,-0.05792359), 
                            mu.history <- c(0.05792359,0.1569446,0,0.08843417,-0.05792359))

Scenario_12_results <- do.call(rbind, lapply(Scenario_12, function(x) as.data.frame(t(x))))

#############################################################################################
addWorksheet(No_borrowing_results_3, "Sheet1")
writeData(No_borrowing_results_3, "Sheet1", x = NULL)
writeData(No_borrowing_results_3, sheet = "Sheet1", c(Scenario_9_results,Scenario_10_results,
                                                      Scenario_11_results,Scenario_12_results))
saveWorkbook(No_borrowing_results_3, "No_borrowing_3_output.xlsx", overwrite = TRUE)
##############################################################################################

#Scenario 13
Scenario_13 <- No_borrowing(trt.effect <- c(0,0,0,0,0),
                            mu.cur <- c(0.1980421,0.1980421,0.1980421,0.1980421,0.1980421), 
                            mu.history <- c(0.2478002,0.2478002,0.2478002,0.2478002,0.2478002))

Scenario_13_results <- do.call(rbind, lapply(Scenario_13, function(x) as.data.frame(t(x))))

#Scenario 14
Scenario_14 <- No_borrowing(trt.effect <- c(0,0,0,0,0),
                            mu.cur <- c(0.1980421,0.1980421,0.1980421,0.1980421,0.1980421), 
                            mu.history <- c(0.1569446,0.1569446,0.1569446,0.1569446,0.1569446))

Scenario_14_results <- do.call(rbind, lapply(Scenario_14, function(x) as.data.frame(t(x))))

#Scenario 15
Scenario_15 <- No_borrowing(trt.effect <- c(0,0,0,0,0),
                            mu.cur <- c(0.1980421,0.1210426,0.05792359,0,-0.05792359),
                            mu.history <- c(0.2478002,0.1569446,0.08843417,0.02866724,-0.02866724))

Scenario_15_results <- do.call(rbind, lapply(Scenario_15, function(x) as.data.frame(t(x))))

#Scenario 16
Scenario_16 <- No_borrowing(trt.effect <- c(0,0,0,0,0),
                            mu.cur <- c(0.1980421,0.1210426,0.05792359,0,-0.05792359), 
                            mu.history <- c(0.1569446,0.08843417,0.02866724,-0.02866724,-0.08843417))

Scenario_16_results <- do.call(rbind, lapply(Scenario_16, function(x) as.data.frame(t(x))))

addWorksheet(No_borrowing_results, "Sheet1")
writeData(No_borrowing_results, "Sheet1", x = NULL)
writeData(No_borrowing_results, sheet = "Sheet1", c(Scenario_1_results,Scenario_2_results,Scenario_3_results,
                                                    Scenario_4_results,Scenario_5_results,Scenario_6_results,
                                                    Scenario_7_results,Scenario_8_results,Scenario_9_results,
                                                    Scenario_10_results,Scenario_11_results,Scenario_12_results,
                                                    Scenario_13_results,Scenario_14_results,Scenario_15_results,
                                                    Scenario_16_results))
saveWorkbook(No_borrowing_results, "No_borrowing_output.xlsx", overwrite = TRUE)




