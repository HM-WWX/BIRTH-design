rm(list=ls())
library(psrwe)
library(RBesT)
library(coda)
library(rjags)
library(R2jags)
library(tictoc)
library(abind)
library(boot)
library(rstan)
library(MASS)
library(lattice)
library(tidyverse)
library(dplyr)
library(doParallel)
library(do)
library(snow)
library(snowfall)
library(patchwork) 
library(PSW)
library(abind)
library(twang)
library(stringr)
library(ggthemr)
library(patchwork)
library(ggplot2)
library(gglabeller)
library(ggrepel)
library(ggannotate)
library(Seurat)
library(scales)
library(invgamma)

setwd("C:/Users/dell/Desktop/Work")

Pooled_MAP.mu.con <- c()
tau.balance <- c()
CI_lower.all <- c()
CI_upper.all <- c()

for (j in 1:60){
  
  set.seed(j+1)
  
  #######################
  ###                 ###
  ###  generate data  ###
  ###                 ###
  #######################

    X_convert = function(X, binary_col) {
      
      
      for (i in 1:length(binary_col)) {
        current_col = X[, i]
        X[which(current_col > 3), i] = 1
        X[which(current_col <= 3), i] = 0
      }
      return(t(X))
    }
    
    
    gen_Y = function(X, beta, var_eps, trt.effect = 0, study.effect = 0) {
      
      
      N = dim(X)[1]
      Y = matrix(0, N)
      for (i in 1:N) {
        xb = t(beta) %*% X[i, ]
        logit_P = xb + rnorm(1, 0, var_eps)
        P = boot::inv.logit(logit_P)
        Y[i, ] = rbinom(1,1,0.5)
      }
      
      return(Y)
    }
    
    
    #Historic
    n.hist = 1 # number of historical trials
    n.hist.total = 200 # total number of patients
    n.hist.samples <- c(200) # number of patients in each trial
    
    
    # Current 
    n.con = 30 #number of patients in control
    
    #Generate covariates X
    p <- 7
    rho <- 0.1 #correlations of the covariates
    eps.var <- 0
    beta_0 <- rep(1,p)
    trt.effect = 0
    
    ##### PSMAP simulation scenario #####
    
    mu.hist <- c(0) 
    
    mu.cur.con = 0 # p=7, muc=0.5
    
    S0 = diag(x = rep(0.5^2, p))
    S0[which(S0 == 0)] = rho * (0.5^2)
    
    S1 = diag(x = rep(0.5^2, p))
    S1[which(S1 == 0)] = rho * (0.5^2)
    
    
    hist.list = list()
    X.hist <- X.hist.init <- list()
    delta = c(0)       #study.effect
    
    #set.seed(888)
    for(i in 1:(n.hist)) {
      X.hist.init = MASS::mvrnorm(n.hist.samples[i], mu = rep(mu.hist[i], p), Sigma = S0)
      X.hist[[i]] = t(X_convert(X.hist.init, c(1, 2)))
      hist.list[[i]] = gen_Y(X.hist[[i]], beta = beta_0, var_eps = eps.var, study.effect = delta[i])
    }
    
    X.hist.combine = abind(X.hist, along=1)
    Y.hist.combine = abind(hist.list, along=1)
    hist.ybar = sapply(hist.list, FUN = sum)
    
    ### Current 
    X.cur.con.init = MASS::mvrnorm(n.con, mu = rep(mu.cur.con, p), Sigma = S1)
    X.cur.con = t(X_convert(X.cur.con.init, c(1, 2)))
    
    #generate outcomes for control
    y.control = gen_Y(X.cur.con, beta = beta_0, var_eps = eps.var, trt.effect = 0)
    ybar.control = sum(y.control)
    
    Y.cur.con=y.control
    S = 1
    Y.combine = rbind(Y.hist.combine, Y.cur.con) 
    dta = rbind(X.hist.combine, X.cur.con)
    covs = NULL
    
    for (k in 1:p) {
      if(k<10){
        covs[k] = paste("V", k, sep = "")
      }else{
        covs[1:9] = paste("V0", 1:9, sep = "")
        covs[10:k] = paste("V", 10:k, sep = "")
      }
    }
    colnames(dta) = covs

    dta = as.data.frame(dta)
    
    Y.hist.combine <- sample(c(rep(1, 100), rep(0, 100)))
    y.control <- sample(c(rep(1, 15), rep(0, 15)))
    ###########################
    ###                     ###
    ###  DAW-MAP vs PS-MAP  ###
    ###                     ###
    ###########################
    
    ###########################################################################DAW-MAP###########################################################
    
    Pooled.fit = function(std_heterogeneity,Ybar.hist, 
                            niter,data.indiv) {
      
      
      # Obtain prior samples when using the PS-MAP prior with the specified target ESS
    
      ess.res = c(0)
      tau.res = c(0)
      theta = list()
    
          dataTemp = list(
            "resp" = sum(Ybar.hist),
            "n" = length(Ybar.hist)
          )
          dataTemp$std_heterogeneity <- std_heterogeneity
          
          
          model <-
            jags.model(
              file = "C:/Users/dell/Desktop/Work/Pooled.txt",
              data = dataTemp,
              n.chains = 1,
              n.adapt = 0.2 * niter,
              quiet = TRUE
            )
          update(model, n.iter = 0.2 * niter, progress.bar = 'none') # burn in
          MAP_model <-
            coda.samples(model,
                         variable.names = "p.pred",
                         thin = 20,
                         n.iter = niter, 
                         progress.bar = 'none',
                         na.rm = TRUE)
          
          theta.pred = c(MAP_model[[1]][,1])
          
        
        mix.res <- automixfit(theta.pred, Nc = 1:5, thresh = 0, type = "beta")
      
      
      #calculate the posterior given prior and current data
      posterior.indiv <- postmix(mix.res, n = length(data.indiv), r = sum(data.indiv))
      results = list( posterior.indiv,theta)
      
      Pooled_MAP_post <- results[[1]]
      mu.con_fat <- rmix(Pooled_MAP_post, n=10000)
      mu.con_p <- mu.con_fat[seq(1, length(mu.con_fat), by=10)]
      mu.con <- mean(mu.con_p)
      
      results = list(posterior.indiv,theta,mu.con_p,mu.con)
      
      return(results)
    }
    
    
    
    Pooled_MAP <- Pooled.fit(std_heterogeneity = 0.4,Ybar.hist=Y.hist.combine,
                          niter = 10000,data.indiv = y.control)  
    
    mu.con.all <- Pooled_MAP[[3]]
    mu.con_mean <- Pooled_MAP[[4]]
    
    Pooled_MAP.mu.con <- c(Pooled_MAP.mu.con,mu.con.all)
    
    
    
    # 样本标准误差
    se <- sd(mu.con.all) / sqrt(length(mu.con.all))
    # 计算95%置信区间
    error_margin <- qt(0.975, df=length(mu.con.all)-1) * se
    CI_lower <- mu.con_mean - error_margin
    CI_upper <- mu.con_mean + error_margin
    CI_lower.all <- c(CI_lower.all,CI_lower)
    CI_upper.all <- c(CI_upper.all,CI_upper)
    
    
    
    print(j)

  
}



# 样本均值
mean_mu.con <- mean(Pooled_MAP.mu.con)

MSE <- mean((Pooled_MAP.mu.con - 0.5)**2)

mean_mu.con-0.5
MSE
# 计算平均95%置信区间
mean(CI_lower.all)
mean(CI_upper.all)




