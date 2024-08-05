

#args <- as.numeric(commandArgs(T))
#set.seed((args[1]*1000))
#job_id <- ((args[1]))
#message("job number ", job_id)
# 

library(dplyr)
library(tidyverse)
library(MASS)
source('functions_withmediation.R')
source('setupmodels_withmediation.R')

reps = 2

res_pcut = data.frame()
results = data.frame()
mvmrres = data.frame()
res_ivw = data.frame()
res_mregger = data.frame()
res_wmedian = data.frame()
res_wmode = data.frame()
res_ivw.s = data.frame()
res_mregger.s = data.frame()
res_wmedian.s = data.frame()
res_wmode.s = data.frame()
results_all = NULL
results_out = NULL
mvmrres <- NULL

for(gm in c(0.5,1)){
for(model in c('A', 'B', 'C')){
  
params <- setup(model)

snps = params[1]      
snpsc = params[2]         
nobs = params[3]
pi = gm


    for(j in 1:reps){
      
      
      dat <- datagenA(snps, snpsc, nobs,0,0.4,pi)
      #(no of snps, snps for confounding var, samplesize, beta1, beta2, snp-confounder effect)
      
      results[1,"model"] <- model
      results[1,"pi"] <- pi
      results[1,"sample.size"] <- nobs
      
      ####OLS regression of X on Y######
      
      ols <- summary(lm(Y ~ X1, data = dat))
      
      results[1,"ols_b"] <- ols$coefficients["X1","Estimate"]
      results[1,"ols_se"] <- ols$coefficients["X1","Std. Error"]
      
      #############################
      ##Generate GWAS data 
      #############################  
      
      MR_dat <- GWASres(dat)
         
      #####################################
      ##Snps associated with each exposure
      #####################################
      
      #select snps associated with the exposure
      res_pcut = NULL
      for (pcut in c(5e-4, 5e-6, 5e-8, 5e-12)) {
        
       
        results[1,"snps_x"] <- snps
        results[1,"snps_c"] <- snpsc
        results[1,"pvalue"] <- pcut
        
        results[1,"obs_X1"] <- summary(lm(Y~X1,data = dat))$coefficients["X1","Estimate"]
        
          MR_dat.X <- MR_dat[MR_dat$X1_p < pcut,]
          results[1,"X1_totalR2"] <- sum(MR_dat.X$X1_r2)            #proportion of X1 explained by snps significantly associated with x1
          results[1,"mr.nsnps"] <- length(MR_dat.X$X1_b)
          
          MR_dat.X$both <- as.numeric(MR_dat.X$X2_p < pcut)
          results[1,"overlap"] <- sum(MR_dat.X$both)
          results[1,"prop_X2"] <- results[1,"overlap"]/results[1,"mr.nsnps"]
         
          MR_dat.Xboth <- MR_dat.X[MR_dat.X$X2_p < pcut,]
          results[1,"X1_R2_X2snps"] <- sum(MR_dat.Xboth$X1_r2)  #proportion of variation in X1 explained by SNPs associated with X2
          
          MR_dat.X2 <- MR_dat[MR_dat$X2_p < pcut,]
          results[1,"X2_totalR2"] <- sum(MR_dat.X2$X2_r2)        #proportion of X2 explained by snps significantly associated with x2
  
          
          #####################################
          ##MR ESTIMATION - using code from the MR base package
          #####################################
        
           
           res.ivw <- mr_ivw(MR_dat.X$X1_b, MR_dat.X$Y_b, MR_dat.X$X1_se, MR_dat.X$Y_se)
           res_ivw[1,"beta_ivw"] <- res.ivw$b
           res_ivw[1,"se_ivw"] <- res.ivw$se
           res_ivw[1,"pval_ivw"] <- res.ivw$pval
           
           res.egger <- mr_egger_regression(MR_dat.X$X1_b, MR_dat.X$Y_b, MR_dat.X$X1_se, MR_dat.X$Y_se)
           res_mregger[1,"beta_egger"] <- res.egger$b
           res_mregger[1,"se_egger"] <- res.egger$se
           res_mregger[1,"pval_egger"] <- res.egger$pval
           
           res.wmed <- mr_weighted_median(MR_dat.X$X1_b, MR_dat.X$Y_b, MR_dat.X$X1_se, MR_dat.X$Y_se)
           res_wmedian[1,"beta_med"] <- res.wmed$b
           res_wmedian[1,"se_med"] <- res.wmed$se
           res_wmedian[1,"pval_med"] <- res.wmed$pval
           
           res.wmode <- mr_weighted_mode(MR_dat.X$X1_b, MR_dat.X$Y_b, MR_dat.X$X1_se, MR_dat.X$Y_se)
           res_wmode[1,"beta_mode"] <- res.wmode$b
           res_wmode[1,"se_mode"] <- res.wmode$se
           res_wmode[1,"pval_mode"] <- res.wmode$pval
          
          #####################################
          ##MVMR ESTIMATION
          #####################################
          
           mvmrres <- MVMR_conf(MR_dat)
           
           
           ####Add in steiger filtering and rerun the univariable MR estimation
           
           MR_dat.S <- MR_dat.X
           MR_dat.S$keep <- as.numeric(MR_dat.S$X2_p > MR_dat.S$X1_p)         #this is currently a very simple steiger filtering (between X1 and X2)- 
                                                                         
           MR_dat.S <- filter(MR_dat.S, MR_dat.S$keep==1)
           
           res.ivw <- mr_ivw(MR_dat.S$X1_b, MR_dat.S$Y_b, MR_dat.S$X1_se, MR_dat.S$Y_se)
           res_ivw.s[1,"beta_ivw.s"] <- res.ivw$b
           res_ivw.s[1,"se_ivw.s"] <- res.ivw$se
           res_ivw.s[1,"pval_ivw.s"] <- res.ivw$pval
           
           res.egger <- mr_egger_regression(MR_dat.S$X1_b, MR_dat.S$Y_b, MR_dat.S$X1_se, MR_dat.S$Y_se)
           res_mregger.s[1,"beta_egger.s"] <- res.egger$b
           res_mregger.s[1,"se_egger.s"] <- res.egger$se
           res_mregger.s[1,"pval_egger.s"] <- res.egger$pval
           
           res.wmed <- mr_weighted_median(MR_dat.S$X1_b, MR_dat.S$Y_b, MR_dat.S$X1_se, MR_dat.S$Y_se)
           res_wmedian.s[1,"beta_med.s"] <- res.wmed$b
           res_wmedian.s[1,"se_med.s"] <- res.wmed$se
           res_wmedian.s[1,"pval_med.s"] <- res.wmed$pval
           
           res.wmode <- mr_weighted_mode(MR_dat.S$X1_b, MR_dat.S$Y_b, MR_dat.S$X1_se, MR_dat.S$Y_se)
           res_wmode.s[1,"beta_mode.s"] <- res.wmode$b
           res_wmode.s[1,"se_mode.s"] <- res.wmode$se
           res_wmode.s[1,"pval_mode.s"] <- res.wmode$pval
           
           results_rep <- bind_cols(results, res_ivw, res_mregger, res_wmedian, res_wmode, mvmrres, 
                                    res_ivw.s, res_mregger.s, res_wmedian.s, res_wmode.s)
           res_pcut <- rbind(res_pcut, results_rep)
           
      
      }
    
      #combine all the results so that the different values of the cut off can be combined into the same dataset.
      results_all <- rbind(results_all, res_pcut)
    }
}
}

#results_out <- rbind(results_out, results_all)

#save(results_all, file=sprintf("results_med_%s.Rda", job_id))
