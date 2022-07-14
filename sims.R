
set.seed(100)
library(dplyr)
library(tidyverse)
library(TwoSampleMR)
source('functions.R')

reps = 1

results = data.frame()
mvmrres = data.frame()
res_ivw = data.frame()
res_mregger = data.frame()
res_wmedian = data.frame()
res_wmode = data.frame()
results_all = NULL

total_snps = 400            #SNPS for X1 and X2 combined
pct = 0.4                   # proportion of snps for X2; no of SNPs for X3 matches the number for X2.
nobs = 100000
beta3 = 0
 
for (pcut in c(5e-4, 5e-6, 5e-8, 5e-12)) {

    for(j in 1:reps){
      
      results[j,"sample size"] <- nobs
      results[j,"pct_c"] <- pct
      results[j,"pvalue"] <- pcut
      
      snps <- round((1-pct)*total_snps)
      snpsc <-round(pct*total_snps)
      
      dat <- datagenA(snps, snpsc, nobs,0,0.4,beta3)
      #(no of snps, snps for confounding var, samplesize, beta1, beta2, beta3)
      
      ##############################
      ##Observational confounded association
      ##############################
      
      results[j,"obs_X1"] <- summary(lm(Y~X1,data = dat))$coefficients["X1","Estimate"]

      #############################
      ##Generate GWAS data 
      #############################  
      
      MR_dat <- GWASres(dat)
         
      #####################################
      ##Snps associated with each exposure
      #####################################
      
      #select snps associated with the exposure
  
          MR_dat.X <- MR_dat[MR_dat$X1_p < pcut,]
          results[j,"X1_totalR2"] <- sum(MR_dat.X$X1_r2)            #proportion of X1 explained by snps significantly associated with x1
          results[j, "mr.nsnps"] <- length(MR_dat.X$X1_b)
          
          MR_dat.X$both <- as.numeric(MR_dat.X$X2_p < pcut)
          results[j,"overlap"] <- sum(MR_dat.X$both)
          results[j,"prop_X2"] <- results[j,"overlap"]/results[j,"mr.nsnps"]
         
          MR_dat.Xboth <- MR_dat.X[MR_dat.X$X2_p < pcut,]
          results[j,"X1_R2_X2snps"] <- sum(MR_dat.Xboth$X1_r2)  #proportion of variation in X1 explained by SNPs associated with X2
          
          MR_dat.X2 <- MR_dat[MR_dat$X2_p < pcut,]
          results[j,"X2_totalR2"] <- sum(MR_dat.X2$X2_r2)        #proportion of X2 explained by snps significantly associated with x2
  
          
          #####################################
          ##MR ESTIMATION - using MR base
          #####################################
          
          #set up the data to use the MR package
          
          exposure_dat <- MR_dat.X[c("X1_b", "X1_se", "X1_p")] 
          exposure_dat <- exposure_dat %>% 
            rename("beta.exposure" = "X1_b") %>%
            rename("se.exposure" = "X1_se") %>%
            rename("pval.exposure" = "X1_p") %>%
            mutate(SNP = row_number()) %>%
            mutate(effect_allele.exposure = "C") %>%
            mutate(other_allele.exposure = "G") %>%
            mutate(eaf.exposure = 0.35) %>%
            mutate(exposure = "X1") %>%
            mutate(id.exposure = "X1")
          
          outcome_dat <- MR_dat.X[c("Y_b", "Y_se", "Y_p")] 
          outcome_dat <- outcome_dat %>% 
            rename("beta.outcome" = "Y_b") %>%
            rename("se.outcome" = "Y_se") %>%
            rename("pval.outcome" = "Y_p") %>%
            mutate(SNP = row_number()) %>%
            mutate(effect_allele.outcome = "C") %>%
            mutate(other_allele.outcome = "G") %>%
            mutate(eaf.outcome = 0.35) %>%
            mutate(outcome = "Y") %>%
            mutate(id.outcome = "Y")
          
           datforMR <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat, action = 1)
           mrres <- mr(datforMR, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_mode", "mr_weighted_median"))  
           
           temp <- subset(mrres, mrres$method == "Inverse variance weighted")
           res_ivw[j,"beta_ivw"] <- temp[, "b"]
           res_ivw[j,"se_ivw"] <- temp[, "se"]
           res_ivw[j,"pval_ivw"] <- temp[, "pval"]
           
           temp <- subset(mrres, mrres$method == "MR Egger")
           res_mregger[j,"beta_egger"] <- temp[, "b"]
           res_mregger[j,"se_egger"] <- temp[, "se"]
           res_mregger[j,"pval_egger"] <- temp[, "pval"]
           
           temp <- subset(mrres, mrres$method == "Weighted median")
           res_wmedian[j,"beta_med"] <- temp[, "b"]
           res_wmedian[j,"se_med"] <- temp[, "se"]
           res_wmedian[j,"pval_med"] <- temp[, "pval"]
           
           temp <- subset(mrres, mrres$method == "Weighted mode")
           res_wmode[j,"beta_mode"] <- temp[, "b"]
           res_wmode[j,"se_mode"] <- temp[, "se"]
           res_wmode[j,"pval_mode"] <- temp[, "pval"]
          
          #####################################
          ##MVMR ESTIMATION
          #####################################
          
          mvmrres[j,] <- MVMR_conf(MR_dat)
      
    }
  #combine all the results so that the different values of the cut off can be combined into the same dataset.

  results_rep <- bind_cols(results, res_ivw, res_mregger, res_wmedian, res_wmode, .name_repair = "universal")

  results_all <- rbind(results_all, results_rep)
  
    }
    

results_all
#save(results_all, file="sim_out.Rda")

  
