#22/5/2024
#updated 24/1/2025

# Summarising the CRP look up results 
#calculating the proportion of the MR studies that show evidence of a causal effect by the different estimators. 
#estimating the effect of CRP on type 2 diabetes and how much it is changed by steiger filtering the results against variables that have evidence of being confounders


library(readxl)
library(tidyverse)
library(TwoSampleMR)


##MR of CRP to Type 2 diabetes


t2d <- 'ebi-a-GCST90018926' ##type 2 diabetes

crp_exp_dat <- extract_instruments(outcomes='ebi-a-GCST90029070')
crp_exp_dat <- clump_data(crp_exp_dat) 

crp_out_t2d <- extract_outcome_data(
  snps = crp_exp_dat$SNP,
  outcomes = t2d
)

dat_crp_t2d <- harmonise_data(
  exposure_dat = crp_exp_dat, 
  outcome_dat = crp_out_t2d
)


res_crp_t2d <- mr(dat_crp_t2d, method_list = c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", 
                                               "mr_egger_regression"))
MR1 <- res_crp_t2d
MR1
p1 <- mr_scatter_plot(res_crp_t2d, dat_crp_t2d)

#######
#summarising look up results - how many MR's are significant etc 

path <- "CRPlookup_MR_DR.xlsx"

tabs <- path %>%
  excel_sheets()

df.names <- paste("tab" , 1:16 ,sep="_")

for(i in 1:16){
 df <- read_excel("CRPlookup_MR_DR.xlsx", sheet = i) 
 df <- df %>%
   mutate_all(funs(str_replace(., "NA","")))
 names(df) <- tolower(names(df))
 assign(df.names[i], df)
}


df <- bind_rows(tab_1, tab_2, tab_3, tab_4, tab_5, tab_6, tab_7, tab_8, tab_9, tab_10, tab_11, tab_12, tab_13, tab_14, tab_15, tab_16)


df <- df %>% 
  drop_na('exposure')


Pvalue_ivw <- df$pvalue_ivw
Pvalue_ivw <- na.omit(Pvalue_ivw)
length(Pvalue_ivw)
mean(as.numeric(Pvalue_ivw < 0.05))

Pvalue_egger <- df$pvalue_egger
Pvalue_egger <- na.omit(Pvalue_egger)
length(Pvalue_egger)
mean(as.numeric(Pvalue_egger < 0.05))

Pvalue_med <- df$pvalue_med
Pvalue_med <- na.omit(Pvalue_med)
length(Pvalue_med)
mean(as.numeric(Pvalue_med < 0.05))

Pvalue_mode <- df$pvalue_mode
Pvalue_mode <- na.omit(Pvalue_mode)
length(Pvalue_mode)
mean(as.numeric(Pvalue_mode < 0.05))

Pvalue_lasso <- df$pvalue_lasso
Pvalue_lasso <- na.omit(Pvalue_lasso)
length(Pvalue_lasso)
mean(as.numeric(Pvalue_lasso < 0.05))

##generating the traits with an MR IVW/wald ratio p value <0.01 to filter against in MR

top_traits <- filter(df, as.numeric(df$pvalue_ivw) <= 0.01)
names(top_traits)<-make.names(names(top_traits),unique = TRUE)

#remove the CRP traits and trait with no valid id 
top_traits <- top_traits %>%
  filter(id.exposure != "ebi-a-GCST90014002") %>%
  filter(id.exposure != "prot-a-670") %>% 
  filter(id.exposure != "Childhood_bmi_1.5years") %>%
  filter(id.exposure != "Childhood_bmi_1year") %>%
  filter(id.exposure != "Childhood_bmi_3months") %>%
  filter(id.exposure != "Childhood_bmi_7years") %>%
  filter(id.exposure != "Childhood_bmi_8months") %>%
  filter(id.exposure != "Childhood_bmi_5years") %>%
  filter(id.exposure != "Childhood_bmi_6weeks") %>%
  filter(id.exposure != "GSCAN_CigDay") %>%
  filter(id.exposure != "GSCAN_SmkInit") %>% 
  filter(id.exposure != "logTG_noUKB") %>%
 filter(id.exposure != "LDL_noUKB")

#save the results
write_csv(top_traits, "data_output/topidentified_traits.csv")
length(top_traits$id.exposure)

##Run MR of the top traits against type 2 diabetes

t2d_results <- NULL

for(i in 1:length(top_traits$id.exposure)){
  exposure <- top_traits$id.exposure[i]
  exp_dat <- extract_instruments(outcomes=exposure)
  
  out_dat <- extract_outcome_data(
    snps = exp_dat$SNP,
    outcomes = t2d
  )
  
  if(!is.null(out_dat)){
  dat_for_MR <- harmonise_data(
    exposure_dat = exp_dat, 
    outcome_dat = out_dat
  )
  
  res <- mr(dat_for_MR, method_list = c("mr_wald_ratio", "mr_ivw"))
  t2d_results <- rbind(t2d_results, res)
}
}


#select those traits that also show evidence of an effect on t2d, also using p-value < 0.01

traits_to_filter <- t2d_results %>% 
  filter(method == "Inverse variance weighted" | method == "Wald ratio") %>%
  filter(pval<=0.01) %>%
  unique()

###Steiger filter the CRP -> t2d results against each of the traits that shows evidence of an effect on CRP and t2d

#identify SNPs that explain more variation in the potential confounder than CRP

results_withsteiger <- NULL 

t2d_outcome_dat <- extract_outcome_data(
  snps = crp_exp_dat$SNP,
  outcomes = t2d
)


for(i in 1:length(traits_to_filter$id.exposure)){
  
  out_dat <- extract_outcome_data(
    snps = crp_exp_dat$SNP,
    outcomes = traits_to_filter$id.exposure[i]
  )

  analysis_dat <- harmonise_data(
    exposure_dat = crp_exp_dat, 
    outcome_dat = out_dat
  )
  
  dat_steiger <- steiger_filtering(analysis_dat)
  
  snps_toremove <- dat_steiger %>%
    #filter(steiger_dir == FALSE)
    filter(steiger_dir == FALSE & steiger_pval < 0.05)

  aa <- as.list(snps_toremove$SNP)
  crp_exp_dat_st <- crp_exp_dat[ ! crp_exp_dat$SNP %in% aa, ]
  t2d_out_dat_st <- t2d_outcome_dat[ ! t2d_outcome_dat$SNP %in% aa, ]

  
  dat_for_MR <- harmonise_data(
    exposure_dat = crp_exp_dat_st, 
    outcome_dat =  t2d_out_dat_st
  )
  
  res <- mr(dat_for_MR, method_list = c("mr_wald_ratio", "mr_ivw", "mr_weighted_median", "mr_weighted_mode", 
                                        "mr_egger_regression"))
  results_withsteiger <- rbind(results_withsteiger, cbind(res, traits_to_filter$id.exposure[i], length(snps_toremove$SNP)))
  
}


write_csv(t2d_results, "data_output/t2d_MR_results.csv")

results_withsteiger <- results_withsteiger %>%
              rename(id.confounder = `traits_to_filter$id.exposure[i]`)

exposure_list <- traits_to_filter[,c("id.exposure", "exposure", 'nsnp')]
exposure_list <- exposure_list %>%
              rename(id.confounder = id.exposure) %>%
              rename(confounder = exposure) %>%
              rename(nsnp.confounder = nsnp)


results <- full_join(results_withsteiger, exposure_list)

write_csv(results, "data_output/steiger_filtered_results.csv")


