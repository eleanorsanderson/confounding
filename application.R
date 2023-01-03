#load the package (you will need to install it before you use it the first time 
#see: https://mrcieu.github.io/TwoSampleMR/index.html)
library(TwoSampleMR)

#select the SNPs associated with age at menarche and participation in strenuous sports.
aam_exp_dat <- extract_instruments(outcomes='ukb-b-3768')
ss_exp_dat <- extract_instruments(outcomes='ukb-b-7663')


#'clump' the data to remove correlated SNPs
aam_exp_dat <- clump_data(aam_exp_dat)
ss_exp_dat <- clump_data(ss_exp_dat)


t2d_out_dat_aam <- extract_outcome_data(
 snps = aam_exp_dat$SNP,
outcomes = 'finn-b-E4_DM2'
)


t2d_out_dat_ss <- extract_outcome_data(
  snps = ss_exp_dat$SNP,
 outcomes = 'finn-b-E4_DM2'
)


#harmonise the data for the exposure and outcome so the same allele is the refrence allele
dat_aam <- harmonise_data(
  exposure_dat = aam_exp_dat, 
  outcome_dat = t2d_out_dat_aam
)

dat_ss <- harmonise_data(
  exposure_dat = ss_exp_dat, 
  outcome_dat = t2d_out_dat_ss
)

#run the main and robust MR analyses
res_aam <- mr(dat_aam)
res_aam
p1_aam <- mr_scatter_plot(res_aam, dat_aam)


#run the main and robust MR analyses
res_ss <- mr(dat_ss)
res_ss
p1_ss <- mr_scatter_plot(res_ss, dat_ss)


###MVMR including BMI

id_exposure_aam <- c('ukb-b-3768', 'ieu-b-40')
id_exposure_ss <- c('ukb-b-7663', 'ieu-b-40')
id_outcome <- 'finn-b-E4_DM2'

exposure_dat_aam <- mv_extract_exposures(id_exposure_aam)
outcome_dat_aam <- extract_outcome_data(exposure_dat_aam$SNP, id_outcome)
mvdat_aam <- mv_harmonise_data(exposure_dat_aam, outcome_dat_aam)
res_aam <- mv_multiple(mvdat_aam)
res_aam

exposure_dat_ss <- mv_extract_exposures(id_exposure_ss)
outcome_dat_ss <- extract_outcome_data(exposure_dat_ss$SNP, id_outcome)
mvdat_ss <- mv_harmonise_data(exposure_dat_ss, outcome_dat_ss)
res_ss <- mv_multiple(mvdat_ss)
res_ss


###filtering to remove SNPs that explain more variation in BMI than the exposure 

#get a list of snps that explain more variation in bmi than aam 
bmi_out_dat_aam <- extract_outcome_data(
  snps = aam_exp_dat$SNP,
  outcomes = 'ieu-b-40'
)

dat_aam_bmi <- harmonise_data(
  exposure_dat = aam_exp_dat, 
  outcome_dat = bmi_out_dat_aam
)

dat_steiger <- steiger_filtering(dat_aam_bmi) 

snps_filtered <- dat_steiger %>% 
  select(SNP, steiger_dir)

snps_toremove <- snps_filtered %>%
            filter(steiger_dir == FALSE)

#remove these snps from the list of aam exposures 
aa <- as.list(snps_toremove$SNP)
aam_exp_dat_st <- aam_exp_dat[ ! aam_exp_dat$SNP %in% aa, ]


#get a list of snps that explain more variation in bmi than ss
bmi_out_dat_ss <- extract_outcome_data(
  snps = ss_exp_dat$SNP,
  outcomes = 'ieu-b-40'
)

dat_ss_bmi <- harmonise_data(
  exposure_dat = ss_exp_dat, 
  outcome_dat = bmi_out_dat_ss
)

dat_steiger_ss <- steiger_filtering(dat_ss_bmi) 

snps_filtered_ss <- dat_steiger_ss %>% 
  select(SNP, steiger_dir)

snps_toremove_ss <- snps_filtered_ss %>%
  filter(steiger_dir == FALSE)

#remove these snps from the list of aam exposures 
ass <- as.list(snps_toremove_ss$SNP)
ss_exp_dat_st <- ss_exp_dat[ ! ss_exp_dat$SNP %in% ass, ]

####################################################
##re estimate the MR with the steiger filtered data

t2d_out_dat_aam <- extract_outcome_data(
  snps = aam_exp_dat_st$SNP,
  outcomes = 'finn-b-E4_DM2'
)


t2d_out_dat_ss <- extract_outcome_data(
  snps = ss_exp_dat_st$SNP,
  outcomes = 'finn-b-E4_DM2'
)

#harmonise the data for the exposure and outcome so the same allele is the refrence allele
dat_aam <- harmonise_data(
  exposure_dat = aam_exp_dat, 
  outcome_dat = t2d_out_dat_aam
)

dat_ss <- harmonise_data(
  exposure_dat = ss_exp_dat, 
  outcome_dat = t2d_out_dat_ss
)

#run the main and robust MR analyses
res_aam_st <- mr(dat_aam)
res_aam_st

#run the main and robust MR analyses
res_ss_st <- mr(dat_ss)
res_ss_st




