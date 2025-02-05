#load the package (you will need to install it before you use it the first time 
#see: https://mrcieu.github.io/TwoSampleMR/index.html)
library(TwoSampleMR)
library(tidyverse)
source('code/MRcorge_edit.R')

#select the SNPs associated with age at menarche
aam_exp_dat <- extract_instruments(outcomes='ukb-b-3768')


#'clump' the data to remove correlated SNPs
aam_exp_dat <- clump_data(aam_exp_dat)


t2d_out_dat_aam <- extract_outcome_data(
 snps = aam_exp_dat$SNP,
outcomes = 'ebi-a-GCST90018926'
)



#harmonise the data for the exposure and outcome so the same allele is the refrence allele
dat_aam <- harmonise_data(
  exposure_dat = aam_exp_dat, 
  outcome_dat = t2d_out_dat_aam
)

#run the main and robust MR analyses
res_aam <- mr(dat_aam, method_list = c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", 
                                       "mr_egger_regression"))
res_aam

p1_aam <- mr_scatter_plot(res_aam, dat_aam)

###MR CORGE plot

res <- mrcorge_ed(dat_aam, rank='beta', K=10, method_list = c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", 
                                                              "mr_egger_regression"))
aam_t2d_mrcorge_output <- res$result %>%
  pivot_wider(names_from = method, values_from = c(b, se, pval))

write_csv(aam_t2d_mrcorge_output, "data_output/mrcorge_aam_t2d.csv")



#IVW plot
res <- mrcorge_ed(dat_aam, rank='beta', K=10, method_list =c("mr_ivw"))
plot_mrcorge(res, scale = 'exp')
ggsave("plots/aam_t2d_ivw_corge.pdf", width = 6.5, height = 5)


#WEighted median plot
res <- mrcorge_ed(dat_aam, rank='beta', K=10, method_list =c("mr_weighted_median"))
plot_mrcorge(res, scale = 'exp')
ggsave("plots/aam_t2d_median_corge.pdf", width = 6.5, height = 5)

#weighted mode plot
res <- mrcorge_ed(dat_aam, rank='beta', K=10, method_list =c("mr_weighted_mode"))
plot_mrcorge(res, scale = 'exp')
ggsave("plots/aam_t2d_mode_corge.pdf", width = 6.5, height = 5)

#MR egger plot 
res <- mrcorge_ed(dat_aam, rank='beta', K=10, method_list =c("mr_egger_regression"))
plot_mrcorge(res, scale = 'exp')
ggsave("plots/aam_t2d_egger_corge.pdf", width = 6.5, height = 5)





###MVMR including BMI

id_exposure_aam <- c('ukb-b-3768', 'ukb-b-4650')

id_outcome <- 'ebi-a-GCST90018926'

exposure_dat_aam <- mv_extract_exposures(id_exposure_aam)
outcome_dat_aam <- extract_outcome_data(exposure_dat_aam$SNP, id_outcome)
mvdat_aam <- mv_harmonise_data(exposure_dat_aam, outcome_dat_aam)
res_aam_mv <- mv_multiple(mvdat_aam)
res_aam_mv


###filtering to remove SNPs that explain more variation in BMI than the exposure 

#get a list of snps that explain more variation in bmi than aam 
bmi_out_dat_aam <- extract_outcome_data(
  snps = aam_exp_dat$SNP,
  outcomes = 'ukb-b-4650'
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



####################################################
##re estimate the MR with the steiger filtered data

t2d_out_dat_aam <- extract_outcome_data(
  snps = aam_exp_dat_st$SNP,
  outcomes = 'ebi-a-GCST90018926'
)

#harmonise the data for the exposure and outcome so the same allele is the refrence allele
dat_aam_steiger <- harmonise_data(
  exposure_dat = aam_exp_dat, 
  outcome_dat = t2d_out_dat_aam
)

#run the main and robust MR analyses
res_aam_st <- mr(dat_aam_steiger)
res_aam_st


####Plot
#add a variable to data to mark the filtered snps
dat_aam$st <- dat_aam$SNP %in% aa

plotdat <- dat_aam %>%
  filter(mr_keep == TRUE) %>%
  select(beta.exposure, beta.outcome, st) %>%
  mutate(beta.outcome = ifelse(beta.exposure<0,-beta.outcome, beta.outcome)) %>%
  mutate(beta.exposure = abs(beta.exposure))

#plot with filtered snps highlighted
plot1 <- ggplot(plotdat, aes(x=beta.exposure)) +
  geom_point(aes(y=beta.outcome,  group=st, colour=st)) +
  geom_abline(intercept = 0, slope = -0.269) +
  scale_colour_brewer(palette = "Set1") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "bottom") +
  labs(x = "SNP-age at menarche association", y = "SNP-T2D association") +
  scale_colour_discrete("Steiger filtered")

ggsave("plots/appliedplot.pdf", width = 6.5, height = 5)



