
##CRP -> T2D MR estimation

#devtools::install_github("zhwm/MRCorge")
library(tidyverse)
library(TwoSampleMR)
library(MRCorge)
source('code/MRcorge_edit.R')

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


###MR Corge and plots 
crp_t2d_mrcorge <- mrcorge_ed(dat_crp_t2d, rank='beta', K=10, method_list =c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", 
                                                                             "mr_egger_regression"))

crp_t2d_mrcorge_output <- crp_t2d_mrcorge$result %>%
  pivot_wider(names_from = method, values_from = c(b, se, pval))

write_csv(crp_t2d_mrcorge_output, "data_output/mrcorge_crp_t2d.csv")



#IVW plot
res <- mrcorge_ed(dat_crp_t2d, rank='beta', K=10, method_list =c("mr_ivw"))
plot_mrcorge(res, scale = 'exp')
ggsave("plots/crp_t2d_ivw_corge.pdf", width = 6.5, height = 5)


#WEighted median plot
res <- mrcorge_ed(dat_crp_t2d, rank='beta', K=10, method_list =c("mr_weighted_median"))
plot_mrcorge(res, scale = 'exp')
ggsave("plots/crp_t2d_median_corge.pdf", width = 6.5, height = 5)

#weighted mode plot
res <- mrcorge_ed(dat_crp_t2d, rank='beta', K=10, method_list =c("mr_weighted_mode"))
plot_mrcorge(res, scale = 'exp')
ggsave("plots/crp_t2d_mode_corge.pdf", width = 6.5, height = 5)

#MR egger plot 
res <- mrcorge_ed(dat_crp_t2d, rank='beta', K=10, method_list =c("mr_egger_regression"))
plot_mrcorge(res, scale = 'exp')
ggsave("plots/crp_t2d_egger_corge.pdf", width = 6.5, height = 5)



##write script to extract the list of SNPs removed for each of top 10 traits
#add results for top 10 traits to paper and full table to supplement 
#scatter plot for CRP on T2D that highlights the SNPs removed by the proteins (with rs429358 again 
#in a different colour so it can be identified from the other SNPs removed)

##table of MR estimates

library(readr)
results <- read_csv("data_output/steiger_filtered_results.csv")

output <- results %>% 
  pivot_wider(names_from = method, values_from = c(b, se, pval))
write_csv(output, "data_output/steiger_filtered_results_pivot.csv")

###Identifying SNPs removed from top 10 confounders

top10 <- output %>% 
  slice_max(`pval_Inverse variance weighted`,n=10)


t2d <- 'ebi-a-GCST90018926'
crp_exp_dat <- extract_instruments(outcomes='ebi-a-GCST90029070')
crp_exp_dat <- clump_data(crp_exp_dat) 

t2d_outcome_dat <- extract_outcome_data(
  snps = crp_exp_dat$SNP,
  outcomes = t2d
)




all_removed_snps <- NULL



for(i in 1:10){
  
  a <- as.character(top10[i,'id.confounder'])

  out_dat <- extract_outcome_data(
    snps = crp_exp_dat$SNP,
    outcomes = a
  )

  analysis_dat <- harmonise_data(
    exposure_dat = crp_exp_dat, 
    outcome_dat = out_dat
  )

  dat_steiger <- steiger_filtering(analysis_dat)

  snps_toremove <- dat_steiger %>%
    filter(steiger_dir == FALSE & steiger_pval < 0.05) %>%
    select(SNP, beta.exposure, beta.outcome, palindromic, ambiguous, id.outcome, outcome, pval.outcome, steiger_pval)

  all_removed_snps <- rbind(all_removed_snps, snps_toremove)
}

write_csv(all_removed_snps, "data_output/confounderSNPS.csv")


#number of unique SNPS identified 
length(unique(all_removed_snps$SNP))

#removing palindromic ambiguous SNPS

df <- all_removed_snps %>%
  filter(palindromic==FALSE & ambiguous==FALSE)
length(unique(df$SNP))


##scatter plot with top 10 confounding SNPs identified


dat_crp_t2d$filtered <- dat_crp_t2d$SNP %in% as.list(unique(df$SNP))

plotdat <- dat_crp_t2d %>%
  filter(mr_keep == TRUE) %>%
  select(beta.exposure, beta.outcome, filtered) %>%
  mutate(beta.outcome = ifelse(beta.exposure<0,-beta.outcome, beta.outcome)) %>%
  mutate(beta.exposure = abs(beta.exposure))

#plot with filtered snps highlighted
plot1 <- ggplot(plotdat, aes(x=beta.exposure)) +
  geom_point(aes(y=beta.outcome,  group=filtered, colour=filtered)) +
  scale_colour_brewer(palette = "Set1") +
  geom_abline(intercept = 0, slope = 0.1245) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "bottom") +
  labs(x = "SNP-CRP association", y = "SNP-T2D association") +
  scale_colour_discrete("Steiger filtered")

ggsave("plots/cprplot_allconfoundingSNPs.pdf", width = 6.5, height = 5)

#scatter plot with rs429358 identified

dat_crp_t2d$filtered <- dat_crp_t2d$SNP %in% as.list('rs429358')

plotdat <- dat_crp_t2d %>%
  filter(mr_keep == TRUE) %>%
  select(beta.exposure, beta.outcome, filtered) %>%
  mutate(beta.outcome = ifelse(beta.exposure<0,-beta.outcome, beta.outcome)) %>%
  mutate(beta.exposure = abs(beta.exposure))

#plot with filtered snps highlighted
plot1 <- ggplot(plotdat, aes(x=beta.exposure)) +
  geom_point(aes(y=beta.outcome,  group=filtered, colour=filtered)) +
  scale_colour_brewer(palette = "Set1") +
  geom_abline(intercept = 0, slope = 0.1245) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "bottom") +
  labs(x = "SNP-CRP association", y = "SNP-T2D association") +
  scale_colour_discrete("Steiger filtered")

ggsave("plots/cprplot_rs429358.pdf", width = 6.5, height = 5)

