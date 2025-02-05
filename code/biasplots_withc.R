# 07/08/2024
#Eleanor Sanderson

#code to make plot of the ols bias against iv bias for the simulation results
#calls the simulation results generated by "ols_2sls_bias.R" (separated to avoid rerunning the sims unnecessarily)
#This version includes an additional non-genetic confounder

##notes:
#bias_iv is the relative size of the pleitropic effect of Z to the effect of Z on X

library(dplyr)
library(ggplot2)

load("relativebias_sim_results_withc.rda")


#positive confounding

results <- results_all %>% 
  filter(pi_zx == 0, beta_cy == 1)

graph_res1 <- results %>% 
  group_by(beta_uy) %>% 
  summarise(across(c(wald, tsls_b, ols_b), mean))


results <- results_all %>% 
  filter(pi_zx == 0.25, beta_cy == 1)

graph_res2 <- results %>% 
  group_by(beta_uy) %>% 
  summarise(across(c(wald, tsls_b, ols_b), mean))

#combined graph

graph_res1$tsls_no <- graph_res1$tsls_b
comb_res <- graph_res1[, c("beta_uy", "tsls_no")]
comb_res <- full_join(comb_res, graph_res2, by="beta_uy")
comb_res$ols_bias <- comb_res$ols_b - 0.4
comb_res$tsls_bias <- comb_res$tsls_b - 0.4
comb_res$tsls_no_bias <- comb_res$tsls_no - 0.4


ggplot(data = comb_res, aes(x=beta_uy)) +
  geom_line(aes(y = ols_bias, colour="Linear regression")) +
  geom_line(aes(y = tsls_bias, colour="IV (direct effect)")) +
  geom_line(aes(y = tsls_no_bias, colour="IV (no direct effect)")) +
  scale_color_brewer(palette="Set1") +
  labs(x = "Effect of heritable confounder on outcome", y = "Bias", colour="Estimator") +
  theme_bw() +
  theme(legend.position = "bottom")



## negative confounding

results <- results_all %>% 
  filter(pi_zx == 0, beta_cy == -1)

graph_res1 <- results %>% 
  group_by(beta_uy) %>% 
  summarise(across(c(wald, tsls_b, ols_b), mean))


results <- results_all %>% 
  filter(pi_zx == 0.25, beta_cy == -1)

graph_res2 <- results %>% 
  group_by(beta_uy) %>% 
  summarise(across(c(wald, tsls_b, ols_b), mean))

#combined graph

graph_res1$tsls_no <- graph_res1$tsls_b
comb_res <- graph_res1[, c("beta_uy", "tsls_no")]
comb_res <- full_join(comb_res, graph_res2, by="beta_uy")
comb_res$ols_bias <- comb_res$ols_b - 0.4
comb_res$tsls_bias <- comb_res$tsls_b - 0.4
comb_res$tsls_no_bias <- comb_res$tsls_no - 0.4


ggplot(data = comb_res, aes(x=beta_uy)) +
  geom_line(aes(y = ols_bias, colour="Linear regression")) +
  geom_line(aes(y = tsls_bias, colour="IV (direct effect)")) +
  geom_line(aes(y = tsls_no_bias, colour="IV (no direct effect)")) +
  scale_color_brewer(palette="Set1") +
  labs(x = "Effect of heritable confounder on outcome", y = "Bias", colour="Estimator") +
  theme_bw() +
  theme(legend.position = "bottom")



