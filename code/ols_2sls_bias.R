# 29/11/2023
#Eleanor Sanderson

#Simulations to show that the bias of the IV estimator relative to the OLS estimator depends on the pleiotropic effect as a proportion 
#of the observed IV exposure association. 
#have a single IV so can illustrate the point and make plots

library(ivreg)
library(dplyr)
library(MASS)

results = data.frame()
results_all = data.frame()

#model parameters
reps = 1000
n = 10000

beta_xy = 0.4
beta_uy = 0.5
beta_ux = 0.5

#data set up
n_IVs <- 1
pi_zx <- 0.25

for(p in seq(0, 0.5, 0.025)){

    for(i in 1:reps){
    results[i, "pi_zu"] <- p
    results[i, "pi_zx"] <- pi_zx
    results[i, "beta_uy"] <- beta_uy
    results[i, "beta_ux"] <- beta_ux

    means <- c(0, 0, 0)                                   
    cov_matrix <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1),
                     ncol = 3)
    errors <- mvrnorm(n = 2*n,
                  mu = means, 
                  Sigma = cov_matrix)

    v_x <- errors[,1]
    v_y <- errors[,2]
    v_u <- errors[,3]

    Z  <- matrix(rbinom(2*n*n_IVs, 2, 0.4), 2*n, n_IVs)
    df <- data.frame(Z)
    df <- df %>% rename_at(vars(starts_with("X")), 
                       funs(str_replace(., "X", "Z")))

    df["U"] <- p*df$Z  + v_u
    df["X1"] <- pi_zx*df$Z + beta_ux*df$U + v_x
    df["Y"] <- beta_xy*df$X + beta_uy*df$U + v_y

    df_1 <- df[1:n,]
    df_2 <- df[(n+1):(2*n),]

    ols <- summary(lm(Y ~ X1, data = df_1))
    tsls <- summary(ivreg(Y ~ X1 | Z, data = df_1))

    zx_est <- lm(X1 ~ Z, data = df_1)$coefficients["Z"]
    zy_est <- lm(Y ~ Z, data = df_2)$coefficients["Z"]

    results[i,"wald"] <- zy_est/zx_est
    results[i,"ols_b"] <- ols$coefficients["X1","Estimate"]
    results[i,"ols_se"] <- ols$coefficients["X1","Std. Error"]
    results[i,"tsls_b"] <- tsls$coefficients["X1","Estimate"]
    results[i,"tsls_se"] <- tsls$coefficients["X1","Std. Error"]
    results[i,"tsls_F"] <- tsls$diagnostics["Weak instruments","statistic"]
    
    }
  results_all <- rbind(results_all, results)
}


results_all[,"bias_ols"] <- results_all[,"beta_uy"]*results_all[,"beta_ux"]
results_all[,"bias_iv"] <- (results_all[,"beta_uy"]*results_all[,"pi_zu"])/(results_all[,"pi_zx"] + results_all[,"beta_ux"]*results_all[,"pi_zu"])


save(results_all, file = "relativebias_sim_results.rda")

