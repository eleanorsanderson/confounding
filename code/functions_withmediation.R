
datagenA <- function(nsnps,snpsc,ss,beta1,beta2,pi){
  
  af = 0.4 
  n=2*ss
  
  G  <- matrix(rbinom(n*nsnps, 2, af), n, nsnps)
  G2 <- matrix(rbinom(n*snpsc, 2, af), n, snpsc)
  
  means <- c(0, 0)                                   
  cov_matrix <- matrix(c(1.5, 0, 0, 1),
                                     ncol = 2)
  
  # create bivariate normal distribution
  errors <- mvrnorm(n = n,
                    mu = means, 
                    Sigma = cov_matrix)
  
  v_x1 <- errors[,1]
  v_y <- errors[,2]
  v_x2 <- rnorm(n,0,1.5)
  v_m <- rnorm(n,0,1.5)
  
  effs_x1 <- abs(rnorm(nsnps,0,0.04))
  effs_x2 <- abs(rnorm(snpsc,0,0.04))
  
  df <- data.frame(cbind(G, G2))
  df <- df %>% rename_at(vars(starts_with("X")), 
                         funs(str_replace(., "X", "G")))
  
  df[,"X2"] <- pi*G2[,]%*%effs_x2 + v_x2
  df[,"M"] <- G[,]%*%effs_x1 + v_m
  df[,"X1"] <- df[,"M"] + df[,"X2"] + v_x1
  df[,"Y"] <- beta1*df[,"X1"] + beta2*df[,"X2"] + v_y  
  
  data <- cbind.data.frame(df)
  return(data)
}



GWASres <- function(dat){
  MR_dat = data.frame()

  dat.1 <- dat[1:nobs,]
  dat.2 <- dat[(nobs+1):(2*nobs),]

  est.snps <- snps + snpsc

    for(i in 1:est.snps){
      a <- summary(lm(dat.1$X1~dat.1[,i]))
      MR_dat[i,"X1_b"] <- a$coefficient[2,1]
      MR_dat[i,"X1_se"] <- a$coefficient[2,2]
      MR_dat[i,"X1_p"] <- a$coefficient[2,4]
      MR_dat[i,"X1_r2"] <- a$r.squared
      b <- summary(lm(dat.1$X2~dat.1[,i]))
      MR_dat[i,"X2_b"] <- b$coefficient[2,1]
      MR_dat[i,"X2_se"] <- b$coefficient[2,2]
      MR_dat[i,"X2_p"] <- b$coefficient[2,4]
      MR_dat[i,"X2_r2"] <- b$r.squared
      c<-summary(lm(dat.2$Y~dat.2[,i]))
      MR_dat[i,"Y_b"] <- c$coefficient[2,1]
      MR_dat[i,"Y_se"] <-c$coefficient[2,2]
      MR_dat[i,"Y_p"] <- c$coefficient[2,4]
  
    }

  return(MR_dat)
}

##MR functions

default_parameters <- function()
{
  list(
    test_dist = "z",
    nboot = 1000,
    Cov = 0,
    penk = 20,
    phi = 1,
    alpha = 0.05,
    Qthresh = 0.05,
    over.dispersion = TRUE,
    loss.function = "huber",
    shrinkage = FALSE
  )
}

#IVW
mr_ivw <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
{
  if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2)
    return(list(b=NA, se=NA, pval=NA, nsnp=NA))
  
  ivw.res <- summary(lm(b_out ~ -1 + b_exp, weights = 1/se_out^2))
  b <- ivw.res$coef["b_exp","Estimate"]
  se <- ivw.res$coef["b_exp","Std. Error"]/min(1,ivw.res$sigma) #sigma is the residual standard error
  pval <- 2 * pnorm(abs(b/se), lower.tail=FALSE)
  Q_df <- length(b_exp) - 1
  Q <- ivw.res$sigma^2 * Q_df
  Q_pval <- pchisq(Q, Q_df, lower.tail=FALSE)
  # from formula phi =  Q/DF rearranged to to Q = phi*DF, where phi is sigma^2
  # Q.ivw<-sum((1/(se_out/b_exp)^2)*(b_out/b_exp-ivw.reg.beta)^2)
  return(list(b = b, se = se, pval = pval, nsnp=length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}

##MR Egger
mr_egger_regression <- function(b_exp, b_out, se_exp, se_out, parameters)
{
  stopifnot(length(b_exp) == length(b_out))
  stopifnot(length(se_exp) == length(se_out))
  stopifnot(length(b_exp) == length(se_out))
  
  # print(b_exp)
  
  nulllist <- list(
    b = NA,
    se = NA,
    pval = NA,
    nsnp = NA,
    b_i = NA,
    se_i = NA,
    pval_i = NA,
    Q = NA,
    Q_df = NA,
    Q_pval = NA,
    mod = NA,
    smod = NA,
    dat = NA
  )
  if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3)
  {
    return(nulllist)
  }
  
  sign0 <- function(x)
  {
    x[x==0] <- 1
    return(sign(x))
  }
  
  to_flip <- sign0(b_exp) == -1
  b_out = b_out*sign0(b_exp)
  b_exp = abs(b_exp)
  dat <- data.frame(b_out=b_out, b_exp=b_exp, se_exp=se_exp, se_out=se_out, flipped=to_flip)
  mod <- lm(b_out ~ b_exp, weights=1/se_out^2)
  smod <- summary(mod)
  if(nrow(coefficients(smod)) > 1)
  {
    b <- coefficients(smod)[2,1]
    se <- coefficients(smod)[2,2] / min(1,smod$sigma)
    pval <- 2 * pt(abs(b / se), length(b_exp) - 2, lower.tail = FALSE)
    b_i <- coefficients(smod)[1,1]
    se_i <- coefficients(smod)[1,2] / min(1,smod$sigma)
    pval_i <- 2 * pt(abs(b_i / se_i), length(b_exp) - 2, lower.tail = FALSE)
    
    Q <- smod$sigma^2 * (length(b_exp) - 2)
    Q_df <- length(b_exp) - 2
    Q_pval <- pchisq(Q, Q_df, lower.tail=FALSE)
  } else {
    warning("Collinearities in MR Egger, try LD pruning the exposure variables.")
    return(nulllist)
  }
  return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), b_i = b_i, se_i = se_i, pval_i = pval_i, Q = Q, Q_df = Q_df, Q_pval = Q_pval, mod = smod, dat = dat))
}



##weighted median
mr_weighted_median <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
{
  if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3)
    return(list(b=NA, se=NA, pval=NA, nsnp=NA))
  
  b_iv <- b_out / b_exp
  VBj <- ((se_out)^2)/(b_exp)^2 + (b_out^2)*((se_exp^2))/(b_exp)^4
  b <- weighted_median(b_iv, 1 / VBj)
  se <- weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, 1 / VBj, parameters$nboot)
  pval <- 2 * pnorm(abs(b/se), lower.tail=FALSE)
  return(list(b=b, se=se, pval=pval, Q=NA, Q_df=NA, Q_pval=NA, nsnp=length(b_exp)))
}

weighted_median_bootstrap <- function(b_exp, b_out, se_exp, se_out, weights, nboot)
{
  med <- rep(0, nboot)
  for(i in 1:nboot){
    b_exp.boot = rnorm(length(b_exp), mean=b_exp, sd=se_exp)
    b_out.boot = rnorm(length(b_out), mean=b_out, sd=se_out)
    betaIV.boot = b_out.boot/b_exp.boot
    med[i] = weighted_median(betaIV.boot, weights)
  }
  return(sd(med))
}

weighted_median <- function(b_iv, weights)
{
  betaIV.order <- b_iv[order(b_iv)]
  weights.order <- weights[order(b_iv)]
  weights.sum <- cumsum(weights.order)-0.5*weights.order
  weights.sum <- weights.sum/sum(weights.order)
  below <- max(which(weights.sum<0.5))
  b = betaIV.order[below] + (betaIV.order[below+1]-betaIV.order[below])*
    (0.5-weights.sum[below])/(weights.sum[below+1]-weights.sum[below])
  return(b)
}


mr_mode <- function(dat, parameters=default_parameters(), mode_method="all")
{
  if("mr_keep" %in% names(dat)) dat <- subset(dat, mr_keep)
  
  if(nrow(dat) < 3) 
  {
    warning("Need at least 3 SNPs")
    return(NULL)
  }
  
  b_exp <- dat$beta.exposure
  b_out <- dat$beta.outcome
  se_exp <- dat$se.exposure
  se_out <- dat$se.outcome
  
  #--------------------------------------#
  #Function to compute the point estimate#
  #--------------------------------------#
  #BetaIV.in: ratio estimates
  #seBetaIV.in: standard errors of ratio estimates
  #' @importFrom stats mad sd
  beta <- function(BetaIV.in, seBetaIV.in, phi)
  {
    #Bandwidth rule - modified Silverman's rule proposed by Bickel (2002)
    s <- 0.9*(min(sd(BetaIV.in), mad(BetaIV.in)))/length(BetaIV.in)^(1/5)
    
    #Standardised weights
    weights <- seBetaIV.in^-2/sum(seBetaIV.in^-2)
    
    beta <- NULL
    
    for(cur_phi in phi)
    {
      #Define the actual bandwidth
      h <- max(0.00000001, s*cur_phi)
      #Compute the smoothed empirical density function
      densityIV <- density(BetaIV.in, weights=weights, bw=h)
      #Extract the point with the highest density as the point estimate 
      beta[length(beta)+1] <- densityIV$x[densityIV$y==max(densityIV$y)]
    }
    return(beta)
  }
  
  #------------------------------------------#
  #Function to estimate SEs through bootstrap#
  #------------------------------------------#
  #BetaIV.in: ratio estimates
  #seBetaIV.in: standard errors of ratio estimates
  #beta_Mode.in: point causal effect estimates
  #' @importFrom stats density pchisq rnorm
  boot <- function(BetaIV.in, seBetaIV.in, beta_Mode.in, nboot)
  {
    #Set up a matrix to store the results from each bootstrap iteration
    beta.boot <- matrix(nrow=nboot, ncol=length(beta_Mode.in))
    
    for(i in 1:nboot) 
    {
      #Re-sample each ratio estimate using SEs derived not assuming NOME
      BetaIV.boot      <- rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in[,1])
      #Re-sample each ratio estimate using SEs derived under NOME
      BetaIV.boot_NOME <- rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in[,2])
      
      #Simple mode, not assuming NOME
      beta.boot[i,1:length(phi)] <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=rep(1, length(BetaIV)), phi=phi)
      #Weighted mode, not assuming NOME
      beta.boot[i,(length(phi)+1):(2*length(phi))] <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=seBetaIV.in[,1], phi=phi)
      #Penalised mode, not assuming NOME
      weights <- 1/seBetaIV.in[,1]^2
      penalty <- pchisq(weights * (BetaIV.boot-beta.boot[i,(length(phi)+1):(2*length(phi))])^2, df=1, lower.tail=FALSE)
      pen.weights <- weights*pmin(1, penalty*parameters$penk)
      
      beta.boot[i,(2*length(phi)+1):(3*length(phi))] <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=sqrt(1/pen.weights), phi=phi)
      
      #Simple mode, assuming NOME
      beta.boot[i,(3*length(phi)+1):(4*length(phi))] <- beta(BetaIV.in=BetaIV.boot_NOME, seBetaIV.in=rep(1, length(BetaIV)), phi=phi)
      #Weighted mode, assuming NOME
      beta.boot[i,(4*length(phi)+1):(5*length(phi))] <- beta(BetaIV.in=BetaIV.boot_NOME, seBetaIV.in=seBetaIV.in[,2], phi=phi)
    }
    return(beta.boot)
  }
  
  # Parameters
  phi <- parameters$phi
  nboot <- parameters$nboot
  alpha <- parameters$alpha
  
  #Ratio estimates
  BetaIV   <- b_out/b_exp
  
  #SEs of ratio estimates
  seBetaIV <- cbind(sqrt((se_out^2)/(b_exp^2) + ((b_out^2)*(se_exp^2))/(b_exp^4)), #SEs NOT assuming NOME
                    se_out/abs(b_exp)) #SEs ASSUMING NOME
  
  #Point causal effect estimate using the simple mode
  beta_SimpleMode <- beta(BetaIV.in=BetaIV, seBetaIV.in=rep(1, length(BetaIV)), phi=phi)
  
  #Point causal effect estimate using the weighted mode (not asusming NOME)
  beta_WeightedMode <- beta(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV[,1], phi=phi)
  weights <- 1/seBetaIV[,1]^2
  penalty <- pchisq(weights * (BetaIV-beta_WeightedMode)^2, df=1, lower.tail=FALSE)
  pen.weights <- weights*pmin(1, penalty*parameters$penk) # penalized 
  
  beta_PenalisedMode <- beta(BetaIV.in=BetaIV, seBetaIV.in=sqrt(1/pen.weights), phi=phi)
  
  #Point causal effect estimate using the weighted mode (asusming NOME)
  beta_WeightedMode_NOME <- beta(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV[,2], phi=phi)
  
  #Combine all point effect estimates in a single vector
  beta_Mode <- rep(c(beta_SimpleMode, beta_WeightedMode, beta_PenalisedMode, beta_SimpleMode, beta_WeightedMode_NOME))
  
  #Compute SEs, confidence intervals and P-value
  beta_Mode.boot <- boot(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV, beta_Mode.in=beta_Mode, nboot=nboot)
  se_Mode <- apply(beta_Mode.boot, 2, mad)
  
  CIlow_Mode <- beta_Mode-qnorm(1-alpha/2)*se_Mode
  CIupp_Mode <- beta_Mode+qnorm(1-alpha/2)*se_Mode
  
  P_Mode <- pt(abs(beta_Mode/se_Mode), df=length(b_exp)-1, lower.tail=FALSE)*2
  
  #Vector to indicate the method referring to each row
  Method <- rep(c('Simple mode', 'Weighted mode', 'Penalised mode', 'Simple mode (NOME)', 'Weighted mode (NOME)'), each=length(phi))
  
  #Return a data frame containing the results
  id.exposure <- ifelse("id.exposure" %in% names(dat), dat$id.exposure[1], "")
  id.outcome <- ifelse("id.outcome" %in% names(dat), dat$id.outcome[1], "")
  Results <- data.frame(
    id.exposure = id.exposure, 
    id.outcome = id.outcome, 
    method = Method, 
    nsnp = length(b_exp), 
    b = beta_Mode, 
    se = se_Mode, 
    ci_low = CIlow_Mode, 
    ci_upp = CIupp_Mode, 
    pval = P_Mode, 
    stringsAsFactors=FALSE
  )
  
  if(mode_method == "all")
  {
    return(Results)
  } else {
    stopifnot(all(mode_method %in% Results$method))
    i <- which(Results$method == mode_method)
    return(list(b = Results$b[i], se = Results$se[i], pval=Results$pval[i], nsnp=length(b_exp)))
  }
  
  return(Results)
}

mr_weighted_mode <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters()) 
{
  index <- !is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)
  if(sum(index) < 3)
    return(list(b=NA, se=NA, pval=NA, nsnp=NA))
  
  b_exp <- b_exp[index]
  b_out <- b_out[index]
  se_exp <- se_exp[index]
  se_out <- se_out[index]
  
  return(mr_mode(data.frame(beta.exposure=b_exp, beta.outcome=b_out, se.exposure=se_exp, se.outcome=se_out), parameters=parameters, mode_method="Weighted mode"))
}




##MVMR estimations
MVMR_conf <- function(MR_dat){
  
  res = data.frame()
  MR_dat$minp <-  apply(MR_dat[,c("X1_p", "X2_p")],1,min)
  MVMR_dat <- MR_dat[MR_dat$minp < pcut,]

  mvmr_out <- summary(lm(Y_b~-1+X1_b+X2_b,weights=1/(Y_se^2), data = MVMR_dat))

  res[1,"mvmr.X1_effect"] <- mvmr_out$coefficients["X1_b","Estimate"]
  res["mvmr.X1_se"] <- mvmr_out$coefficients["X1_b","Std. Error"]
  res["mvmr.X2_effect"] <-mvmr_out$coefficients["X2_b","Estimate"]
  res["mvmr.X2_se"] <- mvmr_out$coefficients["X2_b","Std. Error"]
  res["mvmr.nsnps"] <- length(MVMR_dat$X1_b)

  rho = cor(dat$X1,dat$X2)
  sig12 = as.vector(rho)*MVMR_dat$X1_se*MVMR_dat$X2_se

  delta1 <- lm(X1_b~ -1 + X2_b,data = MVMR_dat)$coefficients["X2_b"]
  delta2 <- lm(X2_b~ -1 + X1_b,data = MVMR_dat)$coefficients["X1_b"]
 
  Qind_1 <- ((MVMR_dat$X1_b - delta1*MVMR_dat$X2_b)^2)/(MVMR_dat$X1_se^2 + (delta1^2)*MVMR_dat$X2_se^2 - 2*delta1*sig12)
  Qind_2 <- ((MVMR_dat$X2_b - delta2*MVMR_dat$X1_b)^2)/(MVMR_dat$X2_se^2 + (delta2^2)*MVMR_dat$X1_se^2 - 2*delta2*sig12)

  res["Cond_F.X1"] <- mean(Qind_1)
  res["Cond_F.X2"] <- mean(Qind_2)

 
# ##MVMR instrument strength testing
# 
  rho = cor(dat$X1,dat$X2)
  sig12 = as.vector(rho)*MR_dat.X$X1_se*MR_dat.X$X2_se

  delta1 <- lm(X1_b~ -1 + X2_b,data = MR_dat.X)$coefficients["X2_b"]
  delta2 <- lm(X2_b~ -1 + X1_b,data = MR_dat.X)$coefficients["X1_b"]
 
  Qind_1 <- ((MR_dat.X$X1_b - delta1*MR_dat.X$X2_b)^2)/(MR_dat.X$X1_se^2 + (delta1^2)*MR_dat.X$X2_se^2 - 2*delta1*sig12)
  Qind_2 <- ((MR_dat.X$X2_b - delta2*MR_dat.X$X1_b)^2)/(MR_dat.X$X2_se^2 + (delta2^2)*MR_dat.X$X1_se^2 - 2*delta2*sig12)
 
  res["CF.X1_uniS"] <- mean(Qind_1)
  res["CF.X2_uniS"] <- mean(Qind_2)
  
  return(res)

}

