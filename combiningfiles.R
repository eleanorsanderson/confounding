library(dplyr)

###Misspecification set up

load("resultsB_1.Rda")

resB <- results_all

for(job_id in 2:20){
  load(file=sprintf("resultsB_%s.Rda", job_id))
  resB <- bind_rows(resB, results_all)  
}

nrow(resB)
res_check <- distinct(resB)
nrow(res_check)

save(resB, file = "results_misspec.Rda")


###Confounder set up 

load("results_1.Rda")

res <- results_all

for(job_id in 2:20){
  load(file=sprintf("results_%s.Rda", job_id))
  res <- bind_rows(res, results_all)  
}

nrow(res)
res_check <- distinct(res)
nrow(res_check)

res_a <- filter(res, res$model=="A")
save(res_a, file = "results_confA.Rda")

res_b <- filter(res, res$model=="B")
save(res_b, file = "results_confB.Rda")

res_c <- filter(res, res$model=="C")
save(res_c, file = "results_confC.Rda")

res_d <- filter(res, res$model=="D")
save(res_d, file = "results_confD.Rda")

res_e <- filter(res, res$model=="E")
save(res_e, file = "results_confE.Rda")

res_f <- filter(res, res$model=="F")
save(res_f, file = "results_confF.Rda")

res_g <- filter(res, res$model=="G")
save(res_g, file = "results_confG.Rda")

res_h <- filter(res, res$model=="H")
save(res_h, file = "results_confH.Rda")