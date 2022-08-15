library(dplyr)

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
save(res_a, file = "results_A.Rda")

res_b <- filter(res, res$model=="B")
save(res_b, file = "results_B.Rda")

res_c <- filter(res, res$model=="C")
save(res_c, file = "results_C.Rda")

res_d <- filter(res, res$model=="D")
save(res_d, file = "results_D.Rda")

res_e <- filter(res, res$model=="E")
save(res_e, file = "results_E.Rda")