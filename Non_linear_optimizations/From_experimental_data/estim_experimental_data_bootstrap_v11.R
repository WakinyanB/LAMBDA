rm(list=ls())

# setwd() # Set working directory

source("Phage_Lambda_functions_v11.R")

B <- 80

data.bs <- read.csv("Outputs_csv/Bootstrapped_data.csv", header=TRUE)
best_parms <- read.csv("Outputs_csv/Best_parms_non_linear_optim_experimental_data.csv", header=TRUE)

t0 <- Sys.time()
tab <- Estim_nmk_parallel_bootstrap(best_parms, bootstrapped.data = data.bs, constants = c("B"=B), maxfeval = 4000, Ncores = 38)
Sys.time()-t0

write.csv(tab, paste0("Outputs_csv/Estim_nmk_parallel_experimental_data_bootstrap_B", B, "_v11.csv"), row.names = FALSE)