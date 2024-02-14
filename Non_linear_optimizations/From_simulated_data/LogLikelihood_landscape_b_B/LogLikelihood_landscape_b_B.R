rm(list=ls())

setwd("C:/Users/wbenhamou/Desktop/LAMBDA/Data_and_codes/")

# Execute functions from external R script
source("Phage_Lambda_functions_v11.R")

## SIMULATED DATA --------------------------------------------------------------

simul_data_FACS <- read.csv("Data/Simulated_data/simul_data_FACS_1.csv", header=TRUE)
simul_data_qPCR <- read.csv("Data/Simulated_data/simul_data_qPCR_1.csv", header=TRUE)

parms <- read.csv("Data/Simulated_data/parms_simul.csv", header=TRUE)
parms$sigma_q <- parms$sigma_g <- parms$sigma_prev <- 0.01

## PARAMETERS OF INTEREST ------------------------------------------------------

n <- 201
constants_grid <- expand.grid("b" = seq(0, 0.2, length.out=n),
                              "B" = seq(0, 100, length.out=n)) %>% as.matrix

## INFERENCE -------------------------------------------------------------------

transf_parms0 <- parm_transformation(parms[!names(parms) %in% c("alpha_w", "alpha_m", "b", "B")])
parm_names <- colnames(transf_parms0)
n_parms <- ncol(transf_parms0)

maxfeval=8000; tol=1e-6; restarts.max=3; solver="lsoda"

# (3) Parallelized optimizations
  
Ncores <- 38

t <- format(Sys.time(), "%Y-%m-%d %Hh%M")
folder_name <- paste0("Lambda_v11_Landscape_estim_nmk_parallel_", stringi::stri_replace_all_fixed(t, c("-"," "),
                                                                                                  "_", vectorize_all=FALSE))
dir.create(folder_name)

cl <- parallel::makeCluster(Ncores)
doParallel::registerDoParallel(cl)

estim_df <- foreach::foreach(i = 1:nrow(constants_grid), .combine = 'rbind', .inorder = FALSE,
                             .packages = c("tidyverse", "deSolve", "boot"), .export = ls(globalenv())) %dopar% {
                               
                               fit <- dfoptim::nmk(par = transf_parms0[1,], fn = negative_log_likelihood,
                                                   control = list(maximize = FALSE, maxfeval=maxfeval,
                                                                  tol=tol, restarts.max=restarts.max),
                                                   parm_names=parm_names, constants=c(constants_grid[i,],"alpha_w"=parms$alpha_w,"alpha_m"=parms$alpha_m),
                                                   data_FACS=simul_data_FACS, data_qPCR=simul_data_qPCR,
                                                   solver=solver, verbose=FALSE, return_all=FALSE)
                               
                               res <- c("Start"=i, constants_grid[i,], "Value"=fit$value, "Convergence"=fit$convergence,
                                        parm_transformation(fit$par %>% setNames(parm_names), reverse=TRUE))
                               saveRDS(res, file = paste0(folder_name, "/result_task_", i, ".rds"))
                               return(res)
                             }
parallel::stopCluster(cl)

# (4) Processing of the outcome

estim_df <- as.data.frame(estim_df)
estim_df <- estim_df[order(estim_df$Start),]

estim_df$Convergence[estim_df$Convergence == 0] <- "Successful completion"
estim_df$Convergence[estim_df$Convergence == 1] <- "Maximum number of evaluations exceeded"
estim_df$Convergence[estim_df$Convergence == 2] <- "Stagnation in Nelder-Mead"
# estim_df$Convergence[estim_df$Convergence == 10] <- "Degeneracy of the Nelder-Mead simplex"

write.csv(estim_df, paste0("Outputs_csv/Lambda_landscape_b_B_v11_n", nrow(constants_grid), "_maxfeval", maxfeval, ".csv"),
          row.names=FALSE)