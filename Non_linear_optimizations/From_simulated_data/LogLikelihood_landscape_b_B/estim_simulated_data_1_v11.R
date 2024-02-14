rm(list=ls())

# setwd() # set working directory

source("Phage_Lambda_functions_v11.R") # Run all functions

## SIMULATED DATA --------------------------------------------------------------

parms <- read.csv("Data/Simulated_data/parms_simul.csv", header= TRUE)
simul_data_FACS <- read.csv("Data/Simulated_data/simul_data_FACS_1.csv", header = TRUE)
simul_data_qPCR <- read.csv("Data/Simulated_data/simul_data_qPCR_1.csv", header = TRUE)

## INFERENCE -------------------------------------------------------------------

# Estimation of alpha_w and alpha_m

estim_alpha_df <- estim_alpha(data=merge.data.frame(simul_data_FACS, simul_data_qPCR, all=TRUE), by_chemostat=FALSE)

# Starting values of estimated parameters

n_starts <- 5*2000
set.seed(123)

parms0 <- data.frame("prev0_epidemic" = n_starts %>% runif(bounds$prev0_epidemic[1], bounds$prev0_epidemic[2]),
                     "prev0_endemic" = n_starts %>% runif(bounds$prev0_endemic[1], bounds$prev0_endemic[2]),
                     "freq0" = n_starts %>% runif(bounds$freq0[1], bounds$freq0[2]),
                     "phi_w" = n_starts %>% runif(bounds$phi_w[1], bounds$phi_w[2]),
                     "phi_m" = n_starts %>% runif(bounds$phi_m[1], bounds$phi_m[2]),
                     "a" = n_starts %>% runif(bounds$a[1], bounds$a[2]),
                     "b" = n_starts %>% runif(bounds$b[1], bounds$b[2]),
                     "r" = n_starts %>% runif(bounds$r[1], bounds$r[2]),
                     "tau" = n_starts %>% runif(bounds$tau[1], bounds$tau[2]),
                     "sigma_prev" = n_starts %>% runif(bounds$sigma_prev[1], bounds$sigma_prev[2]),
                     "sigma_g" = n_starts %>% runif(bounds$sigma_g[1], bounds$sigma_g[2]),
                     "sigma_q" = n_starts %>% runif(bounds$sigma_q[1], bounds$sigma_q[2])) %>%
  as.data.frame

t0 <- Sys.time()

estim_df <- Estim_nmk_parallel(parms0 = parms0, constants = c("alpha_w"=estim_alpha_df$alpha_w,
                                                              "alpha_m"=estim_alpha_df$alpha_m,
                                                              "B"=parms$B),
                               data_FACS = simul_data_FACS, data_qPCR = simul_data_qPCR,
                               maxfeval = 10000, Ncores = 23)
Sys.time()-t0

write.csv(estim_df, "Outputs_csv/Estim_nmk_parallel_simulated_data_sd_noise_001_timestep_01_n5x2000_seed123_v11.csv", row.names = FALSE)