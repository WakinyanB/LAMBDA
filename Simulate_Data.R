rm(list=ls()) # Cleaning objects from the global environment

# Packages (may first require installations: install.packages())
library(tidyverse) # for data manipulation, graphic visualizations, ...
library(plyr) # for data manipulation
library(deSolve) # to solve a system of ordinary differential equations
library(boot) # for the logit and inverse logit function ('logit' and 'inv.logit')

setwd("C:/Users/wbenhamou/Desktop/LAMBDA/Data_and_codes/")

# Execute functions from external R script
source("Phage_Lambda_functions_v11.R")

Sys.setlocale("LC_TIME", "English")

# Simulate data for 2 treatments (epidemic vs endemic) x 4 chemostats each,
# with sampling every 6 min (0.1 h) and without measurement errors.

# Time points
times <- seq(0, 60, 0.1)
n_times <- length(times)

# Parameters

## Values for initial conditions
prev0_epidemic <- 0.01
prev0_endemic <- 0.99
freq0 <- 0.5

# Other parameters
alpha_w <- 7e-03
alpha_m <- 2e-02
phi_w <- 2e-01
phi_m <- 2e-02
a <- 3e-09
b <- 0.1
B <- 80
r <- 1.4
tau <- 1.5 # lysis delay = 40 min
K <- 1e+09
m <- 0.8

parms <- c("prev0_epidemic" = prev0_epidemic, "prev0_endemic" = prev0_endemic, "freq0" = freq0,
           "alpha_w" = alpha_w, "alpha_m" = alpha_m, "phi_w" = phi_w, "phi_m" = phi_m,
           "a" = a, "b" = b, "B" = B, "r" = r, "tau" = tau)

# parms %>% data.frame %>% t %>% write.csv("Data/Simulated_data/parms_simul.csv", row.names = FALSE)

## Numerical simulations

simul <- rbind(lambda_simul(prev0=parms[["prev0_epidemic"]], freq0=parms[["freq0"]], parms=parms, times=times),
               lambda_simul(prev0=parms[["prev0_endemic"]], freq0=parms[["freq0"]], parms=parms, times=times))
simul$prev0 <- c(1, 100) %>% rep(each = n_times)

simul_data_FACS <- data.frame("chemostat"=c("chem_epidemic_", "chem_endemic_") %>%
                                paste0(rep(1:4, each=2)) %>% rep(each=n_times),
                              "time"=simul$time, "prev0"=simul$prev0,
                              "logit_g"=simul$logit_g, "logit_prev"=simul$logit_prev)

simul_data_qPCR <- data.frame("chemostat"=c("chem_epidemic_", "chem_endemic_") %>%
                                paste0(rep(1:4, each=2)) %>% rep(each=n_times),
                              "time"=simul$time, "prev0"=simul$prev0, "logit_q"=simul$logit_q)

## Samplings every 0.1 h with almost no measurement errors (sd = 0.01)

simul_data_FACS.1 <- simul_data_FACS
simul_data_qPCR.1 <- simul_data_qPCR

# Measurement errors (Gaussian noise)

# Standard deviations of measurement errors
sigma_prev <- 0.01
sigma_g <- 0.01
sigma_q <- 0.01

# Parameters in common for .1 and .3
parms.13 <- c(parms, "sigma_prev"=sigma_prev, "sigma_g"=sigma_g, "sigma_q"=sigma_q)

set.seed(123)

simul_data_FACS.1$logit_prev <- simul_data_FACS.1$logit_prev + rnorm(n=nrow(simul_data_FACS.1), mean=0, sd=sigma_prev)
simul_data_FACS.1$logit_g <- simul_data_FACS.1$logit_g + rnorm(n=nrow(simul_data_FACS.1), mean=0, sd=sigma_g)
simul_data_qPCR.1$logit_q <- simul_data_qPCR.1$logit_q + rnorm(n=nrow(simul_data_qPCR.1), mean=0, sd=sigma_q)

# write.csv(simul_data_FACS.1, "Data/Simulated_data/simul_data_FACS_1.csv", row.names = FALSE)
# write.csv(simul_data_qPCR.1, "Data/Simulated_data/simul_data_qPCR_1.csv", row.names = FALSE)

## Samplings every 0.1 h with greater measurement errors (sd = 0.5)

simul_data_FACS.2 <- simul_data_FACS
simul_data_qPCR.2 <- simul_data_qPCR

# Measurement errors (Gaussian noise)

# Standard deviations of measurement errors
sigma_prev <- 0.5
sigma_g <- 0.5
sigma_q <- 0.5

# Parameters in common between .2 and .4
parms.24 <- c(parms, "sigma_prev"=sigma_prev, "sigma_g"=sigma_g, "sigma_q"=sigma_q)

set.seed(123)

simul_data_FACS.2$logit_prev <- simul_data_FACS.2$logit_prev + rnorm(n=nrow(simul_data_FACS.2), mean=0, sd=sigma_prev)
simul_data_FACS.2$logit_g <- simul_data_FACS.2$logit_g + rnorm(n=nrow(simul_data_FACS.2), mean=0, sd=sigma_g)
simul_data_qPCR.2$logit_q <- simul_data_qPCR.2$logit_q + rnorm(n=nrow(simul_data_qPCR.2), mean=0, sd=sigma_q)

# write.csv(simul_data_FACS.2, "Data/Simulated_data/simul_data_FACS_2.csv", row.names = FALSE)
# write.csv(simul_data_qPCR.2, "Data/Simulated_data/simul_data_qPCR_2.csv", row.names = FALSE)

## Samplings every 1 h with almost no measurement errors (sd = 0.01)

times_1h <- 1:60 

simul_data_FACS.3 <- simul_data_FACS[simul_data_FACS$time %in% times_1h,]
simul_data_qPCR.3 <- simul_data_qPCR[simul_data_FACS$time %in% times_1h,]

# Measurement errors (Gaussian noise)

# Standard deviations of measurement errors
sigma_prev <- 0.01
sigma_g <- 0.01
sigma_q <- 0.01

set.seed(123)

simul_data_FACS.3$logit_prev <- simul_data_FACS.3$logit_prev + rnorm(n=nrow(simul_data_FACS.3), mean=0, sd=sigma_prev)
simul_data_FACS.3$logit_g <- simul_data_FACS.3$logit_g + rnorm(n=nrow(simul_data_FACS.3), mean=0, sd=sigma_g)
simul_data_qPCR.3$logit_q <- simul_data_qPCR.3$logit_q + rnorm(n=nrow(simul_data_qPCR.3), mean=0, sd=sigma_q)

# write.csv(simul_data_FACS.3, "Data/Simulated_data/simul_data_FACS_3.csv", row.names = FALSE)
# write.csv(simul_data_qPCR.3, "Data/Simulated_data/simul_data_qPCR_3.csv", row.names = FALSE)

## Samplings every 1 h with greater measurement errors (sd = 0.5)

simul_data_FACS.4 <- simul_data_FACS[simul_data_FACS$time %in% times_1h,]
simul_data_qPCR.4 <- simul_data_qPCR[simul_data_FACS$time %in% times_1h,]

# Measurement errors (Gaussian noise)

# Standard deviations of measurement errors
sigma_prev <- 0.5
sigma_g <- 0.5
sigma_q <- 0.5

set.seed(123)

simul_data_FACS.4$logit_prev <- simul_data_FACS.4$logit_prev + rnorm(n=nrow(simul_data_FACS.4), mean=0, sd=sigma_prev)
simul_data_FACS.4$logit_g <- simul_data_FACS.4$logit_g + rnorm(n=nrow(simul_data_FACS.4), mean=0, sd=sigma_g)
simul_data_qPCR.4$logit_q <- simul_data_qPCR.4$logit_q + rnorm(n=nrow(simul_data_qPCR.4), mean=0, sd=sigma_q)

# write.csv(simul_data_FACS.4, "Data/Simulated_data/simul_data_FACS_4.csv", row.names = FALSE)
# write.csv(simul_data_qPCR.4, "Data/Simulated_data/simul_data_qPCR_4.csv", row.names = FALSE)

# save(simul, parms,
#      simul_data_FACS.1, simul_data_qPCR.1,
#      simul_data_FACS.2, simul_data_qPCR.2,
#      simul_data_FACS.3, simul_data_qPCR.3,
#      simul_data_FACS.4, simul_data_qPCR.4,
#      file = "Data/Simulated_data/simul_data.RData")
