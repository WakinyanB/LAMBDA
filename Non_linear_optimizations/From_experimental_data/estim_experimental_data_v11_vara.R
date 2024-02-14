# Cleaning objects from the workplace
rm(list=ls())

# setwd() # Set working directory

# Packages (may first require installations: install.packages())
library(tidyverse) # for data manipulation
library(plyr) # for data manipulation
library(deSolve) # to solve a system of ordinary differential equations
library(boot) # for the inverse logit function 'inv.logit'

# Execute functions from external R script
source("Phage_Lambda_functions_v11.R")

# Import, processing and graphical visualization of experimental data

## Import

data_FACS <- read.csv("Data/data_cleaned_FACS_v3.csv", header = TRUE)
data_qPCR <- read.csv("Data/data_cleaned_qPCR_v2.csv", header = TRUE)

# Inference

## Estimation of reactivation rates

estim_alpha_exp <- estim_alpha(merge.data.frame(data_FACS[,c(1:3,9:12)], data_qPCR[,c(1:3,6,8)], all=TRUE),
                               by_chemostat = FALSE)

## Non linear optimizations

n_starts <- 2000

# Starting values of estimated parameters
parms0 <- draw_initial_conditions(n=n_starts, seed=123)
parms0 <- parms0[,-which(colnames(parms0)=="a")]

# Fixed burst size (original value)
B <- 80

for(a_var in c(1E-9, 1E-8, 1E-7, 1E-6, 1E-4, 1E-2, 1)){
  
  estim_df <- Estim_nmk_parallel(parms0=parms0, constants=c("alpha_w"=estim_alpha_exp[["alpha_w"]],
                                                            "alpha_m"=estim_alpha_exp[["alpha_m"]],
                                                            "B"=B, "a"=a_var),
                                 data_FACS=data_FACS, data_qPCR=data_qPCR,
                                 maxfeval=4000, Ncores=30)
  
  write.csv(estim_df, paste0("Outputs_csv/Estim_nmk_parallel_experimental_data_nstarts", n_starts, "_a", a_var, "_v11.csv"), row.names = FALSE)
}