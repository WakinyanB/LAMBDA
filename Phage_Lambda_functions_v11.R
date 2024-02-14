#########################################################################
#######                   Bacteriophage lambda                    #######
#######     Numerical simulations and statistical inference       #######
#########################################################################
#######                  R Script with functions                  #######
#########################################################################

# Packages (may first require installations: install.packages())

library(tidyverse) # for data manipulation, graphic visualizations, etc...
library(plyr) # for the function 'ddply'
#library(scales) # to manipulate the internal scaling infrastructure used by ggplot2
#library(RColorBrewer) # for colors and palettes
#library(gridExtra) # to show multiple plots
#library(cowplot) # to show multiple plots
#library(ggpubr) # to show multiple plots
library(deSolve) # to solve a system of ordinary differential equations (ODE integration, function 'ode')
library(boot) # for the inverse logit function 'inv.logit'
library(crayon) # to colour some messages in the terminal
library(stringr) # to handle character strings
library(lme4) # for hierarchical models
library(dfoptim) # for Nelder-Mead algorithm for derivative-free optimization: function 'nmk'
library(parallel) # for parallel calculations
library(doParallel) # for parallel calculations
library(foreach) # for parallel calculations with function '%dopar%'
#library(bigparallelr) # for function nb_cores (recommended number of cores to use)

# Deterministic model using a system of ordinary differential equations (ODE)

lambda_ODE <- function(t, y, parms){ # Model with compartment Y
  # t = current time
  # y = current values of the state variables (at time t)
  # parms = set of parameters
  
  # (1) State variables
  
  # Two viral strains: subscript R -> Resident strain (WT)
  #                    subscript M -> Mutant strain
  
  S <- y[["S"]] # S(t), density of Susceptible bacteria
  
  # L(t), density of bacteria in a lysogenic state
  L_w <- y[["L_w"]]
  L_m <- y[["L_m"]]
  L <- L_w + L_m
  
  # Y(t), densities of infected bacteria prior to lysis
  Y_w <- y[["Y_w"]]
  Y_m <- y[["Y_m"]]
  Y <- Y_w + Y_m
  
  # V(t), density of Virions (free viral particles)
  V_w <- y[["V_w"]]
  V_m <- y[["V_m"]]

  N <- S + L + Y # Total density of bacteria
  
  # (2) Results
  
  # New infections
  abVRS <- parms[["a"]]*parms[["b"]]*V_w*S
  abVMS <- parms[["a"]]*parms[["b"]]*V_m*S
  
  # New reactivations of integrated phages / prophages (i.e., switches to a lytic cycle)
  reactivations_w <- parms[["alpha_w"]]*L_w
  reactivations_m <- parms[["alpha_m"]]*L_m
  
  # New lysis
  lysis_w <- parms[["tau"]]*Y_w
  lysis_m <- parms[["tau"]]*Y_m
  
  # Growth and removal
  growth <- parms[["r"]]*(1-N/parms[["K"]]) # Rate of logistic growth
  m <- parms[["m"]] # Removal rate of the chemostat
  removal.V <- parms[["a"]]*N + m # Removal rate of virions
  # The term a*N represents the rate at which virions adsorbed on bacteria S, L and Y (adsorption in non reversible).
  
  # Temporal derivatives
  
  dS <- growth*S - (abVRS + abVMS) - m*S
  
  dL_w <- growth*L_w + parms[["phi_w"]]*abVRS - reactivations_w - m*L_w
  dL_m <- growth*L_m + parms[["phi_m"]]*abVMS - reactivations_m - m*L_m
  
  dY_w <- (1-parms[["phi_w"]])*abVRS + reactivations_w - lysis_w - m*Y_w
  dY_m <- (1-parms[["phi_m"]])*abVMS + reactivations_m - lysis_m - m*Y_m
  
  dV_w <- lysis_w*parms[["B"]] - removal.V*V_w
  dV_m <- lysis_m*parms[["B"]] - removal.V*V_m
  
  return(list(c(dS, dL_w, dL_m, dY_w, dY_m, dV_w, dV_m)))
}

# Simulating models

remove_01 <- function(x, lim=1e-09){ # Set the values close to 0 or 1 to almost 0 or 1
  
  x[x <= 0] <- lim
  x[x >= 1] <- 1-lim
  
  return(x)
}

lambda_simul <- function(prev0, freq0=0.5, parms, times, solver="lsoda", verbose=TRUE){
  # prev0 = initial value for prevalence
  # freq0 = initial value for frequency of the virulent phage within lysogenized bacteria
  #         (default: 0.5, i.e. initial ratio 1:1)
  # parms = vector of model parameters (with names)
  
  # times = vector of desired time points for the simulation
  # solver = argument 'method' of function 'ode', i.e. the integrator to use
  # verbose = logical indicating whether the messages of the function should be printed (defaut: TRUE)
  
  parm_names <- names(parms)
  
  # (1) Parameters of phenotypic traits

  if(verbose){
    "=" %>% rep(66) %>% c("\n") %>% paste(collapse = '') %>% cat
    cat("  MODEL - No plasticity on any phenotypic traits of phage lambda\n")
    "=" %>% rep(66) %>% c("\n") %>% paste(collapse = '') %>% cat
    " " %>% rep(41) %>% c("Initial prevalence = ", round(100*prev0,2), " %\n", rep(" ", 40),
                          "(with ratio R:M = ", round(1-freq0, 2), ":", round(freq0, 2), ")\n") %>%
      paste(collapse = '') %>% cat
  }
  
  # (2) Remaining parameters
  
  if(!"K" %in% parm_names){parms["K"] <- 1e+09}
  if(!"m" %in% parm_names){parms["m"] <- 0.8}

  # (3) Initial conditions
  
  init <- (c("S" = (1-prev0), "L_w" = (1-freq0)*prev0, "L_m" = freq0*prev0)*parms[["K"]]) %>%
    c("Y_w" = 0, "Y_m" = 0, "V_w" = 0, "V_m" = 0)
  
  # (4) Simulation
  
  simul <- ode(y=init, times=times, parms=parms, func=lambda_ODE, method=solver) %>% as.data.frame
  
  if(verbose){cat("\nSimulation performed using a system of ordinary differential equations\n\n")}
  
  # Densities and prevalence
  
  simul$Y <- simul$Y_w + simul$Y_m
  simul$L <- simul$L_w + simul$L_m
  simul$Infected <- simul$Y + simul$L
  
  simul$Total_density <- simul$S + simul$Infected
  simul$Prevalence <- simul$Infected/simul$Total_density
  
  # Frequencies of virulent phage
  simul$p <- simul$L_m/simul$L # p(t), among lysogenic bacteria (L)
  simul$q <- simul$V_m/(simul$V_w + simul$V_m) # q(t), in free virus stage (V)
  simul$f <- simul$Y_m/simul$Y # f(t), among Y bacteria
  simul$g <- (simul$Y_m + simul$L_m)/simul$Infected # f^{Y+L}(t), among infected bacteria (both Y and L)
  
  # Logit-frequencies of virulent phage
  simul$logit_prev <- simul$Prevalence %>% remove_01 %>% logit
  simul$logit_p <- simul$p %>% remove_01 %>% logit
  simul$logit_q <- simul$q %>% remove_01 %>% logit
  simul$logit_f <- simul$f %>% remove_01 %>% logit
  simul$logit_g <- simul$g %>% remove_01 %>% logit
  # (remove_01 prevents the generation of infinite numbers by the logit function)
  
  return(simul)
}

# Plotting simulations

plot_simul <- function(simul, init_cond="N/D", accuracy=0.1){
  # simul = data.frame of the simulation - e.g. as returned by the above function lambda_simul()
  # col_epidemic (resp. col_endemic) = color for simulations with the epidemic (resp. endemic) treatment
  #                                    by default, col_epidemic is red and col_endemic is blue.
  
  # Only if 'simul' does not already have an 'init_cond' column:
  # init_cond = "epidemic", "endemic" or N/D (not defined, default)
  
  t <- unique(simul$time)
  breaks <- seq(t[1], t[length(t)], 10)
  
  colors <- c("epidemic"='#DA0F0F', "endemic"='#0E5DD6', "N/D"='#000000')
  
  theme_custom <- theme_bw() + theme(axis.text.x = element_text(size=11),
                                     axis.text.y = element_text(size=11),
                                     legend.title = element_blank())
  if(init_cond == "epidemic"){
    simul$init_cond <- "epidemic"
    colors <- colors["epidemic"]
  }else if(init_cond == "endemic"){
    simul$init_cond <- "endemic"
    colors <- colors["endemic"]
  }else{
    if("init_cond" %in% colnames(simul)){
      init_cond <- simul$init_cond %>% unique
      colors <- colors[init_cond %>% match(names(colors))]
    }else{
      simul$init_cond <- "N/D"
      colors <- colors["N/D"]
    }
  }
  if("N/D" %in% init_cond){theme_custom <- theme_custom + theme(legend.position = 'none')}

  return(list(
    "Fig_logit_prevalence" = ggplot(data=simul, aes(x=time, y=logit_prev, color=init_cond)) +
      geom_line(cex = 1.1) +
      scale_color_manual(values = colors) +
      labs(x="time (hours)", y="Logit-prevalence\n") +
      scale_x_continuous(breaks = breaks) +
      scale_y_continuous(labels = scales::label_number(accuracy=accuracy)) +
      theme_custom,
    
    "Fig_logit_inf" = ggplot(data=simul, aes(x=time, y=logit_g, color=init_cond)) +
      geom_line(cex = 1.1) +
      scale_color_manual(values = colors) +
      labs(x="time (hours)", y="Logit-frequency infected by virulent phage (logit scale)\n") +
      scale_x_continuous(breaks = breaks) +
      scale_y_continuous(labels = scales::label_number(accuracy=accuracy)) +
      theme_custom,
    
    "Fig_logit_q" = ggplot(data=simul, aes(x=time, y=logit_q, color=init_cond)) +
      geom_line(cex = 1.1) +
      scale_color_manual(values = colors) +
      labs(x = "time (hours)", y="Logit-frequency free virulent phage\n") +
      scale_x_continuous(breaks = breaks) +
      scale_y_continuous(labels = scales::label_number(accuracy=accuracy)) +
      theme_custom
  ))
}

# Compute basic reproduction number

R0 <- function(alpha, phi, a, b, tau, B, r, S0, K=1E+09, delta=0.8){
  growth <- r*(1-S0/K)
  abStauB <- a*b*S0*tau*B
  A <- growth*(a*S0+delta)*(tau+delta)+abStauB*(alpha+(1-phi)*delta)
  D <- (a*S0+delta)*(alpha+delta)*(tau+delta)
  return((A+sqrt(A^2-4*growth*(1-phi)*abStauB*D))/(2*D))
}

# Inference

estim_alpha <- function(data, threshold_prev=0.95, by_chemostat=TRUE){
  
  tab <- data %>% ddply(~chemostat, function(X){
    tmin <- X$time[which(inv.logit(X$logit_prev)>=threshold_prev)[1]]
    return(subset(X, time>=tmin))})
  tab$lnQ <- tab$logit_q-tab$logit_g
  if(by_chemostat){
    estim <- tab %>% ddply(~chemostat, function(X){
      return(c("prev0"=unique(X$prev0), "Q"=exp(mean(X$lnQ, na.rm=TRUE)), # geometric mean
               "Q_gsem"=exp(sd(X$lnQ, na.rm=TRUE))/sqrt(sum(!is.na(X$lnQ)))))}) # geometric std err. of the mean
    mem <- lmer(logit_g ~ time:chemostat + as.factor(prev0) + (1|chemostat), tab)
    Delta_alpha_estim <- cbind("Delta_alpha" = fixef(mem) %>% {-.[grep("time", names(.))]},
                               # Delta_alpha = alpha_m-alpha_w = -slope
                               "Delta_alpha_sd" = mem %>% vcov %>% diag %>% sqrt %>%
                                 {.[grep("time", names(.))]}) # Standard errors
    estim <- cbind(estim, Delta_alpha_estim[match(rownames(Delta_alpha_estim),
                                                  paste0('time:chemostat', estim$chemostat)),])
    rownames(estim) <- NULL
  }else{
    estim <- data.frame("chemostat"="all", "prev0"=NA, "Q"=exp(mean(tab$lnQ, na.rm=TRUE)), # geometric mean
                        "Q_gsem"=exp(sd(tab$lnQ, na.rm=TRUE))/sqrt(sum(!is.na(tab$lnQ)))) # geometric std err. of the mean
    mem <- lmer(logit_g ~ time + as.factor(prev0) + (1|chemostat), tab)
    Delta_alpha_estim <- cbind("Delta_alpha" = -fixef(mem)[["time"]],
                               # Delta_alpha = alpha_m-alpha_w = -slope
                               "Delta_alpha_sd" = sqrt(diag(vcov(mem)))[["time"]]) # Standard errors
    estim <- cbind(estim, Delta_alpha_estim)
  }
  estim$alpha_w <- estim$Delta_alpha/(estim$Q-1)
  estim$alpha_m <- estim$Q*estim$alpha_w
  return(estim)
}

## Parameter transformations (and reversal)

bounds <- data.frame(
  
  "prev0_epidemic" = c(1e-03, 2e-02),
  "prev0_endemic" = c(0.95, 1),
  "freq0" = c(0.4, 0.6),
  
  "phi_w" = c(0, 1),
  "phi_m" = c(0, 1),
  
  "a" = c(1e-09, 1e-06),
  "b" = c(0, 0.2),
  "B" = c(10, 500),
  "r" = c(0, 5),
  
  "tau" = c(0.5, 4),
  
  "sigma_prev" = c(1e-03, 2),
  "sigma_g" = c(1e-03, 2),
  "sigma_q" = c(1e-03, 2))

rownames(bounds) <- c("lower", "upper")

generalized_logit <- function(x, bounds){ # ]bounds[1], bounds[2][ -> ]-Inf, +Inf[
  # x is a number (or a vector of numbers) between bounds[1] (min) and bounds[2] (max)
  
  if(length(x) == 1){
    if(is.na(x)){
      return(NA)
    }else if(x < bounds[1] | x > bounds[2]){
      stop(paste0("x (= ", x, ") is not between a (= ", bounds[1], ") and b (= ", bounds[2],")..."))
    }else{
      return(log((x-bounds[1])/(bounds[2]-x)))
    }
  }else if(length(x) > 1){
    return(sapply(x, FUN=generalized_logit, bounds=bounds))
  }else{
    return(numeric())
  }
}

inv.generalized_logit <- function(y, bounds){ # R (]-Inf, +Inf[) -> ]bounds[1], bounds[2][
  
  if(length(y) == 1){
    
    num <- bounds[1]+bounds[2]*exp(y)
    denom <- 1+exp(y)
    
    if(is.infinite(num) | is.infinite(denom)){
      if(y %>% names %>% is.null){
        return(bounds[2])
      }else{
        return(setNames(object = bounds[2], nm = names(y)))
      }
    }else{
      return(num/denom)
    }
  }else if(length(y) > 1){
    return(sapply(y, FUN = inv.generalized_logit, bounds = bounds))
  }else{
    return(numeric())
  }
}

parm_transformation <- function(parms, reverse=FALSE){
  
  if(data.class(parms) == "numeric"){
    
    parm_names <- names(parms)
    
    if(reverse == FALSE){
      transf <- c(parms[parm_names == "prev0_epidemic"] %>% generalized_logit(bounds$prev0_epidemic),
                  parms[parm_names == "prev0_endemic"] %>% generalized_logit(bounds$prev0_endemic),
                  parms[parm_names == "freq0"] %>% generalized_logit(bounds$freq0),

                  parms[parm_names == "phi_w"] %>% generalized_logit(bounds$phi_w),
                  parms[parm_names == "phi_m"] %>% generalized_logit(bounds$phi_m),

                  parms[parm_names == "a"] %>% generalized_logit(bounds$a),
                  parms[parm_names == "b"] %>% generalized_logit(bounds$b),
                  parms[parm_names == "B"] %>% generalized_logit(bounds$B),
                  parms[parm_names == "r"] %>% generalized_logit(bounds$r),

                  parms[parm_names == "tau"] %>% generalized_logit(bounds$tau),
                  
                  parms[parm_names == "sigma_prev"] %>% generalized_logit(bounds$sigma_prev),
                  parms[parm_names == "sigma_g"] %>% generalized_logit(bounds$sigma_g),
                  parms[parm_names == "sigma_q"] %>% generalized_logit(bounds$sigma_q))

      return(transf[parm_names %>% match(names(transf))])
      
    }else{
      inv.transf <- c(parms[parm_names == "prev0_epidemic"] %>% inv.generalized_logit(bounds$prev0_epidemic),
                      parms[parm_names == "prev0_endemic"] %>% inv.generalized_logit(bounds$prev0_endemic),
                      parms[parm_names == "freq0"] %>% inv.generalized_logit(bounds$freq0),
                      
                      parms[parm_names == "phi_w"] %>% inv.generalized_logit(bounds$phi_w),
                      parms[parm_names == "phi_m"] %>% inv.generalized_logit(bounds$phi_m),
                      
                      parms[parm_names == "a"] %>% inv.generalized_logit(bounds$a),
                      parms[parm_names == "b"] %>% inv.generalized_logit(bounds$b),
                      parms[parm_names == "B"] %>% inv.generalized_logit(bounds$B),
                      parms[parm_names == "r"] %>% inv.generalized_logit(bounds$r),
                      
                      parms[parm_names == "tau"] %>% inv.generalized_logit(bounds$tau),
                      
                      parms[parm_names == "sigma_prev"] %>% inv.generalized_logit(bounds$sigma_prev),
                      parms[parm_names == "sigma_g"] %>% inv.generalized_logit(bounds$sigma_g),
                      parms[parm_names == "sigma_q"] %>% inv.generalized_logit(bounds$sigma_q))
      
      return(inv.transf[parm_names %>% match(names(inv.transf))])
    }
  }else{ # data.frame or matrix
    return(parms %>% apply(1, parm_transformation, reverse=reverse) %>% t)
  }
}

draw_initial_conditions <- function(n, seed=123){
  set.seed(seed)
  return(data.frame("prev0_epidemic" = n %>% runif(bounds$prev0_epidemic[1], bounds$prev0_epidemic[2]),
                    "prev0_endemic" = n %>% runif(bounds$prev0_endemic[1], bounds$prev0_endemic[2]),
                    "freq0" = n %>% runif(bounds$freq0[1], bounds$freq0[2]),
                    "phi_w" = n %>% runif(bounds$phi_w[1], bounds$phi_w[2]),
                    "phi_m" = n %>% runif(bounds$phi_m[1], bounds$phi_m[2]),
                    "a" = n %>% runif(bounds$a[1], bounds$a[2]),
                    "b" = n %>% runif(bounds$b[1], bounds$b[2]),
                    "r" = n %>% runif(bounds$r[1], bounds$r[2]),
                    "tau" = n %>% runif(bounds$tau[1], bounds$tau[2]),
                    "sigma_prev" = n %>% runif(bounds$sigma_prev[1], bounds$sigma_prev[2]),
                    "sigma_g" = n %>% runif(bounds$sigma_g[1], bounds$sigma_g[2]),
                    "sigma_q" = n %>% runif(bounds$sigma_q[1], bounds$sigma_q[2])))
}

## Negative log-likelihood

negative_log_likelihood <- function(transf_parms, parm_names=NULL, constants=NULL,
                                    data_FACS, data_qPCR, solver="lsoda", verbose=TRUE, return_all=FALSE){
  
  # transf_parms: Transformed parameters to estimate (see parameter transformation)
  #               (NB: this allow to enforce parameter bounds)
  #   parm_names: Vector of names for transf_parms
  #    constants: Vector of named parameters to set as constants (may include initial values)
  #   return_all: (logical) Should all results be returned (TRUE) or only the overall negative log-likelihood (FALSE)?
  #    data_FACS: Data frame with experimental data obtained by FACS. Chemostats, initial condition / treatment
  #               (1, "epidemic", or 100, "endemic"), logit(prevalence) and logit(frequency infected by virulent phage)
  #               must each be specified in separates columns named 'chemostat', 'prev0', 'logit_prev' and 'logit_g',
  #               respectively.
  #    data_qPCR: Data frame with experimental data obtained by qPCR. Chemostats, initial condition / treatment
  #               (1, "epidemic", or 100, "endemic"), and logit(frequency free virulent phage) must each be specified in
  #               separates columns named 'chemostat', 'prev0' and 'logit_q', respectively.
  #      verbose: See function 'lambda_simul' (cf. above)
  
  names(transf_parms) <- parm_names
  
  #   if(is.null(parm_names)){
  #     stop("Parameters in argument 'transf_parms' are unnamed
  #     (note that this can happen with the function 'nmk(b)' which gets rid of these names at each iteration)
  #     To solve this, you may specify parameter names (/!\ in the same order as transf_parms) using the argument 'parm_names'")
  #   }
  
  # (1) Reversal of parameter transformations
  
  parms <- transf_parms %>% parm_transformation(reverse = TRUE) %>% c(constants)
  
  parm_names <- names(parms)
  
  # (2) Numerical simulations
  
  init_cond <- c(data_FACS$prev0, data_qPCR$prev0) %>% unique
  times <- c(data_FACS$time, data_qPCR$time) %>% unique %>% sort
  freq0 <- ifelse("freq0" %in% parm_names, yes = parms[["freq0"]], no = 0.5)
  
  if(1 %in% init_cond){ ## 1. Epidemic
    simul_epidemic <- lambda_simul(prev0=parms[["prev0_epidemic"]], freq0=freq0, parms=parms,
                                   times=times, solver=solver, verbose=verbose)
    if((!times %in% simul_epidemic$time %>% all) | (simul_epidemic[-1,] %>% is.na %>% any)){
      if(return_all){message("Message: simulation failed; detailed results can't be computed.")}
      return(1e+10)
    }
  }
  if(100 %in% init_cond){ ## 2. Endemic
    simul_endemic <- lambda_simul(prev0=parms[["prev0_endemic"]], freq0=freq0, parms=parms,
                                  times=times, solver=solver, verbose=verbose)
    if((!times %in% simul_endemic$time %>% all) | (simul_endemic[-1,] %>% is.na %>% any)){
      if(return_all){message("Message: simulation failed; detailed results can't be computed.")}
      return(1e+10)
    }
  }
  
  # (3) Compute negative log-likelihood (nLL)
  
  # for the logit-prevalence
  nLL_prev <- data_FACS %>% plyr::ddply(~chemostat, function(X){
    prev0 <- unique(X$prev0)
    if(prev0 == 1){
      simul <- simul_epidemic[X$time %>% match(simul_epidemic$time),]
    }else if(prev0 == 100){
      simul <- simul_endemic[X$time %>% match(simul_endemic$time),]
    }else{
      stop("Elements in column 'prev0' of data frame 'data_FACS' are neither 1 nor 100")
    }
    return(data.frame("prev0" = prev0, "time" = X$time,
                      "value" = -dnorm(X$logit_prev, mean=simul$logit_prev, sd=parms[["sigma_prev"]], log=TRUE)))
  })
  # for the logit-frequency of the virulent strain among infected cells
  nLL_g <- data_FACS %>% plyr::ddply(~chemostat, function(X){
    prev0 <- X$prev0 %>% unique
    if(prev0 == 1){
      simul <- simul_epidemic[X$time %>% match(simul_epidemic$time),]
    }else if(prev0 == 100){
      simul <- simul_endemic[X$time %>% match(simul_endemic$time),]
    }else{
      stop("Elements in column 'prev0' of data frame 'data_FACS' are neither 1 nor 100")
    }
    return(data.frame("prev0" = prev0, "time" = X$time,
                      "value" = -dnorm(X$logit_g, mean=simul$logit_g, sd=parms[["sigma_g"]], log=TRUE)))
  })
  # for the logit-frequency of the virulent strain in the free virus stage
  nLL_q <- data_qPCR %>% plyr::ddply(~chemostat, function(X){
    prev0 <- X$prev0 %>% unique
    if(prev0 == 1){
      simul <- simul_epidemic[X$time %>% match(simul_epidemic$time),]
    }else if(prev0 == 100){
      simul <- simul_endemic[X$time %>% match(simul_endemic$time),]
    }else{
      stop("Elements in column prev0 of data frame 'data_FACS' are neither 1 nor 100")
    }
    return(data.frame("prev0" = prev0, "time" = X$time,
                      "value" = -dnorm(X$logit_q, mean=simul$logit_q, sd=parms[["sigma_q"]], log=TRUE)))
  })
  # indicator of bacterial population getting extinct in simulations
  simul_extinctions <- any(simul_epidemic$Total_density < 1e+06) + any(simul_endemic$Total_density < 1e+06)
  
  if(return_all){
    return(list("negative_LogLikelihood_prev" = nLL_prev, "negative_LogLikelihood_g" = nLL_g,
                "negative_LogLikelihood_q" = nLL_q, "Bacterial_extinction_simulation" = simul_extinctions))
  }else{
    # Overall negative log-likelihood
    return(sum(nLL_prev$value, nLL_g$value, nLL_q$value, (1e+10)*simul_extinctions, na.rm = TRUE))
  }
}

# Estimations with non linear optimizations

Estim_nmk_parallel <- function(parms0, constants=NULL, data_FACS, data_qPCR,
                               solver="lsoda", maxfeval=NA, tol=1e-6, restarts.max=3, Ncores=NA){
  # parms0 = a vector or matrix (when multiple starts) with initial values for the parameter to estimate
  # maxfeval, tol and restarts.max = arguments passed to function 'nmk' (see ?dfoptim::nmk)
  # Ncores = number of cores for parallelized optimizations
  #
  # See function 'negative_log_likelihood' for details of the remaining arguments (cf. above)
  
  # (1) Transformations of parameters
  
  transf_parms0 <- parm_transformation(parms0)
  
  # (2) Initializing
  
  if(data.class(transf_parms0) == "numeric"){ # if transf_parms0 is a vector
    transf_parms0 <- transf_parms0 %>% t %>% as.matrix
  }
  parm_names <- colnames(transf_parms0)
  n_starts <- nrow(transf_parms0)
  n_parms <- ncol(transf_parms0)
  
  if(is.na(maxfeval)){maxfeval <- min(5000, max(1500, 20*n_parms^2))} # Default in 'nmk'
  
  data_qPCR <- data_qPCR[-(data_qPCR$logit_q %>% is.na %>% which),]
  
  # (3) Parallelized optimizations
  
  if(is.na(Ncores)){Ncores <- bigparallelr::nb_cores()}
  print(paste("Number of cores:", Ncores))
  
  t <- format(Sys.time(), "%Y-%m-%d %Hh%M")
  folder_name <- paste0("Lambda_v11_Estim_nmk_parallel_", stringi::stri_replace_all_fixed(t, c("-"," "),
                                                                                          "_", vectorize_all=FALSE))
  dir.create(folder_name)
  lapply(list("PROJECT LAMBDA (script v11 - no phenotypic plasticity)\n",
              t, paste0("\nParallel - number of cores: ", Ncores),
              "\nOptimization algorithm: nmk()", "\nControl:",
              paste(c("Max function evaluations", "Tolerance", "Max restarts"),
                    c(maxfeval,tol,restarts.max), sep=" = ", collapse="\n"),
              paste("\nNumber of starts:", n_starts),
              paste0("\nParameters to estimate (", n_parms, "):"), parm_names %>% paste(collapse=", "),
              "\nFixed parameters:", paste(names(constants), constants, sep=" = ", collapse="\n")),
         cat, "\n", file=paste0(folder_name,"/README.txt"), append=TRUE)
  
  cl <- parallel::makeCluster(Ncores)
  doParallel::registerDoParallel(cl)
  
  estim_df <- foreach::foreach(i = 1:n_starts, .combine = 'rbind', .inorder = FALSE,
                               .packages = c("tidyverse", "deSolve", "boot"), .export = ls(globalenv())) %dopar% {
                                 
                                 fit <- dfoptim::nmk(par = transf_parms0[i,], fn = negative_log_likelihood,
                                                     control = list(maximize = FALSE, maxfeval=maxfeval,
                                                                    tol=tol, restarts.max=restarts.max),
                                                     parm_names=parm_names, constants=constants,
                                                     data_FACS=data_FACS, data_qPCR=data_qPCR,
                                                     solver=solver, verbose=FALSE, return_all=FALSE)
                                 
                                 res <- c("Start"=i, "Value"=fit$value, "Convergence"=fit$convergence,
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
  
  return(estim_df)
}

## Create a data frame with results of RDS files (from Estim_nmk_parallel())

data_frame_from_RDS_files.nmk <- function(folder){
  # folder = working directory to the folder containing RDS files
  
  curr_wd <- getwd()
  setwd(folder)
  
  df <- data.frame()
  
  files <- dir()
  
  if(any(files=="README.txt")){files <- files[-which(files=="README.txt")]}
  
  for(file in files){
    df <- rbind(df, readRDS(file))
  }
  colnames(df) <- names(readRDS(file[1]))
  df <- df[order(df$Start),]
  
  df$Convergence[df$Convergence == 0] <- "Successful completion"
  df$Convergence[df$Convergence == 1] <- "Maximum number of evaluations exceeded"
  df$Convergence[df$Convergence == 2] <- "Stagnation in Nelder-Mead"
  # df$Convergence[df$Convergence == 10] <- "Degeneracy of the Nelder-Mead simplex"
  
  setwd(curr_wd)
  
  return(df)
}

# Return best estimates with successful completion and provides graphical vizualisations

process_estim <- function(estim, range=10){
  
  # estim = data frame, as returned by functions Estim_nmk_parallel
  # range = numeric argument to select estimates whose negative log-likelihood are within 'range' of the best one
  
  # Best parameters
  
  estim_success <- estim %>% subset(Convergence == "Successful completion")
  
  best_index <- which.min(estim_success$Value)
  best <- estim_success[best_index,]
  
  estim.range <- estim[abs(estim$Value-best$Value) <= range, ]
  # with successful completion between the best -log-likelihood and the best -log-likelihood + range
  
  # Plot negative log-likelihood for each start with convergence features
  
  color_conv <- c("Successful completion" = "#00BFC4",
                  "Maximum number of evaluations exceeded" = "#F8766D",
                  "Iteration limit had been reached" = "#F8766D",
                  "Stagnation in Nelder-Mead" = "#FF9E00",
                  "Degeneracy of the Nelder-Mead simplex" = "#C77CFF")
  
  # data$Convergence[!data$Convergence %in% names(color_conv)] <- "Other"
  # color_conv[["Other"]] <- "#A2A2A2"
  
  Fig_results_convergence <- ggplot(estim, aes(x = Start, y = Value)) +
    geom_line(cex = 0.1, alpha = 0.1) +
    geom_point(aes(col = Convergence), cex = 1, pch = 19) +
    geom_point(data = best, cex = 5, pch = 1) +
    labs(x = "# start", y = "- log-likelihood\n", col = "Convergence:") +
    scale_color_manual(values = color_conv[estim$Convergence %>% unique %>%
                                             match(names(color_conv))]) +
    # scale_y_log10() +
    theme_classic() + theme(axis.text.x = element_text(size=11),
                            axis.text.y = element_text(size=11),
                            legend.position = "bottom")
  
  # For each parameter, plot estimation against negative log-likelihood
  
  Fig_parms <- (4:ncol(estim)) %>% lapply(FUN = function(i){
    
    parm_name <- colnames(estim)[i]
    
    X <- estim[,c(1:3,i)]
    # X[!X$Start %in% estim_success$Start,4] <- NA
    
    # Fig.1 <- ggplot(X, aes(x = Start, y = X[,4])) +
    #   geom_hline(yintercept = bounds[,parm_name] %>% as.numeric) +
    #   geom_line(aes(col = Value), cex = 0.2, alpha = 0.3) +
    #   geom_point(aes(col = Value)) +
    #   geom_point(data = estim.range, aes(x = Start, y = estim.range[,i]),
    #              cex = 5, pch = 1, col = '#35c3d3') +
    #   geom_point(data = best, aes(x = Start, y = best[,i]), cex = 5, pch = 1) +
    #   scale_color_gradientn(colours = brewer.pal(9, "YlOrRd")[-(1:2)] %>% rev) +
    #   labs(x = "# start", y = paste0(parm_name,"\n"), col = "- log-Likelihood\n") +
    #   theme_classic() + theme(axis.text.x = element_text(size=11),
    #                           axis.text.y = element_text(size=11))
    
    Fig.2 <- ggplot(X, aes(x = X[,4], y = Value)) +
      geom_vline(xintercept = bounds[,parm_name] %>% as.numeric) +
      geom_point(aes(col = Convergence), cex = 0.5, alpha = 0.3) +
      geom_point(data = best, aes(x = best[,i], y = Value), cex = 4, pch = 1) +
      scale_color_manual(values = color_conv[X$Convergence %>% unique %>%
                                               match(names(color_conv))]) +
      labs(x = parm_name, y = "- log-likelihood\n") +
      theme_classic() + theme(axis.text.x = element_text(size=11),
                              axis.text.y = element_text(size=10),
                              legend.position = "none")
    
    Fig.3 <- ggplot(estim.range, aes(x = estim.range[,i], y = Value)) +
      geom_vline(xintercept = bounds[,parm_name] %>% as.numeric) +
      geom_point(aes(col = Convergence)) +
      geom_point(data = best, aes(x = best[,i], y = Value), cex = 5, pch = 1) +
      scale_color_manual(values = color_conv[X$Convergence %>% unique %>%
                                               match(names(color_conv))]) +
      labs(x = parm_name, y = "- log-likelihood\n") +
      theme_classic() + theme(axis.text.x = element_text(size=11),
                              axis.text.y = element_text(size=9),
                              legend.position = "none")
    
    return(list(Fig.2, Fig.3))
  })
  names(Fig_parms) <- colnames(estim)[-(1:3)]
  
  return(list("best" = best, "success" = estim_success,
              "Fig_results_convergence" = Fig_results_convergence, "Fig_parms" = Fig_parms))
}

# Functions for Sieve Bootstrap

# Simulate ARMA model

simul_ARMA <- function(arima.estim, n=100, n.burnin=0, y0=0){
  
  # arima.estim: result of forecast::auto.arima()
  # n: number of returned values
  # n.burnin: length of the burn-in period (values not returned)
  # y0: initial value
  
  coef_arma <- coef(arima.estim)
  s <- sqrt(arima.estim$sigma2)
  
  if(length(coef_arma)==0){
    return(c(y0, rnorm(n.burnin-1+n, mean=0, sd=s))[(n.burnin+1):(n.burnin+n)])
  }
  
  ar <- coef_arma[grep("ar", names(coef_arma))]
  ma <- coef_arma[grep("ma", names(coef_arma))]
  m <- ifelse("intercept" %in% names(coef_arma), yes=coef_arma[["intercept"]], no=0)
  
  p <- length(ar)
  q <- length(ma)
  
  if(p){ # check for stationarity (taken from arima.sim())
    minroots <- min(Mod(polyroot(c(1, -ar))))
    if (minroots <= 1)
      stop("'ar' part of model is not stationary")
  }
  max_order <- max(p,q)
  
  ysim <- c(rep(0, max_order-1), y0)
  e <- c(rep(0, max_order), rnorm(n.burnin-1+n, mean=0, sd=s)) # white noise
  
  for(t in (max_order+1):(max_order+n.burnin+n-1)){
    ysim[t] <- m + sum(ysim[(t-1):(t-p)]*ar) + e[t] + sum(e[(t-1):(t-q)]*ma)
  }
  return(ysim[(max_order+n.burnin):(max_order+n.burnin+n-1)])
}

# Generate bootstrapped data

bootstrapped_data <- function(data, fit, L, seed=123){
  
  # data: experimental data from FACS and qPCR
  # fit: simulations based on the best set of MLE parameter estimates
  # L: a list with, for each chemostat, chemostat's name, mean of residual and output of forecast::auto.arima()
  # seed: seed for set.seed()
  
  data.bs <- data %>%
    
    ddply(~chemostat, function(X){
      
      # Fit
      treatment <- unique(X$prev0)
      fit <- subset(fit, prev0==treatment)
      time_index <- match(X$time, fit$time)
      
      X$logit_prev <- fit$logit_prev[time_index]
      X$logit_g <- fit$logit_g[time_index]
      X$logit_q <- fit$logit_q[time_index]
      
      X$logit_prev[is.na(X$prevalence)] <- NA
      X$logit_g[is.na(X$g)] <- NA
      X$logit_q[is.na(X$q)] <- NA
      
      # Bootstrapped residuals
      arma <- L[names(L)==unique(X$chemostat)][[1]]
      resid_time <- 1:60
      time_index <- match(X$time, resid_time)
      n <- length(resid_time)
      
      set.seed(seed)
      resid.logit_prev <- simul_ARMA(arma$logit_prev, n=n, n.burnin = 10) + arma$means[1]
      resid.logit_g <- simul_ARMA(arma$logit_g, n=n, n.burnin = 10) + arma$means[2]
      resid.logit_q <- simul_ARMA(arma$logit_q, n=n, n.burnin = 10) + arma$means[3]
      
      X$logit_prev <- X$logit_prev + resid.logit_prev[time_index]
      X$logit_g <- X$logit_g + resid.logit_g[time_index]
      X$logit_q <- X$logit_q + resid.logit_q[time_index]
      
      return(X[,colnames(X) %in% c("chemostat", "time", "prev0", "logit_prev", "logit_g", "logit_q")])
    })
  return(data.bs)
}

Estim_nmk_parallel_bootstrap <- function(best_parms, constants=NULL, bootstrapped.data,
                                         solver="lsoda", maxfeval=NA, tol=1e-6, restarts.max=3, Ncores=NA){
  # best_parms = best set of MLE parameter estimates
  # maxfeval, tol and restarts.max = arguments passed to function 'nmk' (see ?dfoptim::nmk)
  # Ncores = number of cores for parallelized optimizations
  #
  # See function 'negative_log_likelihood' for details of the remaining arguments (cf. above)
  
  # (1) Transformations of parameters
  
  transf_parms0 <- parm_transformation(best_parms)
  
  # (2) Initializing
  
  n_datasets <- length(unique(bootstrapped.data$id))
  parm_names <- colnames(transf_parms0)
  n_parms <- ncol(transf_parms0)
  
  transf_parms0 <- transf_parms0 %>% as.vector %>% setNames(parm_names)
  
  if(is.na(maxfeval)){maxfeval <- min(5000, max(1500, 20*n_parms^2))} # Default in 'nmk'
  
  # (3) Parallelized optimizations
  
  if(is.na(Ncores)){Ncores <- bigparallelr::nb_cores()}
  print(paste("Number of cores:", Ncores))
  
  t <- format(Sys.time(), "%Y-%m-%d %Hh%M")
  folder_name <- paste0("Lambda_v11_Estim_nmk_parallel_bootstrap_",
                        stringi::stri_replace_all_fixed(t, c("-"," "), "_", vectorize_all=FALSE))
  dir.create(folder_name)
  lapply(list("PROJECT LAMBDA (script v11 - no phenotypic plasticity)\n",
              t, paste0("\nParallel - number of cores: ", Ncores),
              "\nOptimization algorithm: nmk()", "\nControl:",
              paste(c("Max function evaluations", "Tolerance", "Max restarts"),
                    c(maxfeval,tol,restarts.max), sep=" = ", collapse="\n"),
              paste("\nNumber of bootstrapped datasets:", n_datasets),
              paste0("\nParameters to estimate (", n_parms, "):"), parm_names %>% paste(collapse=", "),
              "\nFixed parameters:", paste(names(constants), constants, sep=" = ", collapse="\n")),
         cat, "\n", file=paste0(folder_name,"/README.txt"), append=TRUE)
  
  cl <- parallel::makeCluster(Ncores)
  doParallel::registerDoParallel(cl)
  
  estim_df <- foreach::foreach(i = 1:n_datasets, .combine = 'rbind', .inorder = FALSE,
                               .packages = c("plyr", "tidyverse", "deSolve", "boot", "lme4"), .export = ls(globalenv())) %dopar% {
                                 
                                 data.i <- bootstrapped.data %>% subset(id==i)
                                 
                                 alpha <- estim_alpha(data.i, by_chemostat = FALSE)
                                 alpha_w <- alpha[["alpha_w"]]
                                 alpha_m <- alpha[["alpha_m"]]
                                 
                                 fit <- dfoptim::nmk(par = transf_parms0, fn = negative_log_likelihood,
                                                     control = list(maximize = FALSE, maxfeval=maxfeval,
                                                                    tol=tol, restarts.max=restarts.max),
                                                     parm_names=parm_names, constants=c("alpha_w"=alpha_w,
                                                                                        "alpha_m"=alpha_m,
                                                                                        constants),
                                                     data_FACS=data.i, data_qPCR=data.i,
                                                     solver=solver, verbose=FALSE, return_all=FALSE)
                                 
                                 res <- c("id"=i, "Value"=fit$value, "Convergence"=fit$convergence,
                                          fit$par %>% setNames(parm_names) %>% parm_transformation(reverse=TRUE) %>%
                                            c("alpha_w"=alpha_w, "alpha_m"=alpha_m, .))
                                 saveRDS(res, file = paste0(folder_name, "/result_task_", i, ".rds"))
                                 return(res)
                               }
  parallel::stopCluster(cl)
  
  estim_df <- as.data.frame(estim_df)
  
  return(estim_df[order(estim_df$id),])
}

# Miscellaneous

## For graphical purpose: Add a '+' before positive numbers ('-' before negative numbers, nothing before 0)

add_plus_sign <- function(x, digits = 2){
  if(length(x) > 1){
    return(sapply(x, add_plus_sign, digits = digits))
  }
  if(isTRUE(all.equal(x,0))){
    return(sprintf(paste0("% .",digits,"f"),x))
  }else{
    return(sprintf(paste0("%+.",digits,"f"),x))
  }
}

scientific_10 <- function(x){
  return(sapply(x, function(x_i){
    if(is.na(x_i)){
      return("")
    }else if(x_i==0){
      return("0")
    }else{
      return(parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x_i))))
    }
  }))
}