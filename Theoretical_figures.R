rm(list=ls())

library(plyr) # for the pipe operator %>%
library(tidyverse) # for data manipulation, graphic visualizations, etc...
library(RColorBrewer) # for colors and palettes
library(cowplot) # to show multiple plots
library(ggpubr) # to show multiple plots
library(deSolve) # to solve a system of ordinary differential equations (ODE integration, function 'ode')
library(boot) # for the inverse logit function 'inv.logit'

# Ordinary differential equations (main model with constant or variable S)

lambda_ODE <- function(t, y, parms, S_cst){
  
  # State variables
  S <- y[["S"]]
  Lw <- y[["Lw"]]
  Lm <- y[["Lm"]]
  Yw <- y[["Yw"]]
  Ym <- y[["Ym"]]
  Vw <- y[["Vw"]]
  Vm <- y[["Vm"]]
  
  # Parameters
  alpha_w <- parms[["alpha_w"]]
  alpha_m <- parms[["alpha_m"]]
  phi_w <- parms[["phi_w"]]
  phi_m <- parms[["phi_m"]]
  tau <- parms[["tau"]]
  r <- parms[["r"]]
  K <- parms[["K"]]
  delta <- parms[["delta"]]
  a <- parms[["a"]]
  b <- parms[["b"]]
  B <- parms[["B"]]
  
  # Infections
  abVwS <- a*b*Vw*S
  abVmS <- a*b*Vm*S
  
  # Reactivations
  react_w <- alpha_w*Lw
  react_m <- alpha_m*Lm
  
  # Lysis
  lysis_w <- tau*Yw
  lysis_m <- tau*Ym
  
  # Total density
  N <- ifelse(S_cst, yes = K, no = S+Lw+Lm+Yw+Ym) 
  
  growth <- r*(1-N/K)

  # Temporal differentiation
  dS <- ifelse(S_cst, yes = 0, no = growth*S - (abVwS+abVmS) - delta*S)
  
  dLw <- growth*Lw + phi_w*abVwS - react_w - delta*Lw
  dLm <- growth*Lm + phi_m*abVmS - react_m - delta*Lm
  
  dYw <- (1-phi_w)*abVwS + react_w - lysis_w - delta*Yw
  dYm <- (1-phi_m)*abVmS + react_m - lysis_m - delta*Ym
  
  dVw <- lysis_w*B - (a*N + delta)*Vw
  dVm <- lysis_m*B - (a*N + delta)*Vm
  
  list(c(dS, dLw, dLm, dYw, dYm, dVw, dVm))
}

lambda_simul <- function(prev0, freq0=0.5, parms, times, S_cst=FALSE){
  
  # prev0 = initial value for prevalence
  # freq0 = initial value for frequency of the virulent phage within lysogenized bacteria
  #         (default: 0.5, i.e. initial ratio 1:1)
  # parms = vector of model parameters (with names)
  # times = vector of desired time points for the simulation
  
  init <- c("S"=(1-prev0)*parms[["K"]], "Lw"=(1-freq0)*prev0*parms[["K"]],
            "Lm"=freq0*prev0*parms[["K"]], "Yw"=0, "Ym"=0, "Vw"=0, "Vm"=0)
  
  simul <- lsoda(y=init, times=times, func=lambda_ODE, parms=parms, S_cst=S_cst) %>% as.data.frame
  
  simul$L <- simul$Lw + simul$Lm
  simul$Y <- simul$Yw + simul$Ym
  
  infected <- simul$L + simul$Y
  
  simul$V <- simul$Vw + simul$Vm
  simul$N <- simul$S + infected
  simul$prev <- infected/simul$N
  
  simul$p <- simul$Lm/simul$L
  simul$q <- simul$Vm/simul$V
  simul$f <- simul$Ym/simul$Y
  simul$g <- (simul$Lm + simul$Ym)/infected
  
  simul$logit_p <- logit(simul$p)
  simul$logit_q <- logit(simul$q)
  simul$logit_f <- logit(simul$f)
  simul$logit_g <- logit(simul$g)
  
  return(simul)
}

# Parameters & initial conditions

alpha_w <- 7e-03
alpha_m <- 2e-02
phi_w <- 2e-01
phi_m <- 2e-02
a <- 3e-09
b <- 0.1
B <- 80
r <- 1.4
tau <- 1.5
K <- 1e+09
delta <- 0.8

parms <- c("alpha_w" = alpha_w, "alpha_m" = alpha_m, "phi_w" = phi_w, "phi_m" = phi_m,
           "a" = a, "b" = b, "B" = B, "r" = r, "tau" = tau, "K" = K, "delta" = delta)

prev0_epidemic <- 0.01
prev0_endemic <- 0.99

# Time points

times <- seq(0,60,0.1)

# Simulation

## The density of susceptible cells decreasing over time

simul <- rbind(lambda_simul(prev0=prev0_epidemic, parms=parms, times=times),
               lambda_simul(prev0=prev0_endemic, parms=parms, times=times)) %>% as.data.frame

simul$Q_YL <- simul$f*(1-simul$p)/((1-simul$f)*simul$p)
simul$Q_VY <- simul$q*(1-simul$f)/((1-simul$q)*simul$f)
simul$Q_VL <- simul$q*(1-simul$p)/((1-simul$q)*simul$p)

simul$treatment <- 1:2 %>% rep(each = length(times)) %>% factor(labels = c("epidemic", "endemic"))

mytheme <- theme_classic() + theme(axis.text.x = element_text(size=10),
                                   axis.text.y = element_text(size=10),
                                   plot.title = element_text(size=9, face=3),
                                   legend.title = element_text(size=10))

mytheme2 <- theme_bw() + theme(axis.text.x = element_text(size=10),
                               axis.text.y = element_text(size=10),
                               plot.title = element_text(size=9, face=3),
                               legend.title = element_text(size=10))

eq <- K*(r-(delta+alpha_w))/(delta+alpha_w+tau)*
  c("L"=(delta+tau)/r, "Y"=alpha_w/r, "V"=B*alpha_w*tau/(a*K*(r-(delta+alpha_w))+delta*r))

scientific_10 <- function(x){
  return(parse(text=gsub("1e[+]", "10^", scales::scientific_format()(x))))
}

(Fig_densities <- ggarrange(
  ggdraw() +
    draw_plot(ggplot(data = simul, aes(x=time, y=S, col = treatment)) +
                geom_line(cex = 0.6) +
                labs(y = "\nS(t)", x = "time (hours)", col="Treatment") +
                scale_y_log10(labels=scientific_10) + annotation_logticks(sides = "l") +
                scale_color_manual(values = c('#DA0F0F', '#0E5DD6')) +
                scale_x_continuous(breaks = seq(times[1], times[length(times)], 10)) +
                mytheme + theme(legend.position = 'none', axis.title.x = element_blank())) +
    draw_plot(ggplot(data = simul, aes(x=time, y=1-S/N, col = treatment)) +
                geom_line(cex = 0.6) +
                geom_hline(yintercept = 1, lty='dashed', cex = 0.2) +
                labs(y = "Prevalence", x = "time (hours)", col="Treatment") +
                scale_color_manual(values = c('#DA0F0F', '#0E5DD6')) +
                scale_x_continuous(breaks = seq(times[1], times[length(times)], 10)) +
                theme_classic() + theme(axis.title.x = element_text(size=6),
                                        axis.text.x = element_text(size=6),
                                        axis.title.y = element_text(size=7),
                                        axis.text.y = element_text(size=6),
                                        legend.position = 'none'),
              x=0.47, y=0.5, width=0.5, height=0.5),
  ggplot(data = simul, aes(x=time, y=L, col = treatment)) +
    geom_line(cex = 0.6) +
    geom_hline(yintercept = eq[["L"]], lty='dashed', cex = 0.5) +
    labs(y = "L(t)", x = "time (hours)", col="Treatment") +
    scale_y_log10(labels=scientific_10) + annotation_logticks(sides = "l") +
    scale_color_manual(values = c('#DA0F0F', '#0E5DD6')) +
    scale_x_continuous(breaks = seq(times[1], times[length(times)], 10)) +
    mytheme + theme(legend.position = 'none', axis.title.x = element_blank()),
  ggplot(data = simul, aes(x=time, y=Y, col = treatment)) +
    geom_line(cex = 0.6) +
    geom_hline(yintercept = eq[["Y"]], lty='dashed', cex = 0.5) +
    labs(y = "\nY(t)", x = "time (hours)", col="Treatment") +
    scale_y_log10(labels=scientific_10) + annotation_logticks(sides = "l") +
    scale_color_manual(values = c('#DA0F0F', '#0E5DD6')) +
    scale_x_continuous(breaks = seq(times[1], times[length(times)], 10)) +
    mytheme + theme(legend.position = 'none'),
  ggplot(data = simul, aes(x=time, y=V, col = treatment)) +
    geom_line(cex = 0.6) +
    geom_hline(yintercept = eq[["V"]], lty='dashed', cex = 0.5) +
    labs(y = "V(t)", x = "time (hours)", col="Treatment") +
    scale_y_log10(labels=scientific_10) + annotation_logticks(sides = "l") +
    scale_color_manual(values = c('#DA0F0F', '#0E5DD6')) +
    scale_x_continuous(breaks = seq(times[1], times[length(times)], 10)) +
    mytheme  + theme(legend.position = 'none'))) # 7 x 5

(Fig_freq <- ggarrange(
  ggplot(data = simul, aes(x=time, y=logit_p, color = treatment)) +
    geom_line(cex = 0.6) +
    scale_color_manual(values = c('#DA0F0F', '#0E5DD6')) +
    labs(title = "among L cells", x = "time (hours)", y = "logit(p(t))", col="Treatment") +
    scale_x_continuous(breaks = seq(times[1], times[length(times)], 10)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(-1, 2)) +
    mytheme2,
  ggplot(data = simul, aes(x=time, y=logit_g, color = treatment)) +
    geom_line(cex = 0.6) +
    scale_color_manual(values = c('#DA0F0F', '#0E5DD6')) +
    labs(title="among infected cells (L or Y)", x="time (hours)", y="logit(g(t))", col="Treatment") +
    scale_x_continuous(breaks = seq(times[1], times[length(times)], 10)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(-1, 2)) +
    mytheme2,
  ggplot(data = simul, aes(x=time, y=logit_f, color = treatment)) +
    geom_line(cex = 0.6) +
    scale_color_manual(values = c('#DA0F0F', '#0E5DD6')) +
    labs(title="among Y cells", x="time (hours)", y="logit(f(t))", col="Treatment") +
    scale_x_continuous(breaks = seq(times[1], times[length(times)], 10)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(0, 3.25)) +
    mytheme2,
  ggplot(data = simul, aes(x=time, y=logit_q, color = treatment)) +
    geom_line(cex = 0.6) +
    scale_color_manual(values = c('#DA0F0F', '#0E5DD6')) +
    labs(title = "among free viral particles (V)", x = "time (hours)", y="logit(q(t))", col="Treatment") +
    scale_x_continuous(breaks = seq(times[1], times[length(times)], 10)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(0, 3.25)) +
    mytheme2,
  nrow = 2, ncol=2, common.legend = TRUE, align = 'v')) # 7 x 5

(Fig_freq.2 <- simul[,which(colnames(simul) %in%
                              c("time","logit_p","logit_q","logit_f","logit_g","treatment"))] %>%
    tidyr::gather("label", "logit_freq", -c(time, treatment)) %>%
    mutate(label = factor(label, levels = c("logit_g","logit_p","logit_f","logit_q"))) %>%
    ggplot(aes(x=time, y=logit_freq, lty=label, col=treatment, alpha=label)) +
    geom_line(cex = 0.6) +
    scale_linetype_manual(values = c(1,4,3,2),
                          labels = c("among infected cells (L or Y)",
                                     "among L cells",
                                     "among Y cells",
                                     "in the free virus stage (V)")) +
    scale_alpha_manual(values=c(0.8,1,1,0.8)) +
    scale_color_manual(values = c('#DA0F0F', '#0E5DD6')) +
    labs(x='time (hours)', y='logit-frequency of virulent phage\n', col="Treatment",
         lty = "Logit-frequency of virulent phage") +
    scale_x_continuous(breaks = seq(times[1], times[length(times)], 10)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    mytheme2 + theme(legend.text = element_text(hjust=0), legend.title=element_text(face='bold')) +
    guides(alpha=FALSE)) # 4 x 7

plot_grid(Fig_densities, Fig_freq.2, rel_heights = c(0.55,0.45),
          labels = paste0(LETTERS[1:2],")"), ncol = 1) # portrait : 8 x 7

logit_g.30 <- subset(simul, treatment=='epidemic' & time==30)$logit_g
logit_g.50 <- subset(simul, treatment=='epidemic' & time==50)$logit_g

predict.0 <- ddply(simul, ~treatment, function(X){
  i <- which(X$prev>=0.955)[1]
  return(c("time"=X$time[i], "logit_g"=X$logit_g[i], "treatment"=X$treatment[i]))
})

predict <- rbind(predict.0,
                 data.frame("time" = 60,
                            "logit_g" = predict.0$logit_g +
                              (alpha_w-alpha_m)*(60-predict.0$time),
                            "treatment" = predict.0$treatment))

predict$treatment <- factor(predict$treatment, labels=c("epidemic", "endemic"))

(Fig_predict_logit_g <- ggplot(simul, aes(x=time, y=logit_g, col=treatment)) +
    geom_line(cex=1.75, alpha=0.2) +
    geom_line(cex=0.7, data=predict) +
    geom_point(data=predict.0, pch=21, col='black', fill='white', size=2.25) +
    # geom_segment(x=30, xend=50, y=logit_g.30, yend=logit_g.30, col='black', lty='dashed', lwd=0.3) +
    # geom_segment(x=50, xend=50, y=logit_g.30, yend=logit_g.50, col='black', lty='dashed', lwd=0.3) +
    scale_color_manual(values = c('#DA0F0F', '#0E5DD6')) +
    labs(x='\ntime (hours)', y="Logit-frequency infected by virulent phage", col="Treatment") +
    scale_x_continuous(breaks = seq(times[1], times[length(times)], 10)) +
    scale_y_continuous(breaks = seq(-1, 2, 0.5)) +
    # annotate(geom='text', label=expression(paste("-", Delta, alpha)),
    #          size = 4.5, x=43, y=logit_g.30+0.3, hjust=0) +
    mytheme)

(Fig_QVL <- ggplot(simul, aes(x=time, y=Q_VL, col=treatment)) +
    geom_line(cex=1.5, alpha=0.35) +
    labs(x='\ntime (hours)', y="Differentiation between virulent\nfree phages and prophages", col="Treatment") +
    # geom_hline(yintercept = alpha_m/alpha_w, lty = 'dashed', lwd=0.3) +
    geom_segment(aes(x=predict.0$time[2], xend=60, y=alpha_m/alpha_w, yend=alpha_m/alpha_w),
                 cex=0.5, col='#0E5DD6') +
    geom_segment(aes(x=predict.0$time[1], xend=60, y=alpha_m/alpha_w, yend=alpha_m/alpha_w),
                 cex=0.5, col='#DA0F0F') +
    scale_color_manual(values = c('#DA0F0F', '#0E5DD6')) +
    scale_x_continuous(breaks = seq(times[1], times[length(times)], 10)) +
    scale_y_continuous(trans='log2', limits=c(2,17)) +
    annotation_logticks(sides = "l") +
    # annotate(geom='text', label=expression(paste("1 + ", frac(paste(Delta,alpha), alpha[w]))),
    #          size = 4.5, x=50, y=1.5*alpha_m/alpha_w, hjust=0) +
    mytheme)

simul[,grep("time|treatment|Q_", colnames(simul))] %>%
  tidyr::gather("id", "Q", -c(time, treatment)) %>%
  mutate(id=factor(id, labels=c("between free viruses (V) and prophages (L)",
                                "between free viruses (V) and Y cells",
                                "between Y cells and prophages (L)"))) %>%
  ggplot(aes(x=time, y=Q, col=treatment)) +
  facet_wrap(~id, ncol=1) +
  geom_line(cex = 0.6) +
  geom_hline(yintercept=c(alpha_m/alpha_w, 1), lty='dashed', cex=0.3) +
  labs(x="time (hours)", y="Differentiation of virulent phage", col="") +
  scale_y_continuous(trans='log2') +
  annotation_logticks(sides = "l") +
  scale_x_continuous(breaks = seq(times[1], times[length(times)], 10)) +
  scale_color_manual(values = c('#DA0F0F', '#0E5DD6')) +
  mytheme2 + theme(legend.position='top') # pdf portrait: 6 x 4

## The density of susceptible cells remaining constant over time

times <- seq(0,9,0.1)

simul2 <- rbind(lambda_simul(prev0=prev0_epidemic, parms=parms, times=times, S_cst=FALSE),
                lambda_simul(prev0=prev0_epidemic, parms=parms, times=times, S_cst=TRUE)) %>% as.data.frame

simul2$model <- 1:2 %>% rep(each = length(times)) %>%
  factor(labels = c("original", "approx"))

ggarrange(
  ggplot(data = simul2, aes(x=time, y=logit_p, lty = model, color = model)) +
    geom_line(cex = 0.6) +
    scale_color_manual(values = c('black', 'red')) +
    labs(title = "among L cells", x = "time (hours)", y = "p(t) (logit scale)\n") +
    scale_x_continuous(breaks = times[1]:times[length(times)]) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    mytheme,
  ggplot(data = simul2, aes(x=time, y=logit_g, color = model, lty = model)) +
    geom_line(cex = 0.6) +
    scale_color_manual(values = c('black', 'red')) +
    labs(title = "among infected cells (L or Y)", x = "time (hours)", y = expression(f^YL*(t)~'(logit'~'scale)')) +
    scale_x_continuous(breaks = times[1]:times[length(times)]) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    mytheme,
  ggplot(data = simul2, aes(x=time, y=logit_f, color = model, lty = model)) +
    geom_line(cex = 0.6) +
    scale_color_manual(values = c('black', 'red')) +
    labs(title = "among Y cells", x = "time (hours)", y = "f(t) (logit scale)\n") +
    scale_x_continuous(breaks = times[1]:times[length(times)]) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    mytheme,
  ggplot(data = simul2, aes(x=time, y=logit_q, color = model, lty = model)) +
    geom_line(cex = 0.6) +
    scale_color_manual(values = c('black', 'red')) +
    labs(title = "among free viral particles (V)", x = "time (hours)", y = "q(t) (logit scale)\n") +
    scale_x_continuous(breaks = times[1]:times[length(times)]) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    mytheme,
  nrow = 2, ncol=2, common.legend = TRUE, align = 'v') # 7 x 5

S0 <- (1-prev0_epidemic)*parms[["K"]]

approx_slope_epidemic <- function(tau, a, b, B, K, phi_w, phi_m, S0){
  Delta_phi <- phi_m-phi_w; tauB <- tau*B
  Z <- (tau-a*K+sqrt((tau-a*K)^2+4*(1-phi_w)*a*b*S0*tau*B))/(2*(1-phi_w)*a*b*S0)
  abS0Z <- a*b*S0*Z
  return(-Delta_phi*abS0Z*((tauB+abS0Z*Z*c(0.5,0)*Delta_phi)/(tauB+abS0Z*Z*(1-phi_w)))) # f=0.5 and f=1
}

(Fig_ratio <- ggarrange(
  ggplot(simul2, aes(x= time, y=L/Y, lty = model)) +
    geom_hline(yintercept = 0, cex = 0.3, col = 'black') +
    geom_line(cex = 0.6, col='#DA0F0F') +
    labs(x = "time (hours)", y="L(t) / Y(t)") +
    scale_linetype_discrete(labels=c("variable S", "constant S")) +
    scale_x_continuous(breaks=0:9) +
    mytheme + theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=9),
                    axis.title.y = element_text(size=10), axis.text.y = element_text(size=9),
                    legend.title = element_blank()),
  ggplot(simul2, aes(x= time, y=V/Y, lty = model)) +
    geom_hline(yintercept = (tau-a*K+sqrt((tau-a*K)^2+4*(1-phi_w)*a*b*S0*tau*B))/(2*(1-phi_w)*a*b*S0),
               cex = 0.3, col = 'black') +
    geom_line(cex = 0.6, col='#DA0F0F') +
    labs(x = "time (hours)", y="V(t) / Y(t)", lty="Availability of susceptible hosts over time") +
    scale_linetype_discrete(labels=c("decreasing", "constant")) +
    scale_x_continuous(breaks=0:9) +
    mytheme + theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=9),
                    axis.title.y = element_text(size=10), axis.text.y = element_text(size=9),
                    legend.title = element_blank()),
  common.legend = TRUE)) # 7 x 3

t0 <- 1
approx_logit_f <- data.frame("time" = times[times>=t0])

approx_slope_logit_f <- approx_slope_epidemic(tau, a, b, B, K, phi_w, phi_m, S0)

t0_index <- which(simul2$time == t0)[2]
approx_logit_f$logit_f.1 <- simul2$logit_f[t0_index]+approx_slope_logit_f[1]*(approx_logit_f$time-t0)
approx_logit_f$logit_f.2 <- simul2$logit_f[t0_index]+approx_slope_logit_f[2]*(approx_logit_f$time-t0)
approx_logit_f$logit_q.1 <- simul2$logit_q[t0_index]+approx_slope_logit_f[1]*(approx_logit_f$time-t0)
approx_logit_f$logit_q.2 <- simul2$logit_q[t0_index]+approx_slope_logit_f[2]*(approx_logit_f$time-t0)

(Fig_approx_slope <- ggarrange(
  ggdraw() + draw_label(label="Logit-frequency of virulent phage", angle=90, size=10),
  ggplot() +
    geom_ribbon(data=approx_logit_f, aes(x=time, ymin=logit_f.1, ymax=logit_f.2),
                fill = '#DA0F0F', col = NA, alpha = 0.1) +
    geom_line(data = simul2, aes(x=time, y=logit_f, lty = model), col='#DA0F0F', cex=0.6) +
    geom_point(data=data.frame("time"=approx_logit_f$time[1], "logit_q"=approx_logit_f$logit_f.1[1]),
               aes(x=time, y=logit_q), pch=21, size=2, col='black', fill='white') +
    labs(x = "time (hours)") +
    scale_x_continuous(breaks=times[1]:times[length(times)]) +
    scale_y_continuous(breaks=seq(1,6,1), limits=c(1,6.5),
                       labels=scales::label_number(accuracy = 0.1)) +
    annotate(geom='text', label='Among Y cells',
             x=0.25, y=6, hjust=0, size=3, fontface=3) +
    mytheme + theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=9),
                    axis.title.y = element_blank(), legend.position = 'none', plot.margin=margin(0,0,0,0)),
  ggplot() +
    geom_ribbon(data=approx_logit_f, aes(x=time, ymin=logit_q.1, ymax=logit_q.2),
                fill = '#DA0F0F', col = NA, alpha = 0.1) +
    geom_line(data = simul2, aes(x=time, y=logit_q, lty = model), color='#DA0F0F', cex=0.6) +
    geom_point(data=data.frame("time"=approx_logit_f$time[1], "logit_q"=approx_logit_f$logit_q.1[1]),
               aes(x=time, y=logit_q), pch=21, size=2, col='black', fill='white') +
    labs(x = "time (hours)") +
    scale_x_continuous(breaks=times[1]:times[length(times)]) +
    scale_y_continuous(breaks=seq(1,6,1), limits=c(1,6.5), labels=scales::label_number(accuracy = 0.1)) +
    annotate(geom='text', label='In the free virus stage (V)',
             x=0.25, y=6, hjust=0, size=3, fontface=3) +
    mytheme + theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=9),
                    axis.title.y = element_blank(), axis.text.y = element_blank(), 
                    legend.position = 'none', plot.margin=margin(0,0,0,0)),
  ncol=3, widths = c(0.1,0.45,0.45))) # 7 x 4

(Fig_approx_vs_Delta_phi <- seq(0,phi_w,0.05) %>% lapply(function(i){
  prev0 <- prev0_epidemic
  slope <- prev0 %>% lapply(function(j){
    return(approx_slope_epidemic(tau=tau, a=a, b=b, B=B, K=K, phi_w=phi_w, phi_m=i, S0=(1-j)*K))
  })
  return(data.frame("Delta_phi"=i-phi_w, "prev0"=prev0,
                    "lower"=sapply(1:length(slope), function(x){return(slope[[x]][1])}),
                    "upper"=sapply(1:length(slope), function(x){return(slope[[x]][2])})))
}) %>% bind_rows %>%
  ggplot(aes(x=Delta_phi, ymin=lower, ymax=upper, fill="Approximation S constant")) +
    geom_hline(yintercept = 0, lwd = 0.5) +
    geom_vline(xintercept = 0, lwd = 1, color = 'grey') +
    geom_ribbon(col=NA, linewidth= 1E-9, alpha=0.3) +
    labs(x=expression(paste(Delta, phi, " = ", phi[m]-phi[w])),
         y="Approximation of the selection gradient\nof the mutant strain at the onset of the epidemic") +
    scale_x_continuous(expand=c(0,0), limits=c(-phi_w, 5E-3)) +
    scale_y_continuous(expand=c(0,0), breaks=seq(0,0.6,0.1)) +
    theme_bw() +
    theme(axis.title.x=element_text(size=12, margin=margin(10,0,0,0)),
          axis.text.x=element_text(size=9),
          axis.title.y=element_text(size=10, margin=margin(0,10,0,0)),
          axis.text.y=element_text(size=9),
          legend.title=element_blank(), legend.position = 'top')) # pdf landscape 5 x 4

plot_grid(Fig_approx_vs_Delta_phi,
          plot_grid(Fig_ratio, Fig_approx_slope + theme(legend.position = 'none'),
                    align='v', ncol=1, rel_heights=c(0.45,0.55), labels=paste0(LETTERS[2:3],")")), # pdf 6 x 6
  ncol = 2, rel_widths=c(0.45,0.55), labels=paste0(LETTERS[1],")")) # pdf landscape 10 x 5

approx_logit_f$logit_q <- simul2$logit_q[which(simul2$time == t0)[1]]+approx_slope_logit_f[2]*(approx_logit_f$time-t0)

(Fig_predict_logit_q <- ggplot() +
    geom_ribbon(data=approx_logit_f, aes(x=time, ymin=logit_q.1, ymax=logit_q.2,
                                         fill="approximation\nS constant "),
                col = NA, alpha = 0.1) +
    geom_line(data = simul2, aes(x=time, y=logit_q, lty = model), color='#DA0F0F', cex=0.7) +
    geom_point(data=data.frame("time"=approx_logit_f$time[1], "logit_q"=approx_logit_f$logit_q.1[1]),
               aes(x=time, y=logit_q), pch=21, size=2.25, col='black', fill='white') +
    labs(x = "\ntime (hours)", y = "Logit-frequency free virulent phage") +
    scale_linetype_discrete(labels=c("variable S", "constant S")) +
    scale_fill_manual(values='#DA0F0F') +
    scale_x_continuous(breaks=times[1]:times[length(times)]) +
    scale_y_continuous(breaks=seq(1,6,1), limits=c(1,6.2)) +
    # geom_segment(data=approx_logit_f, aes(x=time[25], xend=time[50], y=logit_q.2[50], yend=logit_q.2[50]),
    #              col='black', lty='dashed', lwd=0.3) +
    # geom_segment(data=approx_logit_f, aes(x=time[25], xend=time[25], y=logit_q.2[50], yend=logit_q.2[25]),
    #              col='black', lty='dashed', lwd=0.3) +
    # annotate(geom='text', label=expression(paste(" "%prop%-Delta, phi)), size=5, x=3.1, y=4.7, hjust=0) +
    mytheme + theme(legend.title = element_blank()) +
    guides(fill=FALSE))

plot_grid(
  plot_grid(Fig_predict_logit_q + theme(axis.title.x = element_blank(), legend.position = 'none'),
            Fig_predict_logit_g + theme(axis.title.x = element_blank(), legend.position = 'none'),
            Fig_QVL + theme(legend.position = 'none'),
            align='v', ncol=1, labels=paste0(LETTERS[1:3],")"),
            label_x=c(0.2, 0.85,0.85), label_y=0.95),
  plot_grid(get_legend(Fig_predict_logit_q), get_legend(Fig_predict_logit_g), NULL, ncol=1),
  rel_widths = c(0.8,0.2)) # 8.5 x 6

# ggsave("Theoretical_predictions.tiff", width=6, height=8.5, dpi=800)

# Lysis time distributions

## Model without Y stage

lambda_ODE <- function(t, y, parms){ 
  # t = current time
  # y = current values of the state variables (at time t)
  # parms = set of parameters
  
  # (1) State variables
  
  # Two viral strains: subscript w -> Resident strain (WT)
  #                    subscript m -> Mutant strain
  
  S <- y[["S"]] # S(t), density of Susceptible bacteria
  
  # L(t), density of bacteria in a lysogenic state
  L_w <- y[["L_w"]]
  L_m <- y[["L_m"]]
  L <- L_w + L_m
  
  # V(t), density of Virions (free viral particles)
  V_w <- y[["V_w"]]
  V_m <- y[["V_m"]]
  
  N <- S + L # Total density of bacteria
  
  # (2) Results
  
  # New infections
  abVRS <- parms[["a"]]*parms[["b"]]*V_w*S
  abVMS <- parms[["a"]]*parms[["b"]]*V_m*S
  
  # New reactivations of integrated phages / prophages (i.e., switches to a lytic cycle)
  reactivations_w <- parms[["alpha_w"]]*L_w
  reactivations_m <- parms[["alpha_m"]]*L_m
  
  # Growth and removal
  growth <- parms[["r"]]*(1-N/parms[["K"]]) # Rate of logistic growth
  delta <- parms[["delta"]] # Removal rate of the chemostat
  removal.V <- parms[["a"]]*N + delta # Removal rate of virions
  # The term a*N represents the rate at which virions adsorbed on bacteria S, L and Y (adsorption in non reversible).
  
  # Temporal derivatives
  
  dS <- growth*S - (abVRS + abVMS) - delta*S
  
  dL_w <- growth*L_w + parms[["phi_w"]]*abVRS - reactivations_w - delta*L_w
  dL_m <- growth*L_m + parms[["phi_m"]]*abVMS - reactivations_m - delta*L_m
  
  dV_w <- (reactivations_w + (1-parms[["phi_w"]])*abVRS)*parms[["B"]] - removal.V*V_w
  dV_m <- (reactivations_m + (1-parms[["phi_m"]])*abVMS)*parms[["B"]] - removal.V*V_m
  
  return(list(c(dS, dL_w, dL_m, dV_w, dV_m)))
}

# Model with n > 0 Y stage(s)

lambda_ODE.n <- function(t, y, parms){ 
  # t = current time
  # y = current values of the state variables (at time t)
  # parms = set of parameters
  
  # (1) State variables
  
  # Two viral strains: subscript w -> Resident strain (WT)
  #                    subscript m -> Mutant strain
  
  S <- y[["S"]] # S(t), density of Susceptible bacteria
  
  # L(t), density of bacteria in a lysogenic state
  L_w <- y[["L_w"]]
  L_m <- y[["L_m"]]
  L <- L_w + L_m
  
  # Y(t), densities of infected bacteria prior to lysis
  Y_w <- y[grep("Y_w", names(y))]
  Y_m <- y[grep("Y_m", names(y))]
  Y <- sum(Y_w, Y_m)
  
  n <- length(Y_w) # Number of stages within compartment Y
  
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
  lysis_w <- n*parms[["tau"]]*Y_w # vector of length n
  lysis_m <- n*parms[["tau"]]*Y_m # vector of length n
  
  # Growth and removal
  growth <- parms[["r"]]*(1-N/parms[["K"]]) # Rate of logistic growth
  delta <- parms[["delta"]] # Removal rate of the chemostat
  removal.V <- parms[["a"]]*N + delta # Removal rate of virions
  # The term a*N represents the rate at which virions adsorbed on bacteria S, L and Y (adsorption in non reversible).
  
  # Temporal derivatives
  
  dS <- growth*S - (abVRS + abVMS) - delta*S
  
  dL_w <- growth*L_w + parms[["phi_w"]]*abVRS - reactivations_w - delta*L_w
  dL_m <- growth*L_m + parms[["phi_m"]]*abVMS - reactivations_m - delta*L_m
  
  dY_w <- ((1-parms[["phi_w"]])*abVRS + reactivations_w)*c(1, rep(0, n-1)) + c(0, lysis_w[-n]) - lysis_w - delta*Y_w
  dY_m <- ((1-parms[["phi_m"]])*abVMS + reactivations_m)*c(1, rep(0, n-1)) + c(0, lysis_m[-n]) - lysis_m - delta*Y_m
  
  dV_w <- lysis_w[n]*parms[["B"]] - removal.V*V_w
  dV_m <- lysis_m[n]*parms[["B"]] - removal.V*V_m
  
  return(list(c(dS, dL_w, dL_m, dY_w, dY_m, dV_w, dV_m)))
}

lambda_simul <- function(prev0, freq0=0.5, parms, n=1, times, solver="lsoda", verbose=TRUE){
  # prev0 = initial value for prevalence
  # freq0 = initial value for frequency of the virulent phage within lysogenized bacteria
  #         (default: 0.5, i.e. initial ratio 1:1)
  # parms = vector of model parameters (with names)
  # Number of stages within compartment Y
  
  # times = vector of desired time points for the simulation
  # solver = argument 'method' of function 'ode', i.e. the integrator to use
  # verbose = logical indicating whether the messages of the function should be printed (defaut: TRUE)
  
  if(verbose){
    "=" %>% rep(66) %>% c("\n") %>% paste(collapse = '') %>% cat
    cat("  MODEL - No plasticity on any phenotypic traits of phage lambda\n")
    "=" %>% rep(66) %>% c("\n") %>% paste(collapse = '') %>% cat
    " " %>% rep(41) %>% c("Initial prevalence = ", round(100*prev0,2), " %\n", rep(" ", 40),
                          "(with ratio R:M = ", round(1-freq0, 2), ":", round(freq0, 2), ")\n\n",
                          "Number of stage(s) within compartment Y = ", n, "\n") %>%
      paste(collapse = '') %>% cat
  }
  
  # (2) Initial conditions
  
  init <- (c("S" = (1-prev0), "L_w" = (1-freq0)*prev0, "L_m" = freq0*prev0)*parms[["K"]]) %>%
    c("Y_w" = rep(0, n), "Y_m" = rep(0, n), "V_w" = 0, "V_m" = 0)
  
  # (3) Simulation
  
  if(n==0){
    simul <- ode(y=init, times=times, parms=parms, func=lambda_ODE, method=solver) %>% as.data.frame
  }else{
    simul <- ode(y=init, times=times, parms=parms, func=lambda_ODE.n, method=solver) %>% as.data.frame
  }
  if(verbose){cat("\nSimulation performed using a system of ordinary differential equations\n\n")}
  
  # Densities and prevalence
  
  if(n > 1){
    simul$Y_w <- simul[,grep("Y_w", names(simul))] %>% apply(1,sum)
    simul$Y_m <- simul[,grep("Y_m", names(simul))] %>% apply(1,sum)
  }
  
  if(n==0){
    
    simul$L <- simul$L_w + simul$L_m
    simul$Total_density <- simul$S + simul$L
    simul$Prevalence <- simul$L/simul$Total_density
    
    # Frequencies of virulent phage
    simul$p <- simul$L_m/simul$L # p(t), among lysogenic bacteria (L)
    simul$q <- simul$V_m/(simul$V_w + simul$V_m) # q(t), in free virus stage (V)
    simul$g <- simul$p # g(t), among infected bacteria (same as p(t))
  }else{
    simul$Y <- simul$Y_w + simul$Y_m
    simul$L <- simul$L_w + simul$L_m
    simul$Infected <- simul$Y + simul$L
    
    simul$Total_density <- simul$S + simul$Infected
    simul$Prevalence <- simul$Infected/simul$Total_density
    
    # Frequencies of virulent phage
    simul$p <- simul$L_m/simul$L # p(t), among lysogenic bacteria (L)
    simul$q <- simul$V_m/(simul$V_w + simul$V_m) # q(t), in free virus stage (V)
    simul$f <- simul$Y_m/simul$Y # f(t), among Y bacteria
    simul$g <- (simul$Y_m + simul$L_m)/simul$Infected # g(t), among infected bacteria (both Y and L)
  }
  return(simul)
}

n_vec <- 0:10

lysis_time <- seq(0, 1.5, 0.025)

(Fig_PDF <- n_vec[-1] %>% lapply(function(n){
  return(data.frame("n"=n,
                    "lysis_time"=lysis_time,
                    "PDF"=dgamma(lysis_time, shape = n, scale = 1/(n*tau))))
}) %>% bind_rows %>%
    ggplot(aes(x=lysis_time, y=PDF, col=as.factor(n))) +
    geom_vline(xintercept = 1/tau, lty = 'dashed') +
    geom_line(cex=0.8) +
    labs(x="\nLysis time (hours)", y="Probability density\n", col="n, number of stages\nwithin compartment Y") +
    annotate(geom='text', label=expression("Lysis time ~"~Gamma~bgroup("(","n,"~frac(1,n~tau),")")),
             size=3.5, x=1, y=1.7, hjust=0) +
    theme_bw() +
    theme(axis.title.x = element_text(size=11), axis.text.x = element_text(size=11),
          axis.title.y = element_text(size=11), axis.text.y = element_text(size=11),
          legend.position = 'top'))

times <- seq(0,60,0.1)

simul <- n_vec %>% lapply(function(n){
  return(rbind(lambda_simul(prev0=prev0_epidemic, parms=parms, n=n, times=times, verbose=FALSE) %>%
                 cbind("treatment"='epidemic'),
               lambda_simul(prev0=prev0_endemic, parms=parms, n=n, times=times, verbose=FALSE) %>%
                 cbind("treatment"='endemic')) %>%
           cbind("n"=n))
}) %>% bind_rows

simul$treatment <- factor(simul$treatment, levels = c('epidemic', 'endemic'))

Fig_simul.prev <- ggplot(simul, aes(x=time, y=logit(Prevalence), col=interaction(n, treatment))) +
  geom_line(cex=0.8) +
  scale_color_manual("n, number of stages\nwithin compartment Y",
                     labels = c(rep("", length(n_vec)), paste("  ", n_vec)),
                     values = c(colorRampPalette(brewer.pal(9, "YlOrRd"))(length(n_vec)+1)[-1] %>% rev,
                                colorRampPalette(brewer.pal(9, "BuPu"))(length(n_vec)+1)[-1] %>% rev)) +
  labs(x="time (hours)", y="Logit-prevalence\n") +
  scale_x_continuous(breaks=seq(0,60,10)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size=11),
        axis.title.y = element_text(size=10), axis.text.y = element_text(size=11)) +
  guides(col=guide_legend(ncol=2))

Fig_simul.g <- ggplot(simul, aes(x=time, y=logit(g), col=interaction(n, treatment))) +
  geom_line(cex=0.8) +
  scale_color_manual("n, number of stages\nwithin compartment Y\n(epidemic)",
                     labels = c(rep("", length(n_vec)), paste("  ", n_vec)),
                     values = c(colorRampPalette(brewer.pal(9, "YlOrRd"))(length(n_vec)+1)[-1] %>% rev,
                                colorRampPalette(brewer.pal(9, "BuPu"))(length(n_vec)+1)[-1]%>% rev)) +
  labs(x="time (hours)", y="Logit-frequency infected by virulent phage\n") +
  scale_x_continuous(breaks=seq(0,60,10)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size=11),
        axis.title.y = element_text(size=10), axis.text.y = element_text(size=11),
        legend.position = 'none')

Fig_simul.q <- ggplot(simul, aes(x=time, y=logit(q), col=interaction(n, treatment))) +
  geom_line(cex=0.8) +
  scale_color_manual("n, number of stages\nwithin compartment Y",
                     labels = c(rep("", length(n_vec)), paste("  ", n_vec)),
                     values = c(colorRampPalette(brewer.pal(9, "YlOrRd"))(length(n_vec)+1)[-1] %>% rev,
                                colorRampPalette(brewer.pal(9, "BuPu"))(length(n_vec)+1)[-1] %>% rev)) +
  labs(x="time (hours)", y="Logit-frequency free virulent phage\n") +
  scale_x_continuous(breaks=seq(0,60,10)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=11), axis.text.x = element_text(size=11),
        axis.title.y = element_text(size=10), axis.text.y = element_text(size=11),
        legend.position = 'none')

plot_grid(
  ggarrange(Fig_simul.prev + theme(legend.position = 'none'),
            Fig_simul.g, Fig_simul.q, ncol=1, align='v'),
  plot_grid(NULL, Fig_PDF, plot_grid(get_legend(Fig_simul.prev), NULL, rel_widths=c(0.4,0.6)),
            rel_heights=c(0.05,0.5,0.45), ncol=1),
  ncol=2, labels=paste0(LETTERS[1:2], ")"))

# pdf: landscape 10x8

# ggplot(subset(simul, treatment=="epidemic" & n==0), aes(x=time, y=logit(g), col=interaction(n, treatment))) +
#   geom_line(cex=0.8) +
#   scale_color_manual("n, number of stages\nwithin compartment Y\n(epidemic)",
#                      labels = c(rep("", length(n_vec)), paste("  ", n_vec)),
#                      values = c(colorRampPalette(brewer.pal(9, "YlOrRd"))(length(n_vec)+1)[-1] %>% rev,
#                                 colorRampPalette(brewer.pal(9, "BuPu"))(length(n_vec)+1)[-1]%>% rev)) +
#   labs(x="time (hours)", y="Logit-frequency infected by virulent phage\n") +
#   scale_x_continuous(breaks=seq(0,60,10)) +
#   scale_y_continuous(limits=c(-0.7,2)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(size=11),
#         axis.title.y = element_text(size=10), axis.text.y = element_text(size=11),
#         legend.position = 'none')
# 
# ggsave("logit_g_lysis_time_instantaneous.png", width=5, height=3.6, dpi=600)