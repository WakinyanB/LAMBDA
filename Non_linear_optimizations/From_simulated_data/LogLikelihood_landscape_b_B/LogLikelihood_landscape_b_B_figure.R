rm(list=ls())

library(tidyverse)
library(scales)
library(RColorBrewer)
library(viridis)
library(plotly)

setwd("C:/Users/wbenhamou/Desktop/LAMBDA/Data_and_codes/")

a <- 3e-09
b <- 0.1
B <- 80

tab <- read.csv("Outputs_csv/Lambda_landscape_b_B_v11_n40401_maxfeval8000.csv", header = TRUE)
tab2 <- tab
tab2$Value[!tab2$Convergence == "Successful completion"] <- NA

ggplot(subset(tab2, B !=0 & b != 0), aes(x=b, y=B)) +
  geom_tile(aes(fill=-Value)) +
  geom_hline(yintercept = B, lty='dashed', col='white', lwd=0.4) +
  geom_vline(xintercept = b, lty='dashed', col='white', lwd=0.4) +
  scale_fill_viridis_c(option='A') +
  scale_x_continuous(expand=c(0,0), breaks=seq(0,0.2,0.025)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(0,100,10)) +
  labs(fill="Log-likelihood") +
  theme_bw() +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),
        axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 10, hjust = 0),
        legend.text = element_text(size=8)) +
  geom_abline(slope=)

# ggsave("landscape_b_B.png", width = 7, height = 5, dpi = 600)

tab3 <- subset(tab2, B !=0 & b != 0 & Value < -45000)
summary(tab3$b*tab3$B)