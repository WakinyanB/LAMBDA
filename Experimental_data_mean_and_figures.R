rm(list=ls())

setwd("C:/Users/wbenhamou/Desktop/LAMBDA/Data_and_codes/")

data_FACS <- read.csv("Data/data_FACS.csv", header = TRUE)
data_qPCR <- read.csv("Data/data_qPCR.csv", header = TRUE)

data_FACS_mean <- data_FACS %>% plyr::ddply(~prev0+time,
                                            function(subdata){subdata[,-c(1,8)] %>%
                                                apply(2, mean)})
data_FACS_mean <- data.frame("chemostat" = plyr::ddply(data_FACS, ~prev0+time, function(subdata){
  "chem" %>% paste0(gsub("chem", "-", unique(subdata$chemostat)) %>% paste(collapse = ""))})[,3],
  data_FACS_mean)

data_qPCR_mean <- data_qPCR %>% plyr::ddply(~prev0+time,
                                            function(subdata){subdata[,-c(1,9)] %>%
                                                apply(2, mean, na.rm = TRUE)})
data_qPCR_mean <- data.frame("chemostat"=plyr::ddply(data_qPCR, ~prev0+time, function(subdata){
  "chem" %>% paste0(gsub("chem", "-", unique(subdata$chemostat)) %>% paste(collapse = ""))})[,3],
  data_qPCR_mean)

mytheme <- theme_bw()+theme(axis.title.x=element_text(size=10), axis.title.y=element_text(size=10),
                            axis.text.x=element_text(size=9), axis.text.y=element_text(size=9),
                            legend.title = element_text(size=9))

# Plots

Fig_logit_prev <- ggplot(data = data_FACS[!is.na(data_FACS$logit_prev),],
                         aes(x = time, y = logit_prev, group = chemostat, color = as.factor(prev0))) +
  geom_line(linewidth = 0.3) +
  geom_line(data = data_FACS_mean[!is.na(data_FACS_mean$logit_prev),], aes(group = prev0), linewidth = 1.3, alpha = 0.6)+
  scale_color_manual(labels = c("epidemic", "endemic"), values = c('#DA0F0F', '#0E5DD6')) +
  labs(x = "time (hours)", y = "Logit-prevalence\n", col="Treatment:") +
  scale_x_continuous(limits = c(0,60), breaks = seq(0, 60, 10)) +
  scale_y_continuous(limits = c(-9,12), labels = scales::label_number(accuracy = 1)) +
  mytheme

Fig_logit_g <- ggplot(data = data_FACS[!is.na(data_FACS$logit_g),],
                      aes(x = time, y = logit_g, group = chemostat, color = as.factor(prev0))) +
  geom_line(linewidth = 0.3) +
  geom_line(data = data_FACS_mean[!is.na(data_FACS_mean$logit_g),], aes(group = prev0), linewidth = 1.3, alpha = 0.6)+
  scale_color_manual(labels = c("epidemic", "endemic"), values = c('#DA0F0F', '#0E5DD6')) +
  labs(x = "time (hours)", y = "Logit-frequency infected by virulent phage\n", col="Treatment:") +
  scale_x_continuous(limits = c(0,60), breaks = seq(0, 60, 10)) +
  scale_y_continuous(limits = c(-2,6), labels = scales::label_number(accuracy = 1)) +
  mytheme

Fig_logit_q <- ggplot(data = data_qPCR[!is.na(data_qPCR$logit_q),],
                      aes(x = time, y = logit_q, group = chemostat, color = as.factor(prev0))) +
  geom_line(linewidth = 0.3) +
  geom_line(data = data_qPCR_mean[!is.na(data_qPCR_mean$logit_q),], aes(group = prev0), linewidth = 1.3, alpha = 0.6)+
  scale_color_manual(labels = c("epidemic", "endemic"), values = c('#DA0F0F', '#0E5DD6')) +
  labs(x = "time (hours)", y = "Logit-frequency free virulent phage\n", col="Treatment:") +
  scale_x_continuous(limits = c(0,60), breaks = seq(0, 60, 10)) +
  scale_y_continuous(limits = c(-2,9), labels = scales::label_number(accuracy = 1)) +
  mytheme

# With more transparency

Fig_logit_prev.transparent <- ggplot(data = data_FACS[!is.na(data_FACS$logit_prev),],
                                     aes(x = time, y = logit_prev, group = chemostat, color = as.factor(prev0))) +
  geom_line(linewidth = 0.3, alpha = 0.4) +
  geom_line(data = data_FACS_mean[!is.na(data_FACS_mean$logit_prev),], aes(group = prev0), linewidth = 1.3, alpha = 0.3)+
  scale_color_manual(labels = c("epidemic", "endemic"), values = c('#DA0F0F', '#0E5DD6')) +
  labs(x = "time (hours)", y = "Logit-prevalence\n", col="Treatment:") +
  scale_x_continuous(limits = c(0,60), breaks = seq(0, 60, 10)) +
  scale_y_continuous(limits = c(-9,12), labels = scales::label_number(accuracy = 1)) +
  mytheme

Fig_logit_g.transparent <- ggplot(data = data_FACS[!is.na(data_FACS$logit_g),],
                                  aes(x = time, y = logit_g, group = chemostat, color = as.factor(prev0))) +
  geom_line(linewidth = 0.3, alpha = 0.4) +
  geom_line(data = data_FACS_mean[!is.na(data_FACS_mean$logit_g),], aes(group = prev0), linewidth = 1.3, alpha = 0.3)+
  scale_color_manual(labels = c("epidemic", "endemic"), values = c('#DA0F0F', '#0E5DD6')) +
  labs(x = "time (hours)", y = "Logit-frequency infected by virulent phage\n", col="Treatment:") +
  scale_x_continuous(limits = c(0,60), breaks = seq(0, 60, 10)) +
  scale_y_continuous(limits = c(-2,6), labels = scales::label_number(accuracy = 1)) +
  mytheme

Fig_logit_q.transparent <- ggplot(data = data_qPCR[!is.na(data_qPCR$logit_q),],
                                  aes(x = time, y = logit_q, group = chemostat, color = as.factor(prev0))) +
  geom_line(linewidth = 0.3, alpha = 0.4) +
  geom_line(data = data_qPCR_mean[!is.na(data_qPCR_mean$logit_q),], aes(group = prev0), linewidth = 1.3, alpha = 0.3)+
  scale_color_manual(labels = c("epidemic", "endemic"), values = c('#DA0F0F', '#0E5DD6')) +
  labs(x = "time (hours)", y = "Logit-frequency free virulent phage\n", col="Treatment:") +
  scale_x_continuous(limits = c(0,60), breaks = seq(0, 60, 10)) +
  scale_y_continuous(limits = c(-2,9), labels = scales::label_number(accuracy = 1)) +
  mytheme

# save(file = "Data/Experimental_data_mean_and_figures.RData",
#      data_FACS, data_qPCR, data_FACS_mean, data_qPCR_mean,
#      Fig_logit_prev, Fig_logit_g, Fig_logit_q,
#      Fig_logit_prev.transparent, Fig_logit_g.transparent, Fig_logit_q.transparent)