library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(raster)
library(scales)
library(doBy)
library(ggsn)
library(maptools)
library(maps)
library(sp)
library(dplyr)
library(pipeR)
library(plyr)
library(cowplot)
library(lattice)
library(rasterVis)
library(gganimate)
library(reshape2)
library(transformr)
library(ggpubr)
library(gridExtra)
library(rgdal)
library(cowplot)
library(shades)
library(sf)

################################################################################################

####### Plot parameter estimates of final model

# Load models
load("models/mod6_l_pf.R")
load("models/mod6_l_pv.R")

## Create dataframe of values
data <- data.frame(Covariate  = c("Poverty", "Poverty",
                                  "Urban", "Urban",
                                  "Minimum\ntemperature", "Minimum\ntemperature",
                                  "Precipitation", "Precipitation"),
                   
                   mean     = c(mod6_l_pf$summary.fixed$mean[2], mod6_l_pv$summary.fixed$mean[2],
                                mod6_l_pf$summary.fixed$mean[3], mod6_l_pv$summary.fixed$mean[3],
                                mod6_l_pf$summary.fixed$mean[4], mod6_l_pv$summary.fixed$mean[4],
                                mod6_l_pf$summary.fixed$mean[5], mod6_l_pv$summary.fixed$mean[5]),
                   
                   
                   min     = c(mod6_l_pf$summary.fixed$`0.025quant`[2], mod6_l_pv$summary.fixed$`0.025quant`[2],
                               mod6_l_pf$summary.fixed$`0.025quant`[3], mod6_l_pv$summary.fixed$`0.025quant`[3],
                               mod6_l_pf$summary.fixed$`0.025quant`[4], mod6_l_pv$summary.fixed$`0.025quant`[4],
                               mod6_l_pf$summary.fixed$`0.025quant`[5], mod6_l_pv$summary.fixed$`0.025quant`[5]),
                   
                   max     = c(mod6_l_pf$summary.fixed$`0.975quant`[2], mod6_l_pv$summary.fixed$`0.975quant`[2],
                               mod6_l_pf$summary.fixed$`0.975quant`[3], mod6_l_pv$summary.fixed$`0.975quant`[3],
                               mod6_l_pf$summary.fixed$`0.975quant`[4], mod6_l_pv$summary.fixed$`0.975quant`[4],
                               mod6_l_pf$summary.fixed$`0.975quant`[5], mod6_l_pv$summary.fixed$`0.975quant`[5]),
                   
                   Parasite = c("P. falciparum", "P. vivax", 
                                "P. falciparum", "P. vivax",
                                "P. falciparum", "P. vivax",
                                "P. falciparum", "P. vivax"))



data$Covariate <- as.factor(data$Covariate)

l1 <- expression(italic("P. falciparum"), italic("P. vivax"))

### Order for plotting
data$Covariate <- factor(data$Covariate, levels = c("Poverty", 
                                                    "Urban", "Precipitation",
                                                    "Minimum\ntemperature"))

  
  ggplot(data, aes(y = data$mean, x = data$Covariate, colour = data$Parasite)) +
  geom_hline(yintercept = 0, colour = "darkgrey", linetype = "dashed", size = 0.3) +
  geom_linerange(aes(ymin = data$min, ymax = data$max), size = 0.4, position = position_dodge(width = 0.7)) +
  geom_point(aes(x = data$Covariate, y = data$mean, colour = data$Parasite), size = 3,
             position = position_dodge(width = 0.7), shape = 18) +
  theme_classic() + coord_flip() +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size=0.5),
        strip.background = element_blank(),
        legend.position = c(0.18,0.90),
        legend.key = element_blank(),
        legend.background = element_blank(),
        strip.text = element_text(face = "bold.italic")) +
  xlab("") + ylab("Estimate") +
  labs(color = "") +
  scale_colour_manual(labels = l1, values = c("palevioletred", "steelblue"))


################################################################################################
  
####### Plot non-linear climate relationships
### Create df
col1 <- "grey40"
tcol1 <- do.call(rgb,c(as.list(col2rgb(col1)), alpha = 255/4, max = 255))
  
tmin_df_pf <- as.data.frame(mod6_l_pf$summary.random$`inla.group(tmin)`)
tmin_df_pv <- as.data.frame(mod6_l_pv$summary.random$`inla.group(tmin)`)
  
### Unscale temperature values to plot by ddd back mean and times by sd
data <- read.csv("data/inla_input/data.csv")
data_pf <- subset(data, data$parasite == "Falciparum")
data_pv <- subset(data, data$parasite == "Vivax")
  
sd_value <- sd(data_pf$tmin_lag3) 
mean_value <- mean(data_pf$tmin_lag3)
  
## Double check range
#range(data_pf$tmin_lag3)
  
tmin_df_pf$ID_unscaled <- tmin_df_pf$ID*sd_value + mean_value
tmin_df_pv$ID_unscaled <- tmin_df_pv$ID*sd_value + mean_value
  
# Put into df
data_tmin <- data.frame(Parasite = c(rep("P. falciparum", 25),
                                       rep("P. vivax", 25)),
                          x        = c(tmin_df_pf$ID_unscaled,
                                       tmin_df_pv$ID_unscaled),
                          
                          uci        = c(tmin_df_pf$`0.975quant`,
                                         tmin_df_pv$`0.975quant`),
                          
                          lci        = c(tmin_df_pf$`0.025quant`,
                                         tmin_df_pv$`0.025quant`),
                          
                          mean       = c(tmin_df_pf$mean,
                                         tmin_df_pv$mean))
  

ggplot(data_tmin) + 
  geom_ribbon(aes(ymin=exp(data_tmin$lci), ymax=exp(data_tmin$uci), x=data_tmin$x, fill = "95% CI"), alpha = 0.2) +
  geom_hline(yintercept = 1, colour = "grey50", size = 0.3,
            linetype = "dashed") +
    geom_line(aes(x= x, y=exp(mean), colour = "mean")) +
    theme_classic() +
    scale_colour_manual("",values="black") +
    scale_fill_manual("",values=tcol1) +
    xlab("Temperature (°C)") + ylab("Relative risk") +
    theme(axis.line = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          legend.position = "none",
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          legend.key.width = unit(1, "cm"),
          legend.key.height = unit(0.5, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(face = "italic")) +
    facet_wrap(~Parasite, ncol = 1)
  
################################################################################################

####### Look at random effects of models with and without minimum temperature
## Load models
load("models/mod5_l_pf.R")
load("models/mod5_l_pv.R")

load("models/mod6_l_wtmin_pv.R")
load("models/mod6_l_wtmin_pv.R")

t1_df <- data.frame(Model = rep(c("with Tmin",
                                  "without Tmin"), each = 24),
                    Parasite = c(rep("P. falciparum", 12),
                                 rep("P. vivax", 12),
                                 rep("P. falciparum", 12),
                                 rep("P. vivax", 12)),
                    
                    Month = c(rep(1:12, 4)), 
                    mean = c(mod5_l_pf$summary.random$t1$mean,
                             mod5_l_pv$summary.random$t1$mean,
                             mod6_l_wtmin_pf$summary.random$t1$mean,
                             mod6_l_wtmin_pv$summary.random$t1$mean),
                    
                    lci  = c(mod5_l_pf$summary.random$t1$`0.025quant`,
                             mod5_l_pv$summary.random$t1$`0.025quant`,
                             mod6_l_wtmin_pf$summary.random$t1$`0.025quant`,
                             mod6_l_wtmin_pv$summary.random$t1$`0.025quant`),
                    
                    uci  = c(mod5_l_pf$summary.random$t1$`0.975quant`,
                             mod5_l_pv$summary.random$t1$`0.975quant`,
                             mod6_l_wtmin_pf$summary.random$t1$`0.975quant`,
                             mod6_l_wtmin_pv$summary.random$t1$`0.975quant`))

t1_df$Month <- as.factor(t1_df$Month)
levels(t1_df$Month)[levels(t1_df$Month) == "1"]      <- "Jan"
levels(t1_df$Month)[levels(t1_df$Month) == "2"]      <- "Feb"
levels(t1_df$Month)[levels(t1_df$Month) == "3"]      <- "Mar"
levels(t1_df$Month)[levels(t1_df$Month) == "4"]      <- "Apr"
levels(t1_df$Month)[levels(t1_df$Month) == "5"]      <- "May"
levels(t1_df$Month)[levels(t1_df$Month) == "6"]      <- "Jun"
levels(t1_df$Month)[levels(t1_df$Month) == "7"]      <- "Jul"
levels(t1_df$Month)[levels(t1_df$Month) == "8"]      <- "Aug"
levels(t1_df$Month)[levels(t1_df$Month) == "9"]      <- "Sep"
levels(t1_df$Month)[levels(t1_df$Month) == "10"]      <- "Oct"
levels(t1_df$Month)[levels(t1_df$Month) == "11"]      <- "Nov"
levels(t1_df$Month)[levels(t1_df$Month) == "12"]      <- "Dec"

label1 <- bquote(paste("Random effects with T"["min"]))
label2 <- bquote(paste("Random effects without T"["min"]))

ggplot(t1_df, aes(Month, mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin=lci, ymax=uci, colour = Model), position = position_dodge(width = 0.7)) +
  geom_point(aes(fill = Model, colour = Model), shape = 21, position = position_dodge(width = 0.7)) +
  theme_classic() +
  ylab("Relative risk") + 
  theme(axis.line = element_blank(),
        panel.background = element_rect(colour = "black", fill = NA, size=0.5),
        strip.background = element_blank(),
        legend.position = c(0.87,0.90),
        legend.title = element_blank(),
        strip.text.x = element_text(face = "bold.italic", hjust = -0.2),
        strip.text.y = element_text(face = "bold", hjust = -0.2),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(size = 7.5)) +
  facet_wrap(~Parasite, scales = "free", nrow = 1) +
  scale_fill_manual(values = c("salmon", "grey"),
                    labels = c(label1, label2)) +
  scale_colour_manual(values = c("salmon", "grey"),
                      labels = c(label1, label2)) 
  
################################################################################################

####### Compare model posterior predictions with and without climate information
  
  
  