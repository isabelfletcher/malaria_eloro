#############################################################################################################

### Figures for examining the impact of climate and interventions on malaria incidence in El Oro, Ecuador 

#############################################################################################################

## Load libraries
pacman::p_load("raster", "INLA","dplyr", 
               "kableExtra", "reshape2", "ggplot2",
               "scales", "gridExtra", "RColorBrewer",
               "ggsn", "ggpubr")


################################################################################################

## Plot parameter estimates of final model 
load("models/mod6_pf.R")
load("models/mod6_pv.R")

## Create dataframe of values
data <- data.frame(Covariate  = c("Poverty", "Poverty",
                                  "Level of\nurbanization", "Level of\nurbanization",
                                  "Minimum\ntemperature", "Minimum\ntemperature",
                                  "Precipitation", "Precipitation"),
                   
                   mean     = c(mod6_pf$summary.fixed$mean[2], mod6_pv$summary.fixed$mean[2],
                                mod6_pf$summary.fixed$mean[3], mod6_pv$summary.fixed$mean[3],
                                mod6_pf$summary.fixed$mean[4], mod6_pv$summary.fixed$mean[4],
                                mod6_pf$summary.fixed$mean[5], mod6_pv$summary.fixed$mean[5]),
                   
                   
                   min     = c(mod6_pf$summary.fixed$`0.025quant`[2], mod6_pv$summary.fixed$`0.025quant`[2],
                               mod6_pf$summary.fixed$`0.025quant`[3], mod6_pv$summary.fixed$`0.025quant`[3],
                               mod6_pf$summary.fixed$`0.025quant`[4], mod6_pv$summary.fixed$`0.025quant`[4],
                               mod6_pf$summary.fixed$`0.025quant`[5], mod6_pv$summary.fixed$`0.025quant`[5]),
                   
                   max     = c(mod6_pf$summary.fixed$`0.975quant`[2], mod6_pv$summary.fixed$`0.975quant`[2],
                               mod6_pf$summary.fixed$`0.975quant`[3], mod6_pv$summary.fixed$`0.975quant`[3],
                               mod6_pf$summary.fixed$`0.975quant`[4], mod6_pv$summary.fixed$`0.975quant`[4],
                               mod6_pf$summary.fixed$`0.975quant`[5], mod6_pv$summary.fixed$`0.975quant`[5]),
                   
                   Parasite = c("P. falciparum", "P. vivax", 
                                "P. falciparum", "P. vivax",
                                "P. falciparum", "P. vivax",
                                "P. falciparum", "P. vivax"))



data$Covariate <- as.factor(data$Covariate)

l1 <- expression(italic("P. falciparum"), italic("P. vivax"))

### Order for plotting
data$Covariate <- factor(data$Covariate, levels = c("Poverty", 
                                                    "Level of\nurbanization", "Precipitation",
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
        legend.position = c(0.15,0.90),
        legend.key = element_blank(),
        legend.background = element_blank(),
        strip.text = element_text(face = "bold.italic")) +
  xlab("") + ylab("Estimate") +
  labs(color = "") +
  scale_colour_manual(labels = l1, values = c("palevioletred", "steelblue")) +
  scale_y_continuous(limits = c(-1,2),
                     breaks = c(seq(-1,2,0.5)))

################################################################################################

## Non linear relationships of temperature and malaria incidence
load("models/mod5_nl_pf.R")
load("models/mod5_pv.R")

### Create df
col1 <- "grey40"
tcol1 <- do.call(rgb,c(as.list(col2rgb(col1)), alpha = 255/4, max = 255))

tmin_df_pf <- as.data.frame(mod4_pf$summary.random$`inla.group(tmin)`)
tmin_df_pv <- as.data.frame(mod4_pv$summary.random$`inla.group(tmin)`)

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
  geom_ribbon(aes(ymin=data_tmin$lci, ymax=data_tmin$uci, x=data_tmin$x, fill = "95% CI"), alpha = 0.2) +
  geom_hline(yintercept = 0, colour = "grey50", size = 0.3,
             linetype = "dashed") +
  geom_line(aes(x= x, y=mean, colour = "mean")) +
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
  scale_y_continuous(limits = c(-3.5, 3), breaks = c(seq(-3, 3, 1)),
                    labels = c(seq(-3, 3, 1))) +
  scale_x_continuous(limits = c(10, 24), breaks = c(seq(10,24, 2)),
                     labels = c(seq(10,24, 2))) +
  facet_wrap(~Parasite, ncol = 1)

################################################################################################

#### Look at seasonality in malaria incidence explained by temperature
## Load models
# with tmin
load("models/mod6_pf.R")
load("models/mod6_pv.R")

# without tmin
load("models/mod6_wtmin_pf.R")
load("models/mod6_wtmin_pv.R")


## Put into df
t1_df <- data.frame(Model = rep(c("with Tmin",
                                  "without Tmin"), each = 24),
                    Parasite = c(rep("P. falciparum", 12),
                                 rep("P. vivax", 12),
                                 rep("P. falciparum", 12),
                                 rep("P. vivax", 12)),
              
                    Month = c(rep(1:12, 4)), 
                    mean = c(mod6_pf$summary.random$t1$mean,
                             mod6_pv$summary.random$t1$mean,
                             mod6_wtmin_pf$summary.random$t1$mean,
                             mod6_wtmin_pv$summary.random$t1$mean),
                    
                    lci  = c(mod6_pf$summary.random$t1$`0.025quant`,
                             mod6_pv$summary.random$t1$`0.025quant`,
                             mod6_wtmin_pf$summary.random$t1$`0.025quant`,
                             mod6_wtmin_pv$summary.random$t1$`0.025quant`),
                    
                    uci  = c(mod6_pf$summary.random$t1$`0.975quant`,
                             mod6_pv$summary.random$t1$`0.975quant`,
                             mod6_wtmin_pf$summary.random$t1$`0.975quant`,
                             mod6_wtmin_pv$summary.random$t1$`0.975quant`))

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

## Relevel 
levels(t1_df$Parasite)= c("P. falciparum"=expression(paste(bold("A) "), bolditalic("P. falciparum"))),
                          "P. vivax"=expression(paste(bold("B) "), bolditalic("P. vivax"))))

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
        strip.text = element_text(hjust = 0, size = 11),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(size = 7.5)) +
  facet_wrap(~Parasite, scales = "free", nrow = 1,
             labeller = label_parsed) +
  scale_fill_manual(values = c("salmon", "grey"),
                    labels = c(label1, label2)) +
  scale_colour_manual(values = c("salmon", "grey"),
                      labels = c(label1, label2))

  ################################################################################################

  ### Compare model posterior distributions with and without climate information
# Models with climate information
load("models/mod5_l_pf.R")
load("models/mod5_pv.R")

# Models without climate information
load("models/mod3_pf.R") 
load("models/mod3_pv.R") 

## Read in data
data <- read.csv("data/inla_input/data.csv")
# Separate parasite models
data_pf <- subset(data, data$parasite == "Falciparum")
data_pv <- subset(data, data$parasite == "Vivax")


data <- data.frame(Year     = c(rep(data_pf$Year, 2)),
                   
                   Month    = c(rep(data_pf$Month, 2)),
                   
                   Canton  = c(rep(data_pf$Canton, 2)),
                   
                   Population = c(data_pf$Population, data_pv$Population),
                   
                   Parasite = c(rep("P. falciparum", 4872), rep("P. vivax", 4872)),
                   
                   Fit      = c(mod5_l_pf$summary.fitted.values$mean,
                                mod5_pv$summary.fitted.values$mean),
                   
                   Observed = c(data_pf$cases, data_pv$cases),
                   
                   lci      = c(mod5_l_pf$summary.fitted.values[,3],
                                mod5_pv$summary.fitted.values[,3]),
                   
                   uci      =  c(mod5_l_pf$summary.fitted.values[,5],
                                 mod5_pv$summary.fitted.values[,5]))


## Summarise over space
data <- as.data.frame(data %>% group_by(Year, Month, Parasite) %>%
                        dplyr::summarise(Fit = mean(Fit, na.rm = TRUE),
                                         Observed      = mean(Observed, na.rm = TRUE),
                                         lci      = mean(lci, na.rm = TRUE),
                                         uci      = mean(uci, na.rm = TRUE)))

data$Date <- as.Date(with(data, paste(Year, Month, rep(01, nrow(data)), sep = "-")), "%Y-%m-%d")

col1 <- brewer.pal(8, "Set2")[3]
tcol1 <- do.call(rgb,c(as.list(col2rgb(col1)), alpha = 255/4, max = 255))

climate_plot <-
  
  ggplot(data, aes(x = as.Date(Date), y = Observed)) +
  geom_line(aes(colour = "Observed"), linetype = "solid") +
  geom_line(data = data, aes(x = Date, y = Fit, colour = "Fit"), linetype = "dashed") +
  geom_ribbon(data = data, aes(ymin = lci, ymax = uci), alpha = 0.2, colour = tcol1, fill = tcol1) +
  scale_x_date(labels = date_format("%Y"),
               expand = c(0,0),
               breaks = date_breaks("years"),
               date_minor_breaks = "1 month") +
  ylab("") +
  xlab("") +
  theme_classic() +
  facet_grid(~Parasite, scales = "free") +
  scale_colour_manual(name = "legend", values = c("dimgrey", "grey"),
                      breaks = c("Observed", "Fit"),
                      labels = c("Observed", "Modelled")) +
  scale_linetype_manual(name = "legend", values = c("solid", "dashed"),
                        breaks = c("Observed", "Fit"),
                        labels = c("Observed", "Modelled")) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1,2)))) +
  theme(axis.text.x  = element_blank(),
        axis.text.y  = element_text(size = 7),
        axis.ticks.x = element_blank(),
        #legend.text  = element_text(size = 11),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.90, 0.85),
        legend.key.width = unit(2, "line"),
        legend.key.height = unit(1, "line"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.text.x = element_text(face = "bold.italic"),
        strip.background = element_blank()) 


data <- data.frame(Year     = c(rep(data_pf$Year, 2)),
                   
                   Month    = c(rep(data_pf$Month, 2)),
                   
                   Canton  = c(rep(data_pf$Canton, 2)),
                   
                   Parasite = c(rep("P. falciparum", 4872), rep("P. vivax", 4872)),
                   
                   Fit      = c(mod3_pf$summary.fitted.values$mean,
                                mod3_pv$summary.fitted.values$mean),
                   
                   Observed = c(data_pf$cases, data_pv$cases),
                   
                   lci      = c(mod3_pf$summary.fitted.values[,3],
                                mod3_pv$summary.fitted.values[,3]),
                   
                   uci      =  c(mod3_pf$summary.fitted.values[,5],
                                 mod3_pv$summary.fitted.values[,5]))


## Summarise over space
data <- as.data.frame(data %>% group_by(Year, Month, Parasite) %>%
                        dplyr::summarise(Fit = mean(Fit, na.rm = TRUE),
                                         Observed      = mean(Observed, na.rm = TRUE),
                                         lci      = mean(lci, na.rm = TRUE),
                                         uci      = mean(uci, na.rm = TRUE)))

data$Date <- as.Date(with(data, paste(Year, Month, rep(01, nrow(data)), sep = "-")), "%Y-%m-%d")

col1 <- brewer.pal(8, "Set2")[3]
tcol1 <- do.call(rgb,c(as.list(col2rgb(col1)), alpha = 255/4, max = 255))

plot_w_climate <-
  
  ggplot(data, aes(x = as.Date(Date), y = Observed)) +
  geom_line(aes(colour = "Observed"), linetype = "solid") +
  geom_line(data = data, aes(x = Date, y = Fit, colour = "Fit"), linetype = "dashed") +
  geom_ribbon(data = data, aes(ymin = lci, ymax = uci), alpha = 0.2, colour = tcol1, fill = tcol1) +
  scale_x_date(labels = date_format("%Y"),
               expand = c(0,0),
               breaks = date_breaks("years"),
               date_minor_breaks = "1 month") +
  ylab("") +
  xlab("") +
  theme_classic() +
  facet_grid(~Parasite, scales = "free") +
  scale_colour_manual(name = "legend", values = c("dimgrey", "grey"),
                      breaks = c("Observed", "Fit"),
                      labels = c("Observed", "Modelled")) +
  scale_linetype_manual(name = "legend", values = c("solid", "dashed"),
                        breaks = c("Observed", "Fit"),
                        labels = c("Observed", "Modelled")) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1,2)))) +
  theme(axis.text.x  = element_text(angle = 90, hjust = 1, size = 7),
        axis.text.y  = element_text(size = 7),
        #legend.text  = element_text(size = 11),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = 'none',
        legend.key.width = unit(2, "line"),
        legend.key.height = unit(1, "line"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.text = element_blank(),
        strip.background = element_blank()) 

  climate_pred <-
  
  ggarrange(climate_plot,
            plot_w_climate,
            align = "hv",
            labels = c("A) with climate",
                       "B) without climate"),
            font.label = list(size = 11),
            nrow = 2,
            hjust = -0.2)

annotate_figure(climate_pred, 
                left = text_grob("Annual parasite incidence", rot = 90,
                                 vjust = 1.9, hjust = 0.4,
                                 size = 10),
                bottom = text_grob("Time",
                                   vjust = -0.8, hjust = -0.2, size = 10))

################################################################################################

### Compare parameter estimates for 1990-2018 model and intervention model 2001-2015

## 1990-2018 models
load("models/mod6_pf.R")
load("models/mod6_pv.R")

## Intervention models 2001-2015
load("models/mod_int_pf.R")
load("models/mod_int_pv.R")

# Create dataframe of values
data_int <- data.frame(Model = c(rep("1990-2018 model", 14),
                                 rep("2001-2015 model", 14)),
                       
                       Covariate  = c("Minimum\ntemperature", "Minimum\ntemperature",
                                      "Precipitation", "Precipitation",
                                      "Level of\nurbanization", "Level of\nurbanization",
                                      "Poverty", "Poverty",
                                      "Indoor residual\nspraying", "Indoor residual\nspraying",
                                      "Space\nspraying", "Space\nspraying", 
                                      "ULV\nfumigation", "ULV\nfumigation",
                                      
                                      "Minimum\ntemperature", "Minimum\ntemperature",
                                      "Precipitation", "Precipitation",
                                      "Level of\nurbanization", "Level of\nurbanization",
                                      "Poverty", "Poverty",
                                      "Indoor residual\nspraying", "Indoor residual\nspraying",
                                      "Space\nspraying", "Space\nspraying",
                                      "ULV\nfumigation", "ULV\nfumigation"),
                       
                       mean     = c(mod6_pf$summary.fixed$mean[4], mod6_pv$summary.fixed$mean[4],
                                    mod6_pf$summary.fixed$mean[5], mod6_pv$summary.fixed$mean[5],
                                    mod6_pf$summary.fixed$mean[3], mod6_pv$summary.fixed$mean[3],
                                    mod6_pf$summary.fixed$mean[2], mod6_pv$summary.fixed$mean[2],
                                    NA, NA,
                                    NA, NA,
                                    NA, NA,
                                    
                                    mod_int_pf$summary.fixed$mean[2], mod_int_pv$summary.fixed$mean[2],
                                    mod_int_pf$summary.fixed$mean[3], mod_int_pv$summary.fixed$mean[3],
                                    mod_int_pf$summary.fixed$mean[4], mod_int_pv$summary.fixed$mean[4],
                                    mod_int_pf$summary.fixed$mean[5], mod_int_pv$summary.fixed$mean[5],
                                    mod_int_pf$summary.fixed$mean[7], mod_int_pv$summary.fixed$mean[7],
                                    mod_int_pf$summary.fixed$mean[6], mod_int_pv$summary.fixed$mean[6],
                                    mod_int_pf$summary.fixed$mean[8], mod_int_pv$summary.fixed$mean[8]),
                       
                       
                       min     = c(mod6_pf$summary.fixed$`0.025quant`[4], mod6_pv$summary.fixed$`0.025quant`[4],
                                   mod6_pf$summary.fixed$`0.025quant`[5], mod6_pv$summary.fixed$`0.025quant`[5],
                                   mod6_pf$summary.fixed$`0.025quant`[3], mod6_pv$summary.fixed$`0.025quant`[3],
                                   mod6_pf$summary.fixed$`0.025quant`[2], mod6_pv$summary.fixed$`0.025quant`[2],
                                   NA, NA,
                                   NA, NA,
                                   NA, NA,
                                   
                                   mod_int_pf$summary.fixed$`0.025quant`[2], mod_int_pv$summary.fixed$`0.025quant`[2],
                                   mod_int_pf$summary.fixed$`0.025quant`[3], mod_int_pv$summary.fixed$`0.025quant`[3],
                                   mod_int_pf$summary.fixed$`0.025quant`[4], mod_int_pv$summary.fixed$`0.025quant`[4],
                                   mod_int_pf$summary.fixed$`0.025quant`[5], mod_int_pv$summary.fixed$`0.025quant`[5],
                                   mod_int_pf$summary.fixed$`0.025quant`[7], mod_int_pv$summary.fixed$`0.025quant`[7],
                                   mod_int_pf$summary.fixed$`0.025quant`[6], mod_int_pv$summary.fixed$`0.025quant`[6],
                                   mod_int_pf$summary.fixed$`0.025quant`[8], mod_int_pv$summary.fixed$`0.025quant`[8]),
                       
                       max     = c(mod6_pf$summary.fixed$`0.975quant`[4], mod6_pv$summary.fixed$`0.975quant`[4],
                                   mod6_pf$summary.fixed$`0.975quant`[5], mod6_pv$summary.fixed$`0.975quant`[5],
                                   mod6_pf$summary.fixed$`0.975quant`[3], mod6_pv$summary.fixed$`0.975quant`[3],
                                   mod6_pf$summary.fixed$`0.975quant`[2], mod6_pv$summary.fixed$`0.975quant`[2],
                                   NA, NA,
                                   NA, NA,
                                   NA, NA,
                                   
                                   mod_int_pf$summary.fixed$`0.975quant`[2], mod_int_pv$summary.fixed$`0.975quant`[2],
                                   mod_int_pf$summary.fixed$`0.975quant`[3], mod_int_pv$summary.fixed$`0.975quant`[3],
                                   mod_int_pf$summary.fixed$`0.975quant`[4], mod_int_pv$summary.fixed$`0.975quant`[4],
                                   mod_int_pf$summary.fixed$`0.975quant`[5], mod_int_pv$summary.fixed$`0.975quant`[5],
                                   mod_int_pf$summary.fixed$`0.975quant`[7], mod_int_pv$summary.fixed$`0.975quant`[7],
                                   mod_int_pf$summary.fixed$`0.975quant`[6], mod_int_pv$summary.fixed$`0.975quant`[6],
                                   mod_int_pf$summary.fixed$`0.975quant`[8], mod_int_pv$summary.fixed$`0.975quant`[8]),
                       
                       Parasite = c("P. falciparum", "P. vivax", 
                                    "P. falciparum", "P. vivax", 
                                    "P. falciparum", "P. vivax",
                                    "P. falciparum", "P. vivax",
                                    "P. falciparum", "P. vivax", 
                                    "P. falciparum", "P. vivax",
                                    "P. falciparum", "P. vivax",
                                    
                                    "P. falciparum", "P. vivax", 
                                    "P. falciparum", "P. vivax", 
                                    "P. falciparum", "P. vivax",
                                    "P. falciparum", "P. vivax",
                                    "P. falciparum", "P. vivax",
                                    "P. falciparum", "P. vivax",
                                    "P. falciparum", "P. vivax"))

data_int$Covariate <- as.factor(data_int$Covariate)

### Order for plotting
data_int$Covariate <- factor(data_int$Covariate, levels = c("Space\nspraying",
                                                            "ULV\nfumigation",
                                                            "Indoor residual\nspraying",
                                                            "Poverty",
                                                            "Level of\nurbanization",
                                                            "Precipitation",
                                                            "Minimum\ntemperature"))
col1 <- "#08174D"
col2 <- "#339989"

## Relevel 
levels(data_int$Parasite)= c("P. falciparum"=expression(paste(bold("A) "), bolditalic("P. falciparum"))),
                            "P. vivax"=expression(paste(bold("B) "), bolditalic("P. vivax"))))

ggplot(data_int, aes(y = data_int$mean, x = data_int$Covariate, colour = data_int$Model)) +
  geom_hline(yintercept = 0, colour = "darkgrey", linetype = "dashed", size = 0.3) +
  geom_linerange(aes(ymin = data_int$min, ymax = data_int$max), size = 0.4, position = position_dodge(width = 0.7)) +
  geom_point(aes(x = data_int$Covariate, y = data_int$mean, colour = data_int$Model), size = 3,
             position = position_dodge(width = 0.7), shape = 18) +
  theme_classic() + coord_flip() +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_blank(),
        legend.position = c(0.87,0.18),
        legend.key = element_blank(),
        legend.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 11)) +
  xlab("") + ylab("Estimate") +
  labs(color = "") +
  scale_colour_manual(values = c(col1, col2)) + 
  facet_wrap(~Parasite, nrow = 1, labeller = label_parsed) +
  scale_y_continuous(limits = c(-1,1.8),
                     breaks = c(seq(-1,1.8,0.5)))

################################################################################################
### Compare the model improvement for each intervention measure
## Models with all interventions 
load("models/mod_int_pf.R")
load("models/mod_int_nl_pv.R")

## Models without each intervention
load("models/mod_int_w_irs_pf.R")
load("models/mod_int_w_irs_pv.R")
load("models/mod_int_w_fog_pf.R")
load("models/mod_int_w_fog_pv.R")
load("models/mod_int_w_fum_pf.R")
load("models/mod_int_w_fum_pv.R")

## Read in data
data <- read.csv("data/inla_input/data.csv", fileEncoding = "latin1")
## Subset data to intervention period
data <- subset(data, data$Year < 2016 & data$Year > 2000)
# Separate parasite models
data_pf <- subset(data, data$parasite == "Falciparum")
data_pv <- subset(data, data$parasite == "Vivax")

### Models with interventions
data_pf$with_int_pf <- mod_int_pf$summary.fitted.values$`0.5quant`
data_pv$with_int_pv <- mod_int_nl_pv$summary.fitted.values$`0.5quant`
# without IRS
data_pf$without_irs_pf <- mod_int_w_irs_pf$summary.fitted.values$`0.5quant`
data_pv$without_irs_pv <- mod_int_w_irs_pv$summary.fitted.values$`0.5quant`
# without fogging
data_pf$without_fog_pf <- mod_int_w_fog_pf$summary.fitted.values$`0.5quant`
data_pv$without_fog_pv <- mod_int_w_fog_pv$summary.fitted.values$`0.5quant`
# without fumigation
data_pf$without_fum_pf <- mod_int_w_fum_pf$summary.fitted.values$`0.5quant`
data_pv$without_fum_pv <- mod_int_w_fum_pv$summary.fitted.values$`0.5quant`

## Calculate RMSE for each canton
rmse_pf <- as.data.frame(data_pf %>% dplyr::mutate(with_int     = (cases - with_int_pf)^2,
                                                   without_irs  = (cases - without_irs_pf)^2,
                                                   without_fog  = (cases - without_fog_pf)^2,
                                                   without_fum  = (cases - without_fum_pf)^2) %>%
                           group_by(Canton) %>%
                           dplyr::summarise(with_int     = sqrt(mean(with_int, na.rm = T)),
                                            without_irs     = sqrt(mean(without_irs, na.rm = T)),
                                            without_fog     = sqrt(mean(without_fog, na.rm = T)),
                                            without_fum     = sqrt(mean(without_fum, na.rm = T))))

rmse_pf$Parasite <- 'P. falciparum'

rmse_pv <- as.data.frame(data_pv %>% dplyr::mutate(with_int     = (cases - with_int_pv)^2,
                                                   without_irs  = (cases - without_irs_pv)^2,
                                                   without_fog  = (cases - without_fog_pv)^2,
                                                   without_fum  = (cases - without_fum_pv)^2) %>%
                           group_by(Canton) %>%
                           dplyr::summarise(with_int     = sqrt(mean(with_int, na.rm = T)),
                                            without_irs     = sqrt(mean(without_irs, na.rm = T)),
                                            without_fog     = sqrt(mean(without_fog, na.rm = T)),
                                            without_fum     = sqrt(mean(without_fum, na.rm = T))))

rmse_pv$Parasite <- 'P. vivax'

rmse_df <- rbind(rmse_pf, rmse_pv)
rmse_df <- rmse_df %>% dplyr::rename(id = Canton)

## Calculate difference, where positive values indicate improvement in model
rmse_df <- rmse_df %>% dplyr::mutate(difference_irs = without_irs - with_int,
                                     difference_fog = without_fog - with_int,
                                     difference_fum = without_fum - with_int)

# Express as percentage of the original
rmse_df <- rmse_df %>% dplyr::mutate(percent_irs = (difference_irs/without_irs)*100,
                                     percent_fog = (difference_fog/without_fog)*100,
                                     percent_fum = (difference_fum/without_fum)*100)

## Melt to plot
rmse_df <- rmse_df[c(1,6,10:12)]
rmse_df <- melt(rmse_df, id.vars = c("id", "Parasite"))
levels(rmse_df$variable)[levels(rmse_df$variable)=="percent_irs"] <- "Indoor residual spraying"
levels(rmse_df$variable)[levels(rmse_df$variable)=="percent_fog"] <- "Space spraying"
levels(rmse_df$variable)[levels(rmse_df$variable)=="percent_fum"] <- "ULV fumigation"

# Re order for plotting
rmse_df$variable <- factor(rmse_df$variable, levels = c("Indoor residual spraying", 
                                                        "ULV fumigation", 
                                                        "Space spraying"))

rmse_df$value <- replace(rmse_df$value, which(rmse_df$value == 0), NA)

## Replace Chilla canton with NA, as no case data
rmse_df$value[rmse_df$id == "Chilla"] <- NA

ecuador <- getData('GADM', country = "ECU", level = 2)
el_oro  <- subset(ecuador, NAME_1 == "El Oro")
el_oro_f <- fortify(el_oro, region = "NAME_2")
merge_shp_coef <- merge(el_oro_f, rmse_df, by = "id")

## Center colour palette around 0 
limit <- max(abs(merge_shp_coef$value)) * c(-1, 1)

my_palette <- c("#5E083E", "#DDDDDD","#2EA805")

ggplot() + 
  geom_polygon(data = merge_shp_coef, aes(x = long, y = lat, group = group, 
                                          fill = value), color = "black", size = 0.1,
               alpha = 0.6) +
  facet_grid(Parasite~variable, switch = "y") +
  coord_map() +
  theme_void() +
  theme(axis.line = element_blank(),
        panel.background = element_rect(colour = "black", fill = NA, size=0.4),
        strip.background = element_blank(),
        legend.title.align = 0.3,
        legend.text.align = 1,
        strip.text.x = element_text(face = "bold", margin=margin(b=7)),
        strip.text.y = element_text(face = "bold.italic",
                                    vjust = 3.5,
                                    margin=margin(b = 7, l = 7))) +
  scale_fill_gradient2(low = "#5E083E", mid = "#DDDDDD", high = "#2EA805",
                       midpoint = 0, limits = c(-14, max(merge_shp_coef$value, na.rm =T)),
                       name = "Model\nimprovement\n(%)",
                       na.value = "grey43",
                       guide = guide_colourbar(ticks = FALSE, barheight = 5,
                                               barwidth = 1))

################################################################################################

### Compare interannual random effects of 1990-2018 model with 2001-2015 model

# 1990-2018 models
load("models/mod5_l_pf.R")
load("models/mod5_pv.R")

# 2001-2015 model
load("models/mod_int_pf.R")
load("models/mod_int_nl_pv.R")

t2_df <- data.frame(Model = c(rep("Random effects of 1990-2018 model", 58),
                              
                              rep("Random effects of 2001-2015 model", 58)),
                    
                    Parasite = c(rep(c("P. falciparum", "P. vivax"), each = 29),
                                 rep(c("P. falciparum", "P. vivax"), each = 29)),
                    
                    Year = c(1990:2018, 1990:2018,
                             1990:2018, 1990:2018),
                    
                    mean = c(mod5_l_pf$summary.random$t2$mean,
                             mod5_pv$summary.random$t2$mean,
                             rep(NA, 11),
                             mod_int_pf$summary.random$t2$mean,
                             rep(NA, 3),
                             rep(NA, 11),
                             mod_int_nl_pv$summary.random$t2$mean,
                             rep(NA, 3)),
                    
                    lci  = c(mod5_l_pf$summary.random$t2$`0.025quant`,
                             mod5_pv$summary.random$t2$`0.025quant`,
                             rep(NA, 11),
                             mod_int_pf$summary.random$t2$`0.025quant`,
                             rep(NA, 3),
                             rep(NA, 11),
                             mod_int_nl_pv$summary.random$t2$`0.025quant`,
                             rep(NA, 3)),
                    
                    uci  = c(mod5_pf$summary.random$t2$`0.975quant`,
                             mod5_l_pv$summary.random$t2$`0.975quant`,
                             rep(NA, 11),
                             mod_int_pf$summary.random$t2$`0.975quant`,
                             rep(NA, 3),
                             rep(NA, 11),
                             mod_int_nl_pv$summary.random$t2$`0.975quant`,
                             rep(NA, 3)))

col1 <- "#339989"

## Relevel 
levels(t2_df$Parasite)= c("P. falciparum"=expression(paste(bold("A) "), bolditalic("P. falciparum"))),
                          "P. vivax"=expression(paste(bold("B) "), bolditalic("P. vivax"))))

ggplot(t2_df, aes(Year, mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin=lci, ymax=uci, colour = Model), position = position_dodge(width = 0.7)) +
  geom_point(aes(fill = Model, colour = Model), shape = 21, position = position_dodge(width = 0.7)) +
  ylab("Relative risk") + 
  theme_classic() +
  theme(axis.line = element_blank(),
        panel.background = element_rect(colour = "black", fill = NA, size=0.5),
        strip.background = element_blank(),
        legend.position = c(0.67,0.93),
        legend.title = element_blank(),
        strip.text= element_text(hjust = 0, size = 11),
        axis.text.x = element_text(angle = 90, size = 7),
        axis.title = element_text(size = 10),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.4, "cm")) +
  facet_grid(~Parasite, scales = "free", labeller = label_parsed) +
  scale_x_continuous(breaks = 1990:2018) +
  scale_fill_manual(values = c("grey", col1)) +
  scale_colour_manual(values = c("grey", col1)) 
