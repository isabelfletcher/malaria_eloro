## Load libraries
pacman::p_load("dplyr", "ggplot2", "raster",
               "scales", "gridExtra",  "RColorBrewer", 
               "ggpubr", "cowplot", "tidyr", "stringr",
               "mapproj", "reshape2")

##############################################################################################
### Figure S1
data <- read.csv("data.csv")

tiff("supplementary/figure_S1.tif", width = 180, height = 130, units = "mm", res = 520, compression = "lzw")

data %>% subset(Year < 2016 & Year > 2000) %>%
  slice(-seq(0.5 * n())) %>%
  dplyr::group_by(Year, Month) %>%
  dplyr::summarise(fumigation     = sum(blocks_fumigated, na.rm = TRUE),
                   irs            = sum(houses_IRS, na.rm = TRUE),
                   space_spraying = sum(houses_fogged, na.rm = TRUE)) %>%
  melt(id.vars = c("Year", "Month")) %>%
  mutate(Date = str_c(Year, Month, "01", sep = "-"),
         variable = factor(variable, labels = c("A) ULV fumigation", 
                                                "B) Indoor residual spraying", "C) Space spraying"))) %>%
  ggplot(aes(as.Date(Date), value)) +
  geom_line(colour = "darkgrey") +
  theme_classic() +
  xlab("Time") +
  ylab("") +
  scale_x_date(labels = date_format("%Y"),
               expand = c(0,0),
               breaks = date_breaks("years"),
               date_minor_breaks = "1 month") +
  facet_wrap(~variable,
             nrow = 3,
             scales = "free") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.line = element_blank(),
        panel.background = element_rect(colour = "black", fill = NA, size=0.5),
        strip.text = element_text(face="bold", hjust = 0)) +
  scale_y_continuous(labels = comma)

dev.off()

##############################################################################################
### Figure S2
# (nb: best fitting model includes linear climate info for P. falciparum and non-linear for P. vivax)
### 1990-2018 models
load("models/mod5_l_pf.R")
load("models/mod5_nl_pv.R")

### 2001-2015 models
load("models/mod_int_pf.R")
load("models/mod_int_nl_pv.R")

## Full model
data <- read.csv("data.csv")
data <- data %>% mutate(Fit = c(mod5_l_pf$summary.fitted.values$mean, mod5_nl_pv$summary.fitted.values$mean),
                        lci = c(mod5_l_pf$summary.fitted.values[,3], mod5_nl_pv$summary.fitted.values[,3]),
                        uci = c(mod5_l_pf$summary.fitted.values[,5], mod5_nl_pv$summary.fitted.values[,5])) %>%
  # Summarise over space
  dplyr::group_by(Year, Month, parasite) %>%
  dplyr::summarise(Fit = mean(Fit, na.rm = TRUE),
                   Observed   = mean(cases, na.rm = TRUE),
                   lci        = mean(lci, na.rm = TRUE),
                   uci        = mean(uci, na.rm = TRUE),
                   population = mean(Population, na.rm = TRUE)) %>%
  # Calculate API
  mutate(Observed = (Observed*1000)/population,
         Fit      = (Fit*1000)/population,
         lci      = (lci*1000)/population,
         uci      = (uci*1000)/population)

data$parasite <- factor(data$parasite)
levels(data$parasite)[levels(data$parasite)=="Falciparum"] <- "P. falciparum"
levels(data$parasite)[levels(data$parasite)=="Vivax"] <- "P. vivax"

data$Date <- as.Date(with(data, paste(Year, Month, rep(01, nrow(data)), sep = "-")), "%Y-%m-%d")

col1 <- brewer.pal(8, "Set2")[3]
tcol1 <- do.call(rgb,c(as.list(col2rgb(col1)), alpha = 255/4, max = 255))

full_model <-
  
  ggplot(data, aes(x = as.Date(Date), y = Observed)) +
  geom_line(aes(colour = "Observed"), linetype = "solid") +
  geom_line(data = data, aes(x = Date, y = Fit, colour = "Fit"), linetype = "dashed") +
  geom_ribbon(data = data, aes(ymin = lci, ymax = uci), alpha = 0.2, colour = tcol1, fill = tcol1) +
  scale_x_date(labels = date_format("%Y"),
               expand = c(0,0),
               date_breaks = "years",
               date_minor_breaks = "1 month",
               limits = c(as.Date("2001-01-01"), as.Date("2015-12-01"))) +
  ylab("") +
  xlab("") +
  theme_classic() +
  facet_grid(~parasite, scales = "free") +
  scale_colour_manual(name = "legend", values = c("dimgrey", "grey"),
                      breaks = c("Observed", "Fit"),
                      labels = c("Observed", "Modelled")) +
  scale_linetype_manual(name = "legend", values = c("solid", "dashed"),
                        breaks = c("Observed", "Fit"),
                        labels = c("Observed", "Modelled")) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1,2)))) +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.90, 0.85),
        legend.key.width = unit(2, "line"),
        legend.key.height = unit(1, "line"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.text.x = element_text(face = "bold.italic"),
        strip.background = element_blank())


## Intervention model
data <- read.csv("data.csv")
data <- subset(data, data$Year < 2016 & data$Year > 2000)
data <- data %>% mutate(Fit = c(mod_int_pf$summary.fitted.values$mean, mod_int_nl_pv$summary.fitted.values$mean),
                        lci = c(mod_int_pf$summary.fitted.values[,3], mod_int_nl_pv$summary.fitted.values[,3]),
                        uci = c(mod_int_pf$summary.fitted.values[,5], mod_int_nl_pv$summary.fitted.values[,5])) %>%
  # Summarise over space
  dplyr::group_by(Year, Month, parasite) %>%
  dplyr::summarise(Fit = mean(Fit, na.rm = TRUE),
                   Observed   = mean(cases, na.rm = TRUE),
                   lci        = mean(lci, na.rm = TRUE),
                   uci        = mean(uci, na.rm = TRUE),
                   population = mean(Population, na.rm = TRUE)) %>%
  # Calculate API
  mutate(Observed = (Observed*1000)/population,
         Fit      = (Fit*1000)/population,
         lci      = (lci*1000)/population,
         uci      = (uci*1000)/population)

data$parasite <- factor(data$parasite)
levels(data$parasite)[levels(data$parasite)=="Falciparum"] <- "P. falciparum"
levels(data$parasite)[levels(data$parasite)=="Vivax"] <- "P. vivax"

data$Date <- as.Date(with(data, paste(Year, Month, rep(01, nrow(data)), sep = "-")), "%Y-%m-%d")

col1 <- brewer.pal(8, "Set2")[3]
tcol1 <- do.call(rgb,c(as.list(col2rgb(col1)), alpha = 255/4, max = 255))

intervention_model <-
  
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
  facet_grid(~parasite, scales = "free") +
  scale_colour_manual(name = "legend", values = c("dimgrey", "grey"),
                      breaks = c("Observed", "Fit"),
                      labels = c("Observed", "Modelled")) +
  scale_linetype_manual(name = "legend", values = c("solid", "dashed"),
                        breaks = c("Observed", "Fit"),
                        labels = c("Observed", "Modelled")) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1,2)))) +
  theme(axis.text.x  = element_text(angle = 90, hjust = 1, size = 7),
        axis.title = element_text(size = 10),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "none",
        legend.key.width = unit(2, "line"),
        legend.key.height = unit(1, "line"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.text.x = element_text(face = "bold.italic"),
        strip.background = element_blank()) 

plot <-
  
  ggarrange(full_model,
            intervention_model,
            align = "hv",
            labels = c("A) 1990-2018 model",
                       "B) 2001-2015 model"),
            font.label = list(size = 11),
            nrow = 2,
            hjust = -0.1)

tiff("supplementary/figure_S2.tif", height = 130, width = 180, res = 520, units = "mm", compression = "lzw")

annotate_figure(plot, 
                left = text_grob("Annual parasite incidence", rot = 90, size = 10,
                                 vjust = 1.6, hjust = 0.4),
                bottom = text_grob("Time", size = 10,
                                   vjust = -0.8, hjust = -0.2))

dev.off()

##############################################################################################
### Figure S3
data <- read.csv("data.csv", fileEncoding = "latin1")
data <- data %>% dplyr::group_by(Canton) %>%
                 dplyr::summarise(urban = mean(urban, na.rm = TRUE)) %>%
                 dplyr::rename(id = Canton)

data$urban <- cut(data$urban, c(0, 2, 4, 6, 8, 10),
                  labels = c("0-2", "3-4", "5-6", "7-8", "9-10"))

ecuador <- getData('GADM', country = "ECU", level = 2)
el_oro <- ecuador %>% subset(NAME_1 =="El Oro") %>%
  fortify(region = "NAME_2")
el_oro_urban <- merge(el_oro, data, by = "id")

tiff("supplementary/figure_S3.tif", width = 80, height = 80, res = 520, units = "mm", compression = "lzw")

ggplot() + geom_polygon(data = el_oro_urban, aes(x = long, y = lat, group = group, 
                                                   fill = urban), color = "black", size = 0.10) + 
  coord_map() +
  scale_fill_brewer(name = "Urban cover (%)",  palette = "Greys", breaks = rev(levels(el_oro_urban$urban)),
                    drop = FALSE) +
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(colour = "black", fill = NA, size=0.5),
        legend.position = c(0.18, 0.77),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        legend.background = element_blank(),
        legend.key = element_blank()) +
  xlab("") + ylab("") 

dev.off()

##############################################################################################
### Figure S4
data <- read.csv("data.csv", fileEncoding = "latin1")
data <- data %>% subset(Year < 2018) %>%
         dplyr::group_by(Year, parasite, Canton) %>%
         dplyr::summarise(cases      = sum(cases),
                          Population = mean(Population)) %>%
                mutate(API = (cases * 1000)/Population,
                       int = 1) %>%
                mutate(int = ifelse(Year < 2001, 0, 1)) %>% 
        dplyr::group_by(Canton, int, parasite) %>%
        dplyr::summarise(API = mean(API, na.rm = TRUE)) %>%
        dplyr::rename(id = Canton)

int0 <- subset(data, data$int == 0)
int1 <- subset(data, data$int == 1)

int0$API <- cut(int0$API, c(0, 2, 4,  6, 8, 10, 12, 14),
                     labels = c("0-2", "2-4", "4-6", "6-8", "8-10", "10-12",
                                "12-14"))

int0$parasite <- as.factor(int0$parasite)
levels(int0$parasite)[levels(int0$parasite) == "Falciparum"] <- "P. falciparum"
levels(int0$parasite)[levels(int0$parasite) == "Vivax"]      <- "P. vivax"

## Merge data with shapefile to plot
ecuador <- getData('GADM', country = "ECU", level = 2)
el_oro <- ecuador %>% subset(NAME_1 =="El Oro") %>%
  fortify(region = "NAME_2")
el_oro_int0 <- merge(el_oro, int0, by = "id")

int0_map <- 
  
  ggplot() + geom_polygon(data = el_oro_int0, aes(x = long, y = lat, group = group, 
                                                  fill = API), color = "darkgrey", size = 0.10) + 
  coord_map() +
  facet_wrap(~parasite) +
  scale_fill_brewer(name = "Annual parasite \nincidence",  palette = "BuPu", 
                    breaks = rev(levels(el_oro_int0$API)),
                    na.value = "darkgrey") +
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(colour = "black", fill = NA, size=0.5),
        legend.position = c(0.09, 0.74),
        legend.key.size = unit(0.3, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        legend.background = element_blank(),
        legend.key = element_blank()) +
  xlab("") + ylab("")

int1$API <- cut(int1$API, c(0, 1, 2,  3, 4, 5),
                labels = c("0-1", "1-2", "2-3", "3-4", "4-5"))

int1$parasite <- as.factor(int1$parasite)
levels(int1$parasite)[levels(int1$parasite) == "Falciparum"] <- "P. falciparum"
levels(int1$parasite)[levels(int1$parasite) == "Vivax"]      <- "P. vivax"

## Merge data with shapefile to plot
ecuador <- getData('GADM', country = "ECU", level = 2)
el_oro <- ecuador %>% subset(NAME_1 =="El Oro") %>%
  fortify(region = "NAME_2")
el_oro_int1 <- merge(el_oro, int1, by = "id")

int1_map <- 
  
  ggplot() + geom_polygon(data = el_oro_int1, aes(x = long, y = lat, group = group, 
                                                  fill = API), color = "darkgrey", size = 0.10) + 
  coord_map() +
  facet_wrap(~parasite) +
  scale_fill_brewer(name = "Annual parasite \nincidence",  palette = "BuPu", 
                    breaks = rev(levels(el_oro_int1$API)),
                    na.value = "darkgrey",
                    drop = FALSE) +
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(colour = "black", fill = NA, size=0.5),
        legend.position = c(0.09, 0.74),
        legend.key.size = unit(0.3, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        legend.background = element_blank(),
        legend.key = element_blank()) +
  xlab("") + ylab("")

tiff("supplementary/figure_S4.tif", width = 180, height = 180, res = 520, units = "mm", compression = "lzw")

ggarrange(int0_map,
          int1_map,
          align = "hv",
          nrow = 2,
          labels = c("A) 1990-2000",
                     "B) 2001-2015"))


dev.off()

##############################################################################################
### Figure S5
data <- read.csv("data.csv", fileEncoding = "latin1")

data <- data %>% dplyr::group_by(Year, Canton, parasite) %>%
                 dplyr::summarise(cases      = sum(cases, na.rm = TRUE),
                                  Population = mean(Population, na.rm = TRUE),
                                  urban      = mean(urban, na.rm = TRUE)) %>%
                        mutate(API = cases * 1000/Population,
                               urban_category = "Rural") %>%
                        mutate(urban_category = ifelse(urban > 5, "Urban","Rural")) %>% 
                 dplyr::group_by(Year, parasite, urban_category) %>%
                 dplyr::summarise(API = mean(API, na.rm = TRUE))
  
data$parasite <- as.factor(data$parasite)
levels(data$parasite)[levels(data$parasite)=="Falciparum"] <- "P. falciparum"
levels(data$parasite)[levels(data$parasite)=="Vivax"] <- "P. vivax"

col1 <- brewer.pal(8, "Set2")[1]
col2 <- brewer.pal(8, "Set2")[3]

p <- 
ggarrange(data %>% subset(urban_category == "Rural") %>%
          ggplot(aes(Year, API, colour = parasite)) +
            annotate("rect", xmin = 2001, xmax = 2015, ymin = 0, ymax = 12.5, colour = "transparent", fill = 'grey94') +
            geom_line() +
            ylab('') +
            xlab("") +
            theme_classic() +
            theme(strip.background = element_blank(),
                  strip.text = element_text(face = "bold"),
                  legend.title = element_blank(),
                  legend.background = element_blank(),
                  legend.key = element_blank(),
                  legend.text = element_text(face = "italic"),
                  legend.position = "none",
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.line = element_blank(),
                  panel.background = element_rect(colour = "black", fill = NA, size=0.5)) +
            scale_x_continuous(breaks = 1990:2018) +
            scale_colour_manual(values = c("palevioletred", "steelblue")),
          data %>% subset(urban_category == "Urban") %>%
            ggplot(aes(Year, API, colour = parasite)) +
            annotate("rect", xmin = 2001, xmax = 2015, ymin = 0, ymax = 60.5, colour = "transparent", fill = 'grey94') +
            geom_line() +
            ylab('') +
            xlab("Time") +
            theme_classic() +
            theme(strip.background = element_blank(),
                  strip.text = element_text(face = "bold"),
                  legend.title = element_blank(),
                  legend.background = element_blank(),
                  legend.key = element_blank(),
                  legend.text = element_text(face = "italic"),
                  legend.position = c(0.1, 0.85),
                  axis.text.x = element_text(angle = 90),
                  axis.line = element_blank(),
                  panel.background = element_rect(colour = "black", fill = NA, size=0.5)) +
            scale_x_continuous(breaks = 1990:2018) +
            scale_colour_manual(values = c("palevioletred", "steelblue")),
          align = "hv", nrow = 2, ncol = 1)

tiff("supplementary/figure_S5.tif", width = 180, height = 130, units = "mm", res = 520, compression = "lzw")

annotate_figure(p, 
                left = text_grob("Annual parasite incidence", rot = 90,
                                 vjust = 2.5, hjust = 0.4,
                                 size = 11),
                bottom    = text_grob("B) Urban",
                                      vjust = -25, 
                                      hjust = 5.05,
                                      size = 10,
                                      face = "bold"),
                top    = text_grob("A) Rural",
                                   vjust = 1, 
                                   hjust = 5.5,
                                   size = 10,
                                   face = "bold"))

dev.off()

##############################################################################################
### Figure S6
# (nb: best fitting model includes linear climate info for P. falciparum and non-linear for P. vivax)
load("models/mod5_l_pf.R") 
load("models/mod5_nl_pv.R")

## Plot parameter estimates
int0_df_pf <- as.data.frame(mod5_l_pf$marginals.fixed$urban)
int1_df_pf <- as.data.frame(mod5_l_pf$marginals.fixed$urban + abs(mod5_l_pf$marginals.fixed$`urban:int_per1`))

int0_df_pv <- as.data.frame(mod5_nl_pv$marginals.fixed$urban)
int1_df_pv <- as.data.frame(mod5_nl_pv$marginals.fixed$urban + abs(mod5_nl_pv$summary.fixed$mean[5]))

## Remove from y
int1_df_pf$y <-  int1_df_pf$y - abs(as.data.frame(mod5_pf$marginals.fixed$`urban:int_per1`)$y)
int1_df_pv$y <-  int1_df_pv$y - abs(mod5_nl_pv$summary.fixed$mean[5])

## Combine into single df to plot
falciparum_df <- data.frame(int = c(rep("1990-2000", 75),
                                    rep("2001-2015", 75)),
                            x = c(int0_df_pf$x, int1_df_pf$x),
                            y = c(int0_df_pf$y, int1_df_pf$y))

vivax_df <- data.frame(int = c(rep("1990-2000", 75),
                               rep("2001-2015", 75)),
                       x = c(int0_df_pv$x, int1_df_pv$x),
                       y = c(int0_df_pv$y, int1_df_pv$y))

col1 <- brewer.pal(8,"Set2")[3]

param_plot <- 
  
  ggarrange(
    
    ggplot(falciparum_df, aes(x=x, y=y, group = int, colour = "int")) +
      geom_line(aes(colour = int), size = 1) +
      theme_classic() +
      ylab("") +
      xlab("") +
      geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
      theme(axis.line = element_blank(),
            panel.background = element_rect(colour = "black", fill = NA, size=0.5),
            legend.position = c(0.80,0.85),
            legend.background = element_blank(),
            legend.key = element_blank(),
            plot.title = element_text(face = "bold.italic", size = 10,
                                      hjust = 0.5)) +
      scale_colour_manual(values = c("grey", col1),
                          name = "") +
      ggtitle(expression(paste(bold("A) "), bolditalic("P. falciparum")))) ,
    
    
    ggplot(vivax_df, aes(x=x, y=y, group = int, colour = "int")) +
      geom_line(aes(colour = int), size = 1) +
      theme_classic() +
      ylab("") +
      xlab("") +
      geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
      theme(axis.line = element_blank(),
            panel.background = element_rect(colour = "black", fill = NA, size=0.5),
            legend.position = "none",
            axis.title.x = element_text(size = 11),
            plot.title = element_text(face = "bold.italic", size = 10,
                                      hjust = 0.5)) +
      scale_colour_manual(values = c("grey", col1),
                          name = "") +
      ggtitle(expression(paste(bold("B) "), bolditalic("P. vivax")))),
    
    align = "hv", nrow = 1)

tiff("supplementary/figure_S6.tif", width = 180, height = 70, res = 520, units = "mm", compression = "lzw")

annotate_figure(param_plot, 
                left = text_grob("Density", rot = 90,
                                 vjust = 1.9, hjust = 0.4,
                                 size = 12),
                bottom = text_grob("Parameter estimate", 
                                   vjust = -0.9, 
                                   hjust = 0.32,
                                   size = 12))
dev.off()

##############################################################################################
### Figure S7
data <- read.csv("data.csv", fileEncoding = "latin1")
data <- data %>% subset(Year < 2016) %>%
         mutate(incidence = (cases*1000)/Population,
                int = 1) %>%
         mutate(int = ifelse(Year < 2001, 0, 1)) %>% 
  dplyr::group_by(Month, parasite, int) %>%
  dplyr::summarise(incidence = mean(incidence, na.rm = TRUE)) %>%
  mutate(int = factor(int, labels = c("1990 - 2000", "2001 - 2015"))) 

# Labels
data$parasite <- as.factor(data$parasite)
levels(data$parasite)= c("P. falciparum"=expression(paste(bold("A) "), bolditalic("P. falciparum"))),
                         "P. vivax"=expression(paste(bold("B) "), bolditalic("P. vivax"))))  

tiff("supplementary/figure_S7.tif", width = 180, height = 80, units = "mm", res = 520, compression = "lzw")

ggplot(data, aes(x = Month, y = incidence, group = int, linetype = int), colour = "black") +
  geom_line() +
  theme_classic() +
  theme(strip.text = element_text(hjust = 0, size = 11),
        legend.position = c(0.1,0.9),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank()) +
  labs(x = "Month", y = "Incidence (per 1,000)") +
  scale_y_continuous(expand = c(0,0.0),
                     limits = c(0, 0.65)) +
  scale_x_continuous(breaks = c(1:12), labels = month.abb) +
  facet_wrap(~parasite, scales = "free",
             labeller = label_parsed)

dev.off()

##############################################################################################
### Figure S8
## Climate models
load("models/mod5_l_pf.R")
load("models/mod5_nl_pv.R")

## Models without climate
load("models/mod6_wtmin_pf.R")
load("models/mod6_wtmin_pv.R")
load("models/mod6_wprcp_pf.R")
load("models/mod6_wprcp_pv.R")

data <- read.csv("data.csv", fileEncoding = "latin1")

rmse_df <- 
  rbind(mod5_l_pf$summary.fitted.values %>%
                   dplyr::select(`0.5quant`) %>%
                   mutate(observed  = c(subset(data, data$parasite == "Falciparum")$cases),
                          climate   = c(mod5_l_pf$summary.fitted.values$`0.5quant`),
                          parasite  = "P. falciparum",
                          model     = "with_climate",
                          id        = subset(data, data$parasite == "Falciparum")$Canton) %>%
                   tibble::remove_rownames() %>% 
                   dplyr::rename(fit = `0.5quant`),
        mod5_nl_pv$summary.fitted.values %>%
                   dplyr::select(`0.5quant`) %>%
                   mutate(observed  = c(subset(data, data$parasite == "Vivax")$cases),
                          climate   = c(mod5_nl_pv$summary.fitted.values$`0.5quant`),
                          parasite  = "P. vivax",
                          model     = "with_climate",
                          id        = subset(data, data$parasite == "Vivax")$Canton) %>%
                   tibble::remove_rownames() %>% 
                   dplyr::rename(fit = `0.5quant`),
        mod6_wtmin_pf$summary.fitted.values %>%
          dplyr::select(`0.5quant`) %>%
          mutate(observed  = c(subset(data, data$parasite == "Falciparum")$cases),
                 climate   = c(mod5_l_pf$summary.fitted.values$`0.5quant`),
                 parasite  = "P. falciparum",
                 model     = "Minimum temperature",
                 id        = subset(data, data$parasite == "Falciparum")$Canton) %>%
          tibble::remove_rownames() %>% 
          dplyr::rename(fit = `0.5quant`),
        mod6_wtmin_pv$summary.fitted.values %>%
          dplyr::select(`0.5quant`) %>%
          mutate(observed  = c(subset(data, data$parasite == "Vivax")$cases),
                 climate   = c(mod5_nl_pv$summary.fitted.values$`0.5quant`),
                 parasite  = "P. vivax",
                 model     = "Minimum temperature",
                 id        = subset(data, data$parasite == "Vivax")$Canton) %>%
          tibble::remove_rownames() %>% 
          dplyr::rename(fit = `0.5quant`),
        mod6_wprcp_pf$summary.fitted.values %>%
          dplyr::select(`0.5quant`) %>%
          mutate(observed  = c(subset(data, data$parasite == "Falciparum")$cases),
                 climate   = c(mod5_l_pf$summary.fitted.values$`0.5quant`),
                 parasite  = "P. falciparum",
                 model     = "Precipitation",
                 id        = subset(data, data$parasite == "Falciparum")$Canton) %>%
          tibble::remove_rownames() %>% 
          dplyr::rename(fit = `0.5quant`),
        mod6_wprcp_pv$summary.fitted.values %>%
          dplyr::select(`0.5quant`) %>%
          mutate(observed  = c(subset(data, data$parasite == "Vivax")$cases),
                 climate   = c(mod5_nl_pv$summary.fitted.values$`0.5quant`),
                 parasite  = "P. vivax",
                 model     = "Precipitation",
                 id        = subset(data, data$parasite == "Vivax")$Canton) %>%
          tibble::remove_rownames() %>% 
          dplyr::rename(fit = `0.5quant`))

# Calculate rmse
rmse_df <- rmse_df %>% mutate(e         = (observed-fit)^2,
                              e_climate = (observed-climate)^2) %>%
  dplyr::group_by(id, parasite, model) %>%
  dplyr::summarise(rmse         = sqrt(mean(e, na.rm = TRUE)),
                   rmse_climate = sqrt(mean(e_climate, na.rm = TRUE))) %>%
  mutate(rmse_difference = (rmse - rmse_climate)/rmse * 100) %>%
  dplyr::select(id, parasite, model, rmse_difference) %>%
  subset(model != "with_climate")

# No data for canton so replace with NA
rmse_df$rmse_difference[rmse_df$id == "Chilla"] <- NA

## Merge data with shapefile to plot
ecuador <- getData('GADM', country = "ECU", level = 2)
el_oro <- ecuador %>% subset(NAME_1 =="El Oro") %>%
  fortify(region = "NAME_2")
el_oro_rmse <- merge(el_oro, rmse_df, by = "id")

my_palette <- c("#5E083E", "#DDDDDD","#2EA805")

load("functions/align_legend.RData")

p <- ggplot() + 
  geom_polygon(data = el_oro_rmse, aes(x = long, y = lat, group = group, 
                                       fill = rmse_difference), color = "black", size = 0.1,
               alpha = 0.6) +
  facet_grid(parasite~model, switch = "y") +
  theme_void() +
  theme(panel.background = element_rect(colour = "black", fill = NA, size=0.4),
        strip.background = element_blank(),
        legend.title.align = 0.3,
        legend.text.align = 1,
        strip.text.x = element_text(face = "bold", margin=margin(b=7)),
        strip.text.y = element_text(face = "bold.italic",
                                    vjust = 3.5,
                                    margin=margin(b = 7, l = 7))) +
  scale_fill_gradient2(low = "#5E083E", mid = "#DDDDDD", high = "#2EA805",
                       midpoint = 0, limits = c(min(el_oro_rmse$rmse_difference, na.rm =T), 32.5),
                       name = "Model\nimprovement\n(%)",
                       na.value = "grey43",
                       guide = guide_colourbar(ticks = FALSE, barheight = 5,
                                               barwidth = 1))

tiff("supplementary/figure_S8.tif", width = 180, height = 180, units = "mm", res = 520, compression = "lzw")
ggdraw(align_legend(p))
dev.off()


##############################################################################################
### Figure S9
# Models with tmin (nb: best fitting model includes linear climate info for P. falciparum and non-linear for P. vivax)
load("models/mod5_l_pf.R")
load("models/mod5_nl_pv.R")

## Models without tmin
load("models/mod6_wtmin_pf.R")
load("models/mod6_wtmin_pv.R")

t2_df <- rbind(rbind(mod5_l_pf$summary.random$t2 %>%
                       mutate(parasite = "P. falciparum",
                              model    = "Random effects with Tmin") %>%
                       dplyr::select(ID, mean, `0.025quant`, `0.975quant`, parasite, model) %>%
                       dplyr::rename(year = ID,
                                     lci = `0.025quant`,
                                     uci = `0.975quant`),
                     mod5_nl_pv$summary.random$t2 %>%
                       mutate(parasite = "P. vivax",
                              model    = "Random effects with Tmin") %>%
                       dplyr::select(ID, mean, `0.025quant`, `0.975quant`, parasite, model) %>%
                       dplyr::rename(year = ID,
                                     lci = `0.025quant`,
                                     uci = `0.975quant`)),
               rbind(mod6_wtmin_pf$summary.random$t2 %>%
                       mutate(parasite = "P. falciparum",
                              model    = "Random effects without Tmin") %>%
                       dplyr::select(ID, mean, `0.025quant`, `0.975quant`, parasite, model) %>%
                       dplyr::rename(year = ID,
                                     lci = `0.025quant`,
                                     uci = `0.975quant`),
                     mod6_wtmin_pv$summary.random$t2 %>%
                       mutate(parasite = "P. vivax",
                              model    = "Random effects without Tmin") %>%
                       dplyr::select(ID, mean, `0.025quant`, `0.975quant`, parasite, model) %>%
                       dplyr::rename(year = ID,
                                     lci = `0.025quant`,
                                     uci = `0.975quant`)))

# Labels
label1 <- bquote(paste("Random effects with T"["min"]))
label2 <- bquote(paste("Random effects without T"["min"]))

t2_df$parasite <- as.factor(t2_df$parasite)
levels(t2_df$parasite)= c("P. falciparum"=expression(paste(bold("A) "), bolditalic("P. falciparum"))),
                          "P. vivax"=expression(paste(bold("B) "), bolditalic("P. vivax"))))

tiff("supplementary/figure_S9.tif", width = 180, height = 90, res = 520, units = "mm", compression = "lzw")

ggplot(t2_df, aes(year, mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin=lci, ymax=uci, colour = model), position = position_dodge(width = 0.7)) +
  geom_point(aes(fill = model, colour = model), shape = 21, position = position_dodge(width = 0.7)) +
  theme_classic() +
  ylab("Relative risk") + xlab("Year") +
  theme(axis.line = element_blank(),
        panel.background = element_rect(colour = "black", fill = NA, size=0.5),
        strip.background = element_blank(),
        legend.position = c(0.87,0.90),
        legend.title = element_blank(),
        strip.text = element_text(hjust = 0, size = 11),
        axis.text.x = element_text(angle = 90, size = 7),
        axis.title = element_text(size = 10),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(size = 7.5)) +
  facet_grid(~parasite, scales = "free",
             labeller = label_parsed) +
  scale_fill_manual(values = c("salmon", "grey"),
                    labels =c(label1, label2)) +
  scale_colour_manual(values = c("salmon", "grey"),
                      labels =c(label1, label2)) 

dev.off()

##############################################################################################
### Figure S10
data <- read.csv("data.csv")

plot_a <- data %>% dplyr::group_by(Year, Month) %>%
                 dplyr::summarise(tmin = mean(tmin, na.rm = TRUE)) %>%
ggplot(aes(x = Month, y = tmin, group = Year, colour = Year)) +
  geom_line(position = position_dodge(width = 0.8), size = 0.4) +
  theme_classic() +
  theme(strip.text = element_text(face = "bold.italic"),
        legend.position = c(0.82,0.92),
        legend.text  = element_text(size = 7),
        legend.title = element_text(size = 8),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 10)) +
  scale_y_continuous(limits = c(15.3, 20.2)) +
  scale_x_continuous(breaks = c(1:12), labels = month.abb) +
  xlab("Month") + ylab("Minimum temperature (°C)") +
  scale_color_distiller(palette = 'YlOrRd', direction = -1,
                        breaks = c(1990, 2004, 2018),
                        trans = "reverse",
                        guide = guide_colourbar(ticks = FALSE, 
                                                direction = "horizontal",
                                                reverse = TRUE,
                                                title = "Year",
                                                title.position = "top",
                                                title.hjust = 0.5,
                                                barheight = 0.5,
                                                barwidth = 4))

col1 <- brewer.pal(9, "Reds")[7]

plot_b <- data %>% dplyr::group_by(Year, Month) %>%
  dplyr::summarise(tmin = mean(tmin, na.rm = TRUE),
                   prcp = mean(prcp, na.rm = TRUE),
                   prcp_day = mean(prcp_day, na.rm = TRUE)) %>%
  mutate(Date = str_c(Year, Month, "01", sep = "-")) %>%
  ggplot(aes(x = as.Date(Date), y = tmin)) + 
  geom_line(size = 0.4, colour = col1) + 
  theme_classic() +
  xlab("Time") + ylab("Minimum temperature (°C)") +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        axis.text.x  = element_text(angle = 90, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 10),
        axis.line = element_blank(),
        panel.background = element_rect(colour = "black", fill = NA, size=0.5)) +
  scale_y_continuous(limits = c(15.3, 20.2)) +
  scale_x_date(labels = date_format("%Y"),
               expand = c(0,0),
               breaks = date_breaks("years"),
               date_minor_breaks = "1 month") +
  geom_smooth(method = "lm", colour = "grey20", 
              linetype = "dashed",
              size = 0.2)

tiff("supplementary/figure_S10.tif", width = 180, height = 95, units = "mm", res = 520, compression = "lzw")
ggarrange(plot_a,
          plot_b,
          align = "v",
          nrow = 1, ncol = 2,
          labels = c("A)", "B)"))
dev.off()

##############################################################################################
### Table S1
data <- rbind(read.csv("supplementary/table_S1_pf.csv"),
              read.csv("supplementary/table_S1_pv.csv")) %>% 
         mutate(Parasite = rep(c("P. falciparum", "P. vivax"), each = 12))

data$Intervention <- factor(data$Intervention, levels = c("IRS", "Fumigation", "Fogging"),
                            labels = c("Indoor residual spraying", "ULV fumigation", "Space spraying"))
data <- data[order(data$Intervention, data$Parasite, data$Lag),]
data <- data %>% mutate(Mean = round(Mean, 2),
                        LCI  = round(LCI, 2),
                        UCI  = round(UCI, 2),
                        DIC  = round(DIC, 2),
                        WAIC = round(WAIC, 2))


write.csv(data, "supplementary/table_S1.csv")

##############################################################################################
### Table S2
data <- rbind(read.csv("supplementary/table_s2_tmin_lags_l.csv"),
              read.csv("supplementary/table_s2_tmax_lags_l.csv"),
              read.csv("supplementary/table_s2_prcp_lags_l.csv")) %>% 
  mutate(Variable = rep(c("Minimum temperature", "Maximum temperature",
                          "Precipitation"), each = 8)) %>%
  mutate(Estimate = round(Estimate, 2),
         LCI  = round(LCI, 2),
         UCI  = round(UCI, 2),
         DIC  = round(DIC, 2),
         WAIC = round(WAIC, 2))

write.csv(data, "supplementary/table_S2.csv")

##############################################################################################
### Table S3
data <- rbind(read.csv("supplementary/table_s3_tmin_lags_nl.csv"),
              read.csv("supplementary/table_s3_tmax_lags_nl.csv"),
              read.csv("supplementary/table_s3_prcp_lags_nl.csv")) %>% 
  mutate(Variable = rep(c("Minimum temperature", "Maximum temperature",
                          "Precipitation"), each = 8)) %>%
  mutate(DIC  = round(DIC, 2),
         WAIC = round(WAIC, 2))

write.csv(data, "supplementary/table_S3.csv")

##############################################################################################
### Table S4 
# 1990-2018 models
load("models/mod6_pf.R")
load("models/mod6_pv.R")

# Intervention models 2001-2015
load("models/mod_int_pf.R")
load("models/mod_int_pv.R")

estimates_df <- rbind(rbind(mod6_pf$summary.fixed %>%
                           dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
                           mutate(Parasite = "P. falciparum",
                                  Model    = "1990-2018",
                                  Variable = as.factor(rownames(mod6_pf$summary.fixed))) %>%
                           dplyr::slice(-1) %>% 
                           droplevels() %>%
                           mutate(Variable = factor(Variable, labels = c("Precipitation", "Minimum temperature", "Poverty", "Level of urbanization"))) %>%
                           dplyr::rename(LCI = `0.025quant`,
                                         UCI = `0.975quant`,
                                         Estimate = mean),
                         mod6_pv$summary.fixed %>%
                           dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
                           mutate(Parasite = "P. vivax",
                                  Model    = "1990-2018",
                                  Variable = as.factor(rownames(mod6_pv$summary.fixed))) %>%
                           dplyr::slice(-1) %>% 
                           droplevels() %>%
                           mutate(Variable = factor(Variable, labels = c("Precipitation", "Minimum temperature", "Poverty", "Level of urbanization"))) %>%
                           dplyr::rename(LCI = `0.025quant`,
                                         UCI = `0.975quant`,
                                         Estimate = mean)),
                   rbind(mod_int_pf$summary.fixed %>%
                           dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
                           mutate(Parasite = "P. falciparum",
                                  Model    = "2001-2015",
                                  Variable = as.factor(rownames(mod_int_pf$summary.fixed))) %>%
                           dplyr::slice(-1, -6:-8) %>% 
                           droplevels() %>%
                           mutate(Variable = factor(Variable, labels = c("Precipitation", "Minimum temperature", "Poverty", "Level of urbanization"))) %>%
                           dplyr::rename(LCI = `0.025quant`,
                                         UCI = `0.975quant`,
                                         Estimate = mean),
                         mod_int_pv$summary.fixed %>%
                           dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
                           mutate(Parasite = "P. vivax",
                                  Model    = "2001-2015",
                                  Variable = as.factor(rownames(mod_int_pv$summary.fixed))) %>%
                           dplyr::slice(-1, -6:-8) %>% 
                           droplevels() %>%
                           mutate(Variable = factor(Variable, labels = c("Precipitation", "Minimum temperature", "Poverty", "Level of urbanization"))) %>%
                           dplyr::rename(LCI = `0.025quant`,
                                         UCI = `0.975quant`,
                                         Estimate = mean)))

estimates_df <- estimates_df %>% mutate(Estimate = round(Estimate, 2),
                                        LCI      = round(LCI, 2),
                                        UCI      = round(UCI, 2))

estimates_df$Variable <- factor(estimates_df$Variable, levels = c("Minimum temperature", "Precipitation", "Level of urbanization", "Poverty"))
estimates_df <- estimates_df[order(estimates_df$Variable, estimates_df$Parasite, estimates_df$Model),]

write.csv(estimates_df, file = "supplementary/table_S4.csv")

##############################################################################################
### Table S5 
load("models/mod_int_pf.R")
load("models/mod_int_pv.R")

estimates_df <- rbind(mod_int_pf$summary.fixed %>%
                     dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
                     mutate(Parasite = "P. falciparum",
                            Variable = as.factor(rownames(mod_int_pf$summary.fixed))) %>%
                     dplyr::slice(-1:-5) %>% 
                     droplevels() %>%
                     mutate(Variable = factor(Variable, labels = c("ULV fumigation", "Space spraying", "Indoor residual spraying"))) %>%
                     dplyr::rename(Estimate = mean,
                                   LCI      = `0.025quant`,
                                   UCI      = `0.975quant`),
                     mod_int_pv$summary.fixed %>%
                     dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
                     mutate(Parasite = "P. vivax",
                            Variable = as.factor(rownames(mod_int_pf$summary.fixed))) %>%
                     dplyr::slice(-1:-5) %>% 
                     droplevels() %>%
                     mutate(Variable = factor(Variable, labels = c("ULV fumigation", "Space spraying", "Indoor residual spraying"))) %>%
                     dplyr::rename(Estimate = mean,
                                   LCI      = `0.025quant`,
                                   UCI      = `0.975quant`))

estimates_df <- estimates_df %>% mutate(Estimate = round(Estimate, 2),
                                        LCI      = round(LCI, 2),
                                        UCI      = round(UCI, 2))

estimates_df$Variable <- factor(estimates_df$Variable, levels = c("Indoor residual spraying", "ULV fumigation", "Space spraying"))
estimates_df <- estimates_df[order(estimates_df$Variable, estimates_df$Parasite),]


write.csv(estimates_df, file = "supplementary/table_S5.csv")
