#############################################################################################################

### Figures for examining the impact of climate and interventions on spatiotemporal incidence of malaria in El Oro, Ecuador 

#############################################################################################################

## Load libraries
pacman::p_load("raster", "INLA","dplyr", 
               "kableExtra", "reshape2", "ggplot2",
               "scales", "gridExtra", "RColorBrewer",
               "ggsn", "ggpubr", "sf", "ggsn", "shades",
               "cowplot")


################################################################################################

### Figure 1

ecuador_0 <- getData('GADM', country = "ECU", level = 0)
ecuador_0 <- st_as_sf(ecuador_0)

ecuador_1 <- getData('GADM', country = "ECU", level = 1)
ecuador_1 <- st_as_sf(ecuador_1)
el_oro <- subset(ecuador_1, ecuador_1$NAME_1 == "El Oro")

venezuela <- getData('GADM', country = "VEN", level = 0)
venezuela <- st_as_sf(venezuela)

peru <- getData('GADM', country = "PER", level = 0)
peru <- st_as_sf(peru)

brazil <- getData('GADM', country = "BRA", level = 0)
brazil <- st_as_sf(brazil)

colombia <- getData('GADM', country = "COL", level = 0)
colombia <- st_as_sf(colombia)

## Read in data
cases <- read.csv("data/malaria_cases/malaria_cases.csv")

## Total amount of cases
sum_cases <- as.data.frame(cases %>% 
                             group_by(Year, Month) %>%
                             dplyr::summarise(total_cases      = sum(Total_cases, na.rm = TRUE),
                                              total_falciparum = sum(Falciparum, na.rm = TRUE),
                                              total_vivax      = sum(Vivax, na.rm = TRUE)))

# Combine into single variable
sum_cases <- melt(sum_cases, id.vars = c("Year", "Month"), 
                  measure.vars = c("total_cases", "total_falciparum", "total_vivax"),
                  value.name = "Cases", variable.name = "Type")

## Date column
sum_cases$Date <- as.Date(with(sum_cases, paste(Year, Month, rep(01, nrow(sum_cases)), sep = "-")),
                          "%Y-%m-%d")

## Plot time series
l1 <- expression(italic("P. falciparum"), italic("P. vivax"), "Total")

col1 <- saturation("grey", 0.5)

cases_plot <- 
  
  ggplot(sum_cases, aes(x = Date, y = Cases)) + 
  geom_line(aes(colour = Type, group = Type), 
            alpha = 0.7) +
  theme_classic() +
  xlab("Time") +
  theme(legend.title = element_blank(),
        axis.text.x  = element_text(angle = 90, hjust = 1),
        legend.position = c(0.1, 0.85),
        legend.key = element_blank(),
        legend.background = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_colour_manual(name = "legend", values = c("palevioletred", "steelblue","lightgrey"),
                      breaks = c("total_falciparum", "total_vivax", "total_cases"),
                      labels = l1) +
  scale_x_date(labels = date_format("%Y"),
               expand = c(0,0),
               breaks = date_breaks("years"),
               date_minor_breaks = "1 month") +
  guides(colour = guide_legend(override.aes = list(linetype = 1))) +
  guides(col = guide_legend(ncol = 1)) +
  annotate("segment", x = as.Date('2001-01-01'), xend = as.Date('2015-01-01'), y = 950, yend = 950, colour = "black", arrow=arrow(ends = "both")) +
  annotate("text", x = as.Date('2008-01-01'), y = 1050, label = "Period of intensive vector control", size = 3)


#########

## Map
sa <- sf::st_read("data/SA_shp/SouthAmerica.shp")

map <-
  ggplot() + geom_sf(data=sa, fill = "grey95", size = 0.1,
                     colour = "black") +
  ## Add country labels
  annotate('text', x = -9592560, y = -1008478, 
           label = "Peru", size = 3.5) +
  annotate("segment", x = -8400000, xend = -9300000,
           y = -1000000, yend = -1000000,
           size = 0.2) +
  
  annotate('text', x = -9800000, y = 300000, 
           label = "Ecuador", size = 3.5) +
  annotate("segment", x = -8700000, xend = -9500000,
           y = -110000, yend = 100000,
           size = 0.2) +
  
  annotate('text', x = -9450000, y = 1300000, 
           label = "Colombia", size = 3.5) +
  annotate("segment", x = -8250000, xend = -9500000,
           y = 700000, yend = 1100000,
           size = 0.2) +
  
  annotate('text', x = -7550000, y = 1800000, 
           label = "Venezuela", size = 3.5) +
  annotate("segment", x = -7250000, xend = -7650000,
           y = 900000, yend = 1650000,
           size = 0.2) +
  
  scale_y_continuous(limits = c(-4000000, 1856462)) +
  geom_sf(data = el_oro, fill = "salmon", size = 0.05, colour = "salmon") +
  theme_void() +
  ## Add box around location
  annotate("rect", xmin = -8200000, xmax = -9300000,
           ymin = 320000, ymax=-750000, colour = "black",
           fill = "transparent", size = 0.3) +
  
  ## Scalebar and north
  north(sa, symbol = 3, scale = 0.08,
        anchor = c(x=-3550000, y=-2700000)) 

map <- 
  
map +
  scalebar(sa, dist_unit = "km",
           dist = 1000, transform = FALSE,
           dd2km = FALSE,
           model = "WGS84",
           st.dist = 0.0105, 
           height = 0.02, st.size = 3,
           border.size = 0.4,
           anchor = c(x = -3550000, y= -3900000))


### Second map
el_oro <- ggplot() + geom_sf(data=ecuador_0, fill = "grey95", size = 0.1,
                   colour = "black") +
  geom_sf(data = el_oro, fill = "salmon", size = 0.05, colour = "black") +
  geom_sf(data = colombia, fill = "grey95", size = 0.1,
          colour = "black") +
  geom_sf(data = peru, fill = "grey95", size = 0.1,
          colour = "black") +
  geom_sf(data = brazil, fill = "grey95", size = 0.1,
          colour = "black") +
  scale_x_continuous(limits = c(-83.5, -74)) +
  scale_y_continuous(limits = c(-6, 3)) +
  theme_void() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.7))  +
  annotate('text', x = -82.3, y = -2, 
           label = "El Oro", size = 3) +
  annotate("segment", x = -82, xend = -79.8,
           y = -2.8, yend = -3.4,
           size = 0.4) 

############################################
g <- ggplotGrob(el_oro)

map_g <- 
map + annotation_custom(grob = g, 
                        xmin = -13000000,
                        xmax = -5700000,
                        ymin = -2000000,
                        ymax = -4200000)

## Arrange with cases
tiff("figures/figure_1.tif", width = 180, height = 180, res = 360, units = "mm", compression = "lzw")
ggarrange(map_g, cases_plot,
          nrow = 2,
          align = 'v',
          heights = c(2, 1.5),
          labels = c("A)", "B)"),
          hjust = -0.20)
dev.off()

################################################################################################
### Figure 2A

## Plot parameter estimates of final model
load("models/mod6_pf.R")
load("models/mod6_pv.R")

estimates <- rbind(mod6_pf$summary.fixed %>%
                     dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
                     mutate(parasite = "P. falciparum",
                            variable = rownames(mod6_pf$summary.fixed)) %>%
                     dplyr::slice(-1) %>% 
                     dplyr::rename(lci = `0.025quant`,
                                   uci = `0.975quant`),
                   mod6_pv$summary.fixed %>%
                     dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
                     mutate(parasite = "P. vivax",
                            variable = rownames(mod6_pf$summary.fixed)) %>%
                     dplyr::slice(-1) %>% 
                     dplyr::rename(lci = `0.025quant`,
                                   uci = `0.975quant`))


l1 <- expression(italic("P. falciparum"), italic("P. vivax"))

### Order and label for plotting
estimates$variable <- as.factor(estimates$variable)
estimates$variable <- factor(estimates$variable, 
                             levels = c("total_poverty", 
                                        "urban", "prcp",
                                        "tmin"),
                             labels = c("Poverty", 
                                        "Level of\nurbanization", "Precipitation",
                                        "Minimum\ntemperature"))
estimates_plot <- 
  ggplot(estimates, aes(y = estimates$mean, x = estimates$variable, colour = estimates$parasite)) +
  geom_hline(yintercept = 0, colour = "darkgrey", linetype = "dashed", size = 0.3) +
  geom_linerange(aes(ymin = estimates$lci, ymax = estimates$uci), size = 0.4, position = position_dodge(width = 0.7)) +
  geom_point(aes(x = estimates$variable, y = estimates$mean, colour = estimates$parasite), size = 3,
             position = position_dodge(width = 0.7), shape = 18) +
  theme_classic() + coord_flip() +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size=0.5),
        strip.background = element_blank(),
        legend.position = c(0.87,0.93),
        legend.key = element_blank(),
        legend.background = element_blank(),
        strip.text = element_text(face = "bold.italic")) +
  xlab("") + ylab("Estimate") +
  labs(color = "") +
  scale_colour_manual(labels = l1, values = c("palevioletred", "steelblue")) +
  scale_y_continuous(limits = c(-1,2),
                     breaks = c(seq(-1,2,0.5)))

################################################################################################
### Figure 2B

## Non linear relationships of temperature and malaria incidence
load("models/mod5_nl_pf.R")
load("models/mod5_nl_pv.R")

### Create df
col1 <- "grey40"
tcol1 <- do.call(rgb,c(as.list(col2rgb(col1)), alpha = 255/4, max = 255))

### Unscale temperature values to plot by adding back mean and times by sd
data <- read.csv("data/inla_input/data.csv")
# Climate data is repeated for both parasites
unscale <- data %>% subset(parasite == "Falciparum") %>%
  dplyr::summarise(sd   = sd(tmin_lag3),
                   mean = mean(tmin_lag3))

nl_df <- rbind(mod5_nl_pf$summary.random$`inla.group(tmin)` %>%
                 mutate(parasite = "P. falciparum",
                        sd_value       = unscale$sd,
                        mean_value     = unscale$mean) %>%
                 mutate(ID = ID*sd_value + mean_value) %>%
                 dplyr::select(ID, mean, `0.025quant`, `0.975quant`, parasite) %>%
                 dplyr::rename(lci = `0.025quant`,
                               uci = `0.975quant`),
               mod5_pv$summary.random$`inla.group(tmin)` %>%
                 mutate(parasite = "P. vivax",
                        sd_value       = unscale$sd,
                        mean_value     = unscale$mean) %>%
                 mutate(ID = ID*sd_value + mean_value) %>%
                 dplyr::select(ID, mean, `0.025quant`, `0.975quant`, parasite) %>%
                 dplyr::rename(lci = `0.025quant`,
                               uci = `0.975quant`))

nl_plot <- 
  ggplot(nl_df) + 
  geom_ribbon(aes(ymin = nl_df$lci, ymax = nl_df$uci, x = nl_df$ID, fill = ""), alpha = 0.2) +
  geom_hline(yintercept = 0, colour = "grey50", size = 0.3,
             linetype = "dashed") +
  geom_line(aes(x = ID, y = mean, colour = "mean")) +
  theme_classic() +
  scale_colour_manual("", values="black") +
  scale_fill_manual("", values=tcol1) +
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
  facet_wrap(~parasite, ncol = 1)

## Combine 
tiff("figures/paper/figure_2.tif", width = 180, height = 90, units = "mm", res = 520, compression = "lzw")

plot_grid(estimates_plot,
          nl_plot,
          labels = c('A)', 'B)'),
          ncol = 2,
          rel_widths = c(4, 2))

dev.off()

################################################################################################
### Figure 3

#### Look at seasonality in malaria incidence that is explained by temperature
# models with tmin
load("models/mod6_pf.R")
load("models/mod6_pv.R")

# models without tmin
load("models/mod6_wtmin_pf.R")
load("models/mod6_wtmin_pv.R")

t1_df <- rbind(rbind(mod6_pf$summary.random$t1 %>%
                       mutate(parasite = "P. falciparum",
                              month    = month.abb,
                              model    = "with Tmin") %>%
                       dplyr::select(month, mean, `0.025quant`, `0.975quant`, parasite, model) %>%
                       dplyr::rename(lci = `0.025quant`,
                                     uci = `0.975quant`),
                     mod6_pv$summary.random$t1 %>%
                       mutate(parasite = "P. vivax",
                              month    = month.abb,
                              model    = "with Tmin") %>%
                       dplyr::select(month, mean, `0.025quant`, `0.975quant`, parasite, model) %>%
                       dplyr::rename(lci = `0.025quant`,
                                     uci = `0.975quant`)),
               rbind(mod6_wtmin_pf$summary.random$t1 %>%
                       mutate(parasite = "P. falciparum",
                              month    = month.abb,
                              model    = "without Tmin") %>%
                       dplyr::select(month, mean, `0.025quant`, `0.975quant`, parasite, model) %>%
                       dplyr::rename(lci = `0.025quant`,
                                     uci = `0.975quant`),
                     mod6_wtmin_pv$summary.random$t1 %>%
                       mutate(parasite = "P. vivax",
                              month    = month.abb,
                              model    = "without Tmin") %>%
                       dplyr::select(month, mean, `0.025quant`, `0.975quant`, parasite, model) %>%
                       dplyr::rename(lci = `0.025quant`,
                                     uci = `0.975quant`)))

## Relabel for plotting
t1_df$parasite <- as.factor(t1_df$parasite)
levels(t1_df$parasite) = c("P. falciparum" = expression(paste(bold("A) "), bolditalic("P. falciparum"))),
                           "P. vivax"     = expression(paste(bold("B) "), bolditalic("P. vivax"))))
t1_df$month <- factor(t1_df$month, levels = month.abb)

label1 <- bquote(paste("Random effects with T"["min"]))
label2 <- bquote(paste("Random effects without T"["min"]))

tiff("figures/paper/figure_3.tif", width = 180, height = 90, units = "mm", res = 520, compression = "lzw")

ggplot(t1_df, aes(month, mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = lci, ymax = uci, colour = model), position = position_dodge(width = 0.7)) +
  geom_point(aes(fill = model, colour = model), shape = 21, position = position_dodge(width = 0.7)) +
  theme_classic() +
  ylab("Relative risk") + xlab("Month") +
  theme(axis.line = element_blank(),
        panel.background = element_rect(colour = "black", fill = NA, size=0.5),
        strip.background = element_blank(),
        legend.position = c(0.87,0.90),
        legend.title = element_blank(),
        strip.text = element_text(hjust = 0, size = 11),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(size = 7.5)) +
  facet_wrap(~parasite, scales = "free", nrow = 1,
             labeller = label_parsed) +
  scale_fill_manual(values = c("salmon", "grey"),
                    labels = c(label1, label2)) +
  scale_colour_manual(values = c("salmon", "grey"),
                      labels = c(label1, label2))

dev.off()


################################################################################################################
### Figure 4

### Climate suitability for transmission, based on minimum temperatures
data <- read.csv("data/inla_input/data.csv")

suitability <- rbind(data %>% subset(tmin >= 18 & parasite == "Falciparum") %>%
                       group_by(Year) %>% tally() %>%
                       mutate(parasite = "P. falciparum",
                              # divide for number of cantons
                              n = n/14),
                     
                     data %>% subset(tmin >= 15 & parasite == "Vivax") %>%
                       group_by(Year) %>% tally() %>%
                       mutate(parasite = "P. vivax",
                              # divide for number of cantons
                              n = n/14))

tiff("figures/paper/figure_4.tif", width = 180, height = 90, units = "mm", res = 520, compression = "lzw")

ggplot(suitability, aes(x = Year, y = n, group = parasite,
                        colour = "parasite")) +
  geom_smooth(method = "lm", linetype = "dashed",
              fill = "lightgrey",
              size = 0.2, aes(colour = parasite)) +
  geom_line(aes(colour = parasite)) +
  theme_classic() +
  xlab("Time") + ylab("Number of suitable months") +
  theme(axis.text.x  = element_text(angle = 90, hjust = 1),
        axis.line = element_blank(),
        legend.position = "none",
        panel.background = element_rect(colour = "black", fill = NA, size=0.5)) +
  scale_x_continuous(breaks = c(1990:2018)) +
  scale_y_continuous(name="Number of suitable months", limits=c(5, 9.5)) +
  scale_color_manual(values = c("palevioletred","steelblue"),
                     labels = c("T>= 18°C","T>= 15°C")) +
  annotate("text", x = 1992, y = 9.2, label = "T >= 15°C", colour = "steelblue", size = 4) +
  annotate("text", x = 1992, y = 7.55, label = "T >= 18°C", colour = "palevioletred", size = 4) 

dev.off()

################################################################################################
### Figure 5

### Compare model posterior distributions with and without climate information
# Models with climate information (nb: best fitting model includes linear climate info for P. falciparum and non-linear for P. vivax)
load("models/mod5_l_pf.R")
load("models/mod5_pv.R")

# Models without climate information
load("models/mod3_pf.R") 
load("models/mod3_pv.R") 

## Read in data
data <- read.csv("data/inla_input/data.csv")

data <- data %>% mutate(fit = c(mod5_l_pf$summary.fitted.values$mean, mod5_pv$summary.fitted.values$mean),
                        lci = c(mod5_l_pf$summary.fitted.values[,3], mod5_pv$summary.fitted.values[,3]),
                        uci = c(mod5_l_pf$summary.fitted.values[,5], mod5_pv$summary.fitted.values[,5])) %>%
  # Summarise over space
  dplyr::group_by(Year, Month, parasite) %>%
  dplyr::summarise(fit = mean(fit, na.rm = TRUE),
                   observed   = mean(cases, na.rm = TRUE),
                   lci        = mean(lci, na.rm = TRUE),
                   uci        = mean(uci, na.rm = TRUE),
                   population = mean(Population, na.rm = TRUE)) %>%
  # Calculate API
  mutate(observed = (observed*1000)/population,
         fit      = (fit*1000)/population,
         lci      = (lci*1000)/population,
         uci      = (uci*1000)/population)

data$parasite <- factor(data$parasite)
levels(data$parasite)[levels(data$parasite)=="Falciparum"] <- "P. falciparum"
levels(data$parasite)[levels(data$parasite)=="Vivax"] <- "P. vivax"

data$date <- as.Date(with(data, paste(Year, Month, rep(01, nrow(data)), sep = "-")), "%Y-%m-%d")

col1 <- brewer.pal(8, "Set2")[3]
tcol1 <- do.call(rgb,c(as.list(col2rgb(col1)), alpha = 255/4, max = 255))

climate_plot <-
  
  ggplot(data, aes(x = as.Date(date), y = observed)) +
  geom_line(aes(colour = "observed"), linetype = "solid") +
  geom_line(data = data, aes(x = date, y = fit, colour = "fit"), linetype = "dashed") +
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
                      breaks = c("observed", "fit"),
                      labels = c("Observed", "Modelled")) +
  scale_linetype_manual(name = "legend", values = c("solid", "dashed"),
                        breaks = c("observed", "fit"),
                        labels = c("Observed", "Modelled")) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1,2)))) +
  theme(axis.text.x  = element_blank(),
        axis.text.y  = element_text(size = 7),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.90, 0.85),
        legend.key.width = unit(2, "line"),
        legend.key.height = unit(1, "line"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.text.x = element_text(face = "bold.italic"),
        strip.background = element_blank()) 

### Model without climate information
data <- read.csv("data/inla_input/data.csv")

data$Fit <- c(mod3_pf$summary.fitted.values$mean, mod3_pv$summary.fitted.values$mean)
data$lci <- c(mod3_pf$summary.fitted.values[,3], mod3_pv$summary.fitted.values[,3])
data$uci <- c(mod3_pf$summary.fitted.values[,5], mod3_pv$summary.fitted.values[,5])

## Summarise over space
data <- as.data.frame(data %>% group_by(Year, Month, parasite) %>%
                        dplyr::summarise(Fit = mean(Fit, na.rm = TRUE),
                                         Observed      = mean(cases, na.rm = TRUE),
                                         lci      = mean(lci, na.rm = TRUE),
                                         uci      = mean(uci, na.rm = TRUE),
                                         Population = mean(Population, na.rm = TRUE)))

## Calculate API
data <- data %>% mutate(Observed = (Observed*1000)/Population,
                        Fit = (Fit*1000)/Population,
                        lci = (lci*1000)/Population,
                        uci = (uci*1000)/Population)

levels(data$parasite)[levels(data$parasite)=="Falciparum"] <- "P. falciparum"
levels(data$parasite)[levels(data$parasite)=="Vivax"] <- "P. vivax"

data$date <- as.Date(with(data, paste(Year, Month, rep(01, nrow(data)), sep = "-")), "%Y-%m-%d")

plot_w_climate <-
  
  ggplot(data, aes(x = as.Date(date), y = Observed)) +
  geom_line(aes(colour = "Observed"), linetype = "solid") +
  geom_line(data = data, aes(x = date, y = Fit, colour = "Fit"), linetype = "dashed") +
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
        axis.text.y  = element_text(size = 7),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = 'none',
        legend.key.width = unit(2, "line"),
        legend.key.height = unit(1, "line"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.text = element_blank(),
        strip.background = element_blank()) 

tiff("figures/paper/figure_5_API.tif", height = 130, width = 180, res = 460, units = "mm", compression = "lzw")

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

dev.off()

################################################################################################
### Figure 6

### Compare parameter estimates for 1990-2018 model and intervention model 2001-2015
# 1990-2018 models
load("models/mod6_pf.R")
load("models/mod6_pv.R")

# Intervention models 2001-2015
load("models/mod_int_pf.R")
load("models/mod_int_pv.R")

estimates <- rbind(rbind(mod6_pf$summary.fixed %>%
                           dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
                           mutate(parasite = "P. falciparum",
                                  model    = "1990-2018 model",
                                  variable = rownames(mod6_pf$summary.fixed)) %>%
                           dplyr::slice(-1) %>% 
                           dplyr::rename(lci = `0.025quant`,
                                         uci = `0.975quant`),
                         mod6_pv$summary.fixed %>%
                           dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
                           mutate(parasite = "P. vivax",
                                  model    = "1990-2018 model",
                                  variable = rownames(mod6_pf$summary.fixed)) %>%
                           dplyr::slice(-1) %>% 
                           dplyr::rename(lci = `0.025quant`,
                                         uci = `0.975quant`)),
                   rbind(mod_int_pf$summary.fixed %>%
                           dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
                           mutate(parasite = "P. falciparum",
                                  model    = "2001-2015 model",
                                  variable = rownames(mod_int_pf$summary.fixed)) %>%
                           dplyr::slice(-1) %>% 
                           dplyr::rename(lci = `0.025quant`,
                                         uci = `0.975quant`),
                         mod_int_pv$summary.fixed %>%
                           dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
                           mutate(parasite = "P. vivax",
                                  model    = "2001-2015 model",
                                  variable = rownames(mod_int_pv$summary.fixed)) %>%
                           dplyr::slice(-1) %>% 
                           dplyr::rename(lci = `0.025quant`,
                                         uci = `0.975quant`)))

### Order and label for plotting
estimates$variable <- as.factor(estimates$variable)
estimates$variable <- factor(estimates$variable, 
                             levels = c("houses_fogged",
                                        "blocks_fumigated",
                                        "houses_IRS",
                                        "total_poverty", 
                                        "urban", "prcp",
                                        "tmin"),
                             labels = c("Space\nspraying",
                                        "ULV\nfumigation",
                                        "Indoor residual\nspraying",
                                        "Poverty", 
                                        "Level of\nurbanization", "Precipitation",
                                        "Minimum\ntemperature"))

col1 <- "#08174D"
col2 <- "#339989"

## Relevel 
estimates$parasite <- as.factor(estimates$parasite)
levels(estimates$parasite)= c("P. falciparum"=expression(paste(bold("A) "), bolditalic("P. falciparum"))),
                              "P. vivax"=expression(paste(bold("B) "), bolditalic("P. vivax"))))

tiff("figures/paper/figure_6.tif", width = 180, height = 90, units = "mm", res = 520, compression = "lzw")

ggplot(estimates, aes(y = estimates$mean, x = estimates$variable, colour = estimates$model)) +
  geom_hline(yintercept = 0, colour = "darkgrey", linetype = "dashed", size = 0.3) +
  geom_linerange(aes(ymin = estimates$lci, ymax = estimates$uci), size = 0.4, position = position_dodge(width = 0.7)) +
  geom_point(aes(x = estimates$variable, y = estimates$mean, colour = estimates$model), size = 3,
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
  facet_wrap(~parasite, nrow = 1, labeller = label_parsed) +
  scale_y_continuous(limits = c(-1,1.8),
                     breaks = c(seq(-1,1.8,0.5)))

dev.off()

################################################################################################
### Figure 7

### Compare the model improvement for each intervention measure
# Models with all interventions (nb: best fitting model includes linear climate info for P. falciparum and non-linear for P. vivax)
load("models/mod_int_pf.R")
load("models/mod_int_nl_pv.R")

# Models without each intervention
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

rmse_df <- rbind(mod_int_w_irs_pf$summary.fitted.values %>%
                   dplyr::select(`0.5quant`) %>%
                   mutate(observed = c(subset(data, data$parasite == "Falciparum")$cases),
                          int      = c(mod_int_pf$summary.fitted.values$`0.5quant`),
                          parasite = "P. falciparum",
                          model    = "without_irs",
                          id   = subset(data, data$parasite == "Falciparum")$Canton) %>%
                   tibble::remove_rownames() %>% 
                   dplyr::rename(fit = `0.5quant`),
                 mod_int_w_irs_pv$summary.fitted.values %>%
                   dplyr::select(`0.5quant`) %>%
                   mutate(observed = c(subset(data, data$parasite == "Vivax")$cases),
                          int      = c(mod_int_nl_pv$summary.fitted.values$`0.5quant`),
                          parasite = "P. vivax",
                          model    = "without_irs",
                          id   = subset(data, data$parasite == "Vivax")$Canton) %>%
                   tibble::remove_rownames() %>% 
                   dplyr::rename(fit = `0.5quant`),
                 mod_int_w_fog_pf$summary.fitted.values %>%
                   dplyr::select(`0.5quant`) %>%
                   mutate(observed = c(subset(data, data$parasite == "Falciparum")$cases),
                          int      = c(mod_int_pf$summary.fitted.values$`0.5quant`),
                          parasite = "P. falciparum",
                          model    = "without_fog",
                          id   = subset(data, data$parasite == "Falciparum")$Canton) %>%
                   tibble::remove_rownames() %>% 
                   dplyr::rename(fit = `0.5quant`),
                 mod_int_w_fog_pv$summary.fitted.values %>%
                   dplyr::select(`0.5quant`) %>%
                   mutate(observed = c(subset(data, data$parasite == "Vivax")$cases),
                          int      = c(mod_int_nl_pv$summary.fitted.values$`0.5quant`),
                          parasite = "P. vivax",
                          model    = "without_fog",
                          id   = subset(data, data$parasite == "Vivax")$Canton) %>%
                   tibble::remove_rownames() %>% 
                   dplyr::rename(fit = `0.5quant`),
                 mod_int_w_fum_pf$summary.fitted.values %>%
                   dplyr::select(`0.5quant`) %>%
                   mutate(observed = c(subset(data, data$parasite == "Falciparum")$cases),
                          int      = c(mod_int_pf$summary.fitted.values$`0.5quant`),
                          parasite = "P. falciparum",
                          model    = "without_fum",
                          id   = subset(data, data$parasite == "Falciparum")$Canton) %>%
                   tibble::remove_rownames() %>% 
                   dplyr::rename(fit = `0.5quant`),
                 mod_int_w_fum_pv$summary.fitted.values %>%
                   dplyr::select(`0.5quant`) %>%
                   mutate(observed = c(subset(data, data$parasite == "Vivax")$cases),
                          int      = c(mod_int_nl_pv$summary.fitted.values$`0.5quant`),
                          parasite = "P. vivax",
                          model    = "without_fum",
                          id   = subset(data, data$parasite == "Vivax")$Canton) %>%
                   tibble::remove_rownames() %>% 
                   dplyr::rename(fit = `0.5quant`))

# Calculate rmse
rmse_df <- rmse_df %>% mutate(e     = (observed-fit)^2,
                              e_int = (observed-int)^2) %>%
  dplyr::group_by(id, parasite, model) %>%
  dplyr::summarise(rmse     = sqrt(mean(e, na.rm = TRUE)),
                   rmse_int = sqrt(mean(e_int, na.rm = TRUE))) %>%
  mutate(rmse_difference = (rmse - rmse_int)/rmse * 100) %>%
  dplyr::select(id, parasite, model, rmse_difference) 

# No data for canton so replace with NA
rmse_df$rmse_difference[rmse_df$id == "Chilla"] <- NA

# Re order models for plotting
rmse_df$model <- factor(rmse_df$model, levels = c("without_irs", 
                                                  "without_fum", 
                                                  "without_fog"),
                        labels = c("Indoor residual spraying",
                                   "ULV fumigation",
                                   "Space spraying"))

## Merge data with shapefile to plot
ecuador <- getData('GADM', country = "ECU", level = 2)
el_oro <- ecuador %>% subset(NAME_1 =="El Oro") %>%
  fortify(region = "NAME_2")
el_oro_rmse <- merge(el_oro, rmse_df, by = "id")

## Center colour palette around 0 
limit <- max(abs(el_oro_rmse$rmse_difference)) * c(-1, 1)

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
                       midpoint = 0, limits = c(-14, max(el_oro_rmse$rmse_difference, na.rm =T)),
                       name = "Model\nimprovement\n(%)",
                       na.value = "grey43",
                       guide = guide_colourbar(ticks = FALSE, barheight = 5,
                                               barwidth = 1))

tiff("figures/paper/figure_7.tif", width = 180, height = 80, units = "mm", res = 520, compression = "lzw")

ggdraw(align_legend(p))

dev.off()


################################################################################################
### Figure 8

### Compare interannual random effects of 1990-2018 model with 2001-2015 model
# 1990-2018 models (nb: best fitting model includes linear climate info for P. falciparum and non-linear for P. vivax)
load("models/mod5_l_pf.R")
load("models/mod5_pv.R")

# 2001-2015 model
load("models/mod_int_pf.R")
load("models/mod_int_nl_pv.R")

t2_df <- rbind(rbind(mod5_l_pf$summary.random$t2 %>%
                       mutate(parasite = "P. falciparum",
                              model    = "Random effects of 1990-2018 model") %>%
                       dplyr::select(ID, mean, `0.025quant`, `0.975quant`, parasite, model) %>%
                       dplyr::rename(year = ID,
                                     lci = `0.025quant`,
                                     uci = `0.975quant`),
                     mod5_pv$summary.random$t2 %>%
                       mutate(parasite = "P. vivax",
                              model    = "Random effects of 1990-2018 model") %>%
                       dplyr::select(ID, mean, `0.025quant`, `0.975quant`, parasite, model) %>%
                       dplyr::rename(year = ID,
                                     lci = `0.025quant`,
                                     uci = `0.975quant`)),
               rbind(mod_int_pf$summary.random$t2 %>%
                       mutate(parasite = "P. falciparum",
                              model    = "Random effects of 2001-2015 model") %>%
                       dplyr::select(ID, mean, `0.025quant`, `0.975quant`, parasite, model) %>%
                       dplyr::rename(year = ID,
                                     lci = `0.025quant`,
                                     uci = `0.975quant`),
                     mod_int_nl_pv$summary.random$t2 %>%
                       mutate(parasite = "P. vivax",
                              model    = "Random effects of 2001-2015 model") %>%
                       dplyr::select(ID, mean, `0.025quant`, `0.975quant`, parasite, model) %>%
                       dplyr::rename(year = ID,
                                     lci = `0.025quant`,
                                     uci = `0.975quant`)))

col1 <- "#339989"

# Labels
t2_df$parasite <- as.factor(t2_df$parasite)
levels(t2_df$parasite)= c("P. falciparum"=expression(paste(bold("A) "), bolditalic("P. falciparum"))),
                          "P. vivax"=expression(paste(bold("B) "), bolditalic("P. vivax"))))

tiff("figures/paper/figure_8.tif", width = 180, height = 90, res = 520, units = "mm", compression = "lzw")

ggplot(t2_df, aes(year, mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin=lci, ymax=uci, colour = model), position = position_dodge(width = 0.7)) +
  geom_point(aes(fill = model, colour = model), shape = 21, position = position_dodge(width = 0.7)) +
  ylab("Relative risk") + xlab("Year") +
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
  facet_grid(~parasite, scales = "free", labeller = label_parsed) +
  scale_fill_manual(values = c("grey", col1)) +
  scale_colour_manual(values = c("grey", col1)) 

dev.off()
