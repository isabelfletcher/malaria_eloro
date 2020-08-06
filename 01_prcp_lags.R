####################### Precipitation lags #######################

## Here we test for the best monthly time lags for precipitation in spatiotemporal models of malaria incidence

## Load libraries
pacman::p_load("raster", "INLA","dplyr", 
               "kableExtra", "reshape2",
               "spdep", "tidyr")

### Define neighbourhood matrix
ecuador <- getData('GADM', country = "ECU", level = 2)
el_oro  <- subset(ecuador, NAME_1 == "El Oro")

## Create neighbourhood matrix
nb.map <- poly2nb(el_oro)
nb2INLA("map.graph",nb.map)

## Read in data
data <- read.csv("data.csv")
data$cases <- as.numeric(data$cases)

## Add intervention period (2001-2015) before/after
data$int_per <- 0
data$int_per[data$Year > 2000] <- 1

# Separate parasite models
data_pf <- subset(data, data$parasite == "Falciparum")
data_pv <- subset(data, data$parasite == "Vivax")


#######################################################################################################################################################

#### Falciparum models
int_per <- as.factor(data_pf$int_per)

## Fixed effects
y  <- data_pf$cases
n  <- length(y)

## Add API so modify offset to include it 
e  <- (data_pf$Population/12)/1000

## Climate effects, scaling the covariates
prcp <- scale(data_pf$prcp, center = TRUE, scale = TRUE)[,1] 
tmin <- scale(data_pf$tmin, center = TRUE, scale = TRUE)[,1]

urban <- scale(data_pf$urban, center = TRUE, scale = TRUE)[,1]

## Socioeconomic data
total_poverty <- scale(data_pf$total_poverty, center = TRUE, scale = TRUE)[,1]

### Random effects
# Temporal
t1 <- as.factor(data_pf$Month) # Seasonality
t2 <- as.factor(data_pf$Year)  # Interannual 

# Spatial effects: in order to determine the variation due to spatial autocorrelation (structured effects), we specify the besag model
# We specify the iid distribution to take into account unstructured spatial variation (heteorogeneity)
s1 <- rep(1:14, 348) # there are 14 cantons/districts
s2 <- rep(1:14, 348) 


##################################
# Add time lags and scale variables
prcp_lag0 <- scale(data_pf$prcp, center = TRUE, scale = TRUE)[,1]
prcp_lag1 <- scale(data_pf$prcp_lag1, center = TRUE, scale = TRUE)[,1]
prcp_lag2 <- scale(data_pf$prcp_lag2, center = TRUE, scale = TRUE)[,1]
prcp_lag3 <- scale(data_pf$prcp_lag3, center = TRUE, scale = TRUE)[,1]

df_inla_pf <- data.frame(y, tmin, prcp_lag0, prcp_lag1,
                         prcp_lag2, prcp_lag3, 
                         t1, t2, s1, s2, 
                         total_poverty, urban, int_per)


lags <- c("prcp_lag0", "prcp_lag1", "prcp_lag2", "prcp_lag3")

###########################################################################################################

# Loop through each lag
prcp_lags_l <- NULL
prcp_lags_nl <- NULL

for (i in 1:4) {
  
  lag <- extract_numeric(lags[i])
  
  # Get lagged data and scale
  prcp <- scale(df_inla_pf[,grep(lags[i], names(df_inla_pf), value=TRUE)], 
                center = TRUE, scale = TRUE)[,1]
  
  formula_l <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
    f(s2, model = "iid", graph = "map.graph") +
    f(t1, model = "rw1") +
    f(t2, model = "iid") +
    prcp + 
    tmin +
    urban +
    ## Add interaction
    int_per + 
    urban*int_per +
    total_poverty
  
  mod_l <- inla(formula_l, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = FALSE, 
                                       return.marginals = FALSE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))
  
  formula_nl <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
    f(s2, model = "iid", graph = "map.graph") +
    f(t1, model = "rw1") +
    f(t2, model = "iid") +
    f(inla.group(prcp), model = "rw1") + 
    f(inla.group(tmin), model = "rw1") +
    urban +
    ## Add interaction
    int_per + 
    urban*int_per +
    total_poverty
  
  mod_nl <- inla(formula_nl, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                 offset = log(e), verbose = TRUE,
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                        config = FALSE, 
                                        return.marginals = FALSE), 
                 control.predictor = list(link = 1, compute = TRUE), 
                 control.family = list(link = "log"))
  
  ## Extract model DIC, WAIC and model estimates
  df_l <- data.frame(Parasite  = "P. falciparum", 
                     Lag       = paste0(lag),
                     Estimate  = mod_l$summary.fixed$mean[2],
                     LCI       = mod_l$summary.fixed$`0.025quant`[2],
                     UCI       = mod_l$summary.fixed$`0.975quant`[2],
                     DIC       = mod_l$dic$dic,
                     WAIC      = mod_l$waic$waic)
  
  df_nl <- data.frame(Parasite  = "P. falciparum", 
                      Lag       = paste0(lag),
                      DIC       = mod_nl$dic$dic,
                      WAIC      = mod_nl$waic$waic)
  
  prcp_lags_l  <- rbind(prcp_lags_l, df_l)
  prcp_lags_nl <- rbind(prcp_lags_nl, df_nl)
  
}

data_l  <- prcp_lags_l
data_nl <- prcp_lags_nl

#######################################################################################################################################################

#### Vivax models
int_per <- as.factor(data_pv$int_per)

## Fixed effects
y  <- data_pv$cases
n  <- length(y)

## Add API so modify offset to include it 
e  <- (data_pv$Population/12)/1000

## Climate effects, scaling the covariates
prcp <- scale(data_pv$prcp, center = TRUE, scale = TRUE)[,1] 
tmin <- scale(data_pv$tmin, center = TRUE, scale = TRUE)[,1]

urban <- scale(data_pv$urban, center = TRUE, scale = TRUE)[,1]

## Socioeconomic data
total_poverty <- scale(data_pv$total_poverty, center = TRUE, scale = TRUE)[,1]

### Random effects
# Temporal
t1 <- as.factor(data_pv$Month) # Seasonality
t2 <- as.factor(data_pv$Year)  # Interannual 

# Spatial effects: in order to determine the variation due to spatial autocorrelation (structured effects), we specify the besag model
# We specify the iid distribution to take into account unstructured spatial variation (heteorogeneity)
s1 <- rep(1:14, 348) # there are 14 cantons/districts
s2 <- rep(1:14, 348) 


##################################
# Add time lags and scale variables
prcp_lag0 <- scale(data_pv$prcp, center = TRUE, scale = TRUE)[,1]
prcp_lag1 <- scale(data_pv$prcp_lag1, center = TRUE, scale = TRUE)[,1]
prcp_lag2 <- scale(data_pv$prcp_lag2, center = TRUE, scale = TRUE)[,1]
prcp_lag3 <- scale(data_pv$prcp_lag3, center = TRUE, scale = TRUE)[,1]

df_inla_pv <- data.frame(y, tmin, prcp_lag0, prcp_lag1,
                         prcp_lag2, prcp_lag3, 
                         t1, t2, s1, s2, 
                         total_poverty, urban, int_per)


lags <- c("prcp_lag0", "prcp_lag1", "prcp_lag2", "prcp_lag3")

###########################################################################################################

# Loop through each lag
prcp_lags_l <- NULL
prcp_lags_nl <- NULL

for (i in 1:4) {
  
  lag <- extract_numeric(lags[i])
  
  # Get lagged data and scale
  prcp <- scale(df_inla_pv[,grep(lags[i], names(df_inla_pv), value=TRUE)], 
                center = TRUE, scale = TRUE)[,1]
  
  formula_l <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
    f(s2, model = "iid", graph = "map.graph") +
    f(t1, model = "rw1") +
    f(t2, model = "iid") +
    prcp + 
    tmin +
    urban +
    ## Add interaction
    int_per + 
    urban*int_per +
    total_poverty
  
  mod_l <- inla(formula_l, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = FALSE, 
                                       return.marginals = FALSE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))
  
  formula_nl <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
    f(s2, model = "iid", graph = "map.graph") +
    f(t1, model = "rw1") +
    f(t2, model = "iid") +
    f(inla.group(prcp), model = "rw1") + 
    f(inla.group(tmin), model = "rw1") +
    urban +
    ## Add interaction
    int_per + 
    urban*int_per +
    total_poverty
  
  mod_nl <- inla(formula_nl, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                 offset = log(e), verbose = TRUE,
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                        config = FALSE, 
                                        return.marginals = FALSE), 
                 control.predictor = list(link = 1, compute = TRUE), 
                 control.family = list(link = "log"))
  
  ## Extract model DIC, WAIC and model estimates
  df_l <- data.frame(Parasite  = "P. vivax", 
                     Lag       = paste0(lag),
                     Estimate  = mod_l$summary.fixed$mean[2],
                     LCI       = mod_l$summary.fixed$`0.025quant`[2],
                     UCI       = mod_l$summary.fixed$`0.975quant`[2],
                     DIC       = mod_l$dic$dic,
                     WAIC      = mod_l$waic$waic)
  
  df_nl <- data.frame(Parasite  = "P. vivax", 
                      Lag       = paste0(lag),
                      DIC       = mod_nl$dic$dic,
                      WAIC      = mod_nl$waic$waic)
  
  prcp_lags_l  <- rbind(prcp_lags_l, df_l)
  prcp_lags_nl <- rbind(prcp_lags_nl, df_nl)
  
}

## Combine and write output
data_l  <- rbind(data_l, prcp_lags_l)
data_nl <- rbind(data_nl, prcp_lags_nl)

write.csv(data_l,  file = "supplementary/table_s2_prcp_lags_l.csv")
write.csv(data_nl, file = "supplementary/table_s3_prcp_lags_nl.csv")

