########################################################################################################################################

### Bayesian hierarchical models for modelling the effect of interventions on malaria incidence in El Oro

########################################################################################################################################

###########################################################################################################################################

## Model for intervention period only, for when we have data for 2001-2015

##########################################################################################################################################
## Load libraries
pacman::p_load("raster", "INLA", "sf", "sp", "spdep", "dplyr", "ggplot2",
               "scales", "knitr", "kableExtra", "doBy", "ggpubr",
               "RColorBrewer", "Hmisc", "cowplot", "bookdown", "raster",
               "lattice", "rasterVis", "brinla", "miceadds", "gridExtra",
               "ggpubr", "gganimate", "gifski", "transformr", "tidyr",
               "ggregplot", "grid", "MASS", "INLAutils", "reshape2", "zoo")

### Define neighbourhood matrix
ecuador <- getData('GADM', country = "ECU", level = 2)
el_oro  <- subset(ecuador, NAME_1 == "El Oro")

## Create neighbourhood matrix
nb.map <- poly2nb(el_oro)
nb2INLA("map.graph",nb.map)

## Read in data
data <- read.csv("data/data.csv", fileEncoding = "latin1")

## Subset data to intervention period
data <- subset(data, data$Year < 2016 & data$Year > 2000)

## Replace NA values in intervention data to 0
data$houses_fogged[is.na(data$houses_fogged)] <- 0
data$houses_IRS[is.na(data$houses_IRS)]       <- 0
data$blocks_fumigated[is.na(data$blocks_fumigated)]       <- 0

## Time lags
# nb: lag in multiples of 14 for the cantons
data <- data %>% mutate(
  blocks_fumigated_lag1 = lag(blocks_fumigated, 14),
  blocks_fumigated_lag2 = lag(blocks_fumigated, 28),
  blocks_fumigated_lag3 = lag(blocks_fumigated, 42),
  blocks_fumigated_lag4 = lag(blocks_fumigated, 56),
  
  houses_fogged_lag1 = lag(houses_fogged, 14),
  houses_fogged_lag2 = lag(houses_fogged, 28),
  houses_fogged_lag3 = lag(houses_fogged, 42),
  houses_fogged_lag4 = lag(houses_fogged, 56),
  
  houses_IRS_lag1 = lag(houses_IRS, 14),
  houses_IRS_lag2 = lag(houses_IRS, 28),
  houses_IRS_lag3 = lag(houses_IRS, 42),
  houses_IRS_lag4 = lag(houses_IRS, 56))

## Add intervention period before/after
data$int_per <- 0
data$int_per[data$Year > 2000] <- 1

# Separate parasite models
data_pf <- subset(data, data$parasite == "Falciparum")
data_pv <- subset(data, data$parasite == "Vivax")

#######################################################################################################################################################

#### Falciparum models

## Fixed effects
y  <- data_pf$cases
n  <- length(y)

## Add API so modify offset to include it 
e  <- (data_pf$Population/12)/1000

## Climate effects, scaling the covariates
prcp <- scale(data_pf$prcp_lag3, center = TRUE, scale = TRUE)[,1] 
tmin <- scale(data_pf$tmin_lag3, center = TRUE, scale = TRUE)[,1]              

# Land cover variables   
urban <- scale(data_pf$urban, center = TRUE, scale = TRUE)[,1]     

## Socioeconomic data
total_poverty <- scale(data_pf$total_poverty, center = TRUE, scale = TRUE)[,1] 

## Intervention data
blocks_fumigated <- scale(data_pf$blocks_fumigated_lag3, center = TRUE, scale = TRUE)[,1] 
houses_IRS    <- scale(data_pf$houses_IRS_lag3, center = TRUE, scale = TRUE)[,1] 
houses_fogged <- scale(data_pf$houses_fogged_lag2, center = TRUE, scale = TRUE)[,1] 

### Random effects
# Temporal
t1 <- as.factor(data_pf$Month) # Seasonality
t2 <- as.factor(data_pf$Year)  # Interannual 

# In order to include structured and unstructured effects separately in the hierarchical model they need to be specified separately, with different latent models. 
# Structured spatial effects (model = besag)
s1 <- rep(1:14, 180) 

# Unstructured spatial effects (model = iid)
s2 <- rep(1:14, 180) 


df_inla_pf <- data.frame(y, e, prcp, 
                         tmin, urban, 
                         total_poverty, blocks_fumigated, houses_fogged,
                         houses_IRS, t1, t2, s1, s2)


######################################################################
######### Model with interventions 2001-2015
formula <- y ~ 1 + f(s1, model = "iid", graph = "map.graph") +
  f(s2, model = "besag", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  tmin +
  prcp +
  urban +
  total_poverty +
  houses_fogged + houses_IRS + blocks_fumigated 


mod_int_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                    offset = log(e), verbose = TRUE,
                    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                           config = FALSE, 
                                           return.marginals = FALSE), 
                    control.predictor = list(link = 1, compute = TRUE), 
                    control.family = list(link = "log"))

save(mod_int_pf, file = "models/int_mods/mod_int_pf.R")


formula <- y ~ 1 + f(s1, model = "iid", graph = "map.graph") +
  f(s2, model = "besag", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(tmin), model = "rw1") +
  f(inla.group(prcp), model = "rw1") +
  urban +
  total_poverty + 
  houses_fogged + houses_IRS + blocks_fumigated 


mod_int_nl_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                       offset = log(e), verbose = TRUE,
                       control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                              config = FALSE, 
                                              return.marginals = FALSE), 
                       control.predictor = list(link = 1, compute = TRUE), 
                       control.family = list(link = "log"))

save(mod_int_nl_pf, file = "models/int_mods/mod_int_nl_pf.R")


formula <- y ~ 1 + f(s1, model = "iid", graph = "map.graph") +
  f(s2, model = "besag", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(tmin), model = "rw1") +
  f(inla.group(prcp), model = "rw1") +
  urban +
  total_poverty  +
  houses_fogged + 
  #houses_IRS 
  blocks_fumigated 


mod_int_w_irs_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                          offset = log(e), verbose = TRUE,
                          control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                 config = FALSE, 
                                                 return.marginals = FALSE), 
                          control.predictor = list(link = 1, compute = TRUE), 
                          control.family = list(link = "log"))

save(mod_int_w_irs_pf, file = "models/int_mods/mod_int_w_irs_pf.R")

formula <- y ~ 1 + f(s1, model = "iid", graph = "map.graph") +
  f(s2, model = "besag", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(tmin), model = "rw1") +
  f(inla.group(prcp), model = "rw1") +
  urban +
  total_poverty  +
  houses_fogged + 
  houses_IRS 
#blocks_fumigated 


mod_int_w_fum_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                          offset = log(e), verbose = TRUE,
                          control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                 config = FALSE,
                                                 return.marginals = FALSE), 
                          control.predictor = list(link = 1, compute = TRUE), 
                          control.family = list(link = "log"))

save(mod_int_w_fum_pf, file = "models/int_mods/mod_int_w_fum_pf.R")

formula <- y ~ 1 + f(s1, model = "iid", graph = "map.graph") +
  f(s2, model = "besag", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(tmin), model = "rw1") +
  f(inla.group(prcp), model = "rw1") +
  urban +
  total_poverty  +
  #houses_fogged 
  houses_IRS + 
  blocks_fumigated 


mod_int_w_fog_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                          offset = log(e), verbose = TRUE,
                          control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                 config = FALSE, 
                                                 return.marginals = FALSE), 
                          control.predictor = list(link = 1, compute = TRUE), 
                          control.family = list(link = "log"))

save(mod_int_w_fog_pf, file = "models/int_mods/mod_int_w_fog_pf.R")


####### Model without interventions to compare 2001-2015
formula <- y ~ 1 + f(s1, model = "iid", graph = "map.graph") +
  f(s2, model = "besag", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(tmin), model = "rw1") +
  f(inla.group(prcp), model = "rw1") +
  urban +
  total_poverty 

mod_int3_nl_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                       offset = log(e), verbose = TRUE,
                       control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                              config = FALSE, 
                                              return.marginals = FALSE), 
                       control.predictor = list(link = 1, compute = TRUE), 
                       control.family = list(link = "log"))

save(mod_int3_nl_pf, file = "models/int_mods/mod_int3_nl_pf.R")


#######################################################################################################################################################

#### Vivax models

## Fixed effects
y  <- data_pv$cases
n  <- length(y)

## Add API so modify offset to include it 
e  <- (data_pv$Population/12)/1000

## Climate effects, scaling the covariates
prcp <- scale(data_pv$prcp_lag2, center = TRUE, scale = TRUE)[,1]  
tmin <- scale(data_pv$tmin_lag3, center = TRUE, scale = TRUE)[,1]              

# Land cover variables
ndvi <- scale(data_pv$ndvi, center = TRUE, scale = TRUE)[,1]     

## Socioeconomic data
total_poverty <- scale(data_pv$total_poverty, center = TRUE, scale = TRUE)[,1] 

## Intervention data
blocks_fumigated <- scale(data_pv$blocks_fumigated_lag3, center = TRUE, scale = TRUE)[,1] 
houses_IRS    <- scale(data_pv$houses_IRS_lag3, center = TRUE, scale = TRUE)[,1] 
houses_fogged <- scale(data_pv$houses_fogged_lag3, center = TRUE, scale = TRUE)[,1] 

### Random effects
# Temporal
t1 <- as.factor(data_pv$Month) # Seasonality
t2 <- as.factor(data_pv$Year)  # Interannual 

# In order to include structured and unstructured effects separately in the hierarchical model they need to be specified separately, with different latent models. 
# Structured spatial effects (model = besag)
s1 <- rep(1:14, 180) 

# Unstructured spatial effects (model = iid)
s2 <- rep(1:14, 180) 


df_inla_pv <- data.frame(y, e, prcp,  
                         tmin, urban, 
                         total_poverty, blocks_fumigated, houses_fogged,
                         houses_IRS, t1, t2, s1, s2)

######################################################################
######### Model with interventions 2001-2015
formula <- y ~ 1 + f(s1, model = "iid", graph = "map.graph") +
  f(s2, model = "besag", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  tmin +
  prcp +
  urban +
  total_poverty +
  houses_fogged + houses_IRS + blocks_fumigated  


mod_int_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                    offset = log(e), verbose = TRUE,
                    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                           config = FALSE, 
                                           return.marginals = FALSE), 
                    control.predictor = list(link = 1, compute = TRUE), 
                    control.family = list(link = "log"))

save(mod_int_pv, file = "models/int_mods/mod_int_pv.R")


formula <- y ~ 1 + f(s1, model = "iid", graph = "map.graph") +
  f(s2, model = "besag", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(tmin), model = "rw1") +
  f(inla.group(prcp), model = "rw1") +
  urban +
  total_poverty  +
  houses_fogged + houses_IRS + blocks_fumigated 


mod_int_nl_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                       offset = log(e), verbose = TRUE,
                       control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                              config = FALSE, 
                                              return.marginals = FALSE), 
                       control.predictor = list(link = 1, compute = TRUE), 
                       control.family = list(link = "log"))

save(mod_int_nl_pv, file = "models/int_mods/mod_int_nl_pv.R")

formula <- y ~ 1 + f(s1, model = "iid", graph = "map.graph") +
  f(s2, model = "besag", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(tmin), model = "rw1") +
  f(inla.group(prcp), model = "rw1") +
  urban +
  total_poverty  +
  houses_fogged + 
  #houses_IRS 
  blocks_fumigated 


mod_int_w_irs_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                          offset = log(e), verbose = TRUE,
                          control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                 config = FALSE, 
                                                 return.marginals = FALSE), 
                          control.predictor = list(link = 1, compute = TRUE), 
                          control.family = list(link = "log"))

save(mod_int_w_irs_pv, file = "models/int_mods/mod_int_w_irs_pv.R")

formula <- y ~ 1 + f(s1, model = "iid", graph = "map.graph") +
  f(s2, model = "besag", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(tmin), model = "rw1") +
  f(inla.group(prcp), model = "rw1") +
  urban +
  total_poverty  +
  houses_fogged + 
  houses_IRS 
#blocks_fumigated 


mod_int_w_fum_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                          offset = log(e), verbose = TRUE,
                          control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                 config = FALSE, 
                                                 return.marginals = FALSE), 
                          control.predictor = list(link = 1, compute = TRUE), 
                          control.family = list(link = "log"))

save(mod_int_w_fum_pv, file = "models/int_mods/mod_int_w_fum_pv.R")

formula <- y ~ 1 + f(s1, model = "iid", graph = "map.graph") +
  f(s2, model = "besag", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(tmin), model = "rw1") +
  f(inla.group(prcp), model = "rw1") +
  urban +
  total_poverty  +
  #houses_fogged 
  houses_IRS + 
  blocks_fumigated 


mod_int_w_fog_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                          offset = log(e), verbose = TRUE,
                          control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                 config = FALSE, 
                                                 return.marginals = FALSE), 
                          control.predictor = list(link = 1, compute = TRUE), 
                          control.family = list(link = "log"))

save(mod_int_w_fog_pv, file = "models/int_mods/mod_int_w_fog_pv.R")



