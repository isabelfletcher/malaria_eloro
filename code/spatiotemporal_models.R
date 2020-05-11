########################################################################################################################################

### Bayesian hierarchical models for modelling the effect of climate variation and interventions on malaria incidence in El Oro

########################################################################################################################################

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
data <- read.csv("data/inla_input/data.csv")

## Add intervention period before/after
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
prcp <- scale(data_pf$prcp_lag3, center = TRUE, scale = TRUE)[,1]             
tmin <- scale(data_pf$tmin_lag3, center = TRUE, scale = TRUE)[,1] 

urban <- scale(data_pf$urban, center = TRUE, scale = TRUE)[,1]     

## Socioeconomic data
total_poverty <- scale(data_pf$total_poverty, center = TRUE, scale = TRUE)[,1] 

### Random effects
# Temporal
t1 <- as.factor(data_pf$Month) # Seasonality
t2 <- as.factor(data_pf$Year)  # Interannual 

# Spatial
s1 <- rep(1:14, 348) 

df_inla_pf <- data.frame(y, e, prcp, tmin,
                         int_per, urban,  
                         total_poverty, 
                         t1, t2, s1)

########################################################################################

## Baseline model - spatiotemporal effects 

########################################################################################

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") 


mod1_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = FALSE, 
                                       return.marginals = FALSE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))

save(mod1_pf, file = "models/mod1_pf.R")

######################################################################################

## Add t2 random effects

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid")


mod1_2_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                  offset = log(e), verbose = TRUE,
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                         config = FALSE, 
                                         return.marginals = FALSE), 
                  control.predictor = list(link = 1, compute = TRUE), 
                  control.family = list(link = "log"))

save(mod1_2_pf, file = "models/mod1_2_pf.R")

########################################################################################

## Socioeconomic information

########################################################################################

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty 

mod2_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = FALSE, 
                                       return.marginals = FALSE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))

save(mod2_pf, file = "models/mod2_pf.R")

########################################################################################

## Add urban

########################################################################################

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  ## Add interaction
  int_per +
  urban*int_per


mod3_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = FALSE, 
                                       return.marginals = FALSE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))

save(mod3_pf, file = "models/mod3_pf.R")

########################################################################################

## Add temperature

########################################################################################

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  f(inla.group(tmin), model = "rw1") 


mod4_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = FALSE, 
                                       return.marginals = FALSE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))

save(mod4_pf, file = "models/mod4_pf.R")

## Test linear
formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  ## Add interaction with intervention period
  int_per + 
  urban*int_per +
  tmin



mod4_l_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                  offset = log(e), verbose = TRUE,
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                         config = FALSE, 
                                         return.marginals = FALSE), 
                  control.predictor = list(link = 1, compute = TRUE), 
                  control.family = list(link = "log"))

save(mod4_l_pf, file = "models/mod4_l_pf.R")

########################################################################################

## Add precipitation

########################################################################################

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  f(inla.group(tmin), model = "rw1") +
  f(inla.group(prcp), model = "rw1") 

mod5_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = TRUE, 
                                       return.marginals = TRUE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))

save(mod5_pf, file = "models/mod5_pf.R")

## Test linear
formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  ## Add interaction with intervention period
  int_per + 
  urban*int_per +
  tmin +
  prcp


mod5_l_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                  offset = log(e), verbose = TRUE,
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                         config = FALSE, 
                                         return.marginals = FALSE), 
                  control.predictor = list(link = 1, compute = TRUE), 
                  control.family = list(link = "log"))

save(mod5_l_pf, file = "models/mod5_l_pf.R")

### With model without interaction, to plot parameter estimates
formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  tmin +
  prcp


mod6_l_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                  offset = log(e), verbose = TRUE,
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                         config = FALSE, 
                                         return.marginals = FALSE), 
                  control.predictor = list(link = 1, compute = TRUE), 
                  control.family = list(link = "log"))

save(mod6_l_pf, file = "models/mod6_l_pf.R")

########################################################################################

#### Test full model without tmin 

########################################################################################

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  f(inla.group(prcp), model = "rw1") 

mod6_wtmin_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                      offset = log(e), verbose = TRUE,
                      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                             config = FALSE, 
                                             return.marginals = FALSE), 
                      control.predictor = list(link = 1, compute = TRUE), 
                      control.family = list(link = "log"))

save(mod6_wtmin_pf, file = "models/mod6_wtmin_pf.R")

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  prcp

mod6_l_wtmin_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                        offset = log(e), verbose = TRUE,
                        control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                               config = FALSE, 
                                               return.marginals = FALSE), 
                        control.predictor = list(link = 1, compute = TRUE), 
                        control.family = list(link = "log"))

save(mod6_l_wtmin_pf, file = "models/mod6_l_wtmin_pf.R")

########################################################################################

#### Test full model without prcp

########################################################################################

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  f(inla.group(tmin), model = "rw1") 

mod6_wprcp_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                      offset = log(e), verbose = TRUE,
                      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                             config = FALSE, 
                                             return.marginals = FALSE), 
                      control.predictor = list(link = 1, compute = TRUE), 
                      control.family = list(link = "log"))

save(mod6_wprcp_pf, file = "models/mod6_wprcp_pf.R")

#######################################################################################################################################################

#### Vivax models
int_per <- as.factor(data_pv$int_per)

## Fixed effects
y  <- data_pv$cases
n  <- length(y)

## Add API so modify offset to include it 
e  <- (data_pv$Population/12)/1000

## Climate effects, scaling the covariates
prcp <- scale(data_pv$prcp_lag2, center = TRUE, scale = TRUE)[,1]             
tmin <- scale(data_pv$tmin_lag3, center = TRUE, scale = TRUE)[,1] 

urban <- scale(data_pv$urban, center = TRUE, scale = TRUE)[,1]     

## Socioeconomic data
total_poverty <- scale(data_pv$total_poverty, center = TRUE, scale = TRUE)[,1] 

### Random effects
# Temporal
t1 <- as.factor(data_pv$Month) # Seasonality
t2 <- as.factor(data_pv$Year)  # Interannual 

# Spatial
s1 <- rep(1:14, 348) 

df_inla_pv <- data.frame(y, e, prcp, tmin,
                         int_per, urban,  
                         total_poverty, 
                         t1, t2, s1)

########################################################################################

## Baseline model - spatiotemporal effects 

########################################################################################

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") 


mod1_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = FALSE, 
                                       return.marginals = FALSE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))

save(mod1_pv, file = "models/mod1_pv.R")



######################################################################################

## Add t2 random effects

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid")


mod1_2_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                  offset = log(e), verbose = TRUE,
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                         config = FALSE, 
                                         return.marginals = FALSE), 
                  control.predictor = list(link = 1, compute = TRUE), 
                  control.family = list(link = "log"))

save(mod1_2_pv, file = "models/mod1_2_pv.R")

########################################################################################

## Socioeconomic information

########################################################################################

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty 

mod2_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = FALSE, 
                                       return.marginals = FALSE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))

save(mod2_pv, file = "models/mod2_pv.R")

########################################################################################

## Add urban

########################################################################################

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  ## Add interaction
  int_per +
  urban*int_per


mod3_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = FALSE, 
                                       return.marginals = FALSE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))

save(mod3_pv, file = "models/mod3_pv.R")

########################################################################################

## Add temperature

########################################################################################

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  f(inla.group(tmin), model = "rw1") 


mod4_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = FALSE, 
                                       return.marginals = FALSE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))

save(mod4_pv, file = "models/mod4_pv.R")

## Test linear
formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  ## Add interaction with intervention period
  int_per + 
  urban*int_per +
  tmin



mod4_l_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                  offset = log(e), verbose = TRUE,
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                         config = FALSE, 
                                         return.marginals = FALSE), 
                  control.predictor = list(link = 1, compute = TRUE), 
                  control.family = list(link = "log"))

save(mod4_l_pv, file = "models/mod4_l_pv.R")

########################################################################################

## Add precipitation

########################################################################################

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  f(inla.group(tmin), model = "rw1") +
  f(inla.group(prcp), model = "rw1") 

mod5_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = TRUE, 
                                       return.marginals = TRUE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))

save(mod5_pv, file = "models/mod5_pv.R")

## Test linear
formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  ## Add interaction with intervention period
  int_per + 
  urban*int_per +
  tmin +
  prcp


mod5_l_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                  offset = log(e), verbose = TRUE,
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                         config = FALSE, 
                                         return.marginals = FALSE), 
                  control.predictor = list(link = 1, compute = TRUE), 
                  control.family = list(link = "log"))

save(mod5_l_pv, file = "models/mod5_l_pv.R")

### With model without interaction, to plot parameter estimates
formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  tmin +
  prcp


mod6_l_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                  offset = log(e), verbose = TRUE,
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                         config = FALSE, 
                                         return.marginals = FALSE), 
                  control.predictor = list(link = 1, compute = TRUE), 
                  control.family = list(link = "log"))

save(mod6_l_pv, file = "models/mod6_l_pv.R")


########################################################################################

#### Test full model without tmin 

########################################################################################

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  f(inla.group(prcp), model = "rw1") 

mod6_wtmin_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                      offset = log(e), verbose = TRUE,
                      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                             config = FALSE, 
                                             return.marginals = FALSE), 
                      control.predictor = list(link = 1, compute = TRUE), 
                      control.family = list(link = "log"))

save(mod6_wtmin_pv, file = "models/mod6_wtmin_pv.R")

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  prcp

mod6_l_wtmin_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                        offset = log(e), verbose = TRUE,
                        control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                               config = FALSE, 
                                               return.marginals = FALSE), 
                        control.predictor = list(link = 1, compute = TRUE), 
                        control.family = list(link = "log"))

save(mod6_l_wtmin_pv, file = "models/mod6_l_wtmin_pv.R")

########################################################################################

#### Test full model without prcp

########################################################################################

formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  total_poverty +
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  f(inla.group(tmin), model = "rw1") 

mod6_wprcp_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                      offset = log(e), verbose = TRUE,
                      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                             config = FALSE, 
                                             return.marginals = FALSE), 
                      control.predictor = list(link = 1, compute = TRUE), 
                      control.family = list(link = "log"))

save(mod6_wprcp_pv, file = "models/mod6_wprcp_pv.R")
