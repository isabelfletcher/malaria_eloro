####################### Temperature lags #######################

## Here we test for the best monthly time lags for maximum temperature

## Load libraries
pacman::p_load("raster", "INLA","dplyr", 
               "kableExtra", "reshape2",
               "spdep")

### Define neighbourhood matrix
ecuador <- getData('GADM', country = "ECU", level = 2)
el_oro  <- subset(ecuador, NAME_1 == "El Oro")

## Create neighbourhood matrix
nb.map <- poly2nb(el_oro)
nb2INLA("map.graph",nb.map)

## Read in data
data <- read.csv("data/inla_input/data.csv")
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
tmax <- scale(data_pf$tmax, center = TRUE, scale = TRUE)[,1]

urban <- scale(data_pf$urban, center = TRUE, scale = TRUE)[,1]

## Socioeconomic data
total_poverty <- scale(data_pf$total_poverty, center = TRUE, scale = TRUE)[,1]

### Random effects
# Temporal
t1 <- as.factor(data_pf$Month) # Seasonality
t2 <- as.factor(data_pf$Year)  # Interannual 

# Spatial effects: in order to determine the variation due to spatial autocorrelation (structured effects), we specify the besag model
# We specify the iid distribution to take into account unstructured spatial variation (heterogeneity)
s1 <- rep(1:14, 348) # there are 14 cantons/districts
s2 <- rep(1:14, 348) 


##################################
# tmax - add time lags
tmax <- scale(data_pf$tmax, center = TRUE, scale = TRUE)[,1]
tmax_lag1 <- scale(data_pf$tmax_lag1, center = TRUE, scale = TRUE)[,1]
tmax_lag2 <- scale(data_pf$tmax_lag2, center = TRUE, scale = TRUE)[,1]
tmax_lag3 <- scale(data_pf$tmax_lag3, center = TRUE, scale = TRUE)[,1]

df_inla_pf <- data.frame(y, prcp, tmax, tmax_lag1,
                         tmax_lag2, tmax_lag3, 
                         t1, t2, s1, s2, total_poverty, urban, int_per)


###########################################################################################################
## Test with all other covariates
formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
  f(s2, model = "iid", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  tmax +
  prcp + 
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  total_poverty 


mod_tmax_lag0_l_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                           offset = log(e), verbose = TRUE,
                           control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                  config = FALSE, 
                                                  return.marginals = FALSE), 
                           control.predictor = list(link = 1, compute = TRUE), 
                           control.family = list(link = "log"))

formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
  f(s2, model = "iid", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(tmax), model = "rw1") +
  f(inla.group(prcp), model = "rw1") + 
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  total_poverty 

mod_tmax_lag0_nl_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                            offset = log(e), verbose = TRUE,
                            control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                   config = FALSE, 
                                                   return.marginals = FALSE), 
                            control.predictor = list(link = 1, compute = TRUE), 
                            control.family = list(link = "log"))

formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
  f(s2, model = "iid", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  tmax_lag1 +
  prcp + 
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  total_poverty 


mod_tmax_lag1_l_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                           offset = log(e), verbose = TRUE,
                           control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                  config = FALSE, 
                                                  return.marginals = FALSE), 
                           control.predictor = list(link = 1, compute = TRUE), 
                           control.family = list(link = "log"))

formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
  f(s2, model = "iid", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(tmax_lag1), model = "rw1") +
  f(inla.group(prcp), model = "rw1") + 
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  total_poverty 


mod_tmax_lag1_nl_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                            offset = log(e), verbose = TRUE,
                            control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                   config = FALSE, 
                                                   return.marginals = FALSE), 
                            control.predictor = list(link = 1, compute = TRUE), 
                            control.family = list(link = "log"))

formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
  f(s2, model = "iid", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  tmax_lag2 +
  prcp + 
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  total_poverty 


mod_tmax_lag2_l_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                           offset = log(e), verbose = TRUE,
                           control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                  config = FALSE, 
                                                  return.marginals = FALSE), 
                           control.predictor = list(link = 1, compute = TRUE), 
                           control.family = list(link = "log"))

formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
  f(s2, model = "iid", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(tmax_lag2), model = "rw1") +
  f(inla.group(prcp), model = "rw1") + 
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  total_poverty 

mod_tmax_lag2_nl_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                            offset = log(e), verbose = TRUE,
                            control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                   config = FALSE, 
                                                   return.marginals = FALSE), 
                            control.predictor = list(link = 1, compute = TRUE), 
                            control.family = list(link = "log"))


formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
  f(s2, model = "iid", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  tmax_lag3 +
  prcp + 
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  total_poverty 

mod_tmax_lag3_l_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                           offset = log(e), verbose = TRUE,
                           control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                  config = FALSE, 
                                                  return.marginals = FALSE), 
                           control.predictor = list(link = 1, compute = TRUE), 
                           control.family = list(link = "log"))

formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
  f(s2, model = "iid", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(tmax_lag3), model = "rw1") +
  f(inla.group(prcp), model = "rw1") + 
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  total_poverty 


mod_tmax_lag3_nl_pf <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                            offset = log(e), verbose = TRUE,
                            control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                   config = FALSE, 
                                                   return.marginals = FALSE), 
                            control.predictor = list(link = 1, compute = TRUE), 
                            control.family = list(link = "log"))


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
tmax <- scale(data_pv$tmax, center = TRUE, scale = TRUE)[,1]

urban <- scale(data_pv$urban, center = TRUE, scale = TRUE)[,1]

## Socioeconomic data_pv
total_poverty <- scale(data_pv$total_poverty, center = TRUE, scale = TRUE)[,1]

### Random effects
# Temporal
t1 <- as.factor(data_pv$Month) # Seasonality
t2 <- as.factor(data_pv$Year)  # Interannual 

# Spatial effects: in order to determine the variation due to spatial autocorrelation (structured effects), we specify the besag model
# We specify the iid distribution to take into account unstructured spatial variation (heterogeneity)
s1 <- rep(1:14, 348) # there are 14 cantons/districts
s2 <- rep(1:14, 348) 


##################################
# tmax - add time lags
tmax <- scale(data_pv$tmax, center = TRUE, scale = TRUE)[,1]
tmax_lag1 <- scale(data_pv$tmax_lag1, center = TRUE, scale = TRUE)[,1]
tmax_lag2 <- scale(data_pv$tmax_lag2, center = TRUE, scale = TRUE)[,1]
tmax_lag3 <- scale(data_pv$tmax_lag3, center = TRUE, scale = TRUE)[,1]

df_inla_pv <- data.frame(y, prcp, tmax, tmax_lag1,
                         tmax_lag2, tmax_lag3, 
                         t1, t2, s1, s2, urban, total_poverty, int_per)

###########################################################################################################
## Test with all other covariates
formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
  f(s2, model = "iid", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  tmax +
  prcp + 
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  total_poverty 


mod_tmax_lag0_l_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                           offset = log(e), verbose = TRUE,
                           control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                  config = FALSE, 
                                                  return.marginals = FALSE), 
                           control.predictor = list(link = 1, compute = TRUE), 
                           control.family = list(link = "log"))

formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
  f(s2, model = "iid", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(tmax), model = "rw1") +
  f(inla.group(prcp), model = "rw1") + 
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  total_poverty 

mod_tmax_lag0_nl_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                            offset = log(e), verbose = TRUE,
                            control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                   config = FALSE, 
                                                   return.marginals = FALSE), 
                            control.predictor = list(link = 1, compute = TRUE), 
                            control.family = list(link = "log"))

formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
  f(s2, model = "iid", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  tmax_lag1 +
  prcp + 
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  total_poverty 


mod_tmax_lag1_l_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                           offset = log(e), verbose = TRUE,
                           control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                  config = FALSE, 
                                                  return.marginals = FALSE), 
                           control.predictor = list(link = 1, compute = TRUE), 
                           control.family = list(link = "log"))

formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
  f(s2, model = "iid", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(tmax_lag1), model = "rw1") +
  f(inla.group(prcp), model = "rw1") + 
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  total_poverty 


mod_tmax_lag1_nl_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                            offset = log(e), verbose = TRUE,
                            control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                   config = FALSE, 
                                                   return.marginals = FALSE), 
                            control.predictor = list(link = 1, compute = TRUE), 
                            control.family = list(link = "log"))

formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
  f(s2, model = "iid", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  tmax_lag2 +
  prcp + 
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  total_poverty 


mod_tmax_lag2_l_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                           offset = log(e), verbose = TRUE,
                           control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                  config = FALSE, 
                                                  return.marginals = FALSE), 
                           control.predictor = list(link = 1, compute = TRUE), 
                           control.family = list(link = "log"))

formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
  f(s2, model = "iid", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(tmax_lag2), model = "rw1") +
  f(inla.group(prcp), model = "rw1") + 
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  total_poverty 


mod_tmax_lag2_nl_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                            offset = log(e), verbose = TRUE,
                            control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                   config = FALSE, 
                                                   return.marginals = FALSE), 
                            control.predictor = list(link = 1, compute = TRUE), 
                            control.family = list(link = "log"))


formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
  f(s2, model = "iid", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  tmax_lag3 +
  prcp + 
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  total_poverty 


mod_tmax_lag3_l_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                           offset = log(e), verbose = TRUE,
                           control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                  config = FALSE, 
                                                  return.marginals = FALSE), 
                           control.predictor = list(link = 1, compute = TRUE), 
                           control.family = list(link = "log"))

formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +   
  f(s2, model = "iid", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(tmax_lag3), model = "rw1") +
  f(inla.group(prcp), model = "rw1") + 
  urban +
  ## Add interaction
  int_per + 
  urban*int_per +
  total_poverty 


mod_tmax_lag3_nl_pv <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                            offset = log(e), verbose = TRUE,
                            control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                                   config = FALSE, 
                                                   return.marginals = FALSE), 
                            control.predictor = list(link = 1, compute = TRUE), 
                            control.family = list(link = "log"))


########################################################################################################


########################################################################################################

### Model summaries

tmax_table <- data.frame(Parasite = c(rep("P. falciparum", 4),
                                      rep("P. vivax", 4),
                                      rep("P. falciparum", 4),
                                      rep("P. vivax", 4)),
                         
                         Lag = c(0:3, 0:3, 0:3, 0:3),
                         
                         Relationship = c(rep("Linear", 8),
                                          rep("Non-linear", 8)), 
                         
                         DIC      = c(mod_tmax_lag0_l_pf$dic$dic, mod_tmax_lag1_l_pf$dic$dic,
                                      mod_tmax_lag2_l_pf$dic$dic, mod_tmax_lag3_l_pf$dic$dic,
                                      
                                      mod_tmax_lag0_l_pv$dic$dic,mod_tmax_lag1_l_pv$dic$dic,
                                      mod_tmax_lag2_l_pv$dic$dic,mod_tmax_lag3_l_pv$dic$dic,
                                      
                                      mod_tmax_lag0_nl_pf$dic$dic,mod_tmax_lag1_nl_pf$dic$dic,
                                      mod_tmax_lag2_nl_pf$dic$dic,mod_tmax_lag3_nl_pf$dic$dic,
                                      
                                      mod_tmax_lag0_nl_pv$dic$dic,mod_tmax_lag1_nl_pv$dic$dic,
                                      mod_tmax_lag2_nl_pv$dic$dic,mod_tmax_lag3_nl_pv$dic$dic),
                         
                         WAIC      = c(mod_tmax_lag0_l_pf$waic$waic,mod_tmax_lag1_l_pf$waic$waic,
                                       mod_tmax_lag2_l_pf$waic$waic,mod_tmax_lag3_l_pf$waic$waic,
                                       
                                       mod_tmax_lag0_l_pv$waic$waic,mod_tmax_lag1_l_pv$waic$waic,
                                       mod_tmax_lag2_l_pv$waic$waic,mod_tmax_lag3_l_pv$waic$waic,
                                       
                                       mod_tmax_lag0_nl_pf$waic$waic,mod_tmax_lag1_nl_pf$waic$waic,
                                       mod_tmax_lag2_nl_pf$waic$waic,mod_tmax_lag3_nl_pf$waic$waic,
                                       
                                       mod_tmax_lag0_nl_pv$waic$waic,mod_tmax_lag1_nl_pv$waic$waic,
                                       mod_tmax_lag2_nl_pv$waic$waic,mod_tmax_lag3_nl_pv$waic$waic),
                         
                         Estimate = c(mod_tmax_lag0_l_pf$summary.fixed$mean[2], mod_tmax_lag1_l_pf$summary.fixed$mean[2],
                                      mod_tmax_lag2_l_pf$summary.fixed$mean[2], mod_tmax_lag3_l_pf$summary.fixed$mean[2],
                                      
                                      mod_tmax_lag0_l_pv$summary.fixed$mean[2], mod_tmax_lag1_l_pv$summary.fixed$mean[2],
                                      mod_tmax_lag2_l_pv$summary.fixed$mean[2], mod_tmax_lag3_l_pv$summary.fixed$mean[2],
                                      
                                      mod_tmax_lag0_nl_pf$summary.fixed$mean[2], mod_tmax_lag1_nl_pf$summary.fixed$mean[2],
                                      mod_tmax_lag2_nl_pf$summary.fixed$mean[2], mod_tmax_lag3_nl_pf$summary.fixed$mean[2],
                                      
                                      mod_tmax_lag0_nl_pv$summary.fixed$mean[2], mod_tmax_lag1_nl_pv$summary.fixed$mean[2],
                                      mod_tmax_lag2_nl_pv$summary.fixed$mean[2], mod_tmax_lag3_nl_pv$summary.fixed$mean[2]),
                         
                         LCI      = c(mod_tmax_lag0_l_pf$summary.fixed$`0.025quant`[2], mod_tmax_lag1_l_pf$summary.fixed$`0.025quant`[2],
                                      mod_tmax_lag2_l_pf$summary.fixed$`0.025quant`[2], mod_tmax_lag3_l_pf$summary.fixed$`0.025quant`[2],
                                      
                                      mod_tmax_lag0_l_pv$summary.fixed$`0.025quant`[2], mod_tmax_lag1_l_pv$summary.fixed$`0.025quant`[2],
                                      mod_tmax_lag2_l_pv$summary.fixed$`0.025quant`[2], mod_tmax_lag3_l_pv$summary.fixed$`0.025quant`[2],
                                      
                                      mod_tmax_lag0_nl_pf$summary.fixed$`0.025quant`[2], mod_tmax_lag1_nl_pf$summary.fixed$`0.025quant`[2],
                                      mod_tmax_lag2_nl_pf$summary.fixed$`0.025quant`[2], mod_tmax_lag3_nl_pf$summary.fixed$`0.025quant`[2],
                                      
                                      mod_tmax_lag0_nl_pv$summary.fixed$`0.025quant`[2], mod_tmax_lag1_nl_pv$summary.fixed$`0.025quant`[2],
                                      mod_tmax_lag2_nl_pv$summary.fixed$`0.025quant`[2], mod_tmax_lag3_nl_pv$summary.fixed$`0.025quant`[2]),
                         
                         UCI      = c(mod_tmax_lag0_l_pf$summary.fixed$`0.975quant`[2], mod_tmax_lag1_l_pf$summary.fixed$`0.975quant`[2],
                                      mod_tmax_lag2_l_pf$summary.fixed$`0.975quant`[2], mod_tmax_lag3_l_pf$summary.fixed$`0.975quant`[2],
                                      
                                      mod_tmax_lag0_l_pv$summary.fixed$`0.975quant`[2], mod_tmax_lag1_l_pv$summary.fixed$`0.975quant`[2],
                                      mod_tmax_lag2_l_pv$summary.fixed$`0.975quant`[2], mod_tmax_lag3_l_pv$summary.fixed$`0.975quant`[2],
                                      
                                      mod_tmax_lag0_nl_pf$summary.fixed$`0.975quant`[2], mod_tmax_lag1_nl_pf$summary.fixed$`0.975quant`[2],
                                      mod_tmax_lag2_nl_pf$summary.fixed$`0.975quant`[2], mod_tmax_lag3_nl_pf$summary.fixed$`0.975quant`[2],
                                      
                                      mod_tmax_lag0_nl_pv$summary.fixed$`0.975quant`[2], mod_tmax_lag1_nl_pv$summary.fixed$`0.975quant`[2],
                                      mod_tmax_lag2_nl_pv$summary.fixed$`0.975quant`[2], mod_tmax_lag3_nl_pv$summary.fixed$`0.975quant`[2]))

tmax_table <- tmax_table %>% mutate(Estimate = round(Estimate, 2),
                                    LCI      = round(LCI, 2),
                                    UCI      = round(UCI, 2))


kable(tmax_table, caption = " ") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE, 
                font_size = 14) %>%
  collapse_rows(1:2, valign = "top") %>%
  save_kable("model_comparisons/tmax_lags.pdf")

write.csv(tmax_table, file = "model_comparisons/tmax_lags.csv")
