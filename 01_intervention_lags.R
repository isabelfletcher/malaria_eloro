####################### Intervention lags #######################

## Here we test for the best monthly time lags for including the three vector control measures into spatiotemporal models of malaria incidence

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
data <- read.csv("data.csv", fileEncoding = 'latin1')

## Subset data to intervention period 2001-2015
data <- subset(data, data$Year < 2016 & data$Year > 2000)

## Replace NA values in intervention data to 0
data$houses_fogged[is.na(data$houses_fogged)] <- 0
data$houses_IRS[is.na(data$houses_IRS)]       <- 0
data$blocks_fumigated[is.na(data$blocks_fumigated)]       <- 0

## Time lags
# nb: lag in multiples of 14 for the 14 cantons
data <- data %>% mutate(
  blocks_fumigated_lag1 = lag(blocks_fumigated, 14),
  blocks_fumigated_lag2 = lag(blocks_fumigated, 28),
  blocks_fumigated_lag3 = lag(blocks_fumigated, 42),
  blocks_fumigated_lag4 = lag(blocks_fumigated, 56),
  blocks_fumigated_lag5 = lag(blocks_fumigated, 70),
  blocks_fumigated_lag6 = lag(blocks_fumigated, 84),
  
  houses_fogged_lag1 = lag(houses_fogged, 14),
  houses_fogged_lag2 = lag(houses_fogged, 28),
  houses_fogged_lag3 = lag(houses_fogged, 42),
  houses_fogged_lag4 = lag(houses_fogged, 56),
  houses_fogged_lag5 = lag(houses_fogged, 70),
  houses_fogged_lag6 = lag(houses_fogged, 84),
  
  houses_IRS_lag1 = lag(houses_IRS, 14),
  houses_IRS_lag2 = lag(houses_IRS, 28),
  houses_IRS_lag3 = lag(houses_IRS, 42),
  houses_IRS_lag4 = lag(houses_IRS, 56),
  houses_IRS_lag5 = lag(houses_IRS, 70),
  houses_IRS_lag6 = lag(houses_IRS, 84))

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

urban <- scale(data_pf$urban, center = TRUE, scale = TRUE)[,1]

## Socioeconomic data
total_poverty <- scale(data_pf$total_poverty, center = TRUE, scale = TRUE)[,1]

### Random effects
# Temporal
t1 <- as.factor(data_pf$Month) # Seasonality
t2 <- as.factor(data_pf$Year)  # Interannual 

# Spatial effects: in order to determine the variation due to spatial autocorrelation (structured effects), we specify the besag model
# We specify the iid distribution to take into account unstructured spatial variation (heterogeneity)
s1 <- rep(1:14, 180) # there are 14 cantons/districts
s2 <- rep(1:14, 180) 


##### Intervention data - Adding different time lags
blocks_fumigated_lag0 <- scale(data_pf$blocks_fumigated, center = TRUE, scale = TRUE)[,1] 
blocks_fumigated_lag1 <- scale(data_pf$blocks_fumigated_lag1, center = TRUE, scale = TRUE)[,1] 
blocks_fumigated_lag2 <- scale(data_pf$blocks_fumigated_lag2, center = TRUE, scale = TRUE)[,1] 
blocks_fumigated_lag3 <- scale(data_pf$blocks_fumigated_lag3, center = TRUE, scale = TRUE)[,1] 
blocks_fumigated_lag4 <- scale(data_pf$blocks_fumigated_lag4, center = TRUE, scale = TRUE)[,1] 
blocks_fumigated_lag5 <- scale(data_pf$blocks_fumigated_lag5, center = TRUE, scale = TRUE)[,1] 
blocks_fumigated_lag6 <- scale(data_pf$blocks_fumigated_lag6, center = TRUE, scale = TRUE)[,1] 

houses_IRS_lag0 <- scale(data_pf$houses_IRS, center = TRUE, scale = TRUE)[,1] 
houses_IRS_lag1 <- scale(data_pf$houses_IRS_lag1, center = TRUE, scale = TRUE)[,1] 
houses_IRS_lag2 <- scale(data_pf$houses_IRS_lag2, center = TRUE, scale = TRUE)[,1] 
houses_IRS_lag3 <- scale(data_pf$houses_IRS_lag3, center = TRUE, scale = TRUE)[,1] 
houses_IRS_lag4 <- scale(data_pf$houses_IRS_lag4, center = TRUE, scale = TRUE)[,1] 
houses_IRS_lag5 <- scale(data_pf$houses_IRS_lag5, center = TRUE, scale = TRUE)[,1] 
houses_IRS_lag6 <- scale(data_pf$houses_IRS_lag6, center = TRUE, scale = TRUE)[,1] 

houses_fogged_lag0 <- scale(data_pf$houses_fogged, center = TRUE, scale = TRUE)[,1] 
houses_fogged_lag1 <- scale(data_pf$houses_fogged_lag1, center = TRUE, scale = TRUE)[,1] 
houses_fogged_lag2 <- scale(data_pf$houses_fogged_lag2, center = TRUE, scale = TRUE)[,1] 
houses_fogged_lag3 <- scale(data_pf$houses_fogged_lag3, center = TRUE, scale = TRUE)[,1] 
houses_fogged_lag4 <- scale(data_pf$houses_fogged_lag4, center = TRUE, scale = TRUE)[,1] 
houses_fogged_lag5 <- scale(data_pf$houses_fogged_lag5, center = TRUE, scale = TRUE)[,1] 
houses_fogged_lag6 <- scale(data_pf$houses_fogged_lag6, center = TRUE, scale = TRUE)[,1] 

## Assemble into df
df_inla_pf <- data.frame(y, e,
                         houses_fogged_lag0, houses_fogged_lag1, houses_fogged_lag2, houses_fogged_lag3, houses_fogged_lag4, houses_fogged_lag5, houses_fogged_lag6,
                         houses_IRS_lag0, houses_IRS_lag1, houses_IRS_lag2, houses_IRS_lag3, houses_IRS_lag4, houses_IRS_lag5, houses_IRS_lag6,
                         blocks_fumigated_lag0, blocks_fumigated_lag1, blocks_fumigated_lag2, blocks_fumigated_lag3, blocks_fumigated_lag4, blocks_fumigated_lag5, blocks_fumigated_lag6,
                         urban, 
                         total_poverty, 
                         tmin, prcp, 
                         t1, t2, s1, s2)


######### Fogging

fogging_lags <- c("houses_fogged_lag0", "houses_fogged_lag1",
                  "houses_fogged_lag2", "houses_fogged_lag3")

data_pf$houses_fogged_lag0 <- data_pf$houses_fogged
data_pf$houses_fogged <- NULL

df_fogging <- NULL

for (i in c(1:4)) {
  
  # Extract variable combination
  comb_vec <- fogging_lags[i]
  
  ## Create a formula, scaling the covariates
  
  fogging <- scale(data_pf[,grep(comb_vec, names(data_pf), value=TRUE)], 
                   center = TRUE, scale = TRUE)[,1]
  
  formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +      
    f(s2, model = "iid", graph = "map.graph") +
    f(t1, model = "rw1") +
    f(t2, model = "iid") +
    prcp +
    tmin + 
    urban +
    total_poverty +
    houses_IRS_lag0 +
    blocks_fumigated_lag0 +
    fogging 
  
  
  ## Put into model
  model <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = FALSE, 
                                       return.marginals = FALSE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))
  
  ## Table comparison
  fogging_data <- data.frame(Intervention = c("Fogging"),
                             
                             Lag = gsub("houses_fogged_lag", "", comb_vec),
                             
                             Mean = model$summary.fixed$mean[6],
                             
                             LCI = model$summary.fixed$`0.025quant`[6],
                             
                             UCI = model$summary.fixed$`0.975quant`[6],
                             
                             DIC = model$dic$dic,
                             
                             WAIC = model$waic$waic)
  
  df_fogging <- rbind(df_fogging, fogging_data)
  
  rm(model)
  rm(comb_vec)
  
  rm(fogging_data)
  
}

######### Fumigation

fumigation_lags <- c("blocks_fumigated_lag0", "blocks_fumigated_lag1",
                     "blocks_fumigated_lag2", "blocks_fumigated_lag3")

data_pf$blocks_fumigated_lag0 <- data_pf$blocks_fumigated
data_pf$blocks_fumigated <- NULL

df_fumigation <- NULL

for (i in c(1:4)) {
  
  # Extract variable combination
  comb_vec <- fumigation_lags[i]
  
  ## Create a formula, scaling the covariates
  
  fumigation <- scale(data_pf[,grep(comb_vec, names(data_pf), value=TRUE)], 
                      center = TRUE, scale = TRUE)[,1]
  
  formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +      
    f(s2, model = "iid", graph = "map.graph") +
    f(t1, model = "rw1") +
    f(t2, model = "iid") +
    prcp +
    tmin + 
    urban +
    total_poverty +
    houses_IRS_lag0 +
    houses_fogged_lag0 +
    fumigation 
  
  
  ## Put into model
  model <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = FALSE, 
                                       return.marginals = FALSE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))
  
  ## Table comparison
  fumigation_data <- data.frame(Intervention = c("Fumigation"),
                                
                                Lag = gsub("blocks_fumigated_lag", "", comb_vec),
                                
                                Mean = model$summary.fixed$mean[6],
                                
                                LCI = model$summary.fixed$`0.025quant`[6],
                                
                                UCI = model$summary.fixed$`0.975quant`[6],
                                
                                DIC = model$dic$dic,
                                
                                WAIC = model$waic$waic)
  
  df_fumigation <- rbind(df_fumigation, fumigation_data)
  
  rm(model)
  rm(comb_vec)
  
  rm(fumigation_data)
  
}

######### IRS

irs_lags <- c("houses_IRS_lag0", "houses_IRS_lag1",
              "houses_IRS_lag2", "houses_IRS_lag3")

data_pf$houses_IRS_lag0 <- data_pf$houses_IRS
data_pf$houses_IRS <- NULL

df_irs <- NULL

for (i in c(1:4)) {
  
  # Extract variable combination
  comb_vec <- irs_lags[i]
  
  ## Create a formula, scaling the covariates
  
  irs <- scale(data_pf[,grep(comb_vec, names(data_pf), value=TRUE)], 
               center = TRUE, scale = TRUE)[,1]
  
  formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +      
    f(s2, model = "iid", graph = "map.graph") +
    f(t1, model = "rw1") +
    f(t2, model = "iid") +
    prcp +
    tmin + 
    urban +
    total_poverty +
    houses_IRS_lag0 +
    houses_fogged_lag0 +
    irs 
  
  
  ## Put into model
  model <- inla(formula, data = df_inla_pf, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = FALSE, 
                                       return.marginals = FALSE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))
  
  ## Table comparison
  irs_data <- data.frame(       Intervention = c("IRS"),
                                
                                Lag = gsub("houses_IRS_lag", "", comb_vec),
                                
                                Mean = model$summary.fixed$mean[6],
                                
                                LCI = model$summary.fixed$`0.025quant`[6],
                                
                                UCI = model$summary.fixed$`0.975quant`[6],
                                
                                DIC = model$dic$dic,
                                
                                WAIC = model$waic$waic)
  
  df_irs <- rbind(df_irs, irs_data)
  
  rm(model)
  rm(comb_vec)
  
  rm(irs_data)
  
}


data <- rbind(df_irs, df_fogging, df_fumigation)

write.csv(data, "supplementary/table_S1_pf.csv")

#######################################################################################################################################################

#### Vivax models
## Fixed effects
y  <- data_pv$cases
n  <- length(y)

## Add API so modify offset to include it 
e  <- (data_pv$Population/12)/1000

## Climate effects, scaling the covariates
prcp <- scale(data_pv$prcp_lag1, center = TRUE, scale = TRUE)[,1] 
tmin <- scale(data_pv$tmin_lag3, center = TRUE, scale = TRUE)[,1]

urban <- scale(data_pv$urban, center = TRUE, scale = TRUE)[,1]

## Socioeconomic data
total_poverty <- scale(data_pv$total_poverty, center = TRUE, scale = TRUE)[,1]

### Random effects
# Temporal
t1 <- as.factor(data_pv$Month) # Seasonality
t2 <- as.factor(data_pv$Year)  # Interannual 

# Spatial effects: in order to determine the variation due to spatial autocorrelation (structured effects), we specify the besag model
# We specify the iid distribution to take into account unstructured spatial variation (heterogeneity)
s1 <- rep(1:14, 180) # there are 14 cantons/districts
s2 <- rep(1:14, 180) 


##### Intervention data - Adding different time lags
blocks_fumigated_lag0 <- scale(data_pv$blocks_fumigated, center = TRUE, scale = TRUE)[,1] 
blocks_fumigated_lag1 <- scale(data_pv$blocks_fumigated_lag1, center = TRUE, scale = TRUE)[,1] 
blocks_fumigated_lag2 <- scale(data_pv$blocks_fumigated_lag2, center = TRUE, scale = TRUE)[,1] 
blocks_fumigated_lag3 <- scale(data_pv$blocks_fumigated_lag3, center = TRUE, scale = TRUE)[,1] 
blocks_fumigated_lag4 <- scale(data_pv$blocks_fumigated_lag4, center = TRUE, scale = TRUE)[,1] 
blocks_fumigated_lag5 <- scale(data_pv$blocks_fumigated_lag5, center = TRUE, scale = TRUE)[,1] 
blocks_fumigated_lag6 <- scale(data_pv$blocks_fumigated_lag6, center = TRUE, scale = TRUE)[,1] 

houses_IRS_lag0 <- scale(data_pv$houses_IRS, center = TRUE, scale = TRUE)[,1] 
houses_IRS_lag1 <- scale(data_pv$houses_IRS_lag1, center = TRUE, scale = TRUE)[,1] 
houses_IRS_lag2 <- scale(data_pv$houses_IRS_lag2, center = TRUE, scale = TRUE)[,1] 
houses_IRS_lag3 <- scale(data_pv$houses_IRS_lag3, center = TRUE, scale = TRUE)[,1] 
houses_IRS_lag4 <- scale(data_pv$houses_IRS_lag4, center = TRUE, scale = TRUE)[,1] 
houses_IRS_lag5 <- scale(data_pv$houses_IRS_lag5, center = TRUE, scale = TRUE)[,1] 
houses_IRS_lag6 <- scale(data_pv$houses_IRS_lag6, center = TRUE, scale = TRUE)[,1] 

houses_fogged_lag0 <- scale(data_pv$houses_fogged, center = TRUE, scale = TRUE)[,1] 
houses_fogged_lag1 <- scale(data_pv$houses_fogged_lag1, center = TRUE, scale = TRUE)[,1] 
houses_fogged_lag2 <- scale(data_pv$houses_fogged_lag2, center = TRUE, scale = TRUE)[,1] 
houses_fogged_lag3 <- scale(data_pv$houses_fogged_lag3, center = TRUE, scale = TRUE)[,1] 
houses_fogged_lag4 <- scale(data_pv$houses_fogged_lag4, center = TRUE, scale = TRUE)[,1] 
houses_fogged_lag5 <- scale(data_pv$houses_fogged_lag5, center = TRUE, scale = TRUE)[,1] 
houses_fogged_lag6 <- scale(data_pv$houses_fogged_lag6, center = TRUE, scale = TRUE)[,1] 

## Assemble into df
df_inla_pv <- data.frame(y, e,
                         houses_fogged_lag0, houses_fogged_lag1, houses_fogged_lag2, houses_fogged_lag3, houses_fogged_lag4, houses_fogged_lag5, houses_fogged_lag6,
                         houses_IRS_lag0, houses_IRS_lag1, houses_IRS_lag2, houses_IRS_lag3, houses_IRS_lag4, houses_IRS_lag5, houses_IRS_lag6,
                         blocks_fumigated_lag0, blocks_fumigated_lag1, blocks_fumigated_lag2, blocks_fumigated_lag3, blocks_fumigated_lag4, blocks_fumigated_lag5, blocks_fumigated_lag6,
                         urban, 
                         total_poverty, 
                         tmin, prcp, 
                         t1, t2, s1, s2)


######### Fogging

fogging_lags <- c("houses_fogged_lag0", "houses_fogged_lag1",
                  "houses_fogged_lag2", "houses_fogged_lag3")

data_pv$houses_fogged_lag0 <- data_pv$houses_fogged
data_pv$houses_fogged <- NULL

df_fogging <- NULL

for (i in c(1:4)) {
  
  # Extract variable combination
  comb_vec <- fogging_lags[i]
  
  ## Create a formula, scaling the covariates
  
  fogging <- scale(data_pv[,grep(comb_vec, names(data_pv), value=TRUE)], 
                   center = TRUE, scale = TRUE)[,1]
  
  formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +      
    f(s2, model = "iid", graph = "map.graph") +
    f(t1, model = "rw1") +
    f(t2, model = "iid") +
    f(inla.group(prcp), model = "rw1") +
    f(inla.group(tmin), model = "rw1") + 
    urban +
    total_poverty +
    houses_IRS_lag0 +
    blocks_fumigated_lag0 +
    fogging 
  
  
  ## Put into model
  model <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = FALSE, 
                                       return.marginals = FALSE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))
  
  ## Table comparison
  fogging_data <- data.frame(Intervention = c("Fogging"),
                             
                             Lag = gsub("houses_fogged_lag", "", comb_vec),
                             
                             Mean = model$summary.fixed$mean[6],
                             
                             LCI = model$summary.fixed$`0.025quant`[6],
                             
                             UCI = model$summary.fixed$`0.975quant`[6],
                             
                             DIC = model$dic$dic,
                             
                             WAIC = model$waic$waic)
  
  df_fogging <- rbind(df_fogging, fogging_data)
  
  rm(model)
  rm(comb_vec)
  
  rm(fogging_data)
  
}

######### Fumigation

fumigation_lags <- c("blocks_fumigated_lag0", "blocks_fumigated_lag1",
                     "blocks_fumigated_lag2", "blocks_fumigated_lag3")

data_pv$blocks_fumigated_lag0 <- data_pv$blocks_fumigated
data_pv$blocks_fumigated <- NULL

df_fumigation <- NULL

for (i in c(1:4)) {
  
  # Extract variable combination
  comb_vec <- fumigation_lags[i]
  
  ## Create a formula, scaling the covariates
  
  fumigation <- scale(data_pv[,grep(comb_vec, names(data_pv), value=TRUE)], 
                      center = TRUE, scale = TRUE)[,1]
  
  formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +      
    f(s2, model = "iid", graph = "map.graph") +
    f(t1, model = "rw1") +
    f(t2, model = "iid") +
    f(inla.group(prcp), model = "rw1") +
    f(inla.group(tmin), model = "rw1") + 
    urban +
    total_poverty +
    houses_IRS_lag0 +
    houses_fogged_lag0 +
    fumigation 
  
  
  ## Put into model
  model <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = FALSE, 
                                       return.marginals = FALSE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))
  
  ## Table comparison
  fumigation_data <- data.frame(Intervention = c("Fumigation"),
                                
                                Lag = gsub("blocks_fumigated_lag", "", comb_vec),
                                
                                Mean = model$summary.fixed$mean[6],
                                
                                LCI = model$summary.fixed$`0.025quant`[6],
                                
                                UCI = model$summary.fixed$`0.975quant`[6],
                                
                                DIC = model$dic$dic,
                                
                                WAIC = model$waic$waic)
  
  df_fumigation <- rbind(df_fumigation, fumigation_data)
  
  rm(model)
  rm(comb_vec)
  
  rm(fumigation_data)
  
}

######### IRS

irs_lags <- c("houses_IRS_lag0", "houses_IRS_lag1",
              "houses_IRS_lag2", "houses_IRS_lag3")

data_pv$houses_IRS_lag0 <- data_pv$houses_IRS
data_pv$houses_IRS <- NULL

df_irs <- NULL

for (i in c(1:4)) {
  
  # Extract variable combination
  comb_vec <- irs_lags[i]
  
  ## Create a formula, scaling the covariates
  
  irs <- scale(data_pv[,grep(comb_vec, names(data_pv), value=TRUE)], 
               center = TRUE, scale = TRUE)[,1]
  
  formula <- y ~ 1 + f(s1, model = "besag", graph = "map.graph") +      
    f(s2, model = "iid", graph = "map.graph") +
    f(t1, model = "rw1") +
    f(t2, model = "iid") +
    f(inla.group(prcp), model = "rw1") +
    f(inla.group(tmin), model = "rw1") + 
    urban +
    total_poverty +
    houses_IRS_lag0 +
    houses_fogged_lag0 +
    irs 
  
  
  ## Put into model
  model <- inla(formula, data = df_inla_pv, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = FALSE, 
                                       return.marginals = FALSE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))
  
  ## Table comparison
  irs_data <- data.frame(       Intervention = c("IRS"),
                                
                                Lag = gsub("houses_IRS_lag", "", comb_vec),
                                
                                Mean = model$summary.fixed$mean[6],
                                
                                LCI = model$summary.fixed$`0.025quant`[6],
                                
                                UCI = model$summary.fixed$`0.975quant`[6],
                                
                                DIC = model$dic$dic,
                                
                                WAIC = model$waic$waic)
  
  df_irs <- rbind(df_irs, irs_data)
  
  rm(model)
  rm(comb_vec)
  
  rm(irs_data)
  
}


data <- rbind(df_irs, df_fogging, df_fumigation)

write.csv(data, "supplementary/table_S1_pv.csv")

