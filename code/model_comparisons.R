########################################################################################################################################

### Model comparisons

########################################################################################################################################

## Here we create a table of model DIC and WAIC values to select the best fitting model 

library(kableExtra)

### Load models from file
files <- list.files("models", pattern = ".R", full.names = TRUE)

for (i in 1:length(files)) {

load(files[i])

}

## Create df
model_table <- data.frame(Model = c(rep("Baseline", 2),
                                    rep("...+ unstructured yearly random effects", 2),
                                    rep("... + socioeconomic effects", 2),
                                    rep("... + urban", 2),
                                    "... + temperature (linear)", 
                                    "... + temperature (non-linear)", 
                                    "... + temperature (linear)", 
                                    "... + temperature (non-linear)", 
                                    "... + precipitation (linear)", 
                                    "... + precipitation (non-linear)",
                                    "... + precipitation (linear)", 
                                    "... + precipitation (non-linear)"),
                          
                          Parasite = c(rep(c("P. falciparum", "P.vivax"), 4),
                                       rep("P. falciparum",2), 
                                       rep("P. vivax",2), 
                                       rep("P. falciparum",2), 
                                       rep("P. vivax",2)),
                                    
                          
                          DIC   = c(format(mod1_pf$dic$dic,scientific=F), format(mod1_pv$dic$dic,scientific=F),
                                    format(mod1_2_pf$dic$dic,scientific=F), format(mod1_2_pv$dic$dic,scientific=F),
                                    format(mod2_pf$dic$dic,scientific=F), format(mod2_pv$dic$dic,scientific=F),
                                    format(mod3_pf$dic$dic,scientific=F), format(mod3_pv$dic$dic,scientific=F),
                                    
                                    format(mod4_l_pf$dic$dic,scientific=F), format(mod4_pf$dic$dic,scientific=F), 
                                    format(mod4_l_pv$dic$dic,scientific=F), format(mod4_pv$dic$dic,scientific=F),
                                    
                                    format(mod5_l_pf$dic$dic,scientific=F), format(mod5_pf$dic$dic,scientific=F), 
                                    format(mod5_l_pv$dic$dic,scientific=F), format(mod5_pv$dic$dic,scientific=F)),

                          
                          WAIC    = c(format(mod1_pf$waic$waic,scientific=F), format(mod1_pv$waic$waic,scientific=F),
                                      format(mod1_2_pf$waic$waic,scientific=F), format(mod1_2_pv$waic$waic,scientific=F),
                                      format(mod2_pf$waic$waic,scientific=F), format(mod2_pv$waic$waic,scientific=F),
                                      format(mod3_pf$waic$waic,scientific=F), format(mod3_pv$waic$waic,scientific=F),
                                      
                                      format(mod4_l_pf$waic$waic,scientific=F), format(mod4_pf$waic$waic,scientific=F), 
                                      format(mod4_l_pv$waic$waic,scientific=F), format(mod4_pv$waic$waic,scientific=F),
                                      
                                      format(mod5_l_pf$waic$waic,scientific=F), format(mod5_pf$waic$waic,scientific=F), 
                                      format(mod5_l_pv$waic$waic,scientific=F), format(mod5_pv$waic$waic,scientific=F)),
                          
                          ## negative mean natural logarithm of the CPO values  - smaller the better
                          CPO    = c(format(-mean(log(mod1_pf$cpo$cpo), na.rm = T),scientific=F),
                                     format(-mean(log(mod1_pv$cpo$cpo), na.rm = T),scientific=F),
                                     
                                     format(-mean(log(mod1_2_pf$cpo$cpo), na.rm = T),scientific=F),
                                     format(-mean(log(mod1_2_pv$cpo$cpo), na.rm = T),scientific=F),
                                     
                                     format(-mean(log(mod2_pf$cpo$cpo), na.rm = T),scientific=F),
                                     format(-mean(log(mod2_pv$cpo$cpo), na.rm = T),scientific=F),
                                     
                                     format(-mean(log(mod3_pf$cpo$cpo), na.rm = T),scientific=F),
                                     format(-mean(log(mod3_pv$cpo$cpo), na.rm = T),scientific=F),
                                     
                                     
                                     format(-mean(log(mod4_l_pf$cpo$cpo), na.rm = T),scientific=F),
                                     format(-mean(log(mod4_pf$cpo$cpo), na.rm = T),scientific=F),
                                     
                                     format(-mean(log(mod4_l_pv$cpo$cpo), na.rm = T),scientific=F),
                                     format(-mean(log(mod4_pv$cpo$cpo), na.rm = T),scientific=F),
                                     
                                     format(-mean(log(mod5_l_pf$cpo$cpo), na.rm = T),scientific=F),
                                     format(-mean(log(mod5_pf$cpo$cpo), na.rm = T),scientific=F),
                                     
                                     format(-mean(log(mod5_l_pv$cpo$cpo), na.rm = T),scientific=F),
                                     format(-mean(log(mod5_pv$cpo$cpo), na.rm = T),scientific=F)))

model_table$CPO <- round(as.numeric(as.character(model_table$CPO)), 2)

kable(model_table, caption = " ") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE, 
                font_size = 14) %>%
  collapse_rows(1:2) %>%
save_kable("model_comparisons/model_comparisons.pdf")

write.csv(model_table, file = "model_comparisons/model_comparisons.csv")

##############################################################################

### Compare parameter estimates of models for 1990-2018 with intervention models 2001-2015

model_df <- data.frame(Variable = c(rep("Minimum temperature", 4),
                                    
                                    rep("Precipitation", 4),
                                    
                                    rep("Urbanization", 4),
                                    
                                    rep("Poverty", 4)),
                       
                       
                       Parasite = rep(rep(c("P. falciparum", "P. vivax"), each = 2), 4),

                       
                       Model    = rep(c("1990-2018", "2001-2015"), 8),
                       
                       Estimate = c(mod6_pf$summary.fixed$mean[4], mod_int_pf$summary.fixed$mean[2],
                                    mod6_pv$summary.fixed$mean[4], mod_int_pv$summary.fixed$mean[2],
                                    
                                    mod6_pf$summary.fixed$mean[5], mod_int_pf$summary.fixed$mean[3],
                                    mod6_pv$summary.fixed$mean[5], mod_int_pv$summary.fixed$mean[3],
                                    
                                    mod6_pf$summary.fixed$mean[3], mod_int_pf$summary.fixed$mean[4],
                                    mod6_pv$summary.fixed$mean[3], mod_int_pv$summary.fixed$mean[4],
                                    
                                    mod6_pf$summary.fixed$mean[2], mod_int_pf$summary.fixed$mean[5],
                                    mod6_pv$summary.fixed$mean[2], mod_int_pv$summary.fixed$mean[5]),
                       
                       LCI      = c(mod6_pf$summary.fixed$`0.025quant`[4], mod_int_pf$summary.fixed$`0.025quant`[2],
                                    mod6_pv$summary.fixed$`0.025quant`[4], mod_int_pv$summary.fixed$`0.025quant`[2],
                                    
                                    mod6_pf$summary.fixed$`0.025quant`[5], mod_int_pf$summary.fixed$`0.025quant`[3],
                                    mod6_pv$summary.fixed$`0.025quant`[5], mod_int_pv$summary.fixed$`0.025quant`[3],
                                    
                                    mod6_pf$summary.fixed$`0.025quant`[3], mod_int_pf$summary.fixed$`0.025quant`[4],
                                    mod6_pv$summary.fixed$`0.025quant`[3], mod_int_pv$summary.fixed$`0.025quant`[4],
                                    
                                    mod6_pf$summary.fixed$`0.025quant`[2], mod_int_pf$summary.fixed$`0.025quant`[5],
                                    mod6_pv$summary.fixed$`0.025quant`[2], mod_int_pv$summary.fixed$`0.025quant`[5]),
                       
                       UCI      = c(mod6_pf$summary.fixed$`0.975quant`[4], mod_int_pf$summary.fixed$`0.975quant`[2],
                                    mod6_pv$summary.fixed$`0.975quant`[4], mod_int_pv$summary.fixed$`0.975quant`[2],
                                    
                                    mod6_pf$summary.fixed$`0.975quant`[5], mod_int_pf$summary.fixed$`0.975quant`[3],
                                    mod6_pv$summary.fixed$`0.975quant`[5], mod_int_pv$summary.fixed$`0.975quant`[3],
                                    
                                    mod6_pf$summary.fixed$`0.975quant`[3], mod_int_pf$summary.fixed$`0.975quant`[4],
                                    mod6_pv$summary.fixed$`0.975quant`[3], mod_int_pv$summary.fixed$`0.975quant`[4],
                                    
                                    mod6_pf$summary.fixed$`0.975quant`[2], mod_int_pf$summary.fixed$`0.975quant`[5],
                                    mod6_pv$summary.fixed$`0.975quant`[2], mod_int_pv$summary.fixed$`0.975quant`[5]))

model_df <- model_df %>% mutate(Estimate = round(Estimate, 2),
                                LCI      = round(LCI, 2),
                                UCI      = round(UCI, 2))

write.csv(model_df, file = "model_comparisons/model_parameter_estimates.csv")
