########################################################################################################################################

### Table 1

########################################################################################################################################

## Here we create a table of model DIC and WAIC values, and cross-validated log score

library(kableExtra)

### Load models from file
files <- list.files("models", pattern = ".R", full.names = TRUE)

for (i in 1:length(files)) {
  
  load(files[i])
  
}

## Create df
model_table <- data.frame(Model = c(rep("Baseline spatial seasonal", 2),
                                    rep("Unstructured yearly random effects", 2),
                                    rep("Socioeconomic effects", 2),
                                    rep("Urban effects", 2),
                                    "Temperature effects (linear)", 
                                    "Temperature effects (non-linear)", 
                                    "Precipitation effects (linear)", 
                                    "Precipitation effects (non-linear)"),
                          
                          Parasite = c(rep(c("P. falciparum", "P.vivax"), 6)),
                          
                          
                          DIC   = c(format(mod1_pf$dic$dic,scientific=F), format(mod1_pv$dic$dic,scientific=F),
                                    format(mod1_2_pf$dic$dic,scientific=F), format(mod1_2_pv$dic$dic,scientific=F),
                                    format(mod2_pf$dic$dic,scientific=F), format(mod2_pv$dic$dic,scientific=F),
                                    format(mod3_pf$dic$dic,scientific=F), format(mod3_pv$dic$dic,scientific=F),
                                    
                                    format(mod4_l_pf$dic$dic,scientific=F), format(mod4_pv$dic$dic,scientific=F),
                                    
                                    format(mod5_l_pf$dic$dic,scientific=F), format(mod5_nl_pv$dic$dic,scientific=F)),
                          
                          
                          WAIC    = c(format(mod1_pf$waic$waic,scientific=F), format(mod1_pv$waic$waic,scientific=F),
                                      format(mod1_2_pf$waic$waic,scientific=F), format(mod1_2_pv$waic$waic,scientific=F),
                                      format(mod2_pf$waic$waic,scientific=F), format(mod2_pv$waic$waic,scientific=F),
                                      format(mod3_pf$waic$waic,scientific=F), format(mod3_pv$waic$waic,scientific=F),
                                      
                                      format(mod4_l_pf$waic$waic,scientific=F), format(mod4_pv$waic$waic,scientific=F), 
                                      
                                      format(mod5_l_pf$waic$waic,scientific=F), format(mod5_nl_pv$waic$waic,scientific=F)),
                          
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
                                     format(-mean(log(mod4_pv$cpo$cpo), na.rm = T),scientific=F),
                                     
                                     format(-mean(log(mod5_l_pf$cpo$cpo), na.rm = T),scientific=F),
                                     format(-mean(log(mod5_nl_pv$cpo$cpo), na.rm = T),scientific=F)))

model_table$CPO <- round(as.numeric(as.character(model_table$CPO)), 2)

kable(model_table, caption = " ") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE, 
                font_size = 14) %>%
  collapse_rows(1:2) %>%
  save_kable("figures/table_1.pdf")
