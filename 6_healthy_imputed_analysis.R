#imputed analysis using pre 2013 data_cutoff

# Clearing workspace
rm(list = ls())
dev.off()

#set working directory
setwd("C:/Users/bruno/OneDrive/Documents/Code/projects/VAREANZ/Honours")

#source helper functions (also loads in libraries)
source("0_helper_functions.R")

#output directory
output_dir = "outputs/healthy_imputed_analysis/"

#1. import data and cleaning up data ####
an_cohort_ts = read.fst("data/cohorts/an_cohort_ts_2008-01-01_2013-01-01.fst",
                     as.data.table = TRUE)
cohort = copy(an_cohort_ts)

# remove DHB
cohort[,DHB_name := NULL]

#cut factors
cohort = cut_factors(cohort)

# 2. Initial Subsetting and demogrpahics ####
#subetting data - need to center age before adding interaction terms..
raw_subsets = list("Female" = 
                     cohort[gender_code == 0,][,gender_code := NULL],
                   "Male" = 
                     cohort[gender_code == 1,][,gender_code := NULL])
#center age and store results
base_age = c()
for (raw_subset in raw_subsets){
  base_age = c(base_age,raw_subset[,mean(nhi_age)])
  raw_subset[,c_nhi_age := (nhi_age - mean(nhi_age))]
  raw_subset[,nhi_age := NULL]
}
names(base_age) = names(raw_subsets)
base_age


# 3. Creating cohorts *** ######################

# 3a. Comparison cohorts ####
# add prespecified interactions and removes test safe variables
comp_subsets = list()
for (i in 1:length(raw_subsets)){
  
  #copy data table for safety
  dt = copy(raw_subsets[[i]])
  
  # remove TS cols
  TS_vars = c("egfr",
              "hba1c",
              "tchdl",
              "tri")
  cols_to_remove = sapply(names(dt), function(col) {any(sapply(TS_vars,function(x) grepl(x,col)))})
  dt = dt[, !cols_to_remove, with = FALSE]

  #add suneelas interactions
  dt = prespec_ints(dt)
  
  #list,name, append
  list_dt = list(dt)
  names(list_dt) = paste(names(raw_subsets)[i],"comparison",sep = " ")
  comp_subsets = append(comp_subsets,list_dt)
}

# 4. imputation #####################

# dont include VSIMPLE index or outcomes in impuataion as one will crash it and other will bias
imputed_subsets = list()
for (i in 1:length(raw_subsets)){
  
  # copy for safety
  dt = copy(raw_subsets[[i]])
  
  # impute at health values - could consider different imputations for men and women
  dt[, `:=`(
    egfr = fifelse(is.na(egfr), 90, egfr),
    hba1c = fifelse(is.na(hba1c), 40, hba1c),
    tchdl = fifelse(is.na(tchdl), 3.5, tchdl),
    tri = fifelse(is.na(tri), 1.5, tri)
  )]
  
  #list,name, append
  list_dt = list(dt)
  names(list_dt) = names(raw_subsets)[i]
  imputed_subsets = append(imputed_subsets,list_dt)  
}

# 4a. New cohort
subsets = list()
for (i in 1:length(raw_subsets)){
  
  #copy for safety
  dt = copy(imputed_subsets[[i]])
  
  #add old interactions
  dt = prespec_ints(dt)
  
  # transform hba1c and tchdl
  dt = transform_labs_verbose(dt,
                   c_hba1c = 40,
                   s_hba1c = 10,
                   c_tchdl = 3.5,
                   c_tri = 1.5)
  
  #cut up egfr
  dt[,":="(egfr_cat = egfr_groups(c_egfr),
           c_egfr = NULL)]
  dt = cut_factors(dt)
  
  #new interaction
  dt = dt[,":="(age_tchdl = c_nhi_age*c_tchdl)]
  
  #list,name, append
  list_dt = list(dt)
  names(list_dt) = names(raw_subsets)[i]
  subsets = append(subsets,list_dt)  
}

# 5. fitting models ####
# 5a. Comparison models
comp_models = list()
for (i in 1:2){
  #copy for safety
  dt = copy(comp_subsets[[i]])
  
  #remove index
  dt = dt[,-c("VSIMPLE_INDEX_MASTER")]
  
  #fit model
  model = coxph(Surv(total_fu,indicator) ~ .,
                data = dt,
                ties = "breslow",
                model = TRUE)
  
  #list name append
  list_model = list(model)
  names(list_model) = names(comp_subsets)[i]
  comp_models = append(comp_models,list_model)
}

# 5b. models 
models = list()
for (i in 1:2){
  #copy for safety
  dt = copy(subsets[[i]])
  
  #remove index
  dt = dt[,-c("VSIMPLE_INDEX_MASTER")]
  
  #fit model
  model = coxph(Surv(total_fu,indicator) ~ .,
                data = dt,
                ties = "breslow",
                model = TRUE)
  
  #list name append
  list_model = list(model)
  names(list_model) = names(raw_subsets)[i]
  models = append(models,list_model)
}


# MODEL VALIDATION ################################################
###########################



# 6. create prediction dts ####
#reformat models into embedded list
female_models = list("New Model" = copy(models[[1]]),
                     "Baseline Model" = copy(comp_models[[1]]))
male_models = list("New Model" = copy(models[[2]]),
                   "Baseline Model" = copy(comp_models[[2]]))
all_models = list("women" = female_models,
                  "men" = male_models)

pred_dt_list = list()
for (i in 1:2){
  model_list = all_models[[i]]
  new_data = copy(subsets[[i]])
  pred_dt = pred_dt_generator(model_list,
                          new_data)
  
  #list, name, append
  list_pred_dt = list(pred_dt)
  names(list_pred_dt) = names(all_models)[i]
  pred_dt_list = append(pred_dt_list,list_pred_dt)
}

# 7. Demographics tables ####
#simply run the demographics table function on the male and female subsets
source("0_helper_functions.R")
subset_list = list("Women" = list(subsets[[1]][,VSIMPLE_INDEX_MASTER]),
                   "Men" = list(subsets[[2]][,VSIMPLE_INDEX_MASTER]))
demo_table = demo_table_generator(subset_list,
                                  test_rows = TRUE)
demo_table
write.csv(demo_table,
          file = paste0(output_dir,"tables/demo_table_cc_gender.csv"))

# 8. model coefficients

# 9a. coefficients plot ####
source("0_helper_functions.R")
model_list_list = copy(all_models)
coef_plot = coef_plot_generator(model_list_list,
                             return = 1)
coef_plot
ggsave("coef_plot_cc.jpg",
       coef_plot,
       path = paste0(output_dir,"plots/misc"),
       height = 7,
       width = 9)

################################################################################
# 11. Custom SUBSETS ####
#################################################################################


# 11a. flexible subseting
gender_ind = "women" #men or women
cond = "nhi_age >= 30 & nhi_age < 45 & eth_cat == 'NZM'"
cond = "nhi_age >= 60 "
cond = "nhi_age >= 30" #this will select all people
subpop_index_temp = subpop_index(cond)
#subsets[[1]][VSIMPLE_INDEX_MASTER %in% subpop_index_temp]
subpop_pred_dt = pred_dt_list[[gender_ind]][
  VSIMPLE_INDEX_MASTER %in% subpop_index_temp,]
#plot using name of condition
cal_plot = validation(subpop_pred_dt,
                      type = "cal",
                      title = cond)
cal_plot

pred = subpop_pred_dt[,"New Model"][[1]]
obs = subpop_pred_dt[,"indicator"][[1]]

#####################################################################################
# Boostrapping ####
####################################################################################
# Bootstrapping
boot_output = pred_dt_bootstrap_wrapper(subpop_pred_dt)
#model_boot_wrapper(subsets[[1]],comp_subsets[[1]])

################################################################################
# Re_distribution of predictions plot ####
##############################################################################

# selecting which individuals to plot
# random sample of n1 non events and n2 events or just a random sample of n individuals
plot_dt = copy(pred_dt_list[["women"]])
sampled_1 <- plot_dt[indicator==1][sample(.N, 100), ]
sampled_0 <- plot_dt[indicator==0][sample(.N, 250), ]  
# Combine the sampled data back into one data table
plot_dt <- rbind(sampled_1, sampled_0)
#plot_dt = plot_dt[sample(.N,500),]
plot_dt = plot_dt[order(plot_dt$indicator),]

# Create the scatter plot with adjusted data order
plot = ggplot(plot_dt, aes(x = `Baseline Model`, y = `New Model`, color = factor(indicator))) +
  geom_point(aes(shape = factor(indicator)), size = 3, alpha = 0.6) +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  labs(title = "Comparison of Baseline and New Model Predictions",
       x = "Baseline Model Prediction",
       y = "New Model Prediction",
       color = "Indicator",
       shape = "Indicator") +
  theme_minimal() +
  theme(legend.position = "top") +
  coord_fixed(ratio = 1)  # Ensures the plot has an aspect ratio of 1:1


##################################################################################
# 9. Calibration ####  
##################################################################################

# 9. plotting 
#while working so save scrolling up
source("0_helper_functions.R")

# 9a. whole cohort calibration
#plot calibration
whole_cal_plot_list = whole_cohort_plot_generator(pred_dt_list,
                                                  type = "calibration",
                                                  display_n = TRUE,
                                                  display_metrics = TRUE)

# 9b. sub population plots
full_cal_plot_list = cal_sub_pop_plot_generator(pred_dt_list,
                                                all_conds,
                                                display_n = TRUE,
                                                display_metrics = TRUE  )

# 10. Discrimination 
# 8a. whole cohort discrimination
whole_dis_plot_list = whole_cohort_plot_generator(pred_dt_list,
                                                    type = "discrimination",
                                                    display_n = TRUE,
                                                    display_metrics = TRUE)


# # 11b. arranging premade plots
# names(full_cal_plot_list)
# selected_plots = c("Maaori women",
#                    "Maaori men",
#                    "Women 60-74",
#                    "Men 60-74",
#                    "Women Deprivation Q5",
#                    "Men Deprivation Q5")
# selected_plot_list = full_cal_plot_list[selected_plots]
# names(selected_plot_list)
# custom_plot = ggarrange(plotlist = selected_plot_list,
#                         ncol = 2,
#                         nrow = 3)
# custom_plot
# ggsave("custom_calibration_plot.png",
#        custom_plot,
#        path = "plots/validation_plots/calibration",
#        height = 15,
#        width = 9)

# names(full_cal_plot_list)
# selected_plots = c("Women 30-44",
#                    "Men 30-44",
#                    "Women 45-59",
#                    "Men 45-59",
#                    "Women 60-74",
#                    "Men 60-74")
# selected_plot_list = full_cal_plot_list[selected_plots]
# names(selected_plot_list)
# custom_plot = ggarrange(plotlist = selected_plot_list,
#                         ncol = 2,
#                         nrow = 3)
# custom_plot
# ggsave("custom_2_calibration_plot.png",
#        custom_plot,
#        path = "plots/validation_plots/calibration",
#        height = 15,
#        width = 9)


# whole_pop_plots = ggarrange(plotlist = list(whole_cal_list[[1]],
#                                             whole_cal_list[[2]],
#                                             whole_dis_list[[1]],
#                                             whole_dis_list[[2]]),
#             ncol = 2,
#             nrow = 2)
# whole_pop_plots
# ggsave("custom_3_calibration_plot.png",
#        whole_pop_plots,
#        path = "plots/validation_plots/calibration",
#        height = 10,
#        width = 9)

# 13 Net reclassification  ####
reclassifictaion_dt = list(reclassification_testing(pred_dt_list[["women"]]),
                          reclassification_testing(pred_dt_list[["men"]]))
