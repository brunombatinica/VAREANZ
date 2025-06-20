#imputed analysis using pre 2013 data_cutoff

# Clearing workspace
rm(list = ls())
dev.off()

#set working directory
setwd("C:/Users/bruno/OneDrive/Documents/Code/projects/VAREANZ/Honours")

#source helper functions (also loads in libraries)
source("0_helper_functions.R")

# output directory
# Create output directory if it doesn't exist
output_dir = "outputs/cv_rf_imputed_analysis/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir,recursive = TRUE) 
}

###########################################################################################
#1. import data and cleaning up data ####
###########################################################################################

an_cohort_ts = read.fst("data/cohorts/an_cohort_ts_2008-01-01_2013-01-01.fst",
                     as.data.table = TRUE)
cohort = copy(an_cohort_ts)
cohort[,DHB_name := NULL] # remove DHB
cohort = cut_factors(cohort) #cut factors

# 2. Initial Subsetting and demogrpahics ####
#subetting data - need to center age before adding interaction terms..
raw_subsets = list("Female" = 
                     cohort[gender_code == 0,][,gender_code := NULL],
                   "Male" = 
                     cohort[gender_code == 1,][,gender_code := NULL])
#center age and store results
base_age = center_age(raw_subsets) #MODIFIES INPLACE

# 3. Creating cohorts ####
# 3a. Comparison cohorts 
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

# 4. imputation ####
# dont include VSIMPLE index or outcomes in impuataion as one will crash it and other will bias
imputed_subsets = list()
imputed_subsets[["Female"]] = as.data.table(read.csv("data/missforest_imputed_subsets_female.csv"))
imputed_subsets[["Male"]] = as.data.table(read.csv("data/missforest_imputed_subsets_male.csv"))
# for (i in 1:length(raw_subsets)){
  
#   # copy for safety
#   dt = copy(raw_subsets[[i]])
  
#   # remove_index, outcome and follow-up info
#   dt_x = dt[,4:dim(dt)[2]]

#   # convert all ints to factors - all continuous are numerical here
#   fact_cols = sapply(dt_x,function(x) length(unique(x)) < 10)
#   dt_x[,.SD,.SD = fact_cols]

#   # impute - 
#   imputed_dt = missForest(dt_x,
#                           ntree = 10,
#                           maxiter = 10,
#                           verbose=TRUE)$ximp
  
#   #concatenate
#   imputed_dt = cbind(dt[,1:3],imputed_dt)

#   #list,name, append
#   list_dt = list(imputed_dt)
#   names(list_dt) = names(raw_subsets)[i]
#   imputed_subsets = append(imputed_subsets,list_dt)  
# }
#write.csv(imputed_subsets[["Female"]],"data/missforest_imputed_subsets_female.csv")
#write.csv(imputed_subsets[["Male"]],"data/missforest_imputed_subsets_male.csv")

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

###########################################################################################
# 5. fitting models with cross validation ####
###########################################################################################

# fitting and saving model predictions ####
iterations_pred_dt_list = list()
for (i in 1:2){
  # n_cv_iter CV repetitions and average predictions to remove bias from random sampling
  #create pred_dt base
  pred_dt = copy(subsets[[i]][,c("VSIMPLE_INDEX_MASTER","indicator","total_fu"),with=FALSE])
  n_cv_iter = 2

  # n CV repetitions and average predictions to remove bias from random sampling
  for (iteration in 1:n_cv_iter){
    print(paste0("iteration ",iteration," of ",n_cv_iter,Sys.time()))
    # comparison models
    pred_dt[,paste0('Baseline Model_', iteration) := 
              cv_coxph_fit_and_predict(comp_subsets[[i]],
                                      k_folds = 5,
                                      grouped_cv = FALSE,
                                      return_vector = TRUE)]
    
    # new model
    pred_dt[,paste0('New Model_', iteration) := 
              cv_coxph_fit_and_predict(subsets[[i]],
                                      k_folds = 5,
                                      grouped_cv = FALSE,
                                      return_vector = TRUE)]
  }
  
  #list, name, append
  list_pred_dt = list(pred_dt)
  names(list_pred_dt) = names(raw_subsets)[i]
  iterations_pred_dt_list = append(iterations_pred_dt_list,list_pred_dt)
}
iterations_pred_dt_list[[1]]

#average all predictions from the different 
#baseline and new model iterations for each individual
pred_dt_list = list()
for (i in 1:length(iterations_pred_dt_list)){
  temp = iterations_pred_dt_list[[i]][, .(
      VSIMPLE_INDEX_MASTER,
      indicator,  
      total_fu,
      `Baseline Model` = rowMeans(.SD[, grep("^Baseline Model", names(.SD)), with = FALSE]),
      `New Model` = rowMeans(.SD[, grep("^New Model", names(.SD)), with = FALSE])
  )]
  pred_dt_list = append(pred_dt_list,list(temp))
  names(pred_dt_list)[i] = names(iterations_pred_dt_list)[i]
}
pred_dt_list[[1]]

write.csv(pred_dt_list$Female,paste0(output_dir,"pred_dt_list_female_250503.csv"))
write.csv(pred_dt_list$Male,paste0(output_dir,"pred_dt_list_male_250503.csv"))

# output = summary_stats(pred_dt_list)
# write.csv(output,paste0(output_dir,"pred_summary_stats.csv"))


# repeat for lp
iterations_lp_dt_list = list()
for (i in 1:2){
  # n_cv_iterCV repetitions and average predictions to remove bias from random sampling
  #create pred_dt base
  pred_dt = copy(subsets[[i]][,c("VSIMPLE_INDEX_MASTER","indicator","total_fu"),with=FALSE])
  n_cv_iter = 2

  for (iteration in 1:n_cv_iter){
    print(paste0("iteration ",iteration," of ",n_cv_iter,Sys.time()))
    # comparison models
    pred_dt[,paste0('Baseline Model_', iteration) := 
              cv_coxph_fit_and_predict(comp_subsets[[i]],
                                      k_folds = 5,
                                      grouped_cv = FALSE,
                                      return_vector = TRUE,
                                      return_lp = TRUE)]
    
    # new model
    pred_dt[,paste0('New Model_', iteration) := 
              cv_coxph_fit_and_predict(subsets[[i]],
                                      k_folds = 5,
                                      grouped_cv = FALSE,
                                      return_vector = TRUE,
                                      return_lp = TRUE)]
  }
  
  #list, name, append
  list_pred_dt = list(pred_dt)
  names(list_pred_dt) = names(raw_subsets)[i]
  iterations_lp_dt_list = append(iterations_lp_dt_list,list_pred_dt)
}
iterations_lp_dt_list[[1]]

# #average all lps from the different 
# #baseline and new model iterations for each individual
# lp_dt_list = list()
# for (i in 1:length(iterations_lp_dt_list)){
#   temp = iterations_lp_dt_list[[i]][, .(
#       VSIMPLE_INDEX_MASTER,
#       indicator,  
#       total_fu,
#       `Baseline Model` = rowMeans(.SD[, grep("^Baseline Model", names(.SD)), with = FALSE]),
#       `New Model` = rowMeans(.SD[, grep("^New Model", names(.SD)), with = FALSE])
#   )]
#   lp_dt_list = append(lp_dt_list,list(temp))
#   names(lp_dt_list)[i] = names(iterations_pred_dt_list)[i]
# }
# lp_dt_list[[1]]

# output = summary_stats(lp_dt_list)
# write.csv(output,paste0(output_dir,"lp_summary_stats.csv"))

###############################################################################
# 6. global perfomance metrics new test ####
###############################################################################

# your data
pred_dt = copy(pred_dt_list$Female)[1:10000,]
# lp_dt = copy(iterations_lp_dt_list$Female)
# pred_dt$base_lp = lp_dt$`Baseline Model`
# pred_dt$new_lp = lp_dt$`New Model`


pred_dt$risk_base = pred_dt$`Baseline Model`
pred_dt$risk_new = pred_dt$`New Model`
# Your data.frame must have:
#   total_fu       : follow-up time in years
#   indicator      : 1 = event occurred, 0 = censored
#   survprob_base  : predicted S(5y) from baseline model
#   survprob_new   : predicted S(5y) from new model

# Define the 5-year horizon
t0 <- 5

# Precompute risk scores = 1 – survival probability
pred_dt$survprob_base <- 1 - (pred_dt$risk_base / 100)
pred_dt$survprob_new  <- 1 - (pred_dt$risk_new / 100)

#––– Bootstrap –––
B <- 100
n <- nrow(pred_dt)

ci_base    <- numeric(B)
ci_new     <- numeric(B)
r2_base    <- numeric(B)
r2_new     <- numeric(B)
royston_base <- numeric(B)
royston_new <- numeric(B)

set.seed(2025)
for (i in seq_len(B)) {
  print(paste0("bootstrapping ", i, " of ",B, Sys.time()))
  # 1) paired bootstrap resample
  idx <- sample.int(n, size = n, replace = TRUE)
  d   <- copy(pred_dt[idx, ])

  # 2) C-statistic & Somers’ D via Hmisc::rcorr.cens()
  print(paste0("C-statistic base",Sys.time()))
  rc_b       <- rcorr.cens(d$survprob_base, with(d, Surv(total_fu, indicator)))
  Dxy_b      <- unname(rc_b["Dxy"])
  ci_base[i] <- (Dxy_b + 1) / 2
  r2_base[i] <- Dxy_b^2

  print(paste0("C-statistic new",Sys.time()))
  rc_n       <- rcorr.cens(d$survprob_new,  with(d, Surv(total_fu, indicator)))
  Dxy_n      <- unname(rc_n["Dxy"])
  ci_new[i]  <- (Dxy_n + 1) / 2
  r2_new[i]  <- Dxy_n^2

  # # royston
  # print(paste0("royston base",Sys.time()))
  # base_cox = coxph(Surv(total_fu, indicator) ~ survprob_base, data = d)

  # print(paste0("royston new",Sys.time()))
  # new_cox = coxph(Surv(total_fu, indicator) ~ survprob_new, data = d)

  # royston_base[i] = royston(base_cox)["D"]
  # royston_new[i] = royston(new_cox)["D"]

  print(paste0("C-statistic base: ", ci_base[i], " R2 base: ", r2_base[i], " royston base: ", royston_base[i]))
  print(paste0("C-statistic new: ", ci_new[i], " R2 new: ", r2_new[i], " royston new: ", royston_new[i]))
}

base_cox = coxph(Surv(total_fu, indicator) ~ survprob_base, data = pred_dt)
new_cox = coxph(Surv(total_fu, indicator) ~ survprob_new, data = pred_dt)

#––– Original-data estimates –––
orig_rc_b    <- rcorr.cens(pred_dt$survprob_base, with(pred_dt, Surv(total_fu, indicator)))
orig_Dxy_b   <- unname(orig_rc_b["Dxy"])
orig_ci_b    <- (orig_Dxy_b + 1) / 2
orig_r2_b    <- orig_Dxy_b^2
orig_base_royston = royston(base_cox)["D"]

orig_rc_n    <- rcorr.cens(pred_dt$survprob_new, with(pred_dt, Surv(total_fu, indicator)))
orig_Dxy_n   <- unname(orig_rc_n["Dxy"])
orig_ci_n    <- (orig_Dxy_n + 1) / 2
orig_r2_n    <- orig_Dxy_n^2
orig_new_royston = royston(new_cox)["D"]

#––– Compute differences & summaries –––
ci_diff    <- ci_new    - ci_base
r2_diff    <- r2_new    - r2_base
royston_diff = royston_new - royston_base

summ <- function(x) {
  c(mean  = mean(x),
    "lower" = quantile(x, .025),
    "upper" = quantile(x, .975))
}

results <- data.frame(
  Metric            = c("C-statistic", "Explained variance (R²)", "Royston D"),
  Baseline_Mean     = c(summ(ci_base)["mean"],
                        summ(r2_base)["mean"],
                        summ(royston_base)["mean"]),
  Baseline_CI_lower = c(summ(ci_base)["lower.2.5%"],
                        summ(r2_base)["lower.2.5%"],
                        summ(royston_base)["lower.2.5%"]),
  Baseline_CI_upper = c(summ(ci_base)["upper.97.5%"],
                        summ(r2_base)["upper.97.5%"],
                        summ(royston_base)["upper.97.5%"]),
  New_Mean          = c(summ(ci_new)["mean"],
                        summ(r2_new)["mean"],
                        summ(royston_new)["mean"]),
  New_CI_lower      = c(summ(ci_new)["lower.2.5%"],
                        summ(r2_new)["lower.2.5%"],
                        summ(royston_new)["lower.2.5%"]),
  New_CI_upper      = c(summ(ci_new)["upper.97.5%"],
                        summ(r2_new)["upper.97.5%"],
                        summ(royston_new)["upper.97.5%"]),
  Diff_Mean         = c(summ(ci_diff)["mean"],
                        summ(r2_diff)["mean"],
                        summ(royston_diff)["mean"]),
  Diff_CI_lower     = c(summ(ci_diff)["lower.2.5%"],
                        summ(r2_diff)["lower.2.5%"],
                        summ(royston_diff)["lower.2.5%"]),
  Diff_CI_upper     = c(summ(ci_diff)["upper.97.5%"],
                        summ(r2_diff)["upper.97.5%"],
                        summ(royston_diff)["upper.97.5%"]),
  p_value           = c(
    2 * min(mean(ci_diff <= 0),    mean(ci_diff >= 0)),
    2 * min(mean(r2_diff <= 0),    mean(r2_diff >= 0)),
    2 * min(mean(royston_diff <= 0),    mean(royston_diff >= 0))
  )
)

print(results)






###############################################################################
# 6. % events in % population 250503 analysis ####
###############################################################################








################################################################################

#################################################################################

# # 11. Custom SUBSETS ####
# # 11a. flexible subseting
# gender_ind = "women" #men or women
# cond = "nhi_age >= 30 & nhi_age < 45 & eth_cat == 'NZM'"
# cond = "nhi_age >= 60 "
# cond = "nhi_age >= 30" #this will select all people
# subpop_index_temp = subpop_index(cond)
# #subsets[[1]][VSIMPLE_INDEX_MASTER %in% subpop_index_temp]
# subpop_pred_dt = pred_dt_list[[gender_ind]][
#   VSIMPLE_INDEX_MASTER %in% subpop_index_temp,]
# #plot using name of condition
# cal_plot = validation(subpop_pred_dt,
#                       type = "cal",
#                       title = cond)
# cal_plot

# pred = subpop_pred_dt[,"New Model"][[1]]
# obs = subpop_pred_dt[,"indicator"][[1]]

#####################################################################################
# Boostrapping
####################################################################################

# # Bootstrapping
# boot_output = pred_dt_bootstrap_wrapper(subpop_pred_dt)
# boot_output

################################################################################
# Re_distribution of predictions plot
##############################################################################

# # selecting which individuals to plot
# # random sample of n1 non events and n2 events or just a random sample of n individuals
# plot_dt = copy(pred_dt_list[["Female"]])
# sampled_1 <- plot_dt[indicator==1][sample(.N, 100), ]
# sampled_0 <- plot_dt[indicator==0][sample(.N, 250), ]  
# # Combine the sampled data back into one data table
# plot_dt <- rbind(sampled_1, sampled_0)
# #plot_dt = plot_dt[sample(.N,500),]
# plot_dt = plot_dt[order(plot_dt$indicator),]

# # Create the scatter plot with adjusted data order
# plot = ggplot(plot_dt, aes(x = `Baseline Model`, y = `New Model`, color = factor(indicator))) +
#   geom_point(aes(shape = factor(indicator)), size = 3, alpha = 0.6) +
#   scale_color_manual(values = c("0" = "blue", "1" = "red")) +
#   labs(title = "Comparison of Baseline and New Model Predictions",
#        x = "Baseline Model Prediction",
#        y = "New Model Prediction",
#        color = "Indicator",
#        shape = "Indicator") +
#   theme_minimal() +
#   theme(legend.position = "top") +
#   coord_fixed(ratio = 1)  # Ensures the plot has an aspect ratio of 1:1
# plot
# ggsave("re_distribution_plot.png",
#        plot = plot,
#        path = paste0(output_dir,"plots/misc/re_distribution_plots"),
#        height = 5,
#        width = 6.2)

# ##################################################################################

# ##################################################################################

# # 9. plotting ####  
# #while working so save scrolling up
# source("0_helper_functions.R")
# names(pred_dt_list) = c("Women","Men")
# # 9a. whole cohort calibration
# # plot calibration
# whole_cal_plot_list = whole_cohort_plot_generator(pred_dt_list,
#                                                   type = "calibration",
#                                                   display_n = TRUE,
#                                                   display_metrics = FALSE)

# # 9b. sub population plots ####
# full_cal_plot_list = cal_sub_pop_plot_generator(pred_dt_list,
#                                                 all_conds,
#                                                 display_n = TRUE,
#                                                 display_metrics = FALSE)

# # 10. Discrimination #### 
# # 8a. whole cohort discrimination
# whole_dis_plot_list = whole_cohort_plot_generator(pred_dt_list,
#                                                     type = "discrimination",
#                                                     display_n = TRUE,
#                                                     display_metrics = FALSE)


# # # 11b. arranging premade plots

# whole_pop_plots = ggarrange(plotlist = list(whole_cal_plot_list[[1]],
#                                             whole_cal_plot_list[[2]],
#                                             whole_dis_plot_list[[1]],
#                                             whole_dis_plot_list[[2]]),
#             ncol = 2,
#             nrow = 2)
# whole_pop_plots
# ggsave("custom_1_calibration_plot.png",
#        whole_pop_plots,
#        path = paste0(output_dir,"plots/calibration"),
#        height = 10,
#        width = 9)



# # names(full_cal_plot_list)
# # selected_plots = c("Maaori women",
# #                    "Maaori men",
# #                    "Women 60-74",
# #                    "Men 60-74",
# #                    "Women Deprivation Q5",
# #                    "Men Deprivation Q5")
# # selected_plot_list = full_cal_plot_list[selected_plots]
# # names(selected_plot_list)
# # custom_plot = ggarrange(plotlist = selected_plot_list,
# #                         ncol = 2,
# #                         nrow = 3)
# # custom_plot
# # ggsave("custom_calibration_plot.png",
# #        custom_plot,
# #        path = "plots/validation_plots/calibration",
# #        height = 15,
# #        width = 9)

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
#        path = paste0(output_dir,"plots/calibration"),
#        height = 15,
#        width = 9)


# # 13 Net reclassification  ####
# reclassification_dt = list(reclassification_testing(pred_dt_list[["Women"]]),
#                           reclassification_testing(pred_dt_list[["Men"]]))
