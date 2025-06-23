#imputed analysis using pre 2013 data_cutoff

# Clearing workspace
rm(list = ls())
dev.off()

#set working directory
setwd("C:/Users/bruno/OneDrive/Documents/Code/projects/VAREANZ/Honours")

#source helper functions (also loads in libraries)
source("0_helper_functions.R")

#output directory
output_dir = "outputs/rf_imputed_analysis/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir,recursive = TRUE) 
}

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

#sanity check regarding imputations ####
dt = copy(raw_subsets[["Female"]][
   , .SD, .SDcols = c("egfr","hba1c","tchdl","tri")])
dt_imputed = copy(imputed_subsets[["Female"]][ 
   , .SD, .SDcols = c("egfr","hba1c","tchdl","tri")])
#imputed_positions = dt[, lapply(.SD, is.na)]
dt_long = melt(dt, measure.vars = names(dt), variable.name = "Variable", value.name = "Original")
imputed_position = is.na(dt_long[,.SD,.SDcols = "Original"])
dt_imputed_long <- melt(dt_imputed, measure.vars = names(dt_imputed), variable.name = "Variable", value.name = "Imputed")
dt_imputed_long[,imputed_position := imputed_position]
dt_imputed_long[imputed_position == FALSE, Imputed := NA]
dt_imputed_long

x_titles = c(
  egfr="mL/min/1.73m^2",
  hba1c="mmol/L",
  tchdl="mmol/L",
  tri="mmol/L"
)
n_present <- dt_long %>% 
  drop_na(Original) %>%            # keep only rows that were observed
  group_by(Variable) %>% 
  summarise(n_present = n(), .groups = "drop")
n_present
n_imputed <- dt_imputed_long %>% 
  drop_na(Imputed) %>% 
  group_by(Variable) %>% 
  summarise(n_imputed = n(), .groups = "drop")
n_imputed
n_tbl <- full_join(n_present, n_imputed, by = "Variable") %>% 
  replace_na(list(n_present = 0, n_imputed = 0))   # just in case
n_tbl
facet_labs <- n_tbl %>% 
  mutate(label = paste0(
           x_titles[Variable],                       # 1st line
           "\nnumber present: (",  n_present, ")", # 2nd line, part 1
           "\nnumber imputed: (", n_imputed, ")"       # 2nd line, part 2
         )) %>% 
  select(Variable, label) %>%
  tibble::deframe() 
facet_labs

x_titles
imputation_plot = ggplot() +
  geom_histogram(data = dt_long, 
                 aes(x = Original, fill = "Original"), 
                 binwidth = 1, 
                 alpha = 0.5) +
  geom_histogram(data = dt_imputed_long, 
                 aes(x = Imputed, fill = "Imputed"), 
                 binwidth = 1, 
                 alpha = 0.5) +
  facet_wrap(~Variable, 
            scales = "free",
            labeller = labeller(Variable = facet_labs),
            strip.position = "bottom") +
  labs(x = NULL,
       y = "Frequency",
       fill = "Data Type") +  # Legend title
  scale_fill_manual(values = c("Original" = "blue", "Imputed" = "red")) +
  theme_minimal(base_size = 18)+
  theme(legend.position = "bottom",
        strip.placement = "outside"
        )
imputation_plot
ggsave("imputation_sanity_check.png",
       imputation_plot,
       path = paste0(output_dir,"plots/misc"),
       height = 7,
       width = 9)

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

#########################################################################################
# MODEL VALIDATION #####
#########################################################################################

# 6. create prediction dts ####
#reformat models into embedded list
female_models = list("New Model" = copy(models[[1]]),
                     "Baseline Model" = copy(comp_models[[1]]))
male_models = list("New Model" = copy(models[[2]]),
                   "Baseline Model" = copy(comp_models[[2]]))
all_models = list("women" = female_models,
                  "men" = male_models)

female_data = list("New Model" = subsets[[1]],
                   "Baseline Model" = comp_subsets[[1]])
male_data = list("New Model" = subsets[[2]],
                   "Baseline Model" = comp_subsets[[2]])
all_data = list("women" = female_data,
                "men" = male_data)

pred_dt_list = list()
lp_dt_list = list()
for (i in 1:2){
  model_list = all_models[[i]]
  new_data = copy(subsets[[i]])
  lp_dt = pred_dt_generator(model_list,
                          new_data,
                          return_lp = TRUE)
  pred_dt = pred_dt_generator(model_list,
                             new_data)

  #list, name, append
  list_pred_dt = list(pred_dt)
  names(list_pred_dt) = names(all_models)[i]
  pred_dt_list = append(pred_dt_list,list_pred_dt)

  #list, name, append
  list_lp_dt = list(lp_dt)
  names(list_lp_dt) = names(all_models)[i]
  lp_dt_list = append(lp_dt_list,list_lp_dt)
}

test1 = predict(models[[1]],type = "lp",data = subsets[[1]])
test2 = predict(models[[1]],type = "lp",reference = "zero",data = subsets[[2]])



# 7. Demographics tables ####
#simply run the demographics table function on the male and female subsets
source("0_helper_functions.R")
subset_list = list("Women" = list(subsets[[1]][,VSIMPLE_INDEX_MASTER]),
                   "Men" = list(subsets[[2]][,VSIMPLE_INDEX_MASTER]))
demo_table = demo_table_generator(subset_list,
                                  test_rows = TRUE)
demo_table
write.csv(demo_table,
          file = paste0(output_dir,"tables/demo_table_gender.csv"))

# 8. model coefficients

# 9a. coefficients plot ####
source("0_helper_functions.R")
model_list_list = copy(all_models)
coef_plot = coef_plot_generator(model_list_list,
                             return = 1)
coef_plot
ggsave("coef_plot.jpg",
       coef_plot,
       path = paste0(output_dir,"plots/misc"),
       height = 7,
       width = 9)

#alternate background colurs - grey and white to get different colours

# 9b. HR table ####
hr_table = coef_table_generator(model_list_list)
hr_table
write.csv(hr_table,
          paste0(output_dir,"tables/hr_table_imputed_rf.csv"))

#7c Coefficient and P value table
pred = names(summary(models[[1]])$coefficients)
f = as.data.table(summary(models[[1]])$coefficients[,c("coef","Pr(>|z|)")],
                  keep.rownames = TRUE)
m = as.data.table(summary(models[[2]])$coefficients[,c("coef","Pr(>|z|)")],
                  keep.rownames = TRUE)
f[,2:3] = round(f[,2:3],4)
m[,2:3] = round(m[,2:3],4)
output = f[m, on = "rn"]
output
write.csv(output,
          paste0(output_dir,"tables/coef_p_plot.csv"))

# 8. Linear predictor
lp_model_list = list()
for (i in 1:2){
  model = copy(comp_models[[i]])
  data = copy(subsets[[i]])
  lp = predict(model,type = "lp",data = data)
  data = cbind(data,lp)
  #remove all variables except indicator, total_fu, 
  # and any column that contains egfr, hba1c, tchdl, tri
  patterns <- c("indicator","total_fu","lp","egfr", "hba1c", "tchdl", "tri")
  # Select columns that contain any of the patterns or match exact names
  data = data[, .SD, .SDcols = c( 
  names(data)[
    unique(unlist(lapply(patterns, function(p) grep(p, names(data), value = FALSE))))])
  ]
  lp_model = coxph(Surv(total_fu,indicator) ~ .,
                  data = data,
                  ties = "breslow",
                  model = TRUE)
  list_lp_model = list(lp_model)
  names(list_lp_model) = paste0(names(subsets)[i])
  lp_model_list = append(lp_model_list,list_lp_model)
}
temp = list("lp" = lp_model_list, "ignore" = lp_model_list)
lp_table = coef_table_generator(temp)
write.csv(lp_table,
          paste0(output_dir,"tables/coef_p_plot_lp.csv"))

#################################################################################

#################################################################################

###### testing manual stats table ####
stats_table_list = list()
gender_ind = 1 #gender


summary_stats(all_models,all_data)




model$concordance
Surv_ = Surv(time = data$total_fu,event = as.numeric(data$indicator))
c_index = concordance.index(prob[1:1000],data$total_fu[1:1000],data$indicator[1:1000])
string_c = sprintf("%.5f (%.5f - %.5f)",
                     c_index$c.index,
                     c_index$lower,
                     c_index$upper)

coxph(Surv(total_fu,indicator) ~ lp,
      data = data,
      model = TRUE)
sd(lp)



# 10. summary statistics table  ####
#add all models togther
flat_model_list = append(models,comp_models)[c(1,3,2,4)]
summary_stats = c("baseline survival",
                  "Harrels C",
                  "Royston D",
                  "Royston D manual"
                  "Royston R^2")

stats_table = data.frame(matrix(NA,
                                nrow = length(summary_stats),
                                ncol = length(flat_model_list)))
colnames(stats_table) = names(flat_model_list)
#rownames(stats_table) = c(names(model1),names(model2))
rownames(stats_table) = summary_stats
roy_list = list()
for (i in 1:length(flat_model_list)){
  
  bs = base_surv(flat_model_list[[i]])
  stats_table[1,i] = sprintf("%.4f",bs)
  #Harrels C
  model_c = summary(flat_model_list[[i]])$concordance
  print(model_c)
  string_c = sprintf("%.4f (%.4f - %.4f)",
                     model_c[1],
                     model_c[1] - 1.96*model_c[2],
                     model_c[1] + 1.96*model_c[2])
  stats_table[2,i] = string_c

  # Royston D
  roy = royston(flat_model_list[[i]])
  print(roy)
  list_roy = list(roy)
  names(list_roy) = names(flat_model_list)[i]
  roy_list = append(roy_list,list_roy)
  D = roy[1]
  D_se = roy[2]
  R2 = roy[3]
  R2_manual = DtoR2(D)
  R2_ul = DtoR2(D + 1.96*D_se)
  R2_ll = DtoR2(D - 1.96*D_se)



  
  stats_table[3,i] = sprintf("%.4f (%.4f - %.4f)",
                             D,
                             D - 1.96*D_se,
                             D + 1.96*D_se)
  stats_table[4,i] = sprintf("%.4f (%.4f - %.4f)",
                             R2_manual,
                             R2_ll,
                             R2_ul)
  stats_table[5,i] = sprintf("%.4f",roy[3])
}




system.time(royston(all_models[[1]][[1]]))

library(survMisc)
rsq(all_models[[1]][[1]])
summmary(all_models[[1]][[1]])
library(rms)
rsm(all_models[[1]][[1]])
roy_list = append(roy_list,list(roy))
stats_table
write.csv(stats_table,
         paste0(output_dir,"tables/stats_table_cc.csv"))



################################################################################

#################################################################################

# 11. Custom SUBSETS ####
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
# Boostrapping
####################################################################################
# Bootstrapping
boot_output = pred_dt_bootstrap_wrapper(subpop_pred_dt)
#model_boot_wrapper(subsets[[1]],comp_subsets[[1]])
#   #calculate R^2 statistic
# }

# postResample(pred_dt[,"Baseline Model"][[1]]/100,pred_dt[,"indicator"][[1]])
# system.time(royston(model))

################################################################################
# Re_distribution of predictions plot
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
# PLOTTING ####
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

names(full_cal_plot_list)
selected_plots = c("Women 30-44",
                   "Men 30-44",
                   "Women 45-59",
                   "Men 45-59",
                   "Women 60-74",
                   "Men 60-74")
selected_plot_list = full_cal_plot_list[selected_plots]
names(selected_plot_list)
custom_plot = ggarrange(plotlist = selected_plot_list,
                        ncol = 2,
                        nrow = 3)
custom_plot
ggsave("custom_2_calibration_plot.png",
       custom_plot,
       path = paste0(output_dir,"plots/calibration"),
       height = 15,
       width = 9)


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
reclassification_dt = list(reclassification_testing(pred_dt_list[["women"]]),
                          reclassification_testing(pred_dt_list[["men"]]))
write.csv(reclassification_dt[[1]][[1]],
          paste0(output_dir,"tables/reclassification_table_women.csv"))
write.csv(reclassification_dt[[2]][[1]],
          paste0(output_dir,"tables/reclassification_table_men.csv"))



# misc
print(names(comp_subsets[[1]]))
temp = copy(comp_subsets[[1]])
temp[
  ph_bp_lowering_prior_6mths == 1 & 
  ph_lipid_lowering_prior_6mths == 1 & 
  ph_antithrombotic_prior_6mths == 1]
mean(temp$any_med)

# 
mean(cohort[cohort$gender_code == 1,]$tri, na.rm = TRUE)


######################################################################
# proportion captures
####################################################################
library(data.table)
library(ggplot2)
library(scales)
library(ggpubr)     # ggarrange()

#–––––––––––––––––––––#
# Helper to summarise #
#–––––––––––––––––––––#
make_summary <- function(dt, risk_col = "New Model",
                         thr_seq = seq(0, 10, by = 0.1)) {
  n_tot   <- nrow(dt)
  n_ev    <- dt[indicator == 1, .N]

  rbindlist(lapply(thr_seq, function(th) {
    data.table(
      threshold  = th,
      metric     = c("Screen positive", "Events captured"),
      proportion = c(dt[get(risk_col) > th, .N] / n_tot,
                     if (n_ev > 0)
                       dt[get(risk_col) > th & indicator == 1, .N] / n_ev
                     else NA_real_)
    )
  }))
}

#––––––––––––––––#
# Build figures  #
#––––––––––––––––#
plot_threshold_curve <- function(summary_dt, panel_title) {
  ggplot(summary_dt,
         aes(threshold, proportion, colour = metric)) +
    geom_line(linewidth = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(title = panel_title,
         x = "Risk threshold",
         y = "Proportion",
         colour = NULL) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
}

# generate summaries
sum_women <- make_summary(pred_dt_list[["women"]])
sum_men   <- make_summary(pred_dt_list[["men"]])

# create separate plots
p_women <- plot_threshold_curve(sum_women, "Women")
p_men   <- plot_threshold_curve(sum_men,   "Men")

#–––––––––––––––––––––––––#
# Arrange side-by-side    #
#–––––––––––––––––––––––––#
ggarrange(
  p_men, p_women,
  # labels        = c("Men", "Women"),   # panel tags
  ncol          = 2,
  common.legend = TRUE,
  legend        = "bottom"
)



