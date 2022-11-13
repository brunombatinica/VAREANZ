#CC analysis using pre 2015 data cut-off

# Clearing workspace
rm(list = ls())
dev.off()

#set working directory
setwd("//uoa.auckland.ac.nz/Shared/MED/EPBI/ViewData/Users/bbat644/Desktop/code")

#source helper functions (also loads in libraries)
source("0_helper_functions.R")


#1. import data and cleaning up data ####
an_cohort_ts = read.fst("data/cohorts/an_cohort_ts.fst",
                     as.data.table = TRUE)
cohort = an_cohort_ts[complete.cases(an_cohort_ts)]
cohort[,DHB_name := NULL]
dim(cohort)
#cut up factors
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
}
names(base_age) = names(raw_subsets)




# 3. Creating cohorts *** ######################

# 3a. Comparison cohorts ####
# add prespecified interactions and removes test safe variables
TS_vars = c("egfr",
            "c_egfr",
            "egfr_cat_G3a",
            "egfr_cat_G3b",
            "hba1c",
            "c_hba1c",
            "tchdl",
            "c_tchdl",
            "tri",
            "c_tri",
            "log_tri")
comp_subsets = list()
for (i in 1:length(raw_subsets)){
  
  #copy data table for safety
  dt = copy(raw_subsets[[i]])
  
  #remove uncentered age
  dt[,nhi_age := NULL]
  
  #add suneelas interactions
  dt = prespec_ints(dt)
  
  #remove test variables
  dt = dt[,(TS_vars) := NULL]
  
  #list,name, append
  list_dt = list(dt)
  names(list_dt) = paste(names(raw_subsets)[i],"comparison",sep = " ")
  comp_subsets = append(comp_subsets,list_dt)
}


# 3b. New cohort
# olds ints + test safe + new ints
subsets = list()
for (i in 1:length(raw_subsets)){
  
  #copy for safety
  dt = copy(raw_subsets[[i]])
  
  #remove uncentered age
  dt[,nhi_age := NULL]
  
  #add old interactions
  dt = prespec_ints(dt)
  
  # transform hba1c and tchdl
  dt = transform_labs(dt,
                   c_hba1c = 1,
                   c_tchdl = 1,
                   c_tri = 1)
  
  # # #add log triglycerids
  # dt = dt[,":="(log_tri = log2(tri),
  #              tri = NULL)]

  
  #cut up egfr
  dt[,":="(egfr_cat = egfr_groups(egfr),
           egfr = NULL)]
  dt = cut_factors(dt)
  
  
  #new interactions
  dt = dt[,":="(age_tchdl = c_nhi_age*c_tchdl)]
  
  #list,name, append
  list_dt = list(dt)
  names(list_dt) = names(raw_subsets)[i]
  subsets = append(subsets,list_dt)
  
}

# 4. Fitting models
# 4a. Comparison models
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

# 4b. models 
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

# 5. Demographics tables ####
#simply run the demographics table function on the male and female subsets
source("0_helper_functions.R")
subset_list = list("Women" = list(subsets[[1]][,VSIMPLE_INDEX_MASTER]),
                   "Men" = list(subsets[[2]][,VSIMPLE_INDEX_MASTER]))
demo_table = demo_table_generator(subset_list,
                                  test_rows = TRUE)
demo_table
write.csv(demo_table,
          file = "tables/demo_table_cc_gender.csv")


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



# 7. model coefficients

# 7a. coefficients plot
source("0_helper_functions.R")
model_list_list = copy(all_models)
coef_plot = coef_plot_generator(model_list_list,
                             return = 1)
coef_plot
ggsave("coef_plot_cc.png",
       coef_plot,
       path = "plots/misc",
       height = 7,
       width = 9)

#alternate background colurs - grey and white to get different colours

# 7b. HR table ####
hr_table = coef_table_generator(model_list_list)
hr_table
write.csv(hr_table,
          "tables/hr_table_cc.csv")

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
          "tables/coef_p_plot_cc.R")


# 8. summary statistics table ####
#add all models togther
flat_model_list = append(models,comp_models)[c(1,3,2,4)]
summary_stats = c("baseline survival",
                  "Harrels C",
                  "Royston D",
                  "Royston R^2")
stats_table = data.frame(matrix(NA,
                                nrow = length(summary_stats),
                                ncol = length(flat_model_list)))
colnames(stats_table) = names(flat_model_list)
#rownames(stats_table) = c(names(model1),names(model2))
rownames(stats_table) = summary_stats

for (i in 1:length(flat_model_list)){
  
  bs = base_surv(flat_model_list[[i]])
  stats_table[1,i] = sprintf("%.4f",bs)
  
  model_c = summary(flat_model_list[[i]])$concordance
  string_c = sprintf("%.4f (%.4f - %.4f)",
                     model_c[1],
                     model_c[1] - 1.96*model_c[2],
                     model_c[1] + 1.96*model_c[2])
  stats_table[2,i] = string_c
  
  roy = royston(flat_model_list[[i]])
  stats_table[3,i] = sprintf("%.4f (%.4f - %.4f)",
                             roy[1],
                             roy[1] - 1.96*roy[2],
                             roy[1] + 1.96*roy[2])
  stats_table[4,i] = sprintf("%.4f",roy[3])
  
}

stats_table
write.csv(stats_table,
         "tables/stats_table_cc.csv")

# 9. Calibration ####
  
#while working so save scrolling up
source("0_helper_functions.R")

# 9a. whole cohort calibration
#plot calibration
whole_cal_list = list()
for (i in 1:2){
  #whole cohort calibration plots
  cal_plot = validation(pred_dt_list[[i]],
                        type = "cal",
                        str_to_title(names(pred_dt_list)[[i]]),
                        return = 1)
  cal_plot
  
  whole_cal_list = append(whole_cal_list,list(cal_plot))
  
  ggsave(paste(names(pred_dt_list)[[i]],".png",sep=""),
         plot = cal_plot,
         path = "plots/validation_plots/calibration",
         height = 5,
         width = 6.2)
}

whole_cal_plot = ggarrange(plotlist = whole_cal_list,
            ncol = 1)
whole_cal_plot
ggsave("whole_calibration_plot.png",
       plot = whole_cal_plot,
       path = "plots/validation_plots/calibration",
       height = 10,
       width = 6)


# 9b. sub population plots ####
source("0_helper_functions.R")

sub_pop_demo_table_list = list()
full_cal_plot_list = list()
for (i in 1:length(all_conds) ){
  
  conds_name = names(all_conds)[i]
  conds = all_conds[[i]]
  
  index_list = list()
  
  for (gender_ind in 1:length(pred_dt_list)){
    
    gender = names(pred_dt_list)[gender_ind]
    
    cal_plot_list = list()
    
    for (cond_ind in 1:length(conds)){
      
      # 
      subpop_index_temp = subpop_index(conds[[cond_ind]])
      
      
      #subset the relevant pred_dt
      sub_pop_pred_dt = pred_dt_list[[gender_ind]][
        VSIMPLE_INDEX_MASTER %in% subpop_index_temp,]
      
      #extract plot name
      plot_name = eval(parse(text = names(conds)[cond_ind]))
      
      #save to index list
      subpop_index_temp_list = list(subpop_index_temp)
      names(subpop_index_temp_list) = plot_name
      index_list = append(index_list,subpop_index_temp_list)
      
      #plot using name of condition
      cal_plot = validation(sub_pop_pred_dt,
                            type = "cal",
                            title = plot_name)
      
      list_cal_plot = list(cal_plot)
      names(list_cal_plot) = plot_name
      
      
      #append to calplot list ot feed into ggarange and to save for special plotting
      cal_plot_list = append(cal_plot_list,list_cal_plot)
      full_cal_plot_list = append(full_cal_plot_list,list_cal_plot)
      
    }
    
    arrange_plot = ggarrange(plotlist = cal_plot_list,
                             ncol = 2,
                             nrow = ceiling(length(cal_plot_list)/2))
    
    
    arrange_plot
    
    # demo_table = demo_table_generator(index_list,test_rows = TRUE)
    # write.csv(demo_table,
    #           file = paste("tables/cc_subpop_demo/",gender,conds_name,".csv",sep = ""))
    
    ggsave(paste(conds_name,gender,".png",sep = ""),
           arrange_plot,
           path = "plots/validation_plots/calibration/subpopulations",
           height = ceiling(length(cal_plot_list)/2)*5,
           width = 9)
  }
  
}

# 10. Discrimination #### 
source("0_helper_functions.R")
# 8a. whole cohort discrimination

whole_dis_list = list()
for (i in 1:2){
  #whole cohort calibration plots
  dis_plot = validation(pred_dt_list[[i]],
                        type = "dis",
                        str_to_title(names(pred_dt_list)[[i]]),
                        return = 1)
  dis_plot
  
  whole_dis_list = append(whole_dis_list,list(dis_plot))
  
  ggsave(paste(names(pred_dt_list)[[i]],".png",sep=""),
         plot = dis_plot,
         path = "plots/validation_plots/discrimination",
         height = 5,
         width = 6.2)
}

whole_dis_plot = ggarrange(plotlist = whole_dis_list,
                           ncol = 1)
whole_dis_plot
ggsave("whole_discrimination_plot.png",
       plot = whole_dis_plot,
       path = "plots/validation_plots/discrimination",
       height = 10,
       width = 8)
# 
# # 10b. discrimiantion in sub population plots ####
# source("0_helper_functions.R")
# 
# full_dis_plot_list = list()
# 
# for (i in 1:length(all_conds) ){
#   
#   conds_name = names(all_conds)[i]
#   conds = all_conds[[i]]
#   
#   for (gender_ind in 1:length(pred_dt_list)){
#     
#     gender = names(pred_dt_list)[gender_ind]
#     
#     dis_plot_list = list()
#     
#     for (cond_ind in 1:length(conds)){
#       
#       #subset the relevant pred_dt
#       sub_pop_pred_dt = pred_dt_list[[gender_ind]][
#         VSIMPLE_INDEX_MASTER %in% subpop_index(conds[[cond_ind]]),]
#       
#       plot_name = eval(parse(text = names(conds)[cond_ind]))
#       
#       #plot using name of condition
#       dis_plot = validation(sub_pop_pred_dt,
#                             type = "d",
#                             title = plot_name)
#       
#       #list, anme append
#       list_dis_plot = list(dis_plot)
#       names(list_dis_plot) = plot_name
#       full_dis_plot_list = append(full_dis_plot_list,list_dis_plot)
#       
#       #append to calplot list ot feed into ggarange
#       dis_plot_list = append(dis_plot_list,list_dis_plot)
#       
#     }
#     
#     arrange_plot = ggarrange(plotlist = dis_plot_list,
#                              ncol = 2,
#                              nrow = ceiling(length(dis_plot_list)/2))
#     
#     # ggsave(paste(conds_name,gender,".png",sep = ""),
#     #        arrange_plot,
#     #        path = "plots/validation_plots/discrimination/subpopulations",
#     #        height = ceiling(length(dis_plot_list)/2)*5,
#     #        width = 9)
#   }
#   
# }

# 11. Custom plots ####
names(full_cal_plot_list)
selected_plots = c("Maaori women",
                   "Maaori men",
                   "Women 60-74",
                   "Men 60-74",
                   "Women Deprivation Q5",
                   "Men Deprivation Q5")
selected_plot_list = full_cal_plot_list[selected_plots]
names(selected_plot_list)
custom_plot = ggarrange(plotlist = selected_plot_list,
                        ncol = 2,
                        nrow = 3)
custom_plot
ggsave("custom_calibration_plot.png",
       custom_plot,
       path = "plots/validation_plots/calibration",
       height = 15,
       width = 9)

#selected demographics tables

# 12 PI quantification ####
source("0_helper_functions.R")
pi_cohort_list = list()

#create pi cohort
for (i in 1:2){
  temp_cohort = copy(subsets[[i]])
    
  pi_vector = predict(comp_models[[i]],
                      temp_cohort,
                      type = "lp",
                      reference = "zero")

  
  #append pi vector and leave only new variables
  pi_cohort = temp_cohort[,pi := pi_vector]
  pi_cohort = pi_cohort[,c("VSIMPLE_INDEX_MASTER",
                           "total_fu",
                           "indicator",
                           "pi",
                           "egfr_cat_G3a",
                           "egfr_cat_G3b",
                           "c_hba1c",
                           "c_tchdl",
                           "c_tri",
                           "age_tchdl")]
  pi_cohort_list = append(pi_cohort_list,list(pi_cohort))
}

pi_cohort_list

#create pi models
pi_model_list = list()

for (i in 1:2){
  
  temp_cohort = pi_cohort_list[[i]]
  
  pi_model = coxph(Surv(total_fu,indicator) ~ .,
                   data = temp_cohort[,-1],
                   ties = "breslow",
                   model = TRUE)
  
  #list, name, append
  list_pi_model = list(pi_model)
  names(list_pi_model) = names(subsets)[i]
  pi_model_list = append(pi_model_list,list_pi_model)
  
}

#coefficient table #
pi_coef_table = coef_table_generator(list(pi_model_list))
pi_coef_table
write.csv(pi_coef_table,
          "tables/pi_coef_table.csv")




# 13 Net reclassification  ####
dt = pred_dt_list[[2]]

sum(dt$indicator)

old_high_risk = dt[`Old model` > 15]
new_high_risk = dt[`New model` > 15]

dim(old_high_risk)
sum(old_high_risk$indicator)
mean(old_high_risk$indicator)

dim(new_high_risk)
sum(new_high_risk$indicator)
mean(new_high_risk$indicator)

high_gained = new_high_risk[!(VSIMPLE_INDEX_MASTER %in% old_high_risk$VSIMPLE_INDEX_MASTER)]
high_lost = old_high_risk[!(VSIMPLE_INDEX_MASTER %in% new_high_risk$VSIMPLE_INDEX_MASTER)]

mean(high_gained$indicator)
mean(high_lost$indicator)

sum(high_gained$indicator)
sum(high_lost$indicator)

old_low_risk = dt[`Old model` < 5]
new_low_risk = dt[`New model` < 5]

dim(old_low_risk)
sum(old_low_risk$indicator)
mean(old_low_risk$indicator)

dim(new_low_risk)
sum(new_low_risk$indicator)
mean(new_low_risk$indicator)

low_gained = new_low_risk[!(VSIMPLE_INDEX_MASTER %in% old_low_risk$VSIMPLE_INDEX_MASTER)]
low_lost = old_low_risk[!(VSIMPLE_INDEX_MASTER %in% new_low_risk$VSIMPLE_INDEX_MASTER)]

mean(low_gained$indicator)
mean(low_lost$indicator)

sum(low_gained$indicator)
sum(low_lost$indicator)
