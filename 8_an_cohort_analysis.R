#Sub model analysis

# Clearing workspace
rm(list = ls())
dev.off()

#set working directory
setwd("//uoa.auckland.ac.nz/Shared/MED/EPBI/ViewData/Users/bbat644/Desktop/code")

# Loading in necessary libraries and data ###############################
source("0_helper_functions.R")

#import data and cleaning up data ####
#cc_cohort_2008_2015 = read.fst(path = "//uoa.auckland.ac.nz/Shared/MED/EPBI/ViewData/Users/bbat644/Desktop/R Code/CC_Cohort_2008_2015.fst")
an_cohort_ts = read.fst("data/cohorts/an_cohort_ts.fst",
                             as.data.table = TRUE)

cohort = data.table(an_cohort_ts)
names(cohort)
cohort[,c("DHB_name") := NULL]
#cohort = data.table(cc_cohort_2008_2015)
setkey(cohort,VSIMPLE_INDEX_MASTER)

# Initial Subsetting ####
#subetting data - need to center age before adding interaction terms..
raw_subsets = list("women" = 
                     cohort[gender_code == 0,][,gender_code := NULL],
                   "men" = 
                     cohort[gender_code == 1,][,gender_code := NULL])

#center age and store results
base_age = c()
for (raw_subset in raw_subsets){
  base_age = c(base_age,raw_subset[,mean(nhi_age)])
  raw_subset[,c_nhi_age := (nhi_age - mean(nhi_age))]
  #raw_subset[,egfr := NA]
}
names(base_age) = names(raw_subsets)


#*** SUBSET #########################
#Creating subsets #####
subsets = list()
#CREATING SUBSETS WITH ALL THE OLD AND NEW INTERACTIONS!
names = c("Female",
          "Male")
for (i in 1:length(raw_subsets)){
  temp_subset = copy(raw_subsets[[i]])
  
  #remove uncentered age
  temp_subset[,nhi_age := NULL]
  
  #centering labs
  temp_subset = transform_labs(temp_subset,
                                 c_hba1c = 1,
                                 c_tchdl = 1,
                               c_tri = 1)
  
  #group egfr
  temp_subset[,egfr_cat := egfr_groups(egfr)]
  temp_subset[,egfr := NULL]
  
  # #log triglycerides
  # #Mmanually adjust 2 results of 0
  # temp_subset[tri ==0,tri:= 0.1]
  # temp_subset[,log_tri := log2(tri)]
  # temp_subset[,tri := NULL]
  
  #cut factors
  temp_subset = cut_factors(temp_subset)
  
  #add old interactions
  temp_subset = prespec_ints(temp_subset)
  
  #add age_tchdl interaction
  temp_subset = temp_subset[,age_tchdl := c_nhi_age*c_tchdl]
  
  #create list, rename and append
  temp_subset = list(temp_subset)
  names(temp_subset) = names[i]
  subsets = append(subsets,temp_subset)
}

# ****PATTERN SUBMODELS ###############################################

#helper functions ####
# I will label the possible subsets using subset notation - and labelled 
# following binary interpretation of the subset in the order
# { (egfr,1) , (hba1c,2) , (tchdl,3) }
# e.g. if a patient has hba1c and tchdl present = [1,1,0] = 110 = 6
# slght rearranging required so just be careful
# binary relevace

ps_index = function(egfr,hba1c,tchdl,tri){
  # parses 3 inputs missingess to give a binary represenatation of
  # the specific submodel
  # IMPORTANT TO ENTER THE COLUMNS IN THE RIGHT ORDER
  
  return(2^0*as.integer(!is.na(egfr)) +
           2^1*as.integer(!is.na(hba1c)) + 
            2^2*as.integer(!(is.na(tchdl)|is.na(tri))))
}

# cc cohort helper
# which patterns can be used to train a given CC model
man_ccs = list(c(0,1,2,3,4,5,6,7), #0
               c(1,3,5,7), #1 egfr
               c(2,3,6,7), #2 hbac
               c(3,7), #3 egfr +hba1c
               c(4,5,6,7), #4 tchdl
               c(5,7), #5 egfr +tchdl
               c(6,7), #6 hba1c +tchdl
               c(7)) #cc


# actual analaysis ####
missing_index = copy(an_cohort_ts)
#allocating pattern index as specified above
missing_index[,pattern := ps_index(egfr,hba1c,tchdl,tri)]
table(missing_index$pattern)

#creating index cohorts
#pattern submodels cohorts ####
ps_cohorts = list()
ps_indexes = list()
for (gender in 1:2){
  temp_subsets = list()
  temp_indexes = list()
  for (j in 0:7){
    pattern_indexes = missing_index[pattern == j,VSIMPLE_INDEX_MASTER]
    
    temp_subset = subsets[[gender]][
      VSIMPLE_INDEX_MASTER %in% pattern_indexes]
    
    #remove NA columns and turn to list
    temp_subset = temp_subset[
      ,which(colSums(is.na(temp_subset)) != 0) := NULL]
    
    #list,name append cohorts
    list_temp_subset = list(temp_subset)
    names(list_temp_subset) = as.character(j)
    temp_subsets = append(temp_subsets,list_temp_subset)
    
    #
    list_indexes = list(temp_subset$VSIMPLE_INDEX_MASTER)
    names(list_indexes) = as.character(j)
    temp_indexes = append(temp_indexes,list_indexes)
  }
  
  #list, name, append
  temp_subsets = list(temp_subsets)
  names(temp_subsets) = names(raw_subsets)[gender]
  ps_cohorts = append(ps_cohorts,temp_subsets)
  
  temp_indexes = list(temp_indexes)
  names(temp_indexes) = names(raw_subsets)[gender]
  ps_indexes = append(ps_indexes,temp_indexes)
}

#sanity check
for (i in ps_cohorts[[1]]){ print(dim(i))}

#complete case cohorts ####
ccs_cohorts = list()
for (gender in 1:2){
  temp_subsets = list()
  for (j in 0:7){
    #extracting relevant complete case cohorts
    pattern_indexes = missing_index[pattern %in% man_ccs[[j+1]],
                                      VSIMPLE_INDEX_MASTER]
    temp_subset = subsets[[gender]][
      VSIMPLE_INDEX_MASTER %in% pattern_indexes]
    
    #remove NA columns and turn to list
    temp_subset = temp_subset[
      ,which(colSums(is.na(temp_subset)) != 0) := NULL]
    
    #name and append to list
    list_temp_subset = list(temp_subset)
    names(list_temp_subset) = as.character(j)
    temp_subsets = append(temp_subsets,list_temp_subset)
  }
  list_temp_subsets = list(temp_subsets)
  names(list_temp_subsets) = names(raw_subsets)[gender]
  ccs_cohorts = append(ccs_cohorts,list_temp_subsets)
}

# TRAINING MODELS ####
# Trainig complete case model and comparing it to suneelas
ccs_models = list()
ccs_comp_models = list() #ccs comparisons using suneelas equation

for (gender in 1:2){
  gender_ccs_models = list()
  gender_ccs_comp_models = list()
  
  for (j in 0:7){
    training_cohort = ccs_cohorts[[gender]][[j+1]][,-1]
    
    ccs_model = coxph(Surv(total_fu,indicator) ~ .,
                  data = training_cohort,
                  ties = "breslow",
                  model = TRUE)
    
    ccs_comp_model = coxph(Surv(total_fu,indicator) ~
                                      c_en_nzdep_q + 
                                      hx_vdr_diabetes +
                                      hx_af +
                                      ph_bp_lowering_prior_6mths +
                                      ph_lipid_lowering_prior_6mths + 
                                      ph_antithrombotic_prior_6mths +
                                      c_nhi_age + 
                                      eth_cat_Chinese + 
                                      eth_cat_Indian +
                                      eth_cat_NZM + 
                                      eth_cat_Other +
                                      eth_cat_Pacific +
                                      age_bp + 
                                      age_diabetes +
                                      age_af +
                                      bp_diabetes +
                                      at_diabetes + 
                                      bp_ll,
                                    data = training_cohort,
                                    ties = "breslow",
                                    model = TRUE)
    
    
    #list, name and append
    ccs_model  = list(ccs_model)
    names(ccs_model) = as.character(j)
    gender_ccs_models = append(gender_ccs_models,ccs_model)
    
    ccs_comp_model = list(ccs_comp_model)
    names(ccs_comp_model) = as.character(j)
    gender_ccs_comp_models = append(gender_ccs_comp_models,ccs_comp_model)
    
  }
  
  #list name and append
  gender_ccs_models = list(gender_ccs_models)
  names(gender_ccs_models) = names(raw_subsets)[gender]
  ccs_models = append(ccs_models,gender_ccs_models)
  
  gender_ccs_comp_models = list(gender_ccs_comp_models)
  names(gender_ccs_comp_models) = names(raw_subsets)[gender]
  ccs_comp_models = append(ccs_comp_models,gender_ccs_comp_models)
  
}


#summary statistics sanity check
stats_table = list()
#sanity check
for (i in 1:8){
  temp = data.frame(
    "Pattern" = (i-1),
    "CCS C" = summary(ccs_models[[1]][[i]])$concordance[[1]],
    "Comp C" = summary(ccs_comp_models[[1]][[i]])$concordance[[1]])
  stats_table = append(stats_table,list(temp))
}
stats_table = rbindlist(stats_table)
stats_table


# 5. Demographics tables ####
#simply run the demographics table function on the male and female subsets
source("0_helper_functions.R")
subset_list = list("Women" = list(subsets[[1]][,VSIMPLE_INDEX_MASTER]),
                   "Men" = list(subsets[[2]][,VSIMPLE_INDEX_MASTER]))
demo_table = demo_table_generator(subset_list,
                                  test_rows = TRUE)
print.data.frame(demo_table)
write.csv(demo_table,
          file = "tables/demo_table_an_gender.csv")


# patten subset predictions ##############
# 1. CCS 2. Comp CCS 
# use the previously calculated ps_cohorts
source("0_helper_functions.R")
pred_dt_list = list()
for (gender in 1:2){
  
  gender_ps_pred_dts = list()
  
  for (i in 1:8){
    #testing
    temp_cohort = ps_cohorts[[gender]][[i]]
    
    model_list = list("ccs" = ccs_models[[gender]][[i]])
  #                    "ccs_comp" = ccs_comp_models[[gender]][[i]])
    
    #create dataframe of specific pattern
    ps_pred_dt = pred_dt_generator(model_list,
                         temp_cohort)
    
    # datatable to list
    gender_ps_pred_dts = append(gender_ps_pred_dts,list(ps_pred_dt))
    
  }
  
  #rowbind hte list of data.frames
  gender_pred_dt = rbindlist(gender_ps_pred_dts)
  setkey(gender_pred_dt,VSIMPLE_INDEX_MASTER)
  
  #name and append this list of data frame
  list_gender_pred_dt = list(gender_pred_dt)
  names(list_gender_pred_dt) = names(raw_subsets)[gender]
  pred_dt_list = append(pred_dt_list,list_gender_pred_dt)

}


# mean/normal imputation ##################
imp_model_list = list()
for (i in 1:2){

  temp = copy(subsets[[i]])
  
  #set all missing values to 0 (healthy/normal value)
  temp[is.na(c_hba1c),c("c_hba1c"):=0]
  temp[is.na(c_tchdl),c("c_tchdl","age_tchdl"):=0]
  temp[is.na(c_tri),c("c_tri"):=0]
  temp[is.na(egfr_cat_G3a),c("egfr_cat_G3a","egfr_cat_G3b") := 0]
  #hist(temp$c_tchdl)
  
  imp_model = coxph(Surv(total_fu,indicator) ~ .,
                    data = temp[,-1],
                    ties = "breslow",
                    model = TRUE)
  
  #append model
  list_imp_model = list(imp_model)
  names(list_imp_model) = names(subsets)[i]
  imp_model_list = append(imp_model_list,list_imp_model) 
  
  
  #append predictions to pred_dt_list
  imp_pred = pred_dt_generator(list("imp_pred" = imp_model),
                               temp)
  imp_pred[,c("total_fu","indicator"):= NULL]
  
  pred_dt_list[[i]] = pred_dt_list[[i]][
                                imp_pred,
                                on = "VSIMPLE_INDEX_MASTER"]

}




#full comparison models
comp_model_list = list()
for (i in 1:2){
  
  #copy for safety
  temp = copy(subsets[[i]])
  
  comp_model = coxph(Surv(total_fu,indicator) ~
                           c_en_nzdep_q + 
                           hx_vdr_diabetes +
                           hx_af +
                           ph_bp_lowering_prior_6mths +
                           ph_lipid_lowering_prior_6mths + 
                           ph_antithrombotic_prior_6mths +
                           c_nhi_age + 
                           eth_cat_Chinese + 
                           eth_cat_Indian +
                           eth_cat_NZM + 
                           eth_cat_Other +
                           eth_cat_Pacific +
                           age_bp + 
                           age_diabetes +
                           age_af +
                           bp_diabetes +
                           at_diabetes + 
                           bp_ll,
                         data = temp,
                         ties = "breslow",
                         model = TRUE)
  
  #append model
  list_comp_model = list(comp_model)
  names(list_comp_model) = names(subsets)[i]
  comp_model_list = append(comp_model_list,list_comp_model) 
  
  
  #append predictions to pred_dt_list
  comp_pred = pred_dt_generator(list("comp" = comp_model),
                               temp)
  comp_pred[,c("total_fu","indicator"):= NULL]
  
  pred_dt_list[[i]] = pred_dt_list[[i]][
                                comp_pred,
                                on = "VSIMPLE_INDEX_MASTER"]
  
}



# export pred_dt
pred_dt_export = rbindlist(pred_dt_list,idcol = "gender")
pred_dt_export
write.fst(pred_dt_export,
          path = "geographic_analysis/pred_dt_export.fst")



# Coefficient tables ####
model_list = ccs_models
ccs_hr_table = coef_table_generator(model_list)
ccs_hr_table
write.csv(ccs_hr_table,
          "tables/ccs_hr_table.csv")


#normal imputation coef table
model_list = list(list("Women normal imp" = imp_model_list[[1]],
                  "Women CC cohort" = ccs_models$women$`7`,
                  "Women Comparison" = comp_model_list$Female 
                  ),
                  list("Men normal imp" = imp_model_list[[2]],
                       "Men CC cohort" = ccs_models$men$`7`,
                       "Men Comparison" = comp_model_list$Female))


                  
n_imp_hr_table = coef_table_generator(model_list)
n_imp_hr_table
output = n_imp_hr_table[,c(1,4,2,3,7,5,6)]
output
write.csv(output,
          "tables/n_imp_hr_table.csv")


#coefficient plot
female_models = list("Normal Imputation" = copy(imp_model_list[[1]]),
                     "Baseline Model" = copy(comp_model_list$Female))
male_models = list("Normal Imputation" = copy(imp_model_list[[2]]),
                   "Baseline Model" = copy(comp_model_list$Male))
all_models = list("women" = female_models,
                  "men" = male_models)

an_coef_plot = coef_plot_generator(all_models)
an_coef_plot
ggsave("an_coef_plot.png",
       an_coef_plot,
       path = "plots/misc",
       height = 7,
       width = 9)


#normal imputation coefficient plots
model_list = list("Women" = list("Normal imp" = imp_model_list[[1]],
                       "CCS" = ccs_models$women$`7`
),
"Men" = list("Normal imp" = imp_model_list[[2]],
     "Complete case cohort" = ccs_models$men$`7`))
n_imp_coef_plots = coef_plot_generator(model_list,
                                       return = 1)
n_imp_coef_plots

# CALIBRATION AND DISCRIMINATION #######################################
source("0_helper_functions.R")

#actual calibration code #####
i = "Male"
# # validation(pred_dt_list[[1]][VSIMPLE_INDEX_MASTER %in% ps_indexes[[1]][[8]]],
# #            type = "c")
# pred_dt_list_backup = copy(pred_dt_list[[1]])[,c("VSIMPLE_INDEX_MASTER","total_fu","indicator",
#                                          "ccs","comp")]

pred_dt_list_backup = copy(pred_dt_list)
pred_dt_list = copy(pred_dt_list_backup)

#changing pred_dt to focus on the models we want and renaming them
for (i in 1:2){
  # colours based on alphabetical order
  pred_dt_list[[i]][,"CCS" := ccs]
  pred_dt_list[[i]][,"Normal Imputation" := imp_pred]
  pred_dt_list[[i]][,"Baseline Model" := comp]
  pred_dt_list[[i]][,c("ccs","comp","imp_pred") := NULL]
}

#testing
source("0_helper_functions.R")
validation(pred_dt_list[[2]],
           type = "c")

# ggsave(paste("Women",i,".png",sep = ""),
#        plot = last_plot(),
#        height = 10,
#        width = 10,
#        path = "plots/validation_plots/an_calibration")

#PS calibration
for (gender in c("women","men")){
  plot_list = list()
  for (i in 1:8){
    
    subpop_pred_dt = pred_dt_list[[gender]][VSIMPLE_INDEX_MASTER %in% ps_indexes[[gender]][[i]]]
    temp = validation(subpop_pred_dt,
                      type = "c",
                      title = paste("PS index",(i-1),gender),
                      return = 1)
    plot_list = append(plot_list,list(temp))
  }
  
  arrange_plot = ggarrange(plotlist = plot_list,
                           ncol = 2,
                           nrow = ceiling(length(plot_list)/2),
                           common.legend = TRUE)
  
  ggsave(paste("pattern_",gender,".png",sep = ""),
         arrange_plot,
         path = "plots/validation_plots/an_calibration/subpopulations",
         height = ceiling(length(plot_list)/2)*5+2,
         width = 9)
  
  
}


# 9. Calibration

#remove unecessary models

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
         path = "plots/validation_plots/an_calibration",
         height = 5,
         width = 6.2)
}

whole_cal_plot = ggarrange(plotlist = whole_cal_list,
                           ncol = 1)
whole_cal_plot
ggsave("whole_calibration_plot.png",
       plot = whole_cal_plot,
       path = "plots/validation_plots/an_calibration",
       height = 10,
       width = 6)


# 9b. sub population plots ####
source("0_helper_functions.R")


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
    
    demo_table = demo_table_generator(index_list,test_rows = TRUE)
    write.csv(demo_table,
              file = paste("tables/an_subpop_demo/",gender,conds_name,".csv",sep = ""))
    
    ggsave(paste(conds_name,gender,".png",sep = ""),
           arrange_plot,
           path = "plots/validation_plots/an_calibration/subpopulations",
           height = ceiling(length(cal_plot_list)/2)*5,
           width = 9)
  }
  
}

#selected model plot
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
       path = "plots/validation_plots/an_calibration",
       height = 15,
       width = 9)



# 10. Discrimination #### 
source("0_helper_functions.R")
# 8a. whole cohort discrimination

#plot discrimiantion 
for (i in 1:2){
  #whole cohort calibration plots
  dis_plot = validation(pred_dt_list[[i]],
                        type = "d",
                        str_to_title(names(pred_dt_list)[[i]]),
                        return = 1)
  dis_plot
  
  ggsave(paste(names(pred_dt_list)[[i]],".png",sep=""),
         plot = dis_plot,
         path = "plots/validation_plots/an_discrimination",
         height = 5,
         width = 6.2)
}

# 10b. discrimiantion in sub population plots ####
source("0_helper_functions.R")

full_dis_plot_list = list()
for (i in 1:length(all_conds) ){
  
  conds_name = names(all_conds)[i]
  conds = all_conds[[i]]
  
  for (gender_ind in 1:length(pred_dt_list)){
    
    gender = names(pred_dt_list)[gender_ind]
    
    dis_plot_list = list()
    
    for (cond_ind in 1:length(conds)){
      
      #subset the relevant pred_dt
      sub_pop_pred_dt = pred_dt_list[[gender_ind]][
        VSIMPLE_INDEX_MASTER %in% subpop_index(conds[[cond_ind]]),]
      
      plot_name = eval(parse(text = names(conds)[cond_ind]))
      
      #plot using name of condition
      dis_plot = validation(sub_pop_pred_dt,
                            type = "d",
                            title = plot_name)
      
      #list, anme append
      list_dis_plot = list(dis_plot)
      names(list_dis_plot) = plot_name
      full_dis_plot_list = append(full_dis_plot_list,list_dis_plot)
      
      #append to calplot list ot feed into ggarange
      dis_plot_list = append(dis_plot_list,list_dis_plot)
      
    }
    
    arrange_plot = ggarrange(plotlist = dis_plot_list,
                             ncol = 2,
                             nrow = ceiling(length(dis_plot_list)/2))
    
    ggsave(paste(conds_name,gender,".png",sep = ""),
           arrange_plot,
           path = "plots/validation_plots/an_discrimination/subpopulations",
           height = ceiling(length(dis_plot_list)/2)*5,
           width = 9)
  }
  
}



#net reclassification code
dt = pred_dt_export
dt
sum(dt$indicator)

old_high_risk = dt[`comp` > 15]
new_high_risk = dt[`imp_pred` > 15]

dim(old_high_risk)
sum(old_high_risk$indicator)
mean(old_high_risk$indicator)

dim(new_high_risk)
sum(new_high_risk$indicator)
mean(new_high_risk$indicator)


high_new = new_high_risk[!(VSIMPLE_INDEX_MASTER %in% old_high_risk$VSIMPLE_INDEX_MASTER)]
high_old = old_high_risk[!(VSIMPLE_INDEX_MASTER %in% new_high_risk$VSIMPLE_INDEX_MASTER)]
high_both = old_high_risk[(VSIMPLE_INDEX_MASTER %in% new_high_risk$VSIMPLE_INDEX_MASTER)]





dim(high_new)
dim(high_old)
#dim(high_both)

sum(high_new$indicator)
sum(high_old$indicator)
#sum(high_both$indicator)

mean(high_new$indicator)
mean(high_old$indicator)
#mean(high_both$indicator)

#null cox model
high_hr_cohort = rbindlist(list("new" = high_new,"old" = high_old), idcol = TRUE)
high_hr_cohort[,model := factor(.id, levels = c("old","new")) ]
fit = coxph(Surv(total_fu,indicator) ~ model,
            data = high_hr_cohort,
            model = TRUE,
            ties = "breslow")
fit




old_low_risk = dt[`comp` < 5]
new_low_risk = dt[`imp_pred` < 5]

dim(old_low_risk)
sum(old_low_risk$indicator)
mean(old_low_risk$indicator)

dim(new_low_risk)
sum(new_low_risk$indicator)
mean(new_low_risk$indicator)

low_new = new_low_risk[!(VSIMPLE_INDEX_MASTER %in% old_low_risk$VSIMPLE_INDEX_MASTER)]
low_old = old_low_risk[!(VSIMPLE_INDEX_MASTER %in% new_low_risk$VSIMPLE_INDEX_MASTER)]

dim(low_new)
dim(low_old)

sum(low_new$indicator)
sum(low_old$indicator)

mean(low_new$indicator)
mean(low_old$indicator)

low_hr_cohort = rbindlist(list("new" = low_new,"old" = low_old), idcol = TRUE)
low_hr_cohort[,model := factor(.id, levels = c("old","new")) ]
fit = coxph(Surv(total_fu,indicator) ~ model,
            data = low_hr_cohort,
            model = TRUE,
            ties = "breslow")
fit





# old model comparions
# comp_model = coxph(Surv(total_fu,indicator) ~
#                      c_en_nzdep_q + 
#                      hx_vdr_diabetes +
#                      hx_af +
#                      ph_bp_lowering_prior_6mths +
#                      ph_lipid_lowering_prior_6mths + 
#                      ph_antithrombotic_prior_6mths +
#                      c_nhi_age + 
#                      eth_cat_Chinese + 
#                      eth_cat_Indian +
#                      eth_cat_NZM + 
#                      eth_cat_Other +
#                      eth_cat_Pacific +
#                      c_nhi_age*ph_bp_lowering_prior_6mths + 
#                      c_nhi_age*hx_vdr_diabetes +
#                      c_nhi_age*hx_af +
#                      ph_bp_lowering_prior_6mths*hx_vdr_diabetes +
#                      ph_antithrombotic_prior_6mths*hx_vdr_diabetes + 
#                      ph_bp_lowering_prior_6mths*ph_lipid_lowering_prior_6mths,
#                    data = training_cohort,
#                    ties = "breslow",
#                    model = TRUE)



# # misc plot
# pred_dt_export %>%
#   filter(indicator == 1) %>%
#   ggplot() +
#   geom_density(aes(x = imp_pred),colour = "#F8766D") +
#   geom_density(aes(x = comp),color = "#00BFC4") +
#   scale_x_log10()

# # 11. random bar plots ####
# by_eth_dt = an_cohort_ts %>%
#   group_by(eth_cat) %>%
#   summarize(event = mean(indicator)*100)
# 
# by_eth_plot = ggplot(by_eth_dt,aes(x = eth_cat,y = event)) +
#   geom_bar(stat = "identity") +
#   ylab("Event rate (%)") + 
#   xlab("Ethnicity")
# 
# ggsave("by_eth_bar.png",
#        by_eth_plot,
#        path = "plots/misc/",
#        height = 4,
#        width = 4)
# 
# 
# by_dep_dt = an_cohort_ts %>%
#   group_by(c_en_nzdep_q) %>%
#   summarize(event = mean(indicator)*100) %>%
#   mutate(dep = c_en_nzdep_q+3)
# 
# by_dep_plot = ggplot(by_dep_dt,aes(x = dep,y = event)) +
#   geom_bar(stat = "identity") +
#   ylab("Event rate (%)") + 
#   xlab("Deprivation quintile")
# 
# ggsave("by_dep_bar.png",
#        by_dep_plot,
#        path = "plots/misc/",
#        height = 4,
#        width = 4)




