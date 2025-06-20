# 0. Loading in necessary libraries and data ################################
library(Rcpp)
library(fst)
library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(survival)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr)
#library(VIM)
#library(mfp)
#library(survminer)
#library(glmnet)
library(Matrix)
library(stringr)
# library(grid)
# library(Gmisc)
# #library(rms)
# library(patchwork)
#library(missForest)
library(pROC)
library(boot)
library(caret)
library(survAUC)
library(survcomp)
library(Hmisc)
library(sys)

####################################################################################
# 1. Data table helper functions ####
###################################################################################

center_age = function(raw_subsets){
  # centers age inplace and store results
  # returns:
  # base_age: a named vector of the mean age of each subset
  base_age = c()
  for (raw_subset in raw_subsets){
    base_age = c(base_age,raw_subset[,mean(nhi_age)])
    raw_subset[,c_nhi_age := (nhi_age - mean(nhi_age))]
    raw_subset[,nhi_age := NULL]
  }
  names(base_age) = names(raw_subsets)
  return(base_age)
}


cut_factors = function(input_dt){
  # cuts up any categorical (factor) variables into
  # binary variables
  #
  # Inputs:
  # input_dt: datatable: datatable of predictors
  #
  # Outputs:
  # dt: datatable: returns a datatable with 
  #                all the categorical variables
  #               cut into bianry variables
  
  #copy data table for safety
  dt = copy(input_dt)
  
  #extract categorical variables
  vars = copy(names(dt))
  factor_bool = unlist(dt[,lapply(.SD,is.factor)])
  factor_vars = vars[factor_bool]
  for (var in factor_vars){
    
    #extract vector
    fact = copy(dt[,(var),with = FALSE])
    #unlist vector
    vect = unlist(fact)
    #extract levels in vector
    cat_levels = levels(vect)
    
    #just do it with for loops for now
    for (cat in cat_levels[-1]){
      dt[,(paste(var,cat,sep="_")) := as.integer(vect == cat)]
    }
    
    #slightly messy but gets the job done
    #reorder datatable so binary variables are in same place
    #as the categorical was
    new_vars = names(dt)[!(names(dt) %in% vars)]  
    index = which(names(dt) == var)
    new_order = append(vars,new_vars,after = index)
    setcolorder(dt,new_order)
    
    #remove factor variable from data table and var list
    dt[,(var) := NULL]
    vars = copy(names(dt))
  }
  
  #return new data.table
  return(dt)
}

prespec_ints = function(dt){
  #super simple justt adds suneelas 
  #interactions to a data table and returns the data table
  
  #copy for safety
  dt = copy(dt)
  
  #add old interactions
  dt[,":="(age_bp = c_nhi_age*ph_bp_lowering_prior_6mths,
           age_diabetes = c_nhi_age*hx_vdr_diabetes,
           age_af = c_nhi_age*hx_af,
           bp_diabetes = ph_bp_lowering_prior_6mths*hx_vdr_diabetes,
           at_diabetes = ph_antithrombotic_prior_6mths*hx_vdr_diabetes,
           bp_ll = ph_bp_lowering_prior_6mths*ph_lipid_lowering_prior_6mths)]
  
  #returns the new data frame
  return(dt)
  
}


#function to center the lab results 
transform_labs = function(dt,
                       c_egfr = 0,
                       c_hba1c = 0,
                       c_tchdl = 0,
                       c_tri = 0){

  #default does nothing
  #input what value you want to center
  
  #copy dt for safety
  dt = copy(dt)
  
  #centering results according to our prespecified criteria
  if (c_egfr == 1){
    tryCatch({
      dt[,":="(c_egfr = (egfr - 90)/15)]
      #removing uncentered test variables
      dt[,":="(egfr = NULL)]
    },
    error = function(x){
      #deliberately empty
    })
  }
  
  if (c_hba1c == 1){
    tryCatch({
      dt[,":="(c_hba1c = (hba1c - 40)/10)]
      #removing uncentered test variables
      dt[,":="(hba1c = NULL)]
    },
    error = function(x){
      #deliberately empty
    })
  }
  
  if (c_tchdl == 1){
    tryCatch({
      dt[,":="(c_tchdl = (tchdl - 3.5))]
      #removing uncentered test variables
      dt[,":="(tchdl = NULL)]
    },
    error = function(x){
      #deliberately empty
    })
  }
  
  if (c_tri == 1){
    tryCatch({
      dt[,":="(c_tri = (tri - 1.5))]
      #removing uncentered test variables
      dt[,":="(tri = NULL)]
    },
    error = function(x){
      #deliberately empty
    })
  }
  
  return(dt)
  

}

transform_labs_verbose = function(dt,
                       c_egfr = 0,
                       s_egfr = 1,
                       c_hba1c = 0,
                       s_hba1c = 1,
                       c_tchdl = 0,
                       s_tchdl = 1,
                       c_tri = 0,
                       s_tri = 1){

  #default does nothing
  #input what value you want to center
  
  #copy dt for safety
  dt = copy(dt)

  tryCatch({
      dt[,":="(c_egfr = (egfr - c_egfr)/s_egfr)]
      dt[,":="(c_hba1c = (hba1c - c_hba1c)/s_hba1c)]
      dt[,":="(c_tri = (tri - c_tri)/s_tri)]
      dt[,":="(c_tchdl = (tchdl - c_tchdl)/s_tchdl)]
      
      #removing uncentered test variables
      dt[,":="(egfr = NULL)]
      dt[,":="(hba1c = NULL)]
      dt[,":="(tri = NULL)]
      dt[,":="(tchdl = NULL)]
    },error = function(x){
      print("THERE WAS AN ERROR - PLEASE INSPECT")
    })
    
  return(dt)
}


egfr_groups = function(vect){
  # cuts egfr based on KDIGO guidelines
  # technically not perfect as only measures from 1 point in time
  # Input: numerical vector of results
  # output: string vector
  
  groups = cut(-vect,
               breaks = c(-Inf,-60,-45,0),
               labels = c("normal","G3a","G3b"))
  groups = relevel(groups, ref = "normal")
  
  return(groups)
}

###################################################################################
# 2. baseline demographics ####
###################################################################################

demo_table_generator = function(subset_list,
                                test_rows = FALSE){

  # input a list of lists with each elemetn int total list being vectors and returns a baseline demographics 
  # table with a column for each subset
  # name of element in list is columns name
  
  if (!and(exists("total_cohort"),exists("an_cohort_ts_temp")) ){
    total_cohort = read.fst("data/cohorts/hc_cohort_raw.fst",
                            as.data.table = TRUE)
    an_cohort_ts_temp = read.fst("data/cohorts/an_cohort_ts.fst",
                                 as.data.table = TRUE)
  }
  
  setkey(total_cohort, VSIMPLE_INDEX_MASTER)
  
  #some useful numbers to have handy
  total_n = dim(total_cohort[VSIMPLE_INDEX_MASTER %in% unlist((subset_list))])[1]
  
  temp_tibble_table_list = list()
  for (i in 1:length(subset_list)){
    #extract cohort
    indexes = unlist(subset_list[[i]])
    dt = total_cohort[VSIMPLE_INDEX_MASTER %in% indexes,]
    
    #number of participants
    n = matrix(c("Participants",sprintf("%s (%.1f%%)",format(dim(dt)[1],big.mark = ","),dim(dt)[1]*100/total_n)), byrow = TRUE, ncol = 2)
    
    #Age
    age = matrix(c("Age",sprintf("%.1f (%.1f)",mean(dt$nhi_age),sd(dt$nhi_age))),byrow = TRUE, ncol = 2)
    
    #ethnicity
    eth = matrix(c("New Zealand European",sprintf("%s (%.1f%%)",format(sum(dt$eth_cat == "NZE"), big.mark = ","),sum(dt$eth_cat == "NZE")*100/dim(dt)[1]),
                 "Maaori",sprintf("%s (%.1f%%)",format(sum(dt$eth_cat == "NZM"), big.mark = ","),sum(dt$eth_cat == "NZM")*100/dim(dt)[1]),
                 "Pacific",sprintf("%s (%.1f%%)",format(sum(dt$eth_cat == "Pacific"), big.mark = ","),sum(dt$eth_cat == "Pacific")*100/dim(dt)[1]),
                 "Indian",sprintf("%s (%.1f%%)",format(sum(dt$eth_cat == "Indian"), big.mark = ","),sum(dt$eth_cat == "Indian")*100/dim(dt)[1]),
                 "Chinese",sprintf("%s (%.1f%%)",format(sum(dt$eth_cat == "Chinese"), big.mark = ","),sum(dt$eth_cat == "Chinese")*100/dim(dt)[1]),
                 "Other",sprintf("%s (%.1f%%)",format(sum(dt$eth_cat == "Other"), big.mark = ","),sum(dt$eth_cat == "Other")*100/dim(dt)[1])), byrow = TRUE, ncol = 2)
    
    #deprivation
    dep = matrix(c("1",sprintf("%s (%.1f%%)",format(sum(dt$c_en_nzdep_q == "-2"),big.mark = ","),sum(dt$c_en_nzdep_q == "-2")*100/dim(dt)[1]),
                   "2",sprintf("%s (%.1f%%)",format(sum(dt$c_en_nzdep_q == "-1"),big.mark = ","),sum(dt$c_en_nzdep_q == "-1")*100/dim(dt)[1]),
                   "3",sprintf("%s (%.1f%%)",format(sum(dt$c_en_nzdep_q == "0"),big.mark = ","),sum(dt$c_en_nzdep_q == "0")*100/dim(dt)[1]),
                   "4",sprintf("%s (%.1f%%)",format(sum(dt$c_en_nzdep_q == "1"),big.mark = ","),sum(dt$c_en_nzdep_q == "1")*100/dim(dt)[1]),
                   "5",sprintf("%s (%.1f%%)",format(sum(dt$c_en_nzdep_q == "2"),big.mark = ","),sum(dt$c_en_nzdep_q == "2")*100/dim(dt)[1])), byrow = TRUE, ncol = 2)
    
                   
    #conditions
    disease = matrix(c("Diabetes",sprintf("%s (%.1f%%)",format(sum(dt$hx_vdr_diabetes),big.mark = ","),sum(dt$hx_vdr_diabetes)*100/dim(dt)[1]),
                       "Atrial fibrillation",sprintf("%s (%.1f%%)",format(sum(dt$hx_af),big.mark = ","),sum(dt$hx_af)*100/dim(dt)[1])), byrow = TRUE, ncol = 2)
    
    #baseline medication
    bd = matrix(c("Blood-pressure-lowering",
                  sprintf("%s (%.1f%%)",format(sum(dt$ph_bp_lowering_prior_6mths),big.mark = ","),sum(dt$ph_bp_lowering_prior_6mths)*100/dim(dt)[1]),
                  "Lipid-lowering",
                  sprintf("%s (%.1f%%)",format(sum(dt$ph_lipid_lowering_prior_6mths),big.mark = ","),sum(dt$ph_lipid_lowering_prior_6mths)*100/dim(dt)[1]),
                  "Anti-thrombotic",
                  sprintf("%s (%.1f%%)",format(sum(dt$ph_antithrombotic_prior_6mths),big.mark = ","),sum(dt$ph_antithrombotic_prior_6mths)*100/dim(dt)[1])), 
                byrow = TRUE, ncol = 2)
    
    # follow-up
    fu = matrix(c("Total follow-up",sprintf("%s (%.1f)",format(sum(dt$total_fu),big.mark = ","),mean(dt$total_fu))), byrow = TRUE, ncol = 2)
    
    if(test_rows){
      
      dt_tests = an_cohort_ts_temp[VSIMPLE_INDEX_MASTER %in% indexes]
      
      if (dim(dt_tests)[1] != dim(dt)[1]){
        stop("COOHORT ISSUES")
      }
      
      #test completeness
      test_comp = matrix(c("n eGFR",sprintf("%s (%.1f%%)",format(sum(!is.na(dt_tests$egfr)),big.mark = ","),sum(!is.na(dt_tests$egfr))*100/dim(dt_tests)[1]),
                           "n HbA1c",sprintf("%s (%.1f%%)",format(sum(!is.na(dt_tests$hba1c)),big.mark = ","),sum(!is.na(dt_tests$hba1c))*100/dim(dt_tests)[1]),
                           "n TC/HDL ratio",sprintf("%s (%.1f%%)",format(sum(!is.na(dt_tests$tchdl)),big.mark = ","),sum(!is.na(dt_tests$tchdl))*100/dim(dt_tests)[1]),
                           "n TG",sprintf("%s (%.1f%%)",format(sum(!is.na(dt_tests$tri)),big.mark = ","),sum(!is.na(dt_tests$tri))*100/dim(dt_tests)[1])), byrow = TRUE, ncol = 2)
      #test_comp
      
      #test results
      test_res = matrix(c("eGFR",sprintf("%.1f (%.1f)",mean(dt_tests$egfr,na.rm = TRUE), sd(dt_tests$egfr,na.rm = TRUE)),
                          "HbA1c",sprintf("%.1f (%.1f)",mean(dt_tests$hba1c,na.rm = TRUE), sd(dt_tests$hba1c,na.rm = TRUE)),
                          "TC/HDL ratio",sprintf("%.1f (%.1f)",mean(dt_tests$tchdl,na.rm = TRUE), sd(dt_tests$tchdl,na.rm = TRUE)),
                          "TG",sprintf("%.1f (%.1f)",mean(dt_tests$tri,na.rm = TRUE), sd(dt_tests$tri,na.rm = TRUE))) , byrow = TRUE, ncol = 2)
      
      #test_res
      
      # completeness
      #outcomes
      outcome = matrix(c("Events",sprintf("%s (%.1f%%)",format(sum(dt$indicator),big.mark = ","),sum(dt$indicator)*100/dim(dt)[1])), byrow = TRUE, ncol = 2)
      
      #combine into matrix   
      temp_tibble_table = tibble(data.frame(rbind(n,age,eth,dep,disease,bd,fu,test_comp,test_res,outcome)))
      names(temp_tibble_table) = c("Characteristics",names(subset_list)[i])
      temp_tibble_table_list = append(temp_tibble_table_list,list(temp_tibble_table))
      
    } else {
      #outcomes
      outcome = matrix(c("Events",sprintf("%s (%.1f%%)",format(sum(dt$indicator),big.mark = ","),sum(dt$indicator)*100/dim(dt)[1])), byrow = TRUE, ncol = 2)
      
      #combine into matrix   
      temp_tibble_table = tibble(data.frame(rbind(n,age,eth,dep,disease,bd,fu,outcome)))
      names(temp_tibble_table) = c("Characteristics",names(subset_list)[i])
      temp_tibble_table_list = append(temp_tibble_table_list,list(temp_tibble_table))
    }
    
    
    
    
  }
  
  temp_tibble_table_list
  table = Reduce(function(x,y) left_join(x,y,by = "Characteristics"),temp_tibble_table_list)
  
}

###################################################################################
# 3. Visualize model coefficients ####
###################################################################################

#predictor names/order vectors ####
predictor_type =  data.frame(rn = c("age_bp",
                                     "age_diabetes",
                                     "age_af",
                                     "bp_diabetes",                  
                                     "at_diabetes",
                                     "bp_ll",
                                     "ph_antithrombotic_prior_6mths",
                                     "ph_lipid_lowering_prior_6mths",
                                     "ph_bp_lowering_prior_6mths",
                                     "hx_af",
                                     "hx_vdr_diabetes",
                                     "eth_cat_Other",  
                                     "eth_cat_Indian",
                                     "eth_cat_Chinese",
                                     "eth_cat_Pacific",
                                     "eth_cat_NZM",
                                     "c_en_nzdep_q",
                                     "c_nhi_age",
                                     "age_tchdl",
                                     "c_tri",
                                     "c_tchdl", 
                                     "c_hba1c",
                                     "egfr_cat_G3a",
                                     "egfr_cat_G3b"),
                    type = factor(c("Interactions",
                             "Interactions",
                             "Interactions",
                             "Interactions",
                             "Interactions",
                             "Interactions",
                             "Baseline dispensing",
                             "Baseline dispensing",
                             "Baseline dispensing",
                             "Diagnoses",
                             "Diagnoses",
                             "Demographics",
                             "Demographics",
                             "Demographics",
                             "Demographics",
                             "Demographics",
                             "Demographics",
                             "Demographics",
                             "New interactions",
                             "Biochemical",
                             "Biochemical",
                             "Biochemical",
                             "Biochemical",
                             "Biochemical"),
                             levels = c("Biochemical",
                                       "New interactions",
                                       "Demographics",
                                       "Diagnoses",
                                       "Baseline dispensing",
                                       "Interactions")))

#cleaned names 
clean_predictors = data.frame( rn = c("eth_cat_Chinese",
                                     "eth_cat_Indian",
                                     "eth_cat_NZM",
                                     "eth_cat_Other",
                                     "eth_cat_Pacific",
                                     "age_bp",
                                     "age_diabetes",
                                     "age_af",
                                     "bp_diabetes",
                                     "at_diabetes",
                                     "bp_ll",
                                     "ph_bp_lowering_prior_6mths",
                                     "ph_lipid_lowering_prior_6mths",
                                     "ph_antithrombotic_prior_6mths",
                                     "hx_vdr_diabetes",
                                     "hx_af",
                                     "c_en_nzdep_q",
                                     "c_nhi_age",
                                     "c_hba1c",
                                     "c_tchdl",
                                     "age_tchdl",
                                     "c_tri",
                                     "egfr_cat_G3a",
                                     "egfr_cat_G3b"),
                            predictors = factor(c("Chinese",
                                             "Indian",
                                             "Maaori",
                                             "Other",                
                                             "Pacific",
                                             "Age*BP lowering",
                                             "Age*diabetes" ,
                                             "Age*af",
                                             "BP lowering*diabetes",                  
                                             "Antithrombotic*diabetes",
                                             "BP lowering*lipid lowering",
                                             "BP lowering medication",
                                             "Lipid lowering medication",
                                             "Antithrombotic medication",
                                             "Diabetes",
                                             "Atrial fibrillation",
                                             "Deprivation",
                                             "Age",  
                                             "HbA1c",
                                             "TC/HDLC ratio",                      
                                             "Age*TC/HDL ratio",
                                             "Triglycerides",
                                             "G3a CKD",
                                             "G3b CKD"),
                                             levels = c("Age*BP lowering",
                                                        "Age*diabetes" ,
                                                        "Age*af",
                                                        "BP lowering*diabetes",                  
                                                        "Antithrombotic*diabetes",
                                                        "BP lowering*lipid lowering",
                                                        "Antithrombotic medication",
                                                        "Lipid lowering medication",
                                                        "BP lowering medication",
                                                        "Atrial fibrillation",
                                                        "Diabetes",
                                                        "Other", 
                                                        "Chinese",
                                                        "Indian",
                                                        "Pacific",
                                                        "Maaori",
                                                        "Deprivation",
                                                        "Age",  
                                                        "Age*TC/HDL ratio",
                                                        "Triglycerides",
                                                        "TC/HDLC ratio",   
                                                        "HbA1c",
                                                        "G3a CKD",
                                                        "G3b CKD")))
                            
# creates a coefficient dt from a list of a list of cohorts
coef_dt_generator = function(model_list_list){
  
  #copy for safety
  model_list_list = copy(model_list_list)
  
  dt_list = list()
  for (i in 1:length(model_list_list)){
    list_name = str_to_title(names(model_list_list[i]))
    model_list = model_list_list[[i]]
    
    for (i in 1:length(model_list)){
      
      name = str_to_title(names(model_list)[i])
      
      dt = data.table(summary(model_list[[i]])$coef, keep.rownames = TRUE)
      dt[,":="(UL = exp(coef + `se(coef)`),
               LL = exp(coef - `se(coef)`),
               p = `Pr(>|z|)`,
               model = name,
               list_name = list_name)]
      dt_list = append(dt_list,list(dt))
    }
  }
  
  #clean up the plot dt for plotting
  coef_dt = rbindlist(dt_list)
  coef_dt = coef_dt[clean_predictors,on = "rn",nomatch = 0]
  coef_dt = coef_dt[predictor_type,on = "rn", nomatch = 0]
  return(coef_dt)
  
}

#plots the coef dt 
coef_plot_generator = function(model_list_list,
                     return = 0){
  #plots the coefficients of a model
  #Input:
  # model_lists :list of lists : higher degree? list is different subset (e.g. men and women)
  #                second layer of lists = list of model objects
  # return: integer: return model or not
  
  #turn model_list_list into coef dt
  plot_dt = coef_dt_generator(model_list_list)
  
  plot_dt = plot_dt %>% 
    arrange(rn,type)
  
  plot_dt$predictors
  #note coorinate flip to work with geom_point range
  coef_plot = ggplot(plot_dt[order(predictors)]) +
        geom_pointrange(aes(x = predictors,
                        y = `exp(coef)`, 
                        ymin = LL, 
                        ymax = UL,
                        color = model,
                        shape = model),
                    position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 1, color = "black", linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = c("#00BFC4","#F8766D")) +
    labs(color = "Models") +
    guides(shape = "none") +
    xlab("Predictors") +
    facet_grid(type ~ list_name, scales = "free",space = "free") +
    ylab("Hazard ratios") +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.text.y = element_text(angle = 0))

  
  coef_plot
  
  if (return == 1){
    return(coef_plot)
  } else {
    coef_plot
  }
}


# turns the coef dt into a table
coef_table_generator = function(model_list_list){
  # turn a list of model lists into a nice table
  # Inputs:
  # model_list_list - list of models lists
  
  #generate coef_dt
  coef_dt = coef_dt_generator(model_list_list)
  
  #clean up dt
  dt = coef_dt[,c("rn","predictors","list_name","model","exp(coef)","UL","LL","p")]
  dt[,p := round(p,3)]
  dt[,`hazard ratio` := sprintf("%.2f (%.2f, %.2f)",`exp(coef)`,LL,UL)]
  
  #dcast to essentially cast to predictor
  dt_clean = dcast(dt, rn + predictors ~ list_name + model, value.var = c("hazard ratio"))
  
  # #reshape dt_clean specifying 2 fixed columns - only needed if casting multiple columns
  # fixed = 2
  # index = c(1:fixed, rbind( (fixed+1):(fixed+(dim(dt_clean)[2]-2)/2),(fixed+(dim(dt_clean)[2]-2)/2+1):(dim(dt_clean)[2])))
  # dt_clean = dt_clean[,c(index),with = FALSE]

  #quite manually impose order based on clean_predictors vector
  #bind to dt_clean and order
  final_dt = dt_clean[order(-predictors)][,c("rn"):= NULL]
  
  return(final_dt)
}
  
###################################################################################
# 4. prediction functions ####
###############################################################################

#cv_coxmodel fit and predict
cv_coxph_fit_and_predict = function(dt,
                                  k_folds = 5,
                                  grouped_cv = FALSE,
                                  return_vector = FALSE,
                                  return_lp = FALSE){
  # function to fit a coxph model and predict on a cross validated data set
  # inputs:
  # dt = data table to fit model on 
  # k_folds = number of folds
  # grouped_cv = boolean to indicate if cv is grouped
  # return_vector = boolean to indicate if a vector of predictions should be returned
  # Ouputs:
  # a pred_dt datatable

  #set.seed(42) # for reproducibility
  dt = copy(dt)

  if(grouped_cv == TRUE){
    unique_ids = unique(dt$VSIMPLE_INDEX_MASTER)
    id_fold_mapping = data.table(
      VSIMPLE_INDEX_MASTER = unique_ids,
      fold = sample(1:k_folds, size = length(unique_ids), replace = TRUE)
    )
    # Map folds back to original data
    dt = merge(dt, id_fold_mapping, by = "VSIMPLE_INDEX_MASTER")
    fold_indices = dt$fold
  } else {
    fold_indices = sample(1:k_folds, size = nrow(dt), replace = TRUE)
  }

  #pred_dt to store results
  pred_dt = copy(dt)[,c("VSIMPLE_INDEX_MASTER","indicator","total_fu"),with=FALSE]
  pred_dt[,pred_risk := as.double(0)]
  # For each fold
  for (fold in 1:k_folds) {
    # Split into training and validation
    train_data = dt[fold_indices != fold[1],]
    valid_data = dt[fold_indices == fold[1],]
    
    # Remove index from training data
    train_data = train_data[,-c("VSIMPLE_INDEX_MASTER")]
    
    # Fit model on training data
    model = coxph(Surv(total_fu, indicator) ~ .,
                  data = train_data,
                  ties = "breslow",
                  model = TRUE)

    # Predict on validation set
    pred_risk_vector = prediction_helper(model,valid_data,return_lp = return_lp)
    
    # can calculate other metrics here if needed

    # Store predictions
    pred_dt[fold_indices == fold, pred_risk := pred_risk_vector]
  }
  if(return_vector == TRUE){
    return(pred_dt$pred_risk)
  } else {
    return(pred_dt)
  }
}

base_surv = function(model){
  #returns the baseline 5 year survival of the model
  
  bs = exp(-max(basehaz(model,centered = FALSE)$hazard))
  
  return(bs)
}

prediction_helper = function(model,
                             dt,
                             return_lp = FALSE){
  #simple helper function to create prediciton of 5-year absolut CVD risk
  
  #copy data.table for safety
  dt = copy(dt)
  
  #calculating baseline survival
  bs = base_surv(model)
  
  #calculating linear predictor
  lp = predict(model,
               newdata = dt,
               type = "lp",
               reference = "zero")
  
  if(return_lp == TRUE){
    return(lp)
  } else {
    #converting linear predictor into 5 year event estimate
    pred = (1 - (bs^exp(lp)) ) * 100
    #return pred
    return(pred)
  }
}

prediction_wrapper = function(model_data,
                             new_data,
                             name,
                             return_lp = FALSE){
  # wrapper for prediction_helper to enable pattern submodel prediction
  
  #Inputs:
  # models: either a coxph object or a list of coxph objects
  # new_data: either a dt or list of dt's
  # !assume both are single object or both ar list
  
  if (class(model_data) == "coxph"){
    #run simple prediction
    
    pred = prediction_helper(model_data,new_data,return_lp = return_lp)
    output = new_data[,c("VSIMPLE_INDEX_MASTER")][,(name) := pred]
    
    
    
  } else if (class(model_data) == "list"){
    
    outputs = list()
    for (i in 1:length(model_data)){
      pred = prediction_helper(model_data[[i]],new_data[[i]],return_lp = return_lp)
      outputs = append(outputs,new_data[[i]][
        ,c("VSIMPLE_INDEX_MASTER")][
          ,(name):= pred])
    }
    #bind together outputs
    output = rbindlist(outputs)
    
    
  }
  
  #set key of output
  setkey(output,VSIMPLE_INDEX_MASTER)
  return(output)
}

pred_dt_generator = function(model_list,
                          new_data,
                          return_lp = FALSE){
  #creates a prediction dt for a set of models and some input data
  #Inputs: 
  # model_list: list of coxph models (ideally named)
  # new_data: list of dt or 1 dt 
  #   subset(s) of a cohort to predict on 
  #   If only one provided models applied to it
  #   if multiple provided length must match model_list and models applied
  #   to corresponding element and predictions joined
  # 
  
  
  # copy new_data for safety
  dt_list = list()
  if (class(new_data)[1] == "list"){
    #multiple new_data inputted
    #simply copy
    for (i in 1:length(model_list)){
      dt_list = append(dt_list,list(copy(new_data[[i]])))
    }
  } else if (class(new_data)[1] == "data.table") {
    #if only 1 inputted then simply duplicate to length of models
    for (i in 1:length(model_list)){
      dt_list = append(dt_list,list(copy(new_data)))
    }
    
  }
  
  #extract index,total_fu and indicator from first dt
  pred_dt = dt_list[[1]][,c("VSIMPLE_INDEX_MASTER",
                            "total_fu",
                            "indicator")]
  
  #calculate predictions and append to pred_dt
  for (i in 1:length(model_list)){
    
    #extract name
    model_name = names(model_list)[i]
    
    # actual prediction
    #tolerates passing in list of models and dt's
    predictions = prediction_wrapper(model_list[[i]],
                             dt_list[[i]],
                             model_name,
                             return_lp = return_lp)
    
    #joining prediciton to pred_dt
    pred_dt = pred_dt[predictions,on="VSIMPLE_INDEX_MASTER"]
    
  }
  
  #return pred_dt
  return(pred_dt)
}

#####################################################################################
# 5. Boostrapping
####################################################################################

# 11a. i) setting up the bootstrap
pred_dt_bootstrap = function(pred_dt,statistic){
  # Calculates the boostrap CI of a statistic on the pred_dt
  #
  # Inputs:
  # pred_dt = a pred_dt datatable
  # statistic = a bootstrap compatible function
  #
  # Outputs:
  # a sting 

  results_new = boot(pred_dt, statistic = statistic, R = 100)
  results_baseline = boot(pred_dt,statistic = statistic, R = 100, model= "Baseline Model")
  diffs = results_new$t - results_baseline$t
  #output formating
  ci_diffs = quantile(diffs, probs = c(0.025, 0.975))
  #output = sprintf("%.4f (%.4f,%.4f)",mean(diffs),ci_diffs[1],ci_diffs[2])
  output = sprintf("%.3f (%.3f, %.3f), %.3f (%.3f, %.3f),\n diff: %.3f (%.3f, %.3f) ",
  results_new$t0,
  quantile(results_new$t,probs = c(0.025,.975))[1],
  quantile(results_new$t,probs = c(0.025,.975))[2],
  results_baseline$t0,
  quantile(results_baseline$t,probs = c(0.025,.975))[1],
  quantile(results_baseline$t,probs = c(0.025,.975))[2],
  mean(diffs),
  quantile(diffs,probs = c(0.025,.975))[1],
  quantile(diffs,probs = c(0.025,.975))[2])

  print(output)
  return(output)
}

#boostrap functions
auc_bootstrap_func = function(pred_dt,ind,model="New Model"){
  d = pred_dt[ind,] #select the boostrap sample
  #print(d)
  return( auc(d[,"indicator"][[1]],d[,.SD,.SDcols = model][[1]],quiet=TRUE)[1] )
}

brier_bootstrap_func = function(pred_dt,ind,model="New Model"){
  d = pred_dt[ind,] #select the boostrap sample
  prob = d[,.SD,.SDcols = model][[1]] / 100
  obs = d[,"indicator"][[1]]
  #manual
  output = brier_func(prob,obs)
  return(output)
}

brier_func = function(prob,obs){
  #simpe function to calculate MSE
  MSE = mean( (prob - obs)^2 )
  return(MSE)
}

R2_bootstrap_func = function(pred_dt,ind,model="New Model"){
  d = pred_dt[ind,] #select the boostrap sample
  prob = d[,.SD,.SDcols = model][[1]] / 100
  obs = d[,"indicator"][[1]]

  #manual
  SST = mean((mean(obs)-obs)^2)
  SSR = mean((prob - obs)^2)
  R2 = (SST - SSR)/SST

  return(R2)
}

pred_dt_bootstrap_wrapper = function(pred_dt){
  #Computes the auc, brier score and R2 and displays the CI and differences in text
  x = paste0("AUC: ",pred_dt_bootstrap(pred_dt,auc_bootstrap_func),"\n",
  "brier: ",pred_dt_bootstrap(pred_dt,brier_bootstrap_func),"\n",
  "R2: ",pred_dt_bootstrap(pred_dt,R2_bootstrap_func))
  return(x)
}


# 11a model level bootstrappning ###################
# can't be done on subpops onlly on the whole moel
AIC_bootstrap_func =  function(dt_data,ind){
  dt = copy(dt_data[ind,]) #select the boostrap sample
  dt = dt[,-c("VSIMPLE_INDEX_MASTER")]

  #fit model
  model = coxph(Surv(total_fu,indicator) ~ .,
                data = dt,
                ties = "breslow",
                model = TRUE)
  
  AIC_statistic = AIC(model)
  return(AIC_statistic)
}

dt_bootstrap = function(dt1,dt2,statistic,R=100){
  # Calculates the boostrap CI of a statistic on the dt1 and dt2
  #
  # Inputs:
  # dt1 and 2 = a cohort datatable
  # statistic = a bootstrap compatible function
  #
  # Outputs:
  # a sting 

  results_new = boot(dt1, statistic = statistic, R = R)
  results_baseline = boot(dt2,statistic = statistic, R = R)
  diffs = results_new$t - results_baseline$t
  #output formating
  ci_diffs = quantile(diffs, probs = c(0.025, 0.975))
  #output = sprintf("%.4f (%.4f,%.4f)",mean(diffs),ci_diffs[1],ci_diffs[2])
  output = sprintf("New: %.4f (%.4f, %.4f) Baseline: %.4f (%.4f, %.4f), diff: %.4f (%.4f, %.4f) ",
  results_new$t0,
  quantile(results_new$t,probs = c(0.025,.975))[1],
  quantile(results_new$t,probs = c(0.025,.975))[2],
  results_baseline$t0,
  quantile(results_baseline$t,probs = c(0.025,.975))[1],
  quantile(results_baseline$t,probs = c(0.025,.975))[2],
  mean(diffs),
  quantile(diffs,probs = c(0.025,.975))[1],
  quantile(diffs,probs = c(0.025,.975))[2])

  print(output)
  return(output)
}

model_boot_wrapper = function(dt1,dt2){
  x = paste0("AIC: ",dt_bootstrap(dt1,dt2,AIC_bootstrap_func,R=10))
  return(x)
}




# royston_boot_func(data,ind){
#   #remove index
#   dt = data[ind,-c("VSIMPLE_INDEX_MASTER")]
  
#   #fit model
#   model = coxph(Surv(total_fu,indicator) ~ .,
#                 data = dt,
#                 ties = "breslow",
#                 model = TRUE)

#   #calculate R^2 statistic
# }


###################################################################################
###################################################################################
# PLOTTING #####
###################################################################################
###################################################################################



# Calibration plot ####
cal_plot_gg = function(dt_list,
                       plot_title,
                       model_names,
                       plot_subtitle = NULL){
  
  
  #rowbind data
  names(dt_list) = model_names
  plot_data = rbindlist(dt_list, use.names = FALSE, idcol = "model")
  
  #print(max(plot_data[,2:3]))
  #print(plot_data)
  
  if(length(dt_list) == 2){
    
    plot_data$model = factor(plot_data$model, levels = c("New Model","Baseline Model"))
    
    plot = ggplot(plot_data[order(-model)]) +
      geom_abline(slope = 1, 
                  intercept = 0, 
                  linetype = "dashed",
                  alpha = 0.3) +
      geom_point(aes(x = obs, y = pred,shape = model, color = model),size = 1.5) +
      xlim(0,5*(max(plot_data[,2:3]) %/% 5 + 1)) +
      ylim(0,5*(max(plot_data[,2:3]) %/% 5 + 1)) +
      coord_fixed() +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = "Observed events (%)",
           y = "Mean predicted 5 year risk (%)") +
      scale_shape_manual(values = c(17,16)) +
      scale_color_manual(values = c("#F8766D","#00BFC4")) +
      theme_bw() +
      theme(plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 9),
            legend.position = c(0.23,0.87),   #legend formatting options
            legend.title = element_blank(),
            legend.background = element_rect(size = 0,
                                             linetype = "solid",
                                             colour = "white"),
            plot.margin = unit(c(0.02,0,0.02,0),"npc")) 
    
    plot
    
  } else if(length(dt_list) == 3){
    
    print("TEST")
    
    plot_data$model = factor(plot_data$model, levels = c("Normal Imputation","CCS","Baseline Model"))
    
    plot = plot_data[order(-model)] %>%
      ggplot() +
      geom_abline(slope = 1, 
                  intercept = 0, 
                  linetype = "dashed",
                  alpha = 0.3) +
      geom_point(aes(x = obs, y = pred, shape = model,color = model),size = 2,stroke = 1.2) +
      xlim(0,5*(max(plot_data[,2:3]) %/% 5 + 1)) +
      ylim(0,5*(max(plot_data[,2:3]) %/% 5 + 1)) +
      coord_fixed() +
      labs(title = paste("Calibration:",plot_title),
           subtitle = plot_subtitle,
           x = "Observed events (%)",
           y = "Mean predicted 5 year risk (%)") +
      scale_shape_manual(values = c(3,4,20)) +
      scale_color_manual(values = c("#F8766D","#00BA38","dark grey")) +
      theme_bw() +
      theme(plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 9),
            legend.position = c(0.23,0.85),   #legend formatting options
            legend.title = element_blank(),
            legend.background = element_rect(size = 0,
                                             linetype = "solid",
                                             colour = "white"),
            plot.margin = unit(c(0.02,0,0.02,0),"npc"))
    

  }
  
  
  
  
  return(plot)
  
}

# Discrimination plot ####
dis_plot_gg = function(dt_list,
                       plot_title,
                       model_names,
                       plot_subtitle = NULL){
  
  #rowbind data
  names(dt_list) = model_names
  
  #change predicted risk to rank and obsevered rate to proprotion of observed
  plot_data = rbindlist(dt_list, use.names = FALSE, idcol = "model")
  
  #print(max(plot_data[,2:3]))
  #print(plot_data)
  #print("HELLLOO")
  
  
  
  if(length(dt_list) == 2){
    
    plot_data$model = factor(plot_data$model, levels = c("New Model","Baseline Model"))
    
    plot = ggplot(plot_data[order(-model)]) +
      geom_point(aes(x = decile, y = per, color = model,shape = model), size = 1.5) +
      ylim(0,5*(max(plot_data[,5]) %/% 5 + 1)) +
      labs(title = paste("Discrimination:",plot_title),
           subtitle = plot_subtitle) +
      scale_x_continuous(name = "Deciles of predicted risk",
                         limits = c(0,10),
                         breaks = 1:10) +
      ylab("Observed events (%)") +
      scale_shape_manual(values = c(17,16)) +
      scale_color_manual(values = c("#F8766D","#00BFC4")) +
      theme_bw() +
      theme(plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 9),
            legend.position = c(0.23,0.84),   #legend formatting options
            legend.title = element_blank(),
            legend.background = element_rect(size = 0,
                                             linetype = "solid",
                                             colour = "white"),
            plot.margin = unit(c(0.08,0.05,0.08,0.05),"npc")) 
    
  } else if(length(dt_list) == 3){
    
    print("TESTING")
    
    plot_data$model = factor(plot_data$model, levels = c("Normal Imputation","CCS","Baseline Model"))
    print.data.frame(plot_data)
    
    
    
    plot = ggplot(plot_data) +
      geom_point(aes(x = decile, y = per, color = model,shape = model), size = 2, stroke = 1.2) +
      ylim(0,5*(max(plot_data[,5]) %/% 5 + 1)) +
      labs(title = paste("Discrimination:",plot_title),
           subtitle = plot_subtitle) +
      scale_x_continuous(name = "Deciles of predicted risk",
                         limits = c(0,10),
                         breaks = 1:10) +
      ylab("Observed events (%)") +
      scale_shape_manual(values = c(3,4,20)) +
      scale_color_manual(values = c("#F8766D","#00BA38","dark grey")) +
      theme_bw() +
      theme(plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 9),
            legend.position = c(0.18,0.88),   #legend formatting options
            legend.title = element_blank(),
            legend.background = element_rect(size = 0,
                                             linetype = "solid",
                                             colour = "white"),
            plot.margin = unit(c(0.08,0.2,0.08,0.2),"npc")) 
  }

  return(plot)
  ###FIx hte LEGEND!!!
  
}

# wrapper functions ####
validation_plot = function(risk_vectors,
                           indicator_vector,
                           total_fu_vector,
                           type,
                           title = NULL,
                           subtitle = NULL){
  
  # wrapper functions function which takes a vector of risks and a vector of
  # outcomes and creates a risk group df and feeds it into either
  # cal_plot_gg or dis_plot_gg, returning the result
  #
  # Inputs:
  # risk_vector - list of vectors of risk estimates
  # indicator_vector - vectors of events 
  # total_fu vector - vector of fu times
  #     *all vectors should be same length with indexes corresponding to same
  #       person you know what I mean
  # type: string: specify what type ofplot to make
  # title - a overall title for the calibration plot
  # OUTPUTS: a ggplot2 object
  
  #initializing calibration list
  dt_list = list()
  
  #creating the dataframes to plot
  for (i in 1:length(risk_vectors)){
    
    #extracting prediction
    pred = risk_vectors[[i]]
    
    #number of deciles
    n = 10
    
    ##determining which quintile an observation falls into
    #lots of duplicated risk results so quantiles are not of equal length
    deciles = quantile(pred, probs = (0:n)/n)
    dec_group = cut(pred, 
                    breaks = deciles, 
                    labels = 1:n,
                    include.lowest = TRUE,
                    right = TRUE)
    
    #initializing dataframe to populate
    df = data.frame(matrix(nrow = n, ncol = 4))
    names(df) = c("pred","obs","decile","per")
    
    for (j in 1:n){
      
      #group indexes
      group_index = (dec_group == j)
      
      #predicted risk in each decile
      df[j,"pred"] = mean(pred[group_index])
      
      # #observed events use kaplan meier estimate 
      KM_estimate = max(1-survfit(Surv(total_fu_vector[group_index],
                                       indicator_vector[group_index]) ~ 1,
                                  )$surv)
      df[j,"obs"] = KM_estimate*100
      
      # #just raw estimate instead
      # df[j,"obs"] = mean(indicator_vector[group_index]) * 100
      
      #risk decile
      df[j,"decile"] = j
      
      #proportion of total events
      df[j,"per"] = sum(indicator_vector[group_index])/sum(indicator_vector)*100
    }
    
    #debugging print statement
    print(df)
    
    #convert to data.table
    dt = data.table(df)
    
    #add to calibration list
    dt_list = append(dt_list,list(dt))
  }
  
  #extracting n
  n = length(total_fu_vector)
  
  #checking which type of plot to create
  if(and(grepl(type,substring("calibration",0,nchar(type))),
          !grepl(type,substring("discrimination",0,nchar(type))))){
    
    # Running cal_plot_gg function
    plot = cal_plot_gg(dt_list,
                       title,
                       names(risk_vectors),
                       subtitle)
    
  } else if (and(!grepl(type,substring("calibration",0,nchar(type))),
                 grepl(type,substring("discrimination",0,nchar(type))))){
      
    # Running disc_plot_gg function
    plot = dis_plot_gg(dt_list,
                       title,
                       names(risk_vectors),
                       subtitle)
  }

  #return plot output
  return(plot)
  
}


validation = function(pred_dt,
                      type = c("calibration","discrimination"),
                      title = NULL,
                      display_n = FALSE,
                      display_metrics = FALSE,
                      return = 0){
  
  #wrapper function for calibration and discrimination plots
  #
  # Inputs:
  # display_n = ?display number of individuals
  # display_metrics = ?calculate and display boostrap metrics for this pred_dt
  #
  #
  # Outputs:
  # a plot

  #copy pred_dt for safety
  dt = copy(pred_dt)
  
  #remove index
  dt[,VSIMPLE_INDEX_MASTER := NULL]
  
  # extract outcome vectors
  #indexes = pred_dt$VSIMPLE_INDEX_MASTER
  total_fu = dt$total_fu
  dt[,total_fu := NULL]
  indicator = dt$indicator
  dt[,indicator := NULL]
  # extract risk vectors (extract every column)
  risk_vectors = list()
  for (model_name in names(dt)){
    risk_vector = list(dt[,get(model_name)])
    names(risk_vector) = model_name
    risk_vectors = append(risk_vectors,risk_vector)
  }

  # generating subtitle
  subtitle = ""
  if(display_n == TRUE){
    n = dim(dt)[[1]]
    subtitle = paste0(subtitle,"n: ",n)
  } 
  if(display_metrics ==TRUE){
    # calculate boostrap metrics
    boot_output = pred_dt_bootstrap_wrapper(pred_dt)
    subtitle = paste0(subtitle,"\n",boot_output)
  }
  
  #pass vectors into validation_plot function  
  plot = validation_plot(risk_vectors,
                         indicator,
                         total_fu,
                         type,
                         title = title,
                         subtitle = subtitle)
  
  #returning or displaying plot
  if (return == 1){
    return(plot)
  } else {
    #display plot if not returned
    plot
  }
  
}

###########################################################################################
# full plotting function wrappers ####
###########################################################################################

#whole cohort plots
whole_cohort_plot_generator = function(pred_dt_list,
                                      type = "calibration",
                                      display_n = TRUE,
                                      display_metrics = TRUE){
  whole_plot_list = list()
  for (i in 1:2){
    #whole cohort calibration plots
    plot = validation(pred_dt_list[[i]],
                          type = type,
                          str_to_title(names(pred_dt_list)[[i]]),
                          return = 1,
                          display_n = display_n,
                          display_metrics = display_metrics)
    plot
    
    whole_plot_list = append(whole_plot_list,list(plot))
    
    plot_path = paste0(output_dir,"plots/",type,"/")
    if (!dir.exists(plot_path)) {
      dir.create(plot_path,recursive = TRUE) 
    }
    ggsave(paste(names(pred_dt_list)[[i]],".jpg",sep=""),
          plot = plot,
          path = plot_path,
          height = 5,
          width = 6.2)
  }
  whole_plot = ggarrange(plotlist = whole_plot_list,
                           ncol = 2)
  whole_plot
  ggsave(paste0("whole_",type,"_plot.png"),
        plot = whole_plot,
        path = plot_path,
        height = 5,
        width = 15)
  return(whole_plot_list)
}

#demo table and plot list
cal_sub_pop_plot_generator = function(pred_dt_list,
                                  all_conds,
                                  display_n = TRUE,
                                  display_metrics = TRUE){
  # function to generator and save all the subpopulation plots and demo tables
  # Inputs:
  # pred_dt_list: list of pred_dt's for each gender
  # all_conds: list of conditions to loop through
  #
  # Outputs:
  # full_cal_plot_list: list of ggplot2 objects
  
  #initialize lists
  sub_pop_demo_table_list = list()
  full_cal_plot_list = list()
  
  #looping through conditions
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
        
        #demo table indexes
        demo_subpop_index = sub_pop_pred_dt[,VSIMPLE_INDEX_MASTER]
        
        #extract plot name
        plot_name = eval(parse(text = names(conds)[cond_ind]))
        
        #save to index list
        demo_subpop_index_list = list(demo_subpop_index)
        names(demo_subpop_index_list) = plot_name
        index_list = append(index_list,demo_subpop_index_list)
        
        #plot using name of condition
        cal_plot = validation(sub_pop_pred_dt,
                              type = "cal",
                              title = plot_name,
                              display_n = display_n,
                              display_metrics = display_metrics)
        
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
      
      #demo table
      demo_table_path = paste0(output_dir,"tables/cc_subpop_demo/")
      if (!dir.exists(demo_table_path)) {
        dir.create(demo_table_path,recursive = TRUE) 
      }
      demo_table = demo_table_generator(index_list,test_rows = TRUE)
      write.csv(demo_table,
                paste0(demo_table_path,gender,conds_name,".csv"))
      
      #save plot
      output_plot_path = paste0(output_dir,"plots/calibration/subpopulations/")
      if (!dir.exists(output_plot_path)) {
        dir.create(output_plot_path,recursive = TRUE) 
      } 
      ggsave(paste(conds_name,gender,".png",sep = ""),
            arrange_plot,
            path = output_plot_path,
            height = ceiling(length(cal_plot_list)/2)*5,
            width = 9)
    }
  }
  return(full_cal_plot_list)
}

###########################################################################################
# Sub-setting functions #### 
###########################################################################################

# outlining subpopulations ####
#ethnic sub populations
eth_conds = c("sprintf(\"Maaori %s\",gender)" = 
                "eth_cat == \"NZM\"",
              "sprintf(\"Pacific %s\",gender)" = 
                "eth_cat == \"Pacific\"",
              "sprintf(\"Indian %s\",gender)" = 
                "eth_cat == \"Indian\"",
              "sprintf(\"%s of other ethnicities\",str_to_title(gender))" = 
                "eth_cat == \"Other\"",
              "sprintf(\"Chinese %s\",gender)" = 
                "eth_cat == \"Chinese\"",
              "sprintf(\"European %s\",gender)" = 
                "eth_cat == \"NZE\"")

#age based subp opulations
age_conds = c("sprintf(\"%s 30-44\",str_to_title(gender))" = 
                "and(nhi_age >= 30,nhi_age < 45)",
              "sprintf(\"%s 45-59\",str_to_title(gender))" =
                "and(nhi_age >= 45,nhi_age < 60)",
              "sprintf(\"%s 60-74\",str_to_title(gender))" =
                "and(nhi_age >= 60,nhi_age < 75)")

#deprivation sub populations
dep_conds = c("sprintf(\"%s Deprivation Q1\",str_to_title(gender))" = 
                "c_en_nzdep_q == -2",
              "sprintf(\"%s Deprivation Q2\",str_to_title(gender))" =
                "c_en_nzdep_q == -1",
              "sprintf(\"%s Deprivation Q3\",str_to_title(gender))" =
                "c_en_nzdep_q == -0",
              "sprintf(\"%s Deprivation Q4\",str_to_title(gender))" =
                "c_en_nzdep_q == 1",
              "sprintf(\"%s Deprivation Q5\",str_to_title(gender))" =
                "c_en_nzdep_q == 2")   

#baseline_dispensing sub populations
bd_conds = c("sprintf(\"%s with BPL Meds\",str_to_title(gender))" = 
               "ph_bp_lowering_prior_6mths == 1",
             "sprintf(\"%s without BPL Meds\",str_to_title(gender))" = 
               "ph_bp_lowering_prior_6mths == 0",
             "sprintf(\"%s with LL Meds\",str_to_title(gender))" = 
               "ph_lipid_lowering_prior_6mths == 1",
             "sprintf(\"%s without LL Meds\",str_to_title(gender))" = 
               "ph_lipid_lowering_prior_6mths == 0",
             "sprintf(\"%s with AT Meds\",str_to_title(gender))" = 
               "ph_antithrombotic_prior_6mths == 1",
             "sprintf(\"%s without LL Meds\",str_to_title(gender))" = 
               "ph_antithrombotic_prior_6mths == 0")

dhb_conds = c("sprintf(\"%s in Auckland\",str_to_title(gender))" = 
                "DHB_name == \"Auckland\"",
              "sprintf(\"%s in Counties Manukau\",str_to_title(gender))" = 
                "DHB_name == \"Counties Manukau\"",
              "sprintf(\"%s in Waitemata\",str_to_title(gender))" = 
                "DHB_name == \"Waitemata\"",
              "sprintf(\"%s in Northland\",str_to_title(gender))" = 
                "DHB_name == \"Northland\"")

all_conds = list("eth" = eth_conds,
                 "age" = age_conds,
                 "dep" = dep_conds,
                 "bd" = bd_conds,
                 "dhb" = dhb_conds)



subpop_index = function(condition){
  #takes a condition and returns the indexes of all patients in 
  #the an cohort who satisfy this condition
  # Inputs:
  # Condition: string: formula to be passed into i section of an_cohort dt
  # 
  # Outputs:
  # A vector of VSIMPLE_INDEX_MASTER indexes 
  
  if(!exists("an_cohort_indexing")){
    an_cohort_indexing = read.fst("data/cohorts/an_cohort.fst",as.data.table = TRUE)
  }
  
  #actual code
  indexes = an_cohort_indexing[eval(parse(text = condition)),VSIMPLE_INDEX_MASTER]
  
  return(indexes)
}

###########################################################################################
# reclassification testing ####
###########################################################################################

reclassification_testing = function(pred_dt,
                                  high_risk_cutoff = 15,
                                  low_risk_cutoff = 5,
                                  print_results = TRUE){
  # input is a pred_dt with the following columns:
  # VSIMPLE_INDEX_MASTER, total_fu, New Model, Baseline Model
  # output is a list of 3 objectsi:
  # 1. summary of the new and old risk categories
  # 2. cox model between became high risk and lost high risk
  # 3. cox model between became low risk and lost low risk

  # add in follow-up time by joining on VSIMPLE_INDEX_MASTER
  dt = copy(pred_dt[an_cohort_ts,on = .(VSIMPLE_INDEX_MASTER),total_fu := i.total_fu])

  dt[,`:=`(old_high_risk = as.integer(`Baseline Model` > high_risk_cutoff),
          new_high_risk = as.integer(`New Model` > high_risk_cutoff),
          old_low_risk = as.integer(`Baseline Model` < low_risk_cutoff),
          new_low_risk = as.integer(`New Model` < low_risk_cutoff))]
  dt[,`:=`(became_high_risk = as.integer(old_high_risk == 0 & new_high_risk == 1),
          lost_high_risk = as.integer(old_high_risk == 1 & new_high_risk == 0),
          became_low_risk = as.integer(old_low_risk == 0 & new_low_risk == 1),
          lost_low_risk = as.integer(old_low_risk == 1 & new_low_risk == 0))]

  # Calculate mean and sum for each column
  temp_dt = copy(dt[,-c("VSIMPLE_INDEX_MASTER","total_fu","New Model","Baseline Model")])
  summary_dt = data.table(
    column = names(temp_dt),
    mean = sapply(temp_dt, mean),
    sum = sapply(temp_dt, sum),
    n_events = sapply(names(temp_dt), function(x) sum(dt[,"indicator"] == 1 & dt[[x]] == 1)),
    prob_events = format(
      sapply(names(temp_dt), function(x) (sum(dt[,indicator] == 1 & dt[[x]] == 1) / sum(dt[[x]]==1)) ),
      scientific = FALSE,
      digits = 3
    )
  )

  # fit a cox model between the became high risk and the lost high risk to see the hazard ratio of being in 
  # the became high risk vs the lost high risk  
  temp_model_high = coxph(Surv(total_fu, indicator) ~ became_high_risk, 
                    data = dt[became_high_risk == 1 | lost_high_risk == 1])
  #became high risk are ##% higher risk of death than lost high risk

  #repeat for low risk
  temp_model_low = coxph(Surv(total_fu, indicator) ~ became_low_risk, 
                    data = dt[became_low_risk == 1 | lost_low_risk == 1])
  #became low risk are ##% lower risk of death than lost low risk

  if (print_results){
    print(summary_dt)
    print(summary(temp_model_high))
    print(summary(temp_model_low))
  }
  
  output_table_path = paste0(output_dir,"tables/") 
  if (!dir.exists(output_table_path)) {
    dir.create(output_table_path,recursive = TRUE) 
  }
  write.csv(summary_dt,
            paste0(output_table_path,"reclassification_summary.csv"))
  return(list(list(summary_dt),list(summary(temp_model_high)),list(summary(temp_model_low))) )
}


















###################################################################################
# . MISC
###########################################################################
DtoR2 = function(D){
  # kappa = sqrt(8/pi)
  kappa = sqrt(8/pi)
  (D**2 / kappa**2) / ((D**2 / kappa**2) + (pi**2 / 6))
}



summary_stats = function(lp_dt_list){
  stats_table_list = list()
  for (gender_ind in 1:2){
    summary_stats = c("D",
                        "D_se",
                        "string_D",
                        "R2",
                        "R2_se",
                        "string_R2",
                        "C",
                        "C_se",
                        "string_C")
    stats_table = data.frame(matrix(NA,
                                    nrow = length(summary_stats),
                                    ncol = 3))
    colnames(stats_table) = c("New Model",
                              "Baseline Model",
                              "Difference")
    rownames(stats_table) = summary_stats

    lp_dt = lp_dt_list[[gender_ind]]
    model_names = c("New Model","Baseline Model")
    for (i in 1:2){
      model_name = model_names[i]
      dt = copy(lp_dt[,c("total_fu","indicator",model_name),with = FALSE])
      dt[,lp := get(model_name)]
      # remove model_name colum  from dt
      # lp = dt[[model_name]]
      # dt = dt[,.(total_fu,indicator)]
      # dt[,lp := lp]
      dt = dt[order(dt$lp),]
      dt[,prob_level := ((1:nrow(dt)) - (3/8)) / (nrow(dt) + 1/4)]  # in reality the 3/8 and 1/4 are negligible ?unsure why included
      # now for a blom transformation - rank based inverse normal transformation
      dt[,Z_original := qnorm(prob_level) / sqrt(8/pi)] #this is a standard normal distribution
      # now to deal with ties
      dt[,Z := mean(Z_original),by = lp]
      D_model = coxph(Surv(total_fu,indicator) ~ Z,
                    data = dt,
                    ties = "breslow",
                    model = TRUE)
      D = D_model$coef
      D_se = summary(D_model)$coef[,"se(coef)"]
      D_ll = D - 1.96*D_se
      D_ul = D + 1.96*D_se
      string_D = sprintf("%.4f (%.4f - %.4f)",
                          D,
                          D_ll,
                          D_ul)
      stats_table[which(rownames(stats_table) == "D"),i] = D
      stats_table[which(rownames(stats_table) == "D_se"),i] = D_se
      stats_table[which(rownames(stats_table) == "string_D"),i] = string_D

      R2 = DtoR2(D)
      R2_se = sqrt( (64 * D**2 * D_se**2) / ((D**2 + 4)**4) ) #delta method
      R2_ll = DtoR2(D - 1.96*R2_se)
      R2_ul = DtoR2(D + 1.96*R2_se)
      string_R2 = sprintf("%.4f (%.4f - %.4f)",
                          R2,
                          R2_ll,
                          R2_ul)
      stats_table[which(rownames(stats_table) == "R2"),i] = R2
      stats_table[which(rownames(stats_table) == "R2_se"),i] = R2_se
      stats_table[which(rownames(stats_table) == "string_R2"),i] = string_R2

      C_index = concordance.index(dt$lp,dt$total_fu,dt$indicator)
      C = C_index$c.index
      #C_se = C_index$se
      C_ll = C_index$lower
      C_ul = C_index$upper
      C_se = (C_ul - C_ll) / (2 * 1.96)
      string_C = sprintf("%.4f (%.4f - %.4f)",
                          C,
                          C_ll,
                          C_ul)
      stats_table[which(rownames(stats_table) == "C"),i] = C
      stats_table[which(rownames(stats_table) == "C_se"),i] = C_se
      stats_table[which(rownames(stats_table) == "string_C"),i] = string_C
    }
    stats_table

    #calculate difference
    stats_table[which(rownames(stats_table) == "D"),3] = 
      as.numeric(stats_table[which(rownames(stats_table) == "D"),1]) - 
      as.numeric(stats_table[which(rownames(stats_table) == "D"),2])
    stats_table[which(rownames(stats_table) == "R2"),3] = 
      as.numeric(stats_table[which(rownames(stats_table) == "R2"),1]) - 
      as.numeric(stats_table[which(rownames(stats_table) == "R2"),2])
    stats_table[which(rownames(stats_table) == "C"),3] = 
      as.numeric(stats_table[which(rownames(stats_table) == "C"),1]) - 
      as.numeric(stats_table[which(rownames(stats_table) == "C"),2])

    # calculate se(difference)
    stats_table[which(rownames(stats_table) == "D_se"),3] = 
      sqrt(as.numeric(stats_table[which(rownames(stats_table) == "D_se"),1])**2 + 
      as.numeric(stats_table[which(rownames(stats_table) == "D_se"),2])**2)
    stats_table[which(rownames(stats_table) == "R2_se"),3] = 
      sqrt(as.numeric(stats_table[which(rownames(stats_table) == "R2_se"),1])**2 + 
      as.numeric(stats_table[which(rownames(stats_table) == "R2_se"),2])**2)
    stats_table[which(rownames(stats_table) == "C_se"),3] = 
      sqrt(as.numeric(stats_table[which(rownames(stats_table) == "C_se"),1])**2 + 
      as.numeric(stats_table[which(rownames(stats_table) == "C_se"),2])**2)

    # fill in string difference and include p value
    # assuming normal distribution
    stats_table[which(rownames(stats_table) == "string_D"),3] = 
      sprintf("%.4f (%.4f - %.4f) p = %.4f",
              as.numeric(stats_table[which(rownames(stats_table) == "D"),3]),
              as.numeric(stats_table[which(rownames(stats_table) == "D"),3]) - 1.96* as.numeric(stats_table[which(rownames(stats_table) == "D_se"),3]),
              as.numeric(stats_table[which(rownames(stats_table) == "D"),3]) + 1.96* as.numeric(stats_table[which(rownames(stats_table) == "D_se"),3]),
              pnorm(0,as.numeric(stats_table[which(rownames(stats_table) == "D"),3]),as.numeric(stats_table[which(rownames(stats_table) == "D_se"),3])))
    stats_table[which(rownames(stats_table) == "string_R2"),3] = 
      sprintf("%.4f (%.4f - %.4f) p = %.4f",
              as.numeric(stats_table[which(rownames(stats_table) == "R2"),3]),
              as.numeric(stats_table[which(rownames(stats_table) == "R2"),3]) - 1.96*as.numeric(stats_table[which(rownames(stats_table) == "R2_se"),3]),
              as.numeric(stats_table[which(rownames(stats_table) == "R2"),3]) + 1.96*as.numeric(stats_table[which(rownames(stats_table) == "R2_se"),3]),
              pnorm(0,as.numeric(stats_table[which(rownames(stats_table) == "R2"),3]),as.numeric(stats_table[which(rownames(stats_table) == "R2_se"),3])))
    stats_table[which(rownames(stats_table) == "string_C"),3] = 
      sprintf("%.4f (%.4f - %.4f) p = %.4f",
              as.numeric(stats_table[which(rownames(stats_table) == "C"),3]),
              as.numeric(stats_table[which(rownames(stats_table) == "C"),3]) - 1.96*as.numeric(stats_table[which(rownames(stats_table) == "C_se"),3]),
              as.numeric(stats_table[which(rownames(stats_table) == "C"),3]) + 1.96*as.numeric(stats_table[which(rownames(stats_table) == "C_se"),3]),
              pnorm(0,as.numeric(stats_table[which(rownames(stats_table) == "C"),3]),as.numeric(stats_table[which(rownames(stats_table) == "C_se"),3])))

    list_stats_table = list(stats_table)
    names(list_stats_table) = names(lp_dt_list)[gender_ind]
    stats_table_list = append(stats_table_list,list_stats_table)
  }
  return(stats_table_list)
}





# subsetting_helper = function(dt,
#                              conds,
#                              func){
#   
#   # function which subsets a dt according to a "list" of conditions 
#   # and runs a function on each subset returning a list of outputs
#   
#   #copy dt for safety  
#   dt = copy(dt)
#   
#   #initialize output list
#   output_list = list()
#   
#   for (cond_ind in 1:length(conds)){
#     
#     # subsetting the data table
#     subset_dt = dt[VSIMPLE_INDEX_MASTER %in% subpop_index(conds[[cond_ind]]),]
#     
#     # extracting names from conditions
#     name = eval(parse(text = names(conds)[cond_ind]))
#     
#     #running the function on the subset
#     output = func(subset_dt,name)
#     
#     #append function to output_list
#     output_list = append(output_list,list(output))
#     
#     
#   }
#   
#   #return the output list
#   return(output_list)
#   
# }
  
  




# # TestSafe distribution ####
# test_dist_plot = function(input_dt){
#   # takes as inputs a data.table and retursn the test distribution in that data atble
  
  
#   tests = c("egfr","hba1c","tchdl","tri")
#   test_clean = c("eGFR","HbA1c","TC/HDLC ratio","Triglycerides")
#   tests_verbose = c("eGFR (mL/min/1.73m^2)",
#                     "HbA1c (mmol/mol)",
#                     "TC/HDLC ratio",
#                     "Triglycerides (mmol/L)")
#   test_limits = list(c(20,150),
#                      c(20,80),
#                      c(0,10),
#                      c(0,8))
  
  
  
  
#   dt = as_tibble(copy(input_dt))
  
#   egfr_plot = dt %>%
#     filter(!is.na(get(tests[1]))) %>%
#     {ggplot(.) +
#         geom_density(aes_string(x = tests[1]),
#                      fill = "#F8766D", 
#                      size = 0.5,
#                      adjust = 1) + 
#         scale_x_continuous(limits = test_limits[[1]])+
#         labs(x = tests_verbose[1],
#              y = NULL,
#              #title = paste("Distribution of",test_clean[1]),
#              subtitle = paste("n:",as.character(dim(.)[1]))) +
#         theme_bw()}
  
#   hba1c_plot = dt %>%
#     filter(!is.na(get(tests[2]))) %>%
#     {ggplot(.) +
#         geom_density(aes_string(x = tests[2]),
#                      fill = "#619CFF", 
#                      size = 0.5,
#                      adjust = 3) + 
#         scale_x_continuous(limits = test_limits[[2]])+
#         labs(x = tests_verbose[2],
#              y = NULL,
#              #title = paste("Distribution of",test_clean[2]),
#              subtitle = paste("n:",as.character(dim(.)[1]))) +
#         theme_bw()}
  
#   tchdl_plot = dt %>%
#     filter(!is.na(get(tests[3]))) %>%
#     {ggplot(.) +
#         geom_density(aes_string(x = tests[3]),
#                      fill = "#00BA38", 
#                      size = 0.5,
#                      adjust = 3) + 
#         scale_x_continuous(limits = test_limits[[3]])+
#         labs(x = tests_verbose[3],
#              y = NULL,
#              #title = paste("Distribution of",test_clean[3]),
#              subtitle = paste("n:",as.character(dim(.)[1]))) +
#         theme_bw()}
  
#   tri_plot = dt %>%
#     filter(!is.na(get(tests[4]))) %>%
#     {ggplot(.) +
#         geom_density(aes_string(x = tests[4]),
#                      fill = "#C77CFF", 
#                      size = 0.5,
#                      adjust = 3) + 
#         scale_x_continuous(limits = test_limits[[4]])+
#         labs(x = tests_verbose[4],
#              y = NULL,
#              #title = paste("Distribution of",test_clean[4]),
#              subtitle = paste("n:",as.character(dim(.)[1]))) +
#         theme_bw()}
  
#   patchwork = egfr_plot | hba1c_plot | tchdl_plot | tri_plot
#   final_plot = patchwork + 
#     plot_annotation(#title = "Distribution of test results",
#                     theme = theme(plot.title = element_text(
#                       face = "bold",
#                       size = 20)))
  
#   return(final_plot)
  
# }










# star_plot
# dt = copy(pred_dt_list[[1]])
# dt[,`Old model` := NULL]
# dt[,VSIMPLE_INDEX_MASTER := NULL]
# dt[,total_fu := NULL]
# dt[,index := 1:dim(dt)[1]]
# star_plot = ggplot(dt[order(indicator)]) + 
#   geom_point(aes(y = index, x = `New model`, color = factor(indicator), size = factor(indicator),alpha = factor(indicator)),shape = 20) +
#   scale_alpha_manual(values = c(0.5,1)) +
#   scale_size_manual(values = c(1,1)) +
#   scale_color_manual(values = c("grey","red")) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# ggsave("starplot.png",
#        star_plot,
#        path = "plots/misc")
# #now ordered
# dt[order(`New model`),index := 1:dim(dt)[1]]
# dt[order(`New model`)]
# 
# star_plot = ggplot(dt[order(indicator)]) + 
#   geom_point(aes(y = index, x = `New model`, 
#                  color = factor(indicator), 
#                  size = factor(indicator),
#                  alpha = factor(indicator)),shape = 20) +
#   scale_alpha_manual(values = c(0.5,0.5)) +
#   scale_size_manual(values = c(1,1)) +
#   scale_color_manual(values = c("grey","red")) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# star_plot
# ggsave("starplot2.png",
#        star_plot,
#        path = "plots/misc")
# CREATE A DISTRIBUTION STARPLOT!!! - comparing distribution o risk scores for those who get disease and those who dont!


# manual_c_stat = function(){
#   ##########UNFINISHED QUALITY OF LIFE FUNCTION
#   sub_pop_pred_dt[indicator == 1]

#   sub_pop_pred_dt
#   pos_predictions = sub_pop_pred_dt[indicator == 1]
#   neg_predictions = sub_pop_pred_dt[indicator == 0]


#   pred = sub_pop_pred_dt[,"New Model"][[1]]
#   obs = sub_pop_pred_dt[,"indicator"][[1]]

#   pos_predictions = pred[obs == 1]
#   neg_predictions = pred[obs == 0]

#   all_pairs = expand.grid(pos_predictions,neg_predictions)
#   #var1 is pos, var2 is neg
#   concordant_pairs = sum(all_pairs$Var1 > all_pairs$Var2)
#   discordant_pairs = sum(all_pairs$Var1 < all_pairs$Var2)
#   ties = dim(all_pairs)[1] - concordant_pairs - discordant_pairs
#   c = (concordant_pairs + 0.5 * ties)/dim(all_pairs)[1]
# }

