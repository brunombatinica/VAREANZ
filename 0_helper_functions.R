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
library(mfp)
library(survminer)
library(glmnet)
#library(Matrix)
library(stringr)
library(grid)
library(Gmisc)
#library(rms)
library(patchwork)

# 1. Data table helper functions ####

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

# 2. baseline demographics ####
# an_cohort = read.fst("data/cohorts/an_cohort_ts.fst",
#                      as.data.table = TRUE)



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


#hc = read.fst("data/raw_data/VARIANZ2012/VARIANZ_2012.fst")



# 3. Visualize model coefficients ####

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
                             "TestSafe",
                             "TestSafe",
                             "TestSafe",
                             "TestSafe",
                             "TestSafe"),
                             levels = c("TestSafe",
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










# 4. prediction functions ####

# an_cohort_raw
# test = coxph(Surv(total_fu,indicator) ~ hx_af + gender_code,
#              data = an_cohort_raw,
#              ties ="breslow",
#              model = TRUE)
# -max(basehaz(test)$hazard)
# s = survfit(test, type = )
# -max(s$cumhaz)


base_surv = function(model){
  #returns the baseline 5 year survival of the model
  
  bs = exp(-max(basehaz(model,centered = FALSE)$hazard))
  
  return(bs)
}

prediction_helper = function(model,
                             dt){
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
  
  #converting linear predictor into 5 year event estimate
  pred = (1 - (bs^exp(lp)) ) * 100
  
  #return pred
  return(pred)
  
}


prediction_wrapper = function(model_data,
                             new_data,
                             name){
  # wrapper for prediction_helper to enable pattern submodel prediction
  
  #Inputs:
  # models: either a coxph object or a list of coxph objects
  # new_data: either a dt or list of dt's
  # !assume both are single object or both ar list
  
  if (class(model_data) == "coxph"){
    #run simple prediction
    
    pred = prediction_helper(model_data,new_data)
    output = new_data[,c("VSIMPLE_INDEX_MASTER")][,(name) := pred]
    
    
    
  } else if (class(model_data) == "list"){
    
    outputs = list()
    for (i in 1:length(model_data)){
      pred = prediction_helper(model_data[[i]],new_data[[i]])
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
                          new_data){
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
                             model_name)
    
    #joining prediciton to pred_dt
    pred_dt = pred_dt[predictions]
    
  }
  
  #return pred_dt
  return(pred_dt)
}




# Discrimination plot ####

dis_plot_gg = function(dt_list,
                       plot_title,
                       model_names,
                       n){
  
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
           subtitle = paste("n:",n,sep = " ")) +
      scale_x_continuous(name = "Deciles of predicted risk",
                         limits = c(0,10),
                         breaks = 1:10) +
      ylab("Observed events (%)") +
      scale_shape_manual(values = c(17,16)) +
      scale_color_manual(values = c("#F8766D","#00BFC4")) +
      theme_bw() +
      theme(plot.title = element_text(size = 14, face = "bold"),
            legend.position = c(0.23,0.84),   #legend formatting options
            legend.title = element_blank(),
            legend.background = element_rect(size = 0,
                                             linetype = "solid",
                                             colour = "white"),
            plot.margin = unit(c(0.08,0.2,0.08,0.2),"npc")) 
    
  } else if(length(dt_list) == 3){
    
    print("TESTING")
    
    plot_data$model = factor(plot_data$model, levels = c("Normal Imputation","CCS","Baseline Model"))
    print.data.frame(plot_data)
    
    
    
    plot = ggplot(plot_data) +
      geom_point(aes(x = decile, y = per, color = model,shape = model), size = 2, stroke = 1.2) +
      ylim(0,5*(max(plot_data[,5]) %/% 5 + 1)) +
      labs(title = paste("Discrimination:",plot_title),
           subtitle = paste("n:",n,sep = " ")) +
      scale_x_continuous(name = "Deciles of predicted risk",
                         limits = c(0,10),
                         breaks = 1:10) +
      ylab("Observed events (%)") +
      scale_shape_manual(values = c(3,4,20)) +
      scale_color_manual(values = c("#F8766D","#00BA38","dark grey")) +
      theme_bw() +
      theme(plot.title = element_text(size = 14, face = "bold"),
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

# Calibration plot ####

cal_plot_gg = function(dt_list,
                       plot_title,
                       model_names,
                       n){
  
  
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
           subtitle = paste("n:",n,sep = " "),
           x = "Observed events (%)",
           y = "Mean predicted 5 year risk (%)") +
      scale_shape_manual(values = c(17,16)) +
      scale_color_manual(values = c("#F8766D","#00BFC4")) +
      theme_bw() +
      theme(plot.title = element_text(size = 14, face = "bold"),
            legend.position = c(0.23,0.87),   #legend formatting options
            legend.title = element_blank(),
            legend.background = element_rect(size = 0,
                                             linetype = "solid",
                                             colour = "white"),
            plot.margin = unit(c(0.08,0,0.08,0),"npc")) 
    
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
           subtitle = paste("n:",n,sep = " "),
           x = "Observed events (%)",
           y = "Mean predicted 5 year risk (%)") +
      scale_shape_manual(values = c(3,4,20)) +
      scale_color_manual(values = c("#F8766D","#00BA38","dark grey")) +
      theme_bw() +
      theme(plot.title = element_text(size = 14, face = "bold"),
            legend.position = c(0.23,0.85),   #legend formatting options
            legend.title = element_blank(),
            legend.background = element_rect(size = 0,
                                             linetype = "solid",
                                             colour = "white"),
            plot.margin = unit(c(0.08,0,0.08,0),"npc"))
    

  }
  
  
  
  
  return(plot)
  
}

# Validation wrapper functions ####
validation_plot = function(risk_vectors,
                           indicator_vector,
                           total_fu_vector,
                           type,
                           title = NULL){
  
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
                       n)
    
  } else if (and(!grepl(type,substring("calibration",0,nchar(type))),
                 grepl(type,substring("discrimination",0,nchar(type))))){
      
    # Running disc_plot_gg function
    plot = dis_plot_gg(dt_list,
                       title,
                       names(risk_vectors),
                       n)
  }

  #return plot output
  return(plot)
  
}



validation = function(pred_dt,
                      type = c("calibration","discrimination"),
                      title = NULL,
                      return = 0){
  
  #wrapper function for calibration and discrimination plots
  
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
  
  #pass vectors into validation_plot function  
  plot = validation_plot(risk_vectors,
                         indicator,
                         total_fu,
                         type,
                         title = title)
  
  #returning or displaying plot
  if (return == 1){
    return(plot)
  } else {
    #display plot if not returned
    plot
  }
  
}


#Sub-setting functions #### 

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
  # A vector of indexes 
  
  if(!exists("an_cohort_indexing")){
    an_cohort_indexing = read.fst("data/cohorts/an_cohort.fst",as.data.table = TRUE)
  }
  
  #actual code
  indexes = an_cohort_indexing[eval(parse(text = condition)),VSIMPLE_INDEX_MASTER]
  
  return(indexes)
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
  
  




# TestSafe distribution ####


test_dist_plot = function(input_dt){
  # takes as inputs a data.table and retursn the test distribution in that data atble
  
  
  tests = c("egfr","hba1c","tchdl","tri")
  test_clean = c("eGFR","HbA1c","TC/HDLC ratio","Triglycerides")
  tests_verbose = c("eGFR (mL/min/1.73m^2)",
                    "HbA1c (mmol/mol)",
                    "TC/HDLC ratio",
                    "Triglycerides (mmol/L)")
  test_limits = list(c(20,150),
                     c(20,80),
                     c(0,10),
                     c(0,8))
  
  
  
  
  dt = as_tibble(copy(input_dt))
  
  egfr_plot = dt %>%
    filter(!is.na(get(tests[1]))) %>%
    {ggplot(.) +
        geom_density(aes_string(x = tests[1]),
                     fill = "#F8766D", 
                     size = 0.5,
                     adjust = 1) + 
        scale_x_continuous(limits = test_limits[[1]])+
        labs(x = tests_verbose[1],
             y = NULL,
             #title = paste("Distribution of",test_clean[1]),
             subtitle = paste("n:",as.character(dim(.)[1]))) +
        theme_bw()}
  
  hba1c_plot = dt %>%
    filter(!is.na(get(tests[2]))) %>%
    {ggplot(.) +
        geom_density(aes_string(x = tests[2]),
                     fill = "#619CFF", 
                     size = 0.5,
                     adjust = 3) + 
        scale_x_continuous(limits = test_limits[[2]])+
        labs(x = tests_verbose[2],
             y = NULL,
             #title = paste("Distribution of",test_clean[2]),
             subtitle = paste("n:",as.character(dim(.)[1]))) +
        theme_bw()}
  
  tchdl_plot = dt %>%
    filter(!is.na(get(tests[3]))) %>%
    {ggplot(.) +
        geom_density(aes_string(x = tests[3]),
                     fill = "#00BA38", 
                     size = 0.5,
                     adjust = 3) + 
        scale_x_continuous(limits = test_limits[[3]])+
        labs(x = tests_verbose[3],
             y = NULL,
             #title = paste("Distribution of",test_clean[3]),
             subtitle = paste("n:",as.character(dim(.)[1]))) +
        theme_bw()}
  
  tri_plot = dt %>%
    filter(!is.na(get(tests[4]))) %>%
    {ggplot(.) +
        geom_density(aes_string(x = tests[4]),
                     fill = "#C77CFF", 
                     size = 0.5,
                     adjust = 3) + 
        scale_x_continuous(limits = test_limits[[4]])+
        labs(x = tests_verbose[4],
             y = NULL,
             #title = paste("Distribution of",test_clean[4]),
             subtitle = paste("n:",as.character(dim(.)[1]))) +
        theme_bw()}
  
  patchwork = egfr_plot | hba1c_plot | tchdl_plot | tri_plot
  final_plot = patchwork + 
    plot_annotation(#title = "Distribution of test results",
                    theme = theme(plot.title = element_text(
                      face = "bold",
                      size = 20)))
  
  return(final_plot)
  
}


sub_pop_test_dist_plot = function(){
  # TO DO
}
