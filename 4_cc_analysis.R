#CC analysis using pre 2015 data cut-off

# Clearing workspace
rm(list = ls())
dev.off()

# Loading in necessary libraries and data ################################
source("0_helper_functions.R")

#set working directory
setwd("//uoa.auckland.ac.nz/Shared/MED/EPBI/ViewData/Users/bbat644/Desktop/code")

#source helper functions
source("0_helper_functions.R")

#1. import data and cleaning up data ####
an_cohort = read.fst("data/cohorts/an_cohort_ts.fst",
                     as.data.table = TRUE)


setkey(an_cohort,VSIMPLE_INDEX_MASTER)
#remove DHB
an_cohort[,DHB_name := NULL]
cohort = an_cohort[complete.cases(an_cohort)]
dim(cohort)





#removing index
#cohort[,VSIMPLE_INDEX_MASTER := NULL]
#removing result dates
#cohort[,c("egfr_date","hba1c_date","tchdl_date"):= NULL]
#cut up factors

#age distribution of cohort

# predictors disctibution by age




cohort = cut_factors(cohort)

# 2. Initial Subsetting ####
#subetting data - need to center age before adding interaction terms..
raw_subsets = list("female" = 
                     cohort[gender_code == 0,][,gender_code := NULL],
                   "male" = 
                     cohort[gender_code == 1,][,gender_code := NULL])

#center age and store results
base_age = c()
for (raw_subset in raw_subsets){
  base_age = c(base_age,raw_subset[,mean(nhi_age)])
  raw_subset[,c_nhi_age := (nhi_age - mean(nhi_age))]
  raw_subset[,nhi_age := NULL]
}
names(base_age) = names(raw_subsets)

# initial cohort cleaning
subsets = list()
for (i in 1:length(raw_subsets)){
  dt = copy(raw_subsets[[i]])
  
  #add interaction terms
  dt = prespec_ints(dt)
  
  #transform test results
  dt = transform_labs(dt,
                      c_egfr = 1,
                      c_hba1c = 1,
                      c_tchdl = 1,
                      c_tri = 1)
  
  # discretize egfr - DONE AFTER MARTINGALE PLOTS
  #dt[,egfr_cat := egfr_groups(egfr)]
  
  #setnames(dt,"c_egfr","eGFR")
  setnames(dt,"c_hba1c","HbA1c")
  setnames(dt,"c_tchdl","TC/HDL ratio")
  setnames(dt,"c_tri","TG")
  
  #list, name append
  dt_list = list(dt)
  names(dt_list) = names(raw_subsets)[i]
  subsets = append(subsets,dt_list)
}

#* ANALYSIS AND VALIDATION ######


# # ?. Influential observations (not performed) ####
# #fit initial models
# inf_models = list()
# for (i in 1:length(subsets)){
#   #copy or safety
#   dt = copy(subsets[[i]])[,-c("VSIMPLE_INDEX_MASTER")]
# 
#   #fit model
#   inf_model = coxph(Surv(total_fu,indicator) ~ .,
#                 data = dt,
#                 ties = "breslow",
#                 model = TRUE)
# 
#   #list, name, append
#   list_inf_model = list(inf_model)
#   names(list_inf_model) = names(subsets)[i]
#   inf_models = append(inf_models,list_inf_model)
# }
# 
# # dfbetas residuals from survfit
# resid = residuals(models[[1]],
#                   type = "dfbetas",
#                   weighted = TRUE)
# resid_sum = rowSums(resid)
# max(resid_sum)
# hist(resid_sum)
# inf_index = which(resid_sum > 0.05)
# inf_observations = subsets[[1]][inf_index,]
# hist(inf_observations$total_fu)
# hist(inf_observations$c_hba1c)
# 
# hba1c_index = which(subsets[[1]]$c_hba1c > 10)
# hba1c_cohort = subsets[[1]][hba1c_index,]
# hba1c_cohort
# hist(resid_sum[hba1c_index])
# 
# #events most influential

#* Testing linearity ####
#functions
martingale_plot = function(dt,
                          var,
                          var_name = NULL,
                          type = c("full",
                                   "excluded",
                                   "null"),
                          display_resid = 1){
  #generates a crude martingale residual plot

  #var = predictor of interest
  #var_name = name of predictor to use when plotting
  #type = what model's residual to use -
  #           Full = full model
  #           Excluded = Full model minus predictor
  #           Null = null model
  #display_resid = plot includes actual residuals
  
  #if name blank simply use raw predictor
  if (is.null(var_name)){
    var_name = var
  }
  
  

  if (type == "full"){
    form = as.formula("Surv(total_fu,indicator)~.")
  } else if (type == "excluded") {
    form = as.formula(
      paste(
        "Surv(total_fu,indicator) ~ . -",var,sep=" "))
  } else {
    form = as.formula(
      "Surv(total_fu,indicator) ~ 1"
    )
  }
  
  sample_model = coxph(form,
                       ties = "breslow",
                       data = dt[,-c("VSIMPLE_INDEX_MASTER")])
  
  martingale_resid = residuals(sample_model,
                               type = "martingale")
  
  deviance_ = residuals(sample_model,
                        type = "deviance")
  x = unlist(dt[,c(var),with = FALSE])
  ss = smooth.spline(x,martingale_resid,df = 6)
  #ps = pspline(x,martingale_resid,0.5,df = 4)
  #ls = lowess(x,martingale_resid,f = 0.5)
  #plot(x,martingale_resid,lty = 2)
  #lines(ls,col = "red", lty = 2)
  #plot(ls)

  if (display_resid == 1){
    plot = ggplot()+
      geom_point(aes(x = x,
                     y = martingale_resid),
                 color = "lightgrey") +
      geom_line(aes(x = ss$x,
                    y = ss$y),
                color = "red",
                lwd = 2) +
      #labs(title ="Martingale residual plot") +
      xlab(paste("transformed",var_name)) +
      ylab(paste(type,"model residuals")) +
      theme_bw()
  } else {
    plot = ggplot() +
      geom_line(aes(x = ss$x,
                  y = ss$y)) +
      #labs(title ="Martingale residual plot") +
      xlab(paste("transformed",var_name)) +
      ylab(paste(type,"model residuals")) +
      theme_bw()
  }
  return(plot)
}


generate_martingales = function(subset,
                                vars,
                                save = 0){
  
  #extract subset name and subset
  subset_name = names(subset)[1]
  subset = subset[[1]]
  
  
  for (i in 1:length(vars)){
    
    #extract test name
    var_name = names(vars)[i]
    var = vars[[i]]
    
    #initialize plotlist
    plotlist = list()
    for (display in c(1,0)){
      for (type in c("full","excluded","null")){
        plot = martingale_plot(subset,
                               var,
                               var_name = var_name,
                               type = type,
                               display_resid = display)
        plotlist = append(plotlist,list(plot))
      }
    }
    
    fullplot = ggarrange(plotlist = plotlist,
                         ncol = 3,
                         nrow = 2)
    fullplot = annotate_figure(fullplot,
                               top =
                                 text_grob(paste("Martingale residuals of ",
                                                 var_name,
                                                 " in ",
                                                 subset_name,
                                                 "s",
                                                 sep = ""),
                                           face = "bold",
                                           size = 14))
    
    path = "plots/diagnostic_plots/martingale_plots"
    name = paste(var,subset_name,"martingale_residuals_h.png",sep = "_")
    
    if (save == 1){
      ggsave(name,
             fullplot,
             path = path,
             width = 10,
             height = 6)
    } else {
      fullplot
    }
    
  }
}

# Actually running the martingale analysis ####
set.seed(1)
dt_sample = copy(subsets[[1]][
  sample(1:dim(subsets[[i]])[1],1000)])
dt_sample = copy(subsets[[1]])
plot = martingale_plot(subsets[[2]],
            "c_tri",
           type = "full",
           display_resid = 1)
plot

#generating martingale residual plots
for (i in 1:2){
  #copy for safety
  dt = copy(subsets[i])
  
  #use sample for debugging
  #dt_sample = list("temp" = dt[sample(1:100000,10000)])
  generate_martingales(dt,
                       list("age" = "c_nhi_age",
                         "eGFR" = "c_egfr",
                         "HbA1c" = "c_hba1c",
                         "TC/HDLC ratio" = "c_tchdl",
                         "TG" = "c_tri"),
                       save = 1) #alter this to save
}



# 2. transforming c_egfr and c_tri and fitting model ####
for (i in 1:length(subsets)){
  dt = copy(raw_subsets[[i]])
  dt = dt[,c("VSIMPLE_INDEX_MASTER","egfr")]
  dt[,":="(egfr_cat = egfr_groups(egfr),
           egfr = NULL)]
           # log_tri = log2(tri),
           # tri = NULL)]
  dt = cut_factors(dt)

  #append to subsets
  
  subsets[[i]] = subsets[[i]][dt]
  subsets[[i]][,":="(c_egfr = NULL)]
}
class(subsets[[1]])
names(subsets[[1]])

# 
# # OLD JUNK CODE
# # 2. transforming c_egfr and c_tri and fitting model ####
# for (i in 1:length(subsets)){
#   dt = copy(raw_subsets[[i]])
#   dt = dt[,c("VSIMPLE_INDEX_MASTER","egfr","tri")]
#   dt[,":="(egfr_cat = egfr_groups(egfr),
#            egfr = NULL)]
#   # log_tri = log2(tri),
#   # tri = NULL)]
#   dt = cut_factors(dt)
#   
#   #append to subsets
#   
#   subsets[[i]](subsets[[i]],dt)
#   subsets[[i]][,":="(c_egfr = NULL)]
# }
# # 
# # #checking linearity of log transform
# # for (i in 1:2){
# #   #copy for safety
# #   dt = copy(subsets[i])
# # 
# #   #use sample for debugging
# #   #dt_sample = list("temp" = dt[sample(1:100000,10000)])
# #   generate_martingales(dt,
# #                        list("triglycerides" = "log_tri"),
# #                        save = 1) #alter this to save
# # }

#fitting models
models = list()
for (i in 1:length(subsets)){
  #copy or safety
  dt = copy(subsets[[i]])[,-c("VSIMPLE_INDEX_MASTER")]

  #fit model
  model = coxph(Surv(total_fu,indicator) ~ .,
                data = dt,
                ties = "breslow",
                model = TRUE)

  #list, name, append
  list_model = list(model)
  names(list_model) = names(subsets)[i]
  models = append(models,list_model)
}



# 3. Testing PH assumption ####
for (i in 1:2){
  #run cox.zph
  ph_test = cox.zph(models[[i]])
  model_name = names(models)[i]
  
  for (j in 1:4){
    ind = dim(ph_test$table)[1]-j
    
    #set up device
    var_name = row.names(ph_test$table)[ind]
    path = paste("plots/diagnostic_plots/schoenfeld_plots/",
                 gsub("/","",gsub("`","",gsub(" ","",var_name))),
                 "_",
                 model_name,
                 ".png",
                 sep = "")
    png(path)
    plot(ph_test[ind],
         resid = FALSE)
    dev.off()
  } 
}

#constructing global test table

ph_test_f = cox.zph(models[[1]])
ph_test_m = cox.zph(models[[2]])
f = data.frame(ph_test_f$table)[,c(3),drop = FALSE]
m = data.frame(ph_test_m$table)[,c(3),drop = FALSE]
ph_test_table = cbind(f,m)
ph_test_table$rn = rownames(ph_test_table)
names(ph_test_table) = c("Women", "Men","rn")

names(as_tibble(ph_test_table))
names(as_tibble(clean_predictors))
ph_test_table
clean_predictors
ph_table = as_tibble(ph_test_table) %>%
  left_join(as_tibble(clean_predictors), by = "rn")

ph_table[c(4,1,2)]
ph_table$Women = round(ph_table$Women,digits = 3)
ph_table$Men = round(ph_table$Men, digits = 3)
output = ph_table[,c(4,1,2)]
output
write.csv(output,
          "tables/ph_table.csv")


# 4. Interactions ####


# 4a. MFPI #####
#exporting cohorts for MFPI analysis
mfpi_path = "//uoa.auckland.ac.nz/Shared/MED/EPBI/ViewData/Users/bbat644/Desktop/code/MFPI/"
#full_cohort
names(subsets)
for (name in names(subsets)){
  dt = copy(subsets[[name]])
  # remove interactions
  dt[,c("age_bp","age_af","bp_diabetes","at_diabetes","age_diabetes","bp_ll"):= NULL]
  dt
  path = paste(mfpi_path,name,".csv",sep = "")
  write.csv(dt,path,row.names = FALSE)
  set.seed(4)
  sample = dt[sample(dim(dt)[1],100000),]
  path = paste(mfpi_path,name,"sample.csv",sep = "")
  write.csv(sample,path,row.names = FALSE)
}

# 4b. MFPI cont
#fitting model specified by MFPI and seeing which effects are significant
potential_ints = function(dt){
  #super simple justt adds the interactions identified by MFPI 
  #data table and returns the data table
  
  #copy for safety
  dt = copy(dt)
  
  #add old interactions
  dt[,":="(at_age = ph_antithrombotic_prior_6mths*c_nhi_age,
           at_hba1c = ph_antithrombotic_prior_6mths*c_hba1c,
           at_tri = ph_antithrombotic_prior_6mths *log_tri,
           age_hba1c =  c_nhi_age*c_hba1c,
           age_tchdl = c_nhi_age*c_tchdl,
           age_tri = c_nhi_age*log_tri)]
  
  #returns the new data frame
  return(dt)
  
}

MFPI_test_cohorts = list()
MFPI_test_models = list()
for (i in 1:2){
  
  MFPI_test_cohort = potential_ints(subsets[[i]])
  
  #equivalency check
  # MFPI_test_cohort[,c("log_tri","age_tri") := NULL]
  MFPI_test_model = coxph(Surv(total_fu,indicator)~.,
                          data = MFPI_test_cohort[,-1],
                          ties = "breslow",
                          model = TRUE)
  
  MFPI_test_cohorts = append(MFPI_test_cohorts,list(MFPI_test_cohort))
  MFPI_test_models = append(MFPI_test_models,list(MFPI_test_model))
}

#display model and check
MFPI_test_models[[1]]

#SELECT INTERACTION AGE_TCHDL!!!!!




# 4c. GLMNET ####
# # not useful
# dt = copy(subsets[[1]])
# #running on a 10,000 people subset
# set.seed(1)
# ind = sample(dim(dt)[1],10000)
# dt = dt[ind,]
# #specify covariates
# potential_ints = c("hx_vdr_diabetes",
#                    "hx_af",
#                    "ph_bp_lowering_prior_6mths",
#                    "ph_lipid_lowering_prior_6mths",
#                    "ph_antithrombotic_prior_6mths",
#                    "c_nhi_age",
#                    "c_hba1c",
#                    "c_tchdl",
#                    "egfr_cat_G3a",
#                    "egfr_cat_G3b",
#                    "log_tri")
# 
# x_excluded = dt[,-c("VSIMPLE_INDEX_MASTER",
#                     "indicator",
#                     "total_fu",
#                     potential_ints),with= FALSE]
# x = dt[,c(potential_ints),with=FALSE] #what columns to examine for interactions
# 
# #create all possible interactions pf covariates
# int_matrix = data.table() #datatable to store interaction terms
# for (i in (1:(dim(x)[2]-1))){
#   
#   #candidate predictor
#   z2 = as.matrix(x[,i,with = FALSE])[,1]
#   
#   #rest of the matrix
#   x_ = x[,c((i+1):(dim(x)[2])),with = FALSE]
#   
#   #multiply and append
#   int_x = z2*x_
#   colnames(int_x) = unlist(lapply(
#     colnames(int_x), \(input) paste(colnames(x)[i],input,sep = "*")))
#   int_matrix = cbind(int_matrix,int_x)
# }
# names(int_matrix)
# full_x = data.matrix(cbind(x_excluded,x,int_matrix))
# names(full_x)
# 
# #running glmnet
# y = Surv(dt$total_fu,dt$indicator)
# 
# fit = glmnet(full_x,y, family = "cox")
# plot(fit,xvar = "lambda",label = TRUE)
# coef(fit, s = 0.000) #plotting values at a certain level of lambda
# set.seed(1)
# 
# # by C and deviance
# cvfit_C = cv.glmnet(full_x,y,family = "cox",type.measure = "C")
# cvfit_deviance = cv.glmnet(full_x,y,family = "cox",type.measure = "deviance")
# 
# plot(cvfit_C)
# cvfit_C$lambda.1se
# coef(fit, s = cvfit_C$lambda.1se)
# 
# plot(cvfit_deviance)
# cvfit_deviance$lambda.1se
# coef(fit, s = cvfit_deviance$lambda.1se)
# 
# 
# 
# 
