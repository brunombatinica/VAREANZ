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
# 4a. New cohort
subsets = list()
for (i in 1:length(raw_subsets)){
  
  #copy for safety
  dt = copy(raw_subsets[[i]])
  
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

# for computational reasons - while testing
seed = 123
## ---- 2a  Set up imputation spec ----------------------------------------
dt = subsets[[1]]
exclude_vars <- c("VSIMPLE_INDEX_MASTER", "total_fu", "indicator")

pred   <- mice::quickpred(dt, exclude = exclude_vars)
meth   <- mice::make.method(dt)
meth[exclude_vars] <- ""            # "" tells mice to *skip* imputation

## ---- 2b  Run mice -------------------------------------------------------
imp <- mice(
  data             = dt,
  m = 5,              
  method           = meth,          # default is "pmm" for numeric, "polyreg"/"logreg" for factors
  predictorMatrix  = pred,
  maxit            = 5,            # iterations per chain (default 5 is often enough)
  seed             = seed,
  printFlag        = FALSE
  )
fitlist <- with(imp, {
  coxph(Surv(total_fu, indicator) ~ eth_cat_Chinese + eth_cat_Indian + 
    eth_cat_NZM + eth_cat_Other + eth_cat_Pacific + c_en_nzdep_q +
    hx_vdr_diabetes + hx_af + ph_bp_lowering_prior_6mths + ph_lipid_lowering_prior_6mths +
    ph_antithrombotic_prior_6mths + c_nhi_age + age_bp + age_diabetes +
    age_af + bp_diabetes + at_diabetes + bp_ll + c_hba1c + c_tri +
    c_tchdl + egfr_cat_G3a + egfr_cat_G3b + age_tchdl - VSIMPLE_INDEX_MASTER,
        ties = "breslow", model = TRUE)
})
summary(pool(fitlist))

##  keep only rows with no missing values in ANY model variable
orig_data = copy(dt)
vars_used <- c(
  "total_fu", "indicator",
  "eth_cat_Chinese", "eth_cat_Indian", "eth_cat_NZM", "eth_cat_Other",
  "eth_cat_Pacific", "c_en_nzdep_q", "hx_vdr_diabetes", "hx_af",
  "ph_bp_lowering_prior_6mths", "ph_lipid_lowering_prior_6mths",
  "ph_antithrombotic_prior_6mths", "c_nhi_age", "age_bp", "age_diabetes",
  "age_af", "bp_diabetes", "at_diabetes", "bp_ll", "c_hba1c", "c_tri",
  "c_tchdl", "egfr_cat_G3a", "egfr_cat_G3b", "age_tchdl"
)
cc_data <- orig_data %>% filter(complete.cases(select(., all_of(vars_used))))
cat("N imputations  :", imp$m, "\n")
cat("N complete‑case:", nrow(cc_data), "\n")

## ---------------------------------------------------------------
## 2.  Fit the complete‑case Cox model  --------------------------
fit_cc <- coxph(
  Surv(total_fu, indicator) ~
    eth_cat_Chinese + eth_cat_Indian + eth_cat_NZM + eth_cat_Other +
    eth_cat_Pacific + c_en_nzdep_q + hx_vdr_diabetes + hx_af +
    ph_bp_lowering_prior_6mths + ph_lipid_lowering_prior_6mths +
    ph_antithrombotic_prior_6mths + c_nhi_age + age_bp + age_diabetes +
    age_af + bp_diabetes + at_diabetes + bp_ll + c_hba1c + c_tri +
    c_tchdl + egfr_cat_G3a + egfr_cat_G3b + age_tchdl,
  data  = cc_data,
  ties  = "breslow",
  model = TRUE
)
fit_cc



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

