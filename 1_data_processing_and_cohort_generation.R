# Clearing workspace
rm(list = ls())
dev.off()

# Loading in necessary libraries and data ################################
source("0_helper_functions.R")

#set working directory
setwd("//uoa.auckland.ac.nz/Shared/MED/EPBI/ViewData/Users/bbat644/Desktop/code")

# 1. reading in raw data ####
#test results
hba1c_raw = read.fst("data/raw_data/testsafe/HBA1C.fst")
scr_raw = read.fst("data/raw_data/testsafe/SCR.fst")
tchdl_raw = read.fst("data/raw_data/testsafe/TCHDL.fst")
tri_raw = read.fst("data/raw_data/testsafe/TRI.fst")

#extract only results from people in Auckland and Northland

# 2. Creating "initial AN_cohort"
#read in full HC/VARIANZ data
health_contact_raw = read.fst(
  "data/raw_data/varianz2012/VARIANZ_2012.fst",
  as.data.table = TRUE)

#raw data summary
names(health_contact_raw)
head(health_contact_raw)
dim(health_contact_raw)

# 1. initial exclusions #####
hc_dt_raw = data.table(health_contact_raw)

#Remove people with no follow-up 
hc_dt = hc_dt_raw[!(is.na(end_fu_date)),]
dim(hc_dt)

#3x antianginals to make denominators match
hc_dt = hc_dt[!(ph_antianginals_prior_5yrs_3evts == 1)]
dim(hc_dt)


# 2. Cleaning up predictors/outcomes ####
# Refining medication variables ####
#switching bp and ll to integer
num_ind = c("ph_bp_lowering_prior_6mths",
            "ph_lipid_lowering_prior_6mths")
hc_dt[,(num_ind) := lapply(.SD,as.integer),.SDcols = num_ind]


#Combining antiplatelet and anticoagulent
###THERE IS SLIGHT DISCREPENCY BETWEEN THIS DATA 
###THE SUMMARY STATISTICS IN THE PAPER
hc_dt[,ph_antithrombotic_prior_6mths := as.integer(
  or(ph_antiplatelets_prior_6mths,
     ph_anticoagulants_prior_6mths))]

#table(health_contact$en_prtsd_eth)

# ethnicity ####
#Defining ethnicity groups
#ethnicity list
eth = list( "NZE" = 1,
            "NZM" = 2,
            "Pacific" = 3,
            "Indian" = 43,
            "Chinese" = 42,
            "Other" = c(4,5,9)) #other asian,MELAA,other
#creating a new column of ethnicity categories
eth_cat = hc_dt[,en_prtsd_eth]
for (i in names(eth)){
  eth_cat = replace(
    eth_cat,which(eth_cat %in% eth[[i]]),i)  
}
eth_cat = factor(eth_cat)
#Setting NZE to be top level factor
eth_cat = relevel(eth_cat, ref = "NZE")
hc_dt[,eth_cat := eth_cat]

# dhb ####
#add an indicator if in AUckland or northalnd
Auckland_Northland = c("Auckland","Counties Manukau","Northland","Waitemata")
hc_dt[,auckland_northland := as.integer(DHB_name %in% Auckland_Northland)]
#sum(hc_dt$auckland_northland)

# centering deprivation ####
hc_dt[,c_en_nzdep_q := en_nzdep_q - 3]

# Creating outcome variables ####

# indicator
hc_dt[,indicator := as.integer(out_broad_cvd | imp_fatal_cvd)]

# #Non fatal cvd
# health_contact["non_fatal_cvd"] = as.numeric(
#   health_contact$out_broad_cvd &
#     !health_contact$imp_fatal_cvd
# )

# total follow-up
#Varianz_lastcontact does not censor on event
index = as.Date("2012-12-31", "%Y-%m-%d")
hc_dt[,total_fu := as.numeric(
  difftime(
    hc_dt[,end_fu_date],index,"days"),
  units = "days") / 365
]


# 3. clean dataframe
model_vars = c("VSIMPLE_INDEX_MASTER",
               "indicator",
               "total_fu",
               "gender_code",
               "DHB_name",
               "nhi_age",
               "eth_cat",
               "c_en_nzdep_q",
               "hx_vdr_diabetes",
               "hx_af",
               "ph_bp_lowering_prior_6mths",
               "ph_lipid_lowering_prior_6mths",
               "ph_antithrombotic_prior_6mths")

#names(hc_dt)
hc_cohort_raw = hc_dt[,.SD,.SDcols = model_vars]
an_cohort_raw = hc_dt[auckland_northland == 1,
                      .SD,.SDcols = model_vars]
dim(hc_cohort_raw)
dim(an_cohort_raw)


# 4. demographics table 
source("0_helper_functions.R")


# extract subset list
subset_list = 
  list("Auckland/Northland" = 
         list(hc_cohort_raw[DHB_name %in% Auckland_Northland,VSIMPLE_INDEX_MASTER]),
       "Rest of New Zealand" = 
         list(hc_cohort_raw[!(DHB_name %in% Auckland_Northland),VSIMPLE_INDEX_MASTER]),
       "Total New Zealand" = 
         list(hc_cohort_raw[,VSIMPLE_INDEX_MASTER]))

AN_demo_table = demo_table_generator(subset_list)
AN_demo_table
write.csv(AN_demo_table,
          "tables/AN_demo_table.csv")





# 5. testsafe basic analysis

#initial exploration
str(scr_raw)
str(hba1c_raw)
str(tchdl_raw)
str(tri_raw)

#which dates was teh cohort indexed on
min(hba1c_raw$RESULT_DATE)
max(hba1c_raw$RESULT_DATE)


#only look at results from individuals in auckland and northland
an_index = an_cohort_raw$VSIMPLE_INDEX_MASTER

an_hba1c_raw = hba1c_raw %>%
  filter(VSIMPLE_INDEX_MASTER %in% an_index)
an_scr_raw = scr_raw %>% 
  filter(VSIMPLE_INDEX_MASTER %in% an_index)
an_tchdl_raw = tchdl_raw %>% 
  filter(VSIMPLE_INDEX_MASTER %in% an_index)
an_tri_raw = tri_raw %>%
  filter(VSIMPLE_INDEX_MASTER %in% an_index)

str(an_scr_raw)
str(an_hba1c_raw)
str(an_tchdl_raw)
str(an_tri_raw)


dim(an_scr_raw)[1] + dim(an_hba1c_raw)[1] + dim(an_tchdl_raw)[1] + dim(an_tri_raw)[1]

# 6. filter out excluded results ###########
fac_raw = c(an_scr_raw$FACILITY_CODE,an_hba1c_raw$FACILITY_CODE,an_tchdl_raw$FACILITY_CODE,an_tri_raw$FACILITY_CODE)
table(fac_raw)

#Remove out of auckland labs
labs = c("APHDELPHIC",
         "DIAGNOSTIC",
         "inteLAB",
         "LTA",
         "NDHBDELPHIC",
         "NPL",
         "SALAB",
         "WDHB")

sum(!(an_scr_raw$FACILITY_CODE %in% labs))
sum(!(an_hba1c_raw$FACILITY_CODE %in% labs))
sum(!(an_tchdl_raw$FACILITY_CODE %in% labs))
sum(!(an_tri_raw$FACILITY_CODE %in% labs))

#further filter SCR RESULTS
excluded_profile = c("CATHETER URINE",
                     "CLIN. CHEMISTRY URINE",
                     "CLIN.CHEMISTRY MISC FLUID",
                     "TEST NOT FOUND")

table(an_scr_raw$PROFILE_UPPER)
sum((an_scr_raw$PROFILE_UPPER %in% excluded_profile))

#filter out results from out of Auckland labs and from scr_ from excluded profile
hba1c_ = an_hba1c_raw %>%
  filter(FACILITY_CODE %in% labs)
scr_ = an_scr_raw %>% 
  filter(FACILITY_CODE %in% labs &
           !(PROFILE_UPPER %in% excluded_profile))
tchdl_ = an_tchdl_raw %>% 
  filter(FACILITY_CODE %in% labs)
tri_ = an_tri_raw %>%
  filter(FACILITY_CODE %in% labs)

#post exclusions
str(scr_)
str(hba1c_)
str(tchdl_)
str(tri_)

dim(scr_)[1] + dim(hba1c_)[1] + dim(tchdl_)[1] + dim(tri_)[1]

sum(an_cohort_raw$VSIMPLE_INDEX_MASTER %in% scr_$VSIMPLE_INDEX_MASTER) 
mean(an_cohort_raw$VSIMPLE_INDEX_MASTER %in% scr_$VSIMPLE_INDEX_MASTER)

sum(an_cohort_raw$VSIMPLE_INDEX_MASTER %in% hba1c_$VSIMPLE_INDEX_MASTER) 
mean(an_cohort_raw$VSIMPLE_INDEX_MASTER %in% hba1c_$VSIMPLE_INDEX_MASTER)

sum(an_cohort_raw$VSIMPLE_INDEX_MASTER %in% tchdl_$VSIMPLE_INDEX_MASTER) 
mean(an_cohort_raw$VSIMPLE_INDEX_MASTER %in% tchdl_$VSIMPLE_INDEX_MASTER)

sum(an_cohort_raw$VSIMPLE_INDEX_MASTER %in% tri_$VSIMPLE_INDEX_MASTER) 
mean(an_cohort_raw$VSIMPLE_INDEX_MASTER %in% tri_$VSIMPLE_INDEX_MASTER)

#add result variable
hba1c = hba1c_ %>%
  mutate(test = "hba1c",result = EN_OBSR_RESULT_NUM) %>%
  select(VSIMPLE_INDEX_MASTER,FACILITY_CODE,RESULT_DATE,test,result) 
egfr = scr_ %>% 
  mutate(test = "egfr", result = EN_OBSR_RESULT_EGFR) %>%
  select(VSIMPLE_INDEX_MASTER,FACILITY_CODE,RESULT_DATE,test,result)
# scr = scr_ %>% 
#   mutate(test = "scr", result = EN_OBSR_RESULT_NUM) %>%
#   select(VSIMPLE_INDEX_MASTER,FACILITY_CODE,RESULT_DATE,test,result) 
tchdl = tchdl_ %>%
  mutate(test = "tchdl",result = EN_OBSR_RESULT_NUM) %>%
  select(VSIMPLE_INDEX_MASTER,FACILITY_CODE,RESULT_DATE,test,result) 
tri = tri_ %>%
  mutate(test = "tri",result = EN_OBSR_RESULT_NUM) %>%
  select(VSIMPLE_INDEX_MASTER,FACILITY_CODE,RESULT_DATE,test,result)

str(egfr)
str(hba1c)
str(tchdl)
str(tri)
# #remove raw imports
# rm(list = c("hba1c_","hba1c_raw",
#             "scr_","scr_raw",
#             "tchdl_","tchdl_raw",
#             "tri_","tri_raw"))




#bind results together
all_tests = rbind(hba1c,egfr,tchdl,tri)



##############
# Chronologically the analysis now goes to 2_testsafe_analysis
#############





# 7. MoH exclusions ####
#eGFR based exclusions
egfr_exclusions = egfr %>%
  filter(and(RESULT_DATE >= as.Date("2007-12-31"),
             RESULT_DATE < as.Date("2015-01-01"))) %>%
  filter(result < 30.0) %>%
  select(VSIMPLE_INDEX_MASTER,RESULT_DATE,result)
length(unique(egfr_exclusions$VSIMPLE_INDEX_MASTER))




tchdl_exclusions = tchdl %>%
  filter(and(RESULT_DATE >= as.Date("2007-12-31"),
             RESULT_DATE < as.Date("2015-01-01"))) %>%
  filter(result > 8.0) %>%
  select(VSIMPLE_INDEX_MASTER,RESULT_DATE,result)
length(unique(tchdl_exclusions$VSIMPLE_INDEX_MASTER))

tri_exclusions = tri %>%
  filter(and(RESULT_DATE >= as.Date("2007-12-31"),
             RESULT_DATE < as.Date("2015-01-01"))) %>%
  filter(result > 11.0) %>%
  select(VSIMPLE_INDEX_MASTER,RESULT_DATE,result)
length(unique(tri_exclusions$VSIMPLE_INDEX_MASTER))




# export results
write.fst(hba1c,
          path = "data/testsafe/hba1c.fst")
write.fst(egfr,
          path = "data/testsafe/egfr.fst")
# write.fst(scr,
#           path = "data/testsafe/scr.fst")
write.fst(tchdl,
          path = "data/testsafe/tcdl.fst")
write.fst(tri,
          path = "data/testsafe/tri.fst")
write.fst(all_tests,
          path = "data/testsafe/all_tests.fst")

write.fst(egfr_exclusions,
          path = "data/testsafe/egfr_exclusions.fst")
write.fst(tchdl_exclusions,
          path = "data/testsafe/tchdl_exclusions.fst")
write.fst(tri_exclusions,
          path = "data/testsafe/tri_exclusions.fst")





# 4. Testsafe exclusions ONLY DONE ON AN COHORT

an_cohort_raw = as.data.table(an_cohort_raw)
egfr_exc_cohort = an_cohort_raw[
  VSIMPLE_INDEX_MASTER %in% 
    egfr_exclusions$VSIMPLE_INDEX_MASTER]
dim(egfr_exc_cohort)
sum(egfr_exc_cohort$indicator)
mean(egfr_exc_cohort$indicator)

tchdl_exc_cohort = an_cohort_raw[
  VSIMPLE_INDEX_MASTER %in% 
    tchdl_exclusions$VSIMPLE_INDEX_MASTER]
dim(tchdl_exc_cohort)                     
sum(tchdl_exc_cohort$indicator) 
mean(tchdl_exc_cohort$indicator) 

tri_exc_cohort = an_cohort_raw[
  VSIMPLE_INDEX_MASTER %in% 
    tri_exclusions$VSIMPLE_INDEX_MASTER]
dim(tri_exc_cohort)                     
sum(tri_exc_cohort$indicator) 
mean(tri_exc_cohort$indicator) 

an_cohort = an_cohort_raw[!(VSIMPLE_INDEX_MASTER %in% 
                              c(egfr_exclusions$VSIMPLE_INDEX_MASTER,
                                tchdl_exclusions$VSIMPLE_INDEX_MASTER,
                                tri_exclusions$VSIMPLE_INDEX_MASTER))]

dim(an_cohort)

# 5. exporting cox data######
# saving cohorts ####
write.fst(hc_cohort_raw,
          path = "data/cohorts/hc_cohort_raw.fst")
write.fst(an_cohort_raw,
          path = "data/cohorts/an_cohort_raw.fst")
write.fst(an_cohort,
          path = "data/cohorts/an_cohort.fst")






