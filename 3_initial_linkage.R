#Adding test safe variables to the auckland_northland cohort

# Clearing workspace
rm(list = ls())
dev.off()

# Loading in necessary libraries and data ################################
source("0_helper_functions.R")

#set working directory
setwd("C:/Users/bruno/OneDrive/Documents/Code/projects/VAREANZ/Honours")

# 1. Load in data sets ####
# cohort
an_cohort = read.fst(
  path = "data/cohorts/an_cohort.fst",
  as.data.table = TRUE
)

# testsafe results
alltests = read.fst(
  path = "data/testsafe/clean_tests.fst", ####USE clean_tests here for uptodate exclusions
  as.data.table = TRUE)

#clear facility code
alltests[,c("FACILITY_CODE","FACILITY_NAME") := NULL]
table(alltests$test)

#set keys to vsimple_index
setkey(an_cohort,VSIMPLE_INDEX_MASTER)
setkey(alltests,VSIMPLE_INDEX_MASTER)

# 2. Selecting results from ... 
#start = 01-01-2008 to 
#end = 01-01-2013 (Index date) 
# and cleaning ####
start_interval = "2008-01-01"
end_interval = "2013-01-01"
alltests_interval = alltests[and(
  RESULT_DATE < as.Date(start_interval),
  RESULT_DATE >= as.Date(end_interval)),][
    VSIMPLE_INDEX_MASTER %in% an_cohort$VSIMPLE_INDEX_MASTER
  ]


#select result closest to 2013-01-01
#egfr
egfr = alltests_interval[
  test == "eGFR"][
    order(abs(RESULT_DATE - as.Date("2013-01-01"))),
    head(.SD,1), 
    keyby = VSIMPLE_INDEX_MASTER][
      ,.(VSIMPLE_INDEX_MASTER,RESULT_DATE,egfr = result)]

#hba1c
hba1c = alltests_interval[
  test == "HbA1c"][
    order(abs(RESULT_DATE - as.Date("2013-01-01"))),
    head(.SD,1), 
    keyby = VSIMPLE_INDEX_MASTER][
      ,.(VSIMPLE_INDEX_MASTER,RESULT_DATE,hba1c = result)]

#tchdl
tchdl = alltests_interval[
  test =="TCHDL"][
    order(abs(RESULT_DATE - as.Date("2013-01-01"))),
    head(.SD,1),
    keyby = VSIMPLE_INDEX_MASTER][
      ,.(VSIMPLE_INDEX_MASTER,RESULT_DATE,tchdl = result)]

#trig
tri = alltests_interval[
  test == "TG"][
    order(abs(RESULT_DATE - as.Date("2013-01-01"))),
    head(.SD,1), 
    keyby = VSIMPLE_INDEX_MASTER][
      ,.(VSIMPLE_INDEX_MASTER,RESULT_DATE,tri = result)]


# 3. joining and exporting ####
#NA joining the 3 cohorts
# remove result date and join
an_cohort_ts = merge(an_cohort,egfr[,-c("RESULT_DATE")], all.x = TRUE) |> 
  {\(x) merge(x,hba1c[,-c("RESULT_DATE")], all.x = TRUE)}() |> 
  {\(x) merge(x,tchdl[,-c("RESULT_DATE")], all.x = TRUE)}() |>
  {\(x) merge(x,tri[,-c("RESULT_DATE")], all.x = TRUE)}()

sum(complete.cases(an_cohort_ts))
dim(an_cohort_ts[is.na(tchdl)][!is.na(tri)])

#exporting datasets 
write.fst(an_cohort_ts,
          paste0("data/cohorts/an_cohort_ts_",start_interval,"_",end_interval,".fst"))


# ################## FURTHER INVESTIGATIONS #####################
# # 4. number of results in cohort ####
# #summary
# table(alltests_0815$test)
# min(alltests_0815$RESULT_DATE)
# max(alltests_0815$RESULT_DATE)

# sum(!is.na(an_cohort_ts$egfr))
# mean(!is.na(an_cohort_ts$egfr))

# sum(!is.na(an_cohort_ts$hba1c))
# mean(!is.na(an_cohort_ts$hba1c))

# sum(!is.na(an_cohort_ts$tchdl))
# mean(!is.na(an_cohort_ts$tchdl))

# sum(!is.na(an_cohort_ts$tri))
# mean(!is.na(an_cohort_ts$tri))


# # 5. an cohort baseline demogrpahics ####
# source("0_helper_functions.R")
# subset_list = list("Women" = list(an_cohort[gender_code == 0,VSIMPLE_INDEX_MASTER]),
#                    "Men" = list(an_cohort[gender_code == 1,VSIMPLE_INDEX_MASTER]))
# study_cohort_demo_table = demo_table_generator(subset_list,
#                                                test_rows = TRUE)
# print.data.frame(study_cohort_demo_table)
# write.csv(study_cohort_demo_table,
#           "tables/study_cohort_demo_table.csv")

# #hist(scr$RESULT_DATE, breaks = 30)


# # 5a. CC demographics cs AN demographics ####
# subset_list = list(
#   "Complete cases" = 
#     list(an_cohort_ts[complete.cases(an_cohort_ts),VSIMPLE_INDEX_MASTER]),
#   "Incomplete cases" = 
#     list(an_cohort_ts[!complete.cases(an_cohort_ts),VSIMPLE_INDEX_MASTER]),
#   "Study cohort" = 
#     list(an_cohort_ts[,VSIMPLE_INDEX_MASTER]))
  
# missing_data_cohort_demo_table = demo_table_generator(subset_list)
# missing_data_cohort_demo_table
# write.csv(missing_data_cohort_demo_table,
#           "tables/missing_data_cohort_demo_table.csv")

# # portion of missingness by age

# temp = as_tibble(an_cohort_ts)
# temp
# temp = temp %>%
#   mutate(cc = as.integer(!is.na(egfr)&!is.na(hba1c)&!is.na(tchdl)&!is.na(tri)))

# temp %>% 
#   group_by(nhi_age) %>%
#   summarize(pcc = mean(cc)) %>%
#   ggplot() +
#   geom_bar(mapping = aes(x = nhi_age,y = pcc),stat = "identity", width = 1)

# #complete case by age
# cc_age = temp %>%
#   ggplot() +
#   geom_bar(mapping = aes(x = nhi_age, fill = factor(cc)), 
#            position = "fill",
#            width = 1) +
#   labs(#title = "Proportion of complete cases by age",
#        fill = "Complete case",
#        x = "Age (y)",
#        y = "Proportion of individuals with complete test data") +
#   scale_fill_viridis_d(labels = c("No","Yes")) +
#   theme_minimal()

# cc_age

# ggsave("cc_age.png",
#        cc_age,
#        path = "plots/misc",
#        width = 5.5,
#        height = 4)

# #complete case by deprivation
# cc_dep = temp %>%
#   ggplot() +
#   geom_bar(mapping = aes(x = (c_en_nzdep_q + 3), fill = factor(cc)), 
#            position = "fill",
#            width = 1) +
#   labs(title = "Proportion of complete cases by age",
#        fill = "Complete case",
#        x = "Deprivation quintile",
#        y = "Proportion of individuals with complete test data") +
#   scale_fill_viridis_d(labels = c("No","Yes")) +
#   theme_minimal()

# cc_dep
# ggsave("cc_dep.png",
#        cc_dep,
#        path = "plots/misc",
#        width = 5.5,
#        height = 4)

# complete.cases(temp)


# ################

# # 4. predictor distributions (Part 1)
# source("0_helper_functions.R")
# library(patchwork)

# an_cohort_ts = read.fst("data/cohorts/an_cohort_ts.fst",
#           as.data.table = TRUE)


# final_plot = test_dist_plot(an_cohort_ts)
# final_plot
# png("plots/testsafe_plots/an_cohort_test_dist.png",
#     width = 900,
#     height = 300)
#    # res = 1200)
# final_plot
# dev.off()

# # ggsave("an_cohort_test_dist.png",
# #        file = final_plot,
# #        path = "plots/testsafe_plots",
# #        height = 4,
# #        width = 8)


# # # by age predicto distributions
# # ggplot(an_cohort) +
# #   geom_bar(mapping = aes(x = nhi_age, fill = factor(c_en_nzdep_q+3)),position = "fill",width = 1 ) +
# #   labs(title = "Distribution of deprivation status in different ages",
# #        y = "Proportion",
# #        x = "Age (y)",
# #        fill = "Deprviation quintile") +
# #   scale_fill_viridis(option = "magma", end = 0.7, discrete = T) +
# #   theme_bw()
# # 
# # #
# # ggplot(an_cohort) +
# #   geom_bar(mapping = aes(x = nhi_age, fill = factor(eth_cat)),position = "fill",width = 1 ) +
# #   labs(title = "Ethnicity in different ages",
# #        y = "Proportion",
# #        x = "Age (y)",
# #        fill = "Ethnicity") +
# #   theme_bw()
# # 
# # #different predictors
# # an_cohort %>% 
# #   as_tibble() %>%
# #   group_by(nhi_age) %>%
# #   summarize(diabetes = mean(hx_vdr_diabetes),
# #             af = mean(hx_af)) %>%
# #   ggplot() +
# #   geom_smooth(mapping = aes(x = nhi_age, y = diabetes),color = "blue" ) +
# #   geom_smooth(mapping = aes(x = nhi_age, y = af),color = "red" ) +
# #   labs(title = "Ethnicity in different ages",
# #        y = "Proportion",
# #        x = "Age (y)") + 
# #   theme_bw()




