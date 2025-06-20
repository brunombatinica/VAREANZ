# testsafe analysis #####OLD CODE NOT CEHCKED THOUGH# ######

# Clearing workspace
rm(list = ls())
dev.off()

#set working directory
setwd("//uoa.auckland.ac.nz/Shared/MED/EPBI/ViewData/Users/bbat644/Desktop/code")

source("0_helper_functions.R")

#importing testsafe information
all_tests_raw = as_tibble(read.fst("data/testsafe/all_tests.fst"))
# levels(all_tests_raw$test)
# levels(all_tests_raw$test) = list(egfr = "eGFR",
#                 hba1c = "HbA1c",
#                 tchdl = "TC/HDL ratio",
#                 tri = "TG")
all_tests_raw$test[all_tests_raw$test == "egfr"] = "eGFR"
all_tests_raw$test[all_tests_raw$test == "hba1c"] = "HbA1c"
all_tests_raw$test[all_tests_raw$test == "tchdl"] = "TC/HDL ratio"
all_tests_raw$test[all_tests_raw$test == "tri"] = "TG"


#import an_cohort (without the testsafe exclusions)
an_cohort_raw = as_tibble(read.fst("data/cohorts/an_cohort_raw.fst"))

#sanity check
sum(!(all_tests_raw$VSIMPLE_INDEX_MASTER %in% an_cohort_raw$VSIMPLE_INDEX_MASTER))

# 0. plot helper functions #### 
# DATESTUFFF 

#create a list of dates
t0 = as.Date("2012-12-31", "%Y-%m-%d")
x_dates = c()
years = 2003:2025
for (i in years){
  #
  datei = as.Date( paste(as.character(i),"-01-01", sep = ""), "%Y-%m-%d")
  x_dates = append(x_dates, datei)
}
#Creating reusable scales 
x_dates
x_scale_date = scale_x_date(limits = c(min(x_dates),max(x_dates)),
                            breaks = seq.Date(from = min(x_dates), 
                                              to = max(x_dates),
                                              by = "5 years"),
                            minor_breaks = "1 year",
                            date_labels = "%Y")

y_scale = scale_y_continuous(
  label = scales::comma
)
y_scale_log = scale_y_log10(
  label = scales::comma
)

# Predictor missingness ######################

an_tests = all_tests_raw %>%
  filter(VSIMPLE_INDEX_MASTER %in% an_cohort$VSIMPLE_INDEX_MASTER) %>%
  select(-FACILITY_CODE) 

#group and arrange by VSIMPLE_INDEX_MASTER 
by_id = an_tests %>%
  group_by(VSIMPLE_INDEX_MASTER)

included_tests = by_id |>
  filter(RESULT_DATE > as.Date(paste("2008","-01-01",sep = ""),
   RESULT_DATE < as.Date(paste("2013","-01-01",sep = ""))))

temp = included_tests |> 
  reframe(
    "HbA1c" = "HbA1c" %in% test,
    "eGFR" = "eGFR" %in% test,
    "TC/HDL ratio" = "TC/HDL ratio" %in% test,
    "TG" = "TG" %in% test,
    "all" = all("HbA1c" %in% test,"eGFR" %in% test,"TC/HDL ratio" %in% test,"TG" %in% test),
    "Any" = any("HbA1c" %in% test,"eGFR" %in% test,"TC/HDL ratio" %in% test,"TG" %in% test))

temp[-1] |>
  summarize_all(sum) |>
  mutate(across(everything(), ~. / 805807))
 
an_cohort

cc_list = list()
look_back = c("2004","2008","2011","2012")
look_forward = c("2013","2014","2015","2016")
for (i in look_back){
  cc_vector = c() 
  for (j in look_forward){
    
    included_tests = by_id %>%
      filter(RESULT_DATE > as.Date(paste(i,"-01-01",sep = "")), RESULT_DATE < as.Date(paste(j,"-01-01",sep = "")))
    
    cc = included_tests %>%
      summarize(all("HbA1c" %in% test,"eGFR" %in% test,"TC/HDL ratio" %in% test,"TG" %in% test)) #by test id
    
    sum(cc$`all(...)`)
    
    cc_vector = append(cc_vector,sum(cc$`all(...)`))
    
  }
  
  cc_list = append(cc_list,list(cc_vector))
  
}


cc_plot_dt = as.data.table(do.call(rbind,cc_list))/dim(an_cohort)[1]
names(cc_plot_dt) = look_forward
cc_plot_dt[,id :=look_back]
cc_plot_dt
cc_plot_dt = melt(cc_plot_dt, id.vars = c("id"),
                  measure.vars = look_forward)
cc_plot_dt

cc_plot = ggplot(cc_plot_dt) +
  geom_line(aes(x = id, y = value, colour = variable, group = variable), linewidth = 1.5) + 
  labs(x = "Look back year",
       y = "Proportion of cohort with complete laboratory data",
       color = "Look forward year",
       #title = "Proportion of complete cases at \n different cut off for test inclusion"
       ) +
  scale_color_discrete(limits = c("2016","2015","2014","2013")) +
  theme_bw()

ggsave("cc_plot.png",
       cc_plot,
       path = "plots/testsafe_plots",
       width = 5,
       height = 5)


# 3. distributions

an_tests_0815 = an_tests %>%
  filter(and(
    RESULT_DATE < as.Date("2015-01-01"),
    RESULT_DATE >= as.Date("2008-01-01")))

# initial exploration
by_id = an_tests_0815 %>%
  group_by(VSIMPLE_INDEX_MASTER) 

dim(by_id)

#each test
table(by_id$test)

#number per person
test_pp = by_id %>% 
          summarize(egfr = sum(test == "egfr"),
            hba1c = sum(test == "hba1c"),
            tchdl = sum(test == "tchdl"),
            tri = sum(test == "tri"))

summary(test_pp$egfr)
summary(test_pp$hba1c)
summary(test_pp$tchdl)
summary(test_pp$tri)

sum(test_pp$egfr == 0) / 805509
sum(test_pp$hba1c == 0) / 805509
sum(test_pp$tchdl == 0) / 805509
sum(test_pp$tri == 0) / 805509



#Freq poly
ts_tests %>%
  ggplot() +
  geom_freqpoly(mapping = aes(x = RESULT_DATE, col = test),
                breaks = x_dates, size = 1.4) +
  labs(
    title = "Test incidence by year",
    x = "Year",
    y = "Number of tests"
  ) +
  x_scale_date +
  y_scale +
  scale_color_manual(values = c(
    egfr = "orangered",
    hba1c = "royalblue1",
    tchdl = "lawngreen"
  )) 

ggsave("#ts_incidence_line.png")



# 4. Proportion of results missing #############################
by_id = ts_tests %>%
  group_by(VSIMPLE_INDEX_MASTER)

template = data.frame(matrix(NA,4,4))
names(template) = c("Subset","Criteria","n","perc")
template[1,2] = "Any time"
template[2,2] = "Prior to study"
template[3,2] = "prior to 2014"

conditions = c(
  "any" = TRUE,
  "egfr" = expression(("egfr" %in% test)),
  "hba1c" = expression(("hba1c" %in% test)),
  "tchdl" = expression(("tchdl" %in% test)),
  "all" = expression(("egfr" %in% test & "hba1c" %in% test & "tchdl" %in% test))
)

output = data.frame()
for (description in names(conditions)){
  #FILTER TERMS ARE COMPUTATIONALLY HORRIBLE BUT VERY ELEGANT TO WRITE
  results = template
  results[,1] = description
  i = conditions[description]
  
  #any time
  results[1,3] = by_id %>%
    summarize(cond = eval(i)) %>%
    filter(cond) %>%
    count()
  
  #prior to study commencement
  results[2,3] = by_id %>%
    filter(RESULT_DATE < as.Date("2014-01-01") ) %>%
    summarize(cond = eval(i)) %>%
    filter(cond) %>%
    count()
  
  #prior to 1 year post study commencement
  results[3,3] = by_id %>%
    filter(RESULT_DATE < as.Date("2014-01-01") ) %>%
    summarize(cond = eval(i)) %>%
    filter(cond) %>%
    count()
  
  output = rbind(output,results)
}

output[,4] = (output[,3]/dim(ts_cohort)[1] * 100)

# # missing results by age (using bound_ts)
# # try to write it a little nicer
# # slight bodge - age centred at 48.3
# # 5 year age bands
# age_breaks = seq(35,75,5)
# age_label = c(paste("<",as.character(age_breaks[1]),sep = ""))
# for (i in 2:length(age_breaks)){
#   age_label = c(age_label,
#                 paste(age_breaks[i-1],age_breaks[i]-1,sep = "-")
#   )
# }
# age_label = c(age_label,
#               paste(">",age_breaks[length(age_breaks)],sep = ""))
# age_breaks = c(0,age_breaks,100)
# 
# 
# by_id_age = tibble(bound_ts) %>%
#   mutate(age = cut((c_nhi_age + 48.3),age_breaks,label = age_label,right = FALSE,include.lowest = FALSE)) %>%
#   select(VSIMPLE_INDEX_MASTER,age,RESULT_DATE,test) %>%
#   group_by(VSIMPLE_INDEX_MASTER)
# 
# #each id should only have 1 age assocaited with them - can check against TS_cohort_test
# ntests = by_id_age %>%
#   summarize(age = first(age), 
#             n_total = n(),
#             all = )
# 
# #cut test
# cut((test = c(34,35,39,40,55,59,60,74,75,99,1,-1)),age_breaks,label = age_label,right = FALSE,include.lowest = FALSE)
# 
# #plot excluding all people with no tests!
# ggplot(data = by_id_age) +
#   geom_bar(aes(x = age)) +
#   y_scale
# 
# #compared with full ts cohort

#plot quick bargraph
ggplot(data = output, aes(x = Subset,y = perc, fill = Criteria)) +
  geom_col(
    position = "dodge"
  ) +
  labs(
    title = "Percentage of TS cohort with results present",
    x = "Test subsets",
    y = "Percentage %"
  ) +
theme_bw()

ggsave(file = "test_available.png")
output





# #sanity checks
# scr %>% group_by(VSIMPLE_INDEX_MASTER)
# scr[scr$RESULT_DATE < as.Date("2013-01-01"),] %>% group_by(VSIMPLE_INDEX_MASTER)
# by_id %>%
#   filter(RESULT_DATE < as.Date("2013-01-01") ) %>%
#   summarize(cond = eval(i)) %>%
#   filter(cond) %>%
#   count()
# 
# 
# #TEST
# #prior to study commencement
# i = expression(all("egfr" %in% test, "hba1c" %in% test, "tchdl" %in% test))
# by_id %>%
#   filter(RESULT_DATE < as.Date("2013-01-01") ) %>%
#   summarize(cond = eval(i)) %>%
#   filter(cond) %>%
#   count()
# 343844/dim(ts_cohort)[1]





# number of tests per person ####
#Highlightfirst 3 rows!!!!
tests_pp = by_id %>%
  summarize(n_tests = n()) %>%
  group_by(n_tests) %>%
  summarize(n_people = n())

#subset tests per person
stests_pp = alltests %>%
  group_by(test,VSIMPLE_INDEX_MASTER) %>%
  summarize(n_tests = n())


#plotting
ggplot() +
  geom_point(tests_pp, mapping = aes(x = n_tests, y = n_people)) +
  labs(
    title = "Number of tests per person",
    x = "log number of tests",
    y = "number of people"
  ) +
  scale_x_continuous() +
  y_scale

ggsave(file = "testpp.png")




# 5. comparison of baseline demographics! - use old script ####



# 6. Proportion results missing (graphs) ######

by_id = ts_tests %>%
  group_by(VSIMPLE_INDEX_MASTER)

# id_pre_2014 = by_id %>%
#   filter(RESULT_DATE < as.Date("2014-01-01") ) %>%
#   summarize(
#     any = 1,
#     egfr = "egfr" %in% test,
#     hba1c = "hba1c" %in% test,
#     tchdl = "tchdl" %in% test,
#     all = ("egfr" %in% test & "hba1c" %in% test & "tchdl" %in% test))

ntest_ = by_id %>%
  summarize(
    any_ = n(),
    egfr_ = sum(test == "egfr"),
    hba1c_ = sum(test == "hba1c"),
    tchdl_ = sum(test == "tchdl"),
    all_ = ("egfr" %in% test & "hba1c" %in% test & "tchdl" %in% test))

ntest_pre_2013 = by_id %>%
  filter(RESULT_DATE < as.Date("2013-01-01") ) %>%
  summarize(
    any_2013 = n(),
    egfr_2013 = sum(test == "egfr"),
    hba1c_2013 = sum(test == "hba1c"),
    tchdl_2013 = sum(test == "tchdl"),
    all_2013 = ("egfr" %in% test & "hba1c" %in% test & "tchdl" %in% test))

ntest_pre_2015 = by_id %>%
  filter(RESULT_DATE < as.Date("2015-01-01") ) %>%
  summarize(
    any_2015 = n(),
    egfr_2015 = sum(test == "egfr"),
    hba1c_2015 = sum(test == "hba1c"),
    tchdl_2015 = sum(test == "tchdl"),
    all_2015 = ("egfr" %in% test & "hba1c" %in% test & "tchdl" %in% test))

ntest_ #number of tests any time by id
ntest_pre_2013 #number of tests pre 2013 by id
ntest_pre_2015 #number of tests pre 2015 by id

#export these data sets
#write.fst(ntest_,"//uoa.auckland.ac.nz/Shared/MED/EPBI/ViewData/Users/bbat644/Desktop/R Code/TestSafe_copy/ntest_.fst")
#write.fst(ntest_pre_2014,"//uoa.auckland.ac.nz/Shared/MED/EPBI/ViewData/Users/bbat644/Desktop/R Code/TestSafe_copy/ntest_pre_2014.fst")

ts_cohort_tests = ts_cohort %>%
  left_join(ntest_,by = "VSIMPLE_INDEX_MASTER") %>%
  left_join(ntest_pre_2013, by = "VSIMPLE_INDEX_MASTER") %>%
  left_join(ntest_pre_2015, by = "VSIMPLE_INDEX_MASTER") %>%
  replace(is.na(.),0)

dim(ts_cohort_tests)
length(ts_cohort_tests$any_ > 0)
sum(ts_cohort_tests$any_ > 0)
sum(ts_cohort_tests$all_ > 0)

rm( list = c("ntest_","ntest_pre_2013","ntest_pre_2015"))

#7. Missingness by 5 year age groups
#missingness out of total ts cohort

# 5 year age bands
age_breaks = seq(30,75,5)
age_label = c(paste("<",as.character(age_breaks[1]),sep = ""))
for (i in 2:length(age_breaks)){
  age_label = c(age_label,
                paste(age_breaks[i-1],age_breaks[i]-1,sep = "-")
  )
}
age_label = c(age_label,
              paste(">",age_breaks[length(age_breaks)],sep = ""))
age_breaks = c(0,age_breaks,100)

ts_age_group = tibble(ts_cohort_tests) %>%
  mutate(age_group = cut((nhi_age),age_breaks,label = age_label, right = FALSE,include.lowest = FALSE)) 

names(ts_age_group)

# #all results
# ggplot(data = ts_age_group) +
#   geom_bar(mapping = aes(x = age_group, fill = as.factor(all_))) +
#   geom_text(aes(x = age_group, label = ..count..), stat = "count", vjust = 1.5, colour = "white") +
#   y_scale
# 
# #any results
# ggplot(data = ts_age_group) +
#   geom_bar(mapping = aes(x = age_group, fill = as.factor(ifelse(any_ > 0,1,0)))) +
#   geom_text(aes(x = age_group, label = ..count..), stat = "count", vjust = 1.5, colour = "white") +
#   y_scale 

#smarter way of doing it

#automatic colours

plot_data_ts = ts_age_group %>%
  group_by(age_group) %>%
  summarize(n = n(),
            any = sum(any_ > 0),
            all = sum(all_ > 0),
            any2013 = sum(any_2013 > 0),
            all2013 = sum(all_2013 >0),
            any2015 = sum(any_2015 > 0),
            all2015 = sum(all_2015 > 0))

# 
# #full TS cohort
# colours1 = c(
#   "Cohort (Test Safe)" = "#F8766D",
#   "All lab results at any time" = "#619CFF",
#   "All labs results prior to 2014" = "#00BA38"
# )
# 
# ggplot(data = plot_data_ts, aes(x = age_group)) +
#   geom_col( aes( y = n, fill = names(colours1)[1])) +
#   geom_text( aes( y = n, label = n), vjust = -1) +
#   geom_col( aes( y = all, fill = names(colours1)[2])) +
#   geom_text( aes( y = all, label = round(all/n,2), fontface = "bold"),colour = "white", vjust = 1.5) +
#   geom_col( aes(y = all2014, fill = names(colours1)[3])) + 
#   geom_text(aes(y = all2014, label = round(all2014/n,2), fontface = "bold"), colour = "white", vjust = 1.5) +
#   labs(
#     title = "Proportion of whole ts cohort with complete lab results",
#     x = "Age groups",
#     y = "Number of people",
#     fill = "Legend"
#   ) +
#   scale_fill_manual(values = colours1)
# 
# ggsave(file = "#Missing_byage_TS_cohort.png",  width = 12, height = 12)


# any result 2015 (test cohort)SS
colours1 = c(
  "Any lab" = "#F8766D",
  "Any lab prior to 2015" = "#619CFF",
  "All labs prior to 2015" = "#00BA38"
)

ggplot(data = plot_data_ts, aes(x = age_group)) +
  geom_col( aes( y = any, fill = names(colours1)[1])) +
  geom_text( aes( y = any, label = n), vjust = -1) +
  geom_col( aes( y = any2015, fill = names(colours1)[2])) +
  geom_col( aes(y = all2015, fill = names(colours1)[3])) + 
  geom_text( aes( y = any2015, label = round(any2015/n,2), fontface = "bold"),colour = "white", vjust = 1.5) +
  geom_text(aes(y = all2015, label = round(all2015/n,2), fontface = "bold"), colour = "white", vjust = 1.5) +
  labs(
    title = "Proportion of test cohort who have at least 1 lab before 2015",
    x = "Age groups",
    y = "Number of people",
    fill = "Legend"
  ) +
  scale_fill_manual(values = colours1)

ggsave(file = "#any_result_2015.png",  width = 12, height = 12)


# proportion of those with at least 1 test (test_cohort)
colours2 = c(
  "Cohort (those with at least 1 result)" = "#F8766D",
  "All lab results at any time" = "#619CFF",
  "All labs results prior to 2014" = "#00BA38"
)

ggplot(data = plot_data_ts, aes(x = age_group)) +
  geom_col( aes( y = any, fill = names(colours2)[1])) +
  geom_text( aes( y = any, label = any), vjust = -1) +
  geom_col( aes( y = all, fill = names(colours2)[2])) +
  geom_text( aes( y = all, label = round(all/any,2), fontface = "bold"),colour = "white", vjust = 1.5) +
  geom_col( aes(y = all2014, fill = names(colours2)[3])) + 
  geom_text(aes(y = all2014, label = round(all2014/any,2), fontface = "bold"), colour = "white", vjust = 1.5) +
  labs(
    title = "Proportion of those with at least 1 result with complete lab data",
    x = "Age groups",
    y = "Number of people",
    fill = "Legend"
  ) +
  scale_fill_manual(values = colours2)

ggsave(file = "#Missing_byage_test_cohort.png", width = 9, height = 6)

# 2014 vs 2015 ####
colours3 = c(
  "Any lab results" = "#F8766D",
  "All lab results" = "#619CFF",
  "All lab results prior to 2015" = "#C77CFF",
  "All labs results prior to 2013" = "#00BA38"
)

ggplot(data = plot_data_ts, aes(x = age_group)) +
  geom_col( aes( y = any, fill = names(colours3)[1])) +
  geom_col( aes( y = all, fill = names(colours3)[2])) +
  geom_col( aes( y = all2015, fill = names(colours3)[3])) +
  geom_col( aes(y = all2014, fill = names(colours3)[4])) + 
  geom_text( aes( y = all2015, label = round(all2015/any,2), fontface = "bold"),colour = "white", vjust = 1.5) +
  geom_text(aes(y = all2014, label = round(all2014/any,2), fontface = "bold"), colour = "white", vjust = 1.5) +
  labs(
    title = "2014 cut off vs 2015 cut off",
    x = "Age groups",
    y = "Number of people",
    fill = "Legend"
  ) +
  scale_fill_manual(values = colours3)

# ggplot(data = plot_data_ts, aes(x = age_group)) +
#   geom_col( aes( y = any, fill = names(colours3)[1])) +
#   geom_text( aes( y = any, label = any), vjust = -1) +
#   geom_col( aes( y = all, fill = names(colours3)[2])) +
#   geom_col( aes( y = all2015, fill = names(colours3)[3])) +
#   geom_col( aes(y = all2014, fill = names(colours3)[4])) + 
#   geom_text( aes( y = all2015, label = round(all2015/any,2), fontface = "bold"),colour = "white", vjust = 1.5) +
#   geom_text(aes(y = all2014, label = round(all2014/any,2), fontface = "bold"), colour = "white", vjust = 1.5) +
#   labs(
#     title = "2014 cut off vs 2015 cut off",
#     x = "Age groups",
#     y = "Number of people",
#     fill = "Legend"
#   ) +
#   scale_fill_manual(values = colours3)

#ggsave(file = "#2014_vs_2015.png", width = 9, height = 6)

#sanity checks
ages = ts_cohort$nhi_age
sum(ages %in% 35:39)







#8. earliest result graphs ####
# first_hba1c = by_id %>%
#   filter(test == "hba1c") %>%
#   summarize(dhb = first(DHB_name),
#             age = first(nhi_age),
#              first_hba1c = min(RESULT_DATE)) 
# 
# first_egfr = by_id %>%
#   filter(test == "egfr") %>%
#   summarize( dhb = first(DHB_name),
#              age = first(nhi_age),
#              first_egfr = min(RESULT_DATE))
# 
# first_tchdl = by_id %>%
#   filter(test == "tchdl") %>%
#   summarize( dhb = first(DHB_name),
#              age = first(nhi_age),
#              first_tchdl = min(RESULT_DATE))
# 
# #cumulative sum calculations
# sum_hba1c = first_hba1c %>%
#   summarize(date_ = unique(first_hba1c),
#             csum = ecdf(first_hba1c)(unique(first_hba1c)) *n())
# 
# sum_egfr = first_egfr %>%
#   summarize(date_ = unique(first_egfr),
#             csum = ecdf(first_egfr)(unique(first_egfr)) *n())
# 
# sum_tchdl = first_tchdl %>%
#   summarize(date_ = unique(first_tchdl),
#             csum = ecdf(first_tchdl)(unique(first_tchdl)) *n())
# 
# #cummulative sum plots
# colours4 = c(
#   "hba1c" = "#619CFF",
#   "egfr" = "#F8766D",
#   "tchdl" = "#00BA38"
# )
# 
# ggplot() +
#   geom_line(sum_hba1c, mapping = aes(x = date_, y = csum, 
#                                     color = names(colours4)[1]),
#             size = 1.5) +
#   geom_line(sum_egfr, mapping = aes(x = date_, y = csum, 
#                                      color = names(colours4)[2]),
#             size = 1.5) +
#   geom_line(sum_tchdl, mapping = aes(x = date_, y = csum, 
#                                      color = names(colours4)[3]),
#             size = 1.5) +
#   y_scale +
#   x_scale_date +
#   labs( 
#     title = "cumulative some of people with at least 1 test",
#     colour = "legend")
#   scale_colour_manual( values = colours4)
#   
# ggsave(file = "csumtest.png")
# 
# rm(list = c("first_hba1c","first_egfr","first_tchdl","sum_hba1c","sum_egfr","sum_tchdl"))

#cumulative sums
test_cohort = ts_cohort %>%
  filter(VSIMPLE_INDEX_MASTER %in% ts_tests$VSIMPLE_INDEX_MASTER) 

first_hba1c = by_id %>%
  filter(test == "hba1c") %>%
  summarize(first_hba1c_date = min(RESULT_DATE),
            first_hba1c = nth(result,which.min(RESULT_DATE))) 

first_egfr = by_id %>%
  filter(test == "egfr") %>%
  summarize(first_egfr_date = min(RESULT_DATE),
            first_egfr = nth(result,which.min(RESULT_DATE)))

first_tchdl = by_id %>%
  filter(test == "tchdl") %>%
  summarize(first_tchdl_date = min(RESULT_DATE),
            first_tchdl = nth(result,which.min(RESULT_DATE)))

first_results = tibble(test_cohort) %>%
  left_join(first_hba1c, by = "VSIMPLE_INDEX_MASTER") %>%
  left_join(first_egfr, by = "VSIMPLE_INDEX_MASTER") %>%
  left_join(first_tchdl, by = "VSIMPLE_INDEX_MASTER")

rm(list = c("first_hba1c","first_egfr","first_tchdl"))

csum = tibble(first_results) %>%
  summarize(date_ = unique(na.omit(c(first_egfr_date,first_hba1c_date,first_tchdl_date))),
            sum_hba1c = ecdf(first_hba1c_date[!is.na(first_hba1c_date)])(date_) 
            * sum(!is.na(first_hba1c_date)),
            sum_egfr = ecdf(first_egfr_date[!is.na(first_egfr_date)])(date_)
            * sum(!is.na(first_egfr_date)),
            sum_tchdl = ecdf(first_tchdl_date[!is.na(first_tchdl_date)])(date_) 
            * sum(!is.na(first_tchdl_date)))

csum_dhb = tibble(first_results) %>%
  group_by(DHB_name) %>%
  summarize(date_ = unique(na.omit(c(first_egfr_date,first_hba1c_date,first_tchdl_date))),
            sum_hba1c = ecdf(first_hba1c_date[!is.na(first_hba1c_date)])(date_) 
            * sum(!is.na(first_hba1c_date)),
            sum_egfr = ecdf(first_egfr_date[!is.na(first_egfr_date)])(date_)
            * sum(!is.na(first_egfr_date)),
            sum_tchdl = ecdf(first_tchdl_date[!is.na(first_tchdl_date)])(date_) 
            * sum(!is.na(first_tchdl_date)))
#csum by dhb
colours4 = c(
  "hba1c" = "#619CFF",
  "egfr" = "#F8766D",
  "tchdl" = "#00BA38"
)

plots = list()
plot_names = c("Northland","Waitemata","Auckland","Counties Manukau")
for (i in plot_names){
  
  tempplot = csum_dhb %>%
    filter(DHB_name == i) %>% 
    ggplot() +
    geom_line(aes(x = date_, y = sum_hba1c, color = names(colours4)[1]), size = 1.5) +
    geom_line(aes(x = date_, y = sum_egfr, color = names(colours4)[2]), size = 1.5) +
    geom_line(aes(x = date_, y = sum_tchdl,  color = names(colours4)[3]), size = 1.5) +
    geom_hline( yintercept = length(unique(ts_cohort[ts_cohort$DHB_name == i,]$VSIMPLE_INDEX_MASTER)),
                show.legend = T) +
    geom_hline( yintercept = length(unique(test_cohort[test_cohort$DHB_name ==i,]$VSIMPLE_INDEX_MASTER)),
                linetype = "dashed",show.legend = T) +
    x_scale_date +
    y_scale +
    labs(
      title = i, 
      colour = "Legend",
      x = "Date",
      y = "Cummultive sum") +
    scale_color_manual(values = colours4)
  
  plots = append(plots,list(tempplot))
}

ggarrange(plotlist = plots, ncol = 1, nrow = 4, 
          common.legend = TRUE, legend = "bottom") %>%
  annotate_figure(top = text_grob("Cummulative sum of individuals with a given result", 
                                  face = "bold", size = 14))

ggsave(file = "#csum_id_dhb.png", width = 8, height = 10)

#Individual results
ggplot(csum) +
  geom_line(aes(x = date_, y = sum_hba1c, color = names(colours4)[1]), size = 1.5) +
  geom_line(aes(x = date_, y = sum_egfr, color = names(colours4)[2]), size = 1.5) +
  geom_line(aes(x = date_, y = sum_tchdl,  color = names(colours4)[3]), size = 1.5) +
  geom_hline( yintercept = length(unique(ts_cohort$VSIMPLE_INDEX_MASTER)),show.legend = T) +
  geom_hline( yintercept = length(unique(test_cohort$VSIMPLE_INDEX_MASTER)),
              linetype = "dashed",show.legend = T) +
  x_scale_date +
  y_scale +
  labs(
    title = "cumulative some of people with at least 1 test", 
    colour = "legend") +
  scale_color_manual(values = colours4) 
#ggsave(file = "#csum_individual.png")











#9. log reg to predict missingness ####
names(ts_cohort_tests)
ts_cohort_tests = ts_cohort_tests %>%
  mutate(missing_any = as.numeric(all_ == 0),
         missing_all = as.numeric(any_ == 0),
         missing_any_2014 = as.numeric(all_2014 == 0),
         missing_hba1c_2014 = as.numeric(hba1c_2014 == 0),
         missing_egfr_2014 = as.numeric(egfr_2014 == 0),
         missing_tchdl_2014 = as.numeric(tchdl_2014 == 0))

head(ts_cohort_tests)
#logistic regression on what is missing
missing_outcome = "missing_any"
missing_predictors = c("eth_cat",
                       "DHB_name",
                       "gender_code",
                       "I(nhi_age-48)", #Dont need to scale as no interaction terms
                       "c_en_nzdep_q",
                       "ph_bp_lowering_prior_6mths",
                       "ph_lipid_lowering_prior_6mths",
                       "ph_antip_antic_6mths",
                       "indicator") #?using imp_fatala and nonfata cvd here

missing_formula = as.formula(
  paste(missing_outcome," ~ ",
        paste(missing_predictors,collapse = "+"),
        sep = "")) 
missing_formula
missing_fit = glm(missing_formula, family = "binomial", data = ts_cohort_tests)
summary(missing_fit)
exp(summary(missing_fit)$coef)





#10. Variability among results ####
#difference between highest and lowest result
#hba1c
hba1c_var = bound_ts %>%
  filter(test == "hba1c") %>%
  group_by(VSIMPLE_INDEX_MASTER) %>%
  summarize(n = n(),
            min = range(result)[1],
            max = range(result)[2],
            var_ = var(result),
            range_ = (range(result)[2]-range(result)[1]))

hba1c_var %>%
  arrange(desc(range_))

ggplot(hba1c_var) +
  geom_histogram(aes(x = range_),binwidth = 1)

#variability vs number of tests
ggplot(hba1c_var) +
  geom_point(aes(x = hba1c_var$n, y = hba1c_var$range_))

#differences between index and closest possible result
hba1c_temp = bound_ts %>%
  filter(test == "hba1c" ) %>%
  group_by(VSIMPLE_INDEX_MASTER) %>%
  filter(n() > 1) %>% 
  group_split()

n = length(hba1c_temp)
sample_index = sample(1:n,10000)

l = c()
for (i in hba1c_temp){
  i = i %>%
    arrange(RESULT_DATE)
  #print(i)
  if (any(difftime(i$RESULT_DATE,t0, units = "days") < 0)){
    ind = which.min(difftime(i$RESULT_DATE,t0, units = "days"))
    if (ind != 0){
      pre_t0 = i$result[ind]
      post_t0 = tryCatch(i$result[ind+1])
      l = tryCatch(c(l, abs(pre_t0 - post_t0)))
    }
  }
}
l

hist(l)


#11. plot first hba1c results




#11. Lab cov plots ####

#closest result - for each id gives closest test result to t0
closest_hba1c = by_id %>%
  filter(test == "hba1c") %>%
  summarize( hba1c_date_ = RESULT_DATE[which.min(abs(difftime(RESULT_DATE,t0,units = "days")))],
             hba1c_ = result[which.min(abs(difftime(RESULT_DATE,t0,units = "days")))]) 

furthest_hba1c = by_id %>%
  filter(test == "hba1c") %>%
  summarize( hba1c_date_aux = RESULT_DATE[which.max(abs(difftime(RESULT_DATE,t0,units = "days")))],
             hba1c_aux = result[which.max(abs(difftime(RESULT_DATE,t0,units = "days")))]) 

closest_egfr = by_id %>%
  filter(test == "egfr") %>%
  summarize( egfr_date_ = RESULT_DATE[which.min(abs(difftime(RESULT_DATE,t0,units = "days")))],
             egfr_ = result[which.min(abs(difftime(RESULT_DATE,t0,units = "days")))]) 

furthest_egfr = by_id %>%
  filter(test == "egfr") %>%
  summarize( egfr_date_aux = RESULT_DATE[which.max(abs(difftime(RESULT_DATE,t0,units = "days")))],
             egfr_aux = result[which.max(abs(difftime(RESULT_DATE,t0,units = "days")))]) 

closest_tchdl = by_id %>%
  filter(test == "tchdl") %>%
  summarize( tchdl_date_ = RESULT_DATE[which.min(abs(difftime(RESULT_DATE,t0,units = "days")))],
             tchdl_ = result[which.min(abs(difftime(RESULT_DATE,t0,units = "days")))]) 

furthest_tchdl = by_id %>%
  filter(test == "tchdl") %>%
  summarize( tchdl_date_aux = RESULT_DATE[which.max(abs(difftime(RESULT_DATE,t0,units = "days")))],
             tchdl_aux = result[which.max(abs(difftime(RESULT_DATE,t0,units = "days")))]) 

closest_tri = tri %>%
  group_by(VSIMPLE_INDEX_MASTER) %>%
  filter(test == "tri") %>%
  summarize( tri_date_ = RESULT_DATE[which.min(abs(difftime(RESULT_DATE,t0,units = "days")))],
             tri_ = result[which.min(abs(difftime(RESULT_DATE,t0,units = "days")))]) 

furthest_tri = tri %>%
  group_by(VSIMPLE_INDEX_MASTER) %>%
  filter(test == "tri") %>%
  summarize( tri_date_aux = RESULT_DATE[which.max(abs(difftime(RESULT_DATE,t0,units = "days")))],
             tri_aux = result[which.max(abs(difftime(RESULT_DATE,t0,units = "days")))]) 

closest_results = ts_cohort %>%
  left_join(closest_hba1c, by = "VSIMPLE_INDEX_MASTER") %>%
  left_join(furthest_hba1c, by = "VSIMPLE_INDEX_MASTER") %>%
  left_join(closest_egfr, by = "VSIMPLE_INDEX_MASTER") %>%
  left_join(furthest_egfr, by = "VSIMPLE_INDEX_MASTER") %>%
  left_join(closest_tchdl, by = "VSIMPLE_INDEX_MASTER") %>%
  left_join(furthest_tchdl, by = "VSIMPLE_INDEX_MASTER") %>%
  left_join(closest_tri, by = "VSIMPLE_INDEX_MASTER") %>%
  left_join(furthest_tri, by = "VSIMPLE_INDEX_MASTER") 

print(tibble(closest_results), width = Inf)

rm(list = c("closest_hba1c","closest_egfr","closest_tchdl","closest_tri",
            "furthest_hba1c","furthest_egfr","furthest_tchdl","furthest_tri"))

#closest_hba1c
#testing these equations
# test = by_id[1:10000,] %>%
#   arrange(VSIMPLE_INDEX_MASTER)
# test
# 
# test_c = test %>%
#   summarize( resultdate_ = RESULT_DATE[which.min(abs(difftime(RESULT_DATE,t0,units = "days")))],
#              result_ = result[which.min(abs(difftime(RESULT_DATE,t0,units = "days")))])
# 
# test_c
#testing
# by_id_n = by_id %>%
#   mutate(n = n()) %>%
#   arrange(VSIMPLE_INDEX_MASTER)
# 
# print(by_id_n,n = 100)
# hist(closest_hba1c$resultdate_, breaks = "months")
# hist(hba1c_12_14$resultdate_, breaks = "months")
# testing relationship between variables
# hba1c_12_14 %>%
#   mutate(factor = as.factor(gender_code)) %>%
#   ggplot() +
#   geom_density(aes(x = result_, fill = factor), bw = 2)
# 
# hba1c_12_14 %>%
#   mutate(across(gender_code),factor) %>%
#   ggplot() +
#   geom_density(aes(x = result_, fill = factor), bw = 2)

#only take results hba1c from a certain time (2012-2014)
ind = 3
tests = c("hba1c_date_" = "hba1c_",
          "egfr_date_" = "egfr_",
          "tchdl_date_" = "tchdl_",
          "tri_date_" = "tri_")

#Select only results from 2012 - 2014
plot_data = tibble(closest_results) %>%
  filter(.data[[names(tests)[ind]]] > as.Date("2012-01-01") & .data[[names(tests)[ind]]] < as.Date("2014-01-01"))

plot_variables = c("hx_vdr_diabetes",
                   "ph_bp_lowering_prior_6mths",
                   "ph_lipid_lowering_prior_6mths",
                   "ph_antip_antic_6mths",
                   "imp_fatal_cvd",
                   "non_fatal_cvd",
                   "hx_af",
                   "gender_code",
                   "DHB_name",
                   "en_nzdep_q",
                   "eth_cat",
                   "nhi_age",
                   paste(tests[ind],"aux",sep = ""),
                   tests[-ind][1],
                   tests[-ind][2],
                   tests[-ind][3])

plots_ = list()
i = tests[ind]
for (j in plot_variables){
  if (dim(unique(plot_data[,j]))[1] < 8){
    #factor plot
    tempplot = plot_data %>%
      ggplot() +
      geom_density(aes(x = .data[[i]], fill = factor(.data[[j]])), 
                   bw = 0.1, alpha = 0.5, ) + #adjust bw as required
      y_scale + xlim(0,10) +
      labs(
        title = j, 
        colour = "legend"
      ) +
      theme(legend.position = "bottom")
    
    #tempplot
    
  } else {
    #continuous plot
    
    #potentially flip x and y here
    tempplot = plot_data %>%
      ggplot(aes(x = .data[[j]], y = .data[[i]])) +
      stat_bin_2d(binwidth = c(1,1)) +
      y_scale +
      labs(
        title = j
      ) +
      scale_fill_viridis()
  }
  
  plots_ = append(plots_,list(tempplot))
  
}
ggarrange(plotlist = plots_, ncol = min(2,length(plot_variables)), 
          nrow = 1 + (length(plot_variables)-1)%/%2 ) %>%
  annotate_figure(top = text_grob(paste(tests[ind],"covariance plots (results between 2012-2014, 1 result per person)"), 
                                  face = "bold", size = 14))

ggsave(file = paste("#",tests[ind],"_cov_t.png",sep = ""), width = 10, height = 25)

#basic tri values
#dim(tri)
#length(unique(tri$VSIMPLE_INDEX_MASTER))


#Debugging graphs
# plot_data %>%
#   ggplot() +
#   geom_density(aes(x = hba1c_, fill = hx_vdr_diabetes), alpha = 0.5) # +
#   y_scale +
#   xlim(0,10) +
#   labs(
#     title = j, 
#     colour = "legend"
#   ) +
#   theme(legend.position = "bottom")

#testing variables exploration plots
# 
# plot_data %>%
#   ggplot(aes(x= .data[["nhi_age"]], y = result_)) +
#   stat_bin_2d(binwidth = c(1,5)) + 
#   scale_fill_viridis()
# plot_data %>%
#   ggplot() +
#   geom_density(aes(x = result_, fill = factor(.data[["DHB_name"]])), bw = 1, alpha = 0.5) +
#   y_scale +
#   labs(
#     title = i, 
#     colour = "legend"
#   )
# 
# plot_data %>%
#   ggplot(aes(x = .data[["nhi_age"]], y = result_)) +
#   geom_density_2d(stat = "density_2d",
#                   contour_var = "count",
#                   show.legend = T,
#                   bins = 3)
# 
#   stat_density_2d(contour_var = n,
#                   geom = "polygon",
#                   n = 100,
#                   bins = 20)
# 
# plot_data %>%
#   ggplot(aes(x = .data[["nhi_age"]], y = result_)) +
#   stat_density_2d(aes(fill = stat(log(level*1e10))),
#                   geom = "polygon",
#                   n = 100,
#                   breaks = c(1,10,100))
# 
# plot_data %>%
#   ggplot() +
#   geom_point(aes(x = .data[["nhi_age"]], y = result_), position = "jitter") 
# 
# plot_data %>%
#   ggplot() +
#   geom_jitter(aes(x = .data[["nhi_age"]], y = result_)) 
#   
# 
# plot_data %>%
#   ggplot() +
#   geom_density_2d(aes(x = .data[["nhi_age"]], y = result_)) +
#   stat_density_2d()

# 12. simple linear regression to impute labs ####

tests = c("hba1c_date_" = "hba1c_",
          "egfr_date_" = "egfr_",
          "tchdl_date_" = "tchdl_",
          "tri_date_" = "tri_")

ind = 4
full_data = tibble(closest_results) %>%
  filter(.data[[names(tests)[ind]]] > as.Date("2012-01-01") & .data[[names(tests)[ind]]] < as.Date("2014-01-01"))

n = dim(full_data)[1]
set.seed(1)
#test set of 10% can do full CV shebang after 
test_ind = sample(1:n,n/10)

test_data = full_data[test_ind,]
train_data = full_data[-test_ind,]

#just demographics to predict hba1c
names(full_data)
y = tests[ind]
Xp = c("eth_cat", 
       "DHB_name","nhi_age",
       "c_en_nzdep_q",
       "hx_vdr_diabetes",
       "hx_af",
       "ph_bp_lowering_prior_6mths",
       "ph_lipid_lowering_prior_6mths",
       "ph_antip_antic_6mths")
form = as.formula(paste(y,"~",paste(Xp,collapse = " + ")))

#fitting linear regression model #
fit.glm = glm(form, data = train_data, family = "gaussian")
fit.lm = lm(form, data = train_data)

#summary(fit.glm)
#summary(fit.lm)

#make predictions on test set
pred = predict(fit.lm, newdata = test_data, type = "response")

#also including the other lab results
# form2 = as.formula(paste(y,"~",paste(Xp,collapse = " + ")," + ",paste(tests[-ind],collapse = " + ")))
# fit.lm2 = lm(form2, data = train_data)
# summary(fit.lm2)

#using a regression tree
fit.tree = rpart(form, train_data,cp = 0.001) 
print(fit.tree)
pred.tree = predict(fit.tree,test_data)

pred_data = tibble(
  data.frame("test" = test_data[,y], 
             "pred" = pred, 
             "pred.tree" = pred.tree))

#plotting predicted versus real
# names(test_data)
# dim(test_data)
#plot(predict(fit.lm),residuals(fit.lm))
#length(pred)

plot1 = ggplot(data = pred_data, aes(x = pred, y = .data[[y]])) +
  stat_bin_2d(binwidth = c(0.1,0.1)) +
  y_scale +
  scale_fill_viridis()

plot2 = ggplot(data = pred_data, aes(x = pred.tree, y = .data[[y]])) +
  stat_bin_2d(binwidth = c(0.1,0.1)) +
  y_scale +
  scale_fill_viridis()

ggarrange(plot1, plot2, ncol = 1, nrow = 2) %>%
  annotate_figure(top = text_grob(paste(tests[ind],"lin_reg vs tree"), 
                                  face = "bold", size = 14))

ggsave(file = paste("#",tests[ind],"pred.png"), height = 8, width = 5)

#check characteristics of those with high values
# temp = bound_ts %>%
#   filter(test == "hba1c") %>%
#   arrange(desc(result))
# print(temp[1,], width = Inf)
