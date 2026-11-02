# Project title: Metabolic and inflammatory biomarkers trajectories prior to a cancer diagnosis and the risk of subsequent severe infection
# Study population: individuals diagnosed with their first cancer during 1985-2005 (at least three biomarker measurements before the diagnosis)
# Exposure: metabolic and inflammatory biomarkers 
# Outcomes: severe infection, including sepsis and hospitalization for infection as the primary diagnosis (consider both primary and secondary diagnosis in sensitivity analysis)
# covariates: sex, birth year, country of birth, education, employment status, income at the first blood sampling, psychiatric disorders and diabetes
#--------------------------------------------------------------------------------------------------------------------

# The working steps:
# 1. identify the studied biomarkers: at least measured among 20% of cancer patients 
# 2. identify the study population and created the master file with covariates and outcomes
# 3. identify biomarkers trajectories via latent class growth modeling
# 4. analyze the association between biomarker trajectories and the risk of infection using flexible parametric models
# 5. stratified analysis by sex and history of infection before cancer diagnosis

library(dplyr)
library(rmarkdown)
library(markdown)
library(devtools)
library(stringr)
library(tidyr)
library(lcmm)
library(NormPsy)
library(nlme)
library(ggplot2)
library(lattice)
library(haven)
library(survival)
library(survminer)

# step 1: identify the studies biomarkers: at least measured among 20% of cancer patients 

rm(list=ls())
# load the biomarker data
biom <- read.csv(".../biomarkers.csv")
n_distinct(biom$id) 

bios <- c("S_ApoA","S_ApoB","ApoKvot","TC","TG","LDL","HDL","RATHDL","fS_Gluk","S_FAMN","S_Hapt","S_CRP","S_Urat","S_Alb","B_Hb","S_IgA","S_IgG","S_IgM","B_LPKV","S_GLOB")
# subtypes of leukocytes: "Lymfocy","LymfocyP","Monoc","MonocP","Neutrof","NeutrofP","Basof","BasofP","Eos","Eosinof"

results <- NULL
bio <- NULL
for (bio in bios) {
  print(bio)
  tocheck <- NULL
  tocheck <- biom %>%
    filter(!is.na(!!sym(bio)))
  num_ch <- data.frame(
    biomarkers=bio,
    n_test=n_distinct(tocheck$id)
  )
  results <- rbind(results, num_ch)
}

# to check the number of tests of leukocyte subtypes
tocheck <- biom %>%
  filter(!is.na(Lymfocy) | !is.na(LymfocyP))
Lymfocy <- data.frame(
  biomarkers="Lymfocy",
  n_test=n_distinct(tocheck$id)
)
tocheck <- biom %>%
  filter(!is.na(Monoc) | !is.na(MonocP))
Monoc <- data.frame(
  biomarkers="Monoc",
  n_test=n_distinct(tocheck$id)
)
tocheck <- biom %>%
  filter(!is.na(Neutrof) | !is.na(NeutrofP))
Neutrof <- data.frame(
  biomarkers="Neutrof",
  n_test=n_distinct(tocheck$id)
)
tocheck <- biom %>%
  filter(!is.na(Basof) | !is.na(BasofP))
Basof <- data.frame(
  biomarkers="Basof",
  n_test=n_distinct(tocheck$id)
)
tocheck <- biom %>%
  filter(!is.na(Eos) | !is.na(Eosinof))
Eos <- data.frame(
  biomarkers="Eos",
  n_test=n_distinct(tocheck$id)
)

results <- rbind(results, Lymfocy,Monoc,Neutrof,Basof,Eos)

results <- results %>%
  mutate(per=(n_test/665905)*100)
print(results)

# to check if the proportion is similar to that in cancer patients
master <- read_sas(".../cancer_master.sas7bdat") 
cancer <- master %>%
  filter(!is.na(can_date)) %>%
  dplyr::select(id) 
biom_cancer <- cancer %>%
  left_join(biom, by="id")

results_can <- NULL
bio <- NULL
for (bio in bios) {
  print(bio)
  tocheck <- NULL
  tocheck <- biom_cancer %>%
    filter(!is.na(!!sym(bio)))
  num_ch <- data.frame(
    biomarkers=bio,
    n_test_can=n_distinct(tocheck$id)
  )
  results_can <- rbind(results_can, num_ch)
}

# to check the number of tests of leukocyte subtypes
tocheck <- biom_cancer %>%
  filter(!is.na(Lymfocy) | !is.na(LymfocyP))
Lymfocy <- data.frame(
  biomarkers="Lymfocy",
  n_test_can=n_distinct(tocheck$id)
)
tocheck <- biom_cancer %>%
  filter(!is.na(Monoc) | !is.na(MonocP))
Monoc <- data.frame(
  biomarkers="Monoc",
  n_test_can=n_distinct(tocheck$id)
)
tocheck <- biom_cancer %>%
  filter(!is.na(Neutrof) | !is.na(NeutrofP))
Neutrof <- data.frame(
  biomarkers="Neutrof",
  n_test_can=n_distinct(tocheck$id)
)
tocheck <- biom_cancer %>%
  filter(!is.na(Basof) | !is.na(BasofP))
Basof <- data.frame(
  biomarkers="Basof",
  n_test_can=n_distinct(tocheck$id)
)
tocheck <- biom_cancer %>%
  filter(!is.na(Eos) | !is.na(Eosinof))
Eos <- data.frame(
  biomarkers="Eos",
  n_test_can=n_distinct(tocheck$id)
)

results_can <- rbind(results_can, Lymfocy,Monoc,Neutrof,Basof,Eos)
results_can$per_can <- (results_can$n_test_can/169720)*100
results <- full_join(results,results_can,by="biomarkers")
print(results)
# similar trends in the total population and cancer patients
# the following biomarkers (18 biomarkers) are selected for further analyses: "S_ApoA","S_ApoB","ApoKvot","TC","TG","LDL","HDL","RATHDL","fS_Gluk","S_FAMN","S_Hapt","S_CRP","S_Urat","S_Alb","B_Hb","S_IgG","B_LPKV","S_GLOB"
#--------------------------------------------------------------------------------


# step 2: identify the study population and created the master file with covariates and outcomes
rm(list = ls())
#load the master file
master <- read_sas(".../cancer_master.sas7bdat") 
master$id <- as.character(master$id) 
master <- master[,-c(1)]

antidiabetics <- read_sas(".../first_antidiabetecs.sas7bdat")
antidiabetics$id <- as.character(antidiabetics$id) 
antidiabetics <- antidiabetics %>%
  dplyr::select(id, date_first_antidiabetics)

master <- master %>%
  left_join(antidiabetics, by="id") %>%
  mutate(can_date=as.Date(can_date),
         birthdate=as.Date(FODDATn),
         can_year=ifelse(!is.na(can_date),as.numeric(format(can_date, "%Y")), NA),
         first_diabetes=pmin(date_diabetes, date_first_antidiabetics, na.rm = T)) %>%
  mutate(country=ifelse(FodGrEg4=="00", 0, 1),
       age_can=as.integer(can_date - birthdate)/365.25,
       sex=ifelse(Kon=="1",1,2),
       employ_atenroll=ifelse(is.na(employ_atenroll),99,employ_atenroll),
       income_atenroll=ifelse(is.na(income_atenroll),99,income_atenroll),
       education_atenroll=ifelse(is.na(education_atenroll),99,education_atenroll))

# study population: patients diagnosis with their first cancer during 1985-2005
master <- master %>%
  filter(!is.na(can_year) & can_year>1984 & can_year<2006)

# to define subtypes of cancer
master <- master %>%
  mutate(cancer_t1=case_when(
    # by sites
    str_starts(cancer_icd7,"142") ~ 1, # salivary gland
    str_starts(cancer_icd7,"150") ~ 2, # esophagus
    str_starts(cancer_icd7,"151") ~ 3, # stomach
    str_starts(cancer_icd7,"152") ~ 4, # small intestine
    str_starts(cancer_icd7,"153") ~ 5, # colon
    str_starts(cancer_icd7,"154") ~ 6, # rectum
    str_starts(cancer_icd7,"155") ~ 7, # liver
    str_starts(cancer_icd7,"157") ~ 8, # pancreas
    str_starts(cancer_icd7,"162") ~ 9, # lung
    str_starts(cancer_icd7,"170") ~ 10, # breast
    str_starts(cancer_icd7,"171") ~ 11, # cervix
    str_starts(cancer_icd7,"172") | str_starts(cancer_icd7,"173") | str_starts(cancer_icd7,"174") ~ 12, # uterus
    str_starts(cancer_icd7,"175") ~ 13, # ovary
    str_starts(cancer_icd7,"178") ~ 14, # testis
    str_starts(cancer_icd7,"180") ~ 15, # kidney
    str_starts(cancer_icd7,"190") ~ 16, # melanoma
    str_starts(cancer_icd7,"191") ~ 17, # non-melanoma skin
    str_starts(cancer_icd7,"192") ~ 18, # eye
    str_starts(cancer_icd7,"193") ~ 19, # central nervous system (CNS)
    str_starts(cancer_icd7,"194") ~ 20, # thyroid
    str_starts(cancer_icd7,"195") ~ 21, # otherocrine gland
    str_starts(cancer_icd7,"196") ~ 22, # bone
    str_starts(cancer_icd7,"197") ~ 23, # connective tissue
    str_starts(cancer_icd7,"199") ~ 24, # other or unspecified sites
    str_starts(cancer_icd7,"201") ~ 25, # hodgkoin's lymphoma
    str_starts(cancer_icd7,"200") | str_starts(cancer_icd7,"202") | str_starts(cancer_icd7,"2041") ~ 26, # non-hodgkin's lymphoma
    str_starts(cancer_icd7,"2040") | str_starts(cancer_icd7,"2049") ~ 27, # acute lymphoid leukemia (ALL)
    str_starts(cancer_icd7,"2050") | str_starts(cancer_icd7,"2059") | str_starts(cancer_icd7,"2060") | str_starts(cancer_icd7,"2069") ~ 28, # acute myeloid leukemia (AML)
    TRUE ~ NA_real_), 
      
      cancer_t2=case_when(
      #by organ systems
    str_starts(cancer_icd7,"140") | str_starts(cancer_icd7,"141") | str_starts(cancer_icd7,"142") | str_starts(cancer_icd7,"143") | str_starts(cancer_icd7,"144") | str_starts(cancer_icd7,"145") | str_starts(cancer_icd7,"146") | str_starts(cancer_icd7,"147") | str_starts(cancer_icd7,"148") ~ 1, # buccal cavity and pharynx
    str_starts(cancer_icd7,"15") ~ 2, # digestive system
    str_starts(cancer_icd7,"160") | str_starts(cancer_icd7,"161") | str_starts(cancer_icd7,"162") | str_starts(cancer_icd7,"163") | str_starts(cancer_icd7,"164") | str_starts(cancer_icd7,"165") ~ 3, # respiratory system
    str_starts(cancer_icd7,"17") ~ 4, # breast and reproductive system
    str_starts(cancer_icd7,"180") | str_starts(cancer_icd7,"181") ~ 5, # urinary system
    str_starts(cancer_icd7,"20") ~ 6, # hematological malignancy
    str_starts(cancer_icd7,"190") | str_starts(cancer_icd7,"191") | str_starts(cancer_icd7,"192") | str_starts(cancer_icd7,"194") | str_starts(cancer_icd7,"195") | str_starts(cancer_icd7,"196") | str_starts(cancer_icd7,"197") | str_starts(cancer_icd7,"198") | str_starts(cancer_icd7,"199") ~ 7, # other categories
    TRUE ~ NA_real_) 
  )

# to define comorbidities before and after a cancer diagnosis
master <- master %>%
  mutate(cvd_bef_cancer=ifelse(!is.na(date_cvd) & date_cvd < can_date, 1, 0),
         diabete_bef_cancer=ifelse(!is.na(first_diabetes) & first_diabetes < can_date, 1, 0),
         psychi_bef_cancer=ifelse(!is.na(date_psychi) & date_psychi < can_date, 1, 0),
         cvd=ifelse(!is.na(date_cvd) & date_cvd >= can_date, 1, 0),
         ihd=ifelse(!is.na(date_ihd) & date_ihd >= can_date, 1, 0),
         mi=ifelse(!is.na(date_mi) & date_mi >= can_date, 1, 0),
         hf=ifelse(!is.na(date_hf) & date_hf >= can_date, 1, 0),
         af=ifelse(!is.na(date_af) & date_af >= can_date, 1, 0),
         arrhythm=ifelse(!is.na(date_arrhythm) & date_arrhythm >= can_date, 1, 0),
         myocard=ifelse(!is.na(date_myocard) & date_myocard >= can_date, 1, 0),
         cardiomyo=ifelse(!is.na(date_cardiomyo) & date_cardiomyo >= can_date, 1, 0),
         stroke=ifelse(!is.na(date_stroke) & date_stroke >= can_date, 1, 0),
         istroke=ifelse(!is.na(date_istroke) & date_istroke >= can_date, 1, 0),
         subara=ifelse(!is.na(date_subara) & date_subara >= can_date, 1, 0),
         intra=ifelse(!is.na(date_intra) & date_intra >= can_date, 1, 0),
         diabetes=ifelse(!is.na(date_diabetes) & date_diabetes >= can_date, 1, 0))


# to define outcome
infec_both <- read_sas(".../hos_infec_both.sas7bdat")
infec_p <- read_sas(".../hos_infec_p.sas7bdat")
sepsis <- read_sas(".../sepsis.sas7bdat")
death_infec <- read_sas(".../death_infec.sas7bdat")

cancer <- master %>%
  dplyr::select(id,can_date)
infec_both$id <- as.character(infec_both$id) 
infec_both <- infec_both[,2:3]
infec_p$id <- as.character(infec_p$id) 
infec_p <- infec_p[,2:3]
sepsis$id <- as.character(sepsis$id) 
sepsis <- sepsis[,2:3]
death_infec$id <- as.character(death_infec$id) 
death_infec <- death_infec[,2:3]

infec_both_bef <- infec_both %>%
  left_join(cancer, by="id") %>%
  filter(can_date - date_hos_both>0 & can_date - date_hos_both<3650) %>% #infection within 10 years before the cancer diagnosis
  distinct(id) %>%
  mutate(infec_both_bef=1)
infec_both <- infec_both %>%
  left_join(cancer, by="id") %>%
  filter(date_hos_both>=can_date) %>%
  arrange(id,date_hos_both) %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(id, date_hos_both)

infec_p_bef <- infec_p %>%
  left_join(cancer, by="id") %>%
  filter(can_date - date_hos_pri>0 & can_date - date_hos_pri<3650) %>% 
  distinct(id) %>%
  mutate(infec_p_bef=1)
infec_p <- infec_p %>%
  left_join(cancer, by="id") %>%
  filter(date_hos_pri>=can_date) %>%
  arrange(id,date_hos_pri) %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(id, date_hos_pri)

sepsis_bef <- sepsis %>%
  left_join(cancer, by="id") %>%
  filter(can_date - date_sepsis>0 & can_date - date_sepsis<3650) %>% 
  distinct(id) %>%
  mutate(sepsis_bef=1)
sepsis <- sepsis %>%
  left_join(cancer, by="id") %>%
  filter(date_sepsis>=can_date) %>%
  arrange(id,date_sepsis) %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(id, date_sepsis)

master <- master %>%
  left_join(infec_both,by="id") %>%
  left_join(infec_p,by="id") %>%
  left_join(sepsis,by="id") %>%
  left_join(infec_both_bef,by="id") %>%
  left_join(infec_p_bef,by="id") %>%
  left_join(sepsis_bef,by="id") %>%
  left_join(death_infec,by="id") %>%
  mutate(infec_both=ifelse(!is.na(date_hos_both),1,0),
         infec_p=ifelse(!is.na(date_hos_pri),1,0),
         sepsis=ifelse(!is.na(date_sepsis),1,0),
         death_infec=ifelse(!is.na(death_infec),death_infec,0))

saveRDS(master, ".../master_cancer1985_2005.RDS")
#-------------------------------------------------------------------------------


# step 3. identify biomarkers trajectories via latent class growth modeling
rm(list = ls())
biom <- read.csv(".../biomarkers.csv")
biom$id <- as.character(biom$id) 
biom <- biom[,2:45]
master <- readRDS(".../master_cancer1985_2005.RDS")

cancer <- master %>%
  dplyr::select(id, can_date)

cancer_biom <- cancer %>%
  left_join(biom,by="id") %>%
  filter(date_samp <= can_date) %>%
  filter(!is.na(S_ApoA) | !is.na(S_ApoB) | !is.na(ApoKvot) | !is.na(TC) | !is.na(TG) | !is.na(LDL) | !is.na(HDL) | !is.na(RATHDL) | !is.na(fS_Gluk) | !is.na(S_FAMN) | !is.na(S_Hapt) | !is.na(S_CRP) | !is.na(S_Urat) | !is.na(S_Alb) | !is.na(B_Hb) | !is.na(S_IgG) | !is.na(B_LPKV) | !is.na(S_GLOB))
n_distinct(cancer_biom$id)

cancer_biom <- cancer_biom %>%
  mutate(date_samp=as.Date(date_samp),
         can_date=as.Date(can_date)) %>%
  filter(can_date - date_samp>=0 & can_date - date_samp<=3650)
n_distinct(cancer_biom$id)

cancer_biom <- cancer_biom %>%
  group_by(id) %>% 
  filter(n() >=3) %>%
  ungroup()
n_distinct(cancer_biom$id)
# 10,837 individuals with at least 3 blood tests before their first cancer diagnosis 

# to identify latent class of biomarker trajectories
bios <- c("S_ApoA","S_ApoB","ApoKvot","TC","TG","LDL","HDL","RATHDL","fS_Gluk","S_FAMN","S_Hapt","S_CRP","S_Urat","S_Alb","B_Hb","S_IgG","B_LPKV","S_GLOB")
num_tests <- NULL
bm <- NULL
for (bm in bios) {
  print(bm)
  biom <- NULL
  biom <- cancer_biom %>%
    rename(bio=!!bm) %>%
    filter(!is.na(bio)) %>%
    group_by(id) %>% 
    filter(n() >=3) %>%
    ungroup()
  
  biom <- biom %>%
    mutate(date_samp=as.Date(date_samp),
           can_date=as.Date(can_date),
           time=as.integer(date_samp - can_date)/365.25,
           logbio=log(bio)) %>%
    dplyr::select(id, time, logbio) 
  biom$logbio <- as.numeric(biom$logbio)
  biom$id <- as.numeric(biom$id)
  biom$time <- as.numeric(biom$time)
  # number of patients with at least three measurements of the biomarker
  num_t <- NULL
  num_t <- data.frame(
    biomarkers=bm,
    n=n_distinct(biom$id)
  )
  num_tests <- rbind(num_tests,num_t)
  
  
  #linear fixed 
  lf1 <- hlme(logbio ~ time, subject = "id", ng = 1, data = biom)
  lf2 <- hlme(logbio ~ time, subject = "id", ng = 2, data = biom, mixture = ~time, B =lf1)
  lf3 <- hlme(logbio ~ time, subject = "id", ng = 3, data = biom, mixture = ~time, B =lf1)
  lf4 <- hlme(logbio ~ time, subject = "id", ng = 4, data = biom, mixture = ~time, B =lf1)
  lf5 <- hlme(logbio ~ time, subject = "id", ng = 5, data = biom, mixture = ~time, B =lf1)
  
  # to extract coefficients and P-values
  i <- NULL
  for (i in 1:5) {
    model <- NULL
    table <- NULL
    model <- get(paste0("lf",i))
    table <- data.frame(summary(model))
    saveRDS(model,paste0(".../cancer_infection/",bm,"_lf",i,".RDS"))
    write.csv(table,paste0(".../cancer_infection/",bm,"_pvalues_lf",i,".csv"))
  }
  
  # to extract the statistics of different classes
  table1 <- NULL
  table1 <- summarytable(lf1, lf2, lf3, lf4, lf5, which = c("G", "loglik", "conv", "AIC", "BIC", "%class"))
  table1 <- data.frame(round(table1, 2))
  write.csv(table1,paste0(".../cancer_infection/",bm,"_classes_lf.csv"))
  
  # linear random intercept
  lr1 <- hlme(logbio ~ time, subject = "id", ng = 1, random = ~1, data = biom)
  lr2 <- hlme(logbio ~ time, subject = "id", ng = 2, random = ~1, data = biom, mixture = ~time, B =lr1)
  lr3 <- hlme(logbio ~ time, subject = "id", ng = 3, random = ~1, data = biom, mixture = ~time, B =lr1)
  lr4 <- hlme(logbio ~ time, subject = "id", ng = 4, random = ~1, data = biom, mixture = ~time, B =lr1)
  lr5 <- hlme(logbio ~ time, subject = "id", ng = 5, random = ~1, data = biom, mixture = ~time, B =lr1)
  
  # to extract coefficients and P-values
  i <- NULL
  for (i in 1:5) {
    model <- NULL
    table <- NULL
    model <- get(paste0("lr",i))
    table <- data.frame(summary(model))
    saveRDS(model,paste0(".../cancer_infection/",bm,"_lr",i,".RDS"))
    write.csv(table,paste0(".../cancer_infection/",bm,"_pvalues_lr",i,".csv"))
  }
  
  # to extract the statistics of different classes
  table2 <- NULL
  table2 <- summarytable(lr1, lr2, lr3, lr4, lr5, which = c("G", "loglik", "conv", "AIC", "BIC", "%class"))
  table2 <- data.frame(round(table2, 2))
  write.csv(table2,paste0(".../cancer_infection/",bm,"_classes_lr.csv"))
  
  
  # linear random intercept + slope
  lrs1<-  hlme(logbio ~ time, subject = "id", ng = 1, random = ~1 + time, data = biom)
  lrs2 <- hlme(logbio ~ time, subject = "id", ng = 2, random = ~1 + time, data = biom, mixture = ~time, B =lrs1)
  lrs3 <- hlme(logbio ~ time, subject = "id", ng = 3, random = ~1 + time, data = biom, mixture = ~time, B =lrs1)
  lrs4 <- hlme(logbio ~ time, subject = "id", ng = 4, random = ~1 + time, data = biom, mixture = ~time, B =lrs1)
  lrs5 <- hlme(logbio ~ time, subject = "id", ng = 5, random = ~1 + time, data = biom, mixture = ~time, B =lrs1)
  
  # to extract coefficients and P-values
  i <- NULL
  for (i in 1:5) {
    model <- NULL
    table <- NULL
    model <- get(paste0("lrs",i))
    table <- data.frame(summary(model))
    saveRDS(model,paste0(".../cancer_infection/",bm,"_lrs",i,".RDS"))
    write.csv(table,paste0(".../cancer_infection/",bm,"_pvalues_lrs",i,".csv"))
  }
  
  # to extract the statistics of different classes
  table3 <- NULL
  table3 <- summarytable(lrs1, lrs2, lrs3, lrs4, lrs4, which = c("G", "loglik", "conv", "AIC", "BIC", "%class"))
  table3 <- data.frame(round(table3, 2))
  write.csv(table3,paste0(".../cancer_infection/",bm,"_classes_lrs.csv"))
  
  
  # quadratic fixed effect models
  qf1 <- hlme(logbio ~ time + I(time^2), subject = "id", ng = 1, data = biom)
  qf2 <- hlme(logbio ~ time + I(time^2), subject = "id", ng = 2, data = biom, mixture = ~ time + I(time^2), B =qf1)
  qf3 <- hlme(logbio ~ time + I(time^2), subject = "id", ng = 3, data = biom, mixture = ~ time + I(time^2), B =qf1)
  qf4 <- hlme(logbio ~ time + I(time^2), subject = "id", ng = 4, data = biom, mixture = ~ time + I(time^2), B =qf1)
  qf5 <- hlme(logbio ~ time + I(time^2), subject = "id", ng = 5, data = biom, mixture = ~ time + I(time^2), B =qf1)
  
  # to extract coefficients and P-values
  i <- NULL
  for (i in 1:5) {
    model <- NULL
    table <- NULL
    model <- get(paste0("qf",i))
    table <- data.frame(summary(model))
    saveRDS(model,paste0(".../cancer_infection/",bm,"qf",i,".RDS"))
    write.csv(table,paste0(".../cancer_infection/",bm,"_pvalues_qf",i,".csv"))
  }
  
  # to extract the statistics of different classes
  table4 <- NULL
  table4 <- summarytable(qf1, qf2, qf3, qf4, qf5, which = c("G", "loglik", "conv", "AIC", "BIC", "%class"))
  table4 <- data.frame(round(table4, 2))
  write.csv(table4,paste0(".../cancer_infection/",bm,"_classes_qf.csv"))
  
  
  # quadratic random intercept 
  qr1 <- hlme(logbio ~ time + I(time^2), subject = "id", ng = 1, random = ~1, data = biom)
  qr2 <- hlme(logbio ~ time + I(time^2), subject = "id", ng = 2, random = ~1, data = biom, mixture = ~ time + I(time^2), B =qr1)
  qr3 <- hlme(logbio ~ time + I(time^2), subject = "id", ng = 3, random = ~1, data = biom, mixture = ~ time + I(time^2), B =qr1)
  qr4 <- hlme(logbio ~ time + I(time^2), subject = "id", ng = 4, random = ~1, data = biom, mixture = ~ time + I(time^2), B =qr1)
  qr5 <- hlme(logbio ~ time + I(time^2), subject = "id", ng = 5, random = ~1, data = biom, mixture = ~ time + I(time^2), B =qr1)
  
  # to extract coefficients and P-values
  i <- NULL
  for (i in 1:5) {
    model <- NULL
    table <- NULL
    model <- get(paste0("qr",i))
    table <- data.frame(summary(model))
    saveRDS(model,paste0(".../cancer_infection/",bm,"qr",i,".RDS"))
    write.csv(table,paste0(".../cancer_infection/",bm,"_pvalues_qr",i,".csv"))
  }
  
  # to extract the statistics of different classes
  table5 <- NULL
  table5 <- summarytable(qr1, qr2, qr3, qr4, qr5, which = c("G", "loglik", "conv", "AIC", "BIC", "%class"))
  table5 <- data.frame(round(table5, 2))
  write.csv(table5,paste0(".../cancer_infection/",bm,"_classes_qr.csv"))
  
  # quadratic random intercept + slope 
  qrs1 <- hlme(logbio ~ time + I(time^2), subject = "id", ng = 1, random = ~1 + time, data = biom)
  qrs2 <- hlme(logbio ~ time + I(time^2), subject = "id", ng = 2, random = ~1 + time, data = biom, mixture = ~ time + I(time^2), B =qrs1)
  qrs3 <- hlme(logbio ~ time + I(time^2), subject = "id", ng = 3, random = ~1 + time, data = biom, mixture = ~ time + I(time^2), B =qrs1)
  qrs4 <- hlme(logbio ~ time + I(time^2), subject = "id", ng = 4, random = ~1 + time, data = biom, mixture = ~ time + I(time^2), B =qrs1)
  qrs5 <- hlme(logbio ~ time + I(time^2), subject = "id", ng = 5, random = ~1 + time, data = biom, mixture = ~ time + I(time^2), B =qrs1)
  
  # to extract coefficients and P-values
  i <- NULL
  for (i in 1:5) {
    model <- NULL
    table <- NULL
    model <- get(paste0("qrs",i))
    table <- data.frame(summary(model))
    saveRDS(model,paste0(".../cancer_infection/",bm,"qrs",i,".RDS"))
    write.csv(table,paste0(".../cancer_infection/",bm,"_pvalues_qrs",i,".csv"))
  }
  
  # to extract the statistics of different classes
  table6 <- NULL
  table6 <- summarytable(qrs1, qrs2, qrs3, qrs4, qrs5, which = c("G", "loglik", "conv", "AIC", "BIC", "%class"))
  table6 <- data.frame(round(table6, 2))
  write.csv(table6,paste0(".../cancer_infection/",bm,"_classes_qrs.csv"))
  
}

write.csv(num_tests,".../cancer_infection/lcgm_number_individuals.csv")
#--------------------------------------------------------------------------------


# check the outputs from LCGM manually and determine the best-fit models: the lowest BIC, P-values<0.05 for all polynomial terms, and the group probabilities >=5%
rm(list=ls())
# prepare a supplementary table of parameters of LCGM
bios <- c("fS_Gluk","S_FAMN","TC","TG","LDL","HDL","RATHDL","S_ApoA","S_ApoB","ApoKvot","B_Hb","S_GLOB","S_Hapt","S_CRP","S_IgG","B_LPKV","S_Urat","S_Alb")
models <- c("l","q")
bio <- NULL
mo <- NULL
results <- NULL

for (bio in bios) {
  print(bio)
  for (mo in models) {
    print(mo)
    dd <- NULL
    d1 <- NULL
    d2 <- NULL
    d3 <- NULL
    
    d1 <- read.csv(paste0(".../cancer_infection/",bio,"_classes_",mo,"f.csv"), header = T, sep = ",")
    d2 <- read.csv(paste0(".../cancer_infection/",bio,"_classes_",mo,"r.csv"), header = T, sep = ",")
    d3 <- read.csv(paste0(".../cancer_infection/",bio,"_classes_",mo,"rs.csv"), header = T, sep = ",")
    
    if (!"X.class5" %in% colnames(d1)) {
      d1$X.class5 <- NA
    } 
    if (!"X.class5" %in% colnames(d2)) {
      d2$X.class5 <- NA
    } 
    if (!"X.class5" %in% colnames(d3)) {
      d3$X.class5 <- NA
    } 
    
    d1 <- d1 %>%
      filter(G>1) %>%
      rename(C=G,
             logit=loglik) %>%
      mutate(Biomarkers=bio,
             model=mo,
             Class_membership=case_when(
               C==2 ~ paste0(X.class1,"/",X.class2),
               C==3 ~ paste0(X.class1,"/",X.class2,"/",X.class3),
               C==4 ~ paste0(X.class1,"/",X.class2,"/",X.class3,"/",X.class4),
               C==5 ~ paste0(X.class1,"/",X.class2,"/",X.class3,"/",X.class4,"/",X.class5))) %>%
      dplyr::select(Biomarkers,model,logit,BIC,Class_membership)
    d2 <- d2 %>%
      filter(G>1) %>%
      rename(logit_r=loglik,
             BIC_r=BIC) %>%
      mutate(Class_membership_r=case_when(
        G==2 ~ paste0(X.class1,"/",X.class2),
        G==3 ~ paste0(X.class1,"/",X.class2,"/",X.class3),
        G==4 ~ paste0(X.class1,"/",X.class2,"/",X.class3,"/",X.class4),
        G==5 ~ paste0(X.class1,"/",X.class2,"/",X.class3,"/",X.class4,"/",X.class5))) %>%
      dplyr::select(logit_r,BIC_r,Class_membership_r)
    d3 <- d3 %>%
      filter(G>1) %>%
      rename(logit_rs=loglik,
             BIC_rs=BIC) %>%
      mutate(Class_membership_rs=case_when(
        G==2 ~ paste0(X.class1,"/",X.class2),
        G==3 ~ paste0(X.class1,"/",X.class2,"/",X.class3),
        G==4 ~ paste0(X.class1,"/",X.class2,"/",X.class3,"/",X.class4),
        G==5 ~ paste0(X.class1,"/",X.class2,"/",X.class3,"/",X.class4,"/",X.class5))) %>%
      dplyr::select(logit_rs,BIC_rs,Class_membership_rs)
    
    dd <- cbind(d1,d2,d3)
    
    results <- rbind(results,dd)
  }
}

results <- results %>%
  mutate(Biomarkers=case_when(
    Biomarkers=="fS_Gluk" ~ "Glucose", 
    Biomarkers=="S_FAMN" ~ "Fructosamine", 
    Biomarkers=="TC" ~ "TC", 
    Biomarkers=="TG" ~ "TG", 
    Biomarkers=="LDL" ~ "LDL",
    Biomarkers=="HDL" ~ "HDL",
    Biomarkers=="RATHDL" ~ "LDL/HDL ratio", 
    Biomarkers=="S_ApoA" ~ "ApoA1", 
    Biomarkers=="S_ApoB" ~ "ApoB", 
    Biomarkers=="ApoKvot" ~ "ApoA1/ApoB", 
    Biomarkers=="B_Hb" ~ "Hemoglobin",
    Biomarkers=="S_GLOB" ~ "Globulin",
    Biomarkers=="S_Hapt" ~ "Haptoglobin", 
    Biomarkers=="S_CRP" ~ "C-reactive protein", 
    Biomarkers=="S_IgG" ~ "IgG",
    Biomarkers=="B_LPKV" ~ "Leukocyte", 
    Biomarkers=="S_Urat" ~ "Uric acid", 
    Biomarkers=="S_Alb" ~ "Albumin"
  ),
  model=case_when(
    model=="l" ~ "L",
    model=="q" ~ "Q"
  ))

write.csv(results,".../cancer_infection/export/model_parameters.csv")


# best-fit models:
# S_ApoA: qr2
# S_ApoB: lr2
# ApoKvot: lr2
# TC: lf2
# TG: qrs2
# LDL: lf4
# HDL: qrs2
# RATHDL: lr2
# fS_Gluk: lr2
# S_FAMN: lr2
# S_Hapt: qf3
# S_CRP: qr2
# S_Urat: lf3
# S_Alb: qf3
# B_Hb: qf3
# S_IgG: lr2
# B_LPKV: lf3
# S_GLOB: lf2

# "S_ApoA","S_ApoB","ApoKvot","TC","TG","LDL","HDL","RATHDL","fS_Gluk","S_FAMN","S_Hapt","S_CRP","S_Urat","S_Alb","B_Hb","S_IgG","B_LPKV","S_GLOB"

# identified classes of biomarker trajectories for cancer patients
biom <- read.csv(".../biomarkers.csv")
biom$id <- as.character(biom$id) 
biom <- biom[,2:45]
master <- readRDS(".../master_cancer1985_2005.RDS")

cancer <- master %>%
  dplyr::select(id, can_date)

cancer_biom <- cancer %>%
  left_join(biom,by="id") %>%
  filter(date_samp <= can_date) %>%
  filter(!is.na(S_ApoA) | !is.na(S_ApoB) | !is.na(ApoKvot) | !is.na(TC) | !is.na(TG) | !is.na(LDL) | !is.na(HDL) | !is.na(RATHDL) | !is.na(fS_Gluk) | !is.na(S_FAMN) | !is.na(S_Hapt) | !is.na(S_CRP) | !is.na(S_Urat) | !is.na(S_Alb) | !is.na(B_Hb) | !is.na(S_IgG) | !is.na(B_LPKV) | !is.na(S_GLOB))

cancer_biom <- cancer_biom %>%
  mutate(date_samp=as.Date(date_samp),
         can_date=as.Date(can_date)) %>%
  filter(can_date - date_samp>=0 & can_date - date_samp<=3650)

can_biom <- cancer_biom %>%
  group_by(id) %>% 
  filter(n() >=3) %>%
  ungroup() %>%
  distinct(id)
# 10837

bioms <- can_biom %>%
  left_join(cancer_biom, by="id")

S_ApoA <- readRDS(".../cancer_infection/S_ApoAqr2.RDS")
S_ApoB <- readRDS(".../cancer_infection/S_ApoB_lr2.RDS")
ApoKvot <- readRDS(".../cancer_infection/ApoKvot_lr2.RDS")
TC <- readRDS(".../cancer_infection/TC_lf2.RDS")
TG <- readRDS(".../cancer_infection/TGqrs2.RDS")
LDL <- readRDS(".../cancer_infection/LDL_lf4.RDS")
HDL <- readRDS(".../cancer_infection/HDLqrs2.RDS")
RATHDL <- readRDS(".../cancer_infection/RATHDL_lr2.RDS")
fS_Gluk <- readRDS(".../cancer_infection/fS_Gluk_lr2.RDS")
S_FAMN <- readRDS(".../cancer_infection/S_FAMN_lr2.RDS")
S_Hapt <- readRDS(".../cancer_infection/S_Haptqf3.RDS")
S_CRP <- readRDS(".../cancer_infection/S_CRPqr2.RDS")
S_Urat <- readRDS(".../cancer_infection/S_Urat_lf3.RDS")
S_Alb <- readRDS(".../cancer_infection/S_Albqf3.RDS")
B_Hb <- readRDS(".../cancer_infection/B_Hbqf3.RDS")
S_IgG <- readRDS(".../cancer_infection/S_IgG_lr2.RDS")
B_LPKV <- readRDS(".../cancer_infection/B_LPKV_lf3.RDS")
S_GLOB <- readRDS(".../cancer_infection/S_GLOB_lf2.RDS")

bioms <- bioms %>%
  mutate(date_samp=as.Date(date_samp),
         can_date=as.Date(can_date),
         time=as.integer(date_samp - can_date)/365.25) 
bioms$id <- as.numeric(bioms$id)
bioms$time <- as.numeric(bioms$time)

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(S_ApoA))) %>%
  filter(!is.na(S_ApoA)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_S_ApoA <- bio %>%
  left_join(predictClass(S_ApoA, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(S_ApoA=class) 

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(S_ApoB))) %>%
  filter(!is.na(S_ApoB)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_S_ApoB <- bio %>%
  left_join(predictClass(S_ApoB, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(S_ApoB=class) 

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(ApoKvot))) %>%
  filter(!is.na(ApoKvot)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_ApoKvot <- bio %>%
  left_join(predictClass(ApoKvot, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(ApoKvot=class) 

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(TC))) %>%
  filter(!is.na(TC)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_TC <- bio %>%
  left_join(predictClass(TC, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(TC=class) 

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(TG))) %>%
  filter(!is.na(TG)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_TG <- bio %>%
  left_join(predictClass(TG, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(TG=class) 

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(LDL))) %>%
  filter(!is.na(LDL)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_LDL <- bio %>%
  left_join(predictClass(LDL, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(LDL=class) 

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(HDL))) %>%
  filter(!is.na(HDL)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_HDL <- bio %>%
  left_join(predictClass(HDL, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(HDL=class) 

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(RATHDL))) %>%
  filter(!is.na(RATHDL)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_RATHDL <- bio %>%
  left_join(predictClass(RATHDL, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(RATHDL=class) 

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(fS_Gluk))) %>%
  filter(!is.na(fS_Gluk)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_fS_Gluk <- bio %>%
  left_join(predictClass(fS_Gluk, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(fS_Gluk=class) 

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(S_FAMN))) %>%
  filter(!is.na(S_FAMN)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_S_FAMN <- bio %>%
  left_join(predictClass(S_FAMN, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(S_FAMN=class) 

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(S_Hapt))) %>%
  filter(!is.na(S_Hapt)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_S_Hapt <- bio %>%
  left_join(predictClass(S_Hapt, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(S_Hapt=class) 

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(S_CRP))) %>%
  filter(!is.na(S_CRP)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_S_CRP <- bio %>%
  left_join(predictClass(S_CRP, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(S_CRP=class) 

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(S_Urat))) %>%
  filter(!is.na(S_Urat)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_S_Urat <- bio %>%
  left_join(predictClass(S_Urat, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(S_Urat=class) 

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(S_Alb))) %>%
  filter(!is.na(S_Alb)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_S_Alb <- bio %>%
  left_join(predictClass(S_Alb, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(S_Alb=class) 

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(B_Hb))) %>%
  filter(!is.na(B_Hb)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_B_Hb <- bio %>%
  left_join(predictClass(B_Hb, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(B_Hb=class) 

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(S_IgG))) %>%
  filter(!is.na(S_IgG)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_S_IgG <- bio %>%
  left_join(predictClass(S_IgG, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(S_IgG=class) 

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(B_LPKV))) %>%
  filter(!is.na(B_LPKV)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_B_LPKV <- bio %>%
  left_join(predictClass(B_LPKV, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(B_LPKV=class) 

bio <- NULL
bio <- bioms %>% 
  mutate(logbio=as.numeric(log(S_GLOB))) %>%
  filter(!is.na(S_GLOB)) %>%
  dplyr::select(id,time,logbio)
bio$logbio <- as.numeric(bio$logbio)
c_S_GLOB <- bio %>%
  left_join(predictClass(S_GLOB, newdata=bio), by = "id") %>%
  distinct(id, .keep_all = T) %>%
  dplyr::select(id, class)%>%
  rename(S_GLOB=class) 

master$id <- as.numeric(master$id)
can_biom$id <- as.numeric(can_biom$id)
can_biom <- can_biom %>%
  left_join(c_S_ApoA, by="id") %>%
  left_join(c_S_ApoB, by="id") %>%
  left_join(c_ApoKvot, by="id") %>%
  left_join(c_TC, by="id") %>%
  left_join(c_TG, by="id") %>%
  left_join(c_LDL, by="id") %>%
  left_join(c_HDL, by="id") %>%
  left_join(c_RATHDL, by="id") %>%
  left_join(c_fS_Gluk, by="id") %>%
  left_join(c_S_FAMN, by="id") %>%
  left_join(c_S_Hapt, by="id") %>%
  left_join(c_S_CRP, by="id") %>%
  left_join(c_S_Urat, by="id") %>%
  left_join(c_S_Alb, by="id") %>%
  left_join(c_B_Hb, by="id") %>%
  left_join(c_S_IgG, by="id") %>%
  left_join(c_B_LPKV, by="id") %>%
  left_join(c_S_GLOB, by="id") %>%
  left_join(master, by="id")

saveRDS(can_biom, "W:/C6_AmorisGuest/Dang Wei/Bio_com_cancer/Working data/cancer_bio_trajectories.RDS")
#-------------------------------------------------------------------------------

# to make the figure of biomarker trajectories
pdf(".../cancer_infection/export/figure_bio_trajectories.pdf",
    width = 12, height = 12)
par(mfrow=c(5,4)) 

## 1. glu
data_predict <- data.frame(time=seq(-10,0,length.out=160))
predict <- predictY(fS_Gluk, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red","navy"), lty=1, lwd=2, ylab="Glucose (mmol/L)", bty = 'n', shades = TRUE,
     xlab = "Time To Cancer Diagnosis (Years)", legend=NULL, main="Glucose (n=9890)", ylim=c(2,15))

## 2. S_FAMN
predict <- predictY(S_FAMN, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red","navy"), lty=1, lwd=2, ylab="Fructosamine (mmol/L)", bty = 'n', shades = TRUE,
     xlab = "Time To Cancer Diagnosis (Years)", legend=NULL, main="Fructosamine (n=9353)", ylim=c(0,4))

## 3. TC
predict <- predictY(TC, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red","navy"), lty=1, lwd=2, ylab="TC (mmol/L)", bty = 'n', shades = TRUE,
     xlab = "Time To Cancer Diagnosis (Years)", legend=NULL, main="TC (n=10,370)", ylim=c(2,8))

## 4. TG
predict <- predictY(TG, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red","navy"), lty=1, lwd=2, ylab="Triglycerides (mmol/L)", bty = 'n', shades = TRUE,
     xlab = "Time To Cancer Diagnosis (Years)", legend=NULL, main="Triglycerides (n=10,329)", ylim=c(0,3))

## 5. HDL
predict <- predictY(HDL, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red", "navy"), lty=1, lwd=2, ylab="HDL (mmol/L)", bty = 'n', shades = TRUE,
     xlab = "Time To Cancer Diagnosis (Years)", legend=NULL, main="HDL (n=5555)", ylim=c(0,2))

## 6. LDL
predict <- predictY(LDL, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red","navy","orange","green"), lty=1, lwd=2, ylab="LDL (mmol/L)", bty = 'n', shades = TRUE,
     xlab = "Time To Cancer Diagnosis (Years)", legend=NULL, main="LDL (n=5750)", ylim=c(2,7))

## 7. LDL/HDL
predict <- predictY(RATHDL, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red", "navy"), lty=1, lwd=2, ylab="LDL/HDL ratio", bty = 'n', shades = TRUE,
     xlab = "Time To Cancer Diagnosis (Years)", legend=NULL, main="LDL/HDL ratio (n=5555)", ylim=c(0, 5))

## 8. Apoa
predict <- predictY(S_ApoA, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red", "navy"), lty=1, lwd=2, ylab="ApoA1 (g/L)", bty = 'n', shades = TRUE,
     xlab = "Time To Cancer Diagnosis (Years)", legend=NULL, main="ApoA1 (n=5812)", ylim=c(1, 2))

## 9. apob
predict <- predictY(S_ApoB, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red", "navy"), lty=1, lwd=2, ylab="ApoB (g/L)", bty = 'n', shades = TRUE,
     xlab = "Time To Cancer Diagnosis (Years)", legend=NULL, main="ApoB (n=5623)", ylim=c(0.5, 2.5))

## 10. apo ratio
predict <- predictY(ApoKvot, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red","navy"), lty=1, lwd=2, ylab="ApoB/ApoA1 ratio", bty = 'n', shades = TRUE, 
     xlab= "Time To Cancer Diagnosis (Years)", legend=NULL, main= "ApoB/ApoA1 ratio (n=5419)", ylim=c(0,1.5))

## 11. hemoglobin
predict <- predictY(B_Hb, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red","navy", "orange"), lty=1, lwd=2, ylab="Hemoglobin (g/L)", bty = 'n', shades = TRUE,
     xlab = "Time To Cancer Diagnosis (Years)", 
     legend=NULL, main="Hemoglobin (n=6600)", ylim=c(45,200))

## 12. globulin
predict <- predictY(S_GLOB, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red","navy"), lty=1, lwd=2, ylab="Globulin (g/L)", bty = 'n', shades = TRUE, 
     xlab= "Time To Cancer Diagnosis (Years)", legend=NULL, main= "Globulin (n=5688)", ylim=c(20,40))

## 13. haptoglobin
predict <- predictY(S_Hapt, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red","navy","orange"), lty=1, lwd=2, ylab="Haptoglobin (g/L)", bty = 'n', shades = TRUE, 
     xlab= "Time To Cancer Diagnosis (Years)", legend=NULL, main= "Haptoglobin (n=8506)", ylim=c(0,5))

## 14. CRP
predict <- predictY(S_CRP, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red","navy"), lty=1, lwd=2, ylab="CRP (g/L)", bty = 'n', shades = TRUE,
     xlab = "Time To Cancer Diagnosis (Years)", legend=NULL, main="CRP (n=8586)", ylim=c(1,30))

## 15. IgG
predict <- predictY(S_IgG, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red","navy"), lty=1, lwd=2, ylab="IgG (g/L)", bty = 'n', shades = TRUE,
     xlab = "Time To Cancer Diagnosis (Years)", legend=NULL, main="IgG (n=3482)", ylim=c(1,35))

## 16. Leukocyte
predict <- predictY(B_LPKV, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red","navy","orange"), lty=1, lwd=2, ylab="Leukocyte (10^9/L)", bty = 'n', shades = TRUE,
     xlab = "Time To Cancer Diagnosis (Years)", legend=NULL, main="Leukocyte (n=6447)", ylim=c(1,20))

## 17. uric
predict <- predictY(S_Urat, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red", "navy","orange"), lty=1, lwd=2, ylab="Uric acid (µmol/L)", bty = 'n', shades = TRUE,
     xlab = "Time To Cancer Diagnosis (Years)", legend=NULL, main="Uric acid (n=10,109)", ylim=c(100,450))

## 18.albumin
predict <- predictY(S_Alb, data_predict, var.time = "time",  draws = TRUE)
predict[["pred"]] <- exp(predict[["pred"]])
plot(predict, col=c("red","navy", "orange"), lty=1, lwd=2, ylab="Albumin (g/L)", bty = 'n', shades = TRUE,
     xlab = "Time To Cancer Diagnosis (Years)", 
     legend=NULL, main="Albumin (n=10,403)", ylim=c(30,50))

## 19.Empty
plot(1, type="n", axes=FALSE, ann=FALSE)


## 20. Legend
plot(1, type="n", axes=FALSE, ann=FALSE)
# Add the legend
legend("center", 
       legend=c("Class 1", "Class 2", "Class 3", "Class 4"), 
       col=c("red", "navy", "orange","green"), 
       lty=1, lwd=3, cex=1.5, y.intersp=0.8, bty= "n")

dev.off() #turn off develop


#---------------------------------------------------------------------------------


# Step 4: evaluate whether the associations between biomarker trajectories and the risk of infection differed over time using flexible parametric models --> 1. determine the number of knots using crude models, 2. fit multivariable models, 3. plot hazards by groups over time, 4. plot hazard ratios over time
library(biostat3)  
library(rstpm2)  
rm(list = ls())
master <- readRDS("W:/C6_AmorisGuest/Dang Wei/Bio_com_cancer/Working data/cancer_bio_trajectories.RDS")

master <- master %>%
  mutate(can_yearc=ifelse(can_year>1984 & can_year<1996, 0, 1),
         infec_both_bef=ifelse(is.na(infec_both_bef),0,infec_both_bef),
         infec_p_bef=ifelse(is.na(infec_p_bef),0,infec_p_bef),
         sepsis_bef=ifelse(is.na(sepsis_bef),0,sepsis_bef),
         lastDate=ifelse(lastDate<can_date, as.Date("2020-12-31"),lastDate),
         lastDate=as.Date(lastDate),
         end_infec_both=as.Date(pmin(lastDate,date_hos_both,na.rm = TRUE)),
         end_infec_p=as.Date(pmin(lastDate,date_hos_pri,na.rm = TRUE)),
         end_sepsis=as.Date(pmin(lastDate,date_sepsis,na.rm = TRUE)),
         age_end_infec_both=as.integer(end_infec_both - birthdate +1)/365.25,
         age_end_infec_p=as.integer(end_infec_p - birthdate+1)/365.25,
         age_end_sepsis=as.integer(end_sepsis - birthdate+1)/365.25,
         age_canc=ifelse(age_can>0 & age_can<36,0,
                         ifelse(age_can>35 & age_can<61, 1, 2)),
         time_infec_both=age_end_infec_both-age_can,
         time_infec_p=age_end_infec_p-age_can,
         time_sepsis=age_end_sepsis-age_can) %>%
  mutate(education_atenroll=ifelse(education_atenroll=="99",NA,education_atenroll),
         income_atenroll=ifelse(income_atenroll=="99",NA,income_atenroll),
         employ_atenroll=ifelse(employ_atenroll=="99",NA,employ_atenroll))

# To determine the number of knots for each biomarkers: lowest BIC
bios <- c("fS_Gluk","S_FAMN","TC","TG","LDL","HDL","RATHDL","S_ApoA","S_ApoB","ApoKvot","B_Hb","S_GLOB","S_Hapt","S_CRP","S_IgG","B_LPKV","S_Urat","S_Alb")
bm <- NULL

for (bm in bios) {
  print(bm)
  toanalysis <- NULL
  toanalysis <- master %>%
    rename(bio=!!bm) %>%
    filter(!is.na(bio))
  toanalysis$bio <- factor(toanalysis$bio)
  Rbind <- function(object) do.call(rbind,object)
  out <- lapply(1:6, function(i) {
    fitaic <- stpm2(Surv(time_infec_p, infec_p) ~ bio, data=toanalysis, df=i)
    data.frame(
      i,
      AIC=AIC(fitaic),
      BIC=BIC(fitaic),
      beta=as.numeric(coef(fitaic)[2]),
      se=coef(summary(fitaic))[2,2],
      p=coef(summary(fitaic))[2,4])
  })
  out %>% Rbind
  print(out)
}

# degree of freedom
# "fS_Gluk": 3
# "S_FAMN": 3
# "TC": 3
# "TG": 3
# "LDL": 3
# "HDL": 3
# "RATHDL": 3
# "S_ApoA": 3
# "S_ApoB": 3
# "ApoKvot": 3
# "B_Hb": 3
# "S_GLOB": 3
# "S_Hapt": 3
# "S_CRP": 3
# "S_IgG": 3
# "B_LPKV": 3
# "S_Urat": 3
# "S_Alb": 3
#-------------------------------------------------------------------------------

# to plot the hazards acorss letant classes over time using crude models
# col=c("red", "navy", "orange","green")
pdf(".../cancer_infection/export/figure_hazard_infec_p.pdf",
    width = 12, height = 12)
par(mfrow=c(5,4)) 

# fS_Gluk
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(fS_Gluk), data=master, df=3)
plot(crudem,newdata=data.frame(fS_Gluk=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Glucose", ylim=c(0,0.20))
lines(crudem,newdata=data.frame(fS_Gluk=1), type="haz", col="red")
lines(crudem,newdata=data.frame(fS_Gluk=2), type="haz", ci=TRUE, col="navy")

# S_FAMN
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(S_FAMN), data=master, df=3)
plot(crudem,newdata=data.frame(S_FAMN=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Fructosamine", ylim=c(0,0.20))
lines(crudem,newdata=data.frame(S_FAMN=1), type="haz", col="red")
lines(crudem,newdata=data.frame(S_FAMN=2), type="haz", ci=TRUE, col="navy")

# TC
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(TC), data=master, df=3)
plot(crudem,newdata=data.frame(TC=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="TC", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(TC=2), type="haz", ci=TRUE, col="navy")
lines(crudem,newdata=data.frame(TC=1), type="haz", col="red")

# TG
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(TG), data=master, df=3)
plot(crudem,newdata=data.frame(TG=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Triglycerides", ylim=c(0,0.15))
lines(crudem,newdata=data.frame(TG=1), type="haz", col="red")
lines(crudem,newdata=data.frame(TG=2), type="haz", ci=TRUE, col="navy")

# HDL
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(HDL), data=master, df=3)
plot(crudem,newdata=data.frame(HDL=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="HDL", ylim=c(0,0.15))
lines(crudem,newdata=data.frame(HDL=1), type="haz", col="red")
lines(crudem,newdata=data.frame(HDL=2), type="haz", ci=TRUE, col="navy")

# LDL
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(LDL), data=master, df=3)
plot(crudem,newdata=data.frame(LDL=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="LDL", ylim=c(0.02,0.10))
lines(crudem,newdata=data.frame(LDL=2), type="haz", ci=TRUE, col="navy") 
lines(crudem,newdata=data.frame(LDL=1), type="haz", col="red")
lines(crudem,newdata=data.frame(LDL=3), type="haz", col="orange")
lines(crudem,newdata=data.frame(LDL=4), type="haz", col="green")

# RATHDL
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(RATHDL), data=master, df=3)
plot(crudem,newdata=data.frame(RATHDL=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="LDL/HDL ratio", ylim=c(0.02,0.15))
lines(crudem,newdata=data.frame(RATHDL=1), type="haz", col="red")
lines(crudem,newdata=data.frame(RATHDL=2), type="haz", ci=TRUE, col="navy") 

# S_ApoA
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(S_ApoA), data=master, df=3)
plot(crudem,newdata=data.frame(S_ApoA=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="ApoA1", ylim=c(0.02,0.10))
lines(crudem,newdata=data.frame(S_ApoA=1), type="haz", col="red")
lines(crudem,newdata=data.frame(S_ApoA=2), type="haz", col="navy") 

# S_ApoB
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(S_ApoB), data=master, df=3)
plot(crudem,newdata=data.frame(S_ApoB=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="ApoB", ylim=c(0.02,0.15))
lines(crudem,newdata=data.frame(S_ApoB=1), type="haz", col="red")
lines(crudem,newdata=data.frame(S_ApoB=2), type="haz", ci=TRUE, col="navy") 

# ApoKvot
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(ApoKvot), data=master, df=3)
plot(crudem,newdata=data.frame(ApoKvot=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="ApoB/ApoA1 ratio", ylim=c(0.02,0.15))
lines(crudem,newdata=data.frame(ApoKvot=1), type="haz", col="red")
lines(crudem,newdata=data.frame(ApoKvot=2), type="haz", ci=TRUE, col="navy") 

# B_Hb
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(B_Hb), data=master, df=3)
plot(crudem,newdata=data.frame(B_Hb=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Hemoglobin", ylim=c(0.02,0.12))
lines(crudem,newdata=data.frame(B_Hb=3), type="haz", ci=TRUE, col="orange") 
lines(crudem,newdata=data.frame(B_Hb=1), type="haz", col="red")
lines(crudem,newdata=data.frame(B_Hb=2), type="haz", col="navy") 

# S_GLOB
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(S_GLOB), data=master, df=3)
plot(crudem,newdata=data.frame(S_GLOB=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Globulin", ylim=c(0.02,0.10))
lines(crudem,newdata=data.frame(S_GLOB=2), type="haz", ci=TRUE, col="navy") 
lines(crudem,newdata=data.frame(S_GLOB=1), type="haz", col="red")

# S_Hapt
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(S_Hapt), data=master, df=3)
plot(crudem,newdata=data.frame(S_Hapt=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Haptoglobin", ylim=c(0.02,0.10))
lines(crudem,newdata=data.frame(S_Hapt=2), type="haz", ci=TRUE, col="navy") 
lines(crudem,newdata=data.frame(S_Hapt=3), type="haz", ci=TRUE, col="orange") 
lines(crudem,newdata=data.frame(S_Hapt=1), type="haz", col="red")

# S_CRP
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(S_CRP), data=master, df=3)
plot(crudem,newdata=data.frame(S_CRP=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="CRP", ylim=c(0.02,0.15))
lines(crudem,newdata=data.frame(S_CRP=1), type="haz", col="red")
lines(crudem,newdata=data.frame(S_CRP=2), type="haz", ci=TRUE, col="navy") 

# S_IgG
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(S_IgG), data=master, df=3)
plot(crudem,newdata=data.frame(S_IgG=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="IgG", ylim=c(0.02,0.20))
lines(crudem,newdata=data.frame(S_IgG=1), type="haz", col="red")
lines(crudem,newdata=data.frame(S_IgG=2), type="haz", ci=TRUE, col="navy") 

# B_LPKV
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(B_LPKV), data=master, df=3)
plot(crudem,newdata=data.frame(B_LPKV=3), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Leukocyte", ylim=c(0.02,0.12))
lines(crudem,newdata=data.frame(B_LPKV=2), type="haz", ci=TRUE, col="navy") 
lines(crudem,newdata=data.frame(B_LPKV=1), type="haz", ci=TRUE, col="red") 
lines(crudem,newdata=data.frame(B_LPKV=3), type="haz", col="orange")

# S_Urat
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(S_Urat), data=master, df=3)
plot(crudem,newdata=data.frame(S_Urat=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Uric acid", ylim=c(0.02,0.12))
lines(crudem,newdata=data.frame(S_Urat=1), type="haz", col="red")
lines(crudem,newdata=data.frame(S_Urat=2), type="haz", ci=TRUE, col="navy") 
lines(crudem,newdata=data.frame(S_Urat=3), type="haz", ci=TRUE, col="orange") 

# S_Alb
crudem <- stpm2(Surv(time_infec_p, infec_p) ~ factor(S_Alb), data=master, df=3)
plot(crudem,newdata=data.frame(S_Alb=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Albumin", ylim=c(0.02,0.12))
lines(crudem,newdata=data.frame(S_Alb=1), type="haz", col="red")
lines(crudem,newdata=data.frame(S_Alb=2), type="haz", ci=TRUE, col="navy") 
lines(crudem,newdata=data.frame(S_Alb=3), type="haz", ci=TRUE, col="orange") 
lines(crudem,newdata=data.frame(S_Alb=2), type="haz", col="navy") 

plot(1, type="n", axes=FALSE, ann=FALSE)

# add an empty plot
plot(1, type="n", axes=FALSE, ann=FALSE)
# Add the legend
legend("center", 
       legend=c("Class 1", "Class 2", "Class 3", "Class 4"), 
       col=c("red", "navy", "orange", "green"), 
       lty=1, lwd=3, cex=1.5, y.intersp=0.8, bty= "n")

dev.off() #turn off develop


# plot hazards for sepsis
pdf(".../cancer_infection/export/figure_hazard_sepsis.pdf",
    width = 12, height = 12)
par(mfrow=c(5,4)) 

# fS_Gluk
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(fS_Gluk), data=master, df=3)
plot(crudem,newdata=data.frame(fS_Gluk=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Glucose", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(fS_Gluk=1), type="haz", col="red")
lines(crudem,newdata=data.frame(fS_Gluk=2), type="haz", ci=TRUE, col="navy")

# S_FAMN
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(S_FAMN), data=master, df=3)
plot(crudem,newdata=data.frame(S_FAMN=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Fructosamine", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(S_FAMN=1), type="haz", col="red")
lines(crudem,newdata=data.frame(S_FAMN=2), type="haz", ci=TRUE, col="navy")

# TC
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(TC), data=master, df=3)
plot(crudem,newdata=data.frame(TC=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="TC", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(TC=1), type="haz", col="red")
lines(crudem,newdata=data.frame(TC=2), type="haz", col="navy")

# TG
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(TG), data=master, df=3)
plot(crudem,newdata=data.frame(TG=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Triglycerides", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(TG=1), type="haz", col="red")
lines(crudem,newdata=data.frame(TG=2), type="haz", ci=TRUE, col="navy")

# HDL
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(HDL), data=master, df=3)
plot(crudem,newdata=data.frame(HDL=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="HDL", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(HDL=1), type="haz", col="red")
lines(crudem,newdata=data.frame(HDL=2), type="haz", ci=TRUE, col="navy")

# LDL
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(LDL), data=master, df=3)
plot(crudem,newdata=data.frame(LDL=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="LDL", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(LDL=1), type="haz", col="red")
lines(crudem,newdata=data.frame(LDL=4), type="haz", col="green")
lines(crudem,newdata=data.frame(LDL=3), type="haz", col="orange")
lines(crudem,newdata=data.frame(LDL=2), type="haz", col="navy") 

# RATHDL
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(RATHDL), data=master, df=3)
plot(crudem,newdata=data.frame(RATHDL=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="LDL/HDL ratio", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(RATHDL=1), type="haz", col="red")
lines(crudem,newdata=data.frame(RATHDL=2), type="haz", col="navy") 

# S_ApoA
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(S_ApoA), data=master, df=3)
plot(crudem,newdata=data.frame(S_ApoA=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="ApoA1", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(S_ApoA=1), type="haz", col="red")
lines(crudem,newdata=data.frame(S_ApoA=2), type="haz", col="navy") 

# S_ApoB
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(S_ApoB), data=master, df=3)
plot(crudem,newdata=data.frame(S_ApoB=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="ApoB", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(S_ApoB=1), type="haz", col="red")
lines(crudem,newdata=data.frame(S_ApoB=2), type="haz", ci=TRUE, col="navy") 

# ApoKvot
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(ApoKvot), data=master, df=3)
plot(crudem,newdata=data.frame(ApoKvot=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="ApoB/ApoA1 ratio", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(ApoKvot=1), type="haz", col="red")
lines(crudem,newdata=data.frame(ApoKvot=2), type="haz", col="navy") 

# B_Hb
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(B_Hb), data=master, df=3)
plot(crudem,newdata=data.frame(B_Hb=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Hemoglobin", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(B_Hb=1), type="haz", col="red")
lines(crudem,newdata=data.frame(B_Hb=2), type="haz", col="navy") 
lines(crudem,newdata=data.frame(B_Hb=3), type="haz", col="orange") 

# S_GLOB
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(S_GLOB), data=master, df=3)
plot(crudem,newdata=data.frame(S_GLOB=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Globulin", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(S_GLOB=1), type="haz", col="red")
lines(crudem,newdata=data.frame(S_GLOB=2), type="haz", ci=TRUE, col="navy") 

# S_Hapt
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(S_Hapt), data=master, df=3)
plot(crudem,newdata=data.frame(S_Hapt=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Haptoglobin", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(S_Hapt=1), type="haz", col="red")
lines(crudem,newdata=data.frame(S_Hapt=2), type="haz", ci=TRUE, col="navy") 
lines(crudem,newdata=data.frame(S_Hapt=3), type="haz", col="orange") 

# S_CRP
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(S_CRP), data=master, df=3)
plot(crudem,newdata=data.frame(S_CRP=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="CRP", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(S_CRP=1), type="haz", col="red")
lines(crudem,newdata=data.frame(S_CRP=2), type="haz", ci=TRUE, col="navy") 

# S_IgG
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(S_IgG), data=master, df=3)
plot(crudem,newdata=data.frame(S_IgG=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="IgG", ylim=c(0,0.15))
lines(crudem,newdata=data.frame(S_IgG=1), type="haz", col="red")
lines(crudem,newdata=data.frame(S_IgG=2), type="haz", ci=TRUE, col="navy") 

# B_LPKV
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(B_LPKV), data=master, df=3)
plot(crudem,newdata=data.frame(B_LPKV=3), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Leukocyte", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(B_LPKV=3), type="haz", col="orange")
lines(crudem,newdata=data.frame(B_LPKV=2), type="haz", ci=TRUE, col="navy") 
lines(crudem,newdata=data.frame(B_LPKV=1), type="haz", col="red") 

# S_Urat
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(S_Urat), data=master, df=3)
plot(crudem,newdata=data.frame(S_Urat=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Uric acid", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(S_Urat=1), type="haz", col="red")
lines(crudem,newdata=data.frame(S_Urat=2), type="haz", ci=TRUE, col="navy") 
lines(crudem,newdata=data.frame(S_Urat=3), type="haz", ci=TRUE, col="orange") 

# S_Alb
crudem <- stpm2(Surv(time_sepsis, sepsis) ~ factor(S_Alb), data=master, df=3)
plot(crudem,newdata=data.frame(S_Alb=1), type="haz",
     xlab="Time Since Cancer diagnosis (years)", ylab="Hazard (per person-year)", main="Albumin", ylim=c(0,0.10))
lines(crudem,newdata=data.frame(S_Alb=1), type="haz", col="red")
lines(crudem,newdata=data.frame(S_Alb=2), type="haz", ci=TRUE, col="navy") 
lines(crudem,newdata=data.frame(S_Alb=3), type="haz", col="orange") 

plot(1, type="n", axes=FALSE, ann=FALSE)

# add an empty plot
plot(1, type="n", axes=FALSE, ann=FALSE)
# Add the legend
legend("center", 
       legend=c("Class 1", "Class 2", "Class 3", "Class 4"), 
       col=c("red", "navy", "orange", "green"), 
       lty=1, lwd=3, cex=1.5, y.intersp=0.8, bty= "n")

dev.off() #turn off develop
#-------------------------------------------------------------------------------

# to fit the multivariable models
master$sex <- factor(master$sex)
master$age_canc <- factor(master$age_canc)
master$can_yearc <- factor(master$can_yearc)
master$country <- factor(master$country)  
master$education_atenroll <- factor(master$education_atenroll)
master$employ_atenroll <- factor(master$employ_atenroll) 
master$income_atenroll <- factor(master$income_atenroll)
master$diabete_bef_cancer <- factor(master$diabete_bef_cancer)
master$cvd_bef_cancer <- factor(master$cvd_bef_cancer)
master$psychi_bef_cancer <- factor(master$psychi_bef_cancer)
master$infec_p_bef <- factor(master$infec_p_bef)

# fS_Gluk
table(master$fS_Gluk) # class 1 as reference
toanalysis <- master %>%
  filter(!is.na(fS_Gluk)) %>%
  mutate(exposure2=ifelse(fS_Gluk==1,0,1))
fS_Gluk_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
fS_Gluk_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
fS_Gluk_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)

# S_FAMN
table(master$S_FAMN) # class 1 as reference
toanalysis <- master %>%
  filter(!is.na(S_FAMN)) %>%
  mutate(exposure2=ifelse(S_FAMN==1,0,1))
S_FAMN_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
S_FAMN_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
S_FAMN_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)

# TC
table(master$TC) # class 2 as reference
toanalysis <- master %>%
  filter(!is.na(TC)) %>%
  mutate(exposure1=ifelse(TC==2,0,1))
TC_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
TC_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
TC_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)

# TG
table(master$TG) # class 2 as reference
toanalysis <- master %>%
  filter(!is.na(TG)) %>%
  mutate(exposure1=ifelse(TG==2,0,1))
TG_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
TG_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
TG_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)

# HDL
table(master$HDL) # class 1 as reference
toanalysis <- master %>%
  filter(!is.na(HDL)) %>%
  mutate(exposure2=ifelse(HDL==1,0,1))
HDL_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
HDL_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
HDL_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)

# LDL
table(master$LDL) # there are four classes: class 3 as reference
toanalysis <- master %>%
  filter(!is.na(LDL)) %>%
  mutate(exposure0=(LDL==3)+0,
         exposure1=(LDL==1)+0,
         exposure2=(LDL==2)+0,
         exposure4=(LDL==4)+0)
LDL_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure2 + exposure4 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2, exposure4 = 2), df=3)
LDL_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure2 + exposure4 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2, exposure4 = 2), df=3)
LDL_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure2 + exposure4 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2, exposure4 = 2), df=3)

# RATHDL
table(master$RATHDL) # class 2 as reference
toanalysis <- master %>%
  filter(!is.na(RATHDL)) %>%
  mutate(exposure1=ifelse(RATHDL==2,0,1))
RATHDL_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
RATHDL_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
RATHDL_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)

# S_ApoA
table(master$S_ApoA) # class 1 as reference
toanalysis <- master %>%
  filter(!is.na(S_ApoA)) %>%
  mutate(exposure2=ifelse(S_ApoA==1,0,1))
S_ApoA_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
S_ApoA_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
S_ApoA_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)

# S_ApoB
table(master$S_ApoB) # class 2 as reference
toanalysis <- master %>%
  filter(!is.na(S_ApoB)) %>%
  mutate(exposure1=ifelse(S_ApoB==2,0,1))
S_ApoB_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
S_ApoB_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
S_ApoB_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)

# ApoKvot
table(master$ApoKvot) # class 2 as reference
toanalysis <- master %>%
  filter(!is.na(ApoKvot)) %>%
  mutate(exposure1=ifelse(ApoKvot==2,0,1))
ApoKvot_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
ApoKvot_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
ApoKvot_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)

# B_Hb
table(master$B_Hb) # there are three classes: class 2 as reference
toanalysis <- master %>%
  filter(!is.na(B_Hb)) %>%
  mutate(exposure0=(B_Hb==2)+0,
         exposure1=(B_Hb==1)+0,
         exposure3=(B_Hb==3)+0)
B_Hb_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
B_Hb_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
B_Hb_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)

# S_GLOB
table(master$S_GLOB) # class 1 as reference
toanalysis <- master %>%
  filter(!is.na(S_GLOB)) %>%
  mutate(exposure2=ifelse(S_GLOB==1,0,1))
S_GLOB_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
S_GLOB_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
S_GLOB_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)

# S_Hapt
table(master$S_Hapt) # there are three classes: class 3 as reference
toanalysis <- master %>%
  filter(!is.na(S_Hapt)) %>%
  mutate(exposure0=(S_Hapt==3)+0,
         exposure1=(S_Hapt==1)+0,
         exposure2=(S_Hapt==2)+0)
S_Hapt_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
S_Hapt_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
S_Hapt_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)

# S_CRP
table(master$S_CRP) # class 2 as reference
toanalysis <- master %>%
  filter(!is.na(S_CRP)) %>%
  mutate(exposure1=ifelse(S_CRP==2,0,1))
S_CRP_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
S_CRP_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
S_CRP_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)

# S_IgG
table(master$S_IgG) # class 1 as reference
toanalysis <- master %>%
  filter(!is.na(S_IgG)) %>%
  mutate(exposure2=ifelse(S_IgG==1,0,1))
S_IgG_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
S_IgG_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
S_IgG_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)

# B_LPKV
table(master$B_LPKV) # there are three classes: class 1 as reference
toanalysis <- master %>%
  filter(!is.na(B_LPKV)) %>%
  mutate(exposure0=(B_LPKV==1)+0,
         exposure2=(B_LPKV==2)+0,
         exposure3=(B_LPKV==3)+0)
B_LPKV_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure2 = 2, exposure3 = 2), df=3)
B_LPKV_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure2 = 2, exposure3 = 2), df=3)
B_LPKV_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure2 = 2, exposure3 = 2), df=3)

# S_Urat
table(master$S_Urat) # there are three classes: class 2 as reference
toanalysis <- master %>%
  filter(!is.na(S_Urat)) %>%
  mutate(exposure0=(S_Urat==2)+0,
         exposure1=(S_Urat==1)+0,
         exposure3=(S_Urat==3)+0)
S_Urat_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
S_Urat_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
S_Urat_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)

# S_Alb
table(master$S_Alb) # there are three classes: class 3 as reference
toanalysis <- master %>%
  filter(!is.na(S_Alb)) %>%
  mutate(exposure0=(S_Alb==3)+0,
         exposure1=(S_Alb==1)+0,
         exposure2=(S_Alb==2)+0)
S_Alb_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
S_Alb_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
S_Alb_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)


# to plot the hazards across latent classes over time
newdata0 <- data.frame(exposure1=0,sex="1",age_can=0,country="0",can_yearc="0",education_atenroll="0",income_atenroll="0",employ_atenroll="0")
newdata1 <- data.frame(exposure2=0,sex="1",age_can=0,country="0",can_yearc="0",education_atenroll="0",income_atenroll="0",employ_atenroll="0")
newdata2 <- data.frame(exposure1=0,exposure2=0,exposure4=0,sex="1",age_can=0,country="0",can_yearc="0",education_atenroll="0",income_atenroll="0",employ_atenroll="0")
newdata3 <- data.frame(exposure1=0,exposure3=0,sex="1",age_can=0,country="0",can_yearc="0",education_atenroll="0",income_atenroll="0",employ_atenroll="0")
newdata4 <- data.frame(exposure1=0,exposure2=0,sex="1",age_can=0,country="0",can_yearc="0",education_atenroll="0",income_atenroll="0",employ_atenroll="0")
newdata5 <- data.frame(exposure2=0,exposure3=0,sex="1",age_can=0,country="0",can_yearc="0",education_atenroll="0",income_atenroll="0",employ_atenroll="0")

# hospitalizations for infection
pdf(".../cancer_infection/export/figure_hr_infec_p.pdf",
    width = 12, height = 12)
par(mfrow=c(5,4)) 
# fS_Gluk
plot(fS_Gluk_fle1,newdata=newdata1, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure2", main="Glucose", ci=TRUE, ylim=c(0,3))
lines(fS_Gluk_fle1,newdata=newdata1, var="exposure2", type="hr", col="navy")
abline(h=1, col="black", lty=2)
# S_FAMN
plot(S_FAMN_fle1,newdata=newdata1, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure2", main="Fructosamine", ci=TRUE, ylim=c(0,3))
lines(S_FAMN_fle1,newdata=newdata1, var="exposure2", type="hr", col="navy")
abline(h=1, col="black", lty=2)
# TC
plot(TC_fle1,newdata=newdata0, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="TC", ci=TRUE, ylim=c(0,3))
lines(TC_fle1,newdata=newdata0, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)
# TG
plot(TG_fle1,newdata=newdata0, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="Triglycerides", ci=TRUE, ylim=c(0,3))
lines(TG_fle1,newdata=newdata0, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)
# HDL
plot(HDL_fle1,newdata=newdata1, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure2", main="HDL", ci=TRUE, ylim=c(0,3))
lines(HDL_fle1,newdata=newdata1, var="exposure2", type="hr", col="navy")
abline(h=1, col="black", lty=2)
# LDL
plot(LDL_fle1,newdata=newdata2, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="LDL", ci=TRUE, ylim=c(0,3))
lines(LDL_fle1,newdata=newdata2, var="exposure2", type="hr", ci=TRUE, col="navy")
lines(LDL_fle1,newdata=newdata2, var="exposure4", type="hr", ci=TRUE, col="green")
lines(LDL_fle1,newdata=newdata2, var="exposure1", type="hr", col="red")
lines(LDL_fle1,newdata=newdata2, var="exposure2", type="hr", col="navy")
abline(h=1, col="black", lty=2)
# RATHDL
plot(RATHDL_fle1,newdata=newdata0, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="LDL/HDL ratio", ci=TRUE, ylim=c(0,3))
lines(RATHDL_fle1,newdata=newdata0, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)
# S_ApoA
plot(S_ApoA_fle1,newdata=newdata1, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure2", main="ApoA1", ci=TRUE, ylim=c(0,3))
lines(S_ApoA_fle1,newdata=newdata1, var="exposure2", type="hr", col="navy")
abline(h=1, col="black", lty=2)
# S_ApoB
plot(S_ApoB_fle1,newdata=newdata0, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="ApoB", ci=TRUE, ylim=c(0,3))
lines(S_ApoB_fle1,newdata=newdata0, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)
# ApoKvot
plot(ApoKvot_fle1,newdata=newdata0, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="ApoB/ApoA1 ratio", ci=TRUE, ylim=c(0,3))
lines(ApoKvot_fle1,newdata=newdata0, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)
# B_Hb
plot(B_Hb_fle1,newdata=newdata3, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="Hemoglobin", ci=TRUE, ylim=c(0,3))
lines(B_Hb_fle1,newdata=newdata3, var="exposure3", type="hr", ci=TRUE, col="orange")
lines(B_Hb_fle1,newdata=newdata3, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)
# S_GLOB
plot(S_GLOB_fle1,newdata=newdata1, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure2", main="Globulin", ci=TRUE, ylim=c(0,3))
lines(S_GLOB_fle1,newdata=newdata1, var="exposure2", type="hr", col="navy")
abline(h=1, col="black", lty=2)
# S_Hapt
plot(S_Hapt_fle1,newdata=newdata4, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="Haptoglobin", ci=TRUE, ylim=c(0,3))
lines(S_Hapt_fle1,newdata=newdata4, var="exposure2", type="hr", ci=TRUE, col="navy")
lines(S_Hapt_fle1,newdata=newdata4, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)
# S_CRP
plot(S_CRP_fle1,newdata=newdata0, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="CRP", ci=TRUE, ylim=c(0,3))
lines(S_CRP_fle1,newdata=newdata0, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)
# S_IgG
plot(S_IgG_fle1,newdata=newdata1, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure2", main="IgG", ci=TRUE, ylim=c(0,3))
lines(S_IgG_fle1,newdata=newdata1, var="exposure2", type="hr", col="navy")
abline(h=1, col="black", lty=2)
# B_LPKV
plot(B_LPKV_fle1,newdata=newdata5, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure3", main="Leukocyte", ci=TRUE, ylim=c(0,3))
lines(B_LPKV_fle1,newdata=newdata5, var="exposure2", type="hr", ci=TRUE, col="navy")
lines(B_LPKV_fle1,newdata=newdata5, var="exposure3", type="hr", col="orange")
abline(h=1, col="black", lty=2)
# S_Urat
plot(S_Urat_fle1,newdata=newdata3, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="Uric acid", ci=TRUE, ylim=c(0,3))
lines(S_Urat_fle1,newdata=newdata3, var="exposure1", type="hr", col="red")
lines(S_Urat_fle1,newdata=newdata3, var="exposure3", type="hr", ci=TRUE, col="orange")
abline(h=1, col="black", lty=2)
# S_Alb
plot(S_Alb_fle1,newdata=newdata4, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="Albumin", ci=TRUE, ylim=c(0,3))
lines(S_Alb_fle1,newdata=newdata4, var="exposure2", type="hr", ci=TRUE, col="navy")
lines(S_Alb_fle1,newdata=newdata4, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)

# add an empty plot
plot(1, type="n", axes=FALSE, ann=FALSE)

# Add the legend
plot(1, type="n", axes=FALSE, ann=FALSE)
legend("center", 
       legend=c("Class 1 VS. Reference", "Class 2 VS. Reference", "Class 3 VS. Reference", "Class 4 VS. Reference"), 
       col=c("red", "navy", "orange", "green"), 
       lty=1, lwd=3, cex=1.5, y.intersp=0.8, bty= "n")

dev.off() #turn off develop

# sepsis
pdf(".../cancer_infection/export/figure_hr_sepsis.pdf",
    width = 12, height = 12)
par(mfrow=c(5,4)) 
# fS_Gluk
plot(fS_Gluk_fle3,newdata=newdata1, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure2", main="Glucose", ci=TRUE, ylim=c(0,3))
lines(fS_Gluk_fle3,newdata=newdata1, var="exposure2", type="hr", col="navy")
abline(h=1, col="black", lty=2)
# S_FAMN
plot(S_FAMN_fle3,newdata=newdata1, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure2", main="Fructosamine", ci=TRUE, ylim=c(0,3))
lines(S_FAMN_fle3,newdata=newdata1, var="exposure2", type="hr", col="navy")
abline(h=1, col="black", lty=2)
# TC
plot(TC_fle3,newdata=newdata0, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="TC", ci=TRUE, ylim=c(0,3))
lines(TC_fle3,newdata=newdata0, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)
# TG
plot(TG_fle3,newdata=newdata0, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="Triglycerides", ci=TRUE, ylim=c(0,3))
lines(TG_fle3,newdata=newdata0, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)
# HDL
plot(HDL_fle3,newdata=newdata1, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure2", main="HDL", ci=TRUE, ylim=c(0,3))
lines(HDL_fle3,newdata=newdata1, var="exposure2", type="hr", col="navy")
abline(h=1, col="black", lty=2)
# LDL
plot(LDL_fle3,newdata=newdata2, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="LDL", ci=TRUE, ylim=c(0,3))
lines(LDL_fle3,newdata=newdata2, var="exposure2", type="hr", ci=TRUE, col="navy")
lines(LDL_fle3,newdata=newdata2, var="exposure4", type="hr", ci=TRUE, col="green")
lines(LDL_fle3,newdata=newdata2, var="exposure1", type="hr", col="red")
lines(LDL_fle3,newdata=newdata2, var="exposure2", type="hr", col="navy")
abline(h=1, col="black", lty=2)
# RATHDL
plot(RATHDL_fle3,newdata=newdata0, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="LDL/HDL ratio", ci=TRUE, ylim=c(0,3))
lines(RATHDL_fle3,newdata=newdata0, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)
# S_ApoA
plot(S_ApoA_fle3,newdata=newdata1, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure2", main="ApoA1", ci=TRUE, ylim=c(0,3))
lines(S_ApoA_fle3,newdata=newdata1, var="exposure2", type="hr", col="navy")
abline(h=1, col="black", lty=2)
# S_ApoB
plot(S_ApoB_fle3,newdata=newdata0, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="ApoB", ci=TRUE, ylim=c(0,3))
lines(S_ApoB_fle3,newdata=newdata0, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)
# ApoKvot
plot(ApoKvot_fle3,newdata=newdata0, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="ApoB/ApoA1 ratio", ci=TRUE, ylim=c(0,3))
lines(ApoKvot_fle3,newdata=newdata0, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)
# B_Hb
plot(B_Hb_fle3,newdata=newdata3, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="Hemoglobin", ci=TRUE, ylim=c(0,3))
lines(B_Hb_fle3,newdata=newdata3, var="exposure3", type="hr", ci=TRUE, col="orange")
lines(B_Hb_fle3,newdata=newdata3, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)
# S_GLOB
plot(S_GLOB_fle3,newdata=newdata1, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure2", main="Globulin", ci=TRUE, ylim=c(0,3))
lines(S_GLOB_fle3,newdata=newdata1, var="exposure2", type="hr", col="navy")
abline(h=1, col="black", lty=2)
# S_Hapt
plot(S_Hapt_fle3,newdata=newdata4, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="Haptoglobin", ci=TRUE, ylim=c(0,3))
lines(S_Hapt_fle3,newdata=newdata4, var="exposure2", type="hr", ci=TRUE, col="navy")
lines(S_Hapt_fle3,newdata=newdata4, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)
# S_CRP
plot(S_CRP_fle3,newdata=newdata0, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="CRP", ci=TRUE, ylim=c(0,3))
lines(S_CRP_fle3,newdata=newdata0, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)
# S_IgG
plot(S_IgG_fle3,newdata=newdata1, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure2", main="IgG", ci=TRUE, ylim=c(0,4))
lines(S_IgG_fle3,newdata=newdata1, var="exposure2", type="hr", col="navy")
abline(h=1, col="black", lty=2)
# B_LPKV
plot(B_LPKV_fle3,newdata=newdata5, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure3", main="Leukocyte", ci=TRUE, ylim=c(0,3))
lines(B_LPKV_fle3,newdata=newdata5, var="exposure2", type="hr", ci=TRUE, col="navy")
lines(B_LPKV_fle3,newdata=newdata5, var="exposure3", type="hr", col="orange")
abline(h=1, col="black", lty=2)
# S_Urat
plot(S_Urat_fle3,newdata=newdata3, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="Uric acid", ci=TRUE, ylim=c(0,3))
lines(S_Urat_fle3,newdata=newdata3, var="exposure3", type="hr", ci=TRUE, col="orange")
lines(S_Urat_fle3,newdata=newdata3, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)
# S_Alb
plot(S_Alb_fle3,newdata=newdata4, type="hr",
     xlab="Time Since Cancer Daignosis (years)", var="exposure1", main="Albumin", ci=TRUE, ylim=c(0,4))
lines(S_Alb_fle3,newdata=newdata4, var="exposure2", type="hr", ci=TRUE, col="navy")
lines(S_Alb_fle3,newdata=newdata4, var="exposure1", type="hr", col="red")
abline(h=1, col="black", lty=2)

# add an empty plot
plot(1, type="n", axes=FALSE, ann=FALSE)

# Add the legend
plot(1, type="n", axes=FALSE, ann=FALSE)
legend("center", 
       legend=c("Class 1 VS. Reference", "Class 2 VS. Reference", "Class 3 VS. Reference", "Class 4 VS. Reference"), 
       col=c("red", "navy", "orange", "green"), 
       lty=1, lwd=3, cex=1.5, y.intersp=0.8, bty= "n")

dev.off() #turn off develop
#-------------------------------------------------------------------------------

# extract the average hazard ratios for each group
bio1 <- c("TC","TG","LDL","RATHDL","S_ApoB","ApoKvot","B_Hb","S_Hapt","S_CRP","S_Urat","S_Alb")
bio2 <- c("fS_Gluk","S_FAMN","HDL","LDL","S_ApoA","S_GLOB","S_Hapt","S_IgG","B_LPKV","S_Alb")
bio3 <- c("B_Hb","B_LPKV","S_Urat")
bio4 <- c("LDL")
results1 <- NULL
results2 <- NULL
results3 <- NULL
results4 <- NULL
for (i in 1:3) {
  print(i)
  for (bm1 in bio1) {
    print(bm1)
    res <- NULL
    res <- data.frame(
      outcome=i,
      biomarkers=bm1,
      coef=coef(summary(get(paste0(bm1,"_fle",i))))["exposure1",1],
      se=coef(summary(get(paste0(bm1,"_fle",i))))["exposure1",2],
      class1_pvalue=coef(summary(get(paste0(bm1,"_fle",i))))["exposure1",4])
    res <- res %>%
      mutate(lowci=coef - 1.96 * se,
             uppci=coef + 1.96 * se,
             class1_hr=exp(coef),
             class1_lowci=exp(lowci),
             class1_uppci=exp(uppci),
             class1_hrci=paste0(round(class1_hr,2)," (",round(class1_lowci,2),"-",round(class1_uppci,2),")")) %>%
      dplyr::select(-c(coef, se, lowci, uppci))
    results1 <- rbind(results1,res)
  }
  for (bm2 in bio2) {
    print(bm2)
    res <- NULL
    res <- data.frame(
      outcome=i,
      biomarkers=bm2,
      coef=coef(summary(get(paste0(bm2,"_fle",i))))["exposure2",1],
      se=coef(summary(get(paste0(bm2,"_fle",i))))["exposure2",2],
      class2_pvalue=coef(summary(get(paste0(bm2,"_fle",i))))["exposure2",4])
    res <- res %>%
      mutate(lowci=coef - 1.96 * se,
             uppci=coef + 1.96 * se,
             class2_hr=exp(coef),
             class2_lowci=exp(lowci),
             class2_uppci=exp(uppci),
             class2_hrci=paste0(round(class2_hr,2)," (",round(class2_lowci,2),"-",round(class2_uppci,2),")")) %>%
      dplyr::select(-c(coef, se, lowci, uppci))
    results2 <- rbind(results2,res)
  }
  for (bm3 in bio3) {
    print(bm3)
    res <- NULL
    res <- data.frame(
      outcome=i,
      biomarkers=bm3,
      coef=coef(summary(get(paste0(bm3,"_fle",i))))["exposure3",1],
      se=coef(summary(get(paste0(bm3,"_fle",i))))["exposure3",2],
      class3_pvalue=coef(summary(get(paste0(bm3,"_fle",i))))["exposure3",4])
    res <- res %>%
      mutate(lowci=coef - 1.96 * se,
             uppci=coef + 1.96 * se,
             class3_hr=exp(coef),
             class3_lowci=exp(lowci),
             class3_uppci=exp(uppci),
             class3_hrci=paste0(round(class3_hr,2)," (",round(class3_lowci,2),"-",round(class3_uppci,2),")")) %>%
      dplyr::select(-c(coef, se, lowci, uppci))
    results3 <- rbind(results3,res)
  }
  for (bm4 in bio4) {
    print(bm4)
    res <- NULL
    res <- data.frame(
      outcome=i,
      biomarkers=bm4,
      coef=coef(summary(get(paste0(bm4,"_fle",i))))["exposure4",1],
      se=coef(summary(get(paste0(bm4,"_fle",i))))["exposure4",2],
      class4_pvalue=coef(summary(get(paste0(bm4,"_fle",i))))["exposure4",4])
    res <- res %>%
      mutate(lowci=coef - 1.96 * se,
             uppci=coef + 1.96 * se,
             class4_hr=exp(coef),
             class4_lowci=exp(lowci),
             class4_uppci=exp(uppci),
             class4_hrci=paste0(round(class4_hr,2)," (",round(class4_lowci,2),"-",round(class4_uppci,2),")")) %>%
      dplyr::select(-c(coef, se, lowci, uppci))
    results4 <- rbind(results4,res)
  }
}

results <- results1 %>%
  full_join(results2,by=c("outcome","biomarkers")) %>%
  full_join(results3,by=c("outcome","biomarkers")) %>%
  full_join(results4,by=c("outcome","biomarkers"))

results <- results %>%
  mutate(outcomes_d=ifelse(outcome==1,"infec_p",
                           ifelse(outcome==2,"infec_both","sepsis")))
write.csv(results,".../cancer_infection/export/hazard ratios_overall.csv")

# to calculate the number of events and person-years in different classes
bios <- c("fS_Gluk","S_FAMN","TC","TG","LDL","HDL","RATHDL","S_ApoA","S_ApoB","ApoKvot","B_Hb","S_GLOB","S_Hapt","S_CRP","S_IgG","B_LPKV","S_Urat","S_Alb")
bm <- NULL
results <- NULL
for (bm in bios) {
  print(bm)
  res <- NULL
  res <- master %>%
    rename(bio=!!bm) %>%
    group_by(bio) %>%
    summarise(total=n_distinct(id),
              no_infec_p=sum(infec_p==1,na.rm = T),
              ps_infec_p=sum(time_infec_p,na.rm=T),
              no_infec_both=sum(infec_both==1,na.rm = T),
              ps_infec_both=sum(time_infec_both,na.rm=T),
              no_sepsis=sum(sepsis==1,na.rm = T),
              ps_sepsis=sum(time_sepsis,na.rm=T)) %>%
    mutate(biomarkers=bm) %>%
    rename(classes=bio) %>%
    filter(!is.na(classes))
  results <- rbind(results,res)
}
write.csv(results,".../cancer_infection/export/numbers of events and person-years.csv")

# fit multivariable models 2: additionally adjusting for infection, psychiatric disorders, and diabetes before cancer diagnosis
# fS_Gluk
toanalysis <- master %>%
  filter(!is.na(fS_Gluk)) %>%
  mutate(exposure2=ifelse(fS_Gluk==1,0,1))
fS_Gluk_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
fS_Gluk_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
fS_Gluk_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)

# S_FAMN
toanalysis <- master %>%
  filter(!is.na(S_FAMN)) %>%
  mutate(exposure2=ifelse(S_FAMN==1,0,1))
S_FAMN_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
S_FAMN_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
S_FAMN_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)

# TC
toanalysis <- master %>%
  filter(!is.na(TC)) %>%
  mutate(exposure1=ifelse(TC==2,0,1))
TC_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
TC_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
TC_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)

# TG
toanalysis <- master %>%
  filter(!is.na(TG)) %>%
  mutate(exposure1=ifelse(TG==2,0,1))
TG_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
TG_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
TG_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)

# HDL
toanalysis <- master %>%
  filter(!is.na(HDL)) %>%
  mutate(exposure2=ifelse(HDL==1,0,1))
HDL_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
HDL_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
HDL_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)

# LDL
toanalysis <- master %>%
  filter(!is.na(LDL)) %>%
  mutate(exposure0=(LDL==3)+0,
         exposure1=(LDL==1)+0,
         exposure2=(LDL==2)+0,
         exposure4=(LDL==4)+0)
LDL_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure2 + exposure4 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2, exposure4 = 2), df=3)
LDL_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure2 + exposure4 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2, exposure4 = 2), df=3)
LDL_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure2 + exposure4 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2, exposure4 = 2), df=3)

# RATHDL
toanalysis <- master %>%
  filter(!is.na(RATHDL)) %>%
  mutate(exposure1=ifelse(RATHDL==2,0,1))
RATHDL_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
RATHDL_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
RATHDL_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)

# S_ApoA
toanalysis <- master %>%
  filter(!is.na(S_ApoA)) %>%
  mutate(exposure2=ifelse(S_ApoA==1,0,1))
S_ApoA_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
S_ApoA_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
S_ApoA_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)

# S_ApoB
toanalysis <- master %>%
  filter(!is.na(S_ApoB)) %>%
  mutate(exposure1=ifelse(S_ApoB==2,0,1))
S_ApoB_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
S_ApoB_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
S_ApoB_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)

# ApoKvot
toanalysis <- master %>%
  filter(!is.na(ApoKvot)) %>%
  mutate(exposure1=ifelse(ApoKvot==2,0,1))
ApoKvot_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
ApoKvot_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
ApoKvot_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)

# B_Hb
toanalysis <- master %>%
  filter(!is.na(B_Hb)) %>%
  mutate(exposure0=(B_Hb==2)+0,
         exposure1=(B_Hb==1)+0,
         exposure3=(B_Hb==3)+0)
B_Hb_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
B_Hb_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
B_Hb_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)

# S_GLOB
toanalysis <- master %>%
  filter(!is.na(S_GLOB)) %>%
  mutate(exposure2=ifelse(S_GLOB==1,0,1))
S_GLOB_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
S_GLOB_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
S_GLOB_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)

# S_Hapt
toanalysis <- master %>%
  filter(!is.na(S_Hapt)) %>%
  mutate(exposure0=(S_Hapt==3)+0,
         exposure1=(S_Hapt==1)+0,
         exposure2=(S_Hapt==2)+0)
S_Hapt_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
S_Hapt_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
S_Hapt_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)

# S_CRP
toanalysis <- master %>%
  filter(!is.na(S_CRP)) %>%
  mutate(exposure1=ifelse(S_CRP==2,0,1))
S_CRP_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
S_CRP_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
S_CRP_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)

# S_IgG
toanalysis <- master %>%
  filter(!is.na(S_IgG)) %>%
  mutate(exposure2=ifelse(S_IgG==1,0,1))
S_IgG_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
S_IgG_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)
S_IgG_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, df=3)

# B_LPKV
toanalysis <- master %>%
  filter(!is.na(B_LPKV)) %>%
  mutate(exposure0=(B_LPKV==1)+0,
         exposure2=(B_LPKV==2)+0,
         exposure3=(B_LPKV==3)+0)
B_LPKV_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure2 = 2, exposure3 = 2), df=3)
B_LPKV_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure2 = 2, exposure3 = 2), df=3)
B_LPKV_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure2 = 2, exposure3 = 2), df=3)

# S_Urat
toanalysis <- master %>%
  filter(!is.na(S_Urat)) %>%
  mutate(exposure0=(S_Urat==2)+0,
         exposure1=(S_Urat==1)+0,
         exposure3=(S_Urat==3)+0)
S_Urat_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
S_Urat_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
S_Urat_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)

# S_Alb
toanalysis <- master %>%
  filter(!is.na(S_Alb)) %>%
  mutate(exposure0=(S_Alb==3)+0,
         exposure1=(S_Alb==1)+0,
         exposure2=(S_Alb==2)+0)
S_Alb_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
S_Alb_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
S_Alb_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll + diabete_bef_cancer + psychi_bef_cancer + infec_p_bef, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)

# extract the average hazard ratios for each group
bio1 <- c("TC","TG","LDL","RATHDL","S_ApoB","ApoKvot","B_Hb","S_Hapt","S_CRP","S_Urat","S_Alb")
bio2 <- c("fS_Gluk","S_FAMN","HDL","LDL","S_ApoA","S_GLOB","S_Hapt","S_IgG","B_LPKV","S_Alb")
bio3 <- c("B_Hb","B_LPKV","S_Urat")
bio4 <- c("LDL")
results1 <- NULL
results2 <- NULL
results3 <- NULL
results4 <- NULL
for (i in 1:3) {
  print(i)
  for (bm1 in bio1) {
    print(bm1)
    res <- NULL
    res <- data.frame(
      outcome=i,
      biomarkers=bm1,
      coef=coef(summary(get(paste0(bm1,"_fle",i))))["exposure1",1],
      se=coef(summary(get(paste0(bm1,"_fle",i))))["exposure1",2],
      class1_pvalue=coef(summary(get(paste0(bm1,"_fle",i))))["exposure1",4])
    res <- res %>%
      mutate(lowci=coef - 1.96 * se,
             uppci=coef + 1.96 * se,
             class1_hr=exp(coef),
             class1_lowci=exp(lowci),
             class1_uppci=exp(uppci),
             class1_hrci=paste0(round(class1_hr,2)," (",round(class1_lowci,2),"-",round(class1_uppci,2),")")) %>%
      dplyr::select(-c(coef, se, lowci, uppci))
    results1 <- rbind(results1,res)
  }
  for (bm2 in bio2) {
    print(bm2)
    res <- NULL
    res <- data.frame(
      outcome=i,
      biomarkers=bm2,
      coef=coef(summary(get(paste0(bm2,"_fle",i))))["exposure2",1],
      se=coef(summary(get(paste0(bm2,"_fle",i))))["exposure2",2],
      class2_pvalue=coef(summary(get(paste0(bm2,"_fle",i))))["exposure2",4])
    res <- res %>%
      mutate(lowci=coef - 1.96 * se,
             uppci=coef + 1.96 * se,
             class2_hr=exp(coef),
             class2_lowci=exp(lowci),
             class2_uppci=exp(uppci),
             class2_hrci=paste0(round(class2_hr,2)," (",round(class2_lowci,2),"-",round(class2_uppci,2),")")) %>%
      dplyr::select(-c(coef, se, lowci, uppci))
    results2 <- rbind(results2,res)
  }
  for (bm3 in bio3) {
    print(bm3)
    res <- NULL
    res <- data.frame(
      outcome=i,
      biomarkers=bm3,
      coef=coef(summary(get(paste0(bm3,"_fle",i))))["exposure3",1],
      se=coef(summary(get(paste0(bm3,"_fle",i))))["exposure3",2],
      class3_pvalue=coef(summary(get(paste0(bm3,"_fle",i))))["exposure3",4])
    res <- res %>%
      mutate(lowci=coef - 1.96 * se,
             uppci=coef + 1.96 * se,
             class3_hr=exp(coef),
             class3_lowci=exp(lowci),
             class3_uppci=exp(uppci),
             class3_hrci=paste0(round(class3_hr,2)," (",round(class3_lowci,2),"-",round(class3_uppci,2),")")) %>%
      dplyr::select(-c(coef, se, lowci, uppci))
    results3 <- rbind(results3,res)
  }
  for (bm4 in bio4) {
    print(bm4)
    res <- NULL
    res <- data.frame(
      outcome=i,
      biomarkers=bm4,
      coef=coef(summary(get(paste0(bm4,"_fle",i))))["exposure4",1],
      se=coef(summary(get(paste0(bm4,"_fle",i))))["exposure4",2],
      class4_pvalue=coef(summary(get(paste0(bm4,"_fle",i))))["exposure4",4])
    res <- res %>%
      mutate(lowci=coef - 1.96 * se,
             uppci=coef + 1.96 * se,
             class4_hr=exp(coef),
             class4_lowci=exp(lowci),
             class4_uppci=exp(uppci),
             class4_hrci=paste0(round(class4_hr,2)," (",round(class4_lowci,2),"-",round(class4_uppci,2),")")) %>%
      dplyr::select(-c(coef, se, lowci, uppci))
    results4 <- rbind(results4,res)
  }
}

results <- results1 %>%
  full_join(results2,by=c("outcome","biomarkers")) %>%
  full_join(results3,by=c("outcome","biomarkers")) %>%
  full_join(results4,by=c("outcome","biomarkers"))

results <- results %>%
  mutate(outcomes_d=ifelse(outcome==1,"infec_p",
                           ifelse(outcome==2,"infec_both","sepsis")))
write.csv(results,".../cancer_infection/export/hazard ratios_overall_additional adjustment.csv")
#-------------------------------------------------------------------------------

# Step 5. stratified analysis by sex, and history of infection before cancer diagnosis
rm(list=ls())
master <- readRDS(".../cancer_bio_trajectories.RDS")

master <- master %>%
  mutate(can_yearc=ifelse(can_year>1984 & can_year<1996, 0, 1),
         infec_both_bef=ifelse(is.na(infec_both_bef),0,infec_both_bef),
         infec_p_bef=ifelse(is.na(infec_p_bef),0,infec_p_bef),
         sepsis_bef=ifelse(is.na(sepsis_bef),0,sepsis_bef),
         lastDate=ifelse(lastDate<can_date, as.Date("2020-12-31"),lastDate),
         lastDate=as.Date(lastDate),
         end_infec_both=as.Date(pmin(lastDate,date_hos_both,na.rm = TRUE)),
         end_infec_p=as.Date(pmin(lastDate,date_hos_pri,na.rm = TRUE)),
         end_sepsis=as.Date(pmin(lastDate,date_sepsis,na.rm = TRUE)),
         age_end_infec_both=as.integer(end_infec_both - birthdate +1)/365.25,
         age_end_infec_p=as.integer(end_infec_p - birthdate+1)/365.25,
         age_end_sepsis=as.integer(end_sepsis - birthdate+1)/365.25,
         age_canc=ifelse(age_can>0 & age_can<36,0,
                         ifelse(age_can>35 & age_can<61, 1, 2)),
         time_infec_both=age_end_infec_both-age_can,
         time_infec_p=age_end_infec_p-age_can,
         time_sepsis=age_end_sepsis-age_can) %>%
  mutate(education_atenroll=ifelse(education_atenroll=="99",NA,education_atenroll),
         income_atenroll=ifelse(income_atenroll=="99",NA,income_atenroll),
         employ_atenroll=ifelse(employ_atenroll=="99",NA,employ_atenroll))

master$sex <- factor(master$sex)
master$age_canc <- factor(master$age_canc)
master$can_yearc <- factor(master$can_yearc)
master$country <- factor(master$country)  
master$education_atenroll <- factor(master$education_atenroll)
master$employ_atenroll <- factor(master$employ_atenroll) 
master$income_atenroll <- factor(master$income_atenroll)
master$diabete_bef_cancer <- factor(master$diabete_bef_cancer)
master$cvd_bef_cancer <- factor(master$cvd_bef_cancer)
master$psychi_bef_cancer <- factor(master$psychi_bef_cancer)
master$infec_p_bef <- factor(master$infec_p_bef)

bio1 <- c("TC","TG","LDL","RATHDL","S_ApoB","ApoKvot","B_Hb","S_Hapt","S_CRP","S_Urat","S_Alb")
bio2 <- c("fS_Gluk","S_FAMN","HDL","LDL","S_ApoA","S_GLOB","S_Hapt","S_IgG","B_LPKV","S_Alb")
bio3 <- c("B_Hb","B_LPKV","S_Urat")
bio4 <- c("LDL")

results <- NULL

# history of infection before cancer diagnosis
for (n in 0:1) {
  print(n)
  torun <- NULL
  torun <- master %>%
    filter(infec_p_bef==n)
  # fS_Gluk
  toanalysis <- torun %>%
    filter(!is.na(fS_Gluk)) %>%
    mutate(exposure2=ifelse(fS_Gluk==1,0,1))
  fS_Gluk_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  fS_Gluk_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  fS_Gluk_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # S_FAMN
  toanalysis <- torun %>%
    filter(!is.na(S_FAMN)) %>%
    mutate(exposure2=ifelse(S_FAMN==1,0,1))
  S_FAMN_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_FAMN_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_FAMN_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # TC
  toanalysis <- torun %>%
    filter(!is.na(TC)) %>%
    mutate(exposure1=ifelse(TC==2,0,1))
  TC_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  TC_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  TC_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # TG
  toanalysis <- torun %>%
    filter(!is.na(TG)) %>%
    mutate(exposure1=ifelse(TG==2,0,1))
  TG_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  TG_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  TG_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # HDL
  toanalysis <- torun %>%
    filter(!is.na(HDL)) %>%
    mutate(exposure2=ifelse(HDL==1,0,1))
  HDL_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  HDL_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  HDL_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # LDL
  toanalysis <- torun %>%
    filter(!is.na(LDL)) %>%
    mutate(exposure0=(LDL==3)+0,
           exposure1=(LDL==1)+0,
           exposure2=(LDL==2)+0,
           exposure4=(LDL==4)+0)
  LDL_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure2 + exposure4 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2, exposure4 = 2), df=3)
  LDL_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure2 + exposure4 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2, exposure4 = 2), df=3)
  LDL_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure2 + exposure4 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2, exposure4 = 2), df=3)
  # RATHDL
  toanalysis <- torun %>%
    filter(!is.na(RATHDL)) %>%
    mutate(exposure1=ifelse(RATHDL==2,0,1))
  RATHDL_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  RATHDL_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  RATHDL_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # S_ApoA
  toanalysis <- torun %>%
    filter(!is.na(S_ApoA)) %>%
    mutate(exposure2=ifelse(S_ApoA==1,0,1))
  S_ApoA_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_ApoA_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_ApoA_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # S_ApoB
  toanalysis <- torun %>%
    filter(!is.na(S_ApoB)) %>%
    mutate(exposure1=ifelse(S_ApoB==2,0,1))
  S_ApoB_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_ApoB_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_ApoB_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # ApoKvot
  toanalysis <- torun %>%
    filter(!is.na(ApoKvot)) %>%
    mutate(exposure1=ifelse(ApoKvot==2,0,1))
  ApoKvot_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  ApoKvot_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  ApoKvot_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # B_Hb
  toanalysis <- torun %>%
    filter(!is.na(B_Hb)) %>%
    mutate(exposure0=(B_Hb==2)+0,
           exposure1=(B_Hb==1)+0,
           exposure3=(B_Hb==3)+0)
  B_Hb_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
  B_Hb_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
  B_Hb_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
  # S_GLOB
  toanalysis <- torun %>%
    filter(!is.na(S_GLOB)) %>%
    mutate(exposure2=ifelse(S_GLOB==1,0,1))
  S_GLOB_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_GLOB_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_GLOB_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # S_Hapt
  toanalysis <- torun %>%
    filter(!is.na(S_Hapt)) %>%
    mutate(exposure0=(S_Hapt==3)+0,
           exposure1=(S_Hapt==1)+0,
           exposure2=(S_Hapt==2)+0)
  S_Hapt_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
  S_Hapt_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
  S_Hapt_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
  # S_CRP
  toanalysis <- torun %>%
    filter(!is.na(S_CRP)) %>%
    mutate(exposure1=ifelse(S_CRP==2,0,1))
  S_CRP_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_CRP_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_CRP_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # S_IgG
  toanalysis <- torun %>%
    filter(!is.na(S_IgG)) %>%
    mutate(exposure2=ifelse(S_IgG==1,0,1))
  S_IgG_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_IgG_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_IgG_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # B_LPKV
  toanalysis <- torun %>%
    filter(!is.na(B_LPKV)) %>%
    mutate(exposure0=(B_LPKV==1)+0,
           exposure2=(B_LPKV==2)+0,
           exposure3=(B_LPKV==3)+0)
  B_LPKV_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure2 = 2, exposure3 = 2), df=3)
  B_LPKV_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure2 = 2, exposure3 = 2), df=3)
  B_LPKV_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure2 = 2, exposure3 = 2), df=3)
  # S_Urat
  toanalysis <- torun %>%
    filter(!is.na(S_Urat)) %>%
    mutate(exposure0=(S_Urat==2)+0,
           exposure1=(S_Urat==1)+0,
           exposure3=(S_Urat==3)+0)
  S_Urat_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
  S_Urat_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
  S_Urat_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure3 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
  # S_Alb
  toanalysis <- torun %>%
    filter(!is.na(S_Alb)) %>%
    mutate(exposure0=(S_Alb==3)+0,
           exposure1=(S_Alb==1)+0,
           exposure2=(S_Alb==2)+0)
  S_Alb_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
  S_Alb_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
  S_Alb_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure2 + sex + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
  
  # to exract the results
  ress1 <- NULL
  for (i in 1:3) {
    print(i)
    results1 <- NULL
    results2 <- NULL
    results3 <- NULL
    results4 <- NULL
    for (bm1 in bio1) {
      print(bm1)
      res <- NULL
      res <- data.frame(
        outcome=i,
        biomarkers=bm1,
        coef=coef(summary(get(paste0(bm1,"_fle",i))))["exposure1",1],
        se=coef(summary(get(paste0(bm1,"_fle",i))))["exposure1",2],
        class1_pvalue=coef(summary(get(paste0(bm1,"_fle",i))))["exposure1",4])
      res <- res %>%
        mutate(lowci=coef - 1.96 * se,
               uppci=coef + 1.96 * se,
               class1_hr=exp(coef),
               class1_lowci=exp(lowci),
               class1_uppci=exp(uppci),
               class1_hrci=paste0(round(class1_hr,2)," (",round(class1_lowci,2),"-",round(class1_uppci,2),")")) %>%
        dplyr::select(-c(coef, se, lowci, uppci))
      results1 <- rbind(results1,res)
    }
    for (bm2 in bio2) {
      print(bm2)
      res <- NULL
      res <- data.frame(
        outcome=i,
        biomarkers=bm2,
        coef=coef(summary(get(paste0(bm2,"_fle",i))))["exposure2",1],
        se=coef(summary(get(paste0(bm2,"_fle",i))))["exposure2",2],
        class2_pvalue=coef(summary(get(paste0(bm2,"_fle",i))))["exposure2",4])
      res <- res %>%
        mutate(lowci=coef - 1.96 * se,
               uppci=coef + 1.96 * se,
               class2_hr=exp(coef),
               class2_lowci=exp(lowci),
               class2_uppci=exp(uppci),
               class2_hrci=paste0(round(class2_hr,2)," (",round(class2_lowci,2),"-",round(class2_uppci,2),")")) %>%
        dplyr::select(-c(coef, se, lowci, uppci))
      results2 <- rbind(results2,res)
    }
    for (bm3 in bio3) {
      print(bm3)
      res <- NULL
      res <- data.frame(
        outcome=i,
        biomarkers=bm3,
        coef=coef(summary(get(paste0(bm3,"_fle",i))))["exposure3",1],
        se=coef(summary(get(paste0(bm3,"_fle",i))))["exposure3",2],
        class3_pvalue=coef(summary(get(paste0(bm3,"_fle",i))))["exposure3",4])
      res <- res %>%
        mutate(lowci=coef - 1.96 * se,
               uppci=coef + 1.96 * se,
               class3_hr=exp(coef),
               class3_lowci=exp(lowci),
               class3_uppci=exp(uppci),
               class3_hrci=paste0(round(class3_hr,2)," (",round(class3_lowci,2),"-",round(class3_uppci,2),")")) %>%
        dplyr::select(-c(coef, se, lowci, uppci))
      results3 <- rbind(results3,res)
    }
    for (bm4 in bio4) {
      print(bm4)
      res <- NULL
      res <- data.frame(
        outcome=i,
        biomarkers=bm4,
        coef=coef(summary(get(paste0(bm4,"_fle",i))))["exposure4",1],
        se=coef(summary(get(paste0(bm4,"_fle",i))))["exposure4",2],
        class4_pvalue=coef(summary(get(paste0(bm4,"_fle",i))))["exposure4",4])
      res <- res %>%
        mutate(lowci=coef - 1.96 * se,
               uppci=coef + 1.96 * se,
               class4_hr=exp(coef),
               class4_lowci=exp(lowci),
               class4_uppci=exp(uppci),
               class4_hrci=paste0(round(class4_hr,2)," (",round(class4_lowci,2),"-",round(class4_uppci,2),")")) %>%
        dplyr::select(-c(coef, se, lowci, uppci))
      results4 <- rbind(results4,res)
    }
    ress <- results1 %>%
      full_join(results2,by=c("outcome","biomarkers")) %>%
      full_join(results3,by=c("outcome","biomarkers")) %>%
      full_join(results4,by=c("outcome","biomarkers")) 
    ress1 <- rbind(ress1, ress)
  }
  ress1 <- ress1 %>%
    mutate(covar="infec_p_bef",
           covar_c=n)
  results <- rbind(results,ress1)
}

# sex
ress1 <- NULL
for (n in 1:2) {
  print(n)
  torun <- NULL
  torun <- master %>%
    filter(sex==n)
  # fS_Gluk
  toanalysis <- torun %>%
    filter(!is.na(fS_Gluk)) %>%
    mutate(exposure2=ifelse(fS_Gluk==1,0,1))
  fS_Gluk_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  fS_Gluk_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  fS_Gluk_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # S_FAMN
  toanalysis <- torun %>%
    filter(!is.na(S_FAMN)) %>%
    mutate(exposure2=ifelse(S_FAMN==1,0,1))
  S_FAMN_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_FAMN_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_FAMN_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # TC
  toanalysis <- torun %>%
    filter(!is.na(TC)) %>%
    mutate(exposure1=ifelse(TC==2,0,1))
  TC_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  TC_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  TC_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # TG
  toanalysis <- torun %>%
    filter(!is.na(TG)) %>%
    mutate(exposure1=ifelse(TG==2,0,1))
  TG_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  TG_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  TG_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # HDL
  toanalysis <- torun %>%
    filter(!is.na(HDL)) %>%
    mutate(exposure2=ifelse(HDL==1,0,1))
  HDL_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  HDL_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  HDL_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # LDL
  toanalysis <- torun %>%
    filter(!is.na(LDL)) %>%
    mutate(exposure0=(LDL==3)+0,
           exposure1=(LDL==1)+0,
           exposure2=(LDL==2)+0,
           exposure4=(LDL==4)+0)
  LDL_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure2 + exposure4 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2, exposure4 = 2), df=3)
  LDL_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure2 + exposure4 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2, exposure4 = 2), df=3)
  LDL_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure2 + exposure4 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2, exposure4 = 2), df=3)
  # RATHDL
  toanalysis <- torun %>%
    filter(!is.na(RATHDL)) %>%
    mutate(exposure1=ifelse(RATHDL==2,0,1))
  RATHDL_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  RATHDL_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  RATHDL_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # S_ApoA
  toanalysis <- torun %>%
    filter(!is.na(S_ApoA)) %>%
    mutate(exposure2=ifelse(S_ApoA==1,0,1))
  S_ApoA_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_ApoA_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_ApoA_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # S_ApoB
  toanalysis <- torun %>%
    filter(!is.na(S_ApoB)) %>%
    mutate(exposure1=ifelse(S_ApoB==2,0,1))
  S_ApoB_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_ApoB_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_ApoB_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # ApoKvot
  toanalysis <- torun %>%
    filter(!is.na(ApoKvot)) %>%
    mutate(exposure1=ifelse(ApoKvot==2,0,1))
  ApoKvot_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  ApoKvot_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  ApoKvot_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # B_Hb
  toanalysis <- torun %>%
    filter(!is.na(B_Hb)) %>%
    mutate(exposure0=(B_Hb==2)+0,
           exposure1=(B_Hb==1)+0,
           exposure3=(B_Hb==3)+0)
  B_Hb_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure3 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
  B_Hb_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure3 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
  B_Hb_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure3 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
  # S_GLOB
  toanalysis <- torun %>%
    filter(!is.na(S_GLOB)) %>%
    mutate(exposure2=ifelse(S_GLOB==1,0,1))
  S_GLOB_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_GLOB_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_GLOB_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # S_Hapt
  toanalysis <- torun %>%
    filter(!is.na(S_Hapt)) %>%
    mutate(exposure0=(S_Hapt==3)+0,
           exposure1=(S_Hapt==1)+0,
           exposure2=(S_Hapt==2)+0)
  S_Hapt_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
  S_Hapt_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
  S_Hapt_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
  # S_CRP
  toanalysis <- torun %>%
    filter(!is.na(S_CRP)) %>%
    mutate(exposure1=ifelse(S_CRP==2,0,1))
  S_CRP_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_CRP_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_CRP_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # S_IgG
  toanalysis <- torun %>%
    filter(!is.na(S_IgG)) %>%
    mutate(exposure2=ifelse(S_IgG==1,0,1))
  S_IgG_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_IgG_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  S_IgG_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, df=3)
  # B_LPKV
  toanalysis <- torun %>%
    filter(!is.na(B_LPKV)) %>%
    mutate(exposure0=(B_LPKV==1)+0,
           exposure2=(B_LPKV==2)+0,
           exposure3=(B_LPKV==3)+0)
  B_LPKV_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure2 + exposure3 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure2 = 2, exposure3 = 2), df=3)
  B_LPKV_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure2 + exposure3 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure2 = 2, exposure3 = 2), df=3)
  B_LPKV_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure2 + exposure3 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure2 = 2, exposure3 = 2), df=3)
  # S_Urat
  toanalysis <- torun %>%
    filter(!is.na(S_Urat)) %>%
    mutate(exposure0=(S_Urat==2)+0,
           exposure1=(S_Urat==1)+0,
           exposure3=(S_Urat==3)+0)
  S_Urat_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure3 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
  S_Urat_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure3 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
  S_Urat_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure3 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure3 = 2), df=3)
  # S_Alb
  toanalysis <- torun %>%
    filter(!is.na(S_Alb)) %>%
    mutate(exposure0=(S_Alb==3)+0,
           exposure1=(S_Alb==1)+0,
           exposure2=(S_Alb==2)+0)
  S_Alb_fle1 <- stpm2(Surv(time_infec_p, infec_p) ~ exposure1 + exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
  S_Alb_fle2 <- stpm2(Surv(time_infec_both, infec_both) ~ exposure1 + exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
  S_Alb_fle3 <- stpm2(Surv(time_sepsis, sepsis) ~ exposure1 + exposure2 + age_can + country + can_yearc + education_atenroll + employ_atenroll + income_atenroll, data=toanalysis, tvc=list(exposure1 = 2, exposure2 = 2), df=3)
  
  # to exract the results
  ress1 <- NULL
  for (i in 1:3) {
    print(i)
    results1 <- NULL
    results2 <- NULL
    results3 <- NULL
    results4 <- NULL
    for (bm1 in bio1) {
      print(bm1)
      res <- NULL
      res <- data.frame(
        outcome=i,
        biomarkers=bm1,
        coef=coef(summary(get(paste0(bm1,"_fle",i))))["exposure1",1],
        se=coef(summary(get(paste0(bm1,"_fle",i))))["exposure1",2],
        class1_pvalue=coef(summary(get(paste0(bm1,"_fle",i))))["exposure1",4])
      res <- res %>%
        mutate(lowci=coef - 1.96 * se,
               uppci=coef + 1.96 * se,
               class1_hr=exp(coef),
               class1_lowci=exp(lowci),
               class1_uppci=exp(uppci),
               class1_hrci=paste0(round(class1_hr,2)," (",round(class1_lowci,2),"-",round(class1_uppci,2),")")) %>%
        dplyr::select(-c(coef, se, lowci, uppci))
      results1 <- rbind(results1,res)
    }
    for (bm2 in bio2) {
      print(bm2)
      res <- NULL
      res <- data.frame(
        outcome=i,
        biomarkers=bm2,
        coef=coef(summary(get(paste0(bm2,"_fle",i))))["exposure2",1],
        se=coef(summary(get(paste0(bm2,"_fle",i))))["exposure2",2],
        class2_pvalue=coef(summary(get(paste0(bm2,"_fle",i))))["exposure2",4])
      res <- res %>%
        mutate(lowci=coef - 1.96 * se,
               uppci=coef + 1.96 * se,
               class2_hr=exp(coef),
               class2_lowci=exp(lowci),
               class2_uppci=exp(uppci),
               class2_hrci=paste0(round(class2_hr,2)," (",round(class2_lowci,2),"-",round(class2_uppci,2),")")) %>%
        dplyr::select(-c(coef, se, lowci, uppci))
      results2 <- rbind(results2,res)
    }
    for (bm3 in bio3) {
      print(bm3)
      res <- NULL
      res <- data.frame(
        outcome=i,
        biomarkers=bm3,
        coef=coef(summary(get(paste0(bm3,"_fle",i))))["exposure3",1],
        se=coef(summary(get(paste0(bm3,"_fle",i))))["exposure3",2],
        class3_pvalue=coef(summary(get(paste0(bm3,"_fle",i))))["exposure3",4])
      res <- res %>%
        mutate(lowci=coef - 1.96 * se,
               uppci=coef + 1.96 * se,
               class3_hr=exp(coef),
               class3_lowci=exp(lowci),
               class3_uppci=exp(uppci),
               class3_hrci=paste0(round(class3_hr,2)," (",round(class3_lowci,2),"-",round(class3_uppci,2),")")) %>%
        dplyr::select(-c(coef, se, lowci, uppci))
      results3 <- rbind(results3,res)
    }
    for (bm4 in bio4) {
      print(bm4)
      res <- NULL
      res <- data.frame(
        outcome=i,
        biomarkers=bm4,
        coef=coef(summary(get(paste0(bm4,"_fle",i))))["exposure4",1],
        se=coef(summary(get(paste0(bm4,"_fle",i))))["exposure4",2],
        class4_pvalue=coef(summary(get(paste0(bm4,"_fle",i))))["exposure4",4])
      res <- res %>%
        mutate(lowci=coef - 1.96 * se,
               uppci=coef + 1.96 * se,
               class4_hr=exp(coef),
               class4_lowci=exp(lowci),
               class4_uppci=exp(uppci),
               class4_hrci=paste0(round(class4_hr,2)," (",round(class4_lowci,2),"-",round(class4_uppci,2),")")) %>%
        dplyr::select(-c(coef, se, lowci, uppci))
      results4 <- rbind(results4,res)
    }
    ress <- results1 %>%
      full_join(results2,by=c("outcome","biomarkers")) %>%
      full_join(results3,by=c("outcome","biomarkers")) %>%
      full_join(results4,by=c("outcome","biomarkers")) 
    ress1 <- rbind(ress1, ress) 
  }
  ress1 <- ress1 %>%
    mutate(covar="sex",
           covar_c=n)
  results <- rbind(results,ress1) 
}

write.csv(results,".../cancer_infection/export/hazard ratios_stratified analyses.csv")
#-------------------------------------------------------------------------------

