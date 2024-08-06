## Test Code for Stiffness Model Calculations
## Compared to previous calculations from a different software

## Create fake ultrasound and BP data
test_data <- data.frame (
  SBP = seq(100, 200, by=10),
  DBP = seq(50, 100, by=5),
  diam_dia = seq(5, 6, by=0.1),
  diam_sys = seq(6, 7, by=0.1)
  )
## Code to calculate stiffness removed due to patent application from UW-Madison

## Load data and clean for missing variables

## Load necessary packages (and possible some unnecessary ones too)
library(haven)
library("survival")
library("survminer")
library("Hmisc")
library("forestplot")
library("dplyr")
library("ggplot2")
library("tibble")

## Calculate Stiffness Models
setwd("C:/Users/pewow/Documents/MESA/MESA HF and A-Fib")
Raw_stiffness_data <-  read.csv("MESAe1Disten07022014.csv")

## Create stiffness data frame
Stiffness_models <- data.frame(
  ID = Raw_stiffness_data$idno,
  diam_sys = Raw_stiffness_data$sdiamc1,
  diam_dia = with(Raw_stiffness_data,  sdiamc1 - ddiam1),
  SBP = Raw_stiffness_data$sbp1,
  DBP = with(Raw_stiffness_data,  sbp1 - dbp1)
)

## Code to calculate stiffness removed due to patent application from UW-Madison

## Read in Demographics and Co-variates

Raw_Baseline <-  read_sav("MESAE1FinalLabel20240307.sav")

## Create RAAS inhibitor variable for participants taking ARBs, ACEis or ARB/ACEi + diuretic
RAASi <- data.frame(
  ARB = as.numeric(Raw_Baseline$a2a1c),
  ARBd = as.numeric(Raw_Baseline$a2ad1c),
  ACEI = as.numeric(Raw_Baseline$ace1c),
  ACEId = as.numeric(Raw_Baseline$aced1c)
)

RAASi$RAASi_sum <- with(RAASi, ARB + ARBd + ACEI +ACEId) 
RAASi$RAASi_YN <- with(RAASi, as.numeric(RAASi_sum>0)) 

## Create demographics and co-variates data frame
Baseline_data <- data.frame(
  ID =Raw_Baseline$idno,
  age = Raw_Baseline$age1c,
  race = factor(Raw_Baseline$race1c),
  sex = factor(Raw_Baseline$gender1),
  site = factor(Raw_Baseline$site1c),
  educ = ordered(Raw_Baseline$educ1, levels=c(seq(0, 8, 1))),
  bmi = Raw_Baseline$bmi1c,
  smoke = factor(Raw_Baseline$cig1c),
  dm = factor(Raw_Baseline$dm031c),
  SBP_main = Raw_Baseline$sbp1c,
  DBP_main = Raw_Baseline$dbp1c,
  PP_main = Raw_Baseline$spp1c,
  MAP_main = Raw_Baseline$dbp1c + Raw_Baseline$spp1c/3,
  htn = factor(Raw_Baseline$htn1c),
  bp_med = factor(Raw_Baseline$htnmed1c),
  RAASi = factor(RAASi$RAASi_YN ),
  chol = Raw_Baseline$chol1,
  ldl = Raw_Baseline$ldl1,
  hdl = Raw_Baseline$hdl1,
  trig = Raw_Baseline$trig1,
  lipid_med = factor(Raw_Baseline$lipid1c)

)

## Load Events Data
Raw_CVD_Events <-  read_sav("MESAEvThru2019_20220727.sav")
RAW_EF_Events <- read_sav("MESAEvThru2019CHFEF_20220727.sav")
Raw_AF_Events <- read_sav("MESA_AFThru2018_20210614.sav")

## Label HF events as HFrEF
RAW_EF_Events$HFrEF <- as.numeric ( with(RAW_EF_Events,
                            (efclass==1 & efmeas < 50) |
                            (efclass==2 & efmeas < 50) |
                              (efclass==5)
                            ) )
## Label HF events as HFpEF
RAW_EF_Events$HFpEF <- as.numeric ( with(RAW_EF_Events,
                                         (efclass==1 & efmeas >= 50) |
                                           (efclass==3 & efmeas > 50) |
                                           (efclass==4)
                            ) )
idx_indeterminate <- which(RAW_EF_Events$HFpEF + RAW_EF_Events$HFrEF == 0)

RAW_EF_Events$HFrEF[idx_indeterminate ] <- NA
RAW_EF_Events$HFpEF[idx_indeterminate ] <- NA

## Create HF events dataframe
HF_Data <- data.frame(
  ID = Raw_CVD_Events$idno,
  FU_time_HF = Raw_CVD_Events$fuptt,
  HF = Raw_CVD_Events$chf,
  HF_time = Raw_CVD_Events$chftt
)

## Create AF diagnosis dataframe
AF_Data <- data.frame(
  ID = Raw_AF_Events$idno,
  FU_time_AF = Raw_AF_Events$fuptt,
  AF = Raw_AF_Events$af2018,
  AF0 = Raw_AF_Events$pbaf,
  AF_time = Raw_AF_Events$af2018tt
)

HF_Data <- subset(HFAF_Data, !is.na(HFAF_Data$HF))

options(warn=1)

for (i in 1:length(HF_Data$ID) )  {
  
  id = HF_Data$ID[i]
  
  if ( HF_Data$HF[i] == 0 ) { HF_Data$HFrEF[i] <- 0
  HF_Data$HFpEF[i] <- 0 } else {
    
    idx<-which(RAW_EF_Events$idno == id)
    
    if (length(idx) == 0) { 
      HF_Data$HFrEF[i] <- NA
      HF_Data$HFpEF[i] <- NA
    } else {
    
    temp <- which(RAW_EF_Events$ttchf[idx] == min(RAW_EF_Events$ttchf[idx]) )
    
    if (length(temp > 1)) {
      temp_temp <- which(!is.na(RAW_EF_Events$HFrEF [ idx[temp] ]))
      
      if (length(temp_temp) == 1) { temp <- temp[ temp_temp ]
      } else {temp <- temp[ temp_temp [1] ]}
     
    }
    
    idx_event1 <- idx[ temp ]
    HF_Data$HFrEF[i] <- RAW_EF_Events$HFrEF[ idx_event1 ]
    HF_Data$HFpEF[i] <- RAW_EF_Events$HFpEF[ idx_event1 ]
    }
  }
  
}

######## Clean and Merge Datasets
data_merge <- merge(Baseline_data, Stiffness_models)
data_merge <- merge(data_merge, HF_Data)
data_merge <- merge(data_merge, AF_Data)

## Missing demographics
MISSING_Demo <- is.na(data_merge$age) |
  is.na(data_merge$race) |
  is.na(data_merge$sex) |
  is.na(data_merge$site) |
  is.na(data_merge$educ) |
  is.na(data_merge$bmi) |
  is.na(data_merge$smoke) |
  is.na(data_merge$dm) |
  is.na(data_merge$PP_main) |
  is.na(data_merge$MAP_main) |
  is.na(data_merge$htn) |
  is.na(data_merge$bp_med) |
  is.na(data_merge$RAASi) |
  is.na(data_merge$chol) |
  is.na(data_merge$hdl) |
  is.na(data_merge$lipid_med)
sum(MISSING_Demo)

## Missing stiffness
MISSING_Stiff <- is.na(data_merge$Total_PWV)

sum(MISSING_Stiff)

MISSING = MISSING_Stiff | MISSING_Demo

sum(MISSING)

data_clean <- subset(data_merge, subset = !MISSING)
nrow(data_clean)

Q3_struct <- quantile(data_clean$Struct_PWV, probs = 0.75)
Q3_LD <- quantile(data_clean$LD_PWV, probs = 0.75)

data_clean$PWV_cat <- as.numeric(data_clean$Struct_PWV > Q3_struct) +
  2* as.numeric(data_clean$LD_PWV > Q3_LD)


## Created AF data set excluding participants who couldn't be determined at baseline if they had AF
Missing_AF <- is.na(data_clean$AF)
No_Baseline_AF <- (data_clean$AF0==1) | is.na(data_clean$AF0)
Exclude_AF <- Missing_AF | No_Baseline_AF

data_clean_AF <- subset(data_clean, subset = !Exclude_AF )
nrow(data_clean_AF)

###### Make Table 1 - Participant Characteristics
Table1 <- with(data_clean, data.frame(age, sex, race, bmi, educ, 
                                         SBP_main, DBP_main, htn, bp_med, RAASi, 
                                         dm, chol, ldl, hdl, trig, lipid_med,
                                         smoke, Total_PWV, Struct_PWV, LD_PWV
                                        ) )
describe(Table1)

## Make Figure 2 - Unadjusted KM Curves
km_fit_HF <- survfit(Surv(HF_time/365, HF) ~ PWV_cat, data=data_clean)
ggsurvplot(km_fit_HF, 
           title = c("Heart Failure Events"),
           censor=FALSE,
           conf.int=TRUE, 
           ylim = c(0.5, 1),
           legend = c(0.25, 0.25),
           palette = "bmj",
           legend.labs = c("Low Both", "High Structural", "High Load-Dep.", "High Both"),
           ylab = c("Survival"),
           xlab = c("Time (Years)")
)

km_fit_HFrEF <- survfit(Surv(HF_time/365, HFrEF) ~ PWV_cat, data=data_clean)
ggsurvplot(km_fit_HFrEF, 
           title = c("HFrEF Events"),
           censor=FALSE,
           conf.int=TRUE, 
           ylim = c(0.5, 1),
           legend = c(0.25, 0.25),
           palette = "bmj",
           legend.labs = c("Low Both", "High Structural", "High Load-Dep.", "High Both"),
           ylab = c("Survival")
           )

km_fit_HFpEF <- survfit(Surv(HF_time/365, HFpEF) ~ PWV_cat, data=data_clean)
ggsurvplot(km_fit_HFpEF, 
           title = c("HFpEF Events"),
           censor=FALSE,
           conf.int=TRUE, 
           ylim = c(0.5, 1),
           legend = c(0.25, 0.25),
           palette = "bmj",
           legend.labs = c("Low Both", "High Structural", "High Load-Dep.", "High Both"),
           ylab = c("Survival")
)

km_fit_AF <- survfit(Surv(AF_time/365, AF) ~ PWV_cat, data=data_clean_AF)
ggsurvplot(km_fit_AF, 
           title = c("Atrial Fibrillation Events"),
           censor=FALSE,
           conf.int=TRUE, 
           ylim = c(0.5, 1),
           legend = c(0.25, 0.25),
           palette = "bmj",
           legend.labs = c("Low Both", "High Structural", "High Load-Dep.", "High Both"),
           ylab = c("Survival")
)

## Fit Total Effect Cox Models 
## All HF, HFrEF, HFpEF, A-fib
HF_Total.C0 <- coxph(Surv(HF_time, HF) ~ Total_PWV + age + sex + race + educ + site +
                       bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                       lipid_med + MAP_main, data = data_clean)
HF_Struct.C0 <- coxph(Surv(HF_time, HF) ~ Struct_PWV + age + sex + race + educ + site +
                       bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                        lipid_med + MAP_main, data = data_clean)
HF_LD.C0 <- coxph(Surv(HF_time, HF) ~ LD_PWV + age + sex + race + educ + site +
                        bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                    lipid_med + MAP_main, data = data_clean)

HFrEF_Total.C0 <- coxph(Surv(HF_time, HFrEF) ~ Total_PWV + age + sex + race + educ + site +
                       bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                         lipid_med + MAP_main, data = data_clean)
HFrEF_Struct.C0 <- coxph(Surv(HF_time, HFrEF) ~ Struct_PWV + age + sex + race + educ + site +
                        bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                          lipid_med + MAP_main, data = data_clean)
HFrEF_LD.C0 <- coxph(Surv(HF_time, HFrEF) ~ LD_PWV + age + sex + race + educ + site +
                    bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                      lipid_med + MAP_main, data = data_clean)

HFpEF_Total.C0 <- coxph(Surv(HF_time, HFpEF) ~ Total_PWV + age + sex + race + educ + site +
                          bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                          lipid_med + MAP_main, data = data_clean)
HFpEF_Struct.C0 <- coxph(Surv(HF_time, HFpEF) ~ Struct_PWV + age + sex + race + educ + site +
                           bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                           lipid_med + MAP_main, data = data_clean)
HFpEF_LD.C0 <- coxph(Surv(HF_time, HFpEF) ~ LD_PWV + age + sex + race + educ + site +
                       bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                       lipid_med + MAP_main, data = data_clean)

AF_Total.C0 <- coxph(Surv(AF_time, AF) ~ Total_PWV + age + sex + race + educ + site +
                          bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                          lipid_med + MAP_main, data = data_clean_AF)
AF_Struct.C0 <- coxph(Surv(AF_time, AF) ~ Struct_PWV + age + sex + race + educ + site +
                           bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                           lipid_med + MAP_main, data = data_clean_AF)
AF_LD.C0 <- coxph(Surv(AF_time, AF) ~ LD_PWV + age + sex + race + educ + site +
                       bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                       lipid_med + MAP_main, data = data_clean_AF)

## Fit Direct Effect Cox Models (Add PP as co-variate)
## All HF, HFrEF, HFpEF, A-fib
HF_Total.C1 <- coxph(Surv(HF_time, HF) ~ Total_PWV + age + sex + race + educ + site +
                       bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                       lipid_med + MAP_main + PP_main, data = data_clean)
HF_Struct.C1 <- coxph(Surv(HF_time, HF) ~ Struct_PWV + age + sex + race + educ + site +
                        bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                        lipid_med + MAP_main + PP_main, data = data_clean)
HF_LD.C1 <- coxph(Surv(HF_time, HF) ~ LD_PWV + age + sex + race + educ + site +
                    bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                    lipid_med + MAP_main + PP_main, data = data_clean)

HFrEF_Total.C1 <- coxph(Surv(HF_time, HFrEF) ~ Total_PWV + age + sex + race + educ + site +
                          bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                          lipid_med + MAP_main + PP_main, data = data_clean)
HFrEF_Struct.C1 <- coxph(Surv(HF_time, HFrEF) ~ Struct_PWV + age + sex + race + educ + site +
                           bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                           lipid_med + MAP_main + PP_main, data = data_clean)
HFrEF_LD.C1 <- coxph(Surv(HF_time, HFrEF) ~ LD_PWV + age + sex + race + educ + site +
                       bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                       lipid_med + MAP_main + PP_main+ PP_main, data = data_clean)

HFpEF_Total.C1 <- coxph(Surv(HF_time, HFpEF) ~ Total_PWV + age + sex + race + educ + site +
                          bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                          lipid_med + MAP_main + PP_main, data = data_clean)
HFpEF_Struct.C1 <- coxph(Surv(HF_time, HFpEF) ~ Struct_PWV + age + sex + race + educ + site +
                           bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                           lipid_med + MAP_main + PP_main, data = data_clean)
HFpEF_LD.C1 <- coxph(Surv(HF_time, HFpEF) ~ LD_PWV + age + sex + race + educ + site +
                       bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                       lipid_med + MAP_main + PP_main, data = data_clean)

AF_Total.C1 <- coxph(Surv(AF_time, AF) ~ Total_PWV + age + sex + race + educ + site +
                       bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                       lipid_med + MAP_main + PP_main, data = data_clean_AF)
AF_Struct.C1 <- coxph(Surv(AF_time, AF) ~ Struct_PWV + age + sex + race + educ + site +
                        bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                        lipid_med + MAP_main + PP_main, data = data_clean_AF)
AF_LD.C1 <- coxph(Surv(AF_time, AF) ~ LD_PWV + age + sex + race + educ + site +
                    bmi + smoke + dm + htn + bp_med + RAASi + chol + hdl +
                    lipid_med + MAP_main + PP_main, data = data_clean_AF)


## Make Figure 3: Forest Plots for Most Fully Adjusted Models

base_data <- tibble::tibble(mean  = c(10, 1.091, 1.064, 1.234, 
                                      10, 1.109, 1.096, 1.124),
                            lower = c(10, 0.982, 0.938, 1.053, 
                                      10, 1.036, 1.016, 1.020),
                            upper = c(10, 1.223, 1.178, 1.444, 
                                      10, 1.204, 1.161, 1.235),
                            Metric = c("Heart Failure", "Total PWV", "Struct. PWV", "Load-Dep. PWV",
                                       "Atrial Fibrillation", "Total PWV", "Struct. PWV", "Load-Dep. PWV")
                                        )
base_data |>
  forestplot(labeltext = c(Metric),
             boxsize = 0.4,
             clip = c(0.5, 1.5),
             xlab = "HR per 1 SD",
             xticks = c(-0.105, 0, 0.095, 0.182, .262, .336, .405),
             xlog = TRUE,
             xticks.digits = 1,
             lwd.ci = 4,
             lwd.xaxis = 4,
             lwd.zero = 4) |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               zero = "darkgrey",
               txt_gp = fpTxtGp(ticks = gpar(cex = 1.2),
                                xlab = gpar(cex = 1.3),
                                title = gpar(cex = 1.5))
  ) |>
  fp_set_zebra_style("#EFEFEF")

