# AUC calculation from virtual population generated using SMAC data and the LMS method for age vs WT as well as age vs HB


rm(list = ls()) # clear memory
graphics.off()  # clear graphics
cat("\014")     # clear console

# setting the working directory
setwd("XXX")


# load R libraries
library(mStats)
library("ggrepel")                             
library(magrittr) 
library(shiny)
library(shinycssloaders)
library(lhs)
library(VGAM)
library(parallel)
library(dplyr)
library(tictoc)
library(ggplot2)
library(plotly)
library(latex2exp)
library(scales)
library(Hmisc)

# binned data (1000 samples per kg) from virtual population 
VirtPopData=read.csv("VirtPopData_from_SMAC.csv")

#######################################################################################
# Simulated total first-dose exposure levels (AUC0-12h) of DHA after:
# the standard 2.4 mg/kg dosing in children at different body weights and 
# the WHO endorsed dose (3 mg/kg <20 kg or 2.4 mg/kg for >= 20 kg)
#######################################################################################

# DHA PopPK parameters adopted from our current PK study (IV ARS) using the REACH dataset -- Bayesian PK modeling
Parms = c(
  CL = 28.95, # clearance [l/h]
  V = 18.26,  # Volume of the central compartment [l]
  Q= 1.90 ,   # inter-compartmental clearance [l/h]
  Vp = 3.97,  # Volume of the peripheral compartment [l]
  kmat = 0.082, # age-related enzyme-maturation effect 
  betaCL.HB = -0.0356,# fractional change in population mean CL for a unit (g/dl) increase in hemoglobin from the mean 
  sigma2CL = 0.1853, #  variance  
  sigma2V = 0.0656,  # variance 
  sigma2Q2 = 0.8867, # variance
  sigma2Vp =0.4704,  # variance
  covCLV = 0.06205,  # covariance  
  covCLQ2= 0.0753,   # covariance
  covCLVp = -0.0244, # covariance
  covVQ2 =-0.02384,  # covariance
  covVVp = 0.01048,  # covariance
  covQ2Vp = 0.06214, # covariance
  WGTpop = 17.93,    # The reference body weight in the population PK modeling 
  HBpop = 7.91)      # The reference hemoglobin value in the population PK modeling 

# build the covariance matrix 
CovMat=matrix(seq(1,16),nrow=4, byrow = TRUE)

CovMat[1,1]=Parms["sigma2CL"]
CovMat[2,2]=Parms["sigma2V"]
CovMat[3,3]=Parms["sigma2Q2"]
CovMat[4,4]=Parms["sigma2Vp"]

CovMat[2,1]=Parms["covCLV"]
CovMat[3,1]=Parms["covCLQ2"]
CovMat[4,1]=Parms["covCLVp"]

CovMat[3,2]=Parms["covVQ2"]
CovMat[4,2]=Parms["covVVp"]
CovMat[4,3]=Parms["covQ2Vp"]

CovMat_sym <- CovMat                                        # Duplicate matrix
CovMat_sym[upper.tri(CovMat_sym)] <- t(CovMat_sym)[upper.tri(CovMat_sym)] # Insert lower to upper matrix

eta = MASS::mvrnorm(n=nrow(VirtPopData), mu = c(0, 0,0,0), 
                    Sigma = CovMat_sym)

# calculate rounded weight for dose calculation and then the 2 dosage types:
# dose1 is the weight independent dosing (FDA endorsed dosing)
# dose2 is the weight based dosing (WHO endorsed dosing)

VirtPopData_AUC1<- VirtPopData %>% mutate(dose = 2.4*rounded_weight,
                                                  dose_type = "FDA endorsed dosing")
VirtPopData_AUC2<- VirtPopData %>% mutate(dose = case_when(rounded_weight<20 ~ 3*rounded_weight,
                                                                   TRUE ~ 2.4*rounded_weight),
                                                  dose_type = "WHO endorsed dosing")

VirtPopData_AUC <- rbind(VirtPopData_AUC1, VirtPopData_AUC2) %>% 
  dplyr::select(Age, Sex, WT, HB, TEMP, rounded_weight, dose, dose_type)


# calculate DHA CL_pop:
VirtPopData_AUC<- VirtPopData_AUC %>% mutate(CL_Fpop = Parms["CL"] * 
                                                                   exp(Parms["betaCL.HB"]*(HB - Parms["HBpop"]))) #hb covariate
# calculate individual clearance (CL_i)
VirtPopData_AUC<- VirtPopData_AUC %>%
  mutate(CL_Fi = ((CL_Fpop/3.1)*((WT/Parms["WGTpop"])^0.75)*(exp(-Parms["kmat"]*Age)) + (CL_Fpop)*((WT/Parms["WGTpop"])^0.75)*(1 - exp(-Parms["kmat"]*Age)))*exp(eta[,1]))

# calculate individual volume (V_i)
VirtPopData_AUC<- VirtPopData_AUC %>% mutate(V_fi = Parms["V"] * (WT/Parms["WGTpop"])*exp(eta[,2]),
                                                                 Q_i= Parms["Q"]*((WT/Parms["WGTpop"])^0.75)*exp(eta[,3]),
                                                                 Vp_i = Parms["Vp"]*((WT/Parms["WGTpop"])^1.0)*exp(eta[,4])
)

VirtPopData_AUC<- VirtPopData_AUC %>% mutate(beta = 0.5*(Q_i/V_fi+Q_i/Vp_i+CL_Fi/V_fi-sqrt((Q_i/V_fi+Q_i/Vp_i+CL_Fi/V_fi)^2-4*(Q_i/Vp_i)*(CL_Fi/V_fi))),
                                                                 alpha = (Q_i/Vp_i)*(CL_Fi/V_fi)/beta,
                                                                 AA = (alpha-(Q_i/Vp_i))/(V_fi*(alpha-beta)),
                                                                 BB = (beta-(Q_i/Vp_i))/(V_fi*(beta-alpha)))
#calculate drug exposure DHA: 
VirtPopData_AUC<- VirtPopData_AUC %>% 
  mutate(AUC_dose = (dose*(284/384)*1000)*( AA/alpha *(1-exp(-12*alpha)) +  BB/beta *(1-exp(-12*beta)) )) 

# NOTES:
# dose times by (284/384) to account for 
# change in molecular weight between ARS and DHA, assuming 100% conversion; 
# times by 1000 is to convert from hr*mg/L to hr*ng/mL

#######################################################################################
# calculation of Median DHA exposure (h x ng/ml), with 25th and 75th percentiles, 
# after different dose regimens with weights binned into groups for data 
#######################################################################################
VirtPopData_AUC$dose_type = as.factor(VirtPopData_AUC$dose_type)

VirtPopData_AUC<-VirtPopData_AUC %>% group_by(rounded_weight, dose_type) %>% 
  mutate(AUC_med_dose = median(AUC_dose),
         AUC_25_dose = quantile(AUC_dose, probs=0.25, na.rm=T),
         AUC_75_dose = quantile(AUC_dose, probs=0.75,na.rm=T))

pd <- position_dodge(width = 0.8)

result_AUC<-ggplot()+
  geom_point(data=VirtPopData_AUC, aes(x=rounded_weight, y=AUC_med_dose, colour=dose_type), position=pd)+
  geom_errorbar(data=VirtPopData_AUC, aes(x=rounded_weight, y=AUC_med_dose, ymin=AUC_25_dose, ymax=AUC_75_dose, colour=dose_type),
                width=0.4, position=pd)+
  geom_vline(xintercept = 20, linetype='dotted', colour="black") +
  xlab("Bodyweight (kg)")+ylab(expression(DHA~AUC[0-12~h]~(h~x~ng/ml)))+
  theme_bw()+
  scale_x_continuous(breaks=seq(2,40,2))+
  scale_y_continuous(breaks=seq(0,4000,1000), limits = c(0,4000))+
  theme(text = element_text(size = 12),
        legend.justification = c(1,1), legend.position = c(0.99,0.99),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_colour_manual(name = NULL, values = c("red", "black"))

ggsave(plot = result_AUC, filename = paste0("result_AUC_","_PK_from_REACH.jpeg"),
       device = 'jpeg',width = 6, height = 3, dpi = 300)
#######################################################################################

