
rm(list = ls()) # clear memory
graphics.off()  # clear graphics
cat("\014")     # clear console

# setting the working directory
setwd("XXX")


# load R libraries
pacman::p_load(tidyverse, ggplot2,Hmisc,haven, dplyr, psych, gamlss,  car, ggpubr, 
               MASS, stringr)

# source the LMS function  
source("lms_w_GAIC.R")

# set seed
set.seed(31424534)

# load datasets: SMAC NETWORK DATASET

full<-rbind(gambia,ghana,kenya,gabonlam,gabonlib,malawi)

# create id for full dataset for merging purposes
full<-full %>% mutate(id=1:nrow(full))

# change sex and country to factor variables 
full$SEX<-factor(full$SEX, level=c(1,2), labels=c("M", "F"))
full$COUNTRY<-factor(full$COUNTRY, level=c(1,2,3,4,5,6), 
                     labels=c("Gambia", "Ghana", "Kenya","Gabon (Lamberene)",
                              "Gabon (Libreville)", "Malawi"))

########################################################################################
# creating datasets for each type of relationship being looked at, 

full_weight_notclean <- full  %>% dplyr::select(id, AGE, SEX, WEIGHT) 
full_HB_notclean <- full %>% dplyr::select(id, AGE, SEX, HB)
full_TEMP_notclean <- full %>% dplyr::select(id, AGE, SEX, TEMP)

# restricting to 150 months 
full_weight_notclean_below150 <- full_weight_notclean %>% filter(AGE<=150) 
full_HB_notclean_below150  <- full_HB_notclean %>% filter(AGE<=150)
full_TEMP_notclean_below150  <- full_TEMP_notclean %>% filter(AGE<=150) 

# create sex specific complete datasets for weight
complete_weight_boy<-full_weight_notclean_below150 %>% filter(SEX=="M")
complete_weight_girl<-full_weight_notclean_below150 %>% filter(SEX=="F")

#######################################################################################

# LMS calculation WEIGHT vs AGE
m_boy_raw<-lms_w_GAIC(WEIGHT,AGE,families=c("BCCG", "NO"),data=complete_weight_boy,
                      k=4,calibration=F, trans.x=T)

m_girl_raw<-lms_w_GAIC(WEIGHT,AGE,families=c("BCCG", "NO"),data=complete_weight_girl,
                       k=4,calibration=F, trans.x=T,
                       sigma.df = 5, mu.df=10)

# put LMS parameters into dataframe 
complete_weight_boy$M_lms<-m_boy_raw$mu.fv #M
complete_weight_boy$S_lms<-m_boy_raw$sigma.fv #S
complete_weight_boy$L_lms<-m_boy_raw$nu.fv #Lambda

complete_weight_girl$M_lms<-m_girl_raw$mu.fv #M
complete_weight_girl$S_lms<-m_girl_raw$sigma.fv #S
complete_weight_girl$L_lms<-m_girl_raw$nu.fv #Lambda

#######################################################################################

# LMS calculation HB vs AGE
m_HB_raw<-lms_w_GAIC(HB,AGE, data=full_HB,families=c("NO", "BCCG"),
                     k=4,calibration=F, trans.x=T)

#put LMS parameters into dataframe 
full_HB$M_lms<-m_HB_raw$mu.fv #M
full_HB$S_lms<-m_HB_raw$sigma.fv #S
full_HB$L_lms<-m_HB_raw$nu.fv #L

#######################################################################################
# create complete LMS with age dataframe with P1 and P99 as the outer limits: for both weight and HB

weight_boy_raw_LMS_pre<-complete_weight_boy %>% mutate(Sex=1) %>%
  dplyr::select(AGE, WEIGHT, Sex,SEX, M=M_lms, S=S_lms, L=L_lms) %>% group_by(AGE) %>% 
  mutate(Q1 = quantile(WEIGHT, probs=0.01, na.rm=T),
         Q99 = quantile(WEIGHT, probs=0.99, na.rm=T))

weight_boy_raw_LMS<-weight_boy_raw_LMS_pre %>% distinct(AGE, .keep_all=T)

# smooth percentiles to create upper and lower limits of sampling

boyweightsmoothQ1_1<- predict(smooth.spline(x=weight_boy_raw_LMS$AGE,y=weight_boy_raw_LMS$Q1, df = 2))
boyweightsmoothQ99_1<- predict(smooth.spline(x=weight_boy_raw_LMS$AGE,y=weight_boy_raw_LMS$Q99, df = 2))
boyweightsmoothdf<-data.frame(AGE=boyweightsmoothQ1_1$x, P1 = boyweightsmoothQ1_1$y,P99 = boyweightsmoothQ99_1$y)

weight_girl_raw_LMS_pre<-complete_weight_girl %>% mutate(Sex=0) %>%
  dplyr::select(AGE, WEIGHT, Sex,SEX,M=M_lms, S=S_lms, L=L_lms) %>% group_by(AGE) %>% 
  mutate(Q1 = quantile(WEIGHT, probs=0.01, na.rm=T),
         Q99 = quantile(WEIGHT, probs=0.99, na.rm=T)) 

weight_girl_raw_LMS<-weight_girl_raw_LMS_pre %>% distinct(AGE, .keep_all=T)

# smooth percentiles to create upper and lower limits of sampling
girlweightsmoothQ1_1<- predict(smooth.spline(x=weight_girl_raw_LMS$AGE,y=weight_girl_raw_LMS$Q1, df = 5))
girlweightsmoothQ99_1<- predict(smooth.spline(x=weight_girl_raw_LMS$AGE,y=weight_girl_raw_LMS$Q99, df = 5))
girlweightsmoothdf<-data.frame(AGE=girlweightsmoothQ1_1$x, P1 = girlweightsmoothQ1_1$y,P99 = girlweightsmoothQ99_1$y)

# merge back to other dataframe
weight_girl_raw_LMS<-weight_girl_raw_LMS %>% left_join(girlweightsmoothdf) %>% 
  dplyr::select(AGE, Sex,SEX, M, S, L, P1, P99)
weight_boy_raw_LMS<-weight_boy_raw_LMS %>% left_join(boyweightsmoothdf) %>% 
  dplyr::select(AGE, Sex,SEX, M, S, L, P1, P99)

# calculate 1st and 99th percentiles by age for HB
full_HB_lms_pre<-full_HB %>% dplyr::select(AGE, HB, HBM=M_lms, HBS=S_lms, HBL = L_lms) %>% 
  group_by(AGE) %>% mutate(Q1 = quantile(HB, probs=0.01, na.rm=T),
                           Q99 = quantile(HB, probs=0.99, na.rm=T))

full_HB_lms<-full_HB_lms_pre %>% distinct(AGE, .keep_all=T)

# smooth percentiles to create upper and lower limits of sampling
HBsmoothQ1_1<- predict(smooth.spline(x=full_HB_lms$AGE,y=full_HB_lms$Q1, df = 3))
HBsmoothQ99_1<- predict(smooth.spline(x=full_HB_lms$AGE,y=full_HB_lms$Q99, df = 3))
HBsmoothdf<-data.frame(AGE=HBsmoothQ1_1$x, HBP1 = HBsmoothQ1_1$y, HBP99 = HBsmoothQ99_1$y)

HB_raw_LMS<-full_HB_lms %>% left_join(HBsmoothdf) %>% 
  dplyr::select(AGE, HBM, HBS, HBL, HBP1, HBP99)

# merge dataframes to be used for virtual population generation 
raw_combined_LMS_weight<-rbind(weight_boy_raw_LMS, weight_girl_raw_LMS)
raw_combined_LMS<-left_join(raw_combined_LMS_weight,HB_raw_LMS,by="AGE")

#######################################################################################
# create virtual population: Function 
# adopted from Kitabi et al. (2021)* and modified here 

# * Eliford Kitabi and others, Clinical Infectious Diseases, 
# Volume 73, Issue 5, 1 September 2021, Pages 903-906, 
# https://doi.org/10.1093/cid/ciab149

# Function to sample weight and HB given the LMS parameters
# n  is the required sample size per age
# P1 and P99 are the outer limits (data driven) of the sample, anything outside these limits
# wont enter the population
#######################################################################################

sample_wtnhb<-function(Age, Sex, M, L, S, P1, P99, HBM, HBS,HBL, HBP1, HBP99, n, ...){
  sd   <- abs(ifelse(L==0, S, M^L*L*S)) #sd for weight
  HBsd <- abs(ifelse(HBL==0, HBS, HBM^HBL*HBL*HBS))  #sd for hb
  
  mean   <- ifelse(L==0, log(M), M^L) #mean for weight
  HBmean <- ifelse(HBL==0, log(HBM), HBM^HBL)  #mean for hb
  
  zscore_P99 <- ifelse(L==0, ln(P99/M)/S, (P99^L-M^L)/(S*L*M^L))
  zscore_P1 <- ifelse(L==0, ln(P1/M)/S, (P1^L-M^L)/(S*L*M^L))
  HBzscore_P99 <- ifelse(HBL==0, ln(HBP99/HBM)/HBS, (HBP99^HBL-HBM^HBL)/(HBS*HBL*HBM^HBL))
  HBzscore_P1 <- ifelse(HBL==0, ln(HBP1/HBM)/HBS, (HBP1^HBL-HBM^HBL)/(HBS*HBL*HBM^HBL))
  
  wthb   <- MASS::mvrnorm(n=n*1e2, 
                          mu=c(mean, HBmean), 
                          Sigma = matrix(c(sd^2, 0, 0, HBsd^2), nrow = 2, ncol = 2)) %>% as.data.frame(.)
  # mvrnorm produces samples from the specified multivariate normal dist from the MASS package 
  # produces n*100 samples with means being the mean of weight and mean of hb 
  # Sigma covar matrix with variances for weight and hb in diagonal 
  x      <- wthb[,1] #gets weight 
  HBx    <- wthb[,2] #gets hb 
  
  zscore   = (x-mean)/sd #defn of Z-score for weight 
  HBzscore = (HBx-HBmean)/HBsd #defn of Z-score for hb 
  
  zscore    <- zscore[zscore > zscore_P1 & zscore < zscore_P99] #keeps z-scores that are between 1 and 99 percentile from data 
  HBzscore  <- HBzscore[HBzscore > HBzscore_P1 & HBzscore < HBzscore_P99]
  
  zscore    <- sample(zscore, size = n, replace = TRUE) 
  HBzscore  <- sample(HBzscore, size = n, replace = TRUE)
  #gets random sample of zscore of size n with replacement , replacement in case we need to sample a sample
  #larger than the number of zscores
  
  if(L==0) {wt = M*exp(S*zscore)} else {wt = M*(1+L*S*zscore)^(1/L)}
  if(HBL==0) {hb = HBM*exp(HBS*HBzscore)} else {hb = HBM*(1+HBL*HBS*HBzscore)^(1/HBL)}
  
  df <- data.frame(Age=Age, Sex=Sex, M=M, L=L, S=S, P1=P1, P99=P99, HBM=HBM,
                   HBS=HBS, HBP1=HBP1,HBL=HBL, HBP99=HBP99, WT=wt, HB=hb)
  return(df)
}

#######################################################################################
# creating virtual paediatric population with Weights and HB

raw_no_bin_virtual<-raw_combined_LMS %>% 
  dplyr::select(Age=AGE, Sex, L, M, S, P1, P99, HBM, HBS,HBL, HBP1, HBP99) %>% 
  purrr::pmap_dfr(.f=sample_wtnhb, n=30000)

#######################################################################################
# adding in temperature for the complete hb, weight and temp dataset:

# temperature simulated from a normal distribution with mean 38.2C and standard deviation 2.6C
# calculated from the SMAC dataset
temperature <- rnorm(n=10000000, mean = 38.2, sd = 2.6) 

# taking 1000 rows per rounded weight randomly for simulating AUC: 
# keeping pediatric temperatures below 41 and above 35 
 
Raw_no_bin_paed <- raw_no_bin_virtual  %>% 
  mutate(TEMP = sample(x= temperature[between(temperature, 35, 41)], size=n()),
         SEX = ifelse(Sex==1, "M", "F"),
         rounded_weight = round(WT, digits=0)) %>% filter(between(rounded_weight,3,40)) %>% 
  group_by(rounded_weight) %>% slice_sample(n=1000) %>%
  ungroup()

########################################################################################



