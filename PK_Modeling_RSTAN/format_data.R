
rm(list = ls()) # clear memory
graphics.off()  # clear graphics
cat("\014")     # clear console


# Load R libraries
library(dplyr)
library(magrittr)
library(ggplot2)

# Setting the working directory
setwd("XXX")

# load clinical data
df.ARS = read.csv("./Data/REACH_Dataset_ARS_cleaned.csv")
df.DHA = read.csv("./Data/REACH_Dataset_DHA_cleaned.csv")

# FUNCTION to format data for Stan ----------------------------------------
format_data = function(data.keep.ARS, data.keep.DHA, nPar, lwrLim, uprLim)
{
  # Creating the time-index variable
  nARS = data.keep.ARS$TIME.ARS
  nDHA = data.keep.DHA$TIME.DHA
  # Generate all the time-index values (i.e. those when measurement may not carried out too)
  iObsARS = data.keep.ARS %>% group_by(id.ARS) %>% mutate(iObsARS = row_number()) %>% ungroup() %>% dplyr::select(iObsARS) %>% as.matrix() %>% as.double()
  iObsDHA = data.keep.DHA %>% group_by(id.DHA) %>% mutate(iObsDHA = row_number()) %>% ungroup() %>% dplyr::select(iObsDHA) %>% as.matrix() %>% as.double()
  
  nObsARS = length(with(data.keep.ARS, nARS))
  nObsDHA = length(with(data.keep.DHA, nDHA))
  
  nSubjects = data.keep.ARS$id.ARS %>% unique(.) %>% length(.)
  
  iObsMaxArrayARS = tapply(iObsARS, data.keep.ARS$id.ARS, max)
  dimnames(iObsMaxArrayARS) = NULL
  iObsMaxARS = c(iObsMaxArrayARS)
  
  tMaxArrayARS = tapply(nARS, data.keep.ARS$id.ARS, max)
  dimnames(tMaxArrayARS) = NULL
  tMaxARS = c(tMaxArrayARS)
  
  iObsMaxArrayDHA = tapply(iObsDHA, data.keep.DHA$id.DHA, max)
  dimnames(iObsMaxArrayDHA) = NULL
  iObsMaxDHA = c(iObsMaxArrayDHA)
  
  tMaxArrayDHA = tapply(nDHA, data.keep.DHA$id.DHA, max)
  dimnames(tMaxArrayDHA) = NULL
  tMaxDHA = c(tMaxArrayDHA)
  
  ## Start and end indices of samples at actually measured times
  startARS = (1:nObsARS)[!duplicated(data.keep.ARS$id.ARS)]
  endARS = c(startARS[-1] - 1, nObsARS)
  
  startDHA = (1:nObsDHA)[!duplicated(data.keep.DHA$id.DHA)]
  endDHA = c(startDHA[-1] - 1, nObsDHA)
  
  data.stan = 
    list(idARS = data.keep.ARS$id.ARS,
         idDHA = data.keep.DHA$id.DHA,
         timeARS = nARS,
         timeDHA = nDHA,
         tMaxARS = tMaxARS,
         tMaxDHA = tMaxDHA,
         nSubjects = nSubjects,
         nObsARS = nObsARS,
         nObsDHA = nObsDHA,
         iObsARS = iObsARS,
         iObsDHA = iObsDHA,
         iObsMaxARS = iObsMaxARS,
         iObsMaxDHA = iObsMaxDHA,
         startARS = startARS,
         endARS = endARS,
         startDHA = startDHA,
         endDHA = endDHA,
         ConcARS = data.keep.ARS$Conc.ARS,
         ConcDHA = data.keep.DHA$Conc.DHA,
         dose = data.keep.ARS %>% group_by(id.ARS) %>% filter(row_number() == 1) %>% ungroup() %>% dplyr::select(DOSE_nmol) %>% as.matrix() %>% as.numeric(), 
         tDose = 0,
         nDose = 1,
         nPar = nPar,
         a = lwrLim, 
         b = uprLim,
         BQLARS = df.ARS$BQL.ARS,
         BQLDHA = df.DHA$BQL.DHA,
         LLOQARS = 3.1, # nM
         LLOQDHA = 6.9,  # nM
         BWT = data.keep.ARS %>% group_by(id.ARS) %>% filter(row_number() == 1) %>% ungroup() %>% dplyr::select(WT) %>% as.matrix() %>% as.numeric(),
         BWT0 = 17.93, # mean of BW for all considered IDs
         ntPPARS = 101, # Number of time points at which we want to have model simulation + error (posterior predictive)
         ntPPDHA = 101, # Number of time points at which we want to have model simulation + error (posterior predictive)
         Tsim = Tsim, # maximum time for getting posterior predictive
         t_PPARS = seq(0, Tsim, Tsim/100),  # time points at which we want to have model simulation + error (posterior predictive)
         t_PPDHA = seq(0, Tsim, Tsim/100),  # Number of time points at which we want to have model simulation + error (posterior predictive)
         Hb = data.keep.ARS %>% group_by(id.ARS) %>% filter(row_number() == 1) %>% ungroup() %>% dplyr::select(HB) %>% as.matrix() %>% as.numeric(),
         Hb0 = 7.91 # mean of HB for all considered IDs
         )
  
  return(data.stan)
}

# Number of PK parameters:

# Cl_ARS: clearance of parent drug 
# V_ARS: Volume of the central compartment for parent drug
# Cl_DHA: clearance of metabolite 
# V_DHA: Volume of the central compartment for metabolite drug
# Q_ARS: inter-compartmental clearance of parent drug 
# Vp_ARS: Volume of the peripheral compartment for parent drug
# Q_DHA: inter-compartmental clearance of metabolite 
# Vp_DHA: Volume of the peripheral compartment for metabolite drug

nPar = 8 # see above

# Create Stan data 

dataStan = format_data(df.ARS, df.DHA, nPar = nPar)

# Take a look at formatted data
str(dataStan)

# SAVE 
save(dataStan, file = "./StanData/data-stan-REACH-Dataset.RData")
