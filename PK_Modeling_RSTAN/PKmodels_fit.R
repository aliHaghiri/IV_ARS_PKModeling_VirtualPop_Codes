
rm(list = ls()) # clear memory
graphics.off()  # clear graphics
cat("\014")     # clear console

# Load R libraries
library(rstan)
library(cmdstanr)
library(dplyr)
library(parallel)

# Setting the working directory
setwd("XXX")

# Load Stan data
fname = "data-stan-REACH-Dataset"
load(paste0("./StanData/", fname, ".RData"))

# Run Stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Excluding TIME = 0 because ode_rk45 gets these parameters separately
dataStan$t_PPDHA[dataStan$t_PPDHA == 0] = 1e-3

# PK Parameters Initialization

thetaPop_init = c(rnorm(1, log(100), 0.1),  # Cl_ARS
                  rnorm(1, log(40), 0.1),   # V_ARS
                  rnorm(1, log(30), 0.1),   # Cl_DHA
                  rnorm(1, log(15), 0.1),   # V_DHA
                  rnorm(1, log(50), 0.1),   # Q_ARS
                  rnorm(1, log(40), 0.1),   # Vp_ARS
                  rnorm(1, log(30), 0.1),   # Q_DHA
                  rnorm(1, log(15), 0.1))   # Vp_DHA

init_2 <- function(){
  list(thetaPop = thetaPop_init,
       beta_Hb = 0,
       L = matrix(runif(64,0.05,0.15),ncol=8) ,
       omega = matrix(runif(8,0.05,0.15),ncol=8),
       etaStd = matrix(runif(640,0.05,0.15),ncol=80),
       sigmaARS = abs(rnorm(1, 0, 0.2)),
       sigmaDHA = abs(rnorm(1, 0, 0.2))
  )
}

seed = 12
iter = 4000 
adapt_delta = 0.80

mod <- cmdstan_model("PKmodelsODE_REACH.stan")

fit_stan_cmd <- mod$sample(
  data = dataStan, 
  seed = seed, 
  chains = 4, 
  parallel_chains = 4, 
  refresh = iter/10,
  iter_warmup=iter/2,
  iter_sampling=iter/2,
  adapt_delta = adapt_delta,
  init =init_2, 
  thin = 1
)

stan.fsave.name = paste0("_seed",
seed, "-iter", iter, "-adapt_delta", adapt_delta, "_rk45")


fit_stan_cmd$save_object(paste0(fit.dir, "/Results/Fits/", 
                            fname, stan.fsave.name, ".RDS"))

fit_stan <- rstan::read_stan_csv(fit_stan_cmd$output_files())

save(fit_stan, dataStan, file = paste0("fit-",fname, stan.fsave.name, ".RData"))


