// ===================== The function section defines a 2-Compartment Parent-metabolite PK model for ARS-DHA ========================

functions { 
    
  real C_ARS(real t, real D, real Clc, real Vc, real Q1, real Vp2, real BWTi, real BWT0){
    real Cc; 
    real Clc_s; 
    real Vc_s;   
    real Q1_s;   
    real Vp2_s;
    real beta;
    real alpha;
    real AA;
    real BB;
    
    Clc_s = (Clc); //parameters are already scaled in transformed parameters{} chunk by * (BWTi/BWT0)^0.75
    Vc_s = (Vc); // ditto for Vc by * (BWTi/BWT0)^1;
    Q1_s = Q1;
    Vp2_s = Vp2;
    
    beta = 0.5*(Q1_s/Vc_s+Q1_s/Vp2_s+Clc_s/Vc_s-sqrt((Q1_s/Vc_s+Q1_s/Vp2_s+Clc_s/Vc_s)^2-4*(Q1_s/Vp2_s)*(Clc_s/Vc_s)));
    alpha = (Q1_s/Vp2_s)*(Clc_s/Vc_s)/beta;
    AA = (alpha-(Q1_s/Vp2_s))/(Vc_s*(alpha-beta));
    BB = (beta-(Q1_s/Vp2_s))/(Vc_s*(beta-alpha));

    Cc = D*(AA*exp(-alpha*t)+BB*exp(-beta*t));

    return(Cc); 
  }
  
  vector CDHAder(real t,        // time
  vector y,      // state
  vector theta,  // parameters
  vector x_r    // data (real)
  ) {   
  
  real Cm;
  real Clc;
  real Vc;
  real Clm;
  real Vm;
  real Clc_s;
  real Vc_s;
  real Clm_s;
  real Vm_s;
  real BWTi;
  real BWT0;
  real D;
  
  real Q1;   
  real Vp2;
  real Q2;   
  real Vp4;
  real Q1_s;   
  real Vp2_s;
  real Q2_s;   
  real Vp4_s;
  real beta;
  real alpha;
  real AA;
  real BB;  
  
  real Cc_t;
  vector[2] dydt;
  
  Clc = theta[1];  // Cl_ARS
  Vc = theta[2];   // V_ARS
  
  Clm = theta[3];  // Cl_DHA
  Vm = theta[4];   // V_DHA
   
  Q1 = theta[5];   // Q_ARS
  Vp2 = theta[6];  // Vp_ARS
  
  Q2 = theta[7];   // Q_DHA
  Vp4 = theta[8];  // Vp_DHA
  
  D = x_r[1];
  BWTi = x_r[2];
  BWT0 = x_r[3];
  
  //parameters are already (allometrically) scaled in transformed parameters{} chunk
  Clc_s = (Clc);// * (BWTi/BWT0)^0.75;
  Vc_s = (Vc);// * (BWTi/BWT0)^1;
  Clm_s = (Clm);// * (BWTi/BWT0)^0.75;
  Vm_s = (Vm);// * (BWTi/BWT0)^1;
  
  Q1_s = Q1;
  Vp2_s = Vp2;
  Q2_s = Q2;
  Vp4_s = Vp4;
    
  beta = 0.5*(Q1_s/Vc_s+Q1_s/Vp2_s+Clc_s/Vc_s-sqrt((Q1_s/Vc_s+Q1_s/Vp2_s+Clc_s/Vc_s)^2-4*(Q1_s/Vp2_s)*(Clc_s/Vc_s)));
  alpha = (Q1_s/Vp2_s)*(Clc_s/Vc_s)/beta;
  AA = (alpha-(Q1_s/Vp2_s))/(Vc_s*(alpha-beta));
  BB = (beta-(Q1_s/Vp2_s))/(Vc_s*(beta-alpha));

  Cc_t = D*(AA*exp(-alpha*t)+BB*exp(-beta*t));

  dydt[1] = Clc_s/Vm_s * Cc_t - (Clm_s/Vm_s+Q2_s/Vm_s) * y[1] + Q2_s/Vm_s * y[2];
  dydt[2] = Q2_s/Vp4_s * y[1] - Q2_s/Vp4_s * y[2];

  return dydt;
  }
  
  vector PK_model_ARS(int iObsMax, real[] time, real D, real Clc, real Vc, real Q1, real Vp2,  real BWTi, real BWT0){
    vector[iObsMax] CARS;

    for(t in 1:iObsMax){
      CARS[t] = C_ARS(time[t], D, Clc, Vc, Q1, Vp2, BWTi, BWT0);
      
      // Because log(0) is undefined, we'll get errors
      // Sometimes CARS becomes inf, when it becomes really small
      if(CARS[t] <= 1e-5 || is_inf(CARS[t])){CARS[t] = 1e-5;}
    }
    
    return CARS;
  }
  
} 

data { 
  
  int nObsARS;
  int nObsDHA;
  int idARS[nObsARS];
  int idDHA[nObsDHA];
  
  array[nObsARS] real timeARS;
  array[nObsDHA] real timeDHA; 
  int nSubjects;
  
  int iObsARS[nObsARS];
  int iObsDHA[nObsDHA];
  int iObsMaxARS[nSubjects];
  int iObsMaxDHA[nSubjects];
  int startARS[nSubjects];
  int endARS[nSubjects];
  int startDHA[nSubjects];
  int endDHA[nSubjects];
  real ConcARS[nObsARS];
  real ConcDHA[nObsDHA]; 

  
  vector[nSubjects] dose;
  real tDose;
  real nDose;
  int nPar;
  vector[nPar] a;
  vector[nPar] b;
  real BQLARS[nObsARS];
  real BQLDHA[nObsDHA]; 
  real LLOQARS;
  real LLOQDHA;
  
  vector[nSubjects] BWT;
  real BWT0;
  
  int ntPPARS;
  int ntPPDHA;
  real Tsim;
  array[ntPPARS] real t_PPARS;
  array[ntPPDHA] real t_PPDHA;
  
  vector[nSubjects] Hb;
  real Hb0;
  
}

transformed data {
  
  vector[2] C_DHA_init; // intial DHA concentrations for both central and peripheral compartments
  real t0;
  matrix[nSubjects,3]  x_r;
  
  
  C_DHA_init[1] = 0;
  C_DHA_init[2] = 0;
  t0 = 0;
  
  x_r[, 1] = (dose);
  x_r[, 2] = (BWT);
  x_r[, 3] = (rep_vector(BWT0, nSubjects));
  
}

parameters {
  // thetaPop can be negative here because we are sampling on the log scale
  vector[nPar] thetaPop; //<lower = 0>
  real<lower = 0> sigmaARS;
  real<lower = 0> sigmaDHA;
  
  // Inter-Individual variability
  cholesky_factor_corr[nPar] L;
  vector<lower = 0>[nPar] omega;  // the scale of variation of coefficients 
  matrix[nPar, nSubjects] etaStd;
  
  real beta_Hb;
}

transformed parameters {
  
  // thetaInd can be negative here because we are sampling on the log scale
  matrix[nSubjects, nPar] thetaInd;//<lower=0>
  vector<lower=0>[nObsARS] CPredARS;
  vector<lower=0>[2] CPredDHA_2comp[nObsDHA];
  vector<lower=0>[1] CPredDHA[nObsDHA];

  matrix[nPar, nPar] Cor;
  matrix[nPar, nPar] Cov;
  
  vector<lower=0>[8] theta;
  
  real<lower=0> Clc_s;
  real<lower=0> Vc_s;
  
  real<lower=0> Clm_s; 
  real<lower=0> Vm_s;
  
  real<lower=0> Q1_s;
  real<lower=0> Vp2_s;
  
  real<lower=0> Q2_s;
  real<lower=0> Vp4_s;
  
  
  // Individual parameters
  thetaInd = (rep_matrix(thetaPop, nSubjects) + diag_pre_multiply(omega, L) * etaStd)'; // Ali 24-03-2022 diag_pre_multiply(omega, L) * etaStd  
  
  // Calculate correlation and covariance matrices
  Cor = L * L';
  Cov = quad_form_diag(Cor,omega);
  
  // PK Calculations
  for (i in 1:nSubjects) {
    
    // Because parameters are sampled on the log scale, they need to be exponentiated
    Clc_s = exp(thetaInd[i, 1]) * (BWT[i]/BWT0)^0.75;
    Vc_s = exp(thetaInd[i, 2]) * (BWT[i]/BWT0)^1;
    Clm_s = exp(thetaInd[i, 3]) * exp(beta_Hb * (Hb[i] - Hb0)) * (BWT[i]/BWT0)^0.75;
    Vm_s = exp(thetaInd[i, 4]) * (BWT[i]/BWT0)^1;
    
    Q1_s = exp(thetaInd[i, 5]) * (BWT[i]/BWT0)^0.75;
    Vp2_s = exp(thetaInd[i, 6]) * (BWT[i]/BWT0)^1;    
    
    Q2_s = exp(thetaInd[i, 7]) * (BWT[i]/BWT0)^0.75;
    Vp4_s = exp(thetaInd[i, 8]) * (BWT[i]/BWT0)^1;

    CPredARS[startARS[i]:endARS[i]] = PK_model_ARS(
      iObsMaxARS[i],
      timeARS[startARS[i]:endARS[i]],
      dose[i],
      Clc_s,
      Vc_s,
      Q1_s,
      Vp2_s,
      BWT[i],
      BWT0
      );
      
      theta[1] = Clc_s;
      theta[2] = Vc_s;
      theta[3] = Clm_s;
      theta[4] = Vm_s;
      
      theta[5] = Q1_s;
      theta[6] = Vp2_s;
      theta[7] = Q2_s;
      theta[8] = Vp4_s;
      
      CPredDHA_2comp[startDHA[i]:endDHA[i]] = ode_rk45(CDHAder, C_DHA_init, t0, timeDHA[startDHA[i]:endDHA[i]], theta, to_vector(x_r[i]));
      
      for(j in startDHA[i]:endDHA[i]){
        // Make sure CPredDHA != 0 because log(CPredDHA)=nan, hence issues with the likelihood
        if((CPredDHA_2comp[j][1]) <= 1e-5) CPredDHA_2comp[j][1] = 1e-5;
        if((CPredDHA_2comp[j][2]) <= 1e-5) CPredDHA_2comp[j][2] = 1e-5;
        CPredDHA[j][1] = CPredDHA_2comp[j][1]; // central compartment DHA concentration 

      }
      
  }
}  
  
model {
  // Population average
  thetaPop ~ normal(log(100), log(100));

  // Effect of Hb on Cl_DHA
  beta_Hb ~ normal(0, 0.25);
  // Inter-individual variability 
  L ~ lkj_corr_cholesky(2);
  omega ~ normal(0, 1);
  to_vector(etaStd) ~ normal(0, 1);
  
  // Likelihood function 
  // - M3 methods for below quantification limit data (BQL)
  sigmaARS ~ cauchy(0, 2);
  sigmaDHA ~ cauchy(0, 2);
  
  for (i in 1:nObsARS) {
      if (BQLARS[i] == 0) {
        ConcARS[i] ~ lognormal(log(CPredARS[i]), sigmaARS);
      } 
      else {
        target += normal_lcdf(log(LLOQARS) | log(CPredARS[i]), sigmaARS);
      }
      // }
  }
  for (i in 1:nObsDHA) {
      if (BQLDHA[i] == 0) {
        ConcDHA[i] ~ lognormal(log(CPredDHA[i]), sigmaDHA);
      }
      else {
        target += normal_lcdf(log(LLOQDHA) | log(CPredDHA[i]), sigmaDHA);
      }
  }
}

generated quantities {
  vector[ntPPARS*nSubjects] CPredARSErr;
  vector[ntPPDHA*nSubjects] CPredDHAErr;
  vector[ntPPARS*nSubjects] CPredARSErrmean;
  vector[2] CPredDHAErrmean_2comp[ntPPDHA*nSubjects];
  vector[1] CPredDHAErrmean[ntPPDHA*nSubjects];

  vector<lower=0>[8] theta_2;
  
  vector<lower=0>[nObsARS] CPredARSerr;
  array [nObsDHA] real <lower=0> CPredDHAerr;
  
  // Population profiles
  for (i in 1:nSubjects){
    
    theta_2[1] = exp(thetaInd[i, 1]) * (BWT[i]/BWT0)^0.75;
    theta_2[2] = exp(thetaInd[i, 2]) * (BWT[i]/BWT0)^1;
    theta_2[3] = exp(thetaInd[i, 3]) * exp(beta_Hb * (Hb[i] - Hb0)) * (BWT[i]/BWT0)^0.75;
    theta_2[4] = exp(thetaInd[i, 4]) * (BWT[i]/BWT0)^1;
    
    theta_2[5] = exp(thetaInd[i, 5]) * (BWT[i]/BWT0)^0.75;
    theta_2[6] = exp(thetaInd[i, 6]) * (BWT[i]/BWT0)^1;    
    
    theta_2[7] = exp(thetaInd[i, 7]) * (BWT[i]/BWT0)^0.75;
    theta_2[8] = exp(thetaInd[i, 8]) * (BWT[i]/BWT0)^1;
    
    CPredARSErrmean[((i-1)*ntPPARS + 1) : (i*ntPPARS)] = PK_model_ARS(
      ntPPARS,
      t_PPARS,
      dose[i],
      theta_2[1],
      theta_2[2],
      theta_2[5],
      theta_2[6],
      BWT[i],
      BWT0
      );
      
      CPredDHAErrmean_2comp[((i-1)*ntPPDHA + 1) : (i*ntPPDHA)] = ode_rk45(CDHAder, C_DHA_init, t0, t_PPDHA, theta_2, to_vector(x_r[i]));
      
       for(j in ((i-1)*ntPPDHA + 1) : (i*ntPPDHA) ){

        CPredDHAErrmean[j][1] = CPredDHAErrmean_2comp[j][1]; // central compartment DHA concentration 

      }
      
  }
  
  // Posterior predictives
  // Going over a sequence of time points equally spaced with a high resolution
  for (i in 1:(ntPPARS*nSubjects)) {
    CPredARSErr[i] = exp(normal_rng(log(CPredARSErrmean[i]), sigmaARS));
  }
  for (i in 1:(ntPPDHA*nSubjects)) {
    CPredDHAErr[i] = exp(normal_rng(log(CPredDHAErrmean[i][1]), sigmaDHA));
  }
  // Going over a sequence of time points when measurements were carried out in the subjects
  for (i in 1:nObsARS) {
    CPredARSerr[i] = exp(normal_rng(log(CPredARS[i]), sigmaARS));
  }
  for (i in 1:nObsDHA) {
    CPredDHAerr[i] = exp(normal_rng(log(CPredDHA[i][1]), sigmaDHA));
  }
}

