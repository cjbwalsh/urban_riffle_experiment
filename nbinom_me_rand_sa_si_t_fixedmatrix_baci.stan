
data {
  int<lower=0> n_obs;         // Number of observations (sample-units, s-u)
  int<lower=1> n_taxa;        // Number of taxa
  int<lower=1> n_site;        // Number of sites
  int<lower=1> n_sample;      // Number of samples (combinations of t_no, site_no)
  int<lower=1> n_pred;        // Number of predictor variables
  int<lower=1> n_t;           // Number of sampling occasions
  matrix[n_obs,n_pred] u;     // group predictors (model matrix) 
  array[n_obs,n_taxa] int c;  // Count of each taxon in each subsample
  array[n_obs,n_taxa] real s; // Subsample proportion for each s-u
  array[n_obs] int site_no;   // Site number
  array[n_obs] int samp_no;   // Sample number (vectorised [t_no, site_no])
  array[n_obs] int t_no;      // Sampling occasion number
}
parameters {
  vector[n_pred] mu_gamma;         //means of beta parameter hyperdistributions          
  matrix[n_pred,n_taxa] gamma;        // beta parameters of fixed effects in u
  matrix[n_site,n_taxa] a_s_raw;      // raw coefficient of random site effect
  matrix[n_t,n_taxa] a_t_raw;         // raw oefficient of time effect 
  matrix[n_sample,n_taxa] a_st_raw;   // raw coefficient of random sample effect 
                                      //(1 sample per site_no per t_no)
  // see reparameterizations below for explanation of these raw coefficients
  vector<lower=0>[n_taxa] sigma_s; //rsd of hyperdistribution of a_ss among taxa
  vector<lower=0>[n_taxa] sigma_t; //sd of hyperdistribution of a_ts among taxa
  vector<lower=0>[n_taxa] sigma_st;//sd of hyperdistribution of a_sts among taxa
  vector<lower=0>[n_taxa] phi;     //dispersion parameter for neg-binomial distn
  corr_matrix[n_pred] Omega;       // Hyperprior correlation matrix among taxa
  vector<lower=0>[n_pred] tau;     // Scale for correlation matrix
}
transformed parameters {
  matrix[n_obs,n_taxa] log_lambda;  // Log total count
  matrix[n_site,n_taxa] a_s;        // coefficient of random site effect
  matrix[n_t,n_taxa] a_t;           // Coefficient of random time effect 
  matrix[n_sample,n_taxa] a_st;     // coefficient of random sample effect
                         // (sample = 3-4 s-us from each site on each occasion)
  for(j in 1:n_taxa){
  a_s[,j] =  sigma_s[j] * a_s_raw[,j];  
  a_t[,j] =  sigma_t[j] * a_t_raw[,j];  
  a_st[,j] =  sigma_st[j] * a_st_raw[,j];  
  }
// with (e.g.) a_s_raw ~ std_normal(), this implies a_s ~ normal(0, sigma_s)
// See https://mc-stan.org/docs/stan-users-guide/reparameterization.html

  for(i in 1:n_obs){
     for(j in 1:n_taxa){
       //The linear model
      log_lambda[i,j] = a_s[site_no[i],j] +  a_t[t_no[i],j] +  
                       a_st[samp_no[i],j] +  u[i,] * gamma[,j]; 
      }
      }
}

model {
  // Priors
   mu_gamma ~ normal(0,5);
   to_vector(a_s_raw) ~ std_normal();
   to_vector(a_t_raw) ~ std_normal();
   to_vector(a_st_raw) ~ std_normal();
   sigma_s ~ normal(0,1);
   sigma_t ~ normal(0,1);
   sigma_st ~ normal(0,1);
   phi ~ normal(1,1);
   tau ~ exponential(1);
   Omega ~ lkj_corr( 2 );  
   
   // estimation of correlated beta parameters (assembled in matrix gamma)
   for(i in 1:n_taxa){
   target += multi_normal_prec_lpdf(gamma[,i] | mu_gamma, quad_form_diag(Omega, 
                                                                         tau) );
     }
  // Likelihood
  for (i in 1 : n_obs) {
    for(j in 1:n_taxa){
    target += neg_binomial_2_log_lpmf(c[i,j] | log_lambda[i,j] + log(s[i,j]), 
                                      phi[j]);
  // This parameterization adds the marginal log-binomial-probability 
  // resulting from subsampling error to the marginal negative-binomial 
  // probability of the linear model.  It is equivalent to a (50 times) slower 
  // parameterization modelling the marginal binomial and negative-binomial
  // probabilities separately, by looping through all feasible total counts
  // given each count and subsample proportion.
   }
  }
}

// generated quantities {
//   // log-likelihood only used for model comparisons during model development.
//   // The most complex model considered has been used, 
//   // so model comparisons not reported
//   array[n_obs,n_taxa] real log_lik;  
//  for (i in 1 : n_obs) {
//    for(j in 1 : n_taxa){
//    log_lik[i,j] = neg_binomial_2_log_lpmf(c[i,j] | log_lambda[i,j] + 
//                     log(s[i,j]), phi[j]);
//  }
//  }
//  }


