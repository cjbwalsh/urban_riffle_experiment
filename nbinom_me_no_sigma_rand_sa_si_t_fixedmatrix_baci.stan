
data {
  int<lower=0> n_obs;         // Number of observations (sample-units, s-u)
  int<lower=1> n_taxa;        // Number of taxa
  int<lower=1> n_site;        // Number of sites
  int<lower=1> n_sample;      // Number of samples
  int<lower=1> n_pred;        // Number of predictor variables
  int<lower=1> n_t;           // Number of sampling occasions
  matrix[n_obs,n_pred] u;     // group predictors (model matrix) 
  array[n_obs,n_taxa] int c;  // Counts of species in each subsample
  array[n_obs,n_taxa] real s; // Subsample proportion for each observation
  array[n_obs] int site_no;   // Site number
  array[n_obs] int samp_no;   // Sample number
  array[n_obs] int t_no;      // sampling occasion number
}
parameters {
  vector[n_pred] mu_gamma;
  matrix[n_pred,n_taxa] gamma;        // beta parameters of fixed effects in u
  matrix[n_site,n_taxa] a_si;         // coefficient of random site effect
  matrix[n_sample,n_taxa] a_sa;       // coefficient of random sample effect
  matrix[n_t,n_taxa] a_t;             // Coefficient of time effect 
  matrix[n_obs,n_taxa] eps_raw;
  // raw s-u variability (see reparameterization below)   
  vector<lower=0>[n_taxa] sigma_eps; // sd of random s-u error for each taxon
  vector<lower=0>[n_taxa] phi;       // dispersion parameter for neg-binomial distribution
  corr_matrix[n_pred] Omega;         // Hyperprior correlation matrix among taxa
  vector<lower=0>[n_pred] tau;       // Scale for correlation matrix
}
transformed parameters {
  matrix[n_obs,n_taxa] eps;         // Abundance noise
  matrix[n_obs,n_taxa] log_lambda;  // Log total count

//  for(i in 1:n_obs){
    for(j in 1:n_taxa){
      eps[,j] =  sigma_eps[j] * eps_raw[,j];  
   // with eps_raw ~ std_normal(), this implies eps ~ normal(0, sigma_eps)
   // See https://mc-stan.org/docs/stan-users-guide/reparameterization.html
 }
//  }
  
  for(i in 1:n_obs){
     for(j in 1:n_taxa){
       //The linear model
      log_lambda[i,j] = a_si[site_no[i],j] +  a_sa[samp_no[i],j] +  
                        a_t[t_no[i],j] +  u[i,] * gamma[,j] + eps[i,j]; 
      }
      }
}

model {
  // Priors
   mu_gamma ~ normal(0,5);
   to_vector(a_si) ~ normal(0,1);
   to_vector(a_sa) ~ normal(0,1);
   to_vector(a_t) ~ normal(0,1);
   to_vector(eps_raw) ~ std_normal();
   sigma_eps ~ normal(1,1); // sampled poorly if normal(0,1)
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
  // given a each count and subsample proportion.
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

