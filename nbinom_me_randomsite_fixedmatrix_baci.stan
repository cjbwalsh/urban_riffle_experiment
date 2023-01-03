
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
  array[n_obs] int site_no;   // Site number (integer)
  array[n_obs] int samp_no;   // Sample number (integer)
  array[n_obs] int t;         // sampling occasion number (integer)
}
parameters {
  vector[n_pred] mu_gamma;
  matrix[n_pred,n_taxa] gamma;           // beta parameters of fixed effects in u
  matrix[n_site,n_taxa] a_si;            // coefficient of random site effect
  matrix[n_sample,n_taxa] a_sa;          // coefficient of random sample effect
  matrix[n_t,n_taxa] a_t;     // coefficient of random sampling occasion effect
  matrix<upper=5>[n_obs,n_taxa] epi_raw; // raw s-u error term (upper constraint for stability)
  vector<lower=0>[n_taxa] phi;           // dispersion parameter for each taxon
  vector<lower=0>[n_taxa] sd_lam;        // sd of random s-u error
  corr_matrix[n_pred] Omega;             // Hyperprior correlation matrix among taxa
  vector<lower=0>[n_pred] tau;           // Scale for correlation matrix
}
transformed parameters {
  matrix[n_obs,n_taxa] log_lambda;     // Log total count
  matrix[n_obs,n_taxa] epi;            // Abundance noise
   for(i in 1:n_obs){
    for(j in 1:n_taxa){
      epi[i,j] =  sd_lam[j] * epi_raw[i,j];  //stable reparameterization for epi
  }
  }
  for(i in 1:n_obs){
     for(j in 1:n_taxa){
       //The linear model
      log_lambda[i,j] = a_si[site_no[i],j] +  a_sa[samp_no[i],j] +  
                        a_si[t[i],j] +  u[i,] * gamma[,j] + epi[i,j]; 
      }
      }
}
model {
  // Priors
   mu_gamma ~ normal(0 , 5);
   to_vector(a_si) ~ normal(0,2);
   to_vector(a_sa) ~ normal(0,2);
   to_vector(a_t) ~ normal(0,2);
   to_vector(epi_raw) ~ normal(0, 1);
   sd_lam ~ normal(1, 1);
   phi ~ normal(1,1);
   tau ~ exponential(1);
   Omega ~ lkj_corr( 2 );  
   
   // estimation of correlated beta parameters (assembled in matrix gamma)
   for(i in 1:n_taxa){
   target += multi_normal_prec_lpdf(gamma[,i] | mu_gamma , quad_form_diag(Omega , tau) );
     }
  // Likelihood
  for (i in 1 : n_obs) {
    for(j in 1:n_taxa){
  // Note this parameterization adds the marginal log-binomial-probability 
  // resulting from subsampling error to the marginal negative-binomial 
  // probability of the linear model.  It is equivalent to a slower (x ~50!)
  // parameterization modelling the marginal binomial and negative-binomial
  // probabilities separately, by looping through all feasible total counts
  // given a each count and subsample proportion.
    target += neg_binomial_2_log_lpmf(c[i,j] | log_lambda[i,j] + log(s[i,j]), phi[j]);
   }
  }
}

generated quantities {
  // log-likelihood only used for model comparisons during model development.
  // The most complex model considered has been used, so model comparisons not reported
  array[n_obs,n_taxa] real log_lik;  
 for (i in 1 : n_obs) {
   for(j in 1 : n_taxa){
   log_lik[i,j] = neg_binomial_2_log_lpmf(c[i,j] | log_lambda[i,j] + log(s[i,j]), phi[j]);
 }
 }
 }


