data {
  int<lower=0> n_obs;         // Number of observations (sample-units, s-u)
  int<lower=1> n_site;        // Number of sites
  int<lower=1> n_sample;      // Number of samples (no. occasions * no. sites)
  int<lower=1> n_pred;        // Number of predictor variables
  matrix[n_obs,n_pred] u;     // group predictors (model matrix) 
  array[n_obs] real y;        // Normally-distributed response variable
  array[n_obs] int site_no;   // Site number
  array[n_obs] int samp_no;   // Sample number
}
parameters {
  vector[n_pred] gamma;        // beta parameters of fixed effects in u
  vector[n_site] a_si;         // coefficient of random site effect
  vector[n_sample] a_sa;       // coefficient of random sample effect
  real<lower=0> sigma_si;      //sd of hyperdistribution of a_sis among taxa
  real<lower=0> sigma_sa;      //sd of hyperdistribution of a_sas among taxa
  real<lower=0> sigma;         //sd of mu
}
transformed parameters {
  vector[n_obs] mu;  // Log total count
   for(i in 1:n_obs){
       //The linear model
  mu[i] = a_si[site_no[i]] +  a_sa[samp_no[i]] +  u[i,] * gamma; 
      }
}

model {
  // Priors
   a_si ~ normal(0, sigma_si);
   a_sa ~ normal(0, sigma_sa);
   to_vector(gamma) ~ normal(0,5);
   sigma_si ~ student_t(2, 0, 0.3); //exponential(1); //normal(0.5,2); 
   sigma_sa ~ student_t(2, 0, 0.3); //exponential(1); //normal(0,1); 
   sigma ~ normal(0,1);

 // Likelihood
  for (i in 1 : n_obs) {
    target += normal_lpdf(y[i] | mu[i], sigma);
  }
}


