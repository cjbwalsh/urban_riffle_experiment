data {
  int<lower=0> n_obs;         // Number of observations (sample-units, s-u)
  int<lower=1> n_site;        // Number of sites
  int<lower=1> n_sample;      // Number of samples
  int<lower=1> n_pred;        // Number of predictor variables
  int<lower=1> n_t;           // Number of sampling occasions
  matrix[n_obs,n_pred] u;     // group predictors (model matrix) 
  array[n_obs] real y;        // Normally-distributed response variable
  array[n_obs] int site_no;   // Site number
  array[n_obs] int samp_no;   // Sample number (vectorised [t_no, site_no])
  array[n_obs] int t_no;      // sampling occasion number
}
parameters {
  vector[n_pred] gamma;   // beta parameters of fixed effects in u
  vector[n_site] a_s_raw;     // raw coefficient of random site effect 
  vector[n_t] a_t_raw;        // raw coefficient of time effect 
  vector[n_sample] a_st_raw;  // raw coefficient of random sample effect (1 sample per site_no per t_no)
  // see reparameterizations below for explanationof these raw coefficients
  real<lower=0> sigma_s; //sd of hyperdistribution of a_sis among taxa
  real<lower=0> sigma_t;       //sd of hyperdistribution of a_ts among taxa
  real<lower=0> sigma_st;      //sd of hyperdistribution of a_sas among taxa
  real<lower=0> sigma;         //sd of mu
}
transformed parameters {
  vector[n_obs] mu;  // Log total count
  vector[n_site] a_s;     // coefficient of random site effect
  vector[n_t] a_t;        // Coefficient of time effect 
  vector[n_sample] a_st;     // coefficient of random site effect
  
  a_s =  sigma_s * a_s_raw;  
  a_t =  sigma_t * a_t_raw;  
  a_st =  sigma_st * a_st_raw;  
// with (e.g.) a_s_raw ~ std_normal(), this implies a_s ~ normal(0, sigma_s)
// See https://mc-stan.org/docs/stan-users-guide/reparameterization.html

   for(i in 1:n_obs){
       //The linear model
  mu[i] = a_s[site_no[i]] + a_t[t_no[i]]  +  a_st[samp_no[i]] +  
          u[i,] * gamma; 
      }
}

model {
  // Priors
   a_s_raw ~ std_normal();
   a_t_raw ~ std_normal();   ///1
   a_st_raw ~ std_normal();
   to_vector(gamma) ~ normal(0,5);
   sigma_s ~ normal(0,1);  //student_t(2, 0, 0.1); //uniform(.001,10); //exponential(1); //normal(0.5,2); 
   sigma_t ~ normal(0,1);  //student_t(2, 0, 0.1); //uniform(.001,10); //exponential(1); //normal(0.5,2);  
   sigma_st ~ normal(0,1);  //student_t(2, 0, 0.1); //uniform(.001,10); //exponential(1); //
   sigma ~ normal(0,1); //exponential(1); //uniform(.001,10); //

 // Likelihood
  for (i in 1 : n_obs) {
    target += normal_lpdf(y[i] | mu[i], sigma);
  }
}


