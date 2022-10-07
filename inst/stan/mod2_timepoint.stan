data {
  //data to fit
  int N; //number of observations
  int N_t; //number of distinct time points
  array[N] int rank_t; //rank of the time point related to the observations
  array[N] real s; //number of samples by pool
  array[N] int n; //number of pools
  array[N] int k; //number of positive pools

  //data sensitivity specificity
  int<lower = 0> J_sens;
  array[J_sens] int <lower = 0> y_sens;
  array[J_sens] int <lower = 0> n_sens;

  //hyperparameters of the priors
  array[2] real p_sens;
  array[2] real p_intercept;
  real p_phi;
  real spec;

  int inference;
}


parameters {
  array[N_t] real logit_prev_f; //logit-transformed prevalence (done so that we have same priors as for GP model)
  real <lower=0,upper=1> sens; //sensitivity
  real <lower=0>phi; //overdispersion
}

transformed parameters {
  real kappa = phi+2; //overdispersion
  array[N_t] real <lower=0,upper=1> prev_f = inv_logit(logit_prev_f); //true prevalence
  array[N] real <lower=0,upper=1> pool_pos;

  for(i in 1:N){
    pool_pos[i] = 1.0 - pow(1.0 - (prev_f[rank_t[i]] * sens + (1.0-prev_f[rank_t[i]]) * (1.0-spec)),s[i]);
  }
}

model {
  // priors
  logit_prev_f ~ normal(p_intercept[1], p_intercept[2]);
  phi ~ exponential(p_phi);
  sens ~ beta(p_sens[1], p_sens[2]);

  if(inference==1){
    //Pool positivity
    for(i in 1:N){
      target += beta_binomial_lpmf(k[i] | n[i], kappa * pool_pos[i], kappa * (1 - pool_pos[i]));
    }
    //Sensitivity
    for(i in 1:J_sens){
      target += binomial_lpmf(y_sens[i] | n_sens[i], sens);
    }
  }
}

generated quantities{

}
