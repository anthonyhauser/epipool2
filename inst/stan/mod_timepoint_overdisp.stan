data {
  //data to fit
  int N; //number of observations
  int N_t; //number of distinct time points
  int N_pop; //number of population
  array[N] int rank_t; //rank of the time point related to the observations
  array[N] int pop; //rank of the population
  array[N_t] real t; //time points
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
  real spec;

  //Overdispersion
  int N_phi;
  array[N_pop] int use_phi;
  real p_phi;

  int inference;
}

transformed data {
  real delta = 1e-9;
}

parameters {
  array[N_pop] vector[N_t] logit_prev_f; //logit-transformed prevalence (done so that we have same priors as for GP model)
  real <lower=0,upper=1> sens; //sensitivity
  vector<lower=0>[N_phi] phi; //overdispersion parameter
}

transformed parameters {
  vector<lower=0>[N_phi] kappa=phi+2.0; //transformed overdisperion parameter
  array[N_pop] vector[N_t] prev_f;
  array[N] real pool_pos;

  for(i in 1:N_pop){
    prev_f[i] = inv_logit(logit_prev_f[i]);
  }

  for(l in 1:N){
    pool_pos[l] = 1.0 - spec * pow(1.0 - prev_f[pop[l],rank_t[l]] * sens, s[l]);
  }
}

model {
  // priors
  for(i in 1:N_pop){
    logit_prev_f[i] ~ normal(p_intercept[1], p_intercept[2]);
  }
  phi ~ exponential(p_phi);
  sens ~ beta(p_sens[1], p_sens[2]);

  if(inference==1){
    for(i in 1:N){
      if(use_phi[pop[i]]>0){
        target += beta_binomial_lpmf(k[i] | n[i], kappa[use_phi[pop[i]]] * pool_pos[i] + 1e-6,
                                                  kappa[use_phi[pop[i]]] * (1 - pool_pos[i]) + 1e-6);
      }else{
        target += binomial_lpmf(k[i] | n[i], pool_pos[i]);
      }
    }
    for(i in 1:J_sens){
      target += binomial_lpmf(y_sens[i] | n_sens[i], sens);
    }
  }
}

generated quantities{
  array[N_pop, N_t-1] real prev_ratio;
  for(i in 1:N_pop){
    for(j in 1:(N_t-1)){
      prev_ratio[i,j] =prev_f[i,j+1]/prev_f[i,j];
    }
  }
}
