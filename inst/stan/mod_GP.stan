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
  array[2] real p_alpha;
  array[2] real p_lambda;

  int inference;
}

transformed data {
  real delta = 1e-9;
}

parameters {
  //True prevalence
  array[N_pop] real intercept;
  array[N_pop] vector[N_t] beta;
  array[N_pop] real <lower=0> lambda; // lengthscale of f
  array[N_pop] real<lower=0> alpha;

  real <lower=0,upper=1> sens; //sensitivity
}



transformed parameters {
  array[N_pop] vector[N_t] f;
  array[N_pop] vector[N_t] prev_f;
  array[N] real pool_pos;

  array[N_pop] matrix[N_t, N_t] L_K;
  array[N_pop] matrix[N_t, N_t] K;

  for(i in 1:N_pop){
    K[i] = gp_exp_quad_cov(t, alpha[i], lambda[i]);
    for (l in 1:N_t) K[i, l, l] = K[i, l, l] + delta;
    L_K[i] = cholesky_decompose(K[i]);
    f[i] = L_K[i] * beta[i];
    prev_f[i] = inv_logit(intercept[i] + f[i]);
  }

  for(l in 1:N){
    pool_pos[l] = 1.0 - spec * pow(1.0 - prev_f[pop[l],rank_t[l]] * sens, s[l]);
  }
}

model {
  // priors
  for(i in 1:N_pop){
    beta[i] ~ normal(0, 1);
    intercept[i] ~ normal(p_intercept[1], p_intercept[2]);
    lambda[i] ~ normal(p_lambda[1], p_lambda[2]);
    alpha[i] ~ normal(p_alpha[1], p_alpha[2]);
  }

  sens ~ beta(p_sens[1], p_sens[2]);

  if(inference==1){
    for(i in 1:N){
      target += binomial_lpmf(k[i] | n[i], pool_pos[i]);
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
      prev_ratio[i,j] = prev_f[i,j+1]/prev_f[i,j];
    }
  }
}
