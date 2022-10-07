data {
  //data to fit
  int N; //number of observations
  int N_t; //number of distinct time points
  array[N] int rank_t; //rank of the time point related to the observations
  array[N_t] real t; //time points
  array[N] real s; //number of samples by pool
  array[N] int n; //number of pools
  array[N] int k; //number of positive pools


  //hyperparameters of the priors
  array[2] real p_intercept;
  real spec;
  array[2] real p_alpha;
  array[2] real p_lambda;
  real p_phi;

  int inference;
}

transformed data {
  real delta = 1e-9;
  real sens=1;
}

parameters {
  //True prevalence
  real intercept;
  vector[N_t] beta;
  real <lower=0> lambda;// lengthscale of f
  real<lower=0> alpha;

  real <lower=0>phi;
}

transformed parameters {
  real kappa = phi+2;

  vector[N_t] f;
  vector[N_t] prev_f;
  array[N] real pool_pos;

  matrix[N_t, N_t] L_K;
  matrix[N_t, N_t] K;


  K = gp_exp_quad_cov(t, alpha, lambda);
  for (l in 1:N_t) K[l, l] = K[l, l] + delta;
  L_K = cholesky_decompose(K);
  f = L_K * beta;
  prev_f = inv_logit(intercept + f);


  for(l in 1:N){
    pool_pos[l] = 1.0 - pow(1.0 - (prev_f[rank_t[l]] * sens + (1.0-prev_f[rank_t[l]]) * (1.0-spec)),s[l]);
  }
}

model {
  // priors
  beta ~ normal(0, 1);
  intercept ~ normal(p_intercept[1], p_intercept[2]);
  lambda ~ normal(p_lambda[1], p_lambda[2]);
  alpha ~ normal(p_alpha[1], p_alpha[2]);
  phi ~ exponential(p_phi);

  if(inference==1){
    for(i in 1:N){
      target += beta_binomial_lpmf(k[i] | n[i], kappa * pool_pos[i], kappa * (1 - pool_pos[i]));
    }
  }
}

generated quantities{

}
