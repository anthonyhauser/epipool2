data {
  //data to fit
  int N; //number of observations
  int N_t; //number of distinct time points
  array[N] int rank_t; //rank of the time point related to the observations
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
  array[2] real p_alpha;
  array[2] real p_lambda;
  real spec;

  int inference;
}

transformed data {
  real delta = 1e-9;
}

parameters {
  //True prevalence
  real intercept; //intercept of the prevalence
  vector[N_t] beta;
  real <lower=0> lambda; // lengthscale of f
  real<lower=0> alpha; //signal variance
  real <lower=0,upper=1> sens; //sensitivity
}

transformed parameters {
  vector[N_t] f; //logit-transformed prevalence
  vector[N_t] prev_f; //prevalence
  array[N] real pool_pos; //probability of pool being positive

  matrix[N_t, N_t] K; //kernel matrix
  matrix[N_t, N_t] L_K; //cholesky decomposition of K

  //Construct GP
  K = gp_exp_quad_cov(t, alpha, lambda); //define Kernel matrix
  for (l in 1:N_t) K[l, l] = K[l, l] + delta;
  L_K = cholesky_decompose(K); //Cholesky decomposition of Kernel
  f = L_K * beta; //GP
  prev_f = inv_logit(intercept + f);//back logit-tranformed the GP to prevalence

  //Calculate the probability of pool being positive for each observation
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
  sens ~ beta(p_sens[1], p_sens[2]);

  if(inference==1){
    //Pool positivity
    for(i in 1:N){
      target += binomial_lpmf(k[i] | n[i], pool_pos[i]);
    }
    //Sensitivity
    for(i in 1:J_sens){
      target += binomial_lpmf(y_sens[i] | n_sens[i], sens);
    }
  }
}

generated quantities{

}
