data{
  int<lower=0> N;
  int<lower=2> k1;
  int<lower=2> k2;
  vector[N] td;
  vector[N] x1;
  vector[N] x2;
  vector[N] x3;
  vector[N] x4;
  matrix[N,k1-1] c1;
  matrix[N,k2-1] c2;
  vector[N] y;
  // for priors
  real ar_alpha;
  real ar_beta;
}
parameters{
  real beta0;
  real beta1;
  real beta2;
  real beta3;
  real beta4;
  real<lower=0> sigma;
  real<lower=0> sigma_t;
  real error_0;
  vector[k1-1] beta_c1;
  vector[k2-1] beta_c2;
  real<lower=0,upper=1> beta_t;
  vector[N] drift_errors;
}
model {
  vector[N] error_t;
  error_t[1] = error_0;
  for (n in 2:N){
    error_t[n] = error_t[n-1] * beta_t^(td[n]) + drift_errors[n]*
    sqrt((1-beta_t^(2*td[n]))/(1-beta_t^2));
  }
  // estimate
  y ~ normal(
    error_t +
    beta0 + beta1 * x1 + beta2*x2 + beta3*x3 + beta4 * x4 +
     (c1 * beta_c1 + c2 * beta_c2),
    sigma
  );
  // priors to keep values reasonable
  beta0 ~ normal(0,10);
  beta1 ~ normal(0,10);
  beta2 ~ normal(0,10);
  beta3 ~ normal(0,10);
  beta4 ~ normal(0,10);
  beta_c1 ~ normal(0,10);
  beta_c2 ~ normal(0,10);
  beta_t ~ beta(ar_alpha,ar_beta);
  sigma ~ lognormal(0,0.5);
  sigma_t ~ lognormal(0,0.5);
  drift_errors ~ normal(0,sigma_t);
  error_0 ~ normal(0,sigma_t/(1-beta_t));
}

generated quantities {
  vector[N] error_t_hat;
  vector[N] y_hat;
  vector[N] y0_hat;
  vector[N] c1_contribution;
  vector[N] noncat_error;
  // estimated error at each step (same as code above)
  error_t_hat[1] = error_0;
  for (n in 2:N){
    error_t_hat[n] = error_t_hat[n-1] * beta_t^(td[n]) + drift_errors[n]*
    sqrt((1-beta_t^(2*td[n]))/(1-beta_t^2));
  }
  // calculated estimate
  y_hat = error_t_hat + beta0 + beta1 * x1 + beta2*x2 + beta3*x3 + beta4 * x4 +
    (c1 * beta_c1 + c2 * beta_c2);
  // calcualted estimate w/o time error
  y0_hat = y_hat - error_t_hat;
  // calculated estimate without categorical effects
  noncat_error = y - y_hat - (c1 * beta_c1 + c2 * beta_c2);
  // estimated categorical effect at each point
  c1_contribution = c1 * beta_c1;

}
