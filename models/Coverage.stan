data {
  int<lower=0> N;
  real<lower=0, upper=1> Coverage[N];
  real<lower=0> VacT[N];
}

transformed data {
  real y[N];
  
  for (i in 1:N) {
    y[i] = log(1 - Coverage[i]);
  }
}

parameters {
  real<lower=0, upper=1> p_initial;
  real<lower=0, upper=1> p_catch;
  real<lower=0> sigma;
}

transformed parameters {
  real mu[N];
  real beta0 = log(1 - p_initial);
  real beta1 = log(1 - p_catch);
  
  
  for (i in 1:N) {
    mu[i] = beta0 + VacT[i] * beta1;
  }
}

model {
  for (i in 1:N) {
    target += normal_lpdf(y[i] | mu[i], sigma);
  }
}

generated quantities {
  vector[10] Coverage_pred;
  vector[10] Coverage_fitted;
  
  for (i in 1:10) {
    Coverage_pred[i] = 1 - exp(normal_rng(beta0 + (i - 1) * beta1, sigma));
    Coverage_fitted[i] = 1 - exp(beta0 + (i - 1) * beta1);
  }
}

