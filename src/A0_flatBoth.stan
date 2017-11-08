
data {

  // Dimensions
  int<lower=0> nt; // number of observation times
  vector[nt] dA;
  // real logA0_hat;
  // real<lower=0> sigmaA_hat;
}

transformed data {
  real<lower = 0> minA0;
  
  minA0 = 1 - min(dA);
}

parameters {
  
  real<lower = minA0> A0;
  real<lower=0> sigma_A;
}



transformed parameters {
  vector[nt] logA;
  logA = log(A0 + dA);
}



model {
  
  // target += dot_self()
  
  logA ~ normal(log(A0), sigma_A);
  
  target += -logA;

}
