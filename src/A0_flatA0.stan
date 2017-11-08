
data {

  // Dimensions
  int<lower=0> nt; // number of observation times
  vector[nt] dA;
  real<lower=0> sigmalogA_hat;
  real<lower=0> sigmalogA_sd;
}

transformed data {
  real<lower = 0> minA0;
  
  minA0 = 1 - min(dA);
}

parameters {
  
  real<lower = minA0> A0;
  real<lower=0> sigma_logA;
}


transformed parameters {
  vector[nt] logA;
  logA = log(A0 + dA);
}



model {
  
  sigma_logA ~ normal(sigmalogA_hat, sigmalogA_sd);
  logA ~ normal(log(A0), sigma_logA);
  
  target += -logA;

}
