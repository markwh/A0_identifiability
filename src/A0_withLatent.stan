
data {

  // Dimensions
  int<lower=0> nt; // number of observation times
  vector<lower = 0>[nt] dA;
  real logA0_hat;
  real<lower=0> sigmaA_hat;
  // vector<lower=0>[nt] et_est;
  // vector[nt] ds_est;
  // 
  // vector<lower = 0>[nt] p_sd;
  // vector<lower=0>[nt] et_sd;
  // vector<lower=0>[nt] ds_sd;
  
  
}


parameters {
  
  real<lower=0> A0;
  real<lower=0> sigma_A;
  // vector<lower = 0>[nt] p;
  // vector<lower=0>[nt] et;
  // vector[nt] ds;

}



transformed parameters {
  
  vector[nt] logA;
  
  logA = log(A0 + dA);

}



model {
  
  // target += dot_self()
  
  logA ~ normal(log(A0), sigma_A);
  
  target += logA;

}
