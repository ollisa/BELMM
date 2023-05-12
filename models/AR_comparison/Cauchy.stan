data {
  int<lower=0> K; // # series
  int<lower=0> T; // # observations / series
  int<lower=1> H; // # clusters
  matrix[T, K] z; // observations
}
parameters {
  ordered[H] phi0; // AR(1) coefficients
  vector[H] phi1;
  vector<lower=0,upper=40>[H] sigma; 
  simplex[H] theta; 
}
model {
  vector[H] lp;
  vector[H] log_theta;
  
  for (h in 1:H){
  phi0[h] ~ normal(0,40);
  phi1[h] ~ normal(0,1);
  sigma[h] ~ inv_gamma(1,1) T[0,40];
  }
  log_theta = log(theta);
  for (k in 1:K) {
    lp = log_theta;
    for (h in 1:H)
      lp[h] = lp[h] + cauchy_lpdf(z[2:T,k] |phi0[h] + phi1[h] * z[1:(T-1),k],sigma[h]); 
    
    target += log_sum_exp(lp);
  }
}
generated quantities{
real logModel;
vector[H] log_theta =log(theta);

for (k in 1:K) {
 for (h in 1:H)
      log_theta[h] +=  cauchy_lpdf(z[2:T,k] |phi0[h] + phi1[h] * z[1:(T-1),k],sigma[h]); 
}
    logModel = log_sum_exp(log_theta);
 for (h in 1:H)
    logModel += normal_lpdf(phi0[h]|0,40) + normal_lpdf(phi1[h]|0,1) +inv_gamma_lpdf(sigma[h]|1,1); 
}

