data {
  int<lower=0> K; // # series
  int<lower=0> T; // # observations / series
  matrix[K, T] z; // observations
  real<lower=0> sigma_y; // latent 
  real a;
  real b;
}
parameters {
  real<lower=0,upper=20> epsilon[8]; // error to the latent process
  simplex[8] theta; //weights
  vector[T] l[8]; // latent variables
  simplex[9] cht;
}
transformed parameters {
  vector[8] log_theta;
  ordered[8] l_init;
  
  log_theta = log(theta);
  l_init = a + head(cumulative_sum(cht),8)*(b-a); //https://groups.google.com/g/stan-users/c/04GSu-ql3vM
  
} 
model {
  
  vector[8] lp; // for iterating log probability
  
  //Prior for latent process
  l_init ~ normal(-6,3);
  
  for (h in 1:8){
    epsilon[h] ~ inv_gamma(1,1);
  }
  
  for (h in 1:8){
   l[h,1] ~ normal(l_init[h],sigma_y);
    for (t in 2:T){
      l[h,t] ~ normal(l[h,(t-1)],sigma_y);
    }
  }
  
  for (k in 1:K) {
    lp = log_theta;
    for (h in 1:8)
      lp[h] += normal_lpdf(z[k,] | l[h,],epsilon[h]); 
    target += log_sum_exp(lp);
    
  }
}
generated quantities{ 
   real logModel;
   vector[8] log_theta2=log_theta;
   
   for (k in 1:K) {
    for (h in 1:8)
      log_theta2[h] += normal_lpdf(z[k,] | l[h,],epsilon[h]); 
    
    logModel = log_sum_exp(log_theta2)+dirichlet_lpdf(theta | [1,1,1,1,1,1,1,1]');
   }
    for (h in 1:8)
    logModel += inv_gamma_lpdf(epsilon[h] | 1,1)+normal_lpdf(l_init[h] | -6,3) +normal_lpdf(l[h,1]|l_init[h],sigma_y) + normal_lpdf(l[h,2:T] |l[h,1:(T-1)],sigma_y );
   

}

