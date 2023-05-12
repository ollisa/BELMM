data {
  int<lower=0> K; //  series
  int<lower=0> T; // observations / series
  matrix[T, K] z; // observations
  real<lower=0> sigma_y;  
}
parameters {
  vector<lower=0,upper=20>[2] epsilon; 
  simplex[2] theta; 
  vector[T] l[2]; 
  simplex[3] cht;
}
transformed parameters{
  ordered[2] l_init;
  l_init = -320 + head(cumulative_sum(cht),2)*(450); //https://groups.google.com/g/stan-users/c/04GSu-ql3vM
  
}
model {
  
  vector[2] lp; 
  vector[2] log_theta;
  
  log_theta = log(theta);
  
  //Prior for latent process
  for (h in 1:2){
    l_init ~ normal(-30,100);
  }
  for (h in 1:2){
    l[h,1] ~ normal(l_init[h],sigma_y);
    for (t in 2:T){
      l[h,t] ~ normal(l[h,(t-1)],sigma_y);
    }
    epsilon[h] ~ inv_gamma(1,1);
  }
  
  for (k in 1:K) {
    lp = log_theta;
    for (h in 1:2)
      lp[h] = lp[h] + normal_lpdf(z[,k] | l[h,],epsilon[h]); 
    
    target += log_sum_exp(lp);
  }
}
generated quantities{ 
   real logModel; 
   vector[2] log_theta=log(theta); 
   
   for (k in 1:K) {
    for (h in 1:2)
      log_theta[h] += normal_lpdf(z[,k] | l[h,],epsilon[h]); 
   }
         
   logModel = log_sum_exp(log_theta) +normal_lpdf(l[1,1]|l_init[1],20) + normal_lpdf(l[1,2:T] |l[1,1:(T-1)],sigma_y ) +normal_lpdf(l[2,1]|l_init[2],20) + normal_lpdf(l[2,2:T] |l[2,1:(T-1)],sigma_y );
}
