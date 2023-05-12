data {
  int<lower=0> K; //  series
  int<lower=0> T; // observations / series
  matrix[T, K] z; // observations
  real<lower=0> sigma_y;  
}
parameters {
  vector<lower=0,upper=20>[4] epsilon; 
  simplex[4] theta; 
  vector[T] l[4]; 
  simplex[5] cht;
}
transformed parameters{
  ordered[4] l_init;
  l_init = -320 + head(cumulative_sum(cht),4)*(450); //https://groups.google.com/g/stan-users/c/04GSu-ql3vM
  
}
model {
  
  vector[4] lp; 
  vector[4] log_theta;
  
  log_theta = log(theta);
  
  //Prior for latent process
  for (h in 1:4){
    l_init ~ normal(-30,100);
  }
  for (h in 1:4){
    l[h,1] ~ normal(l_init[h],sigma_y);
    for (t in 2:T){
      l[h,t] ~ normal(l[h,(t-1)],sigma_y);
    }
    epsilon[h] ~ inv_gamma(1,1);
  }
  
  for (k in 1:K) {
    lp = log_theta;
    for (h in 1:4)
      lp[h] = lp[h] + normal_lpdf(z[,k] | l[h,],epsilon[h]); 
    
    target += log_sum_exp(lp);
  }
}
generated quantities{ 
   real logModel; 
   vector[4] log_theta=log(theta); 
   
   for (k in 1:K) {
    for (h in 1:4)
      log_theta[h] += normal_lpdf(z[,k] | l[h,],epsilon[h]); 
   }
         
   logModel = log_sum_exp(log_theta) +normal_lpdf(l[1,1]|l_init[1],20) + normal_lpdf(l[1,2:T] |l[1,1:(T-1)],sigma_y ) +normal_lpdf(l[2,1]|l_init[2],20) + normal_lpdf(l[2,2:T] |l[2,1:(T-1)],sigma_y )+normal_lpdf(l[3,1]|l_init[3],20) + normal_lpdf(l[3,2:T] |l[3,1:(T-1)],sigma_y )+normal_lpdf(l[4,1]|l_init[4],20) + normal_lpdf(l[4,2:T] |l[4,1:(T-1)],sigma_y );
}
