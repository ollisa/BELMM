data {
  int<lower=0> K; // # series
  int<lower=0> T; // # observations / series
  int<lower=1> H; // # clusters
  matrix[T, K] z; // observations
  real<lower=0> sigma_y;  
}
parameters {
  vector<lower=0,upper=20>[H] epsilon;
  simplex[H] theta; 
  vector[T] l[H]; 
}
model {
  
  vector[H] lp; 
  vector[H] log_theta;
  
  log_theta = log(theta);
  
  //Prior for latent process
  for (h in 1:H){
    l[h,1] ~ double_exponential(-30,20);
  }
  for (h in 1:H){
    for (t in 2:T){
      l[h,t] ~ double_exponential(l[h,(t-1)],sigma_y);
    }
    epsilon[h] ~ inv_gamma(1,1);
  }
  
  for (k in 1:K) {
    lp = log_theta;
    for (h in 1:H)
    
      lp[h] = lp[h]  + normal_lpdf(z[,k] | l[h,],epsilon[h]); 
    
    target += log_sum_exp(lp);
  }
}
generated quantities{ 
   real logModel; 
   vector[H] log_theta=log(theta); 
   
   for (k in 1:K) {
    for (h in 1:H)
      log_theta[h] += normal_lpdf(z[,k] | l[h,],epsilon[h]); 
   }
         
   logModel = log_sum_exp(log_theta);
   for (h in 1:H)
   logModel += double_exponential_lpdf(l[h,1]|-30,20) + double_exponential_lpdf(l[h,2:T] |l[h,1:(T-1)],sigma_y );
}
