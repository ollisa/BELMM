data {
  int<lower=0> K; // # series
  int<lower=0> T; // # observations / series
  matrix[K, T] z; // observations
  real<lower=0> sigma_y[8]; // latent 
  real a;
  real b;
}
parameters {
  array[3] real<lower=0, upper=20> epsilon; // error to the latent process
  simplex[3] theta; //weights
  array[3] vector[T] l; // latent variables
  simplex[4] cht;
}
transformed parameters {
  vector[3] log_theta;
  ordered[3] l_init;
  
  log_theta = log(theta);
  l_init = a + head(cumulative_sum(cht), 3) * (b - a); //https://groups.google.com/g/stan-users/c/04GSu-ql3vM
}
model {
  vector[3] lp; // for iterating log probability
  
  target += dirichlet_lpdf(theta | [1,1,1]');
  
  for (k in 1 : K) {
    lp = log_theta;
    for (h in 1 : 3) 
      lp[h] += normal_lpdf(z[k,  : ] | l[h,  : ], epsilon[h]);
    target += log_sum_exp(lp);
   }
   
   for (h in 1 : 3){ 
   target += inv_gamma_lpdf(epsilon[h] | 1, 1);
   target += normal_lpdf(l_init[h] | 10, 6);
   target += normal_lpdf(l[h, 1] | l_init[h], sigma_y[h]);
   target += normal_lpdf(l[h, 2 : T] | l[h, 1 : (T - 1)], sigma_y[h]);
   }
  
}
generated quantities {
  real logModel;
  vector[3] log_theta2 = log_theta;
  
  for (k in 1 : K) {
    for (h in 1 : 3) 
      log_theta2[h] += normal_lpdf(z[ k , :] | l[h,  : ], epsilon[h]);
    
    logModel = log_sum_exp(log_theta2) + dirichlet_lpdf(theta | [1,1,1]');
  }
  for(h in 1:3)
    logModel +=  inv_gamma_lpdf(epsilon[h] | 1, 1)
                      + normal_lpdf(l_init[h] | 10 ,6)
                      + normal_lpdf(l[h, 1] | l_init[h], sigma_y[h])
                      + normal_lpdf(l[h, 2 : T] | l[h, 1 : (T - 1)], sigma_y[h]);
}


