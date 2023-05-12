data {
  int<lower=0> K; // # series
  int<lower=0> T; // # observations / series
  matrix[K, T] z; // observations
  real<lower=0> sigma_y;  
  real a;
  real b;
}
parameters {
  array[3] real<lower=0, upper=20> epsilon; 
  simplex[3] theta; 
  array[3] vector[T] l; 
  simplex[4] cht;
}
transformed parameters {
  vector[3] log_theta;
  ordered[3] l_init;
  
  log_theta = log(theta);
  l_init = a + head(cumulative_sum(cht), 3) * (b - a); //https://groups.google.com/g/stan-users/c/04GSu-ql3vM
}
model {
  vector[3] lp; 
    
  for (k in 1 : K) {
    lp = log_theta;
    for (h in 1 : 3) 
      lp[h] += normal_lpdf(z[k,  : ] | l[h,  : ], epsilon[h]);
    target += log_sum_exp(lp);
    for (h in 1 : 3) 
      target += inv_gamma_lpdf(epsilon[h] | 1, 1)
                + normal_lpdf(l_init[h] | -3, 6)
                + normal_lpdf(l[h, 1] | l_init[h], sigma_y)
                + normal_lpdf(l[h, 2 : T] | l[h, 1 : (T - 1)], sigma_y);
  }
}
generated quantities {
  real logModel;
  vector[3] log_theta2 = log_theta;
  
  for (k in 1 : K) {
    for (h in 1 : 3) 
      log_theta2[h] += normal_lpdf(z[ k , :] | l[h,  : ], epsilon[h]);
    
    logModel = log_sum_exp(log_theta2);
  }
  for(h in 1:3)
    logModel +=   inv_gamma_lpdf(epsilon[h] | 1, 1)
                      + normal_lpdf(l_init[h] | -3, 6)
                      + normal_lpdf(l[h, 1] | l_init[h], sigma_y)
                      + normal_lpdf(l[h, 2 : T] | l[h, 1 : (T - 1)], sigma_y);
}


