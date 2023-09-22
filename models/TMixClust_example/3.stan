data {
  int<lower=0> K; //  series
  int<lower=0> T; // observations / series
  matrix[T, K] z; // observations
  real<lower=0> sigma_y;
}
parameters {
  vector<lower=0, upper=20>[3] epsilon;
  simplex[3] theta;
  array[3] vector[T] l;
  simplex[4] cht;
}
transformed parameters {
  ordered[3] l_init;
  l_init = -320 + head(cumulative_sum(cht), 3) * 450; //https://groups.google.com/g/stan-users/c/04GSu-ql3vM
}
model {
  vector[3] lp;
  vector[3] log_theta;
  
  log_theta = log(theta);
  
  //Prior for latent process
  for (h in 1 : 3) {
    l_init ~ normal(-30, 100);
  }
  for (h in 1 : 3) {
    l[h, 1] ~ normal(l_init[h], sigma_y);
    for (t in 2 : T) {
      l[h, t] ~ normal(l[h, t - 1], sigma_y);
    }
    epsilon[h] ~ inv_gamma(1, 1);
  }
  
  for (k in 1 : K) {
    lp = log_theta;
    for (h in 1 : 3) {
      lp[h] = lp[h] + normal_lpdf(z[ : , k] | l[h,  : ], epsilon[h]);
    }
    
    target += log_sum_exp(lp);
  }
}
generated quantities {
  real logModel;
  vector[3] log_theta = log(theta);
  
  for (k in 1 : K) {
    for (h in 1 : 3) {
      log_theta[h] += normal_lpdf(z[ : , k] | l[h,  : ], epsilon[h]);
    }
  }
  
  logModel = log_sum_exp(log_theta) + normal_lpdf(l[1, 1] | l_init[1], 20)
             + normal_lpdf(l[1, 2 : T] | l[1, 1 : (T - 1)], sigma_y)
             + normal_lpdf(l[2, 1] | l_init[2], 20)
             + normal_lpdf(l[2, 2 : T] | l[2, 1 : (T - 1)], sigma_y)
             + normal_lpdf(l[3, 1] | l_init[3], 20)
             + normal_lpdf(l[3, 2 : T] | l[3, 1 : (T - 1)], sigma_y);
}


