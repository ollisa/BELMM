data {
  int<lower=0> K; // # series
  int<lower=0> T; // # observations / series
  matrix[T, K] z; // observations 
  real<lower=0> sigma_y;
  real a;
  real b;
}
parameters {
  array[2] real<lower=0, upper=20> epsilon;
  simplex[2] theta;
  array[2] vector[T] l;
  simplex[3] cht;
}
transformed parameters {
  vector[2] log_theta;
  ordered[2] l_init;
  
  log_theta = log(theta);
  l_init = a + head(cumulative_sum(cht), 2) * (b - a); //https://groups.google.com/g/stan-users/c/04GSu-ql3vM
}
model {
  vector[2] lp;
  
  //Prior for latent process
  l_init ~ normal(3, 6);
  
  for (h in 1 : 2) {
    epsilon[h] ~ inv_gamma(1, 1);
  }
  
  for (h in 1 : 2) {
    l[h, 1] ~ normal(l_init[h], sigma_y);
    for (t in 2 : T) {
      l[h, t] ~ normal(l[h, t - 1], sigma_y);
    }
  }
  
  for (k in 1 : K) {
    lp = log_theta;
    for (h in 1 : 2) {
      lp[h] = lp[h] + normal_lpdf(z[ : , k] | l[h,  : ], epsilon[h]);
    }
    
    target += log_sum_exp(lp);
  }
}
generated quantities {
  real logModel;
  vector[2] lp2;
  vector[2] log_theta2 = log_theta;
  
  for (k in 1 : K) {
    for (h in 1 : 2) {
      log_theta2[h] = log_theta2[h]
                      + normal_lpdf(z[ : , k] | l[h,  : ], epsilon[h]);
      lp2[h] = inv_gamma_lpdf(epsilon[h] | 1, 1)
               + normal_lpdf(l_init[h] | 3, 6)
               + normal_lpdf(l[h, 1] | l_init[h], sigma_y)
               + normal_lpdf(l[h, 2 : T] | l[h, 1 : (T - 1)], sigma_y);
    }
    logModel = log_sum_exp(log_theta2) + sum(lp2);
  }
}


