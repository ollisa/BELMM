data {
  int<lower=0> K; // # series
  int<lower=0> T; // # observations / series
  matrix[K, T] z; // observations 
  array[8] real<lower=0> sigma_y;
  real a;
  real b;
}
parameters {
  real<lower=0, upper=20> epsilon;
  vector[T] l;
  real l_init;
}
transformed parameters {
  
}
model {
  //Prior for latent process
  l_init ~ normal(3, 6);
  
  epsilon ~ inv_gamma(1, 1);
  
  l[1] ~ normal(l_init, sigma_y[1]);
  for (t in 2 : T) {
    l[t] ~ normal(l[t - 1], sigma_y[1]);
  }
  
  for (k in 1 : K) {
    target += normal_lpdf(z[k,  : ] | l, epsilon);
  }
}
generated quantities {
  real logModel;
  
  logModel = inv_gamma_lpdf(epsilon | 1, 1) + normal_lpdf(l_init | 3, 6)
             + normal_lpdf(l[1] | l_init, sigma_y[1])
             + normal_lpdf(l[2 : T] | l[1 : (T - 1)], sigma_y[1]);
  
  for (k in 1 : K) {
    logModel += normal_lpdf(z[k,  : ] | l, epsilon);
  }
}


