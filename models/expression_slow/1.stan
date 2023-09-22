data {
  int<lower=0> K; // # series
  int<lower=0> T; // # observations / series
  matrix[K, T] z; // observations 
  real<lower=0> sigma_y[8]; // latent 
  real a;
  real b;
}
parameters {
  real<lower=0,upper=20> epsilon; // deviance from the latent process
  vector[T] l; // latent variables
  real l_init;
}
transformed parameters {
} 
model {
  
  //Prior for latent process
  

target += inv_gamma_lpdf(epsilon | 1, 1);
target += normal_lpdf(l_init | 10, 6);
target += normal_lpdf(l[1] | l_init, sigma_y[1]);
target += normal_lpdf(l[2 : T] | l[1 : (T - 1)], sigma_y[1]);
   
  for (k in 1:K) {
      target += normal_lpdf(z[k,] | l,epsilon); 
  }
}
generated quantities{ 
  real logModel;
   
  logModel = inv_gamma_lpdf(epsilon | 1,1)+normal_lpdf(l_init | 10,3) +normal_lpdf(l[1]|l_init,sigma_y[1]) + normal_lpdf(l[2:T] |l[1:(T-1)],sigma_y[1] ); 
   
   for (k in 1:K) {
      logModel +=  normal_lpdf(z[k,] | l,epsilon); 
   }

    
}
