data {
 int < lower =0 > N; // number of data points
 vector[N ]  x ; // observation year
 vector[N ] x2 ; // observation year
 vector[N ]  x3 ; // observation year
 vector[N ]  y; // observation life expectancies, Nc countries
 real xpred ; // prediction year
}

parameters {
  real A;
  real B;
  real C;
  real D;
  real < lower =0 > sigma;
}
transformed parameters {
  vector[N ] mu = A*x3 + B*x2 + C*x + D;
}
 
 model {
      //A ~ normal(0,100);
      //B ~ normal(0,100);
      //C ~ normal(0,100);
      //D ~ normal(0,100);
      //sigma ~inv_chi_square(0.01);
      
      y ~ normal (mu , sigma);
   
 }
 
 generated quantities {
  real ypred = normal_rng(A*xpred^3 + B*xpred^2 + C*xpred + D , sigma );
}
