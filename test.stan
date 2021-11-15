data {
 int < lower =0 > N; // number of data points
 vector[N ]  x ; // observation year
 vector[N ]  y; // observation life expectancies, Nc countries
 real xpred ; // prediction year
}

transformed data{
  vector[N ]  x1 = x - x[1]; // transform x to 0,1,2,3,4....
  vector[N ]  x2 ;
  vector[N ]  x3 ; 
  real xpred1 = xpred - x[1];
  for ( j in 1: N){
    x2[j] = x1[j]*x1[j];
    x3[j] = x2[j]*x1[j];
  }
}


parameters {
  real A;
  real B;
  real C;
  real D;
  real < lower =0 > sigma;
}
//transformed parameters {
  //vector[N ] mu = A*x3 + B*x2 + C*x1 + D;
//}
 
 model {
      A ~ normal(0,0.1);
      B ~ normal(0,1);
      C ~ normal(0,10);
      D ~ normal(0,100);
      sigma ~inv_chi_square(0.01);
      
      y ~ normal (A*x3 + B*x2 + C*x1 + D , sigma);
   
 }
 
 generated quantities {
  real ypred = normal_rng(A*xpred1^3 + B*xpred1^2 + C*xpred1 + D , sigma );
}
