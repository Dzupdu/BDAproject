data {
 int < lower =0 > N; // number of data points
 int < lower =0 > Nc; //number of countries/continents
 matrix [Nc,N] y; // observation life expectancies, Nc countries
 int< lower =0 > Npred ; // predict this many years forward
}

transformed data{
  vector[N] x1; // just use years 1,2,...N instead giving x as argument
  for ( i in 1: N){
    x1[i] = i;
  }
  
}

parameters {
  real alpha[Nc];
  real beta[Nc];
  real<lower=.5,upper= 1> lambda[Nc];
  
  real < lower =0 > sigma;
  real bmu0;
  real < lower =0 > bsigma0;
  real amu0;
  real < lower =0 > asigma0;
  real < lower =0 > beta1;
  real < lower =0 > beta2;
}

 
 model {
    bmu0 ~ normal(0,100); // yearly change/slope
    amu0 ~ normal(0,100); // we have observations at x = 1 so lets just use something close to that
    bsigma0 ~ inv_chi_square(0.0001); // deviation around 100
    asigma0 ~ inv_chi_square(0.0001); // deviation around 100
    sigma ~ inv_chi_square(0.01); // deviation around 10 years
    beta1 ~ exponential(0.1); // mean 10
    beta2 ~ exponential(0.1); // mean 10
    lambda ~ beta(beta1, beta2); 
    
   beta ~ normal(bmu0,bsigma0);
   alpha ~ normal(amu0,asigma0);
   
   for ( j in 1: Nc ){
     for ( i in 1: N ){
        y[j, i] ~ normal(alpha[j] - beta[j] * pow(lambda[j], x1[i]) , sigma);
     }
   }
 }
 
 generated quantities {
    matrix [Nc,Npred] ypred;
    vector [Nc*N] log_like;
    vector [Nc*N] y_rep;
    
    for ( j in 1: Nc ){
      for(i in (N+1):(N+Npred)){
        ypred[j,i - N] = normal_rng(alpha[j] - beta[j] * pow(lambda[j], i) , sigma);
      }
    }
    
    for ( j in 1: Nc ){
      for ( i in 1: N){
        log_like[(j-1)*N + i] = normal_lpdf(y[j, i] | (alpha[j] - beta[j] * pow(lambda[j], i)) , sigma);
        y_rep[(j-1)*N + i] = normal_rng(alpha[j] - beta[j] * pow(lambda[j], i), sigma);
      }
    }
}
