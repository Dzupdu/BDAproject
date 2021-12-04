data {
 int < lower =0 > N; // number of data points
 int < lower =0 > Nc; //number of countries/continents
 matrix [Nc,N] y; // observation life expectancies, Nc countries
 int< lower =0 > Npred ; // prediction year
}

transformed data{
  vector[N] x1;
  vector[N] x2;
  vector[N] x3;// just use years 1,2,...N instead giving x as argument
  for ( i in 1: N){
    x1[i] = i;
    x2[i] = i*i;
    x3[i] = i*i*i;
  }
  
}

parameters {
  real A[Nc];
  real B[Nc];
  real C[Nc];
  real D[Nc];
  real < lower =0 > sigma;
  real Amu0;
  real < lower =0 > Asigma0;
  real Bmu0;
  real < lower =0 > Bsigma0;
  real Cmu0;
  real < lower =0 > Csigma0;
  real Dmu0;
  real < lower =0 > Dsigma0;
}
 
 model {
    Amu0 ~ normal(0,0.1);
    Bmu0 ~ normal(0,1);
    Cmu0 ~ normal(0,10);
    Dmu0 ~ normal(0,100);
    Asigma0 ~ inv_chi_square(10);
    Bsigma0 ~ inv_chi_square(1);
    Csigma0 ~ inv_chi_square(0.1);
    Dsigma0 ~ inv_chi_square(0.01);
    sigma ~ inv_chi_square(0.001);
   
    A ~ normal(Amu0,Asigma0);
    B ~ normal(Bmu0,Bsigma0);
    C ~ normal(Cmu0,Csigma0);
    D ~ normal(Dmu0,Dsigma0);
   
    for ( j in 1: Nc ){
      y[j, ] ~ normal (A[j]*x3 + B[j]*x2 + C[j]*x1 + D[j] , sigma);
    }
 }
 
 generated quantities {
    matrix [Nc,Npred] ypred;
    vector [Nc*N] log_like;
    vector [Nc*N] y_rep;
    
    for ( j in 1: Nc ){
      for(i in (N+1):(N+Npred)){
        ypred[j,i - N] = normal_rng(A[j]*i*i*i + B[j]*i*i + C[j]*i + D[j] , sigma);
      }
    }
    
    for ( j in 1: Nc ){
      for ( i in 1: N){
        log_like[(j-1)*N + i] = normal_lpdf(y[j, i] | (A[j]*i*i*i + B[j]*i*i + C[j]*i + D[j]) , sigma);
        y_rep[(j-1)*N + i] = normal_rng(A[j]*i*i*i + B[j]*i*i + C[j]*i + D[j] , sigma);
      }
    }
}