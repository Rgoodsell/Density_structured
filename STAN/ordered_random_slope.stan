data {
int<lower=2> K; // N categories
int<lower=0> N; // N obvs
int<lower=1> D; // N Pop level effects
int<lower=1> Fi; // N fields
int<lower=1,upper=K> y[N]; // Obvs
row_vector[D] x[N]; // Explanatory vars
int<lower=1,upper=Fi> field[N]; // Field index
}

transformed data {
int<lower=0> beta_zero;
row_vector[Fi] gamma_zeroes;

gamma_zeroes = rep_row_vector(0,Fi);
beta_zero = 0;
}

parameters {
vector[D-1] beta_raw;
matrix[D-1,Fi] gamma_raw; // Slope
vector[Fi] alpha; // Intercept
vector<lower=0>[Fi] sigmaAlpha; // sigma for alpha
vector<lower=0>[D-1] sigmaGamma;  // sigma for slope
ordered[K-1] c;


}

transformed parameters{
vector[D] beta;
matrix[D,Fi] gamma;
matrix[D,Fi] eta;
beta = append_row(beta_zero,beta_raw);
gamma = append_row(gamma_zeroes,gamma_raw);

}
model {


sigmaAlpha ~ cauchy(0,20);

// estimate variance for alpha
for(f in 1:Fi){
alpha[f] ~ normal(0,sigmaAlpha);
}


for(d in 1:D-1){
sigmaGamma[d] ~ cauchy(0,20);
beta_raw[d] ~ normal(0,5);
for(f in 1:Fi){
gamma_raw[d,f] ~ normal(0,sigmaGamma[d]);
}}




for (n in 1:N){
y[n] ~ ordered_logistic((x[n]*beta) + (x[n]*gamma[,field[n]]) + alpha[field[n]], c);
}
}