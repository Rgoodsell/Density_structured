

data {
int<lower=2> K;
int<lower=0> N;
int<lower=1> D;
int<lower=1> Fi;
int<lower=1,upper=K> y[N];
row_vector[D] x[N];
int<lower=1,upper=Fi> field[N];
}

transformed data {
int<lower=0> beta_zero;
row_vector[Fi] gamma_zeroes;

gamma_zeroes = rep_row_vector(0,Fi);
beta_zero = 0;
}

parameters {
vector[D-1] beta_raw;
matrix[D-1,Fi] gamma_raw;
vector<lower=0>[D-1] sigmaGamma;  // sigma for intercept
ordered[K-1] c;


}

transformed parameters{
vector[D] beta;
matrix[D,Fi] gamma;
matrix[D,Fi] eta;
beta = append_row(beta_zero,beta_raw);
gamma = append_row(gamma_zeroes,gamma_raw);

for(d in 1:D){
for(f in 1:Fi){
eta[d,f] = beta[d] + gamma[d, f];
}}
}
model {

sigmaGamma ~ cauchy(0,5);



for(d in 1:D-1){
beta_raw[d] ~ normal(0,5);
for(f in 1:Fi){
gamma_raw[d,f] ~ normal(0,sigmaGamma);
}}


for (n in 1:N){
y[n] ~ ordered_logistic(x[n] * eta[,field[n]], c);
}
}
generated quantities{
vector[N] log_lik;

for(n in 1:N)
log_lik[n] = ordered_logistic_lpmf(y[n]| x[n]*eta[,field[n]],c);


}
