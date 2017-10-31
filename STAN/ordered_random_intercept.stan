data {
int<lower=2> K;
int<lower=0> N;
int<lower=1> D;
int<lower=1> Fi; // Number of groups (fields or farms)
int<lower=1,upper=K> y[N];
row_vector[D] x[N];
int<lower=1,upper=Fi> field[N]; // Indexing variable
}

transformed data {
int<lower=0> beta_zero;
beta_zero = 0;
}

parameters {
vector[D-1] beta_raw;
vector[Fi] gamma; // Random intercept
real sigmaGamma;  // Varianace for random intercept
ordered[K-1] c;


}

transformed parameters{
vector[D] beta;
beta = append_row(beta_zero,beta_raw);

}
model {


sigmaGamma ~ cauchy(0,20);

for(d in 1:D-1){
beta_raw[d] ~ normal(0,5);
}

for(f in 1:Fi){
gamma[f] ~ normal(0,sigmaGamma);
}


for (n in 1:N){
y[n] ~ ordered_logistic((x[n] * beta) + gamma[field[n]], c);
}

}




generated quantities{
vector[N] log_lik;



for(n in 1:N)
log_lik[n] = ordered_logistic_lpmf(y[n] |(x[n] * beta) + gamma[field[n]], c);




}
