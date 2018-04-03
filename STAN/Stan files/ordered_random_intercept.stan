data {
int<lower=2> K; // N of Categories
int<lower=0> N; // N Observations
int<lower=1> D; // N Explanatory variables
int<lower=1> Fi; // N groups
int<lower=1,upper=K> y[N]; // Outcome
row_vector[D] x[N]; // Design matrix
int<lower=1,upper=Fi> field[N]; // Indexing variable
}

transformed data { // Transforms for K-1 Parameterisation
int<lower=0> beta_zero;
beta_zero = 0;
}

parameters {
vector[D-1] beta_raw;
vector[Fi] gamma; // group-level intercept
real<lower=0.001> sigmaGamma;  // standard deviation for random intercept
ordered[K-1] c;
}

transformed parameters{
vector[D] beta;
beta = append_row(beta_zero,beta_raw);
}


model {

// Prior specification
sigmaGamma ~ cauchy(0,5);
gamma ~ normal(0,sigmaGamma);
beta_raw ~ normal(0,10);
// Sampling
for (n in 1:N){
y[n] ~ ordered_logistic((x[n] * beta) + gamma[field[n]], c);
}

}




generated quantities{ // For Loo CV 
vector[N] log_lik;
for(n in 1:N)
log_lik[n] = ordered_logistic_lpmf(y[n] |(x[n] * beta) + gamma[field[n]], c);
}
