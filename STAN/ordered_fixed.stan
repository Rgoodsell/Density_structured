
data {
int<lower=2> K; // Categories
int<lower=0> N; // No. Observations
int<lower=1> D; // Explanatory Variables (source state)
int<lower=1,upper=K> y[N]; // Outcome
row_vector[D] x[N]; // Explanatory variable design matrix
}

transformed data {
int<lower=0> beta_zero;
beta_zero = 0;
}

parameters {
vector[D-1] beta_raw; // Source state effect
ordered[K-1] c; // Cutpoints


}

transformed parameters{ // K-1 parameterisation
vector[D] beta;
beta = append_row(beta_zero,beta_raw);
}
model {
for (n in 1:N)
y[n] ~ ordered_logistic(x[n] * beta, c);
}



generated quantities{
vector[N] log_lik;

for(n in 1:N)
log_lik[n] = ordered_logistic_lpmf(y[n]| x[n]*beta,c);

 
}