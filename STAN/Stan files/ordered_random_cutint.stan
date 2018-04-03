data {
int<lower=2> K;// N Categories
int<lower=0> N; // N Observations
int<lower=1> D; // N Explanatory variables (source state)
int<lower=1> Fi; // Number of groups
int<lower=1,upper=K> y[N]; // Outcome
row_vector[D] x[N]; // Design Matrix
int<lower=1,upper=Fi> field[N]; // Indexing variable
}

transformed data { // Transforms for K-1 Parameterisation
int<lower=0> beta_zero;
beta_zero = 0;
}

parameters {
vector[D-1] beta_raw;
ordered[K-1]   c[Fi]; // Field level cutpoints
vector[Fi] gamma; // field effect in linear predictor
real   CutpointsMu; // means for cutpoint differences
real<lower=0> CutpointsSigma; // sds for cutpoint differences
real<lower=0> gammaSigma; // Field effect variance

}

transformed parameters{
vector[D] beta;
beta = append_row(beta_zero,beta_raw);
}

model {
// Prior specification
CutpointsSigma ~ cauchy(0,5);
CutpointsMu ~normal(0,10);
beta_raw ~ normal(0,10);


for(f in 1:Fi){
c[f][1] ~ normal(0,1); // Prior for first cutpoint
// hierarchical cutpoints on differences between cutpoints
(c[f][2] - c[f][1]) ~ normal(CutpointsMu, CutpointsSigma);
(c[f][3] - c[f][2]) ~ normal(CutpointsMu, CutpointsSigma);
(c[f][4] - c[f][3]) ~ normal(CutpointsMu, CutpointsSigma);}

// Sampling
for (n in 1:N){
y[n] ~ ordered_logistic((x[n] * beta) + gamma[field[n]], c[field[n]]);
}
}

generated quantities{ // For Loo CV
vector[N] log_lik;

for(n in 1:N){
log_lik[n] = ordered_logistic_lpmf(y[n] | (x[n] * beta) + gamma[field[n]], c[field[n]]);
}
}
