
data {
int<lower=2> K; // N Categories
int<lower=0> N; // N observations
int<lower=1> D; // N explanatory variables (source state)
int<lower=1> Fi; // N fields/group level effects
int<lower=1,upper=K> y[N]; // Outcome
row_vector[D] x[N]; // Design matrix
int<lower=1,upper=Fi> field[N]; // Group level indexing variable
}

transformed data { // Transforms for K-1 parameterisation
int<lower=0> beta_zero;
vector[D] gammaMu;
beta_zero = 0;
gammaMu = rep_vector(0,D); // Mu for MVN distribution over gamma
}

parameters {
vector[D-1] beta_raw; // Source state effect
matrix[D,Fi] gamma; // Field effect of source state
cholesky_factor_corr[D] L_corr;  // Sigma for intercept
vector<lower=0>[D] sigma_L ; // Standard deviation of correlation matrix components
ordered[K-1] c; // Cutpoints
}

transformed parameters{
vector[D] beta;
matrix[D,Fi] eta;
beta = append_row(beta_zero,beta_raw);

// Construct field level source state effect
for(d in 1:D){
for(f in 1:Fi){
eta[d,f] = beta[d] + gamma[d, f];
}}
}

model {

// Prior specifications
beta_raw ~ normal(0,10);
sigma_L ~ cauchy(0,5);
L_corr ~ lkj_corr_cholesky(1);

for(f in 1:Fi){
gamma[,f] ~ multi_normal_cholesky(gammaMu,diag_pre_multiply(sigma_L,L_corr));
}

// Sampling
for (n in 1:N){
y[n] ~ ordered_logistic(x[n] * eta[,field[n]], c);
}
}


generated quantities{ // For Loo CV
vector[N] log_lik;

for(n in 1:N)
log_lik[n] = ordered_logistic_lpmf(y[n]| x[n]*eta[,field[n]],c);


}
