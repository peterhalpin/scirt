// stan code for running irt

data {

  int<lower=1> J_ind; // n items on individual form
  int<lower=1> K; // n dyads

  // individual data, respondent 1
  int<lower=1> N_ind1; // n non-missing responses on individual form, partner 1

  int<lower=1, upper=J_ind> jj_ind1[N_ind1]; // indexing non-missing data, partner 1

  int<lower=1, upper=K> kk_ind1[N_ind1]; // indexing non-missing data

  int<lower=0, upper=1> y_ind1[N_ind1]; // individual responses

  // item parms (calibrated)
  real<lower=0> alpha_ind[J_ind];  // discrimination parms for individual form
  real beta_ind[J_ind];            // difficulty parms for individual form
}

parameters {
  real theta1[K]; // indvidual abilities
}

model {

  // IRFS
  real p1 = 0;

  // log likelihoods
  real l1 = 0;

  // priors
  for (k in 1:K) {
    theta1[k] ~ normal(0, 1);
  }

  // log likelihoods for individual responses, partner 1
  for (n in 1:N_ind1) {
    p1 = inv_logit(alpha_ind[jj_ind1[n]] * (theta1[kk_ind1[n]] - beta_ind[jj_ind1[n]]));
    l1 = l1 + log(p1) * y_ind1[n] + log(1 - p1) * (1 - y_ind1[n]);
  }
  target += l1;
}
