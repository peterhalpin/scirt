// stan code for running rsc

data {

  int<lower=1> J_ind; // n items on individual form
  int<lower=1> J_col; // n items on group form
  int<lower=1> K; // n dyads

  // individual data, respondent 1
  int<lower=1> N_ind1; // n non-missing responses on individual form, partner 1
  int<lower=1, upper=J_ind> jj_ind1[N_ind1]; // indexing non-missing data, partner 1
  int<lower=1, upper=K> kk_ind1[N_ind1]; // indexing non-missing data
  int<lower=0, upper=1> y_ind1[N_ind1]; // individual responses

  // individual data, respondent 2
  int<lower=1> N_ind2; // n non-missing responses on individual form, partner 2
  int<lower=1, upper=J_ind> jj_ind2[N_ind2]; // indexing non-missing data, partner 2
  int<lower=1, upper=K> kk_ind2[N_ind2]; // indexing non-missing data
  int<lower=0, upper=1> y_ind2[N_ind2]; // individual responses

  // group data
  int<lower=1> N_col; // n non-missing responses on group form
  int<lower=1, upper=J_col> jj_col[N_col]; // indexing non-missing data
  int<lower=1, upper=K> kk_col[N_col]; // indexing non-missing data
  int<lower=0, upper=1> y_col[N_col]; // group responses

  // item parms (calibrated)
  real<lower=0> alpha_ind[J_ind]; // discrimination parms for individual form
  real beta_ind[J_ind]; // difficulty parms for individual form
  real<lower=0> alpha_col[J_col]; // discrimination parms for group form
  real beta_col[J_col]; // difficulty parms for group form
}

parameters {
  real theta1[K]; // indvidual abilities
  real theta2[K]; // indvidual abilities
  real<lower=0, upper=1> w[K]; // beta parameterization weight
}

model {

  // IRFS
  real p1 = 0;
  real p2 = 0;
  real R = 0;

  // log likelihoods
  real l1 = 0;
  real l2 = 0;
  real lR = 0;

  // priors
  for (k in 1:K) {
    theta1[k] ~ normal(0, 1);
    theta2[k] ~ normal(0, 1);
    w[k] ~ beta(1.05, 1.05);
  }

  // log likelihoods for individual responses, partner 1
  for (n in 1:N_ind1) {
    p1 = inv_logit(alpha_ind[jj_ind1[n]] * (theta1[kk_ind1[n]] - beta_ind[jj_ind1[n]]));
    l1 = l1 + log(p1) * y_ind1[n] + log(1 - p1) * (1 - y_ind1[n]);
  }

  // log likelihoods for individual responses, partner 2
  for (n in 1:N_ind2) {
    p2 = inv_logit(alpha_ind[jj_ind2[n]] * (theta2[kk_ind2[n]] - beta_ind[jj_ind2[n]]));
    l2 = l2 + log(p2) * y_ind2[n] + log(1 - p2) * (1 - y_ind2[n]);
  }

  // log likelihood for group responses
  for (n in 1:N_col) {
    p1 = inv_logit(alpha_col[jj_col[n]] * (theta1[kk_col[n]] - beta_col[jj_col[n]]));
    p2 = inv_logit(alpha_col[jj_col[n]] * (theta2[kk_col[n]] - beta_col[jj_col[n]]));
    R  = w[kk_col[n]] * (p1 + p2) + (1 - 2 * w[kk_col[n]]) * p1 * p2;
    lR = lR + log(R)  * y_col[n] + log(1 - R) * (1 - y_col[n]);
  }
  target += l1 + l2 + lR;
}
