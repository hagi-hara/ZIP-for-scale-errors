data {
  int N;
  vector<lower=0>[N] Y;
  vector[N] age;
  // vector[N] age2;
  vector[N] ec;
  vector[N] session;
  vector[N] gender;
}

parameters {
  real beta0;
  real beta1;
  // real beta2;
  real beta3;
  real beta4;
  real beta5;
  real<lower=0> sigma;
}

model {
  Y ~ lognormal(beta0 + beta1 * age + beta3 * ec + beta4 * session + beta5 * gender, sigma);
  beta0 ~ student_t(3, 0, 1);
  beta1 ~ student_t(3, 0, 1);
  // beta2 ~ student_t(3, 0, 1);
  beta3 ~ student_t(3, 0, 1);
  beta4 ~ student_t(3, 0, 1);
  beta5 ~ student_t(3, 0, 1);
  sigma ~ student_t(3, 0, 1);
}

generated quantities {
  vector[N] Y_dist;
  vector[N] log_lik;

  for(i in 1:N) {
    Y_dist[i] = lognormal_rng(beta0 + beta1 * age[i] + beta3 * ec[i] + beta4 * session[i] + beta5 * gender[i], sigma);
  }
  
  for(i in 1:N) {
    log_lik[i] = lognormal_lpdf(Y[i] | beta0 + beta1 * age[i] + beta3 * ec[i] + beta4 * session[i] + beta5 * gender[i], sigma);
  }
}
