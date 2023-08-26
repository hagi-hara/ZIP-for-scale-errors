data {
  int N;
  vector<lower=0>[N] Y;
  vector[N] age;
  vector[N] age2;
  vector[N] session;
  vector[N] obj;
  vector[N] gender;
  int N_age_pred;
  vector[N_age_pred] age_pred;
  vector[N_age_pred] age2_pred;
}

parameters {
  real beta0;
  real beta1;
  real beta2;
  real beta3;
  real beta4;
  real beta5;
  real<lower=0> sigma;
}

model {
  Y ~ lognormal(beta0 + beta1 * age + beta2 * age2 + beta3 * session + beta4 * obj + beta5 * gender, sigma);
  beta0 ~ student_t(3, 0, 1);
  beta1 ~ student_t(3, 0, 1);
  beta2 ~ student_t(3, 0, 1);
  beta3 ~ student_t(3, 0, 1);
  beta4 ~ student_t(3, 0, 1);
  beta5 ~ student_t(3, 0, 1);
  sigma ~ student_t(3, 0, 1);
}

generated quantities {
  vector[N_age_pred] Y_pred;
  vector[N] Y_dist;
  vector[N] log_lik;

  for(i in 1:N_age_pred) {
    Y_pred[i] = exp(beta0 + beta1 * age_pred[i] + beta2 * age2_pred[i] + beta5 * 1/2);
  }
  
  for(i in 1:N) {
    Y_dist[i] = lognormal_rng(beta0 + beta1 * age[i] + beta2 * age2[i] + beta3 * session[i] + beta4 * obj[i] + beta5 * gender[i], sigma);
  }
  
  for(i in 1:N) {
    log_lik[i] = lognormal_lpdf(Y[i] | beta0 + beta1 * age[i] + beta2 * age2[i] + beta3 * session[i] + beta4 * obj[i] + beta5 * gender[i], sigma);
  }
}
