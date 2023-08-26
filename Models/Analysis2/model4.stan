data {
  int N;
  int<lower=0> Y[N];
  vector[N] age;
  vector[N] age2;
  vector[N] countUK;
  vector[N] session;
  vector[N] obj;
  vector[N] dur;
  vector[N] gender;
  int N_age_pred;
  vector[N_age_pred] age_pred;
  vector[N_age_pred] age2_pred;
}

parameters {
  real b_beta0;
  real b_beta1;
  // real b_beta2;
  real b_beta3;
  real b_beta4;
  real b_beta5;
  real b_beta6;
  real b_beta7;
  real p_beta0;
  real p_beta1;
  // real p_beta2;
  real p_beta3;
  real p_beta4;
  real p_beta5;
  real p_beta6;
  real p_beta7;
}

transformed parameters{
  vector[N] theta;
  vector[N] lambda;
  
  for(n in 1:N) {
    theta[n] = inv_logit(b_beta0 + b_beta1 * age[n] + b_beta3 * countUK[n] + b_beta4 * session[n] + b_beta5 * obj[n] + b_beta6 * dur[n] + b_beta7 * gender[n]);
    lambda[n] = exp(p_beta0 + p_beta1 * age[n] + p_beta3 * countUK[n] + p_beta4 * session[n] + p_beta5 * obj[n] + p_beta6 * dur[n] + p_beta7 * gender[n]);
  }
}

model {
  for(n in 1:N) {
    if (Y[n] == 0){
      target += log_sum_exp( bernoulli_lpmf(0|theta[n]), bernoulli_lpmf(1|theta[n]) + poisson_lpmf(0|lambda[n]) );
    } else {
      target += bernoulli_lpmf(1|theta[n]) + poisson_lpmf(Y[n]|lambda[n]);
    }
  }
  b_beta0 ~ student_t(3, 0, 1);
  b_beta1 ~ student_t(3, 0, 1);
  // b_beta2 ~ student_t(3, 0, 1);
  b_beta3 ~ student_t(3, 0, 1);
  b_beta4 ~ student_t(3, 0, 1);
  b_beta5 ~ student_t(3, 0, 1);
  b_beta6 ~ student_t(3, 0, 1);
  b_beta7 ~ student_t(3, 0, 1);
  p_beta0 ~ student_t(3, 0, 1);
  p_beta1 ~ student_t(3, 0, 1);
  // p_beta2 ~ student_t(3, 0, 1);
  p_beta3 ~ student_t(3, 0, 1);
  p_beta4 ~ student_t(3, 0, 1);
  p_beta5 ~ student_t(3, 0, 1);
  p_beta6 ~ student_t(3, 0, 1);
  p_beta7 ~ student_t(3, 0, 1);
}

generated quantities {
  vector[N_age_pred] Y_pred_ber;
  vector[N_age_pred] Y_pred_poi;
  vector[N_age_pred] Y_pred_all;
  vector[N] Y_dist;
  vector[N] log_lik;
  
  for(i in 1:N_age_pred) {
    Y_pred_ber[i] = inv_logit(b_beta0 + b_beta1 * age_pred[i] + b_beta3 * 1/2 + b_beta7 * 1/2);
    Y_pred_poi[i] = exp(p_beta0 + p_beta1 * age_pred[i] + p_beta3 * 1/2 + p_beta7 * 1/2);
    Y_pred_all[i] = Y_pred_ber[i] * Y_pred_poi[i];
  }
  
  for(i in 1:N) {
    Y_dist[i] = bernoulli_rng(theta[i]) * poisson_rng(lambda[i]);
  }
  
  for(i in 1:N) {
    log_lik[i] = (Y[i] == 0) ? 
    (log_sum_exp( bernoulli_lpmf(0|theta[i]), bernoulli_lpmf(1|theta[i]) + poisson_lpmf(0|lambda[i]))):
    (bernoulli_lpmf(1|theta[i]) + poisson_lpmf(Y[i]|lambda[i]));
  }
}
