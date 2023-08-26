data {
  int N;
  int<lower=0> Y[N];
  vector[N] age;
  // vector[N] age2;
  vector[N] countUK;
  vector[N] countUS;
  vector[N] session;
  vector[N] obj;
  vector[N] dur;
  vector[N] gender;
  int N_age_pred;
  vector[N_age_pred] age_pred;
  vector[N_age_pred] age2_pred;
}

parameters {
  real p_beta0;
  real p_beta1;
  // real p_beta2;
  real p_beta3;
  real p_beta4;
  real p_beta5;
  real p_beta6;
  real p_beta7;
  real p_beta8;
}

model {
  Y ~ poisson_log(p_beta0 + p_beta1 * age + p_beta3 * countUK + p_beta4 * countUS + p_beta5 * session + p_beta6 * obj + p_beta7 * dur + p_beta8 * gender);
  p_beta0 ~ student_t(3, 0, 1);
  p_beta1 ~ student_t(3, 0, 1);
  // p_beta2 ~ student_t(3, 0, 1);
  p_beta3 ~ student_t(3, 0, 1);
  p_beta4 ~ student_t(3, 0, 1);
  p_beta5 ~ student_t(3, 0, 1);
  p_beta6 ~ student_t(3, 0, 1);
  p_beta7 ~ student_t(3, 0, 1);
  p_beta8 ~ student_t(3, 0, 1);
}

generated quantities {
  vector[N_age_pred] Y_pred_poi;
  vector[N] Y_dist;
  vector[N] log_lik;
  
  for(i in 1:N_age_pred) {
    Y_pred_poi[i] = exp(p_beta0 + p_beta1 * age_pred[i] + p_beta3 * 1/3 + p_beta4 * 1/3 + p_beta8 * 1/2);
  }
  
  for(i in 1:N) {
    Y_dist[i] = poisson_log_rng(p_beta0 + p_beta1 * age[i] + p_beta3 * countUK[i] + p_beta4 * countUS[i] + p_beta5 * session[i] + p_beta6 * obj[i] + p_beta7 * dur[i] + p_beta8 * gender[i]);
  }
  
  for(i in 1:N){
    log_lik[i] = poisson_log_lpmf(Y[i] | p_beta0 + p_beta1 * age[i] + p_beta3 * countUK[i] + p_beta4 * countUS[i] + p_beta5 * session[i] + p_beta6 * obj[i] + p_beta7 * dur[i] + p_beta8 * gender[i]);
  }
}
