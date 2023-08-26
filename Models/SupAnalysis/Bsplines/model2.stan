data {
  int N_row;
  int N_basis;
  int<lower=0> Y[N_row];
  matrix[N_basis, N_row] B;
  vector[N_row] countUK;
  vector[N_row] countUS;
  vector[N_row] session;
  vector[N_row] obj;
  vector[N_row] dur;
  vector[N_row] gender;
  int N_age_pred;
  matrix[N_basis, N_age_pred] B_pred;
}

parameters {
  row_vector[N_basis] p_beta;
  real p_beta_countUK;
  real p_beta_countUS;
  real p_beta_session;
  real p_beta_obj;
  real p_beta_dur;
  real p_beta_gender;
}

model {
  Y ~ poisson_log(to_vector(p_beta*B) + p_beta_countUK * countUK + p_beta_countUS * countUS + p_beta_session * session + p_beta_obj * obj + p_beta_dur * dur + p_beta_gender * gender);
  p_beta ~ student_t(3, 0, 1);
  p_beta_countUK ~ student_t(3, 0, 1);
  p_beta_countUS ~ student_t(3, 0, 1);
  p_beta_session ~ student_t(3, 0, 1);
  p_beta_obj ~ student_t(3, 0, 1);
  p_beta_dur ~ student_t(3, 0, 1);
  p_beta_gender ~ student_t(3, 0, 1);
}

generated quantities {
  vector[N_age_pred] Y_pred_poi;
  vector[N_row] Y_dist;
  vector[N_row] log_lik;
  
  Y_pred_poi = exp(to_vector(p_beta*B_pred) + p_beta_countUK * 1/3 + p_beta_countUS * 1/3 + p_beta_gender * 1/2);
  
  for(i in 1:N_row) {
    Y_dist[i] = poisson_log_rng(to_vector(p_beta*B)[i] + p_beta_countUK * countUK[i] + p_beta_countUS * countUS[i] + p_beta_session * session[i] + p_beta_obj * obj[i] + p_beta_dur * dur[i] + p_beta_gender * gender[i]);
  }
  
  for(i in 1:N_row) {
    log_lik[i] = poisson_log_lpmf(Y[i] | to_vector(p_beta*B)[i] + p_beta_countUK * countUK[i] + p_beta_countUS * countUS[i] + p_beta_session * session[i] + p_beta_obj * obj[i] + p_beta_dur * dur[i] + p_beta_gender * gender[i]);
  }
}
