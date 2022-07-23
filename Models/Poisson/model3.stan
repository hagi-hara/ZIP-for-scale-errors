data {
  int N_row;
  int N_basis;
  int<lower=0> Y[N_row];
  matrix[N_basis, N_row] B;
  vector[N_row] country;
  vector[N_row] session;
  int N_age_pred;
  matrix[N_basis, N_age_pred] B_pred;
}

parameters {
  row_vector[N_basis] p_beta;
  real p_beta_country;
  real p_beta_session;
}

model {
  Y ~ poisson_log(to_vector(p_beta*B) + p_beta_country * country + p_beta_session * session);
  p_beta ~ student_t(3, 0, 1);
  p_beta_country ~ student_t(3, 0, 1);
  p_beta_session ~ student_t(3, 0, 1);
}

generated quantities {
  vector[N_age_pred] Y_pred_poi;
  vector[N_row] Y_dist;
  vector[N_row] log_lik;
  
  Y_pred_poi = exp(to_vector(p_beta*B_pred) + p_beta_country * 0.5);
  
  for(i in 1:N_row) {
    Y_dist[i] = poisson_log_rng(to_vector(p_beta*B)[i] + p_beta_country * country[i] + p_beta_session * session[i]);
  }
  
  for(i in 1:N_row) {
    log_lik[i] = poisson_log_lpmf(Y[i] | to_vector(p_beta*B)[i] + p_beta_country * country[i] + p_beta_session * session[i]);
  }
}
