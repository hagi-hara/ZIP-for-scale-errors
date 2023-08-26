data {
  int N_row;
  int N_basis;
  int<lower=0> Y[N_row];
  matrix[N_basis, N_row] B;
  vector[N_row] countUK;
  vector[N_row] session;
  vector[N_row] obj;
  vector[N_row] dur;
  vector[N_row] gender;
  int N_age_pred;
  matrix[N_basis, N_age_pred] B_pred;
}

parameters {
  row_vector[N_basis] b_beta;
  row_vector[N_basis] p_beta;
  real b_beta_countUK;
  real b_beta_session;
  real b_beta_obj;
  real b_beta_dur;
  real b_beta_gender;
  real p_beta_countUK;
  real p_beta_session;
  real p_beta_obj;
  real p_beta_dur;
  real p_beta_gender;
}

transformed parameters{
  vector[N_row] theta;
  vector[N_row] lambda;
  
  theta = inv_logit(to_vector(b_beta*B) + b_beta_countUK * countUK + b_beta_session * session + b_beta_obj * obj + b_beta_dur * dur + b_beta_gender * gender);
  lambda = exp(to_vector(p_beta*B) + p_beta_countUK * countUK + p_beta_session * session + p_beta_obj * obj + p_beta_dur * dur + p_beta_gender * gender);
}

model {
  for(n in 1:N_row) {
    if (Y[n] == 0){
      target += log_sum_exp( bernoulli_lpmf(0|theta[n]), bernoulli_lpmf(1|theta[n]) + poisson_lpmf(0|lambda[n]) );
    } else {
      target += bernoulli_lpmf(1|theta[n]) + poisson_lpmf(Y[n]|lambda[n]);
    }
  }
  b_beta ~ student_t(3, 0, 1);
  p_beta ~ student_t(3, 0, 1);
  b_beta_countUK ~ student_t(3, 0, 1);
  b_beta_session ~ student_t(3, 0, 1);
  b_beta_obj ~ student_t(3, 0, 1);
  b_beta_dur ~ student_t(3, 0, 1);
  b_beta_gender ~ student_t(3, 0, 1);
  p_beta_countUK ~ student_t(3, 0, 1);
  p_beta_session ~ student_t(3, 0, 1);
  p_beta_obj ~ student_t(3, 0, 1);
  p_beta_dur ~ student_t(3, 0, 1);
  p_beta_gender ~ student_t(3, 0, 1);
}

generated quantities {
  vector[N_age_pred] Y_pred_ber;
  vector[N_age_pred] Y_pred_poi;
  vector[N_age_pred] Y_pred_all;
  vector[N_row] Y_dist;
  vector[N_row] log_lik;
  
  Y_pred_ber = inv_logit(to_vector(b_beta*B_pred) + b_beta_countUK * 1/2 + b_beta_gender * 1/2);
  Y_pred_poi = exp(to_vector(p_beta*B_pred) + p_beta_countUK * 1/2 + p_beta_gender * 1/2);
  
  for(i in 1:N_age_pred) {
    Y_pred_all[i] = Y_pred_ber[i] * Y_pred_poi[i];
  }
  
  for(i in 1:N_row) {
    Y_dist[i] = bernoulli_rng(theta[i]) * poisson_rng(lambda[i]);
  }
  
  for(i in 1:N_row) {
    log_lik[i] = (Y[i] == 0) ? 
    (log_sum_exp( bernoulli_lpmf(0|theta[i]), bernoulli_lpmf(1|theta[i]) + poisson_lpmf(0|lambda[i]))):
    (bernoulli_lpmf(1|theta[i]) + poisson_lpmf(Y[i]|lambda[i]));
  }
}
