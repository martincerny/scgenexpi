data {
  int num_cells;
  vector[num_cells] regulator_expression;
  vector[num_cells] target_expression;
  real expression_sigma;
  real gp_sigma;
  real gp_length;
  real time_shape_1;
  real time_shape_2;
}


parameters {
  real<lower = 0, upper = 1> time_position[num_cells];
  vector[num_cells] regulator_true_raw;
}

transformed parameters {
  vector[num_cells] regulator_true;
  vector[num_cells] target_estimate;

  {
    //new scope to introduce vars that are not stored in output
    matrix[num_cells, num_cells] cov_m = cov_exp_quad(time_position, gp_sigma, gp_length);
    matrix[num_cells, num_cells] cov_m_chol;

    //Add jitter to ensure positive semidefiniteness
    for(i in 1:num_cells) {
      cov_m[i,i] = cov_m[i,i] + 1e-6;
    }

    cov_m_chol = cholesky_decompose(cov_m);

    regulator_true = cov_m_chol * regulator_true_raw;

    target_estimate = to_vector(time_position);
  }

}

model {
  regulator_expression ~ normal(regulator_true, expression_sigma);
  target_expression ~ normal(target_estimate, expression_sigma);

  regulator_true_raw ~ normal(0,1);

  time_position ~ beta(time_shape_1, time_shape_2);
}

