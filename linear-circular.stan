functions {
  real kernel_estimate(real time, vector value_times, vector values, real bandwidth) {
    int num_values = num_elements(values);
    vector[num_values] weights = exp(-0.5* square((value_times - time)  / bandwidth));
    return sum(weights .* values) / sum(weights);
  }
}

data {
  int num_cells;
  int num_time;

  vector<lower=0>[num_cells] regulator_expression;
  vector<lower=0>[num_cells] target_expression;

  real gp_sigma;
  real gp_length;

  real expression_sigma;
  real w;
  real<lower=0> decay;
  real<lower = 0> initial_condition;

  real bandwidth;

}


transformed data {
  vector[num_cells] zeroes;
    vector[num_time] times;

  for(i in 1:num_cells){
    zeroes[i] = 0;
  }

  for(i in 1:num_time) {
    times[i] = i - 1;
  }

}


parameters {
  vector<lower = 0, upper = 1>[num_cells] time_position;
}

transformed parameters {
  vector[num_time] regulator_profile;
  vector[num_time] target_profile;

  for(t in 1:num_time) {
    regulator_profile[t] = kernel_estimate(t, time_position * (num_time - 1), regulator_expression, inv(num_time - 1));
  }

  //gene RNA synthesis
  { //new scope to make the variables local
    real decay_per_step = exp(-decay);
    real residual;
    vector[num_time] synthesis = w * regulator_profile;

    target_profile[1] = initial_condition;

    //Calculating the integral by trapezoid rule in a single pass for all values
    residual = -0.5 * synthesis[1];
    for (time in 2:num_time)
    {
      residual = (residual + synthesis[time - 1]) * decay_per_step;

      { //new block to let me define new vars
        real integral_value = residual + 0.5 * synthesis[time];
        target_profile[time] = initial_condition * exp(-decay * (time - 1)) + integral_value;
      }
    }
  }



}

model {
  vector[num_cells] predicted_target;

  for(c in 1:num_cells) {
    real time_value = time_position[c] * (num_time - 1);
    predicted_target[c] = kernel_estimate(time_value, times, target_profile, bandwidth);
  }

  target_expression ~ normal_lpdf(predicted_target, expression_sigma);

  {
    matrix[num_cells,num_cells] cov_m;
    //matrix[num_cells,num_cells] cov_m_chol;
    vector[num_cells] inverse_transformed_regulator;
    for(i in 1:num_cells) {
      inverse_transformed_regulator[i] = log_diff_exp(regulator_expression[i], 0); //stands for ln(e^x - 1)
    }

    for(i in 1:num_cells) {
      cov_m[i,i] = square(gp_sigma) + 0.000001;
      for(j in 1:num_cells) {
        real dist = fabs(time_position[i] - time_position[j]);
        real  covariance = square(gp_sigma) * exp(-2* square(sin(pi() * dist)) / square(gp_length));
        cov_m[i,j] = covariance;
        cov_m[j,i] = covariance;
      }
    }
    //cov_m_chol = cholesky_decompose(cov_m);

    target += multi_normal_lpdf(inverse_transformed_regulator | zeroes, cov_m);
    target += regulator_expression - inverse_transformed_regulator; //the log Jacobian of ln(e^x - 1)
  }

  // time_position; is explicitly uniform, no prior here
}

