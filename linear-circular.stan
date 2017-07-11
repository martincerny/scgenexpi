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
  matrix[num_time,num_time] cov_m;
  matrix[num_time,num_time] cov_m_chol;
  vector[num_time] times;

  for(i in 1:num_time) {
    cov_m[i,i] = square(gp_sigma) + 0.000001;
    for(j in 1:num_time) {
      real  covariance = square(gp_sigma) * exp(-2* square(sin(pi() * abs(i - j) / num_time)) / square(gp_length));
      cov_m[i,j] = covariance;
      cov_m[j,i] = covariance;
    }
  }
  cov_m_chol = cholesky_decompose(cov_m);

  for(i in 1:num_time) {
    times[i] = i - 1;
  }
}

parameters {
  vector[num_time] regulator_profile_raw;
  vector<lower = 0, upper = 1>[num_cells] time_position;
}

transformed parameters {
  vector[num_time] regulator_profile = log1p_exp(cov_m_chol * regulator_profile_raw);
  vector[num_time] target_profile;

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
  vector[num_cells] predicted_regulator;
  vector[num_cells] predicted_target;

  for(c in 1:num_cells) {
    real time_value = time_position[c] * (num_time - 1);
    predicted_regulator[c] = kernel_estimate(time_value, times, regulator_profile, bandwidth);
    predicted_target[c] = kernel_estimate(time_value, times, regulator_profile, bandwidth);
  }

  regulator_expression ~ normal_lpdf(predicted_regulator, expression_sigma);
  target_expression ~ normal_lpdf(predicted_target, expression_sigma);
  regulator_profile_raw ~ normal(0,1);
  // time_position; is explicitly uniform, no prior here
}

