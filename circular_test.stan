functions {
  real kernel_estimate(real time, vector value_times, vector values, real bandwidth) {
    int num_values = num_elements(values);
    vector[num_values] weights = exp(-0.5* square((value_times - time)  / bandwidth));
    return sum(weights .* values) / sum(weights);
  }
}

data {
  int num_cells;
  int num_integration_points;
  vector<lower=0>[num_cells] regulator_expression;
  vector<lower=0>[num_cells] target_expression;
  real expression_sigma;
  real regulator_bandwidth;
  real integrated_target_bandwidth;
  int debug;
}

transformed data {
  vector<lower = 0>[num_integration_points] integration_times;
  real integration_step = inv(num_integration_points - 1);
  for(i in 1:num_integration_points) {
    integration_times[i] = (i-1) * integration_step;
  }
  if(debug > 0) {
    print(integration_times);
  }
}

parameters {
  vector<lower = 0, upper = 1>[num_cells] time_position;
  //real<lower=0> time_shape_1;
  //real<lower=0> time_shape_2;
}

transformed parameters {
  vector[num_integration_points] regulator_estimate;
  vector[num_integration_points] integrated_values;
  vector[num_cells] target_estimate;
  real residual = 0;
  real first_term;

  integrated_values[1] = 0;
  regulator_estimate[1] = kernel_estimate(integration_times[1], time_position, regulator_expression, regulator_bandwidth);

  first_term = 0.5 * regulator_estimate[1] * integration_step;

  for(i in 2:num_integration_points) {
    regulator_estimate[i] = kernel_estimate(integration_times[i], time_position, regulator_expression, regulator_bandwidth);
    {
      real current_step = regulator_estimate[i] * integration_step;
      integrated_values[i] = first_term + residual + 0.5 * current_step;
      residual = residual + current_step;
    }
  }

  for(c in 1:num_cells) {
    target_estimate[c] = kernel_estimate(time_position[c], integration_times, integrated_values, integrated_target_bandwidth);
  }
  if(debug > 0) {
    print("Time: ", time_position)
    print("Regulaotr: ", regulator_expression)
    print("RegE: ", regulator_estimate)
    print("Integrated: ", integrated_values)
    print("TargetE: ", target_estimate)
    print("Target: ", target_expression)
  }
}

model {
  target_expression ~ normal(target_estimate, expression_sigma);
  /*
  time_position ~ beta(time_shape_1, time_shape_2);
  time_shape_1 ~ cauchy(1, 2.5);
  time_shape_2 ~ cauchy(1, 2.5);
  */
 //regulator_expression ~ normal(w * time_position, expression_sigma);
}
