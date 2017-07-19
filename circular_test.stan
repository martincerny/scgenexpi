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
}

parameters {
  vector<lower = 0, upper = 1>[num_cells] time_position;
  //real<lower=0> time_shape_1;
  //real<lower=0> time_shape_2;
}

transformed parameters {
  vector<lower=0>[num_cells] regulator_estimate;
  vector<lower=0>[num_cells] target_estimate;

  {
    int ordering[num_cells];
    real regulator_zero;
    real integral;

    for(i in 1:num_cells) {
      ordering[i] = i;
    }
    for(i in 1:num_cells) {
      for(j in 2:num_cells) {
        if(time_position[ordering[j - 1]] > time_position[ordering[j]]) {
          int x = ordering[j];
          ordering[j] = ordering[j - 1];
          ordering[j - 1] = x;
        }
      }
    }


    for(c in 1:num_cells) {
      regulator_estimate[c] = kernel_estimate(time_position[c], time_position, regulator_expression, regulator_bandwidth);
    }
    regulator_zero = 0;

    integral = 0.5 * (regulator_estimate[ordering[1]] + regulator_zero) * time_position[ordering[1]];
    target_estimate[ordering[1]] = integral;

    for(i in 2:num_cells) {
      {
        int cell_i = ordering[i];
        int cell_before = ordering[i-1];
        integral = integral + 0.5 * (regulator_estimate[cell_i] + regulator_estimate[cell_before]) * (time_position[cell_i] - time_position[cell_before]);
        target_estimate[cell_i] = integral;
      }
    }


    if(debug > 0) {
      print("Time: ", time_position)
      print("Ordering: ", ordering)
      print("Regulaotr: ", regulator_expression)
      print("RegE: ", regulator_estimate)
      print("Target: ", target_expression)
      print("TargetE: ", target_estimate)
    }
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
