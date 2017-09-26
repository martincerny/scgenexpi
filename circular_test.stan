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
  //vector<lower=0>[num_cells] regulator_expression;
  //vector<lower=0>[num_cells] target_expression;
  vector[num_cells] regulator_expression;
  vector[num_cells] target_expression;
  real expression_sigma;
  real regulator_bandwidth;
  int debug;
  real gp_sigma;
  real gp_length;
  real gp_mean;
  //real time_shape_1;
  //real time_shape_2;
  real time_position_prior_alpha[num_cells];
  real time_position_prior_beta[num_cells];


}

transformed data {
/*
  matrix[num_cells,num_cells] dist;
  for(i in 1:num_cells) {
    dist[i,i] = 0;
    for(j in (i + 1):num_cells) {
      //real d = sqrt(square(regulator_expression[i] - regulator_expression[j]) + square(target_expression[i] - target_expression[j]));
      //real d = fabs(regulator_expression[i] - regulator_expression[j]) + fabs(target_expression[i] - target_expression[j]);
      real d = fabs(regulator_expression[i] - regulator_expression[j]);
      dist[i,j] = d;
      dist[j,i] = d;
    }
  }
*/


}

parameters {
/*  real<lower=0> dist_bias;
  real<lower=0> dist_coeff;
  real<lower=0.1> dist_sigma;
  real<lower=0> time_shape_1;
  real<lower=0> time_shape_2;

  */

  real<lower = 0, upper = 1> time_position[num_cells];

  vector[num_cells] regulator_true_raw;
}

transformed parameters {
  //vector<lower=0>[num_cells] regulator_estimate;
  //vector<lower=0>[num_cells] target_estimate;
  vector[num_cells] regulator_true;
  //vector[num_cells] regulator_estimate;
  vector[num_cells] target_estimate;

  {
    matrix[num_cells, num_cells] cov_m = cov_exp_quad(time_position, gp_sigma, gp_length);
    matrix[num_cells, num_cells] cov_m_chol;
    //int ordering[num_cells];
    real regulator_zero;
    real integral;
/*
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
*/
    for(i in 1:num_cells) {
      cov_m[i,i] = cov_m[i,i] + 1e-6;
    }

    cov_m_chol = cholesky_decompose(cov_m);

    regulator_true = cov_m_chol * regulator_true_raw;

    target_estimate = to_vector(time_position);
/*
    for(c in 1:num_cells) {
      regulator_estimate[c] = kernel_estimate(time_position[c], to_vector(time_position), regulator_true, regulator_bandwidth);
    }*/
    //regulator_estimate = regulator_true;
    //regulator_zero = 0;

/*
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
*/

    if(debug > 0) {
      print("Time: ", time_position)
      //print("Ordering: ", ordering)
      //print("Regulaotr: ", regulator_expression)
      //print("RegE: ", regulator_estimate)
      //print("Target: ", target_expression)
      //print("TargetE: ", target_estimate)
    }
  }

}

model {
  regulator_expression ~ normal(regulator_true, expression_sigma);
  target_expression ~ normal(target_estimate, expression_sigma);
  /*
  for(i in 1:num_cells) {
    for(j in (i + 1):num_cells) {
      real delta = dist_bias + dist_coeff * fabs(time_position[i] - time_position[j]);
      real sig_sq = square(dist_sigma);
      real alpha = square(delta) / sig_sq;
      real beta = delta/sig_sq;
      dist[i,j] ~ gamma(alpha, beta);
    }
  }*/

  regulator_true_raw ~ normal(0,1);

  //time_position ~ beta(time_shape_1, time_shape_2);
  time_position ~ beta(time_position_prior_alpha, time_position_prior_beta);
  /*
  dist_bias ~ cauchy(0, 2.5);
  dist_coeff ~ cauchy(1, 2.5);
  dist_sigma ~ cauchy(0, 2.5);
  time_position ~ beta(time_shape_1, time_shape_2);
  time_shape_1 ~ cauchy(1, 2.5);
  time_shape_2 ~ cauchy(1, 2.5);
  */
 //regulator_expression ~ normal(w * time_position, expression_sigma);
}
