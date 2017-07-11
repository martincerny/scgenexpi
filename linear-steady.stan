data {
  int num_cells;

  vector<lower=0>[num_cells] regulator_expression;
  vector<lower=0>[num_cells] target_expression;

  real<lower=0> bandwidth;
}

parameters {
  real w;
}

model {
  vector[num_cells] synthesis = (regulator_expression * w);
  vector[num_cells] predicted_target = synthesis + exp(-1) * (target_expression - synthesis);

  for(c in 1:num_cells) {
    vector[num_cells] components = rep_vector(0, num_cells);
    for(c2 in 1:num_cells) {
      components[c2] = normal_lpdf(target_expression[c] | predicted_target[c2] , bandwidth);
    }
    target += log_sum_exp(components);
  }

  w ~ normal(0,3);
}

