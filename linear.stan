data {
  int num_cells;
  int num_samples;
//  int num_genes;
//  int num_regulations;
//  int regulators[num_regulations];
//  int targets[num_regulations];
  matrix<lower=0>[num_cells, num_samples] regulator_expression;
  matrix<lower=0>[num_cells, num_samples] target_expression;

  real bandwidth; //If this should be a parameter, the density should be adjusted to include it
/*  int num_steps;
  int step_starts[num_steps];
  int step_ends[num_steps];*/
}

parameters {
  real<lower=0> w_over_decay;
  real<lower=0> decay;
}

model {
  matrix[num_cells, num_samples] synthesis_over_decay = regulator_expression * w_over_decay;
  real decay_term = exp(-decay);
  matrix[num_cells, num_samples] predicted_target = synthesis_over_decay + exp(-decay) * (target_expression - synthesis_over_decay);

  for(s in 1:(num_samples - 1)) {
    for(c in 1:num_cells) {
      vector[num_cells] components = rep_vector(0, num_cells);
      for(c2 in 1:num_cells) {
        components[c2] = normal_lpdf(target_expression[c, s + 1] | predicted_target[c2,s] , bandwidth);
      }
      //vector[num_cells] components = square((target_expression[c, s + 1] - predicted_target[, s]) / bandwidth) * -0.5;
      target += log_sum_exp(components);
    }
  }

  w_over_decay ~ normal(1, 3);
  decay ~ normal(0,5);
}

generated quantities {
  real w = w_over_decay * decay;
}
