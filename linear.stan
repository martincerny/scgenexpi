data {
  int num_cells;
  int num_samples;
//  int num_genes;
//  int num_regulations;
//  int regulators[num_regulations];
//  int targets[num_regulations];
  matrix<lower=0>[num_cells, num_samples] regulator_expression;
  matrix<lower=0>[num_cells, num_samples] target_expression;

/*  int num_steps;
  int step_starts[num_steps];
  int step_ends[num_steps];*/
}

parameters {
  real w;
  real<lower=0> decay;
}

