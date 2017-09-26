data {
  int num_cells;
  int num_regulators;
  int num_targets;
  int<lower=0> regulator_expression[num_regulators,num_cells];
  int<lower=0> target_expression[num_targets,num_cells];
  real<lower=0> w_prior_sigma;
  //real<lower=0> derivative_tau;

  real<lower=0> regulator_prior_sigma;
  real<lower=0> regulator_prior_mean;
  //real<lower=0> residual_prior_sigma;
}

parameters {
  matrix<lower=0>[num_regulators, num_targets] w;
  matrix<lower=0>[num_cells, num_regulators] regulator_expression_mean;
  //real<lower=0> residual; //intended for 1/dispersion, but maybe not really working
}

transformed parameters {
}


model {
  matrix[num_cells, num_targets] target_expression_mean = regulator_expression_mean * w;// + derivatives_noise;
  //residual ~ normal(0,residual_prior_sigma);
  // int do_print = 0;
  // for(c in 1:num_cells) {
  //   for(t in 1:num_targets) {
  //     if(is_nan(target_expression_mean[c,t])){
  //       print("Nan at",c,":",t);
  //       do_print = 1;
  //     }
  //   }
  // }
  // if(do_print == 1){
  //   print(target_expression_mean);
  //   print(w);
  //   print(regulator_expression_mean);
  // }
  //poisson, because equilibrium over time (and am interested in protein which should smooth time out)
  for(reg in 1:num_regulators) {
    regulator_expression[reg,] ~ poisson(regulator_expression_mean[,reg]);
  }

  //poisson because equilibrium over cells
  for(target in 1:num_targets) {
    target_expression[target,] ~ poisson(target_expression_mean[,target]);//neg_binomial_2(target_expression_mean, 1/residual);
  }

  for(reg in 1:num_regulators) {
    regulator_expression_mean[,reg] ~ lognormal(regulator_prior_mean, regulator_prior_sigma);
  }

  for(target in 1:num_targets){
    w[,target] ~ normal(0, w_prior_sigma);
  }
}
