//Spline code thanks to https://github.com/milkha/Splines_in_Stan/blob/master/b_spline_penalized.stan

functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t))
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}

data {
  int num_cells;
  int num_knots;
  vector[num_knots] knots;
  int spline_degree;
  real regulator_expression[num_cells];
  real target_expression[num_cells];

  real coeff_sigma;
  real expression_sigma;

  real time_position_prior_alpha[num_cells];
  real time_position_prior_beta[num_cells];
}

transformed data {
  int num_basis = num_knots + spline_degree - 1;
}

parameters {
  //row_vector[num_basis] a_raw;
  row_vector[num_basis] coeffs;
  real<lower=0,upper=1> time_position[num_cells];
  /*real a0;
  real<lower=0> sigma;
  real<lower=0> tau;*/
}

transformed parameters {
  //row_vector[num_basis] a;
  vector[num_cells] regulator_true;
  /*a[1] = a_raw[1];
  for (i in 2:num_basis)
    a[i] = a[i-1] + a_raw[i]*tau;
    */

  {
    matrix[num_basis, num_cells] B;

    vector[spline_degree + num_knots] ext_knots_temp;
    vector[2*spline_degree + num_knots] ext_knots;
    ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
    ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
    for (ind in 1:num_basis)
      B[ind,:] = to_row_vector(build_b_spline(time_position, to_array_1d(ext_knots), ind, spline_degree + 1));

    //Hacky adjustment to work with time equal to the last knot
    for(i in 1:num_cells) {
        if(time_position[i] == ext_knots[num_knots + 2*spline_degree - 1]) {
          B[num_basis, i] = 1;
        }
    }

    regulator_true = /*a0*to_vector(X) +*/ to_vector(coeffs * B);
  }

}

model {
  //a_raw ~ normal(0, 1);
  coeffs ~ normal(0, coeff_sigma);
  time_position ~ beta(time_position_prior_alpha, time_position_prior_beta);
  //tau ~ normal(0, 1);
  //sigma ~ cauchy(0, 1);
  regulator_expression ~ normal(regulator_true, expression_sigma);
  target_expression ~ normal(time_position, expression_sigma);
}
