require(deSolve)

generate_random_profile <- function(time, scale, length, mean_func = 0, periodic = FALSE, period = 1) {
  # Construct the squared exponential covariance matrix
  cov_matrix = array(0,c(length(time), length(time)));
  maxTime = max(time)
  minTime = min(time)
  for(i in 1:length(time))  {
    cov_matrix[i,i] = scale^2 + 0.00000001; #Adding jitter to ensure positive semi-definiteness
    if (i < length(time)) {
      for(j in (i+1):length(time)) {
        distance = (time[i] - time[j])
        if(periodic) {
          covariance = (scale^2) * exp(-2*(sin(3.141592 * abs(distance) / period)^2) / (length ^ 2))
        } else {
          covariance = (scale^2) * exp( (-0.5 / (length ^ 2)) * (distance ^ 2) );
        }
        cov_matrix[i,j] = covariance
        cov_matrix[j,i] = covariance
      }
    }
  }
  # Draw from the multinormal distribution using the cholesky decomposition of the cov. matrix
  chol_cov = t(chol(cov_matrix));
  raw_profile = chol_cov %*% rnorm(length(time)) + mean_func;
  # Transform to strictly positive values
  positive_profile = log1p(exp(raw_profile))
  return(t(positive_profile))
}

plot_random_profiles <- function(n, time, scale, length, true_time = time, true_profile = NULL, mean_func = 0, periodic = FALSE, period = 1, ...) {
  profiles = array(0, c(n, length(time)));
  for(i in 1:n) {
    profiles[i,] = generate_random_profile(time, scale, length, mean_func = mean_func, periodic = periodic, period = period);
  }

  ymax = max(profiles)
  ymin = min(profiles)
  if(!is.null(true_profile)) {
    ymax = max(ymax, max(true_profile))
    ymin = min(ymin, min(true_profile))
  }

  matplot(time, t(profiles), type = "l", ylim = c(ymin, ymax), ...)
  if(!is.null(true_profile)) {
    points(true_time, true_profile, pch = 19)
  }
}

linear_regulation_ode_func <- function(time, regulator_expression ) {
  regulator_fun = approxfun(time, regulator_expression,method = "linear", rule = 2);
  return(function(t, state, parms) {
      dY = parms["w"] * regulator_fun(t) - parms["decay"] * state["y"];
      list(dY);
  })
}

plot_linear <- function(init, w, decay, time, regulator_profile = generate_random_profile(time, scale = 3, length = 2)) {
  ode_state = c(y = init)
  ode_func = linear_regulation_ode_func(time - (min(time)), regulator_profile)
  ode_params = c(w = w, decay = decay);
  ode_result = ode(y = ode_state, times = time - (min(time)), func = ode_func, parms = ode_params,method = "ode45")
  plot(time, regulator_profile, type="l", ylim = c(0,max(c(regulator_profile, ode_result[,"y"]))))
  lines(time, ode_result[,"y"], col="green")
  print(cor(regulator_profile, ode_result[,"y"]))
}

simulate_linear <- function(num_cells, num_time, w = runif(1, 0.1, 3), decay = runif(1, 0.1, 1)) {
  regulator_expression = array(0, c(num_cells, num_time));
  target_expression = array(0, c(num_cells, num_time));
  ode_params = c(w = w, decay = decay);
  for(cell in 1:num_cells) {
    for(time in 1:num_time) {
      regulator_profile = generate_random_profile(0:(num_time - 1),scale = 3, length = 2)
      ode_state = c(y = runif(1,0,1))
      ode_func = linear_regulation_ode_func(0:(num_time - 1), regulator_profile)
      ode_result = ode(y = ode_state, times = 0:(num_time - 1), func = ode_func, parms = ode_params,method = "ode45")
      target_expression[cell, time] = ode_result[time,"y"]
      regulator_expression[cell, time] = regulator_profile[time]
    }
  }
  return(list(observed =
                list(
                  num_cells = num_cells,
                  num_samples = num_time,
                  bandwidth = 0.2,
                  regulator_expression = regulator_expression,
                  target_expression = target_expression
                ),
              true =
                list(
                  w = w,
                  decay = decay
                )
        ))
}

simulate_linear_circular <- function(num_cells, min_time = 5, max_time = 10, w = runif(1, 0.1, 3), decay = runif(1, 0.1, 1), init_min = 0, init_max = 1) {
  regulator_expression = array(0, num_cells);
  target_expression = array(0, num_cells);
  ode_params = c(w = w, decay = decay);
  time = seq(0, max_time, by = 0.01)
  for(cell in 1:num_cells) {
      regulator_profile = cos( time ) ^ 2
      ode_state = c(y = runif(1,init_min,init_max))
      ode_func = linear_regulation_ode_func(time, regulator_profile)
      ode_result = ode(y = ode_state, time, func = ode_func, parms = ode_params,method = "ode45")
      time_point = sample(which(time >= min_time), 1)
      target_expression[cell] = ode_result[time_point,"y"]
      regulator_expression[cell] = regulator_profile[time_point]
  }
  return(list(observed =
                list(
                  num_cells = num_cells,
                  bandwidth = 0.2,
                  regulator_expression = regulator_expression,
                  target_expression = target_expression
                ),
              true =
                list(
                  w = w,
                  decay = decay
                )
  ))
}

simulate_linear_circular_2 <- function(num_cells, num_time, scale, length, w = runif(1, 0.1, 3), decay = runif(1, 0.1, 1)) {

  time_points = runif(num_cells,0,1);
  time = time_points * (num_time - 1);

  repeat {
    regulator_profile = array(generate_random_profile(1:num_time, scale, length, periodic = TRUE, period = num_time - 1), num_time);
    if(var(regulator_profile) > 1) {
      break;
    }
  }

  initial_condition = 0.7;

  ode_params = c(w = w, decay = decay);
  ode_state = c(y = initial_condition);
  ode_func = linear_regulation_ode_func(0:(num_time - 1), regulator_profile)
  time_order = order(time)
  ode_result = ode(y = ode_state, time[time_order], func = ode_func, parms = ode_params,method = "ode45")
  ode_result_profile = ode(y = ode_state, 0:(num_time - 1), func = ode_func, parms = ode_params,method = "ode45")

  target_expression_true = array(NA, length(time))
  target_expression_true[time_order] = ode_result[,"y"];

  target_profile = ode_result_profile[,"y"];

  expression_sigma = 0.2;

  target_expression = rnorm(num_cells, target_expression_true , expression_sigma);
  target_expression[target_expression < 0] = 0;

  regulator_expression_true = approx(0:(num_time - 1), regulator_profile, time, method = "linear", rule = 2)$y
  regulator_expression = rnorm(num_cells, regulator_expression_true, expression_sigma);

  return(list(observed =
                list(
                  num_cells = num_cells,
                  num_time = num_time,
                  bandwidth = 0.1,
                  regulator_expression = regulator_expression,
                  target_expression = target_expression,
                  gp_sigma = scale,
                  gp_length = length,
                  expression_sigma = expression_sigma,
                  w = w,
                  decay = decay,
                  initial_condition = initial_condition

                  # int num_cells;
                  # int num_time;
                  #
                  # vector<lower=0>[num_cells] regulator_expression;
                  # vector<lower=0>[num_cells] target_expression;
                  #
                  # real gp_sigma;
                  # real gp_length;
                  #
                  # real expression_sigma;
                  # real w;
                  # real<lower=0> decay;
                  # real<lower = 0> initial_condition;
                  #
                  # real bandwidth;

                ),
              true =
                list(
                  time_position = time_points,
                  regulator_profile = regulator_profile,
                  target_profile = target_profile
                  #vector<lower = 0, upper = 1>[num_cells] time_position;
                  #vector[num_time] regulator_profile
                )
  ))

}

simulate_gp_test <- function(num_cells, num_time, gp_sigma) {
  full_expression = array(0, c(num_time, num_cells));
  time = seq(0,1, by = 1/(num_time - 1))
  for(cell in 1:num_cells) {
    full_expression[, cell] = generate_random_profile(time, gp_sigma, 1);
  }
}
