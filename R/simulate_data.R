require(deSolve)

generate_random_profile <- function(time, scale, length) {
  # Construct the squared exponential covariance matrix
  cov_matrix = array(0,c(length(time), length(time)));
  for(i in 1:length(time))  {
    cov_matrix[i,i] = scale^2 + 0.00001;
    if (i < length(time)) {
      for(j in (i+1):length(time)) {
        covariance = (scale^2) * exp( (-0.5 / (length ^ 2)) * ((time[i] - time[j]) ^ 2) );
        cov_matrix[i,j] = covariance
        cov_matrix[j,i] = covariance
      }
    }
  }
  # Draw from the multinormal distribution using the cholesky decomposition of the cov. matrix
  chol_cov = t(chol(cov_matrix));
  raw_profile = chol_cov %*% rnorm(length(time));
  # Transform to strictly positive values
  positiveProfile = log1p(exp(raw_profile))
  return(t(positiveProfile))
}

plot_random_profiles <- function(n, time, scale, length, true_time = time, true_profile = NULL, ...) {
  profiles = array(0, c(n, length(time)));
  for(i in 1:n) {
    profiles[i,] = generate_random_profile(time, scale, length);
  }

  ymax = max(profiles)
  if(!is.null(true_profile)) {
    ymax = max(ymax, max(true_profile))
  }

  matplot(time, t(profiles), type = "l", ylim = c(0, ymax), ...)
  if(!is.null(true_profile)) {
    points(true_time, true_profile, pch = 19)
  }
}

linear_regulation_ode_func <- function(time, regulator_expression ) {
  regulator_fun = approxfun(time, regulator_expression,method = "linear", rule = 2);
  return(function(t, state, parms) {
      dY = parms["w"] * regulator_fun(t) - parms["decay"] * state["y"];
      list(dX);
  })
}

simulate_linear <- function(num_cells, num_time, w = runif(1, 0.1, 3), decay = runif(1, 0.1, 1)) {
  regulator_expression = array(0, c(num_cells, num_time));
  target_expression = array(0, c(num_cells, num_time));
  ode_params = c(w = w, decay = decay);
  for(cell in 1:num_cells) {
    regulator_expression[cell,] = generate_random_profile(0:(num_time - 1),scale = 3, length = 2)
    ode_state = c(y = runif(1,0,1))
    ode_func = linear_regulation_ode_func(0:(num_time - 1), regulator_expression[cell,])
    target_expression[cell,] = ode(ode_state, times = 0:(num_time - 1), ode_func, method = "ode45")
  }
  return(list(observed =
                list(
                  num_cells = num_cells,
                  num_samples = num_time,
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
