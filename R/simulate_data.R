require(deSolve)

gp_covariance <- function(distance, gp_scale, gp_length, periodic = FALSE, period = 1) {
  if(periodic) {
    return((gp_scale^2) * exp(-2*(sin(3.141592 * abs(distance) / period)^2) / (gp_length ^ 2)))
  } else {
    return((gp_scale^2) * exp( (-0.5 / (gp_length ^ 2)) * (distance ^ 2) ));
  }

}

generate_random_profile <- function(time, scale, length, mean_func = 0, periodic = FALSE, period = 1, positive_transform = TRUE) {
  # Construct the squared exponential covariance matrix
  cov_matrix = array(0,c(length(time), length(time)));
  maxTime = max(time)
  minTime = min(time)
  for(i in 1:length(time))  {
    cov_matrix[i,i] = scale^2 + 0.00000001; #Adding jitter to ensure positive semi-definiteness
    if (i < length(time)) {
      for(j in (i+1):length(time)) {
        distance = (time[i] - time[j])
        covariance = gp_covariance(distance, scale, length, periodic, period)
        cov_matrix[i,j] = covariance
        cov_matrix[j,i] = covariance
      }
    }
  }
  # Draw from the multinormal distribution using the cholesky decomposition of the cov. matrix
  chol_cov = t(chol(cov_matrix));
  raw_profile = chol_cov %*% rnorm(length(time)) + mean_func;
  # Transform to strictly positive values
  if(positive_transform) {
    positive_profile = log1p(exp(raw_profile))
    return(t(positive_profile))
  } else {
    return(t(raw_profile))
  }
}

plot_random_profiles <- function(n, time, scale, length, true_time = time, true_profile = NULL, mean_func = 0, periodic = FALSE, period = 1, positive_transform = TRUE) {
  profiles = array(0, c(n, length(time)));
  for(i in 1:n) {
    profiles[i,] = generate_random_profile(time, scale, length, mean_func = mean_func, periodic = periodic, period = period, positive_transform = positive_transform);
  }


  result = ggmatplot(time, t(profiles))
  if(!is.null(true_profile)) {
    result = result + geom_point(data = data.frame(x = true_time, y = true_profile), aes(x=x, y=y))
  }
  return(result)
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

simulate_linear_circular_2 <- function(num_cells, num_time, scale, gp_length, w = runif(1, 0.1, 3), decay = runif(1, 0.1, 1), min_variance = 1) {

  time_points = runif(num_cells,0,1);
  time = time_points * (num_time - 1);

  repeat {
    regulator_profile = array(generate_random_profile(1:num_time, scale, gp_length, periodic = TRUE, period = num_time - 1), num_time);
    if(var(regulator_profile) > min_variance) {
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
  regulator_expression = regulator_expression_true #rnorm(num_cells, regulator_expression_true, expression_sigma);
  regulator_expression[regulator_expression < 0] = 0;

  return(list(observed =
                list(
                  num_cells = num_cells,
                  num_time = num_time,
                  bandwidth = 0.5,
                  regulator_expression = regulator_expression,
                  target_expression = target_expression,
                  gp_sigma = scale,
                  gp_length = gp_length,
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

generate_random_2d_field <- function(grid_x, grid_y, gp_scale, gp_length) {
  total_points = length(grid_x) * length(grid_y)
  cov_matrix = array(NA, c(total_points,total_points))
  for(x1 in 1:length(grid_x)) {
    for(y1 in 1:length(grid_y)) {
      coord1 = (x1 - 1) * length(grid_y) + y1
      cov_matrix[coord1, coord1] = gp_scale ^ 2 + 1e-6
      for(x2 in 1:length(grid_x)) {
        for(y2 in 1:length(grid_y)) {
          coord2 = (x2 - 1) * length(grid_y) + y2
          if(coord2 > coord1) {
            distance = sqrt((x1-x2) ^ 2 + (y1 - y2) ^ 2)
            covariance = gp_covariance(distance, gp_scale, gp_length)
            cov_matrix[coord1, coord2] = covariance
            cov_matrix[coord2, coord1] = covariance
          }
        }
      }
    }
  }
  chol_cov = t(chol(cov_matrix));
  flat_values = chol_cov %*% rnorm(total_points);

  result = array(NA, c(length(grid_x),length(grid_y)))
  for(x in 1:length(grid_x)) {
    for(y in 1:length(grid_y)) {
      coord = (x - 1) * length(grid_y) + y
      result[x,y] = flat_values[coord]
    }
  }
  return(result)
}

kernel_estimate_2d <- function(xs, ys, values, target_x, target_y, bandwidth) {
  distances = sqrt((xs - target_x) ^ 2 + (ys - target_y) ^ 2)
  weights = exp(-0.5 * (distances / bandwidth) ^ 2)
  return(sum(weights * values) / sum(weights))
}

evolve_2d_field <- function(field, bandwidth, regulator_expression, target_expression, step) {
  next_regulator = array(NA, length(regulator_expression))
  next_target = array(NA, length(target_expression))
  for(i in 1:length(regulator_expression)) {
    regulator_change = kernel_estimate_2d(field$xs,field$ys, field$values, regulator_expression[i], target_expression[i], bandwidth)
    next_regulator[i] = max(regulator_expression[i] + regulator_change, 0)

    line_a = (next_regulator[i] - regulator_expression[i]) / step
    line_b = regulator_expression[i]
    integral_start = (line_b / field$decay) - (line_a / (field$decay ^ 2))
    integral_end = exp(field$decay * step) * (  ( (line_a * step + line_b) / field$decay ) - (line_a / (field$decay ^ 2)))
    next_target[i] = exp(-field$decay * step) * (target_expression[i] + field$w * (integral_end - integral_start))
  }
  return(list(regulator_expression = next_regulator, target_expression = next_target))
}

rotation_field <- function() {
  coords = c(0.5,1.5,2.5)
  list(xs = rep(coords, 3),
       ys = c(rep(coords[1],3),rep(coords[2],3), rep(coords[3],3)),
       values = c(1,    0.5, -0.1,
                  -0.1, -0.5, -1,
                  -0.3,   -0.8, -1.5),
       w = 2,
       decay = 0.8
       )
}

visualize_field <- function(field, initial_pop, bandwidth, step, num_steps) {
  index = 1:length(initial_pop$regulator_expression)
  values = data.frame(reg = initial_pop$regulator_expression, target = initial_pop$target_expression,index = index)
  values$step = 0
  pop = initial_pop;
  for(i in 1:num_steps) {
    pop = evolve_2d_field(field, bandwidth, pop$regulator_expression,pop$target_expression, step)
    new_step = data.frame(reg = pop$regulator_expression, target = pop$target_expression,index = index)
    new_step$step = i
    values = rbind(values, new_step)
  }
  ggplot(values, aes(x=reg, y = target, color = as.factor(index))) + scale_colour_discrete(labels=NULL, guide="none") + geom_path()
}

field_map <- function(field, bandwidth) {
  grid_size = 0.1
  grid_x = seq(min(field$xs) - 0.5, max(field$xs) + 0.5, by = grid_size)
  grid_y = seq(min(field$ys) - 0.5, max(field$ys) + 0.5, by = grid_size)
  regulator = crossing(reg = grid_x, target = grid_y) %>%
    rowwise() %>% mutate(value = kernel_estimate_2d(field$xs,field$ys, field$values, reg, target, bandwidth))



  ggplot(regulator, aes(x = reg, y = target, fill = value)) + geom_tile(size = grid_size) + scale_fill_gradient2() +
  geom_abline(slope = field$w/field$decay , intercept = 0)
}
