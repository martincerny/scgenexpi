empirical_smooth_cdf_logit <- function(samples, slope, x) {
  result = rep_len(0, length(x));
  for(s in samples) {
    entry = x - s
    logistic = 1 / (1 + exp(-slope * entry));
    result = result + logistic;
  }
  return(result / length(samples));
}

empirical_smooth_cdf_normal <- function(samples, bandwidth, x) {
  result = rep_len(0, length(x));
  for(s in samples) {
    result = result + pnorm(x, s, bandwidth);
  }
  return(result / length(samples));
}


empirical_density_logit <- function(samples, slope, x) {
  result = rep_len(0, length(x));
  for(s in samples) {
    entry = x - s
    logistic = 1 / (1 + exp(-slope * entry));
    add = logistic * (1 - logistic);
    result = result + add;
  }
  return((result * slope) / length(samples));
}

empirical_density_normal <- function(samples, bandwidth, x) {
  result = rep_len(0, length(x));
  for(s in samples) {
    result = result + dnorm(x, mean = s, sd = bandwidth);
  }
  return(result / length(samples));

}

kernel_estimate <- function(x_in, y_in, bandwidth, x_out) {
  result = rep_len(0, length(x_out));
  for(i in 1:length(x_out)) {
    weights = exp(-0.5 * ((x_in - x_out[i]) / bandwidth) ^ 2)
    result[i] = sum(weights * y_in) / sum(weights);
  }
  return(result);
}
