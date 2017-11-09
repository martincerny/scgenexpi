density_pcr <- function(x, initial_no, prob_duplicate, n_cycles) {
  if(initial_no <= 0) {
    stop("Initial no. has to be positive")
  }
  if(any(x <= 0)) {
    stop("x has to be positive")
  }

  #computed_densities = array(-1, c(1:x,1:n_cycles))

  if(n_cycles == 0) {
    if(x == initial_no) {
      return(1);
    } else {
      return(0);
    }
  } else {
    if(x < initial_no) {
      return(0)
    } else if(x > initial_no * 2 ^ n_cycles) {
      return(0)
    }
    max_duplicated = floor(x / 2)
    prob = 0
    for(duplicated in 0:max_duplicated) {
      prob = prob + dbinom(duplicated, x - duplicated, prob_duplicate) * density_pcr(x - duplicated, initial_no, prob_duplicate, n_cycles - 1)
    }
    return(prob)
  }
}

rbinom_corr <- function(n, size, prob) {
  return(if_else(size < .Machine$integer.max, as.double(rbinom(n, size, prob)), floor(rnorm(n, size*prob, size*prob*(1-prob)))))
}

simulate_pcr <- function(init, num_samples, num_cycles, amplify_prob ) {
  pcr_state = tibble(init = init) %>%
    mutate(amplified = init)
  result = NULL
  for(s in 1:num_samples) {
    for(i in 1:num_cycles) {
      pcr_state %<>% mutate(amplified = amplified + rbinom_corr(length(amplified), amplified, amplify_prob))
    }
    if(is.null(result)) {
      result = pcr_state
    } else {
      result = rbind(result, pcr_state)
    }
    pcr_state %<>% mutate(amplified = init)
  }
  return(result)
}
