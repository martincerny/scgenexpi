functions {
  int count_num_non_zero(int num_measurements, int[] counts) {
    int num = 0;
    for(i in 1:num_measurements) {
      if(counts[i] > 0) {
        num = num  + 1;
      }
    }
    return(num);
  }

  real non_zero_log_lik(int count, int category, int[] num_non_zeroes, real[] mean_count, real[] dispersion_inv, int ignore_dispersion_threshold)
    if(num_non_zeroes[category] <= ignore_dispersion_threshold) {
      //Ignoring dispersion, modelling as Poisson
      real poisson_prob_0 = poisson_lpmf(0 | mean_count[category]);
      return
        poisson_lpmf(count | mean_count[category])
        - log1m_exp(poisson_prob_0);
    } else {
      real neg_binom_prob_0 = neg_binomial_2_lpmf(0 | mean_count[category], inv(dispersion_inv[category]));
      return
        neg_binomial_2_lpmf(count | mean_count[category], inv(dispersion_inv[category]))
        - log1m_exp(neg_binom_prob_0);
    }

}

data {
  int<lower=0> num_measurements;
  int<lower=0> num_categories;

  int<lower=0> counts[num_measurements];
  int<lower=1, upper = num_categories> categories[num_measurements];

  real prior_zero_alpha;
  real prior_zero_beta;
  real prior_mean_mean;
  real prior_mean_sd;
  //real prior_dispersion_mean;
  //real prior_dispersion_sd;
  real prior_dispersion_rate;
  int ignore_dispersion_threshold;
}

transformed data {
  int<lower=0> num_non_zero = count_num_non_zero(num_measurements, counts);
  int<lower=1> non_zero_counts[num_non_zero];
  int<lower=1, upper = num_categories> non_zero_categories[num_non_zero];
  int<lower=0> num_zeroes[num_categories];
  int<lower=0> num_non_zeroes[num_categories];

  for(i in 1:num_categories) {
    num_zeroes[i] = 0;
    num_non_zeroes[i] = 0;
  }

  {
    int next_id = 1;
    for(i in 1:num_measurements) {
      if(counts[i] == 0){
        num_zeroes[categories[i]] = num_zeroes[categories[i]] + 1;
      } else {
        num_non_zeroes[categories[i]] = num_non_zeroes[categories[i]] + 1;
        non_zero_counts[next_id] = counts[i];
        non_zero_categories[next_id] = categories[i];
        next_id = next_id + 1;
      }
    }
  }
}

parameters {
  real<lower=0, upper=1> zero_prob_ext[num_categories];
  real<lower=0> mean_count[num_categories];
  real<lower=0> dispersion_inv[num_categories];
}

model {
  for(i in 1:num_non_zero) {
    int category = non_zero_categories[i];
    target += non_zero_log_lik(non_zero_counts[i], category, num_non_zeroes, mean_count, dispersion_inv, ignore_dispersion_threshold);
  }

  for(category in 1:num_categories) {
    num_zeroes[category] ~ binomial(num_zeroes[category] + num_non_zeroes[category], zero_prob_ext[category]);
  }

  zero_prob_ext ~ beta(prior_zero_alpha, prior_zero_beta);
  mean_count ~ lognormal(prior_mean_mean, prior_mean_sd);
//  dispersion ~ lognormal(prior_dispersion_mean, prior_dispersion_sd);
  dispersion_inv ~ exponential(prior_dispersion_rate);

}



generated quantities {
  real zero_prob[num_categories];
  vector[num_measurements] log_lik;

  for(category in 1:num_categories) {
    real prob_0;
    if(num_non_zeroes[category] <= ignore_dispersion_threshold) {
      prob_0 = exp(poisson_lpmf(0 | mean_count[category]));
    } else {
      prob_0 = exp(neg_binomial_2_lpmf(0 | mean_count[category], dispersion_inv[category]));
    }

    zero_prob[category] = (zero_prob_ext[category] - prob_0) / (1 - prob_0);
  }

  for(i in 1:num_measurements) {
    int category = categories[i];
    if(counts[i] > 0) {
      log_lik[i] = non_zero_log_lik(counts[i], category, num_non_zeroes, mean_count, dispersion_inv, ignore_dispersion_threshold);
    }
    else {
      log_lik[i] = bernoulli_lpmf(1 | zero_prob_ext[category]);
    }
  }
}

