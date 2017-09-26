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

    if(num_non_zeroes[category] <= ignore_dispersion_threshold) {
      //Ignoring dispersion, modelling as Poisson
      real poisson_prob_0 = poisson_lpmf(0 | mean_count[category]);
      target +=
        poisson_lpmf(non_zero_counts[i] | mean_count[category])
        - log1m_exp(poisson_prob_0);
    } else {
      real neg_binom_prob_0 = neg_binomial_2_lpmf(0 | mean_count[category], inv(dispersion_inv[category]));
      target +=
        neg_binomial_2_lpmf(non_zero_counts[i] | mean_count[category], inv(dispersion_inv[category]))
        - log1m_exp(neg_binom_prob_0);
    }
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

  for(category in 1:num_categories) {
    real prob_0;
    if(num_non_zeroes[category] <= ignore_dispersion_threshold) {
      prob_0 = exp(poisson_lpmf(0 | mean_count[category]));
    } else {
      prob_0 = exp(neg_binomial_2_lpmf(0 | mean_count[category], dispersion_inv[category]));
    }

    zero_prob[category] = (zero_prob_ext[category] - prob_0) / (1 - prob_0);
  }
}

