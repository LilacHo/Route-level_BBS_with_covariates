// Route-varying zero-inflated negative-binomial test model.
//
// The earlier ZINB_test model used one global structural-zero probability for
// the whole species. The PPC results showed that was too blunt: excess zeros
// likely vary by route. This model keeps the same abundance process as the
// improved NB test model, but adds a route-level random effect on the
// zero-inflation logit:
//
//   logit(zi_route[s]) = zi_intercept + sd_zi_route * zi_route_raw[s]

functions {
  real icar_normal_lpdf(vector bb, int ns, array[] int n1, array[] int n2) {
    return -0.5 * dot_self(bb[n1] - bb[n2])
      + normal_lpdf(sum(bb) | 0, 0.001 * ns);
  }

  real zinb2_log_lpmf(int y, real eta, real phi, real zi) {
    if (y == 0) {
      return log_mix(zi, 0, neg_binomial_2_log_lpmf(0 | eta, phi));
    }
    return log1m(zi) + neg_binomial_2_log_lpmf(y | eta, phi);
  }
}

data {
  int<lower=1> nroutes;
  int<lower=1> ncounts;
  int<lower=1> nyears;
  int<lower=1> nobservers;

  array[ncounts] int<lower=0> count;
  array[ncounts] int<lower=1> year;
  array[ncounts] int<lower=1> route;
  array[ncounts] int<lower=0, upper=1> firstyr;
  array[ncounts] int<lower=1> observer;

  array[nroutes] real route_habitat;
  array[nroutes] real route_habitat_slope;

  int<lower=1> fixedyear;

  int<lower=1> N_edges;
  array[N_edges] int<lower=1, upper=nroutes> node1;
  array[N_edges] int<lower=1, upper=nroutes> node2;

  int<lower=0, upper=1> fit_spatial;
}

parameters {
  vector[nroutes] beta_raw;
  real BETA;
  vector[nroutes] rho_beta_raw_hab;
  real rho_BETA_hab;

  vector[nroutes] alpha_raw;
  real ALPHA;
  vector[nroutes] rho_alpha_raw_hab;
  real rho_ALPHA_hab;

  real eta;
  vector[nobservers] obs_raw;

  real zi_intercept;
  vector[nroutes] zi_route_raw;
  real<lower=0> sd_zi_route;

  real<lower=0> sdnoise;
  real<lower=0> sdobs;
  real<lower=0> sdbeta;
  real<lower=0> sdrho_beta_hab;
  real<lower=0> sdalpha;
  real<lower=0> sdrho_alpha_hab;
}

transformed parameters {
  vector[nroutes] beta;
  vector[nroutes] beta_resid;
  vector[nroutes] beta_hab;
  vector[nroutes] alpha;
  vector[nroutes] alpha_hab;
  vector[nroutes] alpha_resid;
  vector[nobservers] obs;
  vector<lower=0, upper=1>[nroutes] zi_route;
  real phi;
  vector[ncounts] E;

  beta_resid = sdbeta * beta_raw + BETA;
  alpha_resid = sdalpha * alpha_raw + ALPHA;

  for (s in 1:nroutes) {
    beta_hab[s] = (sdrho_beta_hab * rho_beta_raw_hab[s] + rho_BETA_hab)
      * route_habitat_slope[s];
    alpha_hab[s] = (sdrho_alpha_hab * rho_alpha_raw_hab[s] + rho_ALPHA_hab)
      * route_habitat[s];
    zi_route[s] = inv_logit(zi_intercept + sd_zi_route * zi_route_raw[s]);
  }

  beta = beta_resid + beta_hab;
  alpha = alpha_resid + alpha_hab;
  obs = sdobs * obs_raw;
  phi = 1 / sqrt(sdnoise);

  for (i in 1:ncounts) {
    E[i] = alpha[route[i]]
      + beta[route[i]] * (year[i] - fixedyear)
      + obs[observer[i]]
      + eta * firstyr[i];
  }
}

model {
  rho_beta_raw_hab ~ normal(0, 1);
  sum(rho_beta_raw_hab) ~ normal(0, 0.001 * nroutes);
  rho_alpha_raw_hab ~ normal(0, 1);
  sum(rho_alpha_raw_hab) ~ normal(0, 0.001 * nroutes);

  beta_raw ~ icar_normal(nroutes, node1, node2);
  if (fit_spatial) {
    alpha_raw ~ icar_normal(nroutes, node1, node2);
  } else {
    alpha_raw ~ normal(0, 1);
    sum(alpha_raw) ~ normal(0, 0.001 * nroutes);
  }

  obs_raw ~ normal(0, 1);
  sum(obs_raw) ~ normal(0, 0.001 * nobservers);

  zi_route_raw ~ normal(0, 1);
  sum(zi_route_raw) ~ normal(0, 0.001 * nroutes);
  zi_intercept ~ normal(-2, 1.5);
  sd_zi_route ~ normal(0, 1);

  sdnoise ~ normal(0, 0.5);
  sdobs ~ normal(0, 0.3);
  sdalpha ~ student_t(3, 0, 1);
  sdrho_alpha_hab ~ student_t(3, 0, 1);
  sdbeta ~ normal(0, 0.1);
  sdrho_beta_hab ~ normal(0, 0.1);

  BETA ~ normal(0, 0.1);
  rho_BETA_hab ~ normal(0, 0.1);
  ALPHA ~ normal(0, 1);
  rho_ALPHA_hab ~ normal(0, 0.5);
  eta ~ normal(0, 1);

  for (i in 1:ncounts) {
    target += zinb2_log_lpmf(count[i] | E[i], phi, zi_route[route[i]]);
  }
}

generated quantities {
  array[ncounts] real log_lik;
  real CH;
  real CH_no_habitat;
  real CH_dif;
  real T;
  real T_no_habitat;
  real T_dif;
  // Per-route change in expected observed count, i.e. the marginal mean
  // (1 - zi_route[s]) * NB_mean. This keeps CH_no_habitat on the same estimand
  // as the NB and ZINB_route_mu models. Here zi_route is time-invariant, so the
  // (1 - zi) factor cancels in each per-route ratio; it is kept explicit for
  // consistency and to document the marginal-mean definition.
  vector[nroutes] CH_route;
  vector[nroutes] CH_no_habitat_route;
  vector[nroutes] CH_dif_route;

  {
    real first_full = 0;
    real last_full = 0;
    real first_no_habitat = 0;
    real last_no_habitat = 0;

    for (i in 1:ncounts) {
      log_lik[i] = zinb2_log_lpmf(count[i] | E[i], phi, zi_route[route[i]]);
    }

    for (s in 1:nroutes) {
      real keep = 1 - zi_route[s];
      real m_first_full = keep * exp(alpha[s] + beta[s] * (1 - fixedyear));
      real m_last_full = keep * exp(alpha[s] + beta[s] * (nyears - fixedyear));
      real m_first_no_habitat = keep * exp(alpha[s] + beta_resid[s] * (1 - fixedyear));
      real m_last_no_habitat = keep * exp(alpha[s] + beta_resid[s] * (nyears - fixedyear));

      first_full += m_first_full;
      last_full += m_last_full;
      first_no_habitat += m_first_no_habitat;
      last_no_habitat += m_last_no_habitat;

      CH_route[s] = 100 * (m_last_full / m_first_full - 1);
      CH_no_habitat_route[s] = 100 * (m_last_no_habitat / m_first_no_habitat - 1);
      CH_dif_route[s] = CH_route[s] - CH_no_habitat_route[s];
    }

    first_full /= nroutes;
    last_full /= nroutes;
    first_no_habitat /= nroutes;
    last_no_habitat /= nroutes;

    CH = 100 * (last_full / first_full - 1);
    CH_no_habitat = 100 * (last_no_habitat / first_no_habitat - 1);
    CH_dif = CH - CH_no_habitat;
    T = 100 * ((last_full / first_full)^(1.0 / nyears) - 1);
    T_no_habitat = 100 * ((last_no_habitat / first_no_habitat)^(1.0 / nyears) - 1);
    T_dif = T - T_no_habitat;
  }
}
