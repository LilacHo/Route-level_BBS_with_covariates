// Habitat-covariate route-level slope model (NB2 count likelihood)
// Merges two prior variants:
//   - Numerically robust parameterization from slope_habitat_route_NB_New.stan:
//     beta_raw, alpha_raw, rho_beta_raw_hab, rho_alpha_raw_hab, and obs_raw
//     use Stan's built-in sum_to_zero_vector type (hard, exact zero-sum
//     constraint) instead of a soft/penalized sum-to-zero constraint. The
//     custom icar_normal_lpdf() helper is dropped; the ICAR pairwise-difference
//     penalty is added directly via target += -0.5 * dot_self(...).
//   - Diagnostic / route-level outputs from slope_habitat_route_NB_test.stan:
//     pointwise log_lik for LOO/PPC, and per-route CH_route /
//     CH_no_habitat_route / CH_dif_route computed directly from first- and
//     last-year expected means (rather than the full nroutes x nyears
//     trajectory arrays), so route-level habitat-excluded trend/change is
//     available without post-processing a saved nsmooth array. Also keeps
//     _test's tighter priors on sdalpha, sdrho_alpha_hab (student_t(3,0,1))
//     and rho_ALPHA_hab (normal(0,0.5)).
//
// Decomposes route intercepts/slopes into a residual (habitat-excluded)
// component and a habitat-driven component:
//   alpha = alpha_resid + alpha_hab,  beta = beta_resid + beta_hab
// beta_resid (i.e. beta with the habitat-slope contribution removed) is the
// route-level habitat-excluded trend; CH_no_habitat_route / T_no_habitat
// (derived from CH_no_habitat_route) give that trend on the percent-change
// scale per route.

data {
  int<lower=1> nroutes;
  int<lower=1> ncounts;
  int<lower=1> nyears;
  int<lower=1> nobservers;

  array [ncounts] int<lower=0> count;   // count observations
  array [ncounts] int<lower=1> year; // year index
  array [ncounts] int<lower=1> route; // route index
  array [ncounts] int<lower=0, upper=1> firstyr; // first year index =1 if observer's first year on route, 0 otherwise
  array [ncounts] int<lower=1> observer;   // observer indicators

  // mean annual habitat suitability on route (centered and scaled)
  array [nroutes] real route_habitat;
  // mean rate of change in annual habitat suitability on route (centered)
  array [nroutes] real route_habitat_slope;

  int<lower=1> fixedyear; // centering value for years

 // spatial neighbourhood information
  int<lower=1> N_edges;
  array [N_edges] int<lower=1, upper=nroutes> node1;  // node1[i] adjacent to node2[i]
  array [N_edges] int<lower=1, upper=nroutes> node2;  // and node1[i] < node2[i]

  int<lower=0, upper=1> fit_spatial;
  // conditional:
  // if 1 then use spatial component to model intercept residual
  // if 0 then use simple exchangeable random effect

}

parameters {

  sum_to_zero_vector[nroutes] beta_raw;// non-centered residual trend, spatial (iCAR)
  real BETA; // hyperparameter mean residual trend
  sum_to_zero_vector[nroutes] rho_beta_raw_hab;// non-centered habitat-based trend
  real rho_BETA_hab; // hyperparameter mean habitat-based trend

  sum_to_zero_vector[nroutes] alpha_raw;// non-centered residual intercept
  real ALPHA; // hyperparameter mean residual intercept
  sum_to_zero_vector[nroutes] rho_alpha_raw_hab;// non-centered habitat-based intercept
  real rho_ALPHA_hab; // hyperparameter mean habitat-based intercept

  real eta; //first-year intercept

  sum_to_zero_vector[nobservers] obs_raw; //observer effects

  real<lower=0> sdnoise;    // inverse of sd of over-dispersion
  real<lower=0> sdobs;    // sd of observer effects
  real<lower=0> sdbeta;    // sd of residual slopes
  real<lower=0> sdrho_beta_hab;    // sd of habitat-change effect on slopes
  real<lower=0> sdalpha;    // sd of residual intercepts
  real<lower=0> sdrho_alpha_hab;    // sd of habitat effect on intercepts

}

transformed parameters{

   vector[nroutes] beta; // full slope
   vector[nroutes] beta_resid; //residual component of slope (route-level habitat-excluded trend)
   vector[nroutes] beta_hab; //habitat component of slope

  vector[nroutes] alpha; // full intercept
  vector[nroutes] alpha_hab; //habitat component intercepts
  vector[nroutes] alpha_resid; // residual component intercepts
  vector[nobservers] obs; // observer effects
  real phi;//dispersion of negative binomial
  vector[ncounts] E;           // predicted log-scale counts (lambda)


// covariate effect on intercepts and slopes

   beta_resid = (sdbeta*beta_raw) + BETA;
   alpha_resid = (sdalpha*alpha_raw) + ALPHA;

   for(s in 1:nroutes){
   beta_hab[s] = ((sdrho_beta_hab*rho_beta_raw_hab[s]) + rho_BETA_hab) * (route_habitat_slope[s]);
   alpha_hab[s] = ((sdrho_alpha_hab*rho_alpha_raw_hab[s]) + rho_ALPHA_hab) * (route_habitat[s]);
   }

   alpha =  alpha_resid + alpha_hab;
   beta = beta_resid  + beta_hab;


   obs = sdobs*obs_raw;

  phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations

  for(i in 1:ncounts){
    E[i] =  alpha[route[i]] + beta[route[i]] * (year[i]-fixedyear)  + obs[observer[i]] + eta*firstyr[i];
  }


}

model {

  // habitat-change effects on slope and intercept (exchangeable, hard sum-to-zero)
  rho_beta_raw_hab ~  normal(0,1);
  rho_alpha_raw_hab ~ normal(0,1);

    // spatially varying residual trend and intercepts
   target += -0.5 * dot_self(beta_raw[node1] - beta_raw[node2]); // ICAR prior
  if(fit_spatial){
   target += -0.5 * dot_self(alpha_raw[node1] - alpha_raw[node2]); // ICAR prior
  }else{
  alpha_raw ~ std_normal(); // simple exchangeable random effect
  }


  sdnoise ~ normal(0,0.5); //prior on scale of extra Poisson log-normal variance

  sdobs ~ normal(0,0.3); //prior on sd of gam hyperparameters

  obs_raw ~ std_normal();//observer effects

  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  rho_BETA_hab ~ normal(0,0.1);// prior on habitat-change effect on slope
  ALPHA ~ normal(0,1);// prior on fixed effect mean intercept
  rho_ALPHA_hab ~ normal(0,0.5);// prior on habitat effect on intercept
  eta ~ normal(0,1);// prior on first-year observer effect

  // Tighter, heavier-tailed priors on the intercept SDs (see
  // slope_habitat_route_NB.stan for the mixing-diagnosis rationale):
  // student_t(3,0,1) concentrates mass at plausible values but keeps heavy
  // tails so genuinely large between-route variance is still reachable.
  sdalpha ~ student_t(3, 0, 1); //prior on sd of intercept
  sdrho_alpha_hab ~ student_t(3, 0, 1); //prior on sd of habitat effect on intercept

  sdbeta ~ normal(0,0.1);// prior on sd of slope spatial variation w mean = 0.04 and 99% < 0.13
  sdrho_beta_hab ~ normal(0,0.1);// prior on sd of habitat-change effect on slope


  // full count likelihood
  count ~ neg_binomial_2_log(E,phi); //vectorized count likelihood with log-transformation

}

generated quantities {
  array[ncounts] real log_lik; // pointwise log-likelihood for LOO / PSIS
  real CH;
  real CH_no_habitat;
  real CH_dif;
  real T;
  real T_no_habitat;
  real T_dif;
  // Per-route change in expected observed count, with and without the
  // habitat-driven component of the slope. CH_no_habitat_route / T_no_habitat
  // (below) is the route-level habitat-excluded trend.
  vector[nroutes] CH_route;
  vector[nroutes] CH_no_habitat_route;
  vector[nroutes] CH_dif_route;
  vector[nroutes] T_route;
  vector[nroutes] T_no_habitat_route;
  vector[nroutes] T_dif_route;

  {
    real first_full = 0;
    real last_full = 0;
    real first_no_habitat = 0;
    real last_no_habitat = 0;

    for (i in 1:ncounts) {
      log_lik[i] = neg_binomial_2_log_lpmf(count[i] | E[i], phi);
    }

    for (s in 1:nroutes) {
      real m_first_full = exp(alpha[s] + beta[s] * (1 - fixedyear));
      real m_last_full = exp(alpha[s] + beta[s] * (nyears - fixedyear));
      real m_first_no_habitat = exp(alpha[s] + beta_resid[s] * (1 - fixedyear));
      real m_last_no_habitat = exp(alpha[s] + beta_resid[s] * (nyears - fixedyear));

      first_full += m_first_full;
      last_full += m_last_full;
      first_no_habitat += m_first_no_habitat;
      last_no_habitat += m_last_no_habitat;

      CH_route[s] = 100 * (m_last_full / m_first_full - 1);
      CH_no_habitat_route[s] = 100 * (m_last_no_habitat / m_first_no_habitat - 1);
      CH_dif_route[s] = CH_route[s] - CH_no_habitat_route[s];

      T_route[s] = 100 * ((m_last_full / m_first_full)^(1.0 / nyears) - 1);
      T_no_habitat_route[s] = 100 * ((m_last_no_habitat / m_first_no_habitat)^(1.0 / nyears) - 1);
      T_dif_route[s] = T_route[s] - T_no_habitat_route[s];
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
