// Endemic-epidemic space-time model of daily cases by ABS: stan/daily.stan (Slater et al. 2025 model 1, daily, fixed
// one-week lag) with extensions, each switched by data so that fits differing by one change can be
// compared. The final model is MODEL in helpers.R. Switches:
// - time_effect: their per-contact probability of infection p (Section 3.1, alpha = c p) made time-varying. p_i,t
//   multiplies the whole epidemic part (AR + NE + MOB) of the destination ABS i; their p_ji,t (Section 3.2) also
//   depends on the origin, ours does not. The coefficients are the values at log p = 0: the geometric mean over weeks,
//   without vaccination.
//   0: p = 1 (their model); 1: log p_t = gamma_t, a weekly RW1 common to all ABS; 2: gamma_t + delta_i,t, an ABS x week
//   interaction whose weekly increments follow a scaled ICAR (Knorr-Held type IV); 3: gamma_t + delta_i,t with delta
//   MacNab's (2022, eqs. 25-26) SVC: delta_i,t = sigma s_i eta_i,t, eta_i a temporal proper CAR with area-specific
//   dependence c_i, log s ~ LCAR (the prior on c_i, the scale and sum-to-zero of log s and the centring of delta are
//   ours). Weekly values interpolated to days.
// - vaccination: p times exp(-tau u_i,t), u = share with a first dose (their Appendix D, which uses u one week earlier;
//   14 days here, ours)
// - testing: 1: the whole mean times (x_t / max x)^beta, x = tests in the ABS's province (their eq. 5); 2: the
//   endemic term times pi_t = (x_t / max x)^beta and the epidemic terms times pi_t / pi_t-s, as the lagged cases are
//   already reported at pi_t-s (pi_t-s: geometric mean over the lag window, weights of the one-week lag)
// - pooled: 0: their independent half-normals; 1: log coefficients ~ normal(mu_k, sd_k) across ABS; 2: as 1 for the
//   endemic coefficient, and for the epidemic coefficients an adaptive Leroux CAR by linear transformation
//   (MacNab 2022, eq. 19; LCAR(f) of her Table 3, from MacNab 2018): log alpha_k,i = mu_k + sd_k s_i z_k,i,
//   z_k ~ LCAR(lambda_k) with unit scale; log s ~ LCAR as her eq. 25 (she gives no prior on s for LCAR(f)); s_i shared
//   by the epidemic coefficients and sum-to-zero on z_k and log s (ours)
// - covariates: static ABS covariates (standardised) in the means of the pooled log coefficients (needs pooled > 0):
//   log alpha_k,i = mu_k + sum_c beta_k,c x_i,c + ..., covariate c on coefficient row cov_row[c] (hhh4 fixed effects
//   with ri() random intercepts, Paul and Held 2011); beta ~ normal(0, 0.5) (ours)
// - negbin: negative binomial in place of the Poisson; phi_area: 0 one size phi; 1 a size phi_i per ABS, log phi_i
//   ~ normal across ABS (their NB had area-specific sizes psi_i)
// - rw_order: 1 or 2, the order of the random walk of gamma; gamma_prior: half-normal scale of its sd
// - npi: log p plus beta' r_t, r_t the share of days t-10..t-4 with each restriction in force (region-wide)
// - mob_closure: the mobility term alone also times exp(beta_mc c_t), c_t that share for the bar and restaurant closure
// - lag_dist: weights u_s of the case lags s = 1..n_lags. 0: fixed, applied in R (equal weights on days 7-13, their
//   weekly Y_t-1 in days); 1: a shifted negative binomial, u_s proportional to NB(s - lag_min | lag_mean, lag_phi) for
//   s >= lag_min, 0 below (Bracher and Held 2022; Douwes-Schultz et al. 2022). lag_min keeps the weights off the
//   first days, where an unshifted lag puts most weight on yesterday's count (reporting, not transmission)
// - mob_lags: 2 or 3 mobility lags (trips of the 7 days up to each lagged case day, the 7 before, and the 7 before
//   those); their castilla "threelags" model has three
// - n_orders: 1: neighbours of first order, weights 1/Ne(j) (theirs); 5: orders o = 1..5 weighted o^-d (power law,
//   Meyer and Held 2014), each order's pressure normalised by the origin's number of neighbours of that order (ours;
//   Meyer and Held normalise over all orders at once); d -> infinity gives theirs
// - within-ABS trips: left in or out of the mobility pressures in R (spacetime_data, no_within)
// - holiday: the reporting factor also times exp(b1 h_t + b2 a_t), h = public holiday, a = the day after one
// - day_effect: the reporting factor also times exp(e_t), a shock of day t shared by all ABS (lab backlogs, batches
//   of results), e_t = sigma_e (z_t - mean of z over its 7-day block), z ~ N(0, 1): sum-to-zero within each block
//   (the cross-validation weeks), so the weekly level stays with gamma and a held-out week has e from its prior (ours);
//   day_effect = 2: one e_g,t per group of ABS (day_group, e.g. province), one sigma_e
// - save_expected: also save the expected cases of each ABS and day (posterior predictive checks)
// - holdout: ABS-days left out of the likelihood (cross-validation); their cases still enter the lagged pressures of
//   later days, as observed history
// No spatial random effect in the endemic part. The lagged pressures are computed in R (spacetime_data), one column
// per lag (lag_dist = 1) or already weighted (one column). Priors not marked as theirs are ours.
functions {
	// Leroux CAR: precision (lambda (D - A) + (1 - lambda) I) / sigma^2, with the eigenvalues of D - A
	real lcar_lpdf(vector x, real lambda, real sigma, array[] int node1, array[] int node2, vector laplacian_eigen) {
		return 0.5 * sum(log(lambda * laplacian_eigen + 1 - lambda)) - rows(x) * log(sigma)
			- 0.5 / square(sigma) * (lambda * dot_self(x[node1] - x[node2]) + (1 - lambda) * dot_self(x));
	}
	// temporal proper CAR per area (one column per area): precision D_T - c_i W_T on the chain of weeks, with
	// week_eigen the eigenvalues of D_T^-1/2 W_T D_T^-1/2 (MacNab 2022, eq. 26)
	real tpcar_lpdf(matrix x, vector c, vector week_n_neighbours, vector week_eigen) {
		int n_weeks = rows(x);
		row_vector[cols(x)] q = week_n_neighbours' * square(x)
			- 2 * c' .* columns_dot_product(x[1:(n_weeks - 1)], x[2:n_weeks]);
		return 0.5 * sum(log1m(week_eigen * c')) - 0.5 * sum(q);
	}
	// p_i,t: weekly log values interpolated linearly between week midpoints, minus the vaccination term, exponentiated
	// npi_lp: the restriction term of each day, common to all ABS
	matrix transmission(vector gamma, matrix delta, array[] real tau, matrix vaccinated, array[] int lo, array[] int hi,
			vector frac, vector npi_lp) {
		matrix[rows(delta), size(lo)] lp;
		for (t in 1:size(lo))
			lp[, t] = (1 - frac[t]) * (gamma[lo[t]] + delta[, lo[t]]) + frac[t] * (gamma[hi[t]] + delta[, hi[t]]) + npi_lp[t];
		if (size(tau) > 0) lp -= tau[1] * vaccinated;
		return exp(lp);
	}
	// shifted negative binomial lag weights on s = 1..n_lags, zero below lag_min, summing to 1
	vector lag_weights(real lag_mean, real lag_phi, int n_lags, int lag_min) {
		vector[n_lags - lag_min + 1] lw;
		for (s in lag_min:n_lags) lw[s - lag_min + 1] = neg_binomial_2_lpmf(s - lag_min | lag_mean, lag_phi);
		return append_row(rep_vector(0, lag_min - 1), softmax(lw));
	}
	// an area x day pressure from its lag columns (column-major area x day rows) and the lag weights
	matrix lagged(matrix columns, vector lag_weight, int n_areas, int n_days) {
		return to_matrix(columns * lag_weight, n_areas, n_days);
	}
	// the neighbour pressure over adjacency orders, weighted by order_weight
	matrix neighbours(array[] matrix columns, vector lag_weight, vector order_weight, int n_areas, int n_days) {
		matrix[n_areas, n_days] m = order_weight[1] * lagged(columns[1], lag_weight, n_areas, n_days);
		for (o in 2:size(columns)) m += order_weight[o] * lagged(columns[o], lag_weight, n_areas, n_days);
		return m;
	}
	// expected cases by component, before the reporting factors: endemic, AR, NE, MOB (summed over the mobility lags);
	// mob_factor: the mobility-only factor of each day
	array[] matrix components(matrix coef, vector pop, matrix p, matrix ar, matrix ne, array[] matrix mob,
			vector mob_factor) {
		array[4] matrix[rows(ar), cols(ar)] x;
		matrix[rows(ar), cols(ar)] m = diag_pre_multiply(coef[4]', mob[1]);
		for (d in 2:size(mob)) m += diag_pre_multiply(coef[3 + d]', mob[d]);
		x[1] = rep_matrix(coef[1]' .* pop, cols(ar));
		x[2] = diag_pre_multiply(coef[2]', ar) .* p;
		x[3] = diag_pre_multiply(coef[3]' .* pop, ne) .* p;    // pop_i factor on the neighbour term, as their model 1
		x[4] = diag_post_multiply(m .* p, mob_factor);
		return x;
	}
}

data {
	int<lower=0> n_areas;
	int<lower=0> n_days;
	array[n_areas, n_days] int<lower=0> cases;                       // Y_i,t
	int<lower=1> n_lags;                                             // case lags s = 1..n_lags
	int<lower=0, upper=1> lag_dist;
	int<lower=1, upper=n_lags> lag_min;
	int<lower=1, upper=5> n_orders;
	int<lower=2, upper=3> mob_lags;
	// lagged pressures, rows = ABS-days (column-major area x day), columns = lags (lag_dist = 1) or one weighted column
	matrix<lower=0>[n_areas * n_days, lag_dist ? n_lags : 1] autoregressive_lags;      // Y_i,t-s
	array[n_orders] matrix<lower=0>[n_areas * n_days, lag_dist ? n_lags : 1] neighbour_lags;   // sum_j v_ji^(o) Y_j,t-s / N_j
	array[mob_lags] matrix<lower=0>[n_areas * n_days, lag_dist ? n_lags : 1] mobility_lags;    // sum_j w_ji^(d) Y_j,t-s / N_j
	vector<lower=0>[n_areas] pop;                                    // N_i
	array[n_days] int<lower=1, upper=7> weekday;                     // Monday = 1
	vector<lower=0>[3 + mob_lags] prior_scale;                       // their half-normal scales: endemic, AR, NE, MOB1, MOB2(, MOB3)
	int<lower=0, upper=3> time_effect;
	int<lower=0, upper=1> vaccination;
	int<lower=0, upper=2> testing;
	int<lower=0, upper=2> pooled;
	int<lower=0> n_cov;                                              // static covariates (0: none)
	matrix[n_areas, n_cov] covariate;                                // standardised over the ABS
	array[n_cov] int<lower=1, upper=3 + mob_lags> cov_row;           // the coefficient each one enters
	int<lower=0, upper=1> negbin;
	int<lower=0, upper=1> phi_area;
	int<lower=1> n_weeks;
	array[n_days] int<lower=1, upper=n_weeks> week_lo;               // week midpoints before and after each day
	array[n_days] int<lower=1, upper=n_weeks> week_hi;
	vector<lower=0, upper=1>[n_days] week_frac;                      // position between them
	vector<lower=1>[n_weeks] week_n_neighbours;                      // D_T of the chain of weeks
	vector[n_weeks] week_eigen;                                      // eigenvalues of D_T^-1/2 W_T D_T^-1/2
	int<lower=0> n_edges;                                            // queen adjacency, i < j
	array[n_edges] int<lower=1, upper=n_areas> node1;
	array[n_edges] int<lower=1, upper=n_areas> node2;
	real<lower=0> icar_scale;                                        // Riebler et al. (2016) scaling factor
	vector[n_areas] laplacian_eigen;                                 // eigenvalues of D - A
	matrix<lower=0, upper=1>[n_areas, n_days] vaccinated;            // u_i,t-14
	matrix<upper=0>[n_areas, n_days] log_tests;                      // log(x_t / max x), province of each ABS
	matrix[n_areas, n_days] log_tests_lag;                           // the same, averaged over the case lags of day t
	int<lower=1, upper=2> rw_order;
	real<lower=0> gamma_prior;
	int<lower=1> n_restrictions;
	matrix<lower=0, upper=1>[n_days, n_restrictions] restrictions;   // share of days t-10..t-4 with each restriction
	int<lower=0, upper=n_restrictions> n_npi;
	array[n_npi] int<lower=1, upper=n_restrictions> npi_index;       // the restrictions in log p
	int<lower=0, upper=1> mob_closure;
	vector<lower=0, upper=1>[n_days] closure;                        // the same share for the bar and restaurant closure
	int<lower=0, upper=1> holiday;
	vector<lower=0, upper=1>[n_days] holiday_day;                    // public holiday (national or Catalan)
	vector<lower=0, upper=1>[n_days] after_holiday;                  // the day after one, not itself a holiday
	int<lower=0, upper=2> day_effect;
	int<lower=1> n_day_groups;
	array[n_areas] int<lower=1, upper=n_day_groups> day_group;       // group of each ABS for day_effect = 2
	int<lower=1> n_blocks;
	array[n_days] int<lower=1, upper=n_blocks> day_block;            // 7-day block of each day (as the CV weeks)
	int<lower=0, upper=1> save_expected;
	array[n_areas, n_days] int<lower=0, upper=1> holdout;
}

transformed data {
	int n_coef = 3 + mob_lags;
	int n_walk = time_effect >= 1 ? n_weeks - 1 : 0;
	int n_interaction = time_effect == 2 ? n_weeks - 1 : 0;
	int fixed = lag_dist == 0 && n_orders == 1;                      // pressures known: the replication's fast path
	matrix[n_days, n_npi] npi = restrictions[, npi_index];
	int n_train = n_areas * n_days - sum(to_array_1d(holdout));
	array[n_train] int train_area;
	array[n_train] int train_day;
	array[n_train] int train_cases;
	matrix[fixed ? n_areas : 0, fixed ? n_days : 0] ar_fixed;
	matrix[fixed ? n_areas : 0, fixed ? n_days : 0] ne_fixed;
	array[mob_lags] matrix[fixed ? n_areas : 0, fixed ? n_days : 0] mob_fixed;
	vector[n_blocks] block_n = rep_vector(0, n_blocks);
	for (t in 1:n_days) block_n[day_block[t]] += 1;
	if (phi_area && !negbin) reject("phi_area needs negbin");
	if (n_cov > 0 && pooled == 0) reject("covariates need pooled coefficients");
	if (fixed) {
		ar_fixed = to_matrix(autoregressive_lags[, 1], n_areas, n_days);
		ne_fixed = to_matrix(neighbour_lags[1][, 1], n_areas, n_days);
		for (d in 1:mob_lags) mob_fixed[d] = to_matrix(mobility_lags[d][, 1], n_areas, n_days);
	}
	{
		int k = 0;
		for (i in 1:n_areas) for (t in 1:n_days) if (holdout[i, t] == 0) {
			k += 1;
			train_area[k] = i;
			train_day[k] = t;
			train_cases[k] = cases[i, t];
		}
	}
}

parameters {
	matrix<lower=0>[pooled ? 0 : n_coef, n_areas] coef_free;
	matrix[pooled ? n_coef : 0, n_areas] coef_z;
	vector[pooled ? n_coef : 0] coef_mu;
	vector<lower=0>[pooled ? n_coef : 0] coef_sd;
	vector<lower=0, upper=1>[pooled == 2 ? n_coef - 1 : 0] coef_lambda;   // LCAR dependence of the epidemic coefficients
	vector[pooled == 2 ? n_areas : 0] coef_log_s;                       // log s_i of the coefficients
	array[pooled == 2] real<lower=0> coef_s_sd;
	array[pooled == 2] real<lower=0, upper=1> coef_s_lambda;
	vector[n_cov] beta_cov;                                             // effect of each covariate on its log coefficient
	vector[6] weekday_raw;
	vector[n_walk] gamma_z;
	array[time_effect >= 1] real<lower=0> gamma_sd;
	matrix[n_areas, n_interaction] delta_z;
	array[time_effect >= 2] real<lower=0> delta_sd;
	matrix[time_effect == 3 ? n_weeks : 0, time_effect == 3 ? n_areas : 0] eta;   // temporal pCAR per ABS (column)
	array[time_effect == 3] real logit_c0;                              // typical temporal dependence
	vector[time_effect == 3 ? n_areas : 0] c_raw;                       // logit c_i - logit_c0, LCAR
	array[time_effect == 3] real<lower=0> c_sd;
	array[time_effect == 3] real<lower=0, upper=1> c_lambda;
	vector[time_effect == 3 ? n_areas : 0] delta_log_s;                 // log s_i of the interaction, LCAR
	array[time_effect == 3] real<lower=0> delta_s_sd;
	array[time_effect == 3] real<lower=0, upper=1> delta_s_lambda;
	array[vaccination] real tau;
	array[testing > 0] real beta_test;
	array[negbin && !phi_area] real<lower=0> inv_sqrt_phi;
	array[phi_area] real log_phi_mu;                                    // typical log phi_i
	array[phi_area] real<lower=0> log_phi_sd;
	vector[phi_area ? n_areas : 0] phi_z;
	vector[n_npi] beta_npi;
	array[mob_closure] real beta_mob_closure;
	array[lag_dist] real<lower=0> lag_mean;                             // mean of s - lag_min
	array[lag_dist] real<lower=0> lag_inv_sqrt_phi;
	array[n_orders > 1] real<lower=0> ne_decay;                         // d of o^-d
	vector[holiday ? 2 : 0] beta_holiday;
	matrix[day_effect ? (day_effect == 2 ? n_day_groups : 1) : 0, n_days] day_z;
	array[day_effect > 0] real<lower=0> day_sd;
}

transformed parameters {
	// rows: endemic, autoregressive, neighbour, mobility lag 1, mobility lag 2(, mobility lag 3)
	matrix<lower=0>[n_coef, n_areas] coef;
	vector[7] weekday_effect = exp(append_row(weekday_raw, -sum(weekday_raw)));
	vector[n_weeks] gamma = rep_vector(0, n_weeks);
	matrix[n_areas, n_weeks] delta = rep_matrix(0, n_areas, n_weeks);
	vector[time_effect == 3 ? n_areas : 0] c_area;                     // c_i
	vector[negbin ? n_areas : 0] phi;                                   // size of each ABS
	vector[lag_dist ? n_lags : 1] lag_weight = rep_vector(1, lag_dist ? n_lags : 1);
	vector[n_orders] order_weight = rep_vector(1, n_orders);
	matrix[day_effect ? (day_effect == 2 ? n_day_groups : 1) : 0, n_days] day_log;   // e_t, or e_g,t by group
	if (pooled == 0) coef = coef_free;
	else {
		matrix[n_coef, n_areas] z = coef_z;
		if (pooled == 2) z[2:n_coef] = diag_post_multiply(z[2:n_coef], exp(coef_log_s));
		matrix[n_coef, n_areas] lm = rep_matrix(coef_mu, n_areas) + diag_pre_multiply(coef_sd, z);
		for (c in 1:n_cov) lm[cov_row[c]] += beta_cov[c] * covariate[, c]';
		coef = exp(lm);
	}
	if (time_effect >= 1) {
		// RW1: weekly steps gamma_sd z; RW2: the steps themselves a random walk (second differences gamma_sd z)
		gamma = gamma_sd[1] * append_row(0, cumulative_sum(rw_order == 1 ? gamma_z : cumulative_sum(gamma_z)));
		gamma -= mean(gamma);
	}
	if (time_effect == 2) {
		matrix[n_areas, n_weeks] walk = rep_matrix(0, n_areas, n_weeks);
		for (w in 2:n_weeks) walk[, w] = walk[, w - 1] + delta_z[, w - 1];
		delta = delta_sd[1] * (walk - rep_matrix(walk * rep_vector(1.0 / n_weeks, n_weeks), n_weeks));
	}
	if (time_effect == 3) {
		// sigma s_i eta_i,t, centred over ABS and over weeks
		matrix[n_weeks, n_areas] e = diag_post_multiply(eta, delta_sd[1] * exp(delta_log_s));
		c_area = inv_logit(logit_c0[1] + c_sd[1] * c_raw);
		e += -rep_matrix(e * rep_vector(1.0 / n_areas, n_areas), n_areas)
			- rep_matrix(rep_row_vector(1.0 / n_weeks, n_weeks) * e, n_weeks) + mean(e);
		delta = e';
	}
	if (negbin) phi = phi_area ? exp(log_phi_mu[1] + log_phi_sd[1] * phi_z) : rep_vector(inv_square(inv_sqrt_phi[1]), n_areas);
	if (lag_dist) lag_weight = lag_weights(lag_mean[1], inv_square(lag_inv_sqrt_phi[1]), n_lags, lag_min);
	if (n_orders > 1) for (o in 2:n_orders) order_weight[o] = pow(o, -ne_decay[1]);
	for (g in 1:rows(day_z)) {
		vector[n_blocks] m = rep_vector(0, n_blocks);
		for (t in 1:n_days) m[day_block[t]] += day_z[g, t];
		day_log[g] = day_sd[1] * (day_z[g] - (m ./ block_n)[day_block]');
	}
}

model {
	{
		matrix[n_areas, n_days] p = transmission(gamma, delta, tau, vaccinated, week_lo, week_hi, week_frac, npi * beta_npi);
		vector[n_days] mob_factor = mob_closure ? exp(beta_mob_closure[1] * closure) : rep_vector(1, n_days);
		vector[n_days] report = weekday_effect[weekday];
		array[4] matrix[n_areas, n_days] x;
		matrix[n_areas, n_days] expected;
		if (fixed) x = components(coef, pop, p, ar_fixed, ne_fixed, mob_fixed, mob_factor);
		else {
			array[mob_lags] matrix[n_areas, n_days] mob;
			for (d in 1:mob_lags) mob[d] = lagged(mobility_lags[d], lag_weight, n_areas, n_days);
			x = components(coef, pop, p, lagged(autoregressive_lags, lag_weight, n_areas, n_days),
				neighbours(neighbour_lags, lag_weight, order_weight, n_areas, n_days), mob, mob_factor);
		}
		if (testing == 2) expected = x[1] .* exp(beta_test[1] * log_tests)
			+ (x[2] + x[3] + x[4]) .* exp(beta_test[1] * (log_tests - log_tests_lag));
		else expected = x[1] + x[2] + x[3] + x[4];
		if (testing == 1) expected = expected .* exp(beta_test[1] * log_tests);
		if (holiday) report = report .* exp(beta_holiday[1] * holiday_day + beta_holiday[2] * after_holiday);
		expected = diag_post_multiply(expected, report);
		if (day_effect) expected = expected .* exp(day_log[day_effect == 2 ? day_group : rep_array(1, n_areas)]);
		if (n_train == n_areas * n_days) {
			if (negbin) for (i in 1:n_areas) target += neg_binomial_2_lpmf(cases[i] | expected[i], phi[i]);
			else for (i in 1:n_areas) target += poisson_lpmf(cases[i] | expected[i]);
		} else {
			vector[n_train] mu;
			for (k in 1:n_train) mu[k] = expected[train_area[k], train_day[k]];
			if (negbin) target += neg_binomial_2_lpmf(train_cases | mu, phi[train_area]);
			else target += poisson_lpmf(train_cases | mu);
		}
	}

	if (pooled == 0) for (k in 1:n_coef) coef_free[k] ~ normal(0, prior_scale[k]);
	else {
		coef_mu ~ normal(log(prior_scale), 3);    // ours: centred on their scales, wide
		coef_sd ~ normal(0, 1);
	}
	beta_cov ~ normal(0, 0.5);                   // ours: up to about 2.7-fold per sd
	if (pooled == 1) to_vector(coef_z) ~ std_normal();
	if (pooled == 2) {
		coef_z[1] ~ std_normal();
		for (k in 2:n_coef) {
			coef_z[k]' ~ lcar(coef_lambda[k - 1], 1, node1, node2, laplacian_eigen);
			target += normal_lpdf(sum(coef_z[k]) | 0, 0.001 * n_areas);   // mu_k is the typical log coefficient
		}
		coef_log_s ~ lcar(coef_s_lambda[1], coef_s_sd[1], node1, node2, laplacian_eigen);
		target += normal_lpdf(sum(coef_log_s) | 0, 0.001 * n_areas);       // geometric mean of s_i is 1
		coef_lambda ~ beta(1, 1);
		coef_s_lambda ~ beta(1, 1);
		coef_s_sd ~ exponential(2);
	}
	weekday_raw ~ normal(0, 0.5);
	gamma_z ~ std_normal();
	gamma_sd ~ normal(0, gamma_prior);           // ours: sd of the weekly change (RW1) or of its change (RW2) in log p
	beta_npi ~ normal(0, 1);                     // ours: log p up to about 2.7-fold with a restriction in force
	beta_mob_closure ~ normal(0, 1);
	if (time_effect == 2) {
		// scaled ICAR on each week's increments, with a soft sum-to-zero over ABS
		target += -0.5 * icar_scale * dot_self(to_vector(delta_z[node1] - delta_z[node2]));
		to_vector(rep_row_vector(1, n_areas) * delta_z) ~ normal(0, 0.001 * n_areas);
	}
	delta_sd ~ normal(0, 0.5);
	if (time_effect == 3) {
		eta ~ tpcar(c_area, week_n_neighbours, week_eigen);
		logit_c0 ~ normal(1, 1);                 // ours: typical c about 0.73
		c_raw ~ lcar(c_lambda[1], 1, node1, node2, laplacian_eigen);
		target += normal_lpdf(sum(c_raw) | 0, 0.001 * n_areas);
		c_sd ~ exponential(2);
		c_lambda ~ beta(1, 1);
		delta_log_s ~ lcar(delta_s_lambda[1], delta_s_sd[1], node1, node2, laplacian_eigen);
		target += normal_lpdf(sum(delta_log_s) | 0, 0.001 * n_areas);
		delta_s_sd ~ exponential(2);
		delta_s_lambda ~ beta(1, 1);
	}
	tau ~ normal(0, 2);
	beta_test ~ normal(0, 1);                    // theirs
	inv_sqrt_phi ~ normal(0, 1);
	log_phi_mu ~ normal(1.5, 1);                 // ours: typical phi about 4.5 (the single phi of the Catalunya fit: 5.4)
	log_phi_sd ~ normal(0, 1);
	phi_z ~ std_normal();
	lag_mean ~ normal(0, 4);                     // ours: half-normal, mean lag lag_min + 0..8 days
	lag_inv_sqrt_phi ~ normal(0, 1);
	ne_decay ~ normal(0, 3);                     // ours: half-normal; 0 = all orders equal, large = first order only
	beta_holiday ~ normal(0, 1);
	to_vector(day_z) ~ std_normal();
	day_sd ~ normal(0, 0.5);                     // ours: half-normal, sd of the daily log factor mostly under 1
}

generated quantities {
	// components without the testing and reporting factors, as daily.stan, summed over ABS for each day and over days
	// for each ABS; expected_day with all factors, for the comparison with observed cases
	vector[n_days] endemic_day;
	vector[n_days] autoregressive_day;
	vector[n_days] neighbour_day;
	vector[n_days] mobility_day;
	vector[n_days] expected_day;
	vector[n_areas] mobility_area;
	vector[n_areas] infections_per_trip_area;
	vector[n_days] infections_per_trip_day;
	real infections_per_trip;
	matrix[save_expected ? n_areas : 0, save_expected ? n_days : 0] expected_cell;
	array[n_areas, n_days] real log_lik;
	{
		matrix[n_areas, n_days] p = transmission(gamma, delta, tau, vaccinated, week_lo, week_hi, week_frac, npi * beta_npi);
		vector[n_days] mob_factor = mob_closure ? exp(beta_mob_closure[1] * closure) : rep_vector(1, n_days);
		vector[n_days] report = weekday_effect[weekday];
		array[mob_lags] matrix[n_areas, n_days] mob;
		array[4] matrix[n_areas, n_days] x;
		matrix[n_areas, n_days] pressure;
		matrix[n_areas, n_days] expected;
		row_vector[n_areas] ones = rep_row_vector(1, n_areas);
		for (d in 1:mob_lags) mob[d] = lagged(mobility_lags[d], lag_weight, n_areas, n_days);
		x = components(coef, pop, p, lagged(autoregressive_lags, lag_weight, n_areas, n_days),
			neighbours(neighbour_lags, lag_weight, order_weight, n_areas, n_days), mob, mob_factor);
		pressure = mob[1];
		for (d in 2:mob_lags) pressure += mob[d];
		if (testing == 2) expected = x[1] .* exp(beta_test[1] * log_tests)
			+ (x[2] + x[3] + x[4]) .* exp(beta_test[1] * (log_tests - log_tests_lag));
		else expected = x[1] + x[2] + x[3] + x[4];
		if (testing == 1) expected = expected .* exp(beta_test[1] * log_tests);
		if (holiday) report = report .* exp(beta_holiday[1] * holiday_day + beta_holiday[2] * after_holiday);
		expected = diag_post_multiply(expected, report);
		if (day_effect) expected = expected .* exp(day_log[day_effect == 2 ? day_group : rep_array(1, n_areas)]);
		endemic_day = (ones * x[1])';
		autoregressive_day = (ones * x[2])';
		neighbour_day = (ones * x[3])';
		mobility_day = (ones * x[4])';
		expected_day = (ones * expected)';
		mobility_area = x[4] * rep_vector(1, n_days);
		infections_per_trip_area = mobility_area ./ (pressure * rep_vector(1, n_days));
		infections_per_trip_day = mobility_day ./ (ones * pressure)';
		infections_per_trip = sum(x[4]) / sum(pressure);
		if (save_expected) expected_cell = expected;
		for (i in 1:n_areas) for (t in 1:n_days)
			log_lik[i, t] = negbin ? neg_binomial_2_lpmf(cases[i, t] | expected[i, t], phi[i])
				: poisson_lpmf(cases[i, t] | expected[i, t]);
	}
}
