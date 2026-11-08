// Endemic-epidemic space-time model of daily cases by ABS. The final model is MODEL in helpers.R; the other
// switches give the variants it was compared with.
functions {
	// Leroux CAR: precision (lambda (D - A) + (1 - lambda) I) / sigma^2, with the eigenvalues of D - A
	real lcar_lpdf(vector x, real lambda, real sigma, array[] int node1, array[] int node2, vector laplacian_eigen) {
		return 0.5 * sum(log(lambda * laplacian_eigen + 1 - lambda)) - rows(x) * log(sigma)
			- 0.5 / square(sigma) * (lambda * dot_self(x[node1] - x[node2]) + (1 - lambda) * dot_self(x));
	}
	// temporal proper CAR per area (MacNab 2022, eq. 26)
	real tpcar_lpdf(matrix x, vector c, vector week_n_neighbours, vector week_eigen) {
		int n_weeks = rows(x);
		row_vector[cols(x)] q = week_n_neighbours' * square(x)
			- 2 * c' .* columns_dot_product(x[1:(n_weeks - 1)], x[2:n_weeks]);
		return 0.5 * sum(log1m(week_eigen * c')) - 0.5 * sum(q);
	}
	// p_i,t: weekly log values interpolated between week midpoints
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
	// neighbour pressure with power-law weights o^-d over orders 1-5, normalised over all orders of each origin (Meyer
	// and Held 2014): W prev, W[i, j] = o_ij^-d / norm_j for the pairs within order 5
	matrix powerlaw_neighbours(matrix prev, real d, array[] int dest, array[] int orig, array[] int ord, matrix n_by_order) {
		int n = rows(prev);
		vector[5] ow;
		for (o in 1:5) ow[o] = pow(o, -d);
		vector[n] r = inv(n_by_order * ow);
		matrix[n, n] w = rep_matrix(0, n, n);
		for (k in 1:size(ord)) w[dest[k], orig[k]] = ow[ord[k]] * r[orig[k]];
		return w * prev;
	}
	// expected cases by component: endemic, AR, NE, MOB
	array[] matrix components(matrix coef, vector pop, matrix p, matrix ar, matrix ne, array[] matrix mob,
			vector mob_factor) {
		array[4] matrix[rows(ar), cols(ar)] x;
		matrix[rows(ar), cols(ar)] m = diag_pre_multiply(coef[4]', mob[1]);
		for (d in 2:size(mob)) m += diag_pre_multiply(coef[3 + d]', mob[d]);
		x[1] = rep_matrix(coef[1]' .* pop, cols(ar));
		x[2] = diag_pre_multiply(coef[2]', ar) .* p;
		x[3] = diag_pre_multiply(coef[3]' .* pop, ne) .* p;
		x[4] = diag_post_multiply(m .* p, mob_factor);
		return x;
	}
}

data {
	int<lower=0> n_areas;
	int<lower=0> n_days;
	array[n_areas, n_days] int<lower=0> cases;
	int<lower=1> n_lags;
	int<lower=0, upper=1> lag_dist;
	int<lower=1, upper=n_lags> lag_min;
	int<lower=1, upper=5> n_orders;
	int<lower=2, upper=3> mob_lags;
	// lagged pressures, rows = ABS-days (column-major area x day), columns = lags (lag_dist = 1) or one weighted column
	matrix<lower=0>[n_areas * n_days, lag_dist ? n_lags : 1] autoregressive_lags;
	array[n_orders] matrix<lower=0>[n_areas * n_days, lag_dist ? n_lags : 1] neighbour_lags;
	int<lower=0, upper=1> ne_pl;
	int<lower=0> pl_nnz;
	array[pl_nnz] int<lower=1, upper=n_areas> pl_dest;
	array[pl_nnz] int<lower=1, upper=n_areas> pl_orig;
	array[pl_nnz] int<lower=1, upper=5> pl_ord;
	matrix<lower=0>[ne_pl ? n_areas : 0, 5] pl_n_by_order;
	matrix<lower=0>[ne_pl ? n_areas : 0, ne_pl ? n_days : 0] pl_prev;
	array[mob_lags] matrix<lower=0>[n_areas * n_days, lag_dist ? n_lags : 1] mobility_lags;
	int<lower=0, upper=1> split_within;
	// the within-ABS part of each mobility pressure, for the split of the mobility channel
	array[mob_lags] matrix<lower=0>[split_within ? n_areas * n_days : 0, lag_dist ? n_lags : 1] within_lags;
	vector<lower=0>[n_areas] pop;
	array[n_days] int<lower=1, upper=7> weekday;
	vector<lower=0>[3 + mob_lags] prior_scale;
	int<lower=0, upper=3> time_effect;
	int<lower=0, upper=1> vaccination;
	int<lower=0, upper=2> testing;
	int<lower=0, upper=2> pooled;
	int<lower=0> n_cov;
	matrix[n_areas, n_cov] covariate;
	array[n_cov] int<lower=1, upper=3 + mob_lags> cov_row;
	int<lower=0, upper=1> negbin;
	int<lower=0, upper=1> phi_area;
	int<lower=1> n_weeks;
	array[n_days] int<lower=1, upper=n_weeks> week_lo;
	array[n_days] int<lower=1, upper=n_weeks> week_hi;
	vector<lower=0, upper=1>[n_days] week_frac;
	vector<lower=1>[n_weeks] week_n_neighbours;
	vector[n_weeks] week_eigen;
	int<lower=0> n_edges;
	array[n_edges] int<lower=1, upper=n_areas> node1;
	array[n_edges] int<lower=1, upper=n_areas> node2;
	real<lower=0> icar_scale;
	vector[n_areas] laplacian_eigen;
	matrix<lower=0, upper=1>[n_areas, n_days] vaccinated;
	matrix<upper=0>[n_areas, n_days] log_tests;
	matrix[n_areas, n_days] log_tests_lag;
	int<lower=1, upper=2> rw_order;
	real<lower=0> gamma_prior;
	int<lower=1> n_restrictions;
	matrix<lower=0, upper=1>[n_days, n_restrictions] restrictions;
	int<lower=0, upper=n_restrictions> n_npi;
	array[n_npi] int<lower=1, upper=n_restrictions> npi_index;
	int<lower=0, upper=1> mob_closure;
	vector<lower=0, upper=1>[n_days] closure;
	int<lower=0, upper=1> holiday;
	vector<lower=0, upper=1>[n_days] holiday_day;
	vector<lower=0, upper=1>[n_days] after_holiday;
	int<lower=0, upper=2> day_effect;
	int<lower=1> n_day_groups;
	array[n_areas] int<lower=1, upper=n_day_groups> day_group;
	int<lower=1> n_blocks;
	array[n_days] int<lower=1, upper=n_blocks> day_block;
	int<lower=0, upper=1> save_expected;
	array[n_areas, n_days] int<lower=0, upper=1> holdout;
}

transformed data {
	int n_coef = 3 + mob_lags;
	int n_walk = time_effect >= 1 ? n_weeks - 1 : 0;
	int n_interaction = time_effect == 2 ? n_weeks - 1 : 0;
	int fixed = lag_dist == 0 && n_orders == 1;
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
	if (ne_pl && (n_orders > 1 || lag_dist)) reject("ne_pl replaces n_orders > 1 and needs the fixed lag");
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
	vector<lower=0, upper=1>[pooled == 2 ? n_coef - 1 : 0] coef_lambda;
	vector[pooled == 2 ? n_areas : 0] coef_log_s;
	array[pooled == 2] real<lower=0> coef_s_sd;
	array[pooled == 2] real<lower=0, upper=1> coef_s_lambda;
	vector[n_cov] beta_cov;
	vector[6] weekday_raw;
	vector[n_walk] gamma_z;
	array[time_effect >= 1] real<lower=0> gamma_sd;
	matrix[n_areas, n_interaction] delta_z;
	array[time_effect >= 2] real<lower=0> delta_sd;
	matrix[time_effect == 3 ? n_weeks : 0, time_effect == 3 ? n_areas : 0] eta;
	array[time_effect == 3] real logit_c0;
	vector[time_effect == 3 ? n_areas : 0] c_raw;
	array[time_effect == 3] real<lower=0> c_sd;
	array[time_effect == 3] real<lower=0, upper=1> c_lambda;
	vector[time_effect == 3 ? n_areas : 0] delta_log_s;
	array[time_effect == 3] real<lower=0> delta_s_sd;
	array[time_effect == 3] real<lower=0, upper=1> delta_s_lambda;
	array[vaccination] real tau;
	array[testing > 0] real beta_test;
	array[negbin && !phi_area] real<lower=0> inv_sqrt_phi;
	array[phi_area] real log_phi_mu;
	array[phi_area] real<lower=0> log_phi_sd;
	vector[phi_area ? n_areas : 0] phi_z;
	vector[n_npi] beta_npi;
	array[mob_closure] real beta_mob_closure;
	array[lag_dist] real<lower=0> lag_mean;
	array[lag_dist] real<lower=0> lag_inv_sqrt_phi;
	array[n_orders > 1 || ne_pl] real<lower=0> ne_decay;
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
	vector[time_effect == 3 ? n_areas : 0] c_area;
	vector[negbin ? n_areas : 0] phi;
	vector[lag_dist ? n_lags : 1] lag_weight = rep_vector(1, lag_dist ? n_lags : 1);
	vector[n_orders] order_weight = rep_vector(1, n_orders);
	matrix[day_effect ? (day_effect == 2 ? n_day_groups : 1) : 0, n_days] day_log;
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
		if (fixed) x = components(coef, pop, p, ar_fixed,
			ne_pl ? powerlaw_neighbours(pl_prev, ne_decay[1], pl_dest, pl_orig, pl_ord, pl_n_by_order) : ne_fixed, mob_fixed, mob_factor);
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
		coef_mu ~ normal(log(prior_scale), 3);
		coef_sd ~ normal(0, 1);
	}
	beta_cov ~ normal(0, 0.5);
	if (pooled == 1) to_vector(coef_z) ~ std_normal();
	if (pooled == 2) {
		coef_z[1] ~ std_normal();
		for (k in 2:n_coef) {
			coef_z[k]' ~ lcar(coef_lambda[k - 1], 1, node1, node2, laplacian_eigen);
			target += normal_lpdf(sum(coef_z[k]) | 0, 0.001 * n_areas);
		}
		coef_log_s ~ lcar(coef_s_lambda[1], coef_s_sd[1], node1, node2, laplacian_eigen);
		target += normal_lpdf(sum(coef_log_s) | 0, 0.001 * n_areas);
		coef_lambda ~ beta(1, 1);
		coef_s_lambda ~ beta(1, 1);
		coef_s_sd ~ exponential(2);
	}
	weekday_raw ~ normal(0, 0.5);
	gamma_z ~ std_normal();
	gamma_sd ~ normal(0, gamma_prior);
	beta_npi ~ normal(0, 1);
	beta_mob_closure ~ normal(0, 1);
	if (time_effect == 2) {
		// scaled ICAR on each week's increments, with a soft sum-to-zero over ABS
		target += -0.5 * icar_scale * dot_self(to_vector(delta_z[node1] - delta_z[node2]));
		to_vector(rep_row_vector(1, n_areas) * delta_z) ~ normal(0, 0.001 * n_areas);
	}
	delta_sd ~ normal(0, 0.5);
	if (time_effect == 3) {
		eta ~ tpcar(c_area, week_n_neighbours, week_eigen);
		logit_c0 ~ normal(1, 1);
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
	beta_test ~ normal(0, 1);
	inv_sqrt_phi ~ normal(0, 1);
	log_phi_mu ~ normal(1.5, 1);
	log_phi_sd ~ normal(0, 1);
	phi_z ~ std_normal();
	lag_mean ~ normal(0, 4);
	lag_inv_sqrt_phi ~ normal(0, 1);
	ne_decay ~ normal(0, 3);
	beta_holiday ~ normal(0, 1);
	to_vector(day_z) ~ std_normal();
	day_sd ~ normal(0, 0.5);
}

generated quantities {
	// components summed over ABS (per day) and over days (per ABS)
	vector[n_days] endemic_day;
	vector[n_days] autoregressive_day;
	vector[n_days] neighbour_day;
	vector[n_days] mobility_day;
	vector[n_days] expected_day;
	vector[n_areas] mobility_area;
	vector[split_within ? n_days : 0] mobility_within_day;
	vector[ne_pl ? 5 : 0] pl_weight;
	vector[split_within ? n_areas : 0] mobility_within_area;
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
			ne_pl ? powerlaw_neighbours(pl_prev, ne_decay[1], pl_dest, pl_orig, pl_ord, pl_n_by_order)
				: neighbours(neighbour_lags, lag_weight, order_weight, n_areas, n_days), mob, mob_factor);
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
		if (ne_pl) for (o in 1:5) pl_weight[o] = pow(o, -ne_decay[1]);
		if (split_within) {
			matrix[n_areas, n_days] w = diag_pre_multiply(coef[4]', lagged(within_lags[1], lag_weight, n_areas, n_days));
			for (d in 2:mob_lags) w += diag_pre_multiply(coef[3 + d]', lagged(within_lags[d], lag_weight, n_areas, n_days));
			w = diag_post_multiply(w .* p, mob_factor);
			mobility_within_day = (ones * w)';
			mobility_within_area = w * rep_vector(1, n_days);
		}
		if (save_expected) expected_cell = expected;
		for (i in 1:n_areas) for (t in 1:n_days)
			log_lik[i, t] = negbin ? neg_binomial_2_lpmf(cases[i, t] | expected[i, t], phi[i])
				: poisson_lpmf(cases[i, t] | expected[i, t]);
	}
}
