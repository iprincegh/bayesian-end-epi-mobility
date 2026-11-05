// Daily version of Slater et al. (2025) model 1 (stan/weekly.stan): the same components, weights and priors. Changes for
// daily counts:
// (1) cases at lags s = 1..S days, weighted by lag_weight u_s, in place of last week's cases: either fixed (equal
//     weights on days 7-13, their weekly Y_t-1 in days) or proportional to Poisson(s - 1 | kappa) (Bracher and Held 2022);
// (2) mobility lag 1 = trips summed over the 7 days up to each case day, lag 2 = the 7 days before them, in place of
//     their weeks t - 1 and t - 2;
// (3) a day-of-week reporting factor on the whole mean (geometric mean 1), entering as their testing factor does.
// The pressure terms are computed in R (helpers.R, daily_inputs) for each lag, as weekly.stan's transformed data block.
functions {
	vector lag_weights(real kappa, int n_lags) {
		vector[n_lags] lw;
		for (s in 1:n_lags) lw[s] = poisson_lpmf(s - 1 | kappa);
		return softmax(lw);
	}
	matrix lagged(array[] matrix x, vector lag_weight) {
		matrix[rows(x[1]), cols(x[1])] m = lag_weight[1] * x[1];
		for (s in 2:size(x)) m += lag_weight[s] * x[s];
		return m;
	}
	// expected cases from the lag-weighted sums: (endemic N_i + autoregressive_i AR + N_i neighbour_i NE
	// + mobility_1,i MOB1 + mobility_2,i MOB2) x weekday factor
	matrix expected_counts(vector endemic_pop, vector autoregressive, vector neighbour_pop, array[] vector mobility,
			matrix ar, matrix ne, matrix mob1, matrix mob2, vector weekday_factor) {
		return diag_post_multiply(rep_matrix(endemic_pop, cols(ar)) + diag_pre_multiply(autoregressive, ar)
			+ diag_pre_multiply(neighbour_pop, ne) + diag_pre_multiply(mobility[1], mob1)
			+ diag_pre_multiply(mobility[2], mob2), weekday_factor);
	}
}

data {
	int<lower=0> n_areas;
	int<lower=0> n_days;
	int<lower=1> n_lags;
	array[n_areas, n_days] int<lower=0> cases;                           // Y_i,t
	array[n_lags] matrix<lower=0>[n_areas, n_days] lagged_cases;         // Y_i,t-s
	array[n_lags] matrix<lower=0>[n_areas, n_days] neighbour_pressure;   // sum_j v_ji Y_j,t-s / N_j
	array[n_lags] matrix<lower=0>[n_areas, n_days] mobility_pressure1;   // sum_j w_ji Y_j,t-s / N_j, trips of days t-s-6..t-s
	array[n_lags] matrix<lower=0>[n_areas, n_days] mobility_pressure2;   // the same, trips of days t-s-13..t-s-7
	array[n_areas] int<lower=0> pop;                                     // N_i
	array[n_days] int<lower=1, upper=7> weekday;                         // Monday = 1
	int<lower=0, upper=1> estimate_lag;                                  // 1: Poisson lag; 0: fixed_lag_weight
	simplex[n_lags] fixed_lag_weight;
}

transformed data {
	vector[n_areas] pop_vector = to_vector(pop);
	// with a fixed lag the lag-weighted sums are data: computed once here, not at every gradient evaluation
	matrix[n_areas, n_days] fixed_ar = lagged(lagged_cases, fixed_lag_weight);
	matrix[n_areas, n_days] fixed_ne = lagged(neighbour_pressure, fixed_lag_weight);
	matrix[n_areas, n_days] fixed_mob1 = lagged(mobility_pressure1, fixed_lag_weight);
	matrix[n_areas, n_days] fixed_mob2 = lagged(mobility_pressure2, fixed_lag_weight);
}

parameters {
	vector<lower=0>[n_areas] endemic;
	array[2] vector<lower=0>[n_areas] mobility;
	vector<lower=0>[n_areas] neighbour;
	vector<lower=0>[n_areas] autoregressive;
	array[estimate_lag] real<lower=0> lag_kappa;
	vector[6] weekday_raw;
}

transformed parameters {
	vector[n_lags] lag_weight = fixed_lag_weight;
	if (estimate_lag) lag_weight = lag_weights(lag_kappa[1], n_lags);
	vector[7] weekday_effect = exp(append_row(weekday_raw, -sum(weekday_raw)));
}

model {
	matrix[n_areas, n_days] expected;                // pop[i] factor on the neighbour term, as weekly.stan
	if (estimate_lag)
		expected = expected_counts(endemic .* pop_vector, autoregressive, pop_vector .* neighbour, mobility,
			lagged(lagged_cases, lag_weight), lagged(neighbour_pressure, lag_weight), lagged(mobility_pressure1, lag_weight),
			lagged(mobility_pressure2, lag_weight), weekday_effect[weekday]);
	else
		expected = expected_counts(endemic .* pop_vector, autoregressive, pop_vector .* neighbour, mobility,
			fixed_ar, fixed_ne, fixed_mob1, fixed_mob2, weekday_effect[weekday]);
	for (i in 1:n_areas) target += poisson_lpmf(cases[i] | expected[i]);

	endemic ~ normal(0, 0.01);
	autoregressive ~ normal(0, 2);
	mobility[1] ~ normal(0, 2);
	mobility[2] ~ normal(0, 2);
	neighbour ~ normal(0, 500);
	lag_kappa ~ normal(4, 2);              // ours: mean lag about 5 days
	weekday_raw ~ normal(0, 0.5);          // ours
}

generated quantities {
	// components without the weekday factor (their components also leave out the testing factor), summed over areas
	// for each day and over days for each area, to keep the output small
	vector[n_days] endemic_day;
	vector[n_days] autoregressive_day;
	vector[n_days] neighbour_day;
	vector[n_days] mobility_day;
	vector[n_areas] mobility_area;
	vector[n_areas] infections_per_trip_area;
	vector[n_days] infections_per_trip_day;
	real infections_per_trip;
	array[n_areas, n_days] real log_lik;
	{
		matrix[n_areas, n_days] pressure = lagged(mobility_pressure1, lag_weight) + lagged(mobility_pressure2, lag_weight);
		matrix[n_areas, n_days] autoregressive_part = diag_pre_multiply(autoregressive, lagged(lagged_cases, lag_weight));
		matrix[n_areas, n_days] neighbour_part = diag_pre_multiply(pop_vector .* neighbour, lagged(neighbour_pressure, lag_weight));
		matrix[n_areas, n_days] mobility_part = diag_pre_multiply(mobility[1], lagged(mobility_pressure1, lag_weight))
			+ diag_pre_multiply(mobility[2], lagged(mobility_pressure2, lag_weight));
		matrix[n_areas, n_days] expected = diag_post_multiply(rep_matrix(endemic .* pop_vector, n_days)
			+ autoregressive_part + neighbour_part + mobility_part, weekday_effect[weekday]);
		endemic_day = rep_vector(sum(endemic .* pop_vector), n_days);
		autoregressive_day = (rep_row_vector(1, n_areas) * autoregressive_part)';
		neighbour_day = (rep_row_vector(1, n_areas) * neighbour_part)';
		mobility_day = (rep_row_vector(1, n_areas) * mobility_part)';
		mobility_area = mobility_part * rep_vector(1, n_days);
		infections_per_trip_area = mobility_area ./ (pressure * rep_vector(1, n_days));
		infections_per_trip_day = mobility_day ./ (rep_row_vector(1, n_areas) * pressure)';
		infections_per_trip = sum(mobility_part) / sum(pressure);
		for (i in 1:n_areas) for (t in 1:n_days) log_lik[i, t] = poisson_lpmf(cases[i, t] | expected[i, t]);
	}
}
