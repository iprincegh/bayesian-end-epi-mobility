// Slater et al. (2025) model 1 (their stan_models/multivariate_castilla_noSI.stan) with our names; same terms, priors and
// outputs. Their names: cases = y, lagged_cases = shifted_y, trips_lag1/2 = shifted_weekly_mobility1/2,
// neighbour_weights = spatial_weights, pop = pops, endemic = alpha_en, autoregressive = alpha_ar, neighbour = alpha_spat,
// mobility = alpha_move, expected = mu, *_part = EN, EP_AR, EP_spat, EP_move, infections_per_trip = inf_per_trip_agg.
// Their EP_spat leaves out the pop[i] factor that their mu has; neighbour_part here matches expected.
data {
	int<lower=0> n_areas;
	int<lower=0> n_weeks;
	array[n_areas, n_weeks] int<lower=0> cases;                 // Y_i,t
	array[n_areas, n_weeks] real<lower=0> lagged_cases;         // Y_i,t-1
	array[n_areas, n_areas, n_weeks] real<lower=0> trips_lag1;  // w_ji,t-1 [origin, destination, week]
	array[n_areas, n_areas, n_weeks] real<lower=0> trips_lag2;  // w_ji,t-2
	array[n_areas, n_areas] real<lower=0> neighbour_weights;    // v_ji = 1/Ne(j)
	array[n_areas] int<lower=0> pop;                            // N_i
}

transformed data {
	matrix[n_areas, n_weeks] mobility_pressure1;   // sum_j w_ji,t-1 Y_j,t-1 / N_j
	matrix[n_areas, n_weeks] mobility_pressure2;   // sum_j w_ji,t-2 Y_j,t-1 / N_j
	matrix[n_areas, n_weeks] neighbour_pressure;   // sum_j v_ji Y_j,t-1 / N_j
	for (t in 1:n_weeks) {
		for (dest in 1:n_areas) {
			mobility_pressure1[dest, t] = 0;
			mobility_pressure2[dest, t] = 0;
			neighbour_pressure[dest, t] = 0;
			for (orig in 1:n_areas) {
				mobility_pressure1[dest, t] += trips_lag1[orig, dest, t] * lagged_cases[orig, t] / pop[orig];
				mobility_pressure2[dest, t] += trips_lag2[orig, dest, t] * lagged_cases[orig, t] / pop[orig];
				neighbour_pressure[dest, t] += neighbour_weights[orig, dest] * lagged_cases[orig, t] / pop[orig];
			}
		}
	}
}

parameters {
	array[n_areas] real<lower=0> endemic;            // alpha^EX_i
	array[2, n_areas] real<lower=0> mobility;        // alpha^mob_i,d
	array[n_areas] real<lower=0> neighbour;          // alpha^spat_i
	array[n_areas] real<lower=0> autoregressive;     // alpha^AR_i
}

transformed parameters {
	array[n_areas, n_weeks] real<lower=0> expected;
	for (t in 1:n_weeks) {
		for (dest in 1:n_areas) {
			expected[dest, t] = endemic[dest] * pop[dest] +
				autoregressive[dest] * lagged_cases[dest, t] +
				pop[dest] * neighbour[dest] * neighbour_pressure[dest, t] +
				mobility[1, dest] * mobility_pressure1[dest, t] +
				mobility[2, dest] * mobility_pressure2[dest, t];
		}
	}
}

model {
	for (i in 1:n_areas) {
		for (t in 1:n_weeks) {
			target += poisson_lpmf(cases[i, t] | expected[i, t]);
		}
	}
	endemic ~ normal(0, 0.01);
	autoregressive ~ normal(0, 2);
	mobility[1, ] ~ normal(0, 2);
	mobility[2, ] ~ normal(0, 2);
	neighbour ~ normal(0, 500);
}

generated quantities {
	array[n_areas] real<lower=0> endemic_part;
	matrix[n_areas, n_weeks] mobility_part;
	array[n_areas, n_weeks] real<lower=0> autoregressive_part;
	array[n_areas, n_weeks] real<lower=0> neighbour_part;
	array[n_areas] real<lower=0> infections_per_trip_area;
	real<lower=0> infections_per_trip;
	array[n_areas, n_weeks] real log_lik;
	for (i in 1:n_areas) {
		endemic_part[i] = endemic[i] * pop[i];
		for (t in 1:n_weeks) {
			mobility_part[i, t] = mobility[1, i] * mobility_pressure1[i, t] + mobility[2, i] * mobility_pressure2[i, t];
			neighbour_part[i, t] = pop[i] * neighbour[i] * neighbour_pressure[i, t];
			autoregressive_part[i, t] = autoregressive[i] * lagged_cases[i, t];
			log_lik[i, t] = poisson_lpmf(cases[i, t] | expected[i, t]);
		}
		infections_per_trip_area[i] = sum(mobility_part[i, ]) / sum(mobility_pressure1[i, ] + mobility_pressure2[i, ]);
	}
	infections_per_trip = sum(mobility_part) / sum(mobility_pressure1 + mobility_pressure2);
}
