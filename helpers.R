# Inputs for the replication of Slater et al. (2025) from the Catalunya ABS data in Spain-data, built as their
# castilla_final.R builds them (lines 14-45 and 135-213), and a report of each fit. Run the scripts from end-epi/slater.
suppressMessages({library(arrow); library(dplyr)})

SPAIN = "../Spain-data"
DAY1 = as.Date("2020-03-01")    # day 1 and the start of week 1, as their time and week = ceiling(time / 7)
FIRST_DATE = as.Date("2020-03-22")   # first modelled day: 1-21 March 2020 are left out, as their weeks 1-3
N_LAGS = 14                     # case lags of the daily model, s = 1..14 days
# lag weights over s = 1..N_LAGS: "week" = equal weights on days 7-13, their weekly Y_t-1 in days;
# "estimated" = Poisson lag with kappa estimated (the fixed weights are then unused)
LAG_WEIGHTS = list(week = c(rep(0, 6), rep(1 / 7, 7), 0), estimated = rep(1 / 14, 14))

spain = function(...) file.path(SPAIN, ...)
out = function(f) {dir.create("outputs", showWarnings = FALSE); file.path("outputs", f)}

# ABS ids, sorted: the 67 ABS of Barcelona city ("Barcelona - 01A" ...) or all 374
area_ids = function(areas = c("barcelona", "catalunya")) {
	a = read.csv(spain("outputs/catalunya/tables/cat_abs_areas.csv"), colClasses = c(id = "character")) |> arrange(id)
	ids = if (match.arg(areas) == "barcelona") a$id[startsWith(a$name, "Barcelona - ")] else a$id
	stopifnot(length(ids) == c(barcelona = 67, catalunya = 374)[areas])
	ids
}

# daily cases (area x day, day 1 = DAY1) and populations
read_cases = function(ids) {
	x = read_parquet(spain("outputs/catalunya/panel/cat_abs_main_noncollinear.parquet"), col_select = c("id", "date", "new_cases", "pop")) |>
		filter(id %in% ids) |> mutate(day = as.numeric(date - DAY1) + 1)
	cases = with(x, tapply(new_cases, list(id, day), sum))
	pop = with(distinct(x, id, pop), setNames(as.integer(round(pop)), id))[ids]
	stopifnot(identical(rownames(cases), ids), !anyNA(cases), ncol(cases) == 427, cases == round(cases))
	storage.mode(cases) = "integer"                   # Stan's int data
	list(cases = cases, pop = pop)
}

# neighbour weights v_ji = 1/Ne(j), [origin j, destination i]: queen adjacency (spdep::poly2nb, as theirs) divided by
# each origin's number of neighbours
neighbour_weights = function(ids) {
	adjacency = readRDS(spain("outputs/catalunya/geometry/cat_abs_topology.rds"))$adjacency[ids, ids]
	diag(adjacency) = 0
	stopifnot(all(rowSums(adjacency) > 0))
	adjacency / rowSums(adjacency)
}

# trips [origin, destination, period] summed by the periods in `when` (date, period), within-ABS trips included
od_sums = function(ids, when) {
	index = \(v) setNames(data.frame(ids, seq_along(ids)), c(v, paste0(v, "_index")))   # ids as integers before collect()
	f = open_dataset(spain("data/cat_od_bha_reprojected")) |>
		inner_join(index("source"), by = "source") |> inner_join(index("target"), by = "target") |>
		inner_join(when, by = "date") |>
		group_by(source_index, target_index, period) |> summarise(trips = sum(trips)) |> collect()
	m = array(0, c(length(ids), length(ids), max(when$period)), dimnames = list(ids, ids, NULL))
	m[cbind(f$source_index, f$target_index, f$period)] = f$trips
	m
}

# window(u, from, to): trips [origin, destination] summed over days u - from .. u - to (from <= 13), for the day
# numbers u in `days`, with sums under 1 set to 0 (their rule for weekly sums)
trip_windows = function(ids, days) {
	first = min(days) - 13
	f = od_sums(ids, data.frame(date = DAY1 + (first:max(days)) - 1, period = seq_along(first:max(days))))
	cum = array(0, dim(f) + c(0, 0, 1))
	for (k in seq_len(dim(f)[3])) cum[, , k + 1] = cum[, , k] + f[, , k]
	function(u, from, to) {
		w = cum[, , u - first + 2 - to] - cum[, , u - first + 1 - from]
		w[w < 1] = 0
		w
	}
}

# daily inputs for the days from .. to: for each case lag s, the sums their transformed data block forms, at day t - s,
# with trips summed over the 7 days up to the case day (lag 1) and the 7 days before them (lag 2), in place of their
# weeks t - 1 and t - 2
daily_inputs = function(areas, from = FIRST_DATE, to = as.Date("2021-05-01")) {
	stopifnot(from >= FIRST_DATE, to <= as.Date("2021-05-01"), from < to)
	ids = area_ids(areas)
	d = read_cases(ids)
	modelled = as.numeric(c(from, to) - DAY1) + 1
	modelled = modelled[1]:modelled[2]
	n = max(modelled)
	days = (min(modelled) - N_LAGS):n
	window = trip_windows(ids, days)
	prevalence = t(d$cases / d$pop)                   # day x origin, Y_j,u / N_j (row = day number)
	trip_pressure = \(from, to) t(vapply(days, \(u) drop(prevalence[u, ] %*% window(u, from, to)), numeric(length(ids))))
	lagged = \(m) lapply(seq_len(N_LAGS), \(s) t(m[match(modelled - s, days), , drop = FALSE]))   # area x day, per lag
	list(ids = ids, dates = DAY1 + modelled - 1, cases = d$cases[, modelled], pop = d$pop,
		lagged_cases = lapply(seq_len(N_LAGS), \(s) d$cases[, modelled - s]),
		neighbour_pressure = lagged(prevalence[days, ] %*% neighbour_weights(ids)),
		mobility_pressure1 = lagged(trip_pressure(6, 0)), mobility_pressure2 = lagged(trip_pressure(13, 7)))
}

# convergence, PSIS-LOO (as theirs: extract_log_lik(...) |> loo()) and trips per new infection
# loo_result: the PSIS-LOO saved by palma/convert.R, for fits stored without log_lik
report = function(fit, loo_result = NULL) {
	s = rstan::summary(fit, pars = c("endemic", "autoregressive", "neighbour", "mobility"))$summary
	cat(sprintf("max Rhat %.3f, min n_eff %.0f, divergences %d, max treedepth hits %d\n", max(s[, "Rhat"], na.rm = TRUE),
		min(s[, "n_eff"], na.rm = TRUE), rstan::get_num_divergent(fit), rstan::get_num_max_treedepth(fit)))
	print(if (is.null(loo_result)) loo::extract_log_lik(fit) |> loo::loo() else loo_result)
	cat("trips per new infection (1 / infections_per_trip):\n")
	print(quantile(1 / rstan::extract(fit)$infections_per_trip, c(0.025, 0.5, 0.975)))
}
