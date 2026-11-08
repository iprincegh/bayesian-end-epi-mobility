# Inputs, switches and helpers of stan/daily.stan and stan/spacetime.stan. Run from end-epi/slater.
suppressMessages({library(arrow); library(dplyr)})

SPAIN = "../Spain-data"
DAY1 = as.Date("2020-03-01")
FIRST_DATE = as.Date("2020-03-22")
N_LAGS = 14
# case-lag weights over s = 1..N_LAGS: equal weights on days 7-13
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
	storage.mode(cases) = "integer"
	list(cases = cases, pop = pop)
}

# neighbour weights 1/Ne(j), [origin, destination]
neighbour_weights = function(ids) {
	adjacency = readRDS(spain("outputs/catalunya/geometry/cat_abs_topology.rds"))$adjacency[ids, ids]
	diag(adjacency) = 0
	stopifnot(all(rowSums(adjacency) > 0))
	adjacency / rowSums(adjacency)
}

# trips [origin, destination, period] summed by period
od_sums = function(ids, when) {
	index = \(v) setNames(data.frame(ids, seq_along(ids)), c(v, paste0(v, "_index")))
	f = open_dataset(spain("data/cat_od_bha_reprojected")) |>
		inner_join(index("source"), by = "source") |> inner_join(index("target"), by = "target") |>
		inner_join(when, by = "date") |>
		group_by(source_index, target_index, period) |> summarise(trips = sum(trips)) |> collect()
	m = array(0, c(length(ids), length(ids), max(when$period)), dimnames = list(ids, ids, NULL))
	m[cbind(f$source_index, f$target_index, f$period)] = f$trips
	m
}

# trips [origin, destination] over days u - from .. u - to, sums under 1 set to 0
trip_windows = function(ids, days, max_from = 13) {
	first = min(days) - max_from
	f = od_sums(ids, data.frame(date = DAY1 + (first:max(days)) - 1, period = seq_along(first:max(days))))
	cum = array(0, dim(f) + c(0, 0, 1))
	for (k in seq_len(dim(f)[3])) cum[, , k + 1] = cum[, , k] + f[, , k]
	function(u, from, to) {
		w = cum[, , u - first + 2 - to] - cum[, , u - first + 1 - from]
		w[w < 1] = 0
		w
	}
}

# daily inputs: lagged cases and the neighbour and mobility pressures for each case lag
daily_inputs = function(areas, from = FIRST_DATE, to = as.Date("2021-05-01")) {
	stopifnot(from >= FIRST_DATE, to <= as.Date("2021-05-01"), from < to)
	ids = area_ids(areas)
	d = read_cases(ids)
	modelled = as.numeric(c(from, to) - DAY1) + 1
	modelled = modelled[1]:modelled[2]
	n = max(modelled)
	days = (min(modelled) - N_LAGS):n
	window = trip_windows(ids, days)
	prevalence = t(d$cases / d$pop)
	trip_pressure = \(from, to) t(vapply(days, \(u) drop(prevalence[u, ] %*% window(u, from, to)), numeric(length(ids))))
	lagged = \(m) lapply(seq_len(N_LAGS), \(s) t(m[match(modelled - s, days), , drop = FALSE]))
	list(ids = ids, dates = DAY1 + modelled - 1, cases = d$cases[, modelled], pop = d$pop,
		lagged_cases = lapply(seq_len(N_LAGS), \(s) d$cases[, modelled - s]),
		neighbour_pressure = lagged(prevalence[days, ] %*% neighbour_weights(ids)),
		mobility_pressure1 = lagged(trip_pressure(6, 0)), mobility_pressure2 = lagged(trip_pressure(13, 7)))
}

# Inputs of stan/spacetime.stan beyond those of daily.stan

# province of each ABS (the one holding most of its population)
abs_province = function(ids) {
	o = readRDS(spain("data/overlay_cat.rds")) |> mutate(province = substr(mitma_id, 1, 2)) |>
		group_by(bha_id, province) |> summarise(pop = sum(pop_piece), .groups = "drop_last") |>
		slice_max(pop, n = 1, with_ties = FALSE)
	p = setNames(o$province, o$bha_id)[ids]
	stopifnot(!anyNA(p), p %in% c("08", "17", "25", "43"))
	p
}

# log(tests / max) in the province of each ABS, 7-day centred mean; lagged: averaged over the case lags
log_tests = function(ids, dates, lagged = FALSE) {
	x = read.delim("mobilityCovid/data/Datos_Pruebas_Realizadas_Historico_11022022.csv", sep = ";", fileEncoding = "latin1")
	code = c(Barcelona = "08", Girona = "17", Lleida = "25", Tarragona = "43")
	x = x[x[[1]] %in% names(code), ]
	day = as.Date(sprintf("%s-%02d-%s", substr(x[[2]], 6, 9), match(substr(x[[2]], 3, 5), toupper(month.abb)), substr(x[[2]], 1, 2)))
	tests = as.numeric(x[[4]]) + as.numeric(sub(",$", "", x[[6]]))
	m = tapply(tests, list(code[x[[1]]], as.character(day)), sum)
	stopifnot(!anyNA(m), all(as.character(c(min(dates) - 3, max(dates) + 3)) %in% colnames(m)))
	s = t(apply(m, 1, \(v) stats::filter(v, rep(1 / 7, 7), sides = 2)))
	dimnames(s) = dimnames(m)
	at = \(u) log(s[, as.character(pmax(u, as.Date(colnames(s)[4]))), drop = FALSE])
	l = if (lagged) Reduce(`+`, Map(\(k, w) w * at(dates - k), seq_len(N_LAGS), LAG_WEIGHTS$week)) else at(dates)
	l = l - apply(at(dates), 1, max)
	`rownames<-`(l[abs_province(ids), , drop = FALSE], ids)
}

# share of each ABS with a first dose 14 days earlier
vaccinated_share = function(ids, dates) {
	x = read_parquet(spain("outputs/catalunya/panel/cat_abs_main_noncollinear.parquet"), col_select = c("id", "date", "vacc_dose1")) |>
		filter(id %in% ids)
	m = with(x, tapply(vacc_dose1, list(id, as.character(date)), sum))[ids, as.character(dates - 14)]
	stopifnot(!anyNA(m), m >= 0, m <= 1)
	m
}

# weekly knots at the week midpoints and the position of each day between them
week_knots = function(n_days) {
	n_weeks = ceiling(n_days / 7)
	mid = 7 * seq_len(n_weeks) - 3
	k = seq_len(n_days)
	lo = pmax(findInterval(k, mid), 1)
	hi = pmin(lo + 1, n_weeks)
	# chain of weeks (temporal proper CAR)
	w = outer(seq_len(n_weeks), seq_len(n_weeks), \(a, b) abs(a - b) == 1) * 1
	n = rowSums(w)
	list(n_weeks = n_weeks, week_lo = lo, week_hi = hi, week_frac = ifelse(k <= mid[1] | k >= mid[n_weeks], 0, (k - mid[lo]) / 7),
		week_n_neighbours = n, week_eigen = eigen(w / sqrt(outer(n, n)), symmetric = TRUE)$values)
}

# queen adjacency: edge list, ICAR scaling factor (Riebler et al. 2016) and Laplacian eigenvalues
icar_graph = function(ids) {
	a = readRDS(spain("outputs/catalunya/geometry/cat_abs_topology.rds"))$adjacency[ids, ids]
	diag(a) = 0
	q = diag(rowSums(a)) - a
	stopifnot("adjacency graph not connected" = qr(q)$rank == length(ids) - 1)
	e = which(upper.tri(a) & a > 0, arr.ind = TRUE)
	list(n_edges = nrow(e), node1 = unname(e[, 1]), node2 = unname(e[, 2]), icar_scale = exp(mean(log(diag(MASS::ginv(q))))),
		laplacian_eigen = pmax(eigen(q, symmetric = TRUE, only.values = TRUE)$values, 0))
}

# share of days t-10..t-4 with each restriction in force
RESTRICTIONS = c("npi_hospitality", "npi_curfew", "npi_alarm", "npi_hibernation")
restriction_share = function(dates, cols = RESTRICTIONS) {
	x = read_parquet(spain("outputs/catalunya/panel/cat_abs_main_noncollinear.parquet"), col_select = all_of(c("id", "date", cols)))
	daily = x |> group_by(date) |> summarise(across(all_of(cols), \(v) {stopifnot(all(v == v[1])); v[1]}))
	m = as.matrix(daily[match(as.character(rep(dates, each = 7) - rep(10:4, length(dates))), as.character(daily$date)), cols])
	stopifnot(!anyNA(m))
	`rownames<-`(rowsum(m, rep(seq_along(dates), each = 7)) / 7, as.character(dates))
}

# restrictions that vary over the days, without duplicates
restriction_columns = function(m) {
	keep = apply(m, 2, \(v) diff(range(v)) > 0)
	for (k in which(keep)) if (any(vapply(which(keep)[which(keep) < k], \(j) isTRUE(all.equal(m[, j], m[, k])), logical(1)))) keep[k] = FALSE
	m[, keep, drop = FALSE]
}

# public holidays and the days after them
holidays = function(dates) {
	h = read.csv(spain("outputs/catalunya/tables/cat_abs_calendar_daily.csv"))
	is_holiday = \(x) h$holiday[match(as.character(x), h$date)]
	now = is_holiday(dates)
	before = is_holiday(dates - 1)
	stopifnot(is.logical(now), !anyNA(now), !anyNA(before))
	list(holiday_day = as.numeric(now), after_holiday = as.numeric(before & !now))
}

# extra pressures, cached per window: within-ABS trips, third mobility lag, adjacency orders 2-5
extra_inputs = function(d) {
	f = out(sprintf("extra_%d_%s_%s.rds", length(d$ids), min(d$dates), max(d$dates)))
	if (file.exists(f)) return(readRDS(f))
	ids = d$ids
	modelled = as.numeric(d$dates - DAY1) + 1
	days = (min(modelled) - N_LAGS):max(modelled)
	k = read_cases(ids)
	prevalence = t(k$cases / k$pop)
	window = trip_windows(ids, days, max_from = 20)
	lagged = \(m) lapply(seq_len(N_LAGS), \(s) t(m[match(modelled - s, days), , drop = FALSE]))
	by_day = \(g) t(vapply(days, g, numeric(length(ids))))
	a = readRDS(spain("outputs/catalunya/geometry/cat_abs_topology.rds"))$adjacency[ids, ids]
	diag(a) = 0
	order = igraph::distances(igraph::graph_from_adjacency_matrix((a > 0) * 1, mode = "undirected"))
	order_weights = \(o) {m = (order == o) * 1; n = rowSums(m); m / ifelse(n > 0, n, 1)}
	x = list(mobility_pressure3 = lagged(by_day(\(u) drop(prevalence[u, ] %*% window(u, 20, 14)))),
		within = lapply(list(c(6, 0), c(13, 7), c(20, 14)), \(w) lagged(by_day(\(u) prevalence[u, ] * diag(window(u, w[1], w[2]))))),
		neighbour_orders = lapply(2:5, \(o) lagged(prevalence[days, ] %*% order_weights(o))))
	saveRDS(x, f)
	x
}

# switches of stan/spacetime.stan and their defaults
SWITCHES = list(time_effect = 0, vaccination = 0, testing = 0, pooled = 0, negbin = 0, rw_order = 1, gamma_prior = 0.5,
	npi = 0, mob_closure = 0, lag_dist = 0, lag_min = 4, phi_area = 0, holiday = 0, no_within = 0, mob_lags = 2,
	ne_powerlaw = 0, day_effect = 0, placebo = 0, no_mobility = 0, no_neighbour = 0, covariates = 0, save_expected = 0)
BASELINE = list(time_effect = 2, negbin = 1, gamma_prior = 0.2)
# CV_REF: reference of the cross-validation rounds; MODEL: the final model
CV_REF = modifyList(BASELINE, list(holiday = 1, phi_area = 1))
MODEL = modifyList(CV_REF, list(pooled = 1, covariates = 1))

# static covariates: coefficient row (1 endemic, 2 autoregressive) and transformation
COVARIATES = list(income_median = list(row = 1, f = log), pct_foreign = list(row = 1, f = identity),
	built_share = list(row = 1, f = log), hh_size = list(row = 2, f = identity), pct_65plus = list(row = 2, f = identity))
abs_covariates = function(ids) {
	x = read_parquet(spain("outputs/catalunya/panel/cat_abs_main_noncollinear.parquet"), col_select = all_of(c("id", names(COVARIATES)))) |>
		filter(id %in% ids) |> distinct()
	stopifnot(nrow(x) == length(ids), !anyNA(x))
	m = vapply(names(COVARIATES), \(v) COVARIATES[[v]]$f(x[[v]][match(ids, x$id)]), numeric(length(ids)))
	scale(m)
}

# Stan data list of stan/spacetime.stan; holdout: ABS x day 0/1 matrix of held-out cells
spacetime_data = function(d, switches, holdout = NULL) {
	s = modifyList(SWITCHES, switches[intersect(names(switches), names(SWITCHES))])
	stopifnot(s$mob_lags %in% 2:3, s$lag_min %in% seq_len(N_LAGS))
	e = if (s$no_within || s$mob_lags == 3 || s$ne_powerlaw) extra_inputs(d)
	mob = list(d$mobility_pressure1, d$mobility_pressure2, e$mobility_pressure3)[seq_len(s$mob_lags)]
	if (s$no_within) mob = Map(\(m, w) Map(\(a, b) pmax(a - b, 0), m, w), mob, e$within[seq_len(s$mob_lags)])
	ne = c(list(d$neighbour_pressure), if (s$ne_powerlaw) e$neighbour_orders)
	if (s$no_mobility) mob = lapply(mob, \(m) lapply(m, \(x) 0 * x))
	if (s$no_neighbour) ne = lapply(ne, \(m) lapply(m, \(x) 0 * x))
	if (s$placebo) {
		n = nrow(d$cases)
		perm = withr::with_seed(1, {p = seq_len(n); while (any(p == seq_len(n))) p = sample(n); p})
		mob = lapply(mob, \(m) lapply(m, \(x) x[perm, , drop = FALSE]))
	}
	columns = \(x) if (s$lag_dist) do.call(cbind, lapply(x, as.vector)) else matrix(as.vector(Reduce(`+`, Map(`*`, x, LAG_WEIGHTS$week))))
	stack = \(l) aperm(array(unlist(l), c(dim(l[[1]]), length(l))), c(3, 1, 2))
	r = restriction_share(d$dates)
	npi_index = if (s$npi) match(colnames(restriction_columns(r)), colnames(r)) else integer(0)
	if (is.null(holdout)) holdout = matrix(0L, nrow(d$cases), ncol(d$cases))
	c(list(n_areas = nrow(d$cases), n_days = ncol(d$cases), cases = d$cases, n_lags = N_LAGS, n_orders = length(ne),
			autoregressive_lags = columns(d$lagged_cases), neighbour_lags = stack(lapply(ne, columns)),
			mobility_lags = stack(lapply(mob, columns)),
			pop = as.numeric(d$pop), weekday = as.integer(format(d$dates, "%u")), prior_scale = c(0.01, 2, 500, rep(2, s$mob_lags))),
		lapply(s[c("time_effect", "vaccination", "testing", "pooled", "negbin", "phi_area", "rw_order", "npi", "mob_closure",
			"lag_dist", "lag_min", "mob_lags", "holiday", "day_effect", "save_expected")], as.integer),
		if (s$testing == 2) {stopifnot("testing = 2 uses the one-week lag weights" = s$lag_dist == 0)
			list(log_tests_lag = log_tests(d$ids, d$dates, lagged = TRUE))} else list(log_tests_lag = log_tests(d$ids, d$dates)),
		with(list(g = if (s$day_effect == 2) as.integer(factor(abs_province(d$ids))) else rep(1L, nrow(d$cases))),
			list(n_day_groups = max(g), day_group = g)),
		if (s$covariates) list(n_cov = length(COVARIATES), covariate = unname(abs_covariates(d$ids)),
			cov_row = array(vapply(COVARIATES, `[[`, 0, "row"))) else list(n_cov = 0L, covariate = matrix(0, nrow(d$cases), 0), cov_row = integer(0)),
		list(gamma_prior = s$gamma_prior, n_blocks = ceiling(ncol(d$cases) / 7), day_block = ceiling(seq_len(ncol(d$cases)) / 7)),
		week_knots(ncol(d$cases)), icar_graph(d$ids), holidays(d$dates),
		list(vaccinated = vaccinated_share(d$ids, d$dates), log_tests = log_tests(d$ids, d$dates),
			n_restrictions = ncol(r), restrictions = unname(r), n_npi = length(npi_index), npi_index = array(npi_index),
			closure = unname(r[, "npi_hospitality"]),
			holdout = array(as.integer(holdout), dim(d$cases))))
}

# variant name of a fit from its switches
spacetime_variant = function(switches) with(modifyList(SWITCHES, switches[intersect(names(switches), names(SWITCHES))]),
	paste0("time", time_effect, if (rw_order == 2) "_rw2", if (gamma_prior != 0.5) paste0("_gsd", gamma_prior),
		if (npi) "_npi", if (mob_closure) "_mobclosure", if (negbin) "_negbin", c("", "_pooled", "_adaptive")[pooled + 1],
		if (testing) c("_testing", "_testing2")[testing], if (vaccination) "_vaccination", if (lag_dist) paste0("_lagnb", lag_min),
		if (phi_area) "_phiarea", if (holiday) "_holiday", if (no_within) "_nowithin", if (mob_lags == 3) "_mob3",
		if (ne_powerlaw) "_powerlaw", if (day_effect) c("_dayeffect", "_dayprov")[day_effect], if (placebo) "_placebo",
		if (no_mobility) "_nomob", if (no_neighbour) "_none", if (covariates) "_cov", if (save_expected) "_ppc"))

# cross-validation variant sets and the fold scheme of each
CV_VARIANTS = c(list(BASELINE), lapply(list(list(lag_dist = 1, lag_min = 4), list(lag_dist = 1, lag_min = 7),
	list(phi_area = 1), list(holiday = 1), list(testing = 1), list(no_within = 1), list(mob_lags = 3), list(ne_powerlaw = 1),
	list(day_effect = 1), list(placebo = 1)),
	\(x) modifyList(BASELINE, x)))
names(CV_VARIANTS) = vapply(CV_VARIANTS, spacetime_variant, character(1))
CV_SETS = list(refinements = CV_VARIANTS,
	heterogeneity = lapply(list(list(), list(pooled = 1), list(pooled = 2), list(testing = 2)), \(x) modifyList(CV_REF, x)))
CV_SETS$coupling = lapply(list(list(), list(ne_powerlaw = 1, no_mobility = 1), list(no_mobility = 1), list(no_neighbour = 1)),
	\(x) modifyList(CV_REF, x))
CV_SETS$daysplit = lapply(list(list(), list(day_effect = 1)), \(x) modifyList(CV_REF, x))
for (k in names(CV_SETS)[-1]) names(CV_SETS[[k]]) = vapply(CV_SETS[[k]], spacetime_variant, character(1))
CV_SCHEME = c(refinements = "week", heterogeneity = "week", coupling = "week", daysplit = "split")

# fold file and holdout matrix of a cross-validation scheme
cv_fold_file = function(window_key, variant, k, K, scheme = "week")
	out(file.path("cv", paste0(window_key, "_", variant, if (scheme == "week") "_K" else "_S", K, "_fold", k, ".rds")))
cv_holdout = function(scheme, n_areas, n_days, k, K) if (scheme == "week") week_fold(n_areas, n_days, k, K) else split_fold(n_areas, n_days, k, K)

# pointwise log predictive density of held-out cells from a draws x cells log-likelihood matrix
lpd = function(ll) apply(ll, 2, \(v) max(v) + log(mean(exp(v - max(v)))))

# week folds: week w held out in fold (w - 1) %% K + 1
week_fold = function(n_areas, n_days, k, K) {
	held = ((ceiling(seq_len(n_days) / 7) - 1) %% K + 1) == k
	matrix(rep(as.integer(held), each = n_areas), n_areas, n_days)
}

# within-day split: the ABS of each day split at random into K folds
split_fold = function(n_areas, n_days, k, K) {
	g = withr::with_seed(1, vapply(seq_len(n_days), \(t) sample(rep_len(seq_len(K), n_areas)), integer(n_areas)))
	(g == k) * 1L
}

# convergence, PSIS-LOO and trips per new infection of a fit
report = function(fit, loo_result = NULL, pars = c("endemic", "autoregressive", "neighbour", "mobility")) {
	s = posterior::summarise_draws(posterior::as_draws_array(as.array(fit, pars = pars)), "rhat", "ess_bulk", "ess_tail")
	cat(sprintf("max R-hat %.3f, min bulk ESS %.0f, min tail ESS %.0f, divergences %d, max treedepth hits %d\n",
		max(s$rhat, na.rm = TRUE), min(s$ess_bulk, na.rm = TRUE), min(s$ess_tail, na.rm = TRUE),
		rstan::get_num_divergent(fit), rstan::get_num_max_treedepth(fit)))
	print(if (is.null(loo_result)) loo::extract_log_lik(fit) |> loo::loo() else loo_result)
	cat("trips per new infection (1 / infections_per_trip):\n")
	print(quantile(1 / rstan::extract(fit)$infections_per_trip, c(0.025, 0.5, 0.975)))
}
