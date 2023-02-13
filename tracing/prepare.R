library(io)
library(magrittr)
library(danitizer)
library(dplyr)
library(tidyr)

hec <- qread("data/hec", na.strings="N/A");
hsc <- qread("data/hsc", na.strings="N/A");

# mouse development:
# P0 is time of birth, which occurs 19 days post conception (dpc)
dpc.birth <- 19;

hec <- lapply(hec, normalize_df_field_names);
hsc <- lapply(hsc, normalize_df_field_names);

fix_group <- function(x) {
	# E and P are uppercase
	x <- toupper(x);
	# unannoated times are E
	x <- gsub("^([0-9.]+)$", "E\\1", x);
	x
}

fix_t_injection <- function(x) {
	# remove p and replace with dpc
	idx <- grep("p|P", x);
	x[idx] <- gsub("p|P", "", x[idx]);
	x <- as.numeric(x);
	x[idx] <- x[idx] + dpc.birth;
	x
}	


# fix group names
hec$mouse$group <- fix_group(hec$mouse$group);
hsc$mouse$group <- fix_group(hsc$mouse$group);

# fix time of injection
hsc$mouse$t_injection <- fix_t_injection(hsc$mouse$t_injection);

# shift reference time point to birth (t = 0)
# t_analysis already uses this reference
hsc$mouse$t_injection <- hsc$mouse$t_injection - dpc.birth;
hec$mouse$t_injection <- hec$mouse$t_injection - dpc.birth;

# select key fields
hec$mouse <- select(hec$mouse, id, t_injection, t_analysis, sex, dose_tam, group)
hsc$mouse <- select(hsc$mouse, id, t_injection, t_analysis, sex, dose_tam, group)

remove_id <- function(xs) {
	with(xs,
		list(
			mouse = mouse,
			bm_fl = select(bm_fl, -id),
			neonatal_spleen = select(neonatal_spleen, -id),
			spleen = select(spleen, -id),
			thymus = select(thymus, -id),
			peritoneal_cavity = select(peritoneal_cavity, -id)
		)
	)
}

devel_stage_factor <- function(groups) {
	levels <- unique(groups);
	stage <- ifelse(grepl("E", levels), 0, 1);
	days <- as.numeric(gsub("(E|P)(\\d+).*", "\\2", levels));
	days.end <- suppressWarnings(as.numeric(gsub(".*-(\\d+)", "\\1", levels)));
	idx <- order(stage, days, days.end);
	factor(groups, levels[idx])
}

fix_devel_stages <- function(xs) {
	# remove problematic stages
	idx <- xs$mouse$group != "E0";
	ys <- lapply(xs, function(x) x[idx, ]);

	ys$mouse$group <- devel_stage_factor(ys$mouse$group);

	ys
}

hec <- fix_devel_stages(hec);
hsc <- fix_devel_stages(hsc);

qwrite(hec, "data/hec.rds");
qwrite(hsc, "data/hsc.rds");

# ---

# combine tables column-wise
combine_tables <- function(xs) {
	with(remove_id(xs), data.frame(
		mouse,
		bm_fl = bm_fl,
		neonatal_spleen = neonatal_spleen,
		spleen = spleen,
		thymus = thymus,
		peritoneal = peritoneal_cavity
	))
}


hec.w <- combine_tables(hec);
hsc.w <- combine_tables(hsc);


qwrite(hec.w, "data/hec_wide.rds");
qwrite(hsc.w, "data/hsc_wide.rds");


# stack tables
# main: main table name
stack_tables <- function(xs, main, values_to="count") {
	others <- setdiff(names(xs), main);
	names(others) <- others;
	ys <- do.call(rbind, lapply(others,
		function(other) {
			pivot_longer(xs[[other]], !id, names_to = "cell", values_to=values_to, values_drop_na = TRUE) %>%
				mutate(site = other)
		}
	));
	
	zs <- left_join(xs[[main]], ys, by="id");
	zs[!is.na(zs[[values_to]]), ]
}

hec.l <- stack_tables(hec, "mouse");
hsc.l <- stack_tables(hsc, "mouse");

qwrite(hec.l, "data/hec_long.rds");
qwrite(hsc.l, "data/hsc_long.rds");



normalize_wrt_hsc <- function(xs, min.val=1) {
	xs$bm_fl$hsc[which(xs$bm_fl$hsc < min.val)] <- NA;

	c(
		xs["mouse"],
		lapply(xs[! names(xs) %in% "mouse"],
			function(x) {
				as.data.frame(lapply(x, function(z) {
					if (is.numeric(z)) {
						z / xs$bm_fl$hsc
					} else {
						z
					}
				}))
			}
		)
	)
}

hec.n <- normalize_wrt_hsc(hec);
hsc.n <- normalize_wrt_hsc(hsc);

hec.wn <- combine_tables(hec.n);
hsc.wn <- combine_tables(hsc.n);

qwrite(hec.wn, "data/normalized/hec_wide_normalized.rds");
qwrite(hsc.wn, "data/normalized/hsc_wide_normalized.rds");


hec.ln <- stack_tables(hec.n, main="mouse", values_to="ratio");
hsc.ln <- stack_tables(hsc.n, main="mouse", values_to="ratio");

qwrite(hec.ln, "data/normalized/hec_long_normalized.rds");
qwrite(hsc.ln, "data/normalized/hsc_long_normalized.rds");

