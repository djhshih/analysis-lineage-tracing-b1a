library(io)
library(magrittr)
library(danitizer)
library(dplyr)

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


# join tables

hec.d <- with(hec, data.frame(
	mouse,
	bm_fl = select(bm_fl, -id),
	neonatal_spleen = select(neonatal_spleen, -id),
	spleen = select(spleen, -id),
	thymus = select(thymus, -id),
	peritoneal = select(peritoneal_cavity, -id)
));

hsc.d <- with(hsc, data.frame(
	mouse,
	bm_fl = select(bm_fl, -id),
	neonatal_spleen = select(neonatal_spleen, -id),
	spleen = select(spleen, -id),
	thymus = select(thymus, -id),
	peritoneal = select(peritoneal_cavity, -id)
));

qwrite(hec.d, "data/hec.rds");
qwrite(hsc.d, "data/hsc.rds");

