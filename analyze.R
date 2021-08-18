library(io)
library(ggplot2)
library(dplyr)

out.fname <- filename("lineage", path="plot", date=NA);
pdf.fname <- insert(out.fname, ext="pdf");


hec.w <- qread("data/hec_wide.rds");
hsc.w <- qread("data/hsc_wide.rds");


ggplot(hec.w, aes(t_injection, bm_fl.hsc)) + 
	theme_classic() +
	geom_point()

ggplot(hec.w, aes(t_analysis, bm_fl.hsc)) +
	theme_classic() +
	geom_point()

ggplot(hec.w, aes(t_injection, t_analysis, size=bm_fl.hsc)) +
	theme_classic() +
	geom_point()


# ---

hec.wn <- qread("data/normalized/hec_wide_normalized.rds");
hsc.wn <- qread("data/normalized/hsc_wide_normalized.rds");

# @param d  table in wide format
# compartment: location.cell (e.g. bm_fl.hsc)
ratio_trace_plot <- function(d, cell) {
	d$ratio <- d[[cell]];
	g <- ggplot(d, aes(t_analysis, ratio)) +
		theme_classic() +
		geom_hline(yintercept=1, colour="grey60") +
		stat_smooth(method="glm") +
		geom_point() +
		ggtitle(cell)

	if (length(unique(d$group)) > 1) {
		g <- g + facet_grid(group ~ ., scales="free_y");
	}

	g
}

ratio_trace_plot(hec.wn, "bm_fl.hsc")
ratio_trace_plot(hec.wn, "bm_fl.st_hsc")
ratio_trace_plot(hec.wn, "bm_fl.mpp")
ratio_trace_plot(hec.wn, "bm_fl.pro_b")
ratio_trace_plot(hec.wn, "bm_fl.pre_b")

grep("b1a", names(hec.wn), value=TRUE)
ratio_trace_plot(hec.wn, "peritoneal.b1a")
ratio_trace_plot(hec.wn, "spleen.b1a")

# examine E7.5
hec.wn.e7.5 <- filter(hec.wn, group == "E7.5");
ratio_trace_plot(hec.wn.e7.5, "bm_fl.hsc")
ratio_trace_plot(hec.wn.e7.5, "bm_fl.mpp")
ratio_trace_plot(hec.wn.e7.5, "bm_fl.pro_b")
ratio_trace_plot(hec.wn.e7.5, "peritoneal.b1a")

# re-analyze peritoneal.b1a after removing the lone outlier
outlier <- filter(hec.wn, group == "E7.5", peritoneal.b1a > 5)$id;
filter(hec.wn, id %in% outlier);

qdraw(
	ratio_trace_plot(hec.wn.e7.5, "peritoneal.b1a"),
	file = insert(pdf.fname, c("hec", "ratio", "grid", "peritoneal", "b1a", "e7-5"))
);

qdraw(
	ratio_trace_plot(filter(hec.wn.e7.5, ! id %in% outlier), "peritoneal.b1a"),
	file = insert(pdf.fname, c("hec", "ratio", "grid", "peritoneal", "b1a", "e7-5", "outlier-rm"))
);


ratio_trace_plot(hsc.wn, "bm_fl.hsc")
ratio_trace_plot(hsc.wn, "bm_fl.st_hsc")
ratio_trace_plot(hsc.wn, "bm_fl.mpp")
ratio_trace_plot(hsc.wn, "bm_fl.pro_b")
ratio_trace_plot(hsc.wn, "bm_fl.pre_b")


# ---

hec.l <- qread("data/hec_long.rds");
hsc.l <- qread("data/hsc_long.rds");

# @param d  table in long format
count_trace_plots <- function(d) {
	ggplot(d, aes(t_analysis, count, colour=cell)) +
		theme_classic() +
		facet_grid(group ~ cell) +
		ylim(0, 100) +
		guides(colour="none") +
		stat_smooth(method="glm") +
		geom_point() +
		scale_x_continuous(n.breaks=3) +
		theme(panel.spacing = unit(0.7, "lines"))
}

qdraw(
	count_trace_plots(filter(hec.l, site == "bm_fl")),
	width = 15, height = 5,
	file = insert(pdf.fname, c("hec", "count", "grid", "bm_fl"))
)

qdraw(
	count_trace_plots(filter(hec.l, site == "spleen")),
	width = 9, height = 5,
	file = insert(pdf.fname, c("hec", "count", "grid", "spleen"))
)

qdraw(
	count_trace_plots(filter(hec.l, site == "thymus")),
	width = 6.5, height = 5,
	file = insert(pdf.fname, c("hec", "count", "grid", "thymus"))
)

qdraw(
	count_trace_plots(filter(hec.l, site == "peritoneal_cavity")),
	width = 6.5, height = 5,
	file = insert(pdf.fname, c("hec", "count", "grid", "peritoneal"))
)


qdraw(
	count_trace_plots(filter(hsc.l, site == "bm_fl")),
	width = 15, height = 5,
	file = insert(pdf.fname, c("hsc", "count", "grid", "bm_fl"))
)

qdraw(
	count_trace_plots(filter(hsc.l, site == "spleen")),
	width = 9, height = 5,
	file = insert(pdf.fname, c("hsc", "count", "grid", "spleen"))
)

qdraw(
	count_trace_plots(filter(hsc.l, site == "thymus")),
	width = 6.5, height = 5,
	file = insert(pdf.fname, c("hsc", "count", "grid", "thymus"))
)

qdraw(
	count_trace_plots(filter(hsc.l, site == "peritoneal_cavity")),
	width = 6.5, height = 5,
	file = insert(pdf.fname, c("hsc", "count", "grid", "peritonael"))
)


# ---

hec.ln <- qread("data/normalized/hec_long_normalized.rds");
hsc.ln <- qread("data/normalized/hsc_long_normalized.rds");

unique(hec.ln$site)

# @param d  table in long format
ratio_trace_plots <- function(d, ylim=NULL) {
	ggplot(d, aes(t_analysis, ratio, colour=cell)) +
		theme_classic() +
		geom_hline(yintercept=1, colour="grey60") +
		guides(colour="none") +
		coord_cartesian(ylim = ylim) +
		stat_smooth(method="glm") +
		geom_point() +
		scale_x_continuous(n.breaks=3) +
		theme(panel.spacing = unit(0.7, "lines"))
}

qdraw(
	ratio_trace_plots(filter(hec.ln, site == "bm_fl"), ylim=c(0, 4)),
	width = 15, height = 5,
	file = insert(pdf.fname, c("hec", "ratio", "grid", "bm_fl"))
)

qdraw(
	ratio_trace_plots(filter(hec.ln, site == "spleen")),
	width = 9, height = 5,
	file = insert(pdf.fname, c("hec", "ratio", "grid", "spleen"))
)

qdraw(
	ratio_trace_plots(filter(hec.ln, site == "thymus")),
	width = 6.5, height = 5,
	file = insert(pdf.fname, c("hec", "ratio", "grid", "thymus"))
)

qdraw(
	ratio_trace_plots(filter(hec.ln, site == "peritoneal_cavity"), ylim=c(0, 8)),
	width = 6.5, height = 5,
	file = insert(pdf.fname, c("hec", "ratio", "grid", "peritoneal"))
)

ratio_trace_plots(filter(hec.ln, site == "neonatal_spleen"));


qdraw(
	ratio_trace_plots(filter(hsc.ln, site == "bm_fl")),
	width = 15, height = 5,
	file = insert(pdf.fname, c("hsc", "ratio", "grid", "bm_fl"))
)

qdraw(
	ratio_trace_plots(filter(hsc.ln, site == "spleen")),
	width = 9, height = 5,
	file = insert(pdf.fname, c("hsc", "ratio", "grid", "spleen"))
)

qdraw(
	ratio_trace_plots(filter(hsc.ln, site == "thymus")),
	width = 6.5, height = 5,
	file = insert(pdf.fname, c("hsc", "ratio", "grid", "thymus"))
)

qdraw(
	ratio_trace_plots(filter(hsc.ln, site == "peritoneal_cavity")),
	width = 6.5, height = 5,
	file = insert(pdf.fname, c("hsc", "ratio", "grid", "peritoneal"))
)

ratio_trace_plots(filter(hsc.ln, site == "neonatal_spleen"));


