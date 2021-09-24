library(ggplot2)
library(reshape2)
library(io)
library(dplyr)
library(ggsci)

out.fname <- filename("lineage-sim");
pdf.fname <- insert(out.fname, ext="pdf");


# --- Parameters

# number of cell types
J <- 6;

# number of samples
N <- 2;

# number of analysis time points
T <- 500;

# unlabeled cells and labeled cells
L <- 2;     

# time offset in order to account for unobservable period
s0 <- 5;


# labelling times
ss <- c(10, 35, 60, 85, 110) + s0;
#ss <- seq(7, 12) + s;

# shift time s.t. birth occurs t = 0
t.shift <- -200;
# pre-natal time is compressed
t.neg.scale <- 10;


names(ss) <- c("E7.5", "E8.5", "E9.5", "E10.5", "E11.5");

stopifnot(s0 > 1 && s0 < T)


# --- Random variables

set.seed(6)

compart.names <- c("unlabeled", "labeled");

cell.types <- c("HEC", "fetal HSC", "adult HSC", "MPP", "B-1 pro", "B-1");
stopifnot(length(cell.types) == J);

cell.types.obs <- c("HEC", "HSC", "MPP", "B-1 pro", "B-1");

# unknown initial common progenitor size
n0 <- 1e6;

# unknown labeling efficiency
#lambda <- runif(N);
lambda <- seq(0.2, 0.8, length.out=N);

# differentiation rates
# rows: from;  columns: to

A0 <- matrix(
	c(
		0, 0.01, 0, 0, 0, 0,
		0, 0, 0.01, 0.2, 0, 0,
		0, 0, 0, 0.1, 0, 0,
		0, 0, 0, 0, 0.2, 0,
		0, 0, 0, 0, 0, 0.2,
		0, 0, 0, 0, 0, 0
	),
	nrow=J, byrow=TRUE
);

A1 <- A0;
A1[1, 4] <- 0.05;

A2 <- A1;
A2[1, 5] <- 0.05;

As <- list(H0=A0, H1=A1, H2=A2);

lapply(As,
	function(A) {
		# ensure that the diagonals are 0
		stopifnot(diag(A) == rep(0, nrow(A)))
	}
)

# net proliferation rate
# assume that beta is invariant across time
#beta <- runif(J, -0.5, 0.5);  # arbitrary upperbound
beta <- c(-0.01, -0.01, 0.2, 0.5, 0.1, 0);
#beta <- c(0, 0.1, 0.3, 0.05, 0);
#beta <- c(0, 0.01, 0.04, 1, 0);

# relative limiting capacities
nl <- 1e5;
kappa <- c(10, 10, 1, 9, 52, 100);

params0 <- list(
	J = J, N = N, T = T, L = 2,
	n0 = n0, nl = nl,
	lambda = lambda, beta = beta,
	kappa = kappa
);


# --- Functions

cell_factor <- function(j, cell.types) {
	factor(j, levels=1:length(cell.types), labels=cell.types)
}

# net change = net proliferation + differentiation into - differentiation out of
update_ntl_exp <- function(n.tm1.l, n.tm1, params) {
	pmax(
		0,
		n.tm1.l + t(t(n.tm1.l) * beta) + n.tm1.l %*% params$A  - t(t(n.tm1.l) * params$alpha.out
	))
}

# @param n.tm1.l   compartment size for selected label component at t - 1
# @param n.tm1     total size across both labeled and unlabeled compartments
#                  at t - 1
update_ntl_logistic <- function(n.tm1.l, n.tm1, params) {
	pmax(
		0,
		n.tm1.l + t(t(n.tm1.l) * beta * pmax(0, 1 - t(n.tm1) / params$kappa.p)) + n.tm1.l %*% params$A  - t(t(n.tm1.l) * params$alpha.out
	))
}


# simulate cell differentiation trajectory
simulate_trajectory <- function(params, update_ntl) {

	# Setup initial conditions

	params <- within(params,
		{
			alpha.out <- rowSums(A);
			kappa.p <- nl * kappa;
		}
	);

	n <- with(params, array(0, dim = c(n=N, j=J, t=T, l=L)));

	# only unlabeled common progenitor is present at t = 0
	n[, 1, 1, 1] <- params$n0;


	# Iterate until labeling

	for (t in 2:params$s) {
		# all cells are unlabeled right now,
		# so we only update the unlabeled compartment
		n[, , t, 1] <- update_ntl(
			n[, , t - 1, 1],
			apply(n[, , t - 1, ], c(1, 2), sum),
			params
		);
	}


	# Apply labeling

	# label the first common progenitor
	m <- n[, 1, params$s, 1];

	# unlabeled
	n[, 1, params$s, 1] <- with(params, m * (1 - lambda));

	# labeled
	n[, 1, params$s, 2] <- with(params, m * lambda);

	# all other cells remain unlabeled


	# Iterate until end of time

	for (t in (params$s+1):params$T) {
		for (l in 1:params$L) {
			n[, , t, l] <- update_ntl(
				n[, , t - 1, l],
				apply(n[, , t - 1, ], c(1, 2), sum),
				params
			);
		}
	}

	# Calculate observed variable

	n
}

# --- 

params <- within(params0, {
	ss <- ss;
	As <- As;
});

qwrite(params, insert(out.fname, tag="params", ext="rds"));

nss <- lapply(As,
	function(A) {
		ns <- lapply(ss,
			function(s) {
				params <- within(params0, {
					s <- s;
					A <- A;
				});
				simulate_trajectory(params, update_ntl_logistic)
			}
		);
		
		ns
	}
);

calculate_fractions <- function(n) {
	n[, , , dim(n)[4]] / pmax(1, apply(n, c(1, 2, 3), sum));
}

fss <- lapply(nss, function(ns) lapply(ns, calculate_fractions));


# --- Generate plots

reshape_fractions <- function(fs, cell.types) {
	f.d <- melt(fs, varnames=names(dim(fs[[1]])));

	f.d$s <- factor(f.d$L1, levels=names(ss));
	f.d$L1 <- NULL;

	f.d$j <- cell_factor(f.d$j, cell.types);

	# scale time
	f.d$t <- f.d$t + t.shift;
	idx <- f.d$t < 0;
	f.d$t[idx] <- f.d$t[idx] / t.neg.scale;

	f.d
}

# normalize fractions against the jth component,
# and average over samples
normalize_fractions <- function(f, j) {
	fn <- f / f[, rep(j, dim(f)[2]), ];
	apply(fn, c(2, 3), mean)
}

# combine two populations together
combine_populations <- function(n, from, to) {
	n[, to, , ] <- n[, to, , ] + n[, from, , ];
	n[, -from, , ]
}


graphics.off()

for (h in names(fss)) {

	fs <- fss[[h]];
	ns <- nss[[h]];

	n.d <- reshape_fractions(ns, cell.types);
	n.d$l <- factor(n.d$l, labels=compart.names);

	g <- ggplot(filter(n.d, n == params$N), aes(x = t, y = value, colour = j, linetype=l)) +
		theme_classic() +
		geom_line() + facet_grid(j ~ s, scales="free_y") +
		guides(colour = "none") +
		scale_colour_npg() +
		theme(legend.title=element_blank()) +
		xlab("time") + ylab("compartment size")
	qdraw(g, width = 7, file = insert(pdf.fname, c(tolower(h), "latent", "label-n")))

	f.d <- reshape_fractions(fs, cell.types);

	g <- ggplot(filter(f.d, n == 1), aes(x = t, y = value, colour = j)) +
		theme_classic() +
		geom_line(linetype=2) + facet_grid(s ~ j) +
		guides(colour = "none") +
		scale_colour_npg() +
		xlab("analysis time") + ylab("% labelled")
	qdraw(g, width = 6, file = insert(pdf.fname, c(tolower(h), "latent", "label-prop")))

	if (FALSE) {
		# normalize against the common progenitor
		fns <- lapply(fs, function(f) normalize_fractions(f, 1));

		fn.d <- reshape_fractions(fns, cell.types);

		g <- ggplot(fn.d, aes(x = t, y = value, colour = j)) +
			theme_classic() +
			geom_hline(yintercept = 1, linetype=3, colour="grey30") +
			geom_line(linetype=2) + 
			facet_grid(s ~ j, scales="free_y") +
			guides(colour = "none") +
			scale_colour_npg() +
			xlab("analysis time") + ylab("label ratio vs. HEC")
		qdraw(g, file = insert(pdf.fname, c(tolower(h), "latent", "label-ratio")))
	}

	# combine second and third populations because
	# they are not distinguishable during observation
	ns.obs <- lapply(ns, combine_populations, from=3, to=2);

	fs.obs <- lapply(ns.obs, calculate_fractions);

	# normalize against the second progenitor
	fns.obs <- lapply(fs.obs, function(f) normalize_fractions(f, 2));

	fn.obs.d <- reshape_fractions(fns.obs, cell.types.obs);

	g <- ggplot(fn.obs.d[fn.obs.d$j != "HEC", ], aes(x = t, y = value, colour = j)) +
		theme_classic() +
		geom_hline(yintercept = 1, linetype=3, colour="grey30") +
		geom_line() + facet_grid(s ~ j, scales="free_y") +
		guides(colour = "none") +
		scale_colour_npg() +
		coord_cartesian(ylim=c(0, 2)) +
		xlab("analysis time") + ylab("label ratio vs. HSC")
	qdraw(g, file = insert(pdf.fname, c(tolower(h), "observed", "label-ratio")))

	fn.obs.f.d <- filter(fn.obs.d, t == max(fn.obs.d$t), j != "HEC");

	g <- ggplot(fn.obs.f.d, aes(x = j, fill = j, y = value)) +
		theme_classic() +
		geom_col() + facet_wrap(~ s, ncol=1) +
		geom_hline(yintercept = 1, linetype=3, colour="grey30") +
		scale_fill_npg() +
		guides(fill = "none") +
		theme(
			axis.text.x = element_text(angle = 40, hjust=1),
			strip.background = element_blank()
		) +
		xlab("") + ylab("label ratio vs. HSC")
	qdraw(g, width = 1.5, file = insert(pdf.fname, c(tolower(h), "observed", "label-ratio", "end")))

}

