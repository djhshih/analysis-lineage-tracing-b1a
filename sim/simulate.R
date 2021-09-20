library(ggplot2)
library(reshape2)
library(io)
library(dplyr)

out.fname <- filename("lineage-sim");
pdf.fname <- insert(out.fname, ext="pdf");


# --- Parameters

# number of cell types
J <- 6;

# number of samples
N <- 4;

# number of analysis time points
T <- 300;

# unlabeled cells and labeled cells
L <- 2;     

# time offset in order to account for unobservable period
s <- 5;

# labelling times
ss <- c(10, 20, 30, 40, 50, 60) + s;
#ss <- seq(7, 12) + s;

stopifnot(s > 1 && s < T)


# --- Random variables

set.seed(6)

cell.types <- c("HEC", "fetal HSC", "adult HSC", "MPP", "B-1 pro", "B-1");
stopifnot(length(cell.types) == J);

cell.types.obs <- c("HEC", "HSC", "MPP", "B-1 pro", "B-1");

# unknown initial common progenitor size
n0 <- 1e11;

# unknown labeling efficiency
#lambda <- runif(N);
lambda <- seq(0.1, 0.9, length.out=N);

# differentiation rates
# rows: from;  columns: to

A0 <- matrix(
	c(
		0, 0.05, 0, 0, 0, 0,
		0, 0, 0.1, 0.2, 0, 0,
		0, 0, 0, 0.05, 0, 0,
		0, 0, 0, 0, 0.2, 0,
		0, 0, 0, 0, 0, 0.2,
		0, 0, 0, 0, 0, 0
	),
	nrow=J, byrow=TRUE
);

A1 <- A0;
A1[1, 4] <- 0.1;

A2 <- A1;
A2[1, 5] <- 0.2;


A <- A0;
#A <- A1;
#A <- A2;

# ensure that the diagonals are 0
stopifnot(diag(A) == rep(0, nrow(A)))

# net proliferation rate
# assume that beta is invariant across time
#beta <- runif(J, -0.5, 0.5);  # arbitrary upperbound
beta <- c(0, 0, 0.05, 0.5, 0.1, 0);
#beta <- c(0, 0.1, 0.3, 0.05, 0);
#beta <- c(0, 0.01, 0.04, 1, 0);

# relative limiting capacities
nl <- 1e9;
kappa <- c(1, 1, 0.1, 9, 52);

params <- list(
	J = J, N = N, T = T, L = 2, s = s,
	n0 = n0, nl = nl,
	lambda = lambda, A = A, beta = beta,
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
		n.tm1.l + t(t(n.tm1.l) * beta) + n.tm1.l %*% A  - t(t(n.tm1.l) * params$alpha.out
	))
}

# @param n.tm1.l   compartment size for selected label component at t - 1
# @param n.tm1     total size across both labeled and unlabeled compartments
#                  at t - 1
update_ntl_logistic <- function(n.tm1.l, n.tm1, params) {
	pmax(
		0,
		n.tm1.l + t(t(n.tm1.l) * beta * (1 - t(n.tm1) / params$kappa.p)) + n.tm1.l %*% A  - t(t(n.tm1.l) * params$alpha.out
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


ns <- lapply(ss,
	function(s) {
		params2 <- params;
		params2$s <- s;
		simulate_trajectory(params2, update_ntl_logistic)
		#simulate_trajectory(params2, update_ntl_exp)
	}
);
names(ns) <- ss;

calculate_fractions <- function(n) {
	n[, , , dim(n)[4]] / pmax(1, apply(n, c(1, 2, 3), sum));
}

fs <- lapply(ns, calculate_fractions);

# --- Generate plots

graphics.off()


reshape_fractions <- function(fs, cell.types) {
	f.d <- melt(fs, varnames=names(dim(fs[[1]])));
	f.d$s <- as.integer(f.d$L1); f.d$L1 <- NULL;
	f.d$j <- cell_factor(f.d$j, cell.types);

	f.d
}

f.d <- reshape_fractions(fs, cell.types);

g <- ggplot(filter(f.d, n == params$N), aes(x = t, y = value, colour = j)) +
	theme_classic() +
	geom_line(linetype=2) + facet_grid(s ~ j) +
	guides(colour = "none") +
	xlab("analysis time") + ylab("% labelled")
qdraw(g, width = 6, file = insert(pdf.fname, c("latent", "label-prop")))

# normalize fractions against the jth component,
# and average over samples
normalize_fractions <- function(f, j) {
	fn <- f / f[, rep(j, dim(f)[2]), ];
	apply(fn, c(2, 3), mean)
}

# normalize against the common progenitor
fns <- lapply(fs, function(f) normalize_fractions(f, 1));

fn.d <- reshape_fractions(fns, cell.types);

g <- ggplot(fn.d, aes(x = t, y = value, colour = j)) +
	theme_classic() +
	geom_hline(yintercept = 1, linetype=3, colour="grey30") +
	geom_line(linetype=2) + 
	facet_grid(s ~ j, scales="free_y") +
	guides(colour = "none") +
	xlab("analysis time") + ylab("label ratio vs. HEC")
#qdraw(g, file = insert(pdf.fname, c("latent", "label-ratio")))

combine_populations <- function(n, from, to) {
	n[, to, , ] <- n[, to, , ] + n[, from, , ];
	n[, -from, , ]
}


# combine second and third populations because
# they are not distinguishable during observation
ns.obs <- lapply(ns, combine_populations, from=3, to=2);

plot(ns[[1]][1, 1, , 1], type="l")
plot(ns[[1]][1, 2, , 1], type="l")
plot(ns[[1]][1, 3, , 1], type="l")
plot(ns[[1]][1, 4, , 1], type="l")
plot(ns[[1]][1, 5, , 1], type="l")
plot(ns[[1]][1, 6, , 1], type="l")

plot(ns[[1]][1, 2, , 2] + ns[[1]][1, 3, , 2], type="l")
#plot(ns.obs[[1]][1, 2, , 1], type="l")

fs.obs <- lapply(ns.obs, calculate_fractions);

# normalize against the second progenitor
fns.obs <- lapply(fs.obs, function(f) normalize_fractions(f, 2));

fn.obs.d <- reshape_fractions(fns.obs, cell.types.obs);

g <- ggplot(fn.obs.d[fn.obs.d$j != "HEC", ], aes(x = t, y = value, colour = factor(j))) +
	theme_classic() +
	geom_hline(yintercept = 1, linetype=3, colour="grey30") +
	geom_line() + facet_grid(s ~ j, scales="free_y") +
	guides(colour = "none") +
	coord_cartesian(ylim=c(0, 2)) +
	xlab("analysis time") + ylab("label ratio vs. HSC")
qdraw(g, file = insert(pdf.fname, c("observed", "label-ratio")))

qwrite(params, insert(out.fname, tag="params", ext="rds"));

