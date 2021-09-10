library(ggplot2)
library(reshape2)
library(io)

out.fname <- filename("lineage-sim");
pdf.fname <- insert(out.fname, ext="pdf");


# --- Parameters

# number of cell types
J <- 5;

# number of samples
N <- 9;

# number of analysis time points
T <- 300;

# unlabeled cells and labeled cells
L <- 2;     

# labeling time point 
s <- 10;
#s <- 150;

stopifnot(s > 1 && s < T)


# --- Random variables

set.seed(6)

# unknown initial common progenitor size
n0 <- 1e3;

# unknown labeling efficiency
#lambda <- runif(N);
lambda <- seq(0.1, 0.9, by=0.1);

# differentiation rate
# rows: from;  columns: to
A0 <- matrix(c(0, 0.05, 0.00, 0, 0,  0, 0, 0.009, 0, 0,   0, 0, 0, 0.045, 0,  0, 0, 0, 0, 1,  0, 0, 0, 0, 0), nrow=J, byrow=TRUE);
A1 <- matrix(c(0, 0.05, 0.2, 0, 0,  0, 0, 0.009, 0, 0,   0, 0, 0, 0.045, 0,  0, 0, 0, 0, 1,  0, 0, 0, 0, 0), nrow=J, byrow=TRUE);
A2 <- matrix(c(0, 0.05, 0, 0.2, 0,  0, 0, 0.009, 0, 0,   0, 0, 0, 0.045, 0,  0, 0, 0, 0, 1,  0, 0, 0, 0, 0), nrow=J, byrow=TRUE);

#A <- A0;
A <- A1;
#A <- A2;

# ensure that the diagonals are 0
stopifnot(diag(A) == rep(0, nrow(A)))

# net proliferation rate
# assume that beta is invariant across time
#beta <- runif(J, -0.5, 0.5);  # arbitrary upperbound
beta <- c(0.3, 0.2, 0.1, 0.05, 0);
#beta <- c(0, 0.1, 0.3, 0.05, 0);
#beta <- c(0, 0.01, 0.04, 1, 0);

# relative limiting capacities
nl <- 1e9;
kappa <- c(10, 1, 2.9, 9, 52);
#kappa <- c(10, 1, 3, 9, 50);

params <- list(
	J = J, N = N, T = T, L = 2, s = s,
	n0 = n0, nl = nl,
	lambda = lambda, A = A, beta = beta,
	kappa = kappa
);


# --- Functions

cell_factor <- function(j) {
	factor(j, levels=1:J, labels=c("HEC", "HSC", "ST-HSC", "MPP", "other"))
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

	f <- n[, , , params$L] / pmax(1, apply(n, c(1, 2, 3), sum));

	list(n = n, f = f)
}

# --- 

out <- simulate_trajectory(params, update_ntl_logistic);

# --- Generate plots

graphics.off()

f <- out$f;

f.d <- melt(f, varnames=names(dim(f)));
f.d$j <- cell_factor(f.d$j);

g <- ggplot(f.d, aes(x = t, y = value, colour = factor(j))) +
	theme_classic() +
	geom_line(linetype=2) + facet_grid(n ~ j) +
	guides(colour = "none") +
	xlab("analysis time") + ylab("% labelled")
qdraw(g, width = 6, file = insert(pdf.fname, c("latent", "label-prop")))


# normalize against the common progenitor
fn <- f / f[1:N, rep(1, J), 1:T];
fn <- apply(fn, c(2, 3), mean);

fn.d <- melt(fn, varnames=names(dim(f)[-1]));
fn.d$j <- cell_factor(fn.d$j);

g <- ggplot(fn.d, aes(x = t, y = value, colour = factor(j))) +
	theme_classic() +
	geom_hline(yintercept = 1, linetype=3, colour="grey30") +
	geom_line(linetype=2) + 
	facet_grid(j ~ ., scales="free_y") +
	guides(colour = "none") +
	xlab("analysis time") + ylab("label ratio vs. HEC")
qdraw(g, file = insert(pdf.fname, c("latent", "label-ratio")))


# normalize against the second progenitor
fn2 <- f / f[1:N, rep(2, J), 1:T];
fn2 <- apply(fn2, c(2, 3), mean);

fn2.d <- melt(fn2, varnames=names(dim(f)[-1]));
fn2.d$j <- cell_factor(fn2.d$j);

g <- ggplot(fn2.d[fn2.d$j != "HEC", ], aes(x = t, y = value, colour = factor(j))) +
	theme_classic() +
	geom_hline(yintercept = 1, linetype=3, colour="grey30") +
	geom_line() + facet_grid(j ~ ., scales="free_y") +
	guides(colour = "none") +
	xlab("analysis time") + ylab("label ratio vs. HSC")
qdraw(g, file = insert(pdf.fname, c("observed", "label-ratio")))

qwrite(params, insert(out.fname, tag="params", ext="rds"));

