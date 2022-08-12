library(io)

# times are in hours

delta <- 0.20;

# number of samples
N <- 1e8;

# number of cell types
J <- 5;

# HSC	  ST HSC	MPP	  CLP	  CD19 B-pro
# 4291	10562	  99220	44683	1083

# HSC, ST-HSC, MMP, CMP, CLP
target <- c(4291, 10562, 99220, NA, 44683);
#target <- c(4291, NA, NA, NA, NA);

# differentiation and net proliferation rates from Busch et al. 2015
A.lower <- matrix(
	c(
		# HSC
		0,  0,  0,     0,     0,
		# ST-HSC
		0,     0,  0,  0,     0,
		# MMP
		0,     0,     0,  0,  0,
		# CMP
		0,     0,     0,     0,     0,
		# CLP
		0,     0,     0,     0,     0
	),
	ncol=J, byrow=TRUE
);

A.upper <- matrix(
	c(
		# HSC
		4/24,  4/24,  0,     0,     0,
		# ST-HSC
		0,     4/24,  4/24,  0,     0,
		# MMP
		0,     0,     4/24,  4/24,  4/24,
		# CMP
		0,     0,     0,     0,     0,
		# CLP
		0,     0,     0,     0,     0
	),
	ncol=J, byrow=TRUE
);

n0.lower <- 1;
n0.upper <- 100;

# number of time points
T.lower <- 5 * 24;
T.upper <- 6 * 24;

sample_params <- function() {
	A0 <- matrix(runif(prod(dim(A.lower)), min=A.lower, max=A.upper), nrow=nrow(A.lower));

	# net proliferation rate
	beta <- diag(A0);

	# differentiation rates
	A <- A0;
	diag(A) <- 0;

	# total efflux rates
	alpha.out <- rowSums(A);

	n0 <- floor(runif(1, min=n0.lower, max=n0.upper));
	
	T <- floor(runif(1, min=T.lower, max=T.upper));

	list(
		n0 = n0, beta = beta, alpha.out = alpha.out, A = A, T = T
	)
}

update_n_exp <- function(n, params) {
	n.net.prolif <- rpois(length(n), n * params$beta);
	# n.efflux <= n
	n.efflux <- pmin(n, rpois(length(n), n * params$alpha.out));

	# normalize probability within each row
	probs <- params$A / params$alpha.out;
	probs[is.nan(probs)] <- 0;

	n.influx <- integer(length(n));
	for (i in 1:length(n.efflux)) {
		if (sum(probs[i, ]) > 0) {
			n.influx <- n.influx + rmultinom(1, n.efflux[i], probs[i, ]);
		}
	}

	n + n.net.prolif - n.efflux + n.influx
}

sim <- function(params) {
	# initialize
	n <- integer(J);
	n[1] <- params$n0;

	# iterate
	for (t in 2:params$T) {
		n <- update_n_exp(n, params);
	}

	n
}

accepted <- list();
for (i in 1:N) {

	if (i %% 1000 == 0) {
		message("round ", i)
	}

	params <- sample_params();
	n <- sim(params);
	score <- mean(abs(target - n) / n, na.rm=TRUE)

	if (score < delta) {
		message("acceptance")
		params$n <- n;
		accepted[[length(accepted) + 1]] <- params
	}
}

n0 <- unlist(lapply(accepted, function(x) x$n0));
beta1 <- unlist(lapply(accepted, function(x) x$beta[1])) * 24;
alpha.out1 <- unlist(lapply(accepted, function(x) x$alpha.out[1])) * 24;

length(accepted)

hist(beta1)
hist(n0)
range(n0)

qwrite(accepted, "accepted.rds")

