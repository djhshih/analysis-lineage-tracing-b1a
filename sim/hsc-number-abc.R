library(io)

# times are in days
tsize <- 24;

# number of samples
N <- 5e7;
#N <- 1e7;

# number of cell types
J <- 5;

cells <- c("HSC", "ST-HSC", "MMP", NA, "CLP");
# HSC, ST-HSC, MMP, CMP, CLP
target <- c(4291, 10562, 99220, NA, 44683);
#target <- c(4291, NA, NA, NA, NA);
target.sd <- c(791, 695, 25019, NA, 30116);
alpha <- 0.01;
delta <- qnorm(1 - alpha/2);

resume <- TRUE;
accepted.fn <- "accepted.rds";

thresholds <- delta * target.sd / target;

# differentiation and net proliferation rates from Busch et al. 2015
A.lower <- matrix(
	c(
		# HSC
		0,     0,     0,     0,     0,
		# ST-HSC
		0,     0,     0,     0,     0,
		# MMP
		0,     0,     0,     0,     0,
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
		4/tsize,  4/tsize,  0,     0,     0,
		# ST-HSC
		0,     4/tsize,  4/tsize,  0,     0,
		# MMP
		0,     0,     4/tsize,  4/tsize,  4/tsize,
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
T.lower <- 5 * tsize;
T.upper <- 6 * tsize;

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

if (resume) {
	set.seed();
	accepted <- qread(accepted.fn);
} else {
	accepted <- list();
}

for (i in 1:N) {

	if (i %% 1000 == 0) {
		message("round ", i)
	}

	params <- sample_params();
	pred <- sim(params);
	scores <- abs(pred - target) / target;

	if (all(scores < thresholds, na.rm=TRUE)) {
		message("acceptance")
		params$n <- pred;
		accepted[[length(accepted) + 1]] <- params;
	}
}

length(accepted)

n0 <- unlist(lapply(accepted, function(x) x$n0));
beta1 <- unlist(lapply(accepted, function(x) x$beta[1])) * tsize;
alpha.out1 <- unlist(lapply(accepted, function(x) x$alpha.out[1])) * tsize;

betas <- matrix(unlist(lapply(accepted, function(x) x$beta)), ncol=length(accepted)) * tsize;

rowMeans(betas)

plot(betas[1, ], betas[2, ])

ns <- matrix(unlist(lapply(accepted, function(x) x$n)), ncol=length(accepted));
pred <- rowMeans(ns)
pred.sd <- apply(ns, 1, sd);

library(ggplot2)
library(ggsci)

d <- rbind(
	data.frame(
		x = factor(cells, levels=cells[!is.na(cells)]),
		y = target,
		ymin = target - target.sd,
		ymax = target + target.sd,
		group = "observed"
	),
	data.frame(
		x = factor(cells, levels=cells[!is.na(cells)]),
		y = pred,
		ymin = pred - pred.sd,
		ymax = pred + pred.sd,
		group = "predicted"
	)
);

d <- d[complete.cases(d), ];

qdraw(
	ggplot(d, aes(x=x, y=y, ymin=ymin, ymax=ymax, fill=group)) + theme_classic() +
		geom_col(position=position_dodge()) + 
		geom_errorbar(width=0.3, position=position_dodge(width=0.9)) +
		scale_y_continuous() +
		scale_fill_npg() +
		theme(legend.title=element_blank()) +
		xlab("") + ylab("cell number")
	,
	file = "hsc-number_abc-fit.pdf"
)

mean(n0)
print(quantile(n0, c(0.025, 0.975)))

qdraw(
	ggplot(data.frame(x=n0), aes(x=x)) + theme_classic() +
		geom_histogram(binwidth=3) +
		#geom_vline(xintercept=quantile(n0, c(0.025, 0.975)), col="firebrick") +
		xlab("initial HSC number")
	,
	file="hsc-number_n0_hist.pdf"
)

hist(n0, breaks=50)
plot(density(n0, from=1, to=100, bw=4))
hist(beta1, breaks=50)
plot(n0, beta1)

smoothScatter(n0, beta1)

d <- MASS::kde2d(n0, beta1, n=60);
qdraw(
	{
		contour(d,
			col=hcl.colors(12, "viridis"),
			las=1,
			ylab="net prolif. rate of HSC",
			xlab="initial HSC number",
			ylim=c(0, 4)
		)
		#points(n0, beta1, pch=20, col="#00000020")
	},
	file = "hsc-number_contour_n0-beta1.pdf"
);

hist(n0[beta1 < 3], breaks=100)

quantile(n0, c(0.05))
mean(n0 <= 1)
summary(n0)

idx <- n0 == 1;
beta1[idx]

hist(alpha.out1)
plot(n0, alpha.out1)
dev.off()

qwrite(accepted, "accepted.rds")

