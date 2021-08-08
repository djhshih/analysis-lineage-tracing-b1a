library(io)
library(ggplot2)

hec <- qread("data/hec.rds");
hsc <- qread("data/hsc.rds");

ggplot(hec, aes(t_injection, bm_fl.hsc)) +
	geom_point()

