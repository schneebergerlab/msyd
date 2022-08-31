#!/usr/bin/env Rscript

library(ggplot2)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)
lensdf <- read.table(args[1], header=TRUE)
lensdf$index <- 1:nrow(lensdf)
coredf <- lensdf[apply((lensdf > 0), 1, all), ]

print(lensdf[lensdf$ref > 170000, ])



ggplot(data=pivot_longer(lensdf, -c(index, chr), names_to='key', values_to='value'), aes(x = key, y=value, group=index, color=chr)) +
	#scale_y_log10() +
	ylab("length on genome") + 
	xlab("genome name") +
	geom_line(alpha=0.3)

ggplot() +
	xlim(0, 100000) +
	xlab("length on reference") +
	geom_histogram(data=lensdf, binwidth=2000, aes(x=ref, fill='crosssyn')) +
	geom_histogram(data=coredf, binwidth=2000, aes(x=ref, fill='coresyn'))

ggplot() +
	xlim(0, 10000) +
	xlab("length on reference") +
	geom_histogram(data=lensdf, binwidth=200, aes(x=ref, fill='crosssyn')) +
	geom_histogram(data=coredf, binwidth=200, aes(x=ref, fill='coresyn'))

ggplot() +
	xlim(100000, NA) +
	xlab("length on reference") +
	#scale_y_log10() +
	geom_histogram(data=lensdf, aes(x=ref, fill='crosssyn')) +
	geom_histogram(data=coredf, aes(x=ref, fill='coresyn'))

ggplot() +
	scale_y_log10() +
	xlab("length on reference") +
	xlim(0, 100000) +
	geom_freqpoly(data=lensdf, aes(x=ref, color=chr))

ggplot() +
	xlim(0, 5000) +
	xlab("length on reference") +
	geom_freqpoly(data=lensdf, aes(x=ref, color=chr))




