#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

data <- read.table(args[1])

png(file = args[2], width=15, height=10, units="cm", res=1200)

barplot(data[,2], names.arg = data[,1], xlab = "Chromosome", ylab = "Count", col = "blue", main = paste0("Chrosome distribution_", args[3]))

dev.off()