#######################################
#                                     #
#       Simulate with Polyester       #
#                                     #
#######################################

# Author: Jorge Mestre
# Last modified: 18/03/2022 by Jorge Mestre

suppressMessages(library(polyester))
suppressMessages(library(Biostrings))

args <- commandArgs(trailingOnly = TRUE)
ref.trans <- args[1]
count.mat <- args[2]
out.dir <- args[3]
in.seed <- args[4]

count.mat <- read.table(count.mat, header=F)
count.mat <- as.matrix(count.mat)

simulate_experiment_countmat(ref.trans, readmat=count.mat, outdir=out.dir,
                             paired = TRUE, readlen=100, seed=in.seed) 