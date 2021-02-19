library(qtl2)
library(ggplot2)
library(cowplot)
library(dplyr)
library(optparse)
library(tidyr)

option_list = list(
  make_option(c("-j", "--json"), type="character", default=NULL),
  make_option(c("-p", "--phenotype_file"), type="character", default=NULL),
  make_option(c("-o", "--out_prefix"), type="character", default=NULL))
      
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# read in the JSON that directs R/qtl2 to sample genotypes,
# phenotypes, and covariates
bxd <- read_cross2(opt$json)

# insert pseudomarkers into the genotype map
# following lines are from https://kbroman.org/pages/teaching.html
gmap <- insert_pseudomarkers(bxd$gmap, step=0.2, stepwidth='max')
pmap <- interp_map(gmap, bxd$gmap, bxd$pmap)

# calculate QTL genotype probabilities
# (error_prob taken directly from R/qtl2 user guide)
pr <- calc_genoprob(bxd, gmap, error_prob=0.002, map_function="c-f")

# read in the phenotype values for each BXD strain
phen_df = read.csv(opt$phenotype_file, header=T)

# calculate kinship between strains using the
# "leave one chromosome out" method
k = calc_kinship(pr, 'loco')

# get special covariates for the X
Xcovar <- get_x_covar(bxd)

# get the phenotype as a log10-transformed fraction...
phen_matrix_rate = as.matrix(phen_df$rate)

phenotype_rate = as.matrix(phen_matrix_rate[,1])
strain_names = phen_df$bxd_strain_conv
rownames(phenotype_rate) = strain_names

# convert all covariates to numeric
phen_df$is_ail[phen_df$n_intercross_gens == 0] <- 0
phen_df$is_ail[phen_df$n_intercross_gens > 0] <- 1

# get covariates
covariate_cols = c("n_intercross_gens", "is_ail", "epoch")
covariate_matrix = as.matrix(phen_df[covariate_cols])
rownames(covariate_matrix) = phen_df$bxd_strain_conv

# perform a genome scan, accounting for kinship and
# epoch as an additive covarirate
out <- scan1(pr, phenotype_rate, kinship=k, 
             addcovar=covariate_matrix)

# perform a permutation test to assess significance
operm <- scan1perm(pr, phenotype_rate, kinship=k, 
                   addcovar=covariate_matrix, n_perm=100)

# get the LOD threshold for a < 0.05
lod_cutoff_sig = summary(operm, alpha=0.05 / 15)[1]

# plot peaks and LOD threshold
ymx <- max(c(maxlod(out), lod_cutoff_sig))

# plot LOD scores genome-wide for fraction phenotype
setEPS()
fname = "supp_figure_3a.eps"
outfile = sprintf("%s/%s", opt$out_prefix, fname)
postscript(outfile, width=7, height=4)
par(mar=c(4.1, 4.1, 1.6, 1.1))
color <- "green3"
plot(out, pmap, lodcolumn=1, col=color, ylim=c(0, ymx*1.05))
abline(h=lod_cutoff_sig, col='green3', lwd=2, lty=2)

legend("topright", lwd=2, col=color, "overall rate", bg="gray90", lty=c(1,1,2))

dev.off()