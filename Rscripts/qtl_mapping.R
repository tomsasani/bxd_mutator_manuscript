library(qtl2)
library(ggplot2)
library(cowplot)
library(dplyr)
library(optparse)
library(tidyr)
library(RNOmni)

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

# and further subset to only include the relevant mutation type
phen_df_sub = subset(phen_df, base_mut == "C>A")

# calculate kinship between strains using the
# "leave one chromosome out" method
k = calc_kinship(pr, 'loco')

# get special covariates for the X
Xcovar <- get_x_covar(bxd)

# get the phenotype as a log10-transformed fraction...
phen_df_sub_frac = subset(phen_df_sub, estimate_type == "fraction")
phen_matrix_frac = as.matrix(log2(subset(phen_df_sub_frac, bxd_strain_conv != "BXD68_RwwJ_0462")$estimate))
#phen_matrix_frac = as.matrix(subset(phen_df_sub_frac, bxd_strain_conv != "BXD68_RwwJ_0462")$estimate)
phen_matrix_frac = as.matrix(RankNorm(subset(phen_df_sub_frac, bxd_strain_conv != "BXD68_RwwJ_0462")$estimate))

phenotype_frac = as.matrix(phen_matrix_frac[,1])
strain_names = subset(phen_df_sub_frac, bxd_strain_conv != "BXD68_RwwJ_0462")$bxd_strain_conv
rownames(phenotype_frac) = strain_names

# and as a rate
phen_df_sub_rate = subset(phen_df_sub, estimate_type == "rate")
phen_matrix_rate = as.matrix(subset(phen_df_sub_rate, bxd_strain_conv != "BXD68_RwwJ_0462")$estimate)
phen_matrix_rate = as.matrix(RankNorm(subset(phen_df_sub_rate, bxd_strain_conv != "BXD68_RwwJ_0462")$estimate))

strain_names = subset(phen_df_sub_rate, bxd_strain_conv != "BXD68_RwwJ_0462")$bxd_strain_conv
phenotype_rate = as.matrix(phen_matrix_rate[,1])
rownames(phenotype_rate) = strain_names

# get covariates to include
phen_df_sub_frac$is_ail[phen_df_sub_frac$n_intercross_gens == 0] <- 0
phen_df_sub_frac$is_ail[phen_df_sub_frac$n_intercross_gens > 0] <- 1


#count_covar = subset(phen_df_sub, estimate_type == "count" & base_mut == "C>A")$estimate

# get covariates
covariate_cols = c("n_intercross_gens", "is_ail", "epoch")
covariate_matrix = as.matrix(phen_df_sub_frac[covariate_cols])
#covariate_matrix_with_count = cbind(covariate_matrix, count_covar)
rownames(covariate_matrix) = phen_df_sub_frac$bxd_strain_conv
#rownames(covariate_matrix_with_count) = phen_df_sub_frac$bxd_strain_conv

Xcovar = get_x_covar(bxd)

# perform a genome scan, accounting for kinship and
# epoch as an additive covarirate
out_rate <- scan1(pr, phenotype_rate, kinship=k, 
             addcovar=covariate_matrix, Xcovar=Xcovar)

out_frac <- scan1(pr, phenotype_frac, kinship=k, 
                  addcovar=covariate_matrix, Xcovar=Xcovar)

print (dim(phenotype_rate))
#print (phen_df_sub_rate$bxd_strain_conv)

# perform a permutation test to assess significance
operm_rate <- scan1perm(pr, phenotype_rate, kinship=k, 
                   addcovar=covariate_matrix, n_perm=100)

operm_frac <- scan1perm(pr, phenotype_frac, kinship=k, 
                        addcovar=covariate_matrix, n_perm=100)

# get the LOD threshold for a < 0.05
lod_cutoff_sig_rate = summary(operm_rate, alpha=0.05 / 15)[1]
lod_cutoff_sig_frac = summary(operm_frac, alpha=0.05 / 15)[1]

print (lod_cutoff_sig_rate)
print (lod_cutoff_sig_frac)

print (bayes_int(out_frac, pmap, lodcolumn=1, chr=4, prob=0.95))
print (bayes_int(out_rate, pmap, lodcolumn=1, chr=4, prob=0.95))

print(find_peaks(out_frac, pmap, threshold=lod_cutoff_sig_frac))
print(find_peaks(out_rate, pmap, threshold=lod_cutoff_sig_rate))

# plot peaks and LOD threshold
ymx_rate <- maxlod(out_rate)
ymx_frac <- maxlod(out_frac)

# plot LOD scores genome-wide for fraction phenotype
setEPS()
fname = "figure_2a.eps"
outfile = sprintf("%s/%s", opt$out_prefix, fname)
postscript(outfile, width=7, height=4)
par(mar=c(4.1, 4.1, 1.6, 1.1))
color <- c("green3", "slateblue")
plot(out_rate, pmap, lodcolumn=1, col=color[1], ylim=c(0, ymx_frac*1.05))
plot(out_frac, pmap, lodcolumn=1, col=color[2], ylim=c(0, ymx_frac*1.05), add=T)
abline(h=lod_cutoff_sig_frac, col='slateblue', lwd=2, lty=2)
abline(h=lod_cutoff_sig_rate, col='green3', lwd=2, lty=2)

legend("topright", lwd=2, col=color, c("C>A rate", "C>A fraction"), bg="gray90", lty=c(1,1,2))

dev.off()

setEPS()
fname = "figure_2b.eps"
outfile = sprintf("%s/%s", opt$out_prefix, fname)
postscript(outfile, width=4, height=4)
par(mar=c(4.1, 4.1, 1.6, 1.1))
color <- c("green3", "slateblue")
plot(out_rate, pmap$`4`, lodcolumn=1, col=color[1], ylim=c(0, ymx_frac*1.05))
plot(out_frac, pmap$`4`, lodcolumn=1, col=color[2], ylim=c(0, ymx_frac*1.05), add=T)
abline(h=lod_cutoff_sig_frac, col='slateblue', lwd=2, lty=2)
abline(h=lod_cutoff_sig_rate, col='green3', lwd=2, lty=2)
legend("topleft", lwd=2, col=color, c("C>A rate", "C>A fraction"), bg="gray90", lty=c(1,1,2))
dev.off()

print ('done with plotting QTL')

# find the maximum LOD peak
max_peak_frac = find_peaks(out_frac, pmap, threshold = lod_cutoff_sig_frac)[1,]

# below is some file formatting to be able to plot
# a phenotype x genotype plot 
g_frac <- maxmarg(pr, pmap, chr=max_peak_frac$chr, 
                  pos=max_peak_frac$pos, return_char=TRUE)

g_frac = setNames(stack(g_frac)[2:1], c('strain','haplotype'))

g_frac = g_frac %>% replace_na(list(haplotype = "H"))

vars_to_include = c("bxd_strain_conv", "estimate")

p_frac = phen_df_sub_frac[vars_to_include]
colnames(p_frac) <- c("strain", "fraction")

g_new_frac = inner_join(p_frac, g_frac)$haplotype
names(g_new_frac) = inner_join(p_frac, g_frac)$strain

p_new_frac = inner_join(p_frac, g_frac)$fraction
names(p_new_frac) = inner_join(p_frac, g_frac)$strain

setEPS()
fname = "figure_2c.eps"
outfile = sprintf("%s/%s", opt$out_prefix, fname)
postscript(outfile, width=4, height=4)
plot_pxg(g_new_frac, p_new_frac, SEmult=2, 
         ylab="C>A singleton mutation fraction",
         xlab="Founder haplotype background")
dev.off()


