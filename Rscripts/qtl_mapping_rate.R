library(qtl2)
library(ggplot2)
library(cowplot)
library(dplyr)
library(optparse)
library(RNOmni)
library(stringr)

option_list = list(
  make_option(c("-j", "--json"), type="character", default=NULL),
  make_option(c("-p", "--phenotype_file"), type="character", default=NULL))

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

print (head(phen_df))

phen_df = subset(phen_df, bxd_strain_conv != "BXD68_RwwJ_0462")

# calculate kinship between strains using the
# "leave one chromosome out" method
k = calc_kinship(pr, 'loco')

# get special covariates for the X
Xcovar <- get_x_covar(bxd)


for (mut_type in c("C>A", "C>T", "C>G", "A>T", "A>G", "A>C", "CpG>TpG"))
  {
    phen_df_sub = subset(phen_df, base_mut == mut_type)
    
    # get the phenotype as a log10-transformed fraction...
    phen_df_sub_frac = subset(phen_df_sub, estimate_type == "fraction")
    phen_matrix_frac = as.matrix(log10(phen_df_sub_frac$estimate))
    #phen_matrix_frac = as.matrix(RankNorm(phen_df_sub_frac$estimate))
    
    phenotype_frac = as.matrix(phen_matrix_frac[,1])
    rownames(phenotype_frac) = phen_df_sub_frac$bxd_strain_conv
    
    # and as a rate
    phen_df_sub_rate = subset(phen_df_sub, estimate_type == "rate")
    phen_matrix_rate = as.matrix(phen_df_sub_rate$estimate)
    phenotype_rate = as.matrix(phen_matrix_rate[,1])
    rownames(phenotype_rate) = phen_df_sub_rate$bxd_strain_conv
    
    # get covariates to include
    covariate_cols = c("epoch", "n_intercross_gens")
    covariate_matrix = as.matrix(phen_df_sub_frac[covariate_cols])
    rownames(covariate_matrix) = phen_df_sub_frac$bxd_strain_conv
    
    # perform a genome scan, accounting for kinship and
    # epoch as an additive covarirate
    out_rate <- scan1(pr, phenotype_rate, kinship=k, 
                      addcovar=covariate_matrix, Xcovar = Xcovar)
    
    out_frac <- scan1(pr, phenotype_frac, kinship=k, 
                      addcovar=covariate_matrix, Xcovar = Xcovar)
    
    # perform a permutation test to assess significance
    operm_rate <- scan1perm(pr, phenotype_rate, kinship=k, 
                            addcovar=covariate_matrix, Xcovar=Xcovar, n_perm=100)
    
    operm_frac <- scan1perm(pr, phenotype_frac, kinship=k, 
                            addcovar=covariate_matrix, Xcovar=Xcovar, n_perm=100)
    
    # get the LOD threshold for a < 0.05
    lod_cutoff_sig_rate = summary(operm_rate, alpha=0.05 / 15)[1]
    lod_cutoff_sig_frac = summary(operm_frac, alpha=0.05 / 15)[1]
    
    # plot peaks and LOD threshold
    ymx_rate <- maxlod(out_rate)
    ymx_frac <- maxlod(out_frac)
    
    ylim_val = max(c(ymx_frac, ymx_rate, lod_cutoff_sig_frac, lod_cutoff_sig_rate))
    
    formatted_mut_type = str_replace(mut_type, ">", ".")
    
    # plot LOD scores genome-wide for fraction phenotype
    setEPS()
    suff = ".eps"
    outfile = sprintf("/Users/tomsasani/harrislab/bxd_mutator_ms/plots/all_qtl_maps/supp_figure_3_%s%s", formatted_mut_type, suff)
    print (outfile)
    
    postscript(outfile, width=10, height=5)
    par(mar=c(4.1, 4.1, 1.6, 1.1))
    color <- c("green3", "slateblue")
    plot(out_rate, pmap, lodcolumn=1, col=color[1], ylim=c(0, ylim_val*1.15))
    plot(out_frac, pmap, lodcolumn=1, col=color[2], ylim=c(0, ylim_val*1.15), add=T)
    abline(h=lod_cutoff_sig_frac, col='slateblue', lwd=2, lty=2)
    abline(h=lod_cutoff_sig_rate, col='green3', lwd=2, lty=2)
    legend_name_rate = sprintf("%s rate", mut_type)
    legend_name_frac = sprintf("%s fraction", mut_type)
    
    legend("topleft", lwd=2, col=color, c(legend_name_rate, legend_name_frac), bg="gray90", lty=c(1,1,2))
    dev.off()
    
  }



