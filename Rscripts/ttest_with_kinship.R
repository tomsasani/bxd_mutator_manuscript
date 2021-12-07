library(qtl2)
#library(ggplot2)
#library(cowplot)
library(dplyr)
library(optparse)
#library(gaston)
library(coxme)

option_list = list(
  make_option(c("-j", "--json"), type="character", default=NULL),
  make_option(c("-p", "--phenotype_file"), type="character", default=NULL))
      
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# read in the JSON that directs R/qtl2 to sample genotypes,
# phenotypes, and covariates
bxd <- read_cross2(opt$json)

# read in the phenotype values for each BXD strain
phen_df = read.csv(opt$phenotype_file, header=T)

# and further subset to only include the relevant mutation type
phen_df_sub = subset(phen_df, base_mut == "C>A" & estimate_type == "fraction")

# subset cross2 to relevant BXDs
bxd = bxd[phen_df_sub$bxd_strain_conv, ]

# insert pseudomarkers into the genotype map
# following lines are from https://kbroman.org/pages/teaching.html
gmap <- insert_pseudomarkers(bxd$gmap, step=0.2, stepwidth='max')
pmap <- interp_map(gmap, bxd$gmap, bxd$pmap)

# calculate QTL genotype probabilities
# (error_prob taken directly from R/qtl2 user guide)
pr <- calc_genoprob(bxd, gmap, error_prob=0.002, map_function="c-f")

# calculate kinship between strains using the
# "overall" method
k = calc_kinship(pr, 'overall')

# convert haplotype column to binary
phen_df_sub <- phen_df_sub %>% mutate(haplotype_at_qtl_binary = case_when(
    haplotype_at_qtl == "D" ~ 1 ,
    haplotype_at_qtl == "B"  ~ 0,
  )
)

x1 <- cbind(1, as.matrix(phen_df_sub$haplotype_at_qtl_binary))
Y <- phen_df_sub$estimate

dat <- data.frame(Y, x=x1[,2], id = 1:length(Y))

k = as.matrix(k)
colnames(k) <- c(1:length(Y))
rownames(k) <- c(1:length(Y))
# https://sahirbhatnagar.com/blog/2017/10/12/mixed-models-with-kinship-in-r/
print (k[1:5,1:5])
gfit1 <- lmekin(Y ~ x + (1|id), data=dat, varlist=k)

print (gfit1)