library(cowplot)
library(optparse)
library(ggplot2)

option_list = list(
  make_option(c("-s", "--tidy_spectra"), type="character", default=NULL))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# read in tidy dataframes
tidy_spectra = read.csv(opt$tidy_spectra)

tidy_spectra = subset(tidy_spectra, estimate_type == "fraction" & base_mut == "C>A")

# make figure
setEPS()
fname = "plots/figure_2b.eps"
postscript(fname, width=4, height=3.5)
ggplot(tidy_spectra, aes(x=haplotype_at_qtl, y=estimate)) +
    # plot boxplot
    geom_boxplot(outlier.shape=NA, fill='white', col='lightgrey') +
    # plot jitter
    geom_jitter(aes(fill=haplotype_at_qtl), width=0.25, pch=21, col='white', size=3) +
    theme_cowplot() +
    theme(legend.position = "none")  +
    labs(x="Haplotype at QTL",
         y="C>A singleton fraction") +
    scale_fill_manual(values=c("cornflowerblue", "firebrick"))

dev.off()