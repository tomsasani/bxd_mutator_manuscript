library(ggplot2)
library(cowplot)
library(optparse)

cbPalette <- cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

option_list = list(
  make_option(c("-r", "--tidy_rates"), type="character", default=NULL),
  make_option(c("-s", "--tidy_spectra"), type="character", default=NULL))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# read in tidy dataframes
tidy_rates = read.csv(opt$tidy_rates)
tidy_spectra = read.csv(opt$tidy_spectra)

# make sure epoch is treated as factor
tidy_rates$epoch = as.factor(tidy_rates$epoch)
tidy_spectra$epoch = as.factor(tidy_spectra$epoch)

names(tidy_rates)[names(tidy_rates) == 'epoch'] <- 'Epoch'
names(tidy_spectra)[names(tidy_spectra) == 'epoch'] <- 'Epoch'

# make figure 1a
f1a <- ggplot(tidy_rates, aes(x=Epoch, y=total_muts, col=Epoch)) +
    # plot boxplot
    geom_boxplot(outlier.shape=NA) +
    # plot jitter
    geom_jitter(width=0.25) +
    theme_cowplot() +
    theme(legend.position = "none")  +
    labs(x="Epoch",
         y="Number of homozgyous singletons") +
    scale_fill_manual(values=cbPalette)

ggsave("plots/ashbrook/figure_1a.eps", f1a, width=6, height=4)

# subset spectra to only include fractions
tidy_spectra = subset(tidy_spectra, estimate_type == "fraction")

# make figure 1b
f1b <- ggplot(tidy_spectra, aes(x=base_mut, y=estimate, col=base_mut)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(width=0.25) + 
    theme_cowplot() +
    theme(legend.position = "none")  +
    labs(x="Mutation type", 
         y="Fraction of singletons") +
    scale_fill_manual(values=cbPalette)

ggsave("plots/ashbrook/figure_1b.eps", f1b, width=6, height=4)
