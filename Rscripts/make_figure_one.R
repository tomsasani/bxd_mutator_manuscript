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

# fit Poisson regression modeling numbers of singletons as function
# of inbreeding time
m = glm(total_muts ~ n_inbreeding_gens, family=poisson(link="identity"), data=tidy_rates)

print (summary(m))

# use model to predict Y values given X
preds = predict(m, type='response', se.fit=TRUE)

# manually get 95% confidence bands around slope
# use (1.96 * SEM)
tidy_rates$ci_lo = preds$fit - (1.96 * preds$se.fit)
tidy_rates$ci_hi = preds$fit + (1.96 * preds$se.fit)

# make figure 1a
f1a <- ggplot(tidy_rates) +
    # plot regression line
    geom_line(data=cbind(tidy_rates, pred=preds$fit), aes(x=n_inbreeding_gens, y=pred), col='black') +
    # plot confidence bands
    geom_ribbon(aes(x=n_inbreeding_gens, ymin=ci_lo, ymax=ci_hi), alpha=1, fill='grey') +
    geom_point(aes(x=n_inbreeding_gens, y=total_muts, fill=Epoch), pch=21, col='black', size=3) +
    theme_cowplot() +
    labs(x="Number of generations of inbreeding",
            y="Number of homozygous singletons") +
    scale_fill_manual(values=cbPalette)

ggsave("plots/figure_1a.eps", f1a, width=6, height=4)

# subset spectra to only include fractions
tidy_spectra = subset(tidy_spectra, estimate_type == "fraction")
m = aov(estimate ~ base_mut + Epoch, data=tidy_spectra)
print (anova(m, test="Chisq"))

# make figure 1b
f1b <- ggplot(tidy_spectra, aes(x=base_mut, y=estimate, fill=Epoch)) +
    geom_boxplot() +
    theme_cowplot() +
    labs(x="Mutation type", 
         y="Fraction of singletons") +
    scale_fill_manual(values=cbPalette)

ggsave("plots/figure_1b.eps", f1b, width=6, height=4)
