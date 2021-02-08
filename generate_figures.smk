# set name of working directory
WORKDIR = "/Users/tomsasani/harrislab/bxd_mutator_ms/"

# set the path to the CONDA YAML
CONDA_YAML="/Users/tomsasani/harrislab/bxd_mutator_ms/env.yaml"

muts = ["C>A", "C>T", "C>G", "A>T", "A>G", "A>C"]

# include separate "rule" files for each figure
include: WORKDIR + "rules/make_figure_one.smk"
include: WORKDIR + "rules/make_figure_two.smk"
include: WORKDIR + "rules/make_figure_three.smk"
include: WORKDIR + "rules/make_figure_four.smk"
#include: WORKDIR + "rules/make_figure_five.smk"
include: WORKDIR + "rules/make_supp_figures.smk"

# pseudo-rule to collect all output figures
main_figures = ["1a", "1b", "1c", "1d", 
				   "2a", "2b", "2c", "2d",
				   "3a", "3b",
				   "4a", "4bc"]

supp_figures = ["1", "4", "5a", "5b",
				"6a", "6b", "7"]

rule all:
	input:
		# generate all of the main and supplementary figures
		expand(WORKDIR + "plots/figure_{fig_num}.eps", fig_num=main_figures),
		expand(WORKDIR + "plots/supp_figure_{fig_num}.eps", fig_num=supp_figures),
		# we generate supplementary figure 3 a little differently, since it
		# comprises a sub-panel for every mutation type
		#expand(WORKDIR + "plots/all_qtl_maps/supp_figure_3_{mut_type}.eps", mut_type = [m.replace('>', '.') for m in muts])

# annotate a "BED-like" file of variants
rule annotate_vars:
	input: 
		var_dir = WORKDIR + "data/{var_type}/",
		strain_metadata = WORKDIR + "data/bam_names_to_metadata.xlsx",
		strain_callable_bp = WORKDIR + "data/smp2denominator.csv",
		annotate_py = WORKDIR + "scripts/annotate_vars.py",
		qtl_geno = WORKDIR + "Rqtl_data/bxd.geno.new.updated"
	output:
		WORKDIR + "csv/annotated_{var_type}.csv"
	conda:
		CONDA_YAML
	shell:
		"""
		python {input.annotate_py} \
			   --var_dir {input.var_dir} \
			   --strain_genotypes {input.qtl_geno} \
			   --strain_metadata {input.strain_metadata} \
			   --callable_bp {input.strain_callable_bp} \
			   --out {output}
		"""

# ---
# convert the "BED-like" file of annotated variants into
# a "tidy" dataframe 
# --- 
rule generate_tidy_data:
	input:
		annotated_singletons = WORKDIR + "csv/annotated_singletons.csv",
		tidy_data_py = WORKDIR + "scripts/make_tidy_data.py"
	output:
		mut_rates = WORKDIR + "csv/tidy_mutation_rates.csv",
		mut_spectra = WORKDIR + "csv/tidy_mutation_spectra.csv"
	shell:
		"""
		python {input.tidy_data_py} \
			   --annotated_singleton {input.annotated_singletons} \
			   --out_pref /Users/tomsasani/harrislab/bxd_mutator_ms/csv/
		"""


rule make_epoch_sharing_fig:
	input:
		annotated_vars = WORKDIR + "csv/annotated_bxd_private.csv",
		heatmap_py = WORKDIR + "scripts/make_epoch_sharing_heatmap.py"
	output:
		WORKDIR + "plots/epoch_heatmap.eps"
	shell:
		"""
		python {input.heatmap_py} --annotated_vars {input.annotated_vars} \
									 --out {output}
		"""
