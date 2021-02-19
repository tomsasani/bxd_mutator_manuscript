muts = ["C>A", "C>T", "C>G", "A>T", "A>G", "A>C"]

chroms = list(map(str, range(1, 20)))
chroms = ['chr' + c for c in chroms]

# include separate "rule" files for each figure
include: "rules/make_figure_one.smk"
include: "rules/make_figure_two.smk"
include: "rules/make_figure_three.smk"
include: "rules/make_figure_four.smk"
include: "rules/make_supp_figures.smk"

# pseudo-rule to collect all output figures
main_figures = ["1a", "1b", "1c", "1d", 
				   #"2a", "2b", "2c", 
				   "2d",
				   "3a", "3b",
				   "4a", "4bc", "4d"]

supp_figures = ["2", 
#"3a", 
"4a", "4b", "5a", "5b",
				"6a", "6b", "6c", "6d", "7"]

rule all:
	input:
		# generate all of the main and supplementary figures
		expand("plots/figure_{fig_num}.eps", fig_num=main_figures),
		expand("plots/supp_figure_{fig_num}.eps", fig_num=supp_figures),
		# we generate supplementary figure 3 a little differently, since it
		# comprises a sub-panel for every mutation type
		#expand("plots/all_qtl_maps/supp_figure_3_{mut_type}.eps", 
		#				mut_type = [m.replace('>', '.') for m in muts])

rule annotate_vars:
	input: 
		var_list = expand("data/{{var_type}}_vars/{chrom}_{{var_type}}_vars.exclude.csv", chrom=chroms),
		strain_metadata = "data/bam_names_to_metadata.xlsx",
		strain_callable_bp = "data/smp2denominator.csv",
		annotate_py = "scripts/annotate_vars.py",
		qtl_geno = "Rqtl_data/bxd.geno.new"
	output:
		"csv/annotated_{var_type}_vars.csv"
	shell:
		"""
		python {input.annotate_py} \
			   --var_list {input.var_list} \
			   --strain_genotypes {input.qtl_geno} \
			   --strain_metadata {input.strain_metadata} \
			   --callable_bp {input.strain_callable_bp} \
			   --out {output}
		"""

rule generate_tidy_data:
	input:
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		tidy_data_py = "scripts/make_tidy_data.py"
	output:
		mut_rates = "csv/tidy_mutation_rates.csv",
		mut_spectra = "csv/tidy_mutation_spectra.csv"
	shell:
		"""
		python {input.tidy_data_py} \
			   --annotated_singleton {input.annotated_singletons} \
			   --out_pref /Users/tomsasani/harrislab/bxd_mutator_ms/csv/
		"""


rule make_epoch_sharing_fig:
	input:
		annotated_vars = "csv/annotated_bxd_private.csv",
		heatmap_py = "scripts/make_epoch_sharing_heatmap.py"
	output:
		"plots/epoch_heatmap.eps"
	shell:
		"""
		python {input.heatmap_py} --annotated_vars {input.annotated_vars} \
									 --out {output}
		"""
