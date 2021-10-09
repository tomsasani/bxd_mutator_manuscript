# NOTE: all rules in `generate_figures.smk` use relative paths,
# since I assume you're running this pipeline locally. 

muts = ["C>A", "C>T", "C>G", "A>T", "A>G", "A>C"]

chroms = list(map(str, range(1, 20)))
chroms = ['chr' + c for c in chroms]

samples, = glob_wildcards("data/coverage_files/{sample}.bam.thresholds.bed")

# include separate "rule" files for each figure
include: "rules/make_figure_one.smk"
include: "rules/make_figure_two.smk"
include: "rules/make_figure_three.smk"
include: "rules/make_figure_four.smk"
include: "rules/make_figure_four_b.smk"
include: "rules/make_figure_five.smk"
include: "rules/make_supp_figures.smk"
include: "rules/make_main_bxd_paper_figs.smk"


# pseudo-rule to collect all output figures
main_figures = ["1a", "1b", 
				   "2a", "2b", "2c", 
				   "2d",
				   "3a", "3b",
				   "4a", "4c", "4d",
				   "5a", "5b", "5c"]

supp_figures = ["2", 
"3a", 
"4",
"5a", "5b",
				"6a", "6b", "6c", "6d", "7",]


rule all:
	input:
		# generate all of the main and supplementary figures
		expand("plots/figure_{fig_num}.eps", fig_num=main_figures),
		expand("plots/supp_figure_{fig_num}.eps", fig_num=supp_figures),
		"plots/ashbrook/epoch_sharing_heatmap.eps",
		"plots/ashbrook/figure_1a.eps",
		"plots/ashbrook/figure_1b.eps",
		"plots/figure_4b.pdf",
		# we generate supplementary figure 3 a little differently, since it
		# comprises a sub-panel for every mutation type
		expand("plots/all_qtl_maps/supp_figure_3_{mut_type}.eps", 
						mut_type = [m.replace('>', '.') for m in muts])

rule annotate_vars:
	input: 
		var_list = expand("data/{{var_type}}_vars/{chrom}_{{var_type}}_vars.exclude.csv", chrom=chroms),
		strain_metadata = "data/bam_names_to_metadata.xlsx",
		strain_callable_bp = "data/callable_bp.per_sample.csv",
		annotate_py = "py_scripts/annotate_vars.py",
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
		tidy_data_py = "py_scripts/make_tidy_data.py"
	output:
		mut_rates = "csv/tidy_mutation_rates.csv",
		mut_spectra = "csv/tidy_mutation_spectra.csv"
	shell:
		"""
		python {input.tidy_data_py} \
			   --annotated_singleton {input.annotated_singletons} 
		"""

rule calculate_callable_bp:
	input:
		py_script = "py_scripts/calculate_callable_bp.py",
		coverage_fhs = expand("data/coverage_files/{sample}.bam.thresholds.bed", sample=samples),
		exclude = "data/mm10.seg_dups.simple_repeats.merged.bed.gz"
	output:
		"data/callable_bp.per_sample.csv"
	shell:
		"""
		python {input.py_script} --coverage_files {input.coverage_fhs} \
								 --exclude {input.exclude} \
								 --out {output}
		"""
