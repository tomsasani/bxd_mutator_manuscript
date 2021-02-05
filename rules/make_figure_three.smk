# set name of working directory
WORKDIR = "/Users/tomsasani/harrislab/bxd_mutator_ms/"

rule make_figure_three_a:
	input:
		annotated_singletons = WORKDIR + "csv/annotated_singletons.csv",
		py_script = WORKDIR + "scripts/mutation_comparison.py"
	output:
		WORKDIR + "plots/figure_3a.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
							   --out {output} \
							   -subset_key haplotype_at_qtl \
							   -plot_type heatmap
		"""

rule make_figure_three_b:
	input:
		annotated_singletons = WORKDIR + "csv/annotated_singletons.csv",
		ohno_data = WORKDIR + "data/ohno_etal.xls",
		py_script = WORKDIR + "scripts/compare_signatures_TOYKO.py"
	output:
		WORKDIR + "plots/figure_3b.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
						 	   --ohno_data {input.ohno_data} \
							   --out {output}
		"""