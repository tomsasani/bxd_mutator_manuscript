
rule make_figure_three_a:
	input:
		annotated_singletons = "csv/annotated_singletons.csv",
		py_script = "scripts/mutation_comparison.py"
	output:
		"plots/figure_3a.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
							   --out {output} \
							   -subset_key haplotype_at_qtl \
							   -plot_type heatmap
		"""

rule make_figure_three_b:
	input:
		annotated_singletons = "csv/annotated_singletons.csv",
		ohno_data = "data/41598_2014_BFsrep04689_MOESM2_ESM.xls",
		py_script = "scripts/compare_signatures_TOYKO.py"
	output:
		"plots/figure_3b.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
						 	   --ohno_data {input.ohno_data} \
							   --out {output}
		"""