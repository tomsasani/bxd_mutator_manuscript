
rule make_figure_one_abd: 
	input:
		py_script = "scripts/make_boxplots.py",
		mut_rates = "csv/tidy_mutation_rates.csv",
		mut_spectra = "csv/tidy_mutation_spectra.csv"
	output:
		"plots/figure_1a.eps",
		"plots/figure_1b.eps",
		"plots/figure_1d.eps"
	shell:
		"""
		python {input.py_script} --tidy_rates {input.mut_rates} \
								  --tidy_spectra {input.mut_spectra} \
								  --outdir /Users/tomsasani/harrislab/bxd_mutator_ms/plots/
		"""

rule make_figure_one_c: 
	input:
		py_script = "scripts/plot_conservation.py",
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		annotated_fixed = "csv/annotated_fixed_vars.csv",
	output:
		"plots/figure_1c.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
							   --annotated_fixed {input.annotated_fixed} \
							   --out {output} 
		"""