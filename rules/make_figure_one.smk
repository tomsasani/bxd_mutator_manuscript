# set name of working directory
WORKDIR = "/Users/tomsasani/harrislab/bxd_mutator_ms/"

rule make_figure_one_abd: 
	input:
		py_script = WORKDIR + "scripts/make_boxplots.py",
		mut_rates = WORKDIR + "csv/tidy_mutation_rates.csv",
		mut_spectra = WORKDIR + "csv/tidy_mutation_spectra.csv"
	output:
		WORKDIR + "plots/figure_1a.eps",
		WORKDIR + "plots/figure_1b.eps",
		WORKDIR + "plots/figure_1d.eps"
	shell:
		"""
		python {input.py_script} --tidy_rates {input.mut_rates} \
								  --tidy_spectra {input.mut_spectra} \
								  --outdir /Users/tomsasani/harrislab/bxd_mutator_ms/plots/
		"""

rule make_figure_one_c: 
	input:
		py_script = WORKDIR + "scripts/plot_conservation.py",
		annotated_singletons = WORKDIR + "csv/annotated_singletons.csv",
		annotated_common = WORKDIR + "csv/annotated_common.csv",
	output:
		WORKDIR + "plots/figure_1c.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
							   --annotated_common {input.annotated_common} \
							   --out {output} 
		"""