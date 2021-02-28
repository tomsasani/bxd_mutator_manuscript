rule make_figure_one_ab: 
	input:
		R_script = "Rscripts/make_figure_one.R",
		mut_rates = "csv/tidy_mutation_rates.csv",
		mut_spectra = "csv/tidy_mutation_spectra.csv"
	output:
		"plots/figure_1a.eps",
		"plots/figure_1b.eps",
	shell:
		"""
		Rscript {input.R_script} -r {input.mut_rates} \
								  -s {input.mut_spectra} 
		"""