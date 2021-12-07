rule make_epoch_sharing_heatmap:
	input:
		py_script = "py_scripts/make_epoch_sharing_heatmap.py",
		annotated_shared_vars = "csv/annotated_shared_vars.csv"
	output:
		"plots/ashbrook/epoch_sharing_heatmap.eps"
	shell:
		"""
		python {input.py_script} --annotated_vars {input.annotated_shared_vars} \
							   --out {output} 
		"""

rule make_boxplots:
	input:
		R_script = "Rscripts/make_figure_one_ashbrook.R",
		tidy_rate = "csv/tidy_mutation_rates.csv",
		tidy_spectra = "csv/tidy_mutation_spectra.csv"
	output:
		"plots/ashbrook/figure_1a.eps",
		"plots/ashbrook/figure_1b.eps"
	shell:
		"""
		Rscript {input.R_script} -r {input.tidy_rate} \
								 -s {input.tidy_spectra}
		"""