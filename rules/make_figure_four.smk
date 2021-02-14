rule make_figure_four_a:
	input:
		mgp_mutyh_vars ="data/mgp.mutyh.csv",
		py_script = "scripts/make_mutyh_grid.py"
	output:
		"plots/figure_4a.eps"
	shell:
		"""
		python {input.py_script} --strain_vars {input.mgp_mutyh_vars} \
									--out {output}
		"""

rule make_figure_four_bc:
	input:
		dumont_xls = "data/SuppTables_concat.xlsx",
		py_script = "scripts/compare_mgp_spectra.py"
	output:
		"plots/figure_4bc.eps"
	shell:
		"""
		python {input.py_script} --dumont_xls {input.dumont_xls} \
									--out {output}
		"""
rule make_figure_four_d:
	input:
		wild_mutyh_vars = "data/wild.mutyh.csv",
		py_script = "scripts/make_mutyh_grid.py"
	output:
		"plots/figure_4d.eps",
	shell:
		"""
		python {input.py_script} --strain_vars {input.wild_mutyh_vars} \
									--out {output} \
									-is_wild
		"""