rule make_figure_four_c:
	input:
		mgp_mutyh_vars ="data/mgp.mutyh_genotypes.csv",
		py_script = "py_scripts/make_mutyh_grid.py"
	output:
		"plots/figure_4c.eps"
	shell:
		"""
		python {input.py_script} --strain_vars {input.mgp_mutyh_vars} \
									--out {output}
		"""

rule make_figure_four_d:
	input:
		dumont_xls = "data/SuppTables_concat.xlsx",
		py_script = "py_scripts/compare_mgp_spectra.py"
	output:
		"plots/figure_4d.eps"
	shell:
		"""
		python {input.py_script} --dumont_xls {input.dumont_xls} \
									--out {output}
		"""
		
rule make_figure_four_a:
	input:
		wild_mutyh_vars = "data/wild.mutyh_genotypes.csv",
		py_script = "py_scripts/make_mutyh_grid.py"
	output:
		"plots/figure_4a.eps",
	shell:
		"""
		python {input.py_script} --strain_vars {input.wild_mutyh_vars} \
									--out {output} \
									-is_wild
		"""