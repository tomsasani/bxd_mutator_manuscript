# set name of working directory
WORKDIR = "/Users/tomsasani/harrislab/bxd_mutator_ms/"

rule make_figure_four_a:
	input:
		mgp_var_info = WORKDIR + "data/mgp.mutyh.csv",
		py_script = WORKDIR + "scripts/mgp_mouse_mutyh.py"
	output:
		WORKDIR + "plots/figure_4a.eps"
	shell:
		"""
		python {input.py_script} --strain_vars {input.mgp_var_info} \
									--out {output}
		"""

rule make_figure_four_bc:
	input:
		dumont_xls = WORKDIR + "data/dumont_singletons.xlsx",
		py_script = WORKDIR + "scripts/compare_mgp_spectra.py"
	output:
		WORKDIR + "plots/figure_4bc.eps"
	shell:
		"""
		python {input.py_script} --dumont_xls {input.dumont_xls} \
									--out {output}
		"""
rule make_figure_four_d:
	input:
		wild_var_info = WORKDIR + "data/wild.mutyh.csv",
		py_script = WORKDIR + "scripts/wild_mouse_mutyh.py"
	output:
		out_a = WORKDIR + "plots/figure_4da.eps",
		out_b = WORKDIR + "plots/figure_4db.eps"
	shell:
		"""
		python {input.py_script} --strain_vars {input.wild_var_info} \
									--out_a {output.out_a} \
									--out_b {output.out_b}
		"""