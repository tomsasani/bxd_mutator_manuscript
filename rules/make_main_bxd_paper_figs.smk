rule make_epoch_sharing_heatmap:
	input:
		py_script = "py_scripts/make_epoch_sharing_heatmap.py",
		annotated_shared_vars = "csv/annotated_shared_vars.csv"
	output:
		"plots/epoch_sharing_heatmap.eps"
	shell:
		"""
		python {input.py_script} --annotated_vars {input.annotated_shared_vars} \
							   --out {output} 
		"""
