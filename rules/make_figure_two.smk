
# generate a `.geno` file that R/qtl2 uses to map sample IDs
# to genotypes at each of ~7,000 autosomal markers
rule generate_bxd_geno:
	input:
		gts = "data/bxd_genotypes_at_markers.csv",
		py_script = "scripts/make_bxd_geno_map.py"
	output:
		"Rqtl_data/bxd.geno.from_vcf"
	shell:
		"""
		python {input.py_script} --geno_file {input.gts} \
									--output {output} 
		"""

rule make_figure_two_ab: 
	input:
		mut_spectra = "csv/tidy_mutation_spectra.csv",
		qtl_rscript = "Rscripts/qtl_mapping.R",
		qtl_json = "Rqtl_data/bxd.json",
		qtl_geno = "Rqtl_data/bxd.geno.from_vcf",
		qtl_gmap = "Rqtl_data/bxd.gmap",
		qtl_pmap = "Rqtl_data/bxd.pmap",
		outpref = "plots"
	output:
		"plots/figure_2a.eps",
		"plots/figure_2b.eps",
		"plots/figure_2c.eps"
	shell:
		"""
		Rscript {input.qtl_rscript} -j {input.qtl_json} \
									-p {input.mut_spectra} \
									-o {input.outpref} 
		"""

rule make_figure_two_d:
	input:
		py_script = "scripts/pca_projection.py",
		mut_spectra = "csv/tidy_mutation_spectra.csv",
		dumont_xls = "data/SuppTables_concat.xlsx"
	output:
		"plots/figure_2d.eps"
	shell:
		"""
		python {input.py_script} --tidy_spectra {input.mut_spectra} \
							  --dumont_xls {input.dumont_xls} \
							  --out {output}
		"""