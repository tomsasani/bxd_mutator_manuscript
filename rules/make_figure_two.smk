# set name of working directory
WORKDIR = "/Users/tomsasani/harrislab/bxd_mutator_ms/"

# generate a `.geno` file that R/qtl2 uses to map sample IDs
# to genotypes at each of ~7,000 autosomal markers
rule generate_bxd_geno:
	input:
		gts = WORKDIR + "data/bxd_genotypes_at_rsids.csv",
		py_script = WORKDIR + "scripts/make_bxd_geno_map.py"
	output:
		WORKDIR + "Rqtl_data/bxd.geno.new.updated"
	shell:
		"""
		python {input.py_script} --geno_file {input.gts} \
									--output {output} 
		"""

rule make_figure_two_ab: 
	input:
		mut_spectra = WORKDIR + "csv/tidy_mutation_spectra.csv",
		qtl_rscript = WORKDIR + "Rscripts/qtl_mapping.R",
		qtl_json = WORKDIR + "Rqtl_data/bxd.json",
		#qtl_geno = WORKDIR + "Rqtl_data/bxd.geno.new.updated",
		qtl_geno = WORKDIR + "Rqtl_data/bxd.geno.broman.conv",
		qtl_gmap = WORKDIR + "Rqtl_data/bxd.gmap",
		qtl_pmap = WORKDIR + "Rqtl_data/bxd.pmap",
		outpref = WORKDIR + "plots"
	output:
		WORKDIR + "plots/figure_2a.eps",
		WORKDIR + "plots/figure_2b.eps",
		WORKDIR + "plots/figure_2c.eps"
	shell:
		"""
		Rscript {input.qtl_rscript} -j {input.qtl_json} \
									-p {input.mut_spectra} \
									-o {input.outpref} 
		"""

rule make_figure_two_d:
	input:
		py_script = WORKDIR + "scripts/pca_projection.py",
		mut_spectra = WORKDIR + "csv/tidy_mutation_spectra.csv",
		dumont_xls = "/Users/tomsasani/Downloads/msz026_supp/SuppTables_concat.xlsx"
	output:
		WORKDIR + "plots/figure_2d.eps"
	shell:
		"""
		python {input.py_script} --tidy_spectra {input.mut_spectra} \
							  --dumont_xls {input.dumont_xls} \
							  --out {output}
		"""