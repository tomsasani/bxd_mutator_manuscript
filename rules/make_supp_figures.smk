# set name of working directory
WORKDIR = "/Users/tomsasani/harrislab/bxd_mutator_ms/"

rule make_supp_figure_one:
	input:
		annotated_singletons = WORKDIR + "csv/annotated_singletons.csv",
		py_script = WORKDIR + "scripts/mutation_comparison.py"
	output:
		WORKDIR + "plots/supp_figure_1.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
							   --out {output} \
							   -subset_key haplotype \
							   -plot_type scatter
		"""

rule make_supp_figure_three:
	input:
		mut_spectra = WORKDIR + "csv/tidy_mutation_spectra.csv",
		qtl_rscript = WORKDIR + "Rscripts/qtl_mapping_ind_muts.R",
		qtl_json = WORKDIR + "Rqtl_data/bxd.json",
		#qtl_geno = WORKDIR + "Rqtl_data/bxd.geno.new.updated",
		qtl_geno = WORKDIR + "Rqtl_data/bxd.geno.broman.conv",
		qtl_gmap = WORKDIR + "Rqtl_data/bxd.gmap",
		qtl_pmap = WORKDIR + "Rqtl_data/bxd.pmap",
		outpref = WORKDIR + "plots/all_qtl_maps"
	output:
		WORKDIR + "plots/all_qtl_maps/supp_figure_3_C.T.eps",
        WORKDIR + "plots/all_qtl_maps/supp_figure_3_C.A.eps",

		WORKDIR + "plots/all_qtl_maps/supp_figure_3_C.G.eps",
		WORKDIR + "plots/all_qtl_maps/supp_figure_3_A.T.eps",
		WORKDIR + "plots/all_qtl_maps/supp_figure_3_A.G.eps",
		WORKDIR + "plots/all_qtl_maps/supp_figure_3_A.C.eps",
		WORKDIR + "plots/all_qtl_maps/supp_figure_3_CpG.TpG.eps" 
	shell:
		"""
		Rscript {input.qtl_rscript} -j {input.qtl_json} \
									-p {input.mut_spectra} \
									-o {input.outpref} 
		""" 

rule make_supp_figure_four:
	input:
		py_script = WORKDIR + "scripts/mutation_rates_by_haplotype.py",
		mut_spectra = WORKDIR + "csv/tidy_mutation_spectra.csv"
	output:
		WORKDIR + "plots/supp_figure_4.eps"
	shell:
		"""
		python {input.py_script} --tidy_spectra {input.mut_spectra} \
							   --out {output} 
		"""

rule make_supp_figure_five_a:
	input:
		annotated_singletons = WORKDIR + "csv/annotated_singletons.csv",
		cosmic_sig = WORKDIR + "data/sbs36_cosmic_signatures.csv",
		py_script = WORKDIR + "scripts/compare_signatures_regression.py"
	output:
		WORKDIR + "plots/supp_figure_5a.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
						 	   --cosmic_signature {input.cosmic_sig} \
							   --out {output} \
							   --sig_name SBS36_mm10
		"""

rule make_supp_figure_five_b:
	input:
		annotated_singletons = WORKDIR + "csv/annotated_singletons.csv",
		cosmic_sig = WORKDIR + "data/sbs18_cosmic_signatures.csv",
		py_script = WORKDIR + "scripts/compare_signatures_regression.py"
	output:
		WORKDIR + "plots/supp_figure_5b.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
						 	   --cosmic_signature {input.cosmic_sig} \
							   --out {output} \
							   --sig_name SBS18_mm10
		"""

rule make_supp_figure_six_a:
	input:
		py_script = WORKDIR + "scripts/bxd_68_vs_cosmic.py",
		annotated_singletons = WORKDIR + "csv/annotated_singletons.csv",
		cosmic_sig = WORKDIR + "data/sbs36_cosmic_signatures.csv",
	output:
		WORKDIR + "plots/supp_figure_6a.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
						 	   --cosmic_signature {input.cosmic_sig} \
							   --out {output}
		"""
rule make_supp_figure_six_b:
	input:
		py_script = WORKDIR + "scripts/bxd_68_vs_cosmic.py",
		annotated_singletons = WORKDIR + "csv/annotated_singletons.csv",
		cosmic_sig = WORKDIR + "data/sbs18_cosmic_signatures.csv",
	output:
		WORKDIR + "plots/supp_figure_6b.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
						 	   --cosmic_signature {input.cosmic_sig} \
							   --out {output}
		"""

rule make_supp_figure_seven:
	input:
		script_py = WORKDIR + "scripts/mutation_rates_by_inbreeding_time.py",
		mut_spectra = WORKDIR + "csv/tidy_mutation_spectra.csv"
	output:
		WORKDIR + "plots/supp_figure_7.eps"
	shell:
		"""
		python {input.script_py} --tidy_spectra {input.mut_spectra} \
							   --out {output} 
		"""
