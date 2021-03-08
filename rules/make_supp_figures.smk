chroms = list(map(str, range(1, 20)))
chroms = ['chr' + c for c in chroms]

rule make_supp_figure_two: 
	input:
		py_script = "py_scripts/plot_conservation.py",
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		annotated_fixed = "csv/annotated_fixed_vars.csv",
	output:
		"plots/supp_figure_2.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
							   --annotated_fixed {input.annotated_fixed} \
							   --out {output} 
		"""

rule make_supp_figure_three_a:
	input:
		mut_rates = "csv/tidy_mutation_rates.csv",
		qtl_rscript = "Rscripts/qtl_mapping_overall_rate.R",
		qtl_geno = "Rqtl_data/bxd.geno.new",
		qtl_json = "Rqtl_data/bxd.json",
		qtl_gmap = "Rqtl_data/bxd.gmap",
		qtl_pmap = "Rqtl_data/bxd.pmap",
		outpref = "plots/"
	output:
		"plots/supp_figure_3a.eps"
	shell:
		"""
		Rscript {input.qtl_rscript} -j {input.qtl_json} \
									-p {input.mut_rates} \
									-o {input.outpref} 
		""" 

rule make_supp_figure_three_b:
	input:
		mut_spectra = "csv/tidy_mutation_spectra.csv",
		qtl_rscript = "Rscripts/qtl_mapping_ind_muts.R",
		qtl_geno = "Rqtl_data/bxd.geno.new",
		qtl_json = "Rqtl_data/bxd.json",
		qtl_gmap = "Rqtl_data/bxd.gmap",
		qtl_pmap = "Rqtl_data/bxd.pmap",
		outpref = "plots/all_qtl_maps"
	output:
		"plots/all_qtl_maps/supp_figure_3_C.T.eps",
        "plots/all_qtl_maps/supp_figure_3_C.A.eps",
		"plots/all_qtl_maps/supp_figure_3_C.G.eps",
		"plots/all_qtl_maps/supp_figure_3_A.T.eps",
		"plots/all_qtl_maps/supp_figure_3_A.G.eps",
		"plots/all_qtl_maps/supp_figure_3_A.C.eps",
		"plots/all_qtl_maps/supp_figure_3_CpG.TpG.eps" 
	shell:
		"""
		Rscript {input.qtl_rscript} -j {input.qtl_json} \
									-p {input.mut_spectra} \
									-o {input.outpref} 
		""" 

rule make_supp_figure_four_a:
	input:
		py_script = "py_scripts/mutation_rates_by_haplotype.py",
		mut_spectra = "csv/tidy_mutation_spectra.csv"
	output:
		"plots/supp_figure_4a.eps"
	shell:
		"""
		python {input.py_script} --tidy_spectra {input.mut_spectra} \
							   --out {output} 
		"""

rule make_supp_figure_four_b:
	input:
		py_script = "py_scripts/mutation_rates_by_inbreeding_time.py",
		mut_spectra = "csv/tidy_mutation_spectra.csv"
	output:
		"plots/supp_figure_4b.eps"
	shell:
		"""
		python {input.py_script} --tidy_spectra {input.mut_spectra} \
							   --out {output} 
		"""
tissues = ["amygdala", "hematopoietic_stem_cells", "kidney",
		   "liver", "retina", "spleen"]

rule make_supp_figure_five:
	input:
		py_script = "py_scripts/mutyh_expression_in_tissues.py",
		strain_metadata = "data/bam_names_to_metadata.xlsx",
		mut_spectra = "csv/tidy_mutation_spectra.csv",
		expression_data = expand("data/gene_network_expression/{tissue}_rnaseq.csv", tissue=tissues)
	output:
		"plots/supp_figure_5.eps"
	shell:
		"""
		python {input.py_script} --strain_metadata {input.strain_metadata} \
								 --tidy_spectra {input.mut_spectra} \
								 --out {output}
		"""


rule make_supp_figure_six_a:
	input:
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		cosmic_sig = "data/sigProfiler_SBS_signatures_SBS36.csv",
		py_script = "py_scripts/compare_signatures_OR_COSMIC.py"
	output:
		"plots/supp_figure_6a.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
						 	   --cosmic_signature {input.cosmic_sig} \
							   --out {output} \
							   --sig_name SBS36_mm10
		"""

rule make_supp_figure_six_b:
	input:
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		cosmic_sig = "data/sigProfiler_SBS_signatures_SBS18.csv",
		py_script = "py_scripts/compare_signatures_OR_COSMIC.py"
	output:
		"plots/supp_figure_6b.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
						 	   --cosmic_signature {input.cosmic_sig} \
							   --out {output} \
							   --sig_name SBS18_mm10
		"""

rule make_supp_figure_seven_a:
	input:
		py_script = "py_scripts/compare_signatures_BXD68_COSMIC.py",
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		cosmic_sig = "data/sigProfiler_SBS_signatures_SBS36.csv",
	output:
		"plots/supp_figure_7a.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
						 	   --cosmic_signature {input.cosmic_sig} \
							   --out {output} \
							   --sig_name SBS36_mm10

		"""
rule make_supp_figure_seven_b:
	input:
		py_script = "py_scripts/compare_signatures_BXD68_COSMIC.py",
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		cosmic_sig = "data/sigProfiler_SBS_signatures_SBS18.csv",
	output:
		"plots/supp_figure_7b.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
						 	   --cosmic_signature {input.cosmic_sig} \
							   --out {output} \
							   --sig_name SBS18_mm10
		"""

rule make_supp_figure_seven_c:
	input:
		py_script = "py_scripts/compare_signatures_BXD68_TOYKO.py",
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		ohno_data = "data/41598_2014_BFsrep04689_MOESM2_ESM.xls",
	output:
		"plots/supp_figure_7c.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
						 	   --ohno_data {input.ohno_data} \
							   --out {output} 
		"""
rule make_supp_figure_seven_d:
	input:
		py_script = "py_scripts/compare_signatures_BXD68_OR.py",
		annotated_singletons = "csv/annotated_singleton_vars.csv",
	output:
		"plots/supp_figure_7d.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
							   --out {output}
		"""

rule make_supp_figure_eight:
	input:
		dumont_xls = "data/SuppTables_concat.xlsx",
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		py_script = "py_scripts/compare_mgp_spectra_3mer.py"
	output:
		"plots/supp_figure_8.eps"
	shell:
		"""
		python {input.py_script} --dumont_xls {input.dumont_xls} \
								 --annotated_singletons {input.annotated_singletons} \
									--out {output}
		"""