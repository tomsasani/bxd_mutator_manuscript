samples, = glob_wildcards("data/nucleotide_composition/{sample}_nucleotide_composition.csv")

chroms = list(map(str, range(1, 20)))
chroms = ['chr' + c for c in chroms]

rule make_supp_figure_two:
	input:
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		py_script = "scripts/mutation_comparison.py",
		nuc_comp = expand("data/nucleotide_composition/{sample}_nucleotide_composition.csv", sample=samples)
	output:
		"plots/supp_figure_2.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
							   --out {output} \
							   -nmers_for_normalization {input.nuc_comp} \
							   -subset_key haplotype \
							   -plot_type scatter
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
		py_script = "scripts/mutation_rates_by_haplotype.py",
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
		py_script = "scripts/mutation_rates_by_inbreeding_time.py",
		mut_spectra = "csv/tidy_mutation_spectra.csv"
	output:
		"plots/supp_figure_4b.eps"
	shell:
		"""
		python {input.py_script} --tidy_spectra {input.mut_spectra} \
							   --out {output} 
		"""

rule make_supp_figure_five_a:
	input:
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		cosmic_sig = "data/sigProfiler_SBS_signatures_SBS36.csv",
		py_script = "scripts/compare_signatures_OR_COSMIC.py"
	output:
		"plots/supp_figure_5a.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
						 	   --cosmic_signature {input.cosmic_sig} \
							   --out {output} \
							   --sig_name SBS36_mm10
		"""

rule make_supp_figure_five_b:
	input:
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		cosmic_sig = "data/sigProfiler_SBS_signatures_SBS18.csv",
		py_script = "scripts/compare_signatures_OR_COSMIC.py"
	output:
		"plots/supp_figure_5b.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
						 	   --cosmic_signature {input.cosmic_sig} \
							   --out {output} \
							   --sig_name SBS18_mm10
		"""

rule make_supp_figure_six_a:
	input:
		py_script = "scripts/compare_signatures_BXD68_COSMIC.py",
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		cosmic_sig = "data/sigProfiler_SBS_signatures_SBS36.csv",
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
		py_script = "scripts/compare_signatures_BXD68_COSMIC.py",
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		cosmic_sig = "data/sigProfiler_SBS_signatures_SBS18.csv",
	output:
		"plots/supp_figure_6b.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
						 	   --cosmic_signature {input.cosmic_sig} \
							   --out {output} \
							   --sig_name SBS18_mm10
		"""

rule make_supp_figure_six_c:
	input:
		py_script = "scripts/compare_signatures_BXD68_TOYKO.py",
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		ohno_data = "data/41598_2014_BFsrep04689_MOESM2_ESM.xls",
	output:
		"plots/supp_figure_6c.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
						 	   --ohno_data {input.ohno_data} \
							   --out {output} 
		"""
rule make_supp_figure_six_d:
	input:
		py_script = "scripts/compare_signatures_BXD68_OR.py",
		annotated_singletons = "csv/annotated_singleton_vars.csv",
	output:
		"plots/supp_figure_6d.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
							   --out {output}
		"""

rule make_supp_figure_seven:
	input:
		dumont_xls = "data/SuppTables_concat.xlsx",
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		py_script = "scripts/compare_mgp_spectra_3mer.py"
	output:
		"plots/supp_figure_7.eps"
	shell:
		"""
		python {input.py_script} --dumont_xls {input.dumont_xls} \
								 --annotated_singletons {input.annotated_singletons} \
									--out {output}
		"""

rule make_supp_figure_eight:
	input:
		wild_singletons = expand("data/wild_singleton_vars/{chrom}_singleton_spectrum.csv", chrom=chroms),
		py_script = "scripts/pca_projection.wild_mice.py"
	output:
		"plots/supp_figure_8.eps"
	shell:
		"""
		python {input.py_script}  --wild_singleton_vars {input.wild_singletons} \
									--out {output}
		"""
