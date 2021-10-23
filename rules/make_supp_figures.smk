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
	output:
		"plots/supp_figure_3a.eps"
	shell:
		"""
		Rscript {input.qtl_rscript} -j {input.qtl_json} \
									-p {input.mut_rates} 
		""" 

rule make_supp_figure_three_b:
	input:
		mut_spectra = "csv/tidy_mutation_spectra.csv",
		qtl_rscript = "Rscripts/qtl_mapping_ind_muts.R",
		qtl_geno = "Rqtl_data/bxd.geno.new",
		qtl_json = "Rqtl_data/bxd.json",
		qtl_gmap = "Rqtl_data/bxd.gmap",
		qtl_pmap = "Rqtl_data/bxd.pmap",
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
									-p {input.mut_spectra}
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

rule make_supp_figure_four:
	input:
		py_script = "py_scripts/mutyh_expression_in_tissues.py",
		strain_metadata = "data/bam_names_to_metadata.xlsx",
		mut_spectra = "csv/tidy_mutation_spectra.csv",
		expression_data = expand("data/gene_network_expression/{tissue}_rnaseq.csv", tissue=tissues)
	output:
		"plots/supp_figure_4.eps"
	shell:
		"""
		python {input.py_script} --strain_metadata {input.strain_metadata} \
								 --tidy_spectra {input.mut_spectra} \
								 --out {output}
		"""


rule make_supp_figure_five_a:
	input:
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		cosmic_sig = "data/sigProfiler_SBS_signatures_SBS36.csv",
		py_script = "py_scripts/compare_signatures_OR_COSMIC.py"
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
		py_script = "py_scripts/compare_signatures_OR_COSMIC.py"
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
		py_script = "py_scripts/compare_signatures_BXD68_COSMIC.py",
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
		py_script = "py_scripts/compare_signatures_BXD68_COSMIC.py",
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
		py_script = "py_scripts/compare_signatures_BXD68_TOYKO.py",
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
		py_script = "py_scripts/compare_signatures_BXD68_OR.py",
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
		py_script = "py_scripts/compare_mgp_spectra_3mer.py"
	output:
		"plots/supp_figure_7.eps"
	shell:
		"""
		python {input.py_script} --dumont_xls {input.dumont_xls} \
								 --annotated_singletons {input.annotated_singletons} \
									--out {output}
		"""

rule make_wild_popfiles:
	input:
		py_script = "py_scripts/make_wild_popfile.py",
	output:
		"data/mutyper_input/wild.popfile.{species}.txt"
	shell:
		"""
		python {input.py_script} --out {output} --species {wildcards.species}
		"""

rule mask_ancestral_genome:
	input:
		ancestral_ref = "data/ref/mm10.ancestral.{chrom}.fa",
		exclude = "data/BEDFILEOFEXCLUDEDSITES.amibiguousPolarizing.0based.bed.gz"
	output:
		"data/ref/mm10.ancestral.{chrom}.ambig_polarized_masked.fa"
	shell:
		"""
		bedtools maskfasta -fi {input.ancestral_ref} -bed {input.exclude} -fo {output}
		"""


rule generate_mutyper_sfs_input:
	input: 
		ancestral_ref = "data/ref/mm10.ancestral.{chrom}.ambig_polarized_masked.fa",
		exclude = "data/mm10.seg_dups.simple_repeats.merged.bed.gz",
		bcftools = "/Users/tomsasani/bin/bcftools",
		bgzip = "/Users/tomsasani/bin/bgzip",
	output:
		temp("data/mutyper_input/mutyper_in.{chrom}.vcf.gz")
	shell:
		"""
		{input.bcftools} view -c 1:minor \
							  -T ^{input.exclude} \
							  https://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz \
							  {wildcards.chrom} | \
							  grep -v "\./\." \
							  | \
							  mutyper variants --k 1 \
							  				   --chrom_pos 0 \
											   {input.ancestral_ref} - \
							  | \
							  {input.bgzip} -cf > {output}
		"""

rule tabix_mutyper_input:
	input:
		vcf = "data/mutyper_input/mutyper_in.{chrom}.vcf.gz",
		tabix = "/Users/tomsasani/bin/tabix"
	output:
		"data/mutyper_input/mutyper_in.{chrom}.vcf.gz.tbi"
	shell:
		"{input.tabix} -p vcf {input.vcf}"


rule generate_mutyper_sfs:
	input:
		bcftools = "/Users/tomsasani/bin/bcftools",
		mutyper_in = "data/mutyper_input/mutyper_in.{chrom}.vcf.gz",
		mutyper_in_index = "data/mutyper_input/mutyper_in.{chrom}.vcf.gz.tbi",
		popfile = "data/mutyper_input/wild.popfile.{species}.txt"
	output:
		"data/mutyper_output/ksfs.{chrom}.{species}.txt"
	shell:
		"""
		{input.bcftools} view -S {input.popfile} {input.mutyper_in} | \
		 {input.bcftools} view -c 1:minor | mutyper ksfs - > {output}
		 """


rule plot_wild_sfs:
	input:
		py_script = "py_scripts/plot_ca_sfs.py",
		ksfs_files = expand("data/mutyper_output/ksfs.{chrom}.{species}.txt", 
							species=["Mmd", "Ms", "Mmm", "Mmc",], chrom=["chr4"],)
	output:
		"plots/supp_figure_8.eps"
	shell:
		"""
		python {input.py_script} --wild_sfs {input.ksfs_files} --out {output}
		"""


rule identify_fixed_differences:
	input:
		exclude = "data/mm10.seg_dups.simple_repeats.merged.bed.gz",
		py_script = "py_scripts/identify_fixed_diffs_b6_d2.py",
		ref = "data/ref/mm10.fa"
	output:
		"csv/fixed_diffs/fixed_differences_d2_b6_{chrom}.csv"
	shell:
		"""
		python {input.py_script} --ref {input.ref} \
								 --chrom {wildcards.chrom} \
								 --out {output} \
								 -exclude {input.exclude}
		"""

# rule concatenate_fixed_vars:
# 	input:
# 		expand("csv/fixed_diffs/fixed_differences_d2_b6_{chrom}.csv")#, chrom=["4"])#chrom=list(map(str, range(1, 20))))
# 	output:
# 		"csv/fixed_diffs/fixed_differences_d2_b6_combined.csv"
# 	shell:
# 		"""
# 		cat {input} > {output}
# 		"""

rule count_rare_ca:
	input:
		exclude = "data/mm10.seg_dups.simple_repeats.merged.bed.gz",
		py_script = "py_scripts/count_rare_ca_vars_wild_mice.py",
		ref = "data/ref/mm10.fa",
		fixed_vars = "csv/fixed_diffs/fixed_differences_d2_b6_{chrom}.csv"

	output:
		"csv/rare_ca_muts/lenient/rare_ca_mutations_wild.{chrom}.csv"
	shell:
		"""
		python py_scripts/count_rare_ca_vars_wild_mice.py --ref {input.ref} \
								 --fixed_vars {input.fixed_vars} \
								 --out {output} \
								 --chrom {wildcards.chrom} \
								 -nmer 1 \
								 -min_dp 1 \
								 -min_gq 1 \
								 -exclude {input.exclude}
		"""


rule make_supp_figure_eight_b:
	input:
		counts = expand("csv/rare_ca_muts/lenient/rare_ca_mutations_wild.{chrom}.csv", chrom=list(map(str, range(1, 20)))),
		py_script = "py_scripts/plot_rare_ca_variants.py"
	output:
		"plots/supp_figure_8b.eps"
	shell:
		"""
		python {input.py_script} --variants {input.counts} \
								 --out {output}
		"""
