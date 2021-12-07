chroms = list(map(str, range(1, 20)))
chroms = ['chr' + c for c in chroms]

EXCLUDE = "data/mm10.seg_dups.simple_repeats.merged.bed.gz"

rule make_supp_figure_two:
	input: 
		strain_metadata = "data/bam_names_to_metadata.xlsx",
		mut_spectra = "csv/tidy_mutation_spectra.csv",
		py_script = "py_scripts/bxd_generation_times.py"
	output:
		"plots/supp_figure_2a.eps",
		"plots/supp_figure_2b.eps"
	shell:
		"""
		python {input.py_script} \
			   --strain_metadata {input.strain_metadata} \
			   --tidy_spectra {input.mut_spectra} \
		"""

rule make_supp_figure_three: 
	input:
		py_script = "py_scripts/plot_conservation.py",
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		annotated_fixed = "csv/annotated_fixed_vars.csv",
	output:
		"plots/supp_figure_3.eps"
	shell:
		"""
		python {input.py_script} --annotated_singletons {input.annotated_singletons} \
							   --annotated_fixed {input.annotated_fixed} \
							   --out {output} 
		"""

rule make_supp_figure_four_a:
	input:
		mut_rates = "csv/tidy_mutation_rates.csv",
		qtl_rscript = "Rscripts/qtl_mapping_overall_rate.R",
		qtl_geno = "Rqtl_data/bxd.geno.new",
		qtl_json = "Rqtl_data/bxd.json",
		qtl_gmap = "Rqtl_data/bxd.gmap",
		qtl_pmap = "Rqtl_data/bxd.pmap",
	output:
		"plots/supp_figure_4a.eps"
	shell:
		"""
		Rscript {input.qtl_rscript} -j {input.qtl_json} \
									-p {input.mut_rates} 
		""" 

rule make_supp_figure_four_b:
	input:
		mut_spectra = "csv/tidy_mutation_spectra.csv",
		qtl_rscript = "Rscripts/qtl_mapping_ind_muts.R",
		qtl_geno = "Rqtl_data/bxd.geno.new",
		qtl_json = "Rqtl_data/bxd.json",
		qtl_gmap = "Rqtl_data/bxd.gmap",
		qtl_pmap = "Rqtl_data/bxd.pmap",
	output:
		"plots/all_qtl_maps/supp_figure_4_C.T.eps",
		"plots/all_qtl_maps/supp_figure_4_C.A.eps",
		"plots/all_qtl_maps/supp_figure_4_C.G.eps",
		"plots/all_qtl_maps/supp_figure_4_A.T.eps",
		"plots/all_qtl_maps/supp_figure_4_A.G.eps",
		"plots/all_qtl_maps/supp_figure_4_A.C.eps",
		"plots/all_qtl_maps/supp_figure_4_CpG.TpG.eps" 
	shell:
		"""
		Rscript {input.qtl_rscript} -j {input.qtl_json} \
									-p {input.mut_spectra}
		""" 

rule make_supp_figure_five:
	input:
		dumont_xls = "data/SuppTables_concat.xlsx",
		annotated_singletons = "csv/annotated_singleton_vars.csv",
		py_script = "py_scripts/compare_mgp_spectra_3mer.py"
	output:
		"plots/supp_figure_5.eps"
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

rule unzip_ancestral_genome:
	input:
		"data/ref/mm10.ancestral.{chrom}.fa.gz"
	output:
		"data/ref/mm10.ancestral.{chrom}.fa"
	shell:
		"""
		gunzip {input}
		"""

rule mask_ancestral_genome:
	input:
		ancestral_ref = "data/ref/mm10.ancestral.{chrom}.fa",
		exclude = "data/unpolarizable_sites.mm10.wild.bed.gz"
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


rule make_supp_figure_six:
	input:
		py_script = "py_scripts/plot_ca_sfs.py",
		ksfs_files = expand("data/mutyper_output/ksfs.{chrom}.{species}.txt", 
							species=["Mmd", "Ms", "Mmm", "Mmc",], chrom=["chr4"],)
	output:
		"plots/supp_figure_6.eps"
	shell:
		"""
		python {input.py_script} --wild_sfs {input.ksfs_files} --out {output}
		"""

rule get_sigprofiler_spectra:
	input:
		py_script = "py_scripts/make_spectra_file_for_sig_profiler.py",
		singleton_vars = "csv/annotated_singleton_vars.csv",
	output:
		"csv/matrix4sigprofiler.txt"
	shell:
		"""
		python {input.py_script} --singleton_vars {input.singleton_vars} \
								 --out {output}
		"""

rule initialize_sigprofiler:
	input:
	output: "data/sigprof_initialized_confirmation.txt"
	run:
		from SigProfilerMatrixGenerator import install as genInstall
		genInstall.install('mm10')

		shell("echo done >> {output}")
	
rule run_sigprofiler:
	input:
		py_script = "py_scripts/run_sigprofiler.py",
		spectra = "csv/matrix4sigprofiler.txt",
		conf = "data/sigprof_initialized_confirmation.txt"
	output:
		"data/sigprofiler_outdir/"
	shell:
		"""
		python {input.py_script} --spectra {input.spectra} \
								 --outdir {output}
		"""

rule count_windowed_ca_wild_mice:
	input:
		exclude = EXCLUDE,
		ref = "data/ref/mm10.fa",
		py_script = "py_scripts/identify_windowed_ca.wild_mice.py",
	output:
		"data/wild.windowed_singletons.csv"
	shell:
		"""
		python {input.py_script} --ref {input.ref} \
								--region chr4:114800000-118300000 \
								--out {output} \
								-exclude {input.exclude} \
								-nmer 1 \
								-min_dp 10 \
								-min_gq 20
		"""

rule make_supp_figure_seven:
	input:
		py_script = "py_scripts/compare_spectra_wild.py",
		singleton_vars = expand("data/singleton_vars/{chrom}_singleton_vars.wild_mice.exclude.csv", chrom=chroms),
	output:
		"plots/supp_figure_7a.eps",
		"plots/supp_figure_7b.eps",
		"plots/supp_figure_7c.eps"
	shell:
		"""
		python {input.py_script} --wild_singleton_vars {input.singleton_vars}
		"""

rule make_supp_figure_eight:
	input:
		py_script = "py_scripts/quantify_variability_cast_dom.py",
		csv = "data/wild.windowed_singletons.csv"
	output:
		"plots/supp_figure_8.eps"
	shell:
		"""
		python {input.py_script} --windowed_vars {input.csv} \
								 --out {output} 
		"""

rule intersect_fixed_svs:
	input:
		svs = "data/sv_files/structural_variants.fixed_d2_b6_differences.mgp.bed",
		gtf = "data/sv_files/gencode.mm10.genes.gtf"
	output:
		"data/sv_files/exons_overlapping_svs.bed"
	shell:
		"""
		grep exon {input.gtf} | \
		bedtools intersect -a - -b {input.svs} | \
		cut -f 9 | cut -d ';' -f 2 | cut -d '"' -f 2 > {output}
		"""
