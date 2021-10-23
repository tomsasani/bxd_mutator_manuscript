# pseudo-rule to collect all output figures
figures = ["1", "2", "3", "4"]
				   
rule all:
	input:
		# generate all of the main and supplementary figures
		expand("plots/response_figures/response_figure_{fig_num}.eps", fig_num=figures),
        "data/sigprofiler_outdir/"

rule make_figure_one_two:
	input: 
		strain_metadata = "data/bam_names_to_metadata.xlsx",
		mut_spectra = "csv/tidy_mutation_spectra.csv",
		py_script = "py_scripts/bxd_generation_times.py"
	output:
		"plots/response_figures/response_figure_1.eps",
		"plots/response_figures/response_figure_2.eps"
	shell:
		"""
		python {input.py_script} \
			   --strain_metadata {input.strain_metadata} \
			   --tidy_spectra {input.mut_spectra} \
		"""

rule make_figure_three:
	input:
		py_script = "py_scripts/mutyh_expression_gastrointestinal.py",
		strain_metadata = "data/bam_names_to_metadata.xlsx",
		mut_spectra = "csv/tidy_mutation_spectra.csv",
		expression_data = "data/gene_network_expression/gastrointestinal_rnaseq.csv"
	output:
		"plots/response_figures/response_figure_3.eps"
	shell:
		"""
		python {input.py_script} --strain_metadata {input.strain_metadata} \
								 --tidy_spectra {input.mut_spectra} \
								 --expression_data {input.expression_data} \
								 --out {output}
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

EXCLUDE = "data/mm10.seg_dups.simple_repeats.merged.bed.gz"

rule count_windowed_ca_wild_mice:
	input:
		exclude = EXCLUDE,
		ref = "data/ref/mm10.fa",
		py_script = "py_scripts/identify_windowed_ca.wild_mice.py",
	output:
		"csv/wild.windowed_singletons.csv"
	shell:
		"""
        python {input.py_script} --ref {input.ref} \
                                --region chr4:110000000-120000000 \
                                --out {output} \
                                -exclude {input.exclude} \
								-nmer 1 \
								-min_dp 10 \
								-min_gq 20
		"""

rule make_figure_four:
	input:
		py_script = "py_scripts/wild_ca_distribs.py",
		csv = "csv/wild.windowed_singletons.csv"
	output:
		"plots/response_figures/response_figure_4.eps"
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