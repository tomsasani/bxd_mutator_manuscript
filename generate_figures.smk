WORKDIR = "/Users/tomsasani/harrislab/bxd_mutator_ms/"

CONDA_YAML="/Users/tomsasani/harrislab/bxd_mutator_ms/env.yaml"

rule all:
	input:
		WORKDIR + "plots/figure_1a.eps",
		WORKDIR + "plots/figure_1b.eps",
		WORKDIR + "plots/figure_1c.eps",
		WORKDIR + "plots/figure_1d.eps",
		WORKDIR + "plots/figure_2a.eps",
		WORKDIR + "plots/figure_2b.eps",
		WORKDIR + "plots/figure_2c.eps",
		WORKDIR + "plots/figure_3a.eps",
		WORKDIR + "plots/figure_3b.eps",
		WORKDIR + "plots/supp_figure_1.eps",
		WORKDIR + "plots/supp_figure_2a.eps",
#		WORKDIR + "plots/supp_figure_3.eps",
		WORKDIR + "plots/epoch_heatmap.eps"

rule annotate_vars:
	input: 
		var_dir = WORKDIR + "data/{var_type}/",
		strain_metadata = WORKDIR + "data/bam_names_to_metadata.xlsx",
		strain_callable_bp = WORKDIR + "data/smp2denominator.csv",
		annotate_py = WORKDIR + "scripts/annotate_vars.py"
	output:
		WORKDIR + "csv/annotated_{var_type}.csv"
	conda:
		CONDA_YAML
	shell:
		"""
		python {input.annotate_py} \
			   --var_dir {input.var_dir} \
			   --strain_metadata {input.strain_metadata} \
			   --callable_bp {input.strain_callable_bp} \
			   --out {output}
		"""
rule generate_tidy_data:
	input:
		annotated_singletons = WORKDIR + "csv/annotated_singletons.csv",
		tidy_data_py = WORKDIR + "scripts/make_tidy_data.py"
	output:
		mut_rates = WORKDIR + "csv/tidy_mutation_rates.csv",
		mut_spectra = WORKDIR + "csv/tidy_mutation_spectra.csv"
	shell:
		"""
		python {input.tidy_data_py} \
			   --annotated_singleton {input.annotated_singletons} \
			   --out_pref /Users/tomsasani/harrislab/bxd_mutator_ms/csv/
		"""

rule make_figure_one_abd: 
	input:
		boxplot_py = WORKDIR + "scripts/make_boxplots.py",
		mut_rates = WORKDIR + "csv/tidy_mutation_rates.csv",
		mut_spectra = WORKDIR + "csv/tidy_mutation_spectra.csv"
	output:
		WORKDIR + "plots/figure_1a.eps",
		WORKDIR + "plots/figure_1b.eps",
		WORKDIR + "plots/figure_1d.eps"
	shell:
		"""
		python {input.boxplot_py} --tidy_rates {input.mut_rates} \
								  --tidy_spectra {input.mut_spectra} \
								  --outdir /Users/tomsasani/harrislab/bxd_mutator_ms/plots/
		"""

rule make_figure_one_c: 
	input:
		cons_py = WORKDIR + "scripts/plot_conservation.py",
		annotated_singletons = WORKDIR + "csv/annotated_singletons.csv",
		annotated_common = WORKDIR + "csv/annotated_common.csv",
	output:
		WORKDIR + "plots/figure_1c.eps"
	shell:
		"""
		python {input.cons_py} --annotated_singletons {input.annotated_singletons} \
							   --annotated_common {input.annotated_common} \
							   --out {output} 
		"""

rule make_figure_two_abc: 
	input:
		mut_spectra = WORKDIR + "csv/tidy_mutation_spectra.csv",
		qtl_rscript = WORKDIR + "Rscripts/qtl_mapping.R",
		qtl_json = WORKDIR + "Rqtl_data/bxd.json",
		qtl_geno = WORKDIR + "Rqtl_data/bxd.geno.new.updated",
		qtl_gmap = WORKDIR + "Rqtl_data/bxd.gmap",
		qtl_pmap = WORKDIR + "Rqtl_data/bxd.pmap",
		outdir = WORKDIR + "plots/"
	output:
		WORKDIR + "plots/figure_2a.eps",
		WORKDIR + "plots/figure_2b.eps",
		WORKDIR + "plots/figure_2c.eps",
		WORKDIR + "plots/supp_figure_2a.eps"
	shell:
		"""
		Rscript {input.qtl_rscript} -j {input.qtl_json} \
									-p {input.mut_spectra} \
									-o {input.outdir}
		"""

rule make_figure_three_a:
	input:
		annotated_singletons = WORKDIR + "csv/annotated_singletons.csv",
		comp_py = WORKDIR + "scripts/mutation_comparison.py"
	output:
		WORKDIR + "plots/figure_3a.eps"
	shell:
		"""
		python {input.comp_py} --annotated_singletons {input.annotated_singletons} \
							   --out {output} \
							   -subset_key haplotype_at_qtl \
							   -plot_type heatmap
		"""

rule make_figure_three_b:
	input:
		annotated_singletons = WORKDIR + "csv/annotated_singletons.csv",
		cosmic_sig = WORKDIR + "data/sbs36_cosmic_signatures.csv",
		comp_py = WORKDIR + "scripts/compare_signatures.py"
	output:
		WORKDIR + "plots/figure_3b.eps"
	shell:
		"""
		python {input.comp_py} --annotated_singletons {input.annotated_singletons} \
						 	   --cosmic_signature {input.cosmic_sig} \
							   --out {output}
		"""

rule make_supp_figure_one:
	input:
		annotated_singletons = WORKDIR + "csv/annotated_singletons.csv",
		comp_py = WORKDIR + "scripts/mutation_comparison.py"
	output:
		WORKDIR + "plots/supp_figure_1.eps"
	shell:
		"""
		python {input.comp_py} --annotated_singletons {input.annotated_singletons} \
							   --out {output} \
							   -subset_key haplotype \
							   -plot_type scatter
		"""

rule make_epoch_sharing_fig:
	input:
		annotated_vars = WORKDIR + "csv/annotated_bxd_private.csv",
		heatmap_py = WORKDIR + "scripts/make_epoch_sharing_heatmap.py"
	output:
		WORKDIR + "plots/epoch_heatmap.eps"
	shell:
		"""
		python {input.heatmap_py} --annotated_vars {input.annotated_vars} \
							   	  --out {output}
		"""
