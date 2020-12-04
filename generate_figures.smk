WORKDIR = "/Users/tomsasani/harrislab/bxd_mutator_ms/"

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
		WORKDIR + "plots/supp_figure_2a.eps"

rule annotate_vars:
	input: 
		var_dir = "/Users/tomsasani/harrislab/bxd/csvs/{var_type}/",
		strain_metadata = WORKDIR + "data/bam_names_to_metadata.xlsx",
		strain_callable_bp = WORKDIR + "data/smp2denominator.csv",
		annotate_py = WORKDIR + "scripts/annotate_vars.py"
	output:
		WORKDIR + "csv/annotated_{var_type}.csv"
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
		annotated_singleton = WORKDIR + "csv/annotated_singleton.csv",
		tidy_data_py = WORKDIR + "scripts/make_tidy_data.py"
	output:
		combined_mut_summary = WORKDIR + "csv/tidy_mutation_rates.csv",
		disjoint_mut_summary = WORKDIR + "csv/tidy_mutation_fractions.csv"
	shell:
		"""
		python {input.tidy_data_py} \
			   --annotated_singleton {input.annotated_singleton} \
			   --outdir /Users/tomsasani/harrislab/bxd_mutator_ms/csv/
		"""

rule make_figure_one_abd: 
	input:
		boxplot_py = WORKDIR + "scripts/make_boxplots.py",
		tidy_singleton_combined = WORKDIR + "csv/tidy_mutation_rates.csv",
		tidy_singleton_disjoint = WORKDIR + "csv/tidy_mutation_fractions.csv"
	output:
		WORKDIR + "plots/figure_1a.eps",
		WORKDIR + "plots/figure_1b.eps",
		WORKDIR + "plots/figure_1d.eps"
	shell:
		"""
		python {input.boxplot_py} --tidy_rates {input.tidy_singleton_combined} \
								  --tidy_fracs {input.tidy_singleton_disjoint} \
								  --outdir /Users/tomsasani/harrislab/bxd_mutator_ms/plots/
		"""

rule make_figure_one_c: 
	input:
		cons_py = WORKDIR + "scripts/plot_conservation.py",
		annotated_singleton = WORKDIR + "csv/annotated_singleton.csv",
		annotated_common = WORKDIR + "csv/annotated_common.csv",
	output:
		WORKDIR + "plots/figure_1c.eps"
	shell:
		"""
		python {input.cons_py} --annotated_singletons {input.annotated_singleton} \
								  --annotated_common {input.annotated_common} \
								  --out {output} 
		"""

rule make_figure_two_abc: 
	input:
		tidy_singleton_combined = WORKDIR + "csv/tidy_mutation_rates.csv",
		tidy_singleton_disjoint = WORKDIR + "csv/tidy_mutation_fractions.csv",
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
									-p {input.tidy_singleton_disjoint} \
									-o {input.outdir}
		"""

rule make_figure_three_a:
	input:
		annotated_singleton = WORKDIR + "csv/annotated_singleton.csv",
		comp_py = WORKDIR + "scripts/mutation_comparison.py"
	output:
		WORKDIR + "plots/figure_3a.eps"
	shell:
		"""
		python {input.comp_py} --annotated_singletons {input.annotated_singleton} \
							   --out {output} \
							   -subset_key haplotype_at_qtl \
							   -plot_type heatmap
		"""

rule make_figure_three_b:
	input:
		annotated_singleton = WORKDIR + "csv/annotated_singleton.csv",
		cosmic_sig = WORKDIR + "data/sbs36_cosmic_signatures.csv",
		comp_py = WORKDIR + "scripts/compare_signatures.py"
	output:
		WORKDIR + "plots/figure_3b.eps"
	shell:
		"""
		python {input.comp_py} --annotated_singletons {input.annotated_singleton} \
						 	   --cosmic_signature {input.cosmic_sig} \
							   --out {output}
		"""
