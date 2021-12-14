rule make_figure_three_a:
    input:
        annotated_singletons = "csv/annotated_singleton_vars.csv",
        py_script = "py_scripts/mutation_comparison.py"
    output:
        "plots/figure_3a.eps"
    shell:
        """
        python {input.py_script} --annotated_singletons {input.annotated_singletons} \
                               --out {output} \
                               -subset_key haplotype_at_qtl \
                               -plot_type heatmap
        """

rule make_figure_three_b:
    input:
        annotated_singletons = "csv/annotated_singleton_vars.csv",
        ohno_data = "data/41598_2014_BFsrep04689_MOESM2_ESM.xls",
        py_script = "py_scripts/compare_signatures_OR_TOYKO.py"
    output:
        "plots/figure_3b.eps"
    shell:
        """
        python {input.py_script} --annotated_singletons {input.annotated_singletons} \
                                --ohno_data {input.ohno_data} \
                               --out {output}
        """

rule make_figure_three_c:
    input:
        annotated_singletons = "csv/annotated_singleton_vars.csv",
        cosmic_sig = "data/sigProfiler_SBS_signatures_SBS18.csv",
        py_script = "py_scripts/compare_signatures_OR_COSMIC.py"
    output:
        "plots/figure_3c.eps"
    shell:
        """
        python {input.py_script} --annotated_singletons {input.annotated_singletons} \
                                --cosmic_signature {input.cosmic_sig} \
                               --out {output} \
                               --sig_name SBS18_mm10
        """
