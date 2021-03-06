rule make_geno_file:
    input:
        strain_metadata = "data/bam_names_to_metadata.xlsx",
        geno = "Rqtl_data/bxd.geno",
        py_script = "py_scripts/orig_to_new_geno.py"
    output:
        "Rqtl_data/bxd.geno.new"
    shell:
        """
        python {input.py_script} --geno {input.geno} \
                                 --strain_metadata {input.strain_metadata} \
                                 --out {output}
        """

rule make_figure_two_a: 
    input:
        mut_spectra = "csv/tidy_mutation_spectra.csv",
        qtl_rscript = "Rscripts/qtl_mapping.R",
        qtl_geno = "Rqtl_data/bxd.geno.new",
        qtl_json = "Rqtl_data/bxd.json",
        qtl_gmap = "Rqtl_data/bxd.gmap",
        qtl_pmap = "Rqtl_data/bxd.pmap",
    output:
        "plots/figure_2a.eps",
        "plots/figure_2a_inset.eps"
    shell:
        """
        Rscript {input.qtl_rscript} -j {input.qtl_json} \
                                    -p {input.mut_spectra}
        """

rule make_figure_two_b:
    input:
        mut_spectra = "csv/tidy_mutation_spectra.csv",
        qtl_rscript = "Rscripts/make_figure_2b.R",
    output:
        "plots/figure_2b.eps",
    shell:
        """
        Rscript {input.qtl_rscript} -s {input.mut_spectra}
        """

rule make_figure_two_c:
    input:
        py_script = "py_scripts/pca_projection.py",
        mut_spectra = "csv/tidy_mutation_spectra.csv",
        dumont_xls = "data/SuppTables_concat.xlsx"
    output:
        "plots/figure_2c.eps"
    shell:
        """
        python {input.py_script} --tidy_spectra {input.mut_spectra} \
                              --dumont_xls {input.dumont_xls} \
                              --out {output}
        """
