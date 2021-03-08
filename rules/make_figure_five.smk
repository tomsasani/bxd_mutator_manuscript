chroms = list(map(str, range(1, 20)))
chroms = ['chr' + c for c in chroms]

rule make_figure_five:
    input:
        py_script = "scripts/compare_spectra_wild.py",
        singleton_vars = expand("data/singleton_vars/{chrom}_singleton_vars.wild_mice.exclude.csv", chrom=chroms),
    output:
        "plots/figure_5a.eps",
        "plots/figure_5b.eps",
        "plots/figure_5c.eps"
    shell:
        """
        python {input.py_script} --wild_singleton_vars {input.singleton_vars}
        """