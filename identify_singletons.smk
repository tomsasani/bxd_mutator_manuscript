from cyvcf2 import VCF

# NOTE: if running `identify_singletons.smk` on a cluster (using
# something like SGE), absolute paths to all files need to be provided

# ---
# SET PATH TO WORKING DIRECTORY HERE
# ---
WORKDIR = '/net/harris/vol1/project/bxd_manuscript_fresh/'

# --- 
# SET PATH TO VCF FILE HERE
# ---
VCFFILE = "/net/harris/vol1/project/bxd_mutators/bxd.2020_10.snp_eff.GRCm38.86.vcf.gz"

# ---
# SET PATH TO BEDOPS WIG2BED BINARY HERE
# ---
WIG2BED = "/net/harris/vol1/home/sasanit/bin/wig2bed"

# make a list of chromosomes
chroms = list(map(str, range(1, 20)))
chroms = ['chr' + c for c in chroms]

# read in VCF and get a list of sample IDs
VCF_H = VCF(VCFFILE)
samples = VCF_H.samples
samples = [s for s in samples if 'BXD' in s]

EXCLUDE = WORKDIR + "data/mm10.seg_dups.simple_repeats.merged.bed.gz"

rule all:
    input:
        expand(WORKDIR + "data/singleton_vars/{chrom}_singleton_vars.exclude.csv", chrom=chroms),
        expand(WORKDIR + "data/fixed_vars/{chrom}_fixed_vars.exclude.csv", chrom=chroms),
        expand(WORKDIR + "data/singleton_vars/{chrom}_singleton_vars.wild_mice.exclude.csv", chrom=chroms),
        WORKDIR + "data/mgp.mutyh_genotypes.csv",
        WORKDIR + "data/wild.mutyh_genotypes.csv",

rule download_reference:
    input:
    output:
        WORKDIR + "data/ref/mm10.fa.gz"
    shell:
        """
        wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz -O {output}
        """

rule unzip_reference:
    input:
        WORKDIR + "data/ref/mm10.fa.gz"
    output:
        WORKDIR + "data/ref/mm10.fa"
    shell:
        """
        gunzip {input}
        """

rule download_conservation:
    input:
    output:
        temp(WORKDIR + "data/conservation_scores/{chrom}.phastCons60way.placental.wigFix.gz")
    shell:
        """
        wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/phastCons60way/placental/{wildcards.chrom}.phastCons60way.placental.wigFix.gz \
                -O {output}
        """

rule unzip_conservation:
    input:
        WORKDIR + "data/conservation_scores/{chrom}.phastCons60way.placental.wigFix.gz"
    output:
        temp(WORKDIR + "data/conservation_scores/{chrom}.phastCons60way.placental.wigFix")
    shell:
        """
        gunzip {input}
        """

rule convert_cons_to_bed:
    input:
        wig2bed_binary = WIG2BED,
        wig = WORKDIR + "data/conservation_scores/{chrom}.phastCons60way.placental.wigFix",
    output:
        temp(WORKDIR + "data/conservation_scores/{chrom}.phastCons60way.placental.bed")
    shell:
        """
        {input.wig2bed_binary} --zero-indexed < {input.wig} > {output}
        """

rule compress_and_index_cons:
    input:
        WORKDIR + "data/conservation_scores/{chrom}.phastCons60way.placental.bed"
    output:
        WORKDIR + "data/conservation_scores/{chrom}.phastCons60way.placental.bed.gz",
        WORKDIR + "data/conservation_scores/{chrom}.phastCons60way.placental.bed.gz.tbi"
    shell:
        """
        bgzip {input}
        tabix -p bed {input}.gz
        """

rule identify_singleton_vars:
    input:
        vcf = VCFFILE,
        ref = WORKDIR + 'data/ref/mm10.fa',
        exclude = EXCLUDE,
        pcons = WORKDIR + "data/conservation_scores/{chrom}.phastCons60way.placental.bed.gz",
        py_script = WORKDIR + "py_scripts/identify_singletons.py",
        haps = expand(WORKDIR + 'data/haplotypes/{sample}_haplotypes.csv', sample=samples),
    output:
        WORKDIR + "data/singleton_vars/{chrom}_singleton_vars.exclude.csv",
    shell:
        """
        python {input.py_script} --vcf {input.vcf} \
                                 --chrom {wildcards.chrom} \
                                 --out {output} \
                                 --haplotypes {input.haps} \
                                 --ref {input.ref} \
                                 --pcons {input.pcons} \
                                 -min_dp 10 \
                                 -min_gq 20 \
                                 -nmer 3 \
                                 -exclude {input.exclude}
        """

rule generate_windows:
    input:
        WORKDIR + "data/mm10.chrom.sizes"
    output:
        WORKDIR + "data/mm10.50kbp.windows"
    shell:
        """
        module load bedtools/2.29.2
        bedtools makewindows -g {input} -w 50000 > {output}
        """

rule identify_fixed_vars:
    input:
        vcf = VCFFILE,
        haps = expand(WORKDIR + 'data/haplotypes/{sample}_haplotypes.csv', sample=samples),
        ref = WORKDIR + 'data/ref/mm10.fa',
        exclude = EXCLUDE,
        pcons = WORKDIR + "data/conservation_scores/{chrom}.phastCons60way.placental.bed.gz",
        intervals = WORKDIR + "data/mm10.50kbp.windows",
        py_script = WORKDIR + "py_scripts/identify_fixed_vars.py"
    output:
        WORKDIR + "data/fixed_vars/{chrom}_fixed_vars.exclude.csv"
    shell:
        """
        python {input.py_script} --vcf {input.vcf} \
                                 --chrom {wildcards.chrom} \
                                 --intervals {input.intervals} \
                                 --haplotypes {input.haps} \
                                 --ref {input.ref} \
                                 --pcons {input.pcons} \
                                 --out {output} \
                                 -min_dp 10 \
                                 -min_gq 20 \
                                 -nmer 3 \
                                 -exclude {input.exclude}
        """

rule identify_singletons_wild_mice:
    input:
        ref = WORKDIR + 'data/ref/mm10.fa',
        exclude = EXCLUDE,
        pcons = WORKDIR + "data/conservation_scores/{chrom}.phastCons60way.placental.bed.gz",
        py_script = WORKDIR + "py_scripts/identify_singletons.wild_mice.py"
    output:
        WORKDIR + "data/singleton_vars/{chrom}_singleton_vars.wild_mice.exclude.csv",
    shell:
        """
        python {input.py_script} --chrom {wildcards.chrom} \
                                 --out {output} \
                                 --ref {input.ref} \
                                 --pcons {input.pcons} \
                                 -min_dp 10 \
                                 -min_gq 20 \
                                 -nmer 3 \
                                 -exclude {input.exclude}
        """

rule find_inf_sites:
    input:
        vcf = VCFFILE,
        exclude = EXCLUDE,
        py_script = WORKDIR + "py_scripts/catalog_inf_sites.py"
    output:
        fh_a = WORKDIR + "data/inf_sites_for_hmm/{chrom}_inf_site_positions.csv",
        fh_b = WORKDIR + "data/inf_sites_for_hmm/{chrom}_inf_site_states.csv"
    shell:
        """
        python {input.py_script} --vcf {input.vcf} \
                                  --chrom {wildcards.chrom} \
                                  --fh_a {output.fh_a} \
                                  --fh_b {output.fh_b} \
                                  -exclude {input.exclude}
        """

rule extract_haplotypes:
    input:
        inf_site_positions = expand(WORKDIR + "data/inf_sites_for_hmm/{chrom}_inf_site_positions.csv", chrom=chroms),
        inf_site_states = expand(WORKDIR + "data/inf_sites_for_hmm/{chrom}_inf_site_states.csv", chrom=chroms),
        py_script = WORKDIR + "py_scripts/extract_haplotypes.py"
    output:
        haps = WORKDIR + 'data/haplotypes/{sample}_haplotypes.csv'
    shell:
        """
        python {input.py_script} \
                --inf_site_positions {input.inf_site_positions} \
                --inf_site_states {input.inf_site_states} \
                --sample {wildcards.sample} \
                --out {output}
        """


rule get_wild_mutyh_alleles:
    input:
        py_script = "py_scripts/mutyh_mutations_in_mice.py"
    output:
        WORKDIR + "data/wild.mutyh_genotypes.csv"
    shell:
        """
        python {input.py_script} --out {output} \
                                 -is_wild
        """

rule get_mgp_mutyh_alleles:
    input:
        py_script = "py_scripts/mutyh_mutations_in_mice.py"
    output:
        WORKDIR + "data/mgp.mutyh_genotypes.csv"
    shell:
        """
        python {input.py_script} --out {output}
        """
