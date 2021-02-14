# Discovery of a quantitative trait locus for the germline mutation rate in mice

Thomas A. Sasani, David G. Ashbrook, Abraham A. Palmer, Robert W. Williams, Jonathan K. Pritchard, Kelley Harris

The code in this repository is sufficient to reproduce the entire manuscript from "top to bottom." This includes everything from downloading a reference genome to generating supplementary figures. However, it's also possible to
simply generate the figures in the manuscript using pre-computed data, since the latter require very large input files and quite a bit of time to create.

The basic outline of the pipeline is as follows:

1) Identify high-quality singleton mutations from the BXD VCF using `count_mutations.smk`.

2) Annotate singleton calls and generate figures using `generate_figures.smk`.

>Identifying singletons is much, much easier if the `count_mutations.smk` pipeline is run on a high-performance computing system. The `count_mutations.smk` pipeline involves downloading the BXD VCF (~70 Gbp), the mm10 reference genome (~3 Gbp), and many Gbps of phastCons scores.

>For this reason, the `data/` directory already contains the files output by the `count_mutations.smk` pipeline for those users that cannot or don't want to run all of the steps in that first pipeline.

## Table of Contents

1. [Dependencies](#dependencies)
2. [Directory structure](#directory-structure)
3. [Usage for generating manuscript figures using pre-computed data](#usage-figures)
4. [Usage for generating all raw data](#usage-data)
5. [How to generate raw data on a HPC](#usage-hpc)

## Dependencies
Make sure that these are installed and in your system `$PATH`!

### Required for all pipelines
* [conda](https://docs.conda.io/en/latest/)
* [snakemake](https://snakemake.readthedocs.io/en/stable/)

### Required only if you want to generate singleton data from scratch
* [bedtools (v2.29.2)](https://bedtools.readthedocs.io/en/latest/)
* [bedops (v2.4.38)](https://bedops.readthedocs.io/en/latest/)
* [tabix (v1.10.2-125-g4162046)](http://www.htslib.org/doc/tabix.html)
* [mutyper](https://harrispopgen.github.io/mutyper/install.html)

All other `python` dependencies will be handled by `conda` when a pipeline is executed. 

## Directory structure

* `rules/`
    * individual `snakemake` "rule files" that are imported by the main pipelines
* `py_scripts/`
    * all of the `python` scripts called by `snakemake` rules
* `R_scripts/`
    * all of the `R` scripts called by `snakemake` rules
* `Rqtl_data/`
    * data used by R/qtl2 for QTL mapping
* `data/`
    * raw data files output by the `identify_singletons.smk` pipeline (e.g., singleton calls, HMM-inferred haplotypes, etc.)
* `misc/`
    * miscellaneous data files needed to count or annotate singletons
* `figure_generation.yaml`
    * `conda` YAML file containing all of the dependencies required to generate figures.
* `singleton_calling.yaml`
    * `conda` YAML file containing all of the dependencies required to call singletons.
* `generate_figures.smk`
    * main `snakemake` pipeline that generates main and supplementary figures
* `identify_singletons.smk`
    * main `snakemake` pipeline that identifies singletons using the BXD VCF

## Usage for generating manuscript figures using pre-computed data

```
# enter a directory of your choice (make sure
# you have at least 1 Gbp of free space here)
cd $WORKDIR

# clone the BXD analysis GitHub repository
git clone $name

# create a new conda environment with all of the dependencies
# required for running the pipeline
conda env create --name bxd_analysis --file figure_generation.yaml

# activate the new environment
conda activate bxd_analysis

# run the pipeline
snakemake \
        -j 1 \ # number of simultaneous jobs to run
        -s generate_figures.smk # name of the snakemake pipeline
```

This will produce plots from every main and supplementary figure in the manuscript, accessible in the `plots/` directory.

## Usage for generating all raw data

#### What steps are involved?

1) Download mm10 reference assemblies for each chromosome.
2) Download phastCons 60-way placental mammal WIG files for each chromosome, and convert to BED format.
3) Use an HMM to identify D2 and B6 haplotype tracts in each BXD.
4) Determine the 3-mer nucleotide composition of D2 and B6 haplotypes in each BXD.
5) Identify singletons in each BXD.
6) Identify "fixed" variants in D2 that we'll use for conservation score comparisons.
7) Get genotypes of every BXD at each of ~7,000 markers we'll use for QTL mapping.

To generate all of the **raw data in the manuscript**, it's highly recommended that you run the following on a machine with the ability to run many simultaneous processes, and with at least 100 Gbp of free space in the working directory.

```
# enter a directory of your choice (make sure
# you have at least 100 Gbp of free space here)
cd $WORKDIR

# clone the BXD analysis GitHub repository
git clone $name

# create a new conda environment with all of the dependencies
# required for running the pipeline
conda env create --name singleton_calling --file singleton_calling.yaml

# activate the new environment
conda activate singleton_calling

# ---
# before running the pipeline, use your text editor of choice
# and make sure that the WORKDIR variable at the top of 
# `singleton_calling.smk` points to the current directory you
# `cd`'ed into at the beginning of these steps
# ---

# run the pipeline
snakemake \
        -j 20 \ # number of simultaneous jobs to run
        -s singleton_calling.smk # name of the snakemake pipeline
```

## How to generate raw data on a HPC

The Department of Genome Sciences at UW uses the Sun Grid Engine (SGE) for job management on the cluster. As an example, I've included a config file that can be used in conjunction with the following command to run the `singleton_calling.smk` pipeline on a system with SGE. If you use a different job management system (e.g., SLURM), there is documentation on how to submit jobs using `snakemake` [here](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

```
# run the pipeline using SGE
snakemake \
        -j 20 \
        -s singleton_calling.smk \
        --cluster-config singleton_calling_cluster_config.json \
        --cluster "qsub -l centos={cluster.os} \
                        -l mfree={cluster.memory} \
                        -l h_rt={cluster.time} \
                        -l s_rt={cluster.time} \
                        -o {cluster.erroutdir} \
                        -e {cluster.erroutdir}"
```

After generating raw data, you can run the figure generation pipeline as described above. 