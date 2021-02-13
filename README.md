## Discovery of a quantitative trait locus for the germline mutation rate in mice

The code in this repository is sufficient to reproduce the entire manuscript from "top to bottom." This includes everything from
downloading a reference genome to generating supplementary figures. To do so, **there are a few dependencies you'll need**:

### Software dependencies
Make sure that these are installed prior to running any pipelines!

#### Only required if you want to generate figures using pre-computed singleton data
* [conda](https://docs.conda.io/en/latest/)
* [snakemake](https://snakemake.readthedocs.io/en/stable/)

#### Required if you want to generate singleton data from scratch
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [bedops](https://bedops.readthedocs.io/en/latest/)
* [tabix](http://www.htslib.org/doc/tabix.html)

All `python` dependencies will be handled by `snakemake` and `conda` when a pipeline is executed. 

### Directory structure

* `rules/`
    * individual `snakemake` "rule files" that are imported by the main pipelines
* `py_scripts/`
    * all of the `python` scripts called by `snakemake` rules
* `Rscripts/`
    * all of the `R` scripts called by `snakemake` rules
* `Rqtl_data/`
    * data used by R/qtl2 for QTL mapping
* `data/`
    * raw data files output by the `count_mutations.smk` pipeline (e.g., singleton calls, HMM-inferred haplotypes, etc.)
* `misc/`
    * miscellaneous data files needed to count singletons
* `figure_generation.yaml`
    * `conda` YAML file containing all of the dependencies required to generate figures.
* `singleton_calling.yaml`
    * `conda` YAML file containing all of the dependencies required to call singletons.

### Usage

The basic outline of the pipeline is as follows:

1) Identify high-quality singleton mutations from the BXD VCF using `count_mutations.smk`

2) Annotate singleton calls and generate figures using `generate_figures.smk`

##### NOTE!

**Identifying singletons is much, much easier if the `count_mutations.smk` pipeline is run on a high-performance computing system. The `count_mutations.smk` pipeline involves downloading the BXD VCF (~70 Gbp), the mm10 reference genome (~3 Gbp), and many Gbps of phastCons scores.**

**For this reason, the `data/` directory already contains the files output by the `count_mutations.smk` pipeline for those users that cannot or don't want to run all of the steps in that first pipeline.**

To generate all of the figures in the manuscript **using pre-computed files in the `data/` directory**, run the following steps.

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