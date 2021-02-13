### Discovery of a quantitative trait locus for the germline mutation rate in mice

The code in this repository is sufficient to reproduce the entire manuscript from "top to bottom." This includes everything from
downloading a reference genome to generating supplementary figures. To do so, **there are a few dependencies you'll need**:

#### Software dependencies
Make sure that these are installed prior to running any pipelines!
* [conda](https://docs.conda.io/en/latest/)
* [snakemake](https://snakemake.readthedocs.io/en/stable/)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [bedops](https://bedops.readthedocs.io/en/latest/)
* [tabix](http://www.htslib.org/doc/tabix.html)

All `python` dependencies will be handled by `snakemake` and `conda` when a pipeline is executed. 

#### Directory structure

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
* `env.yaml`
    * `conda` YAML file containing all of the dependencies required to run the pipelines. these dependencies are installed into a `conda` environment "on the fly" when a pipeline is executed

#### Usage

The basic outline of the pipeline is as follows:

1) Identify high-quality singleton mutations from the BXD VCF using `count_mutations.smk`

2) Annotated singleton calls and generate figures using `generate_figures.smk`

##### NOTE!!!

**Identifying singletons is much, much easier if the `count_mutations.smk` pipeline is run on a high-performance computing cluster. The `count_mutations.smk` pipeline involves downloading the BXD VCF (~70 Gbp), the mm10 reference genome (~3 Gbp), and many Gbps of phastCons scores.**

**For this reason, the `data/` directory already contains the files output by the `count_mutations.smk` pipeline for those users that cannot or don't want to run all of the steps in that pipeline.**

To generate all of the figures in the manuscript **using pre-computed files in the `data/` directory**, run the following steps.

```
# enter a directory of your choice (make sure)
# you have at least 1 Gbp of free space here
cd $WORKDIR

# clone the BXD analysis GitHub repository
git clone $name

# create a new conda environment with all of the dependencies
# required for running the pipeline
conda env create --name bxd_analysis --file env.yaml

# activate the new environment
conda activate bxd_analysis

# run the pipeline
snakemake \
        -j 1 \ # number of simultaneous jobs to run
        -s generate_figures.smk # name of the snakemake pipeline
```




