# Discovery of a quantitative trait locus for the germline mutation rate in mice

Thomas A. Sasani, David G. Ashbrook, Abraham A. Palmer, Robert W. Williams, Jonathan K. Pritchard, Kelley Harris

The code in this repository is sufficient to reproduce the entire manuscript from "top to bottom." This includes everything from downloading a reference genome to generating supplementary figures. However, it's also possible to
simply generate the figures in the manuscript using pre-computed data, since the latter require very large input files and quite a bit of time to create.

**The basic outline of the pipeline is as follows:**

1) Identify high-quality singleton mutations from the BXD VCF using `identify_singletons.smk`.

2) Annotate singleton calls and generate figures using `generate_figures.smk`.

**IMPORTANT NOTE:**

>Identifying singletons is much, much easier if the `identify_singletons.smk` pipeline is run on a high-performance computing system. The `identify_singletons.smk` pipeline involves downloading the BXD VCF (~70 Gb), the mm10 reference genome (~1 Gb), and many Gbs of phastCons scores. It also involves hundreds of individual steps, so it will finish much more quickly if those steps are run in parallel.

>If you want to generate all raw data from scratch, see the section entitled [Usage for generating all raw data](#usage-data). 

For this reason, the `data/` directory already contains the files output by the `identify_singletons.smk` pipeline for those users that cannot or don't want to run all of the steps in that first pipeline.

## Table of Contents

1. [Dependencies](#dependencies)
2. [Directory structure](#directory-structure)
3. [Usage for generating manuscript figures using precomputed data](#usage-for-generating-manuscript-figures-using-precomputed-data)
4. [Usage for generating all raw data](#usage-for-generating-all-raw-data)
5. [How to generate raw data on a HPC](#how-to-generate-raw-data-on-a-hpc)

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
```
.
|__rules                        # individual `snakemake` "rule files" that are imported by the main pipelines
|__py_scripts                   # all of the `python` scripts called by `snakemake` rules
|__R_scripts                    # all of the `R` scripts called by `snakemake` rules
|__Rqtl_data                    # data used by R/qtl2 for QTL mapping
|__data                         # raw data files output by the `identify_singletons.smk` pipeline
|__misc                         # miscellaneous data files needed to count or annotate singletons
|__figure_generation.yaml       # YAML file containing all of the dependencies required to generate figures
|__singleton_calling.yaml       # YAML file containing all of the dependencies required to call singletons
|__generate_figures.smk         # main `snakemake` pipeline that generates main and supplementary figures
|__identify_singletons.smk      # main `snakemake` pipeline that identifies singletons using the BXD VCF
```
## Usage for generating manuscript figures using precomputed data

If you'd rather not spend the time (and compute) needed to generate all of the singleton data from scratch, I've included the output of `identify_singletons.smk` in the `data/` directory. These files include raw singleton calls and HMM-derived haplotype tracts in each BXD, among others. The command-line instructions below will reproduce every figure in the manuscript using these data.

#### What steps are involved?

1) Annotate singletons and fixed variants with various metadata.
2) Construct "tidy" dataframes containing summary information about mutation rates and spectra in each BXD.
3) Perform QTL scans for phenotypes related to the mutation rate.
4) Make plots that constitute all figures.

In a number of analyses, we compare mutation spectra between the BXD singletons and previously published datasets. **Before moving on to the command-line instructions below, first obtain these third-party datasets as follows:** 

* singleton data from [Dumont 2019](https://academic.oup.com/mbe/article/36/5/865/5315518)
    * download ZIP file from [Supplementary data](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/mbe/36/5/10.1093_molbev_msz026/2/msz026_supp.zip?Expires=1616302324&Signature=u8neUFiV~0aBABDNG-ZPeMwd~usDZPmIO5TVjLHqKVcjHXrUWBm7MnR1ZJpSkMmDmQhMGrcdK~G7hySKLp79xgpQnj-SCFD09Hj7e9uCi9oYvVT-guMav1JY6qEMzSCubzlChpHfItUKJt15lXbxmuT2FxTibIs2gSrXvHCexmwGLxQYCoIAZJHY1nOjOfSDDlIejE-aGrPFozB86PXTZz~uM9JuAnmfZ5wmARxwuEzOHMfYZWh7WnWzeXEaNwKqYzrYHYDhej5sq~LSOfsQTzSPI-nrtn~KOV7x9ckk0RzqJ0kIhmF0uBLMkF8grDXTTRTHEjYV1dBALvxO1ZMVmA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA) associated with the manuscript
    * uncompress file and move `SuppTables_concat.xlsx` into the `data` directory of the cloned version of this repository
* *de novo* germline mutation data from [Ohno et al. 2014](https://www.nature.com/articles/srep04689)
    * download [Supplementary Data 1](https://static-content.springer.com/esm/art%3A10.1038%2Fsrep04689/MediaObjects/41598_2014_BFsrep04689_MOESM2_ESM.xls)
    * move `41598_2014_BFsrep04689_MOESM2_ESM.xls` into the `data` directory of the cloned version of this repository
* mutation signatures from [COSMIC](https://cancer.sanger.ac.uk/cosmic)
    * download [SBS36](https://cancer.sanger.ac.uk/sigs-assets-20/SBS_vignettes/sigProfiler_SBS_signatures_SBS36.csv)
    * download [SBS18](https://cancer.sanger.ac.uk/sigs-assets-20/SBS_vignettes/sigProfiler_SBS_signatures_SBS18.csv)
    * move these two files (`sigProfiler_SBS_signatures_SBS36.csv` and `sigProfiler_SBS_signatures_SBS18.csv`) into the `data` directory of the cloned version of this repository

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

This will produce plots from every main and supplementary figure in the manuscript, which will be added to a new `plots/` directory.

## Usage for generating all raw data

#### What steps are involved?

1) Download mm10 reference assembly.
2) Download phastCons 60-way placental mammal WIG files for each chromosome, and convert to BED format.
3) Use an HMM to identify D2 and B6 haplotype tracts in each BXD.
4) Determine the 3-mer nucleotide composition of D2 and B6 haplotypes in each BXD.
5) Identify singletons in each BXD.
6) Identify "fixed" variants in D2 that we'll use for conservation score comparisons.

To generate all of the **raw data in the manuscript**, it's highly recommended that you run the following on a machine with the ability to run many simultaneous processes, and **with at least 100 Gbp of free space in the working directory.**

Step #4 above is by far the most time-consuming -- the method used to get 3-mer nucleotide compositions in haplotype in each sample on each chromosome is slow and generates a unique file for every combination of sample and chromosome (almost 3,000 individual files). Since the `data/` directory already contains precomputed versions of these files, I've removed the Snakemake rule that generates the files from `identify_singletons.smk`. If you want to recreate these files, just uncomment the part of the `rule all:` input that generates the nucleotide composition files and re-run the pipeline.

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