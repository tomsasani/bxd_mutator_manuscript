# A natural mutator allele shapes mutation spectrum variation in mice

*Thomas A. Sasani, David G. Ashbrook, Annabel C. Beichman, Lu Lu, Abraham A. Palmer, Robert W. Williams, Jonathan K. Pritchard, Kelley Harris*

:mouse: --> :dna: --> :bar_chart:

## Brief overview of the manuscript and associated code

We recently published a [manuscript](https://www.biorxiv.org/content/10.1101/2021.03.12.435196v1) describing our discovery of a genetic modifier of the germline mutation rate in mice. As part of this analysis, we identified *de novo* germline mutations in approximately 100 recombinant inbred mice, known as the [BXD family](https://www.cell.com/cell-systems/fulltext/S2405-4712(20)30503-2?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2405471220305032%3Fshowall%3Dtrue).

Using these germline mutations (which we referred to as "singletons" in the manuscript), we then mapped quantitative trait loci (QTL) for various phenotypes related to the mutation spectrum. Notably, we found a QTL for the rate of C>A germline mutations. 

This repository includes all of the code necessary to reproduce the main figures from our paper. This code is packaged into a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline, enabling users to generate all of these figures using a handful of simple commands, described below. 

## Table of Contents

1. [Dependencies](#dependencies)
2. [Directory structure](#directory-structure)
3. [Usage for generating manuscript figures](#usage-for-generating-manuscript-figures)
4. [Running tests](#running-tests)
5. [Dash application](#dash-application)

## Dependencies

1. `mamba` (by way of `conda`)

First, make sure that [`conda`](https://docs.conda.io/en/latest/) is installed and in your system `$PATH`. I used `conda v4.9.2`, but more recent versions will likely work, as well. 

Then, install [`mamba`](https://github.com/mamba-org/mamba) using the following command:

```
conda install mamba -n base -c conda-forge
```

> Why `mamba`?
> These days, `conda` can be *extremely* slow, and will often hang when building a new Python environment. I recommend using `mamba` as a 1-to-1 replacement for `conda` when building the Python env for this repository.

2. `perl`

Download [`perl`](https://www.perl.org) and make sure that the `perl` executable is in your system `$PATH`. I used `perl v5.32.1`, but more recent versions may work.

All other `python` and `R` dependencies will be handled by `mamba` before executing a pipeline. 

## Directory structure
```
.
|__rules                                # individual `snakemake` "rule files" that are imported by the main pipelines
|__py_scripts                           # all of the `python` scripts called by `snakemake` rules
|__R_scripts                            # all of the `R` scripts called by `snakemake` rules
|__Rqtl_data                            # data used by R/qtl2 for QTL mapping
|__data                                 # raw singleton data files, plus third-party datasets
|__tests                                # `pytest` tests for various utility functions
|__figure_generation.yaml               # YAML file containing all of the dependencies required to generate figures
|__generate_figures.smk                 # main `snakemake` pipeline that generates main and supplementary figures
```
## Usage for generating manuscript figures

All of the BXD singleton data (in addition to relevant third-party data) are included in the `data/` directory. The command-line instructions below will reproduce all of the main figures in the paper using these data.

```
# enter a directory of your choice
cd $WORKDIR

# clone the BXD analysis GitHub repository
git clone https://github.com/tomsasani/bxd_mutator_manuscript.git

# enter the directory
cd bxd_mutator_manuscript

# check out this simplified branch of the repo
git checkout simplified_code_sample

# create a new conda environment with all of the dependencies
# required for running the pipeline
mamba env create --name figure_generation --file figure_generation.yaml

# activate the new environment
conda activate figure_generation

# run the pipeline
snakemake \
        -j 1 \ # number of simultaneous jobs to run
        -s generate_figures.smk # name of the snakemake pipeline
```

This will produce plots from every main figure in the manuscript, which will be added to a new `plots/` directory.

#### What steps are involved?

1) Annotate singletons with various metadata.
2) Construct "tidy" dataframes containing summary information about mutation rates and spectra in each BXD.
3) Perform QTL scans for phenotypes related to the mutation rate.
4) Make plots.

In a number of analyses, we compare mutation spectra between the BXD singletons and previously published datasets. **I've included the files below in the `data/` directory, but here are instructions for downloading if necessary:** 

* singleton data from [Dumont 2019](https://academic.oup.com/mbe/article/36/5/865/5315518)
    * download ZIP file from [Supplementary data](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/mbe/36/5/10.1093_molbev_msz026/2/msz026_supp.zip?Expires=1616302324&Signature=u8neUFiV~0aBABDNG-ZPeMwd~usDZPmIO5TVjLHqKVcjHXrUWBm7MnR1ZJpSkMmDmQhMGrcdK~G7hySKLp79xgpQnj-SCFD09Hj7e9uCi9oYvVT-guMav1JY6qEMzSCubzlChpHfItUKJt15lXbxmuT2FxTibIs2gSrXvHCexmwGLxQYCoIAZJHY1nOjOfSDDlIejE-aGrPFozB86PXTZz~uM9JuAnmfZ5wmARxwuEzOHMfYZWh7WnWzeXEaNwKqYzrYHYDhej5sq~LSOfsQTzSPI-nrtn~KOV7x9ckk0RzqJ0kIhmF0uBLMkF8grDXTTRTHEjYV1dBALvxO1ZMVmA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA) associated with the manuscript
    * uncompress file and move `SuppTables_concat.xlsx` into the `data` directory of this repository
* *de novo* germline mutation data from [Ohno et al. 2014](https://www.nature.com/articles/srep04689)
    * download [Supplementary Data 1](https://static-content.springer.com/esm/art%3A10.1038%2Fsrep04689/MediaObjects/41598_2014_BFsrep04689_MOESM2_ESM.xls)
    * move `41598_2014_BFsrep04689_MOESM2_ESM.xls` into the `data` directory of this repository
* mutation signatures from [COSMIC](https://cancer.sanger.ac.uk/cosmic)
    * download [SBS36](https://cancer.sanger.ac.uk/sigs-assets-20/SBS_vignettes/sigProfiler_SBS_signatures_SBS36.csv)
    * download [SBS18](https://cancer.sanger.ac.uk/sigs-assets-20/SBS_vignettes/sigProfiler_SBS_signatures_SBS18.csv)
    * move these two files (`sigProfiler_SBS_signatures_SBS36.csv` and `sigProfiler_SBS_signatures_SBS18.csv`) into the `data` directory of this repository

## Running tests

There are a handful of `pytest` tests for various utility functions that get used by the singleton calling and figure generation pipelines. 

To run these tests, you'll first need to install the repository as a package using `pip`. Make sure you're in the top level of this directory, then run:

```
# activate previously created env
conda activate figure_generation

pip install -e .
```

Then, you can simply run `pytest .`, and the tests will run.

## Dash application

We've also made a Dash app that enables interactive exploration of some results from the manuscript. [Check it out!](https://bxd-mutator-exploration.herokuapp.com)