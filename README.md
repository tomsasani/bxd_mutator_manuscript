# A wild-derived antimutator drives germline mutation spectrum differences in a genetically diverse murine family

*Thomas A. Sasani, David G. Ashbrook, Annabel C. Beichman, Lu Lu, Abraham A. Palmer, Robert W. Williams, Jonathan K. Pritchard, Kelley Harris*

:mouse: --> :dna: --> :bar_chart:

The code in this repository uses [Snakemake](https://snakemake.readthedocs.io/en/stable/) to reproduce the entire [manuscript](https://www.biorxiv.org/content/10.1101/2021.03.12.435196v1) from "top to bottom." This includes everything from downloading a reference genome to generating supplementary figures. However, it's also possible to simply generate the figures in the manuscript (step #2 below), since step #1 requires very large input files and quite a bit of time to execute.

**The basic outline is as follows:**

1) Identify high-quality singleton mutations from the BXD VCF using `identify_singletons.smk`.

2) Annotate singleton mutations and generate figures using `generate_figures.smk`.

**If you just want to generate all of the manuscript figures, skip to the section entitled [Usage for generating manuscript figures using precomputed data](#usage-for-generating-manuscript-figures-using-precomputed-data).**

**IMPORTANT NOTE:**

>I highly recommend that the `identify_singletons.smk` pipeline is run on a high-performance computing system or cluster. The `identify_singletons.smk` pipeline assumes you've downloaded the BXD VCF (~50 Gb), and the hundreds of individual steps in the pipeline include downloading the mm10 reference genome (~1 Gb) and many Gbs of phastCons scores.

>Depending on the particular steps of the pipeline you want to run, some of these steps/downloads can be avoided by editing the `identify_singletons.smk` pipeline directly. For example, if you don't want to get singleton data from the wild mouse genomes, simply comment out those output files in `rule: all`. Or if you already have a copy of the mm10 reference, just add it to the `data/ref` directory.

## Table of Contents

1. [Dependencies](#dependencies)
2. [Directory structure](#directory-structure)
3. [Usage for generating manuscript figures using precomputed data](#usage-for-generating-manuscript-figures-using-precomputed-data)
4. [Usage for generating all raw data](#usage-for-generating-all-raw-data)
5. [How to generate raw data on a HPC](#how-to-generate-raw-data-on-a-hpc)

## Dependencies
Make sure that these are installed and in your system `$PATH`! Versions in parentheses are the ones I used at the time of manuscript posting. I haven't experimented with other versions, so YMMV.

### Common requirements for all pipelines
* [conda (v4.9.2)](https://docs.conda.io/en/latest/)

### Requirements if you want to generate singleton data from scratch
* [bedtools (v2.29.2)](https://bedtools.readthedocs.io/en/latest/)
* [bedops (v2.4.38)](https://bedops.readthedocs.io/en/latest/)
* [tabix (v1.10.2-125-g4162046)](http://www.htslib.org/doc/tabix.html)
* [mutyper (v0.5.0)](https://harrispopgen.github.io/mutyper/install.html)
* [bcftools (v1.12)](https://samtools.github.io/bcftools/bcftools.html)
* [bgzip (v1.12)](http://www.htslib.org/doc/bgzip.html)

### Requirements if you want to reproduce manuscript figures + analysis
* [perl (v5.32.1)](https://www.perl.org)
* [pal2nal (v14)](http://www.bork.embl.de/pal2nal/#Download)
* [SigProfilerExtractor (v1.1.3)](https://github.com/AlexandrovLab/SigProfilerExtractor)
* [codeml](http://abacus.gene.ucl.ac.uk/software/paml.html)
    * `codeml` is used by `ete3` for the PAML selection scans. Rather than installing `codeml` directly, I recommend using `ete3 upgrade-external-tools` to install, followed by copying the `codeml` binary from `/Users/your_username/.etetoolkit/ext_apps-latest/bin/codeml` (the default installation directory for `ete3 upgrade-external-tools`) into `/Users/your_username/opt/anaconda3/envs/figure_generation/bin/` (i.e., the `bin` directory of the `conda` installation directory created after initializing an environment below).

All other `python` and `R` dependencies will be handled by `conda` before executing a pipeline. 

## Directory structure
```
.
|__rules                                # individual `snakemake` "rule files" that are imported by the main pipelines
|__py_scripts                           # all of the `python` scripts called by `snakemake` rules
|__R_scripts                            # all of the `R` scripts called by `snakemake` rules
|__Rqtl_data                            # data used by R/qtl2 for QTL mapping
|__data                                 # raw data files output by the `identify_singletons.smk` pipeline, plus third-party datasets
|__figure_generation.yaml               # YAML file containing all of the dependencies required to generate figures
|__singleton_calling.yaml               # YAML file containing all of the dependencies required to call singletons
|__generate_figures.smk                 # main `snakemake` pipeline that generates main and supplementary figures
|__identify_singletons.smk              # main `snakemake` pipeline that identifies singletons using the BXD VCF
|__singleton_calling_hpc_config.json    # JSON config file template for running pipelines on cluster like SGE 
```
## Usage for generating manuscript figures using precomputed data

If you'd rather not spend the time (and compute) needed to generate all of the singleton data from scratch, I've included the output of `identify_singletons.smk` in the `data/` directory. The command-line instructions below will reproduce all but one or two figures in the manuscript (like Supp. Fig. 1), which were created by hand in Illustrator.

#### What steps are involved?

1) Annotate singletons and fixed variants with various metadata.
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

```
# enter a directory of your choice
cd $WORKDIR

# clone the BXD analysis GitHub repository
git clone https://github.com/harrispopgen/bxd_mutator_manuscript.git

# enter the directory
cd bxd_mutator_manuscript

# create a new conda environment with all of the dependencies
# required for running the pipeline
conda env create --name figure_generation --file figure_generation.yaml

# activate the new environment
conda activate figure_generation

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
3) Identify singletons in each BXD.
4) Identify "fixed" variants in D2 that we'll use for conservation score comparisons.

To generate all of the **raw data in the manuscript**, it's highly recommended that you run the following on a machine with the ability to run many simultaneous processes, and **with at least 20 Gb of free space in the working directory.**

**The pipeline below also assumes that you've downloaded the BXD VCF (about 50 Gb). The path to this VCF must be edited in the Snakemake file.**

```
# enter a directory of your choice (make sure
# you have at least 20 Gb of free space here)
cd $WORKDIR

# clone the BXD analysis GitHub repository
git clone https://github.com/harrispopgen/bxd_mutator_manuscript.git

# enter the directory
cd bxd_mutator_manuscript

# create a new conda environment with all of the dependencies
# required for running the pipeline
conda env create --name singleton_calling --file singleton_calling.yaml

# activate the new environment
conda activate singleton_calling

# install mutyper
pip install mutyper==0.5.0

# ---
# before running the pipeline, use your text editor of choice
# and make sure that the WORKDIR variable at the top of 
# `identify_singletons.smk` is the absolute path to the directory
# that you're currently in, and that all of
# the other paths at the top of the pipeline point to the right
# binaries and directories on your machine.
# ---

# run the pipeline
snakemake \
        -j 20 \ # number of simultaneous jobs to run (change if you want)
        -s identify_singletons.smk # name of the snakemake pipeline
```

## How to generate raw data on a HPC

The Department of Genome Sciences at UW uses the Sun Grid Engine (SGE) for job management on the cluster. As an example, I've included a config file that can be used in conjunction with the following command to run the `identify_singletons.smk` pipeline on a system with SGE. Just edit the various parameters in the file accordingly. If you use a different job management system (e.g., SLURM), there is documentation on how to submit jobs using `snakemake` [here](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

```
# run the pipeline using SGE
snakemake \
        -j 20 \
        -s identify_singletons.smk \
        --cluster-config singleton_calling_hpc_config.json \
        --cluster "qsub -l centos={cluster.os} \
                        -l mfree={cluster.memory} \
                        -l h_rt={cluster.time} \
                        -l s_rt={cluster.time} \
                        -o {cluster.erroutdir} \
                        -e {cluster.erroutdir}" \
        --latency-wait 30 \
        --rerun-incomplete
```

After generating raw data, you can run the figure generation pipeline as described above. 