import pandas as pd
import argparse
import doctest
from utils import *

doctest.testmod()


p = argparse.ArgumentParser()
p.add_argument("--strain_metadata", required=True, 
                    help="""metadata about strains in Excel format.""")
p.add_argument("--strain_genotypes", required=True,
                    help="""genotypes for each strain at each marker used in QTL scans""")
p.add_argument("--var_list", required=True, nargs="*",
                    help="""paths to per-chromosome BED files of variants""")
p.add_argument("--callable_bp", required=True,
                    help="""file containing a column with sample IDs and a \
                            column with the number of haploid callable base \
                            pairs in those samples.""")
p.add_argument("--out", required=True,
                    help="""name of the output extended BED file containing \
                            annotated variants""")
args = p.parse_args()

# read in strain metadata and construct relevant dictionaries that
# map sample names to metadata
summary = pd.read_excel(args.strain_metadata)

# convert verbose filial generations to integer numbers of generations in the file
summary['gen_at_seq'] = summary['Generation at sequencing'].apply(get_generation)

# map expanded strain names to generation numbers
strain2inbreed_gen = dict(zip(summary['bam_name'], summary['gen_at_seq']))

# map expanded strain names to epochs
strain2epoch = dict(zip(summary['bam_name'], summary['Epoch']))

# map expanded strain names to genenetwork names
strain2genenet = dict(zip(summary['bam_name'], summary['GeneNetwork name']))

# map expanded strain names to the number of generations
# each strain was intercrossed prior to inbreeding (0 for F2-derived strains)
strain2intercross_gens = dict(zip(summary['bam_name'],
                                  summary['n_advanced_intercross_gens']))

# map strain names to callable base pairs (referred to as "denominators")
denominators = pd.read_csv(args.callable_bp)
strain2denom = dict(zip(denominators['sample'], denominators['autosomal_callable_bp']))

# read in the genotypes of every strain at each QTL marker,
# and subset to the rsids that correspond to the length of the 
# 95% QTL credible interval. we use this data to query
# the strains' haplotype data to figure out which strains
# were B or D at the QTL
genos = pd.read_csv(args.strain_genotypes)

rsids = ["rs47460195",
         "rs52263933",
         "rs32445859",
         "rs13477933",
         "rs28272806",
         "rs28256540",
         "rs27509845"]

genos_at_markers = genos[genos['marker'].isin(rsids)]

# read in the file of variants, containing one BED-format line per variant
variants = combine_chr_df(args.var_list)

# add a column describing the 1-mer mutation type corresponding to all
# 3-mer "kmers" in the dataframe
variants['base_mut'] = variants['kmer'].apply(to_base_mut, cpg=True)

# remove any potential variants identified in founder genomes. we're
# only interested in variantss observed in BXD RILs. not sure if I've
# ever seen a singleton in a founder, but better safe than sorry
founder_samps = ['4512-JFI-0333_C57BL_6J_two_lanes_phased_possorted_bam',
                 '4512-JFI-0334_DBA_2J_three_lanes_phased_possorted_bam']

variants = variants[~variants['bxd_strain'].isin(founder_samps)]

# reformat BXD strain names from BAM notation
variants['bxd_strain_conv'] = variants['bxd_strain'].apply(lambda x: convert_bxd_name(x))

# add columns to the dataframe with relevant metadata
variants['epoch'] = variants['bxd_strain'].apply(lambda s: strain2epoch[s])
variants['n_inbreeding_gens'] = variants['bxd_strain'].apply(lambda s: strain2inbreed_gen[s])
# remove strains if they have "bad" values in the `n_inbreeding_gens` column,
# due to them being backcrossed during inbreeding 
variants = variants[~variants['n_inbreeding_gens'].isin([-1, "NA"])]
variants['n_intercross_gens'] = variants['bxd_strain'].apply(lambda s: strain2intercross_gens[s])
variants['n_callable_bp'] = variants['bxd_strain'].apply(lambda s: strain2denom[s] \
                                                                if s in strain2denom else "NA")
variants['haplotype_at_qtl'] = variants['bxd_strain_conv'].apply(lambda s: find_haplotype(genos_at_markers, 
                                                                                            s) if s in list(genos_at_markers) else "NA")
variants = variants[variants['haplotype_at_qtl'] != "NA"]
variants = variants.query('n_inbreeding_gens >= 20')

print ("Total of {} mutations in {} strains".format(variants.shape[0], len(pd.unique(variants['bxd_strain_conv']))))

variants.to_csv(args.out, index=False)
