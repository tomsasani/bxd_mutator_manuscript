import pandas as pd
from figure_gen_utils import to_base_mut, revcomp
from collections import defaultdict
import argparse

def convert_kmer_to_sig_profiler_format(k):
    """
    convert the mutation types from in the BXD data
    to a format acceptable to SigProfilerExtractor. namely,
    convert A>T mutations to their reverse complements and
    format mutations as a A[B>D]C format, where A is the five-prime
    nucleotide, B>D is the mutation, and C is the three-prime nucleotide.
    """
    base_mut = to_base_mut(k)
    fprime, tprime = k[0], k[2]
    if base_mut[0] == "A":
        base_mut = "{}>{}".format(revcomp(base_mut.split('>')[0]),
                                  revcomp(base_mut.split('>')[1]))
        fprime, tprime = revcomp(k[0]), revcomp(k[2])
    new_mutation_type = "{}[{}]{}".format(fprime, base_mut, tprime)
    return new_mutation_type

p = argparse.ArgumentParser()
p.add_argument("--singleton_vars", required=True,
            help="""annotated variants in extended BED format""")
p.add_argument("--out", required=True)
args = p.parse_args()

singleton = pd.read_csv(args.singleton_vars)

group_cols = [
    "bxd_strain_conv",
    'kmer',
    "haplotype_at_qtl",
]

# convert to wide-form dataframe
singleton_tidy = singleton.groupby(group_cols).count().add_suffix('_count').reset_index()
singleton_tidy = singleton_tidy[group_cols + ["chrom_count"]]
singleton_tidy['base_mut'] = singleton_tidy['kmer'].apply(to_base_mut, cpg=False)

strain2mut_count = defaultdict(lambda: defaultdict(int))

# make a new column in the dataframe with the reformatted
# SigProfilerExtractor mutation type
singleton_tidy['Mutation Types'] = singleton_tidy['kmer'].apply(lambda k: convert_kmer_to_sig_profiler_format(k))

# store the singleton counts of every every mutation
# type in every strain
for strain, sub_df in singleton_tidy.groupby([
        "bxd_strain_conv",
        "Mutation Types",
]):
    s, m = strain
    strain2mut_count[m][s] = sum(sub_df.chrom_count)

out_fh = open(args.out, "w")

strains = singleton_tidy.bxd_strain_conv.unique()

# generate the header
header = ["Mutation Types"]
header.extend(list(map(str, strains)))
print ('\t'.join(header), file=out_fh)

for m, m_df in singleton_tidy.groupby('Mutation Types'):
    counts = []
    for s in strains:
        if s not in strain2mut_count[m].keys(): counts.append(0)
        else: counts.append(strain2mut_count[m][s])

    base_mut = m_df.base_mut.unique()[0]
    out_vals = [m]
    out_vals.extend(list(map(str, counts)))

    print ('\t'.join(out_vals), file=out_fh)
