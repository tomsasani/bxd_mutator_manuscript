import argparse
import pandas as pd
from collections import defaultdict
from figure_gen_utils import convert_bxd_name

p = argparse.ArgumentParser()
p.add_argument(
    "--geno",
    required=True,
    help=
    """original geno file containing genotypes for all BXDs at each of ~7300 markers""",
)
p.add_argument(
    "--strain_metadata",
    required=True,
    help="""Excel file containing strain metadata""",
)
p.add_argument(
    "--out",
    required=True,
    help="""name of output geno file""",
)
args = p.parse_args()

outfh = open(args.out, "w")

# read in strain metadata
strain_metadata = pd.read_excel(args.strain_metadata)

# read in genotype file
geno = pd.read_csv(args.geno)

# map the GeneNetwork names of strains to their names as
# they appear in the BAM files/VCF sample info
gn2orig = dict(zip(strain_metadata['GeneNetwork name'], strain_metadata['bam_name']))

# map GeneNetwork names to the converted format used in all
# annotated/tidy data files
gn2conv = defaultdict()
for gn in gn2orig:
    orig = gn2orig[gn]
    if type(orig) != str: continue
    conv = convert_bxd_name(orig)
    gn2conv[gn.split(' ')[0]] = conv

def column_converter(cname, gn2conv):
    """
    helper function to apply to each of the
    column names in the geno dataframe (i.e., the
    gene network names of each BXD
    """
    if cname == "marker": return cname
    else:
        if cname in gn2conv: return gn2conv[cname]
        else: return "NOT_IN_METADATA"

# rename columns so that they match the "converted" BAM names
geno.columns = geno.columns.to_series().apply(lambda n: column_converter(n, gn2conv))

# output the formatted geno file
cols2keep = [c for c in list(geno) if c != "NOT_IN_METADATA"]
geno = geno[cols2keep]

geno.to_csv(args.out, index=False)
