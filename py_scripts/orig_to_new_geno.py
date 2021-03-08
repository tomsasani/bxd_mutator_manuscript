import csv
import argparse
import pandas as pd
from collections import defaultdict
from figure_gen_utils import convert_bxd_name

p = argparse.ArgumentParser()
p.add_argument("--geno", required=True, 
                    help="""original geno file containing genotypes
                            for all BXDs at each of ~7300 markers""")
p.add_argument("--strain_metadata", required=True,
                    help="""Excel file containing strain metadata""")
p.add_argument("--out", required=True,
                    help="""name of output geno file""")
args = p.parse_args()

outfh = open(args.out, "w")

# read in strain metadata
strain_metadata = pd.read_excel(args.strain_metadata)

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

# read in the original geno file
with open(args.geno, "r") as csvfh:
    f = csv.reader(csvfh)
    bad_indices = []
    for i,l in enumerate(f):
        # skip the header
        if l[0].startswith("#"): continue
        # generate a new header with updated
        # sample names
        elif l[0] == "marker":
            new_header = ['marker']
            # loop over samples in the order they appear in the 
            # original geno header
            for s_i,s in enumerate(l[1:]):
                # keep track of samples that are in the original geno
                # file but that I leave out of my own analyses
                if s not in gn2conv: bad_indices.append(s_i)
                else: new_header.append(gn2conv[s])
            print (','.join(new_header), file=outfh)
        else:
            new_out = [l[0]]
            # loop over haplotypes in the order they appear
            for v_i, v in enumerate(l[1:]):
                if v_i in bad_indices: continue
                new_out.append(v)
            print (','.join(new_out), file=outfh)
