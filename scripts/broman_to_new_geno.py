import csv
import argparse
import pandas as pd
from collections import defaultdict

p = argparse.ArgumentParser()
p.add_argument("--csv")
p.add_argument("--xls")
p.add_argument("--out")
def convert_bxd_name(name: str) -> str:
    """
    depending on the VCF file we're using we may
    have to convert BXD strain names for portability
    """
    bxd_line_name = '_'.join(name.split('_phased')[0].split('_')[1:])
    bxd_line_num = name.split('_phased')[0].split('_')[0].split('-')[-1]

    bxd_line_new = bxd_line_name + '_' + bxd_line_num

    return bxd_line_new

args = p.parse_args()

outfh = open(args.out, "w")

info = pd.read_excel(args.xls)

gn2orig = dict(zip(info['GeneNetwork name'], info['bam_name']))

gn2conv = defaultdict()
for gn in gn2orig:
    orig = gn2orig[gn]
    if type(orig) != str: continue
    conv = convert_bxd_name(orig)
    gn2conv[gn.split(' ')[0]] = conv

with open(args.csv, "r") as csvfh:
    f = csv.reader(csvfh)
    bad_indices = []
    for i,l in enumerate(f):
        if i < 3: continue
            #print (' '.join(l), file=outfh)
        elif i == 3:
            new_head = ['marker']
            for s_i,s in enumerate(l[1:]):
                if s not in gn2conv: bad_indices.append(s_i)
                else:
                    new_head.append(gn2conv[s])
            print (','.join(new_head), file=outfh)
        else:
            new_out = [l[0]]
            for v_i, v in enumerate(l[1:]):
                if v_i in bad_indices: continue
                new_out.append(v)
            print (','.join(new_out), file=outfh)
