import argparse
import matplotlib.pyplot as plt
import glob
import csv
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('--haplotype_dir')
args = p.parse_args()

b6, d2 = 0, 0

for f in glob.glob(args.haplotype_dir + '*haplotypes.csv'):
    with open(f, 'r') as fh:
        reader = csv.reader(fh)
        for l in reader:
            c, s, e, hap = l
            length = int(e) - int(s)
            if hap == "1": d2 += length
            else: b6 += length
print (b6, d2)
