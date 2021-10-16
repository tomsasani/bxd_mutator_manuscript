from cyvcf2 import VCF
from collections import defaultdict
import argparse

p = argparse.ArgumentParser()
p.add_argument("--out")
p.add_argument("--species")
args = p.parse_args()

vcfh = "http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz"

vcf = VCF(vcfh, gts012=True)

# map sample indices to sample IDs and vice versa
smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))
# map "ancestries" (i.e., mouse subspecies) to sample indexes
anc2idx = defaultdict(list)
anc = [s.split('_')[0] for s in vcf.samples]
smp2anc = dict(zip(vcf.samples, anc))
for smp in smp2idx:
    anc = smp2anc[smp]

outfh = open(args.out, "w")
for s in smp2anc:
    if smp2anc[s] != args.species: continue
    print (s, file=outfh)