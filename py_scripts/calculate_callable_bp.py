import csv
import argparse
from figure_gen_utils import make_interval_tree

p = argparse.ArgumentParser()
p.add_argument(
    '--coverage_files',
    required=True,
    nargs="*",
    help="""list of paths to sample threshold BED files""",
)
p.add_argument(
    '--exclude',
    required=True,
    help="""path to file containing regions we want to mask""",
)
p.add_argument(
    '--out',
    required=True,
    help="""path to output file""",
)
args = p.parse_args()

outfh = open(args.out, "w")
print(",".join(["sample", "autosomal_callable_bp"]), file=outfh)

exclude = make_interval_tree(args.exclude)

chroms = list(map(str, range(1, 20)))
chroms = ['chr' + c for c in chroms]

for fh in args.coverage_files:
    # read in file with coverage information from mosdepth
    sample = fh.split('/')[-1].split('.')[0]

    total_bp_callable = 0

    with open(fh, "r") as f:
        csvf = csv.reader(f, delimiter='\t')
        for l in csvf:
            # skip header line
            if l[0] == "#chrom": continue
            # make sure we only look at autosomes
            if l[0] not in chroms: continue
            chrom, start, end, _, _, _, tenx, _ = l

            # figure out how much of this interval overlaps the masked regions
            exclude_overlaps = exclude[chrom].search(int(start), int(end))
            exclude_overlap_sum = sum(
                [e.end - e.start for e in exclude_overlaps])

            # get the fraction of bases masked in this interval
            # (i.e., in seg dups or simple repeats)
            interval_size = int(end) - int(start)
            masked_frac = exclude_overlap_sum / interval_size

            # we'll assume that if X% of bases in the interval were masked by
            # the exclude file, then X% of the bases covered by at least 10 reads
            # were also masked by the exclude file
            adj_tenx_num = int(tenx) * (1 - masked_frac)
            total_bp_callable += adj_tenx_num

            #if sample == "4512-JFI-0354_BXD199_phased_possorted_bam":
            #    print (chrom, start, end, exclude_overlap_sum, interval_size, tenx)

    print(",".join([sample, str(total_bp_callable)]), file=outfh)
