from SigProfilerExtractor import sigpro as sig
import argparse

p = argparse.ArgumentParser()
p.add_argument(
    "--spectra",
    required=True,
    help=
    """path to file containing mutation spectra counts in samples in SigProfiler format""",
)
p.add_argument(
    "--outdir",
    required=True,
    help="""name of output directory to store results""",
)
args = p.parse_args()

if __name__ == "__main__":

    sig.sigProfilerExtractor(
        'matrix',
        args.outdir,
        args.spectra,
        maximum_signatures=10,
        nmf_replicates=100,
        reference_genome="mm10",
        opportunity_genome="mm10",
    )
