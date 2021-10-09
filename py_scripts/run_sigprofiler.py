from SigProfilerExtractor import sigpro as sig
import argparse

p = argparse.ArgumentParser()
p.add_argument("--spectra")
p.add_argument("--outdir")
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
