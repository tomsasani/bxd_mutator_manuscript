import pandas as pd
import argparse

p = argparse.ArgumentParser()
p.add_argument("--tidy_mutation_spectra")
p.add_argument(
    "-sig_profiler_activities",
    default=
    "data/sigprofiler_outdir/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities_refit.txt",
)
args = p.parse_args()

spectra = pd.read_csv(args.tidy_mutation_spectra)

activ = pd.read_csv(args.sig_profiler_activities, sep='\t')
activ['has_sbs18'] = activ['SBS18'].apply(lambda s: int(s) > 0)
activ['has_sbs1'] = activ['SBS1'].apply(lambda s: int(s) > 0)
activ['has_sbs5'] = activ['SBS5'].apply(lambda s: int(s) > 0)
activ['has_sbs30'] = activ['SBS30'].apply(lambda s: int(s) > 0)


merged = activ.merge(spectra, left_on="Samples", right_on="bxd_strain_conv")
merged = merged.query("estimate_type == 'fraction' and base_mut == 'C>A'")
merged['predicted_to_have_sbs18'] = merged['haplotype_at_qtl'].apply(lambda h: h == "D")

print (merged.groupby('has_sbs18').size()) 

print (merged.groupby(['haplotype_at_qtl', 'has_sbs18']).size())

print (merged.query('predicted_to_have_sbs18 != has_sbs18')[['epoch', 'n_inbreeding_gens']])

for c in ['has_sbs30', 'has_sbs5', 'has_sbs1', 'has_sbs18']:
    print (merged.groupby(c).size())


groupby_cols = ["Samples", "SBS1", "SBS5", "SBS18", "SBS30"]
merged = merged[groupby_cols]
merged_tidy = merged.melt(id_vars=["Samples"], var_name="signature", value_name="sig_count")

bxd68_data = merged_tidy[merged_tidy['Samples'].str.contains("BXD68")]
print ("{} of bxd68 mutations are sbs18".format(bxd68_data.query('signature == "SBS18"')['sig_count'].values[0] / sum(bxd68_data['sig_count'])))
mutation_total = sum(merged_tidy.sig_count.values)

for sig in groupby_cols[1:]:
    print (sig, sum(merged_tidy.query(f"signature == '{sig}'").sig_count.values) / mutation_total)