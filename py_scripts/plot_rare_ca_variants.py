import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse

p = argparse.ArgumentParser()
p.add_argument("--variants", nargs="*")
p.add_argument("--out")
args = p.parse_args()

out_df = []

for fh in args.variants:
    df = pd.read_csv(fh)
    out_df.append(df)

df = pd.concat(out_df).groupby(["sample", "ancestry"]).sum().reset_index()

print (df)

df['d2_frac'] = df['d2_count'] / (df['d2_count'] + df['b6_count'])
df['ca_frac'] = df['ca_count'] / df['mut_count']

f, ax = plt.subplots()

sns.scatterplot(
    x="d2_frac",
    y="ca_frac",
    hue="ancestry",
    data=df,
    ax=ax,
)

ax.set_ylabel("Fraction of rare variants that are C>A")
ax.set_xlabel("Proportion of fixed differences between\nB6 and D2 that are D2 in the focal strain")

f.tight_layout()

f.savefig(args.out)
