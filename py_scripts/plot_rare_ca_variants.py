import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse

p = argparse.ArgumentParser()
p.add_argument("--variants")
p.add_argument("--out")
args = p.parse_args()


df = pd.read_csv(args.variants)

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

f.savefig(args.out, dpi=300)
