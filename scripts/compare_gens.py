import pandas as pd 
import matplotlib.pyplot as plt

df = pd.read_csv("/Users/tomsasani/harrislab/bxd_mutator_ms/csv/tidy_mutation_rates.csv")

f, ax = plt.subplots()

ax.scatter(df['l_n'], df['n_inbreeding_gens'])

f.savefig('o.png')