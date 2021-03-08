import csv
import sys
import argparse
from pomegranate import *
import pandas as pd
from collections import OrderedDict
import numpy as np

def bake_hmm():
	"""
	generate a simple HMM with pre-set emission and
	transition probabilities.
	"""

	# emission probabilities for B6 and D2 regions
	# logic here is that given the noise in genotyping
	# B6 vs D2 sites (any given region will likely contain
	# a mix), we can model the D2 state has having a little
	# bit of B6 in it, and vice versa. choice of a 90/10 split
	# is totally arbitrary, though.
	d1 = DiscreteDistribution({0: 0.9, 2: 0.1})
	d2 = DiscreteDistribution({2: 0.9, 0: 0.1})
	s1 = State(d1, name='B6')
	s2 = State(d2, name='D2')
   
	# initialize model  
	model = HiddenMarkovModel()
	model.add_states(s1, s2)
   
	# add transition probabilities -- assume a 
	# "real" transition from a B6->D2 or D2->B6 haplotype
	# is exceedingly rare
	model.add_transition(model.start, s1, 0.5)
	model.add_transition(model.start, s2, 0.5)
	model.add_transition(s1, s1, 0.9999)
	model.add_transition(s2, s2, 0.9999)
	
	model.add_transition(s1, s2, 0.0001)
	model.add_transition(s2, s1, 0.0001)

	model.bake()
	
	return model

p = argparse.ArgumentParser()
p.add_argument("--inf_site_positions", required=True, nargs="*",
					help="""positions of each informative site""")
p.add_argument("--inf_site_states", required=True, nargs="*",
					help="""BXD genotypes at each informative site""")
p.add_argument("--sample", required=True,
					help="""name of sample we want to get haplotypes for""")
p.add_argument("--out", required=True,
					help="""name of output file""")
args = p.parse_args()

# initialize the HMM
hmm_model = bake_hmm()

# --
# output full haplotype blocks to a new file
# --
# assign HMM predictions to this strain's states
# we'll generate a block file for each strain
hap_fh = open(args.out, 'w')

chroms = range(1, 20)
chroms = list(map(str, chroms))

# sort inputs by chromosome
inf_pos = list(sorted(args.inf_site_positions, 
						key=lambda p: p.split('/')[-1].split('_')[0]))
inf_gts = list(sorted(args.inf_site_states, 
						key=lambda p: p.split('/')[-1].split('_')[0]))

for states_fh, sites_fh in zip(inf_gts, inf_pos):

	chrom = states_fh.split('/')[-1].split('_')[0]

	# first, read in the states file, which is just a
	# matrix of sites x samples. we can store as int8
	# since values must be in (0,1,2)
	states_matrix = pd.read_csv(states_fh, header=2, dtype=np.int8)

	# get states for sample of interest
	states = states_matrix[args.sample]

	# then, read in the positions files, which is just
	# a single column of positions that we can access
	# by corresponding index in the states matrix. these
	# positions can be in the hundreds of millions, so we'll 
	# use int32
	pos_matrix = pd.read_csv(sites_fh, header=None, dtype=np.int32)

	# pomegranate can handle missing values, so we'll
	# pretend that UNK and HET sites are missing
	states[states == 1] = np.nan

	# get our HMM predictions
	preds = np.array(hmm_model.predict(states))

	prev_pos, prev_state = 0, -1
	prev_xo_pos = 0
	# loop over every position across the chromosome
	for i,pred in enumerate(preds):
		if i == 0:
			# if we're at the first state, record the position 
			# and state. this will be the start of our first hap block.
			# we can grab the position by accessing the corresponding
			# position in `pos_matrix`
			prev_pos, prev_state = pos_matrix.at[i, 0], pred 
			prev_xo_pos = pos_matrix.at[i, 0]
		else:
			# if we reach a state that doesn't match our previous state,
			# output the full block
			if pred != prev_state:
				print (','.join(list(map(str, [chrom, prev_pos, 
									pos_matrix.at[i, 0], prev_state]))), file=hap_fh)
				# then, reset the last seen state and its position
				prev_pos, prev_state = pos_matrix.at[i, 0], pred 
				prev_xo_pos = pos_matrix.at[i, 0]
			# if we hit the last line of the dataframe,
			# output whatever the last block is.
			elif i == preds.shape[0] - 1:
				print (','.join(list(map(str, [chrom, prev_pos, 
									pos_matrix.at[i, 0], prev_state]))), file=hap_fh)
			# otherwise, do nothing
			else: 
				prev_xo_pos = pos_matrix.at[i, 0]
				continue
