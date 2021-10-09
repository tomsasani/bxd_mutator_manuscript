from skbio.alignment import TabularMSA
from skbio import DNA, Protein
import skbio.io
from io import StringIO
import argparse
import numpy as np
from collections import Counter

p = argparse.ArgumentParser()
p.add_argument("--msa")
p.add_argument("-split_mmd", action="store_true")
args = p.parse_args()

align = TabularMSA.read(args.msa, constructor=Protein, format="fasta" )
   
def change_seq(seq, bxd68=False):
    seq_changes = {5: "R",
                   24: "C",
                   69: "R",
                   153: "Q",
                   312: "P",
                   313: "P" }

    codons_to_change = [5, 24, 69, 312, 313]
    if bxd68: codons_to_change = [5, 24, 69, 153, 312, 313]

    new_seq = []
    non_gap_codons = 1
    for codon_i, codon in enumerate(list(seq)):
        if codon == '-':
            new_seq.append(codon)
        elif non_gap_codons in codons_to_change:
            new_seq.append(seq_changes[non_gap_codons])
            non_gap_codons += 1
        else:
            new_seq.append(codon)
            non_gap_codons += 1


    return ''.join(new_seq)

with StringIO() as fh:
    info =  align.write(fh).getvalue().split('\n')

    prev_name = None
    for i,e in enumerate(info):
        if i % 2 == 0:
            if len(e) == 0: continue
            organism = e.split('organism')[1].split(']')[0].split('=')[1].split(' ')
            prev_name = '_'.join(organism)

        if i % 2 != 0:
            if prev_name == "Mus_musculus" and args.split_mmd:

                protein_seq = e.replace('-', '')

                print ('>{}\n{}'.format("C57BL/6J", e))
                print ('>{}\n{}'.format("DBA/2J", change_seq(e)))
                print ('>{}\n{}'.format("BXD68", change_seq(e, bxd68=True)))

            else:
                print ('>{}\n{}'.format(prev_name, e))


