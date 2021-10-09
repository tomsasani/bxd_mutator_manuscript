import sys
from Bio.Seq import Seq
from Bio import SeqIO
import argparse

def change_mm_seq(seq, bxd68=False):
    """
    generate a new transcript sequence for DBA/2J
    and C57BL/6J using the Mus musculus as the baselin
    """
    
    seq_changes = {5: (1, 'A', 'G'),
                   24: (0, 'C', 'T'),
                   69: (2, 'C', 'G'),
                   153: (1, 'G', 'A'),
                   312: (0, 'A', 'C'),
                   313: (0, 'T', 'C') }

    new_seq = ''

    codon_nums = [5, 24, 69, 312, 313]
    if bxd68: codon_nums = [5, 24, 69, 153, 312, 313]

    prev_nuc = 0
    for codon_num in codon_nums:
        # figure out the position of the nucleotide
        # that gets changed in the transcript to produce
        # the amino acid change of interest
        nuc_loc, nuc_ref, nuc_alt = seq_changes[codon_num]

        nuc_num = (codon_num * 3) - 3 + nuc_loc

        new_seq += seq[prev_nuc:nuc_num] + nuc_alt

        prev_nuc = nuc_num + 1

    new_seq += seq[prev_nuc:]

    return new_seq

p = argparse.ArgumentParser()
p.add_argument("--fa")
p.add_argument("-split_mmd", action="store_true")
args = p.parse_args()

record_dict = SeqIO.to_dict(SeqIO.parse(args.fa, "fasta"))

for r in record_dict:
    seq_record = record_dict[r]

    full_name = []

    desc = seq_record.description.split(' ')
    for i,d in enumerate(desc):
        if i == 0: continue
        if "PREDICTED" in d: continue
        if "mutY" in d or "lymphocyte" in d or "serotransferrin" in d or "transferrin" in d: break
        full_name.append(d)

    if '_'.join(full_name) == "Mus_musculus" and args.split_mmd:

        print (">{}\n{}".format("C57BL/6J",
                                seq_record.seq))
        dba_seq = change_mm_seq(seq_record.seq)

        print (">{}\n{}".format("DBA/2J",
                                dba_seq))
        bxd68_seq = change_mm_seq(seq_record.seq, bxd68=True)
        print (">{}\n{}".format("BXD68",
                                bxd68_seq))
    else: 
        print (">{}\n{}".format('_'.join(full_name),
                                seq_record.seq))

