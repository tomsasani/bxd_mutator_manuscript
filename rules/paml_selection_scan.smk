PAL2NAL = "/Users/tomsasani/src/pal2nal.v14/pal2nal.pl"

FASTA = "data/ete3_data/multiple_murine/rodent_refseq_cds.fa"
PROTEIN = "data/ete3_data/multiple_murine/rodent_protein_alignment.fa"
TREE = "data/ete3_data/multiple_murine/rodent_tree.formatted.nwk"

rule all:
	input: "data/ete3_data/multiple_murine/selection_scan_results.csv"

rule codon_aware_aln:
	input: 
		protein = "data/ete3_data/multiple_murine/rodent_protein_alignment.fa",
		msa = "data/ete3_data/multiple_murine/rodent_refseq_cds.fa",
		pal2nal = PAL2NAL
	output: temp("data/ete3_data/multiple_murine/rodent_refseq_cds.codon_aligned.fa")
	shell:
		"""
		perl {input.pal2nal} -output fasta {input.protein} {input.msa} > {output} 
		"""

rule trim_aln:
	input: 
		msa = "data/ete3_data/multiple_murine/rodent_refseq_cds.codon_aligned.fa",
		script = "py_scripts/trim_sequences.py"
	output: temp("data/ete3_data/multiple_murine/rodent_refseq_cds.codon_aligned.trimmed.fa")
	shell:
		"""
		python {input.script} --msa {input.msa} > {output} 
		"""

rule run_selection_scans:
	input:
		msa = "data/ete3_data/multiple_murine/rodent_refseq_cds.codon_aligned.trimmed.fa",
		tree = "data/ete3_data/multiple_murine/rodent_tree.formatted.nwk"
	output:
		"data/ete3_data/multiple_murine/selection_scan_results.csv"
	run:
		from ete3 import EvolTree, PhyloTree
		from ete3.treeview.layouts import evol_clean_layout

		tree = EvolTree(input.tree) 
		tree.link_to_alignment(alignment=input.msa, alg_format="fasta")
		
		def compare_models(tree, m1="M8a", m2="M8"):
			tree.run_model (m1)
			tree.run_model (m2)
			
			pval = tree.get_most_likely (m2,m1)
		
			model2 = tree.get_evol_model(m2)

			if pval < 0.05:
				print (f'{m2} model wins.')
				for s in range(len(model2.sites['BEB']['aa'])):
					if model2.sites['BEB']['p2'][s] > 0.95:
						print ('positively selected site %s at position: %s, with probability: %s') % (model2.sites['BEB']['aa'][s], s+1, model2.sites['BEB']['p2'][s])
			else:
				print (f'{m1} model is not rejected')

			return pval

		m1_v_m2 = compare_models(tree, m1="M1", m2="M2")
		m7_v_m8 = compare_models(tree, m1="M7", m2="M8")

		with open(output[0], "w") as fh:
			print (f"M1 vs. M2 p-value: {m1_v_m2}", file=fh)
			print (f"M7 vs. M8 p-value: {m7_v_m8}", file=fh)
		



