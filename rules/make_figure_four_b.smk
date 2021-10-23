PAL2NAL = "/Users/tomsasani/src/pal2nal.v14/pal2nal.pl"

FASTA = "data/ete3_data/multi_organism_mutyh_cds.fa"
PROTEIN = "data/ete3_data/multi_organism_mutyh_protein_alignment.fa"
TREE = "data/ete3_data/multi_organism_mutyh_tree.nwk"

rule format_tree:
	input: TREE
	output: "data/ete3_data/multi_organism_mutyh_tree.formatted.nwk"
	shell:
		"""
		sed -e 's/[0-9]*//g' {input} | sed -e 's/rodents//g' | \
									   sed -e 's/frogs_&_toads_and_birds//g' | \
									   sed -e 's/primates//g' | \
									   sed -e 's/rodents_and_primates//g' | \
									   sed -e 's/Multiple_organisms//g' | \
									   sed -e 's/_and_//g' | \
									   sed -e "s/\.//g" | sed -e "s/://g" | sed -e "s/e-//g" > {output}
		"""
		
rule format_protein_alignment:
	input: 
		protein = PROTEIN,
		script = "py_scripts/reformat_protein_alignment.py"
	output: 
		temp("data/ete3_data/multi_organism_mutyh_protein_alignment.reformatted.aln")
	shell:
		"""
		python {input.script} --msa {input.protein} -split_mmd > {output}
		"""

rule format_msa:
	input: 
		fasta = FASTA,
		script = "py_scripts/reformat_msa.py"
	output: 
		temp("data/ete3_data/multi_organism_mutyh_cds.reformatted.fa")
	shell:
		"""
		python {input.script} --fa {input.fasta} -split_mmd > {output}
		"""

rule codon_aware_aln:
	input: 
		protein = "data/ete3_data/multi_organism_mutyh_protein_alignment.reformatted.aln",
		msa = "data/ete3_data/multi_organism_mutyh_cds.reformatted.fa",
		pal2nal = PAL2NAL
	output: "data/ete3_data/multi_organism_mutyh_cds.codon_aligned.fa"
	shell:
		"""
		perl {input.pal2nal} -output fasta {input.protein} {input.msa} > {output} 
		"""

rule trim_aln:
	input: 
		msa = "data/ete3_data/multi_organism_mutyh_cds.codon_aligned.fa",
		script = "py_scripts/trim_sequences.py"
	output: "data/ete3_data/multi_organism_mutyh_cds.codon_aligned.trimmed.fa"
	shell:
		"""
		python {input.script} --msa {input.msa} -only_plot_mutyh > {output} 
		"""

rule visualize_phylo:
	input:
		msa = "data/ete3_data/multi_organism_mutyh_cds.codon_aligned.trimmed.fa",
		tree = "data/ete3_data/multi_organism_mutyh_tree.formatted.separate_mus.nwk"
	output:
		"plots/figure_4b.pdf"
	run:
		from ete3 import EvolTree, PhyloTree
		from ete3.treeview.layouts import evol_clean_layout

		tree = EvolTree(input.tree) 
		tree.link_to_alignment(alignment=input.msa, alg_format="fasta")
		
		tree.render(output[0], layout=evol_clean_layout)
		



