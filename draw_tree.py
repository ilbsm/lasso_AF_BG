from Bio import Phylo

tree = list(Phylo.parse('align/_representatives.phylotree.ph', 'newick'))[0]
Phylo.draw(tree)

