from Bio import Phylo

tree = list(Phylo.parse('align/_representatives.tree.dnd', 'newick'))[0]
Phylo.draw(tree)
