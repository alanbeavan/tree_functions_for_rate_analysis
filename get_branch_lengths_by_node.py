#!usr/bin/env python3.6
"""
Write a table of the average lineage branch length either side of a node.

Get the node info from trees_to_analyse
Use functions from tree module to get the species-like side and the
deletion side.
Write each side to a table with 3 columns - tree, species like average branch
length, deletion average branch length.

Annoyingly, because I've written is_species_like using newick,
I need to write a convertion function
"""

import re
import sys
import statistics
import ete3
import newick
import tree_module as tm
import my_module as mod

def get_args():
    """Get user arguments."""
    if len(sys.argv) == 5:
        return sys.argv[1:]
    else:
        print("\nUSAGE: python get_branch_lengths_by_node.py treedir nodefile species_tree outfile\n")
        exit()

def extract_node(node_info, treedir):
    """Find the node to analyse in ete3."""
    treefile, nodelab = node_info.split()
    tree = ete3.Tree(treedir + "/" + treefile, format = 1)
    node = tree.search_nodes(name = nodelab)    
    if not node:   #Here there is only one outgroup taxon.
        for child in tree.children:
            if not re.search("_", child.name):
                node = child
    else:
        node = node[0]

    return node

def ete3_2_newick(tree):
    """Convert an ete3 tree to newick format."""
    string = tree.write()
    tree = newick.loads(string)[0]
    return tree

def get_ave_branch_length(tree):
    """Find the average root-tip distance for all leaves."""
    leaves = tree.get_leaves()
    distances = []
    for leaf in leaves:
        distance = 0
        node = leaf
        while node:
            distance += node.dist
            if node == tree:
                break
            node = node.up
        distances.append(distance)

    return statistics.mean(distances)


def write_table(outfile, species_like, with_deletion):
    """Write the results to a table."""
    out = open(outfile, "w")
    out.write("gene_tree\tspecies_like_length\tlength_with_deletion\n")
    for key in species_like.keys():
        out.write(key + "\t" + str(species_like[key]) + "\t" + str(with_deletion[key]) + "\n")
    out.close()

def main():
    """
    Use funcitons to:
    
    get user arguments.
    extract the node of interest.
    assess which node follows the species tree.
    get the average branch length.
    write to a table.
    """
    treedir, nodefile, species_tree, outfile = get_args()
    species_tree = tm.read_tree(species_tree)
    sp_nodes = tm.get_all_nodes(species_tree, tm.calculate_nnodes(species_tree))
    species_like = {}
    with_deletion = {}
    for line in mod.get_file_data(nodefile):
        try:
            treefile = line.split()[0]
            node = extract_node(line, treedir)
            comparison_node = ete3_2_newick(node.children[1])
            flag = 0
            for child in node.children:
                newick_tree = ete3_2_newick(child)
                if flag ==1:
                    with_deletion[treefile] = get_ave_branch_length(child)
                elif tm.node_with_deletion(comparison_node, newick_tree):
                    with_deletion[treefile] = get_ave_branch_length(child)
                else:
                    species_like[treefile] = get_ave_branch_length(child)
                    flag = 1
        except:
            continue


    write_table(outfile, species_like, with_deletion)

if __name__ == "__main__":
    main()
