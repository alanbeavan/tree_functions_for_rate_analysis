#!usr/bin/env python3.6
"""
Extract the subtrees of interest and their sequences for analysis.

From the list of nodes with the motif looked for:
Read in the tree.
Take the subtree (node) that is the parent of the duplication node we want to
    analyse.
Write this to a file with a unique name using indexing.
Also read the sequences with names identical to those in the tree and write
    them to a sequence file with the same index.
"""

import sys
import newick
import tree_module as tm
import my_module as mod

def get_args():
    """Read user arguments, else remind user of input."""
    if len(sys.argv) == 6:
        return sys.argv[1:]
    print("USAGE: python extract_nodes_and_sequences.py nodeslist " + \
          "input_treedir input_seqdir output_treedir output_seqdir\n")
    exit()

def extract_treenode(node_info, treedir):
    """
    Take the node from a tree above that specified in the node_info.

    ARGUMENTS:
    node_info - the name of the orthgroup and the node label
    treedir - the location of the tree files

    RETURNS:
    node_with_outgroup - the node of interest with the smallest possible
        outgroup attached.

    REQUIRES:
    newick
    tree_module
    """
    og, node_lab = node_info.split()
    tree = tm.read_tree(treedir + "/" + og + "_tree.txt")
    print(treedir + "/" + og + "_tree.txt")
    nodes = []
    nodes = tm.get_all_nodes(tree, tm.calculate_nnodes(tree), [], [])
    print(nodes[0:3])
    for node in nodes:
        for des in node.descendants:
            if des.name == node_lab:
                return node
    print(og + " " + node_lab + " didn't work for some reason")
    return 0

def write_tree(tree, treedir, index):
    """
    Write the tree to a file.

    ARGUMENTS:
    tree - the tree to be written
    treedir - the location of the output file
    index - a unique identifier to tag the filename with.

    REQUIRES:
    newick
    """
    newick.write(tree, treedir + "/tree_" + str(index) + ".txt")


def write_seqs(tree, in_seqdir, out_seqdir, index, node_info):
    """
    Write sequences of proteins from tree to a file.

    ARGUMENTS:
    tree - the tree with tips that need sequences written
    seqdir - the location of the output file
    index - a unique identifier to tag the filename with.
    """
    og = node_info.split()[0]
    print(tree)
    tips = tree.get_leaf_names()
    print(tips)
    in_seqs = mod.read_fasta(in_seqdir + "/" + og + ".fa")
    newseqs = {}
    for tip in tips:
        print(og + " " + tip)
        newseqs[str(tip)] = in_seqs[str(tip)]
    out = open(out_seqdir + "/seqs_" + str(index) + ".fa", "w")
    for key, value in newseqs.items():
        out.write(">" + key + "\n" + value + "\n")
    out.close()


def main():
    """Using tree reading and writing funciutons, perform the above."""
    sys.setrecursionlimit(20000)
    nodesfile, in_treedir, in_seqdir, out_treedir, out_seqdir = get_args()
    nodeslines = mod.get_file_data(nodesfile)
    i = 0
    for line in nodeslines:
        print(line)
        for local in locals():
            print(local)
        if 'subtree' in locals():
            print(subtree)
        subtree = extract_treenode(line, in_treedir)
        print(subtree)
        if subtree:
            write_tree(subtree, out_treedir, i)
            write_seqs(subtree, in_seqdir, out_seqdir, i, line)
            print()
        i += 1


if __name__ == "__main__":
    main()
