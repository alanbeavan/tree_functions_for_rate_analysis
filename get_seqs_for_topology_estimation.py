#!usr/bin/env python3.6
"""
Write the sequences needed for my analysis.

Read trees and find the node of interest.
Take the seqeuences of those tips plus one sequence from the nearest outgroup.
Write to a file
"""

import sys
import newick
import tree_module as tm
import my_module as mod

def get_args():
    """Read user arguments, else remind user of input."""
    if len(sys.argv) == 5:
        return sys.argv[1:]
    print("USAGE: python get_seqs_for_topology_estimation nodeslist " + \
          "input_treedir input_seqdir output_seqdir\n")
    exit()

def write_seqs(tree, in_seqdir, out_seqdir, index, node_info):
    """
    Write sequences of proteins from tree to a file.
    These are the sequences involved in the duplication node and one from
    the nearest outgroup.

    ARGUMENTS:
    tree - the tree with tips that need sequences written
    seqdir - the location of the output file
    index - a unique identifier to tag the filename with.
    """
    og, node_lab = node_info.split()
    print(og, node_lab)
    print(tree)
    tips = []
    for node in tree.descendants:
        if node.name == node_lab:
            tips.extend(node.get_leaf_names())
        else:
            tips.append(node.get_leaf_names()[0])
            outgroup = node.get_leaf_names()[0]
    in_seqs = mod.read_fasta(in_seqdir + "/" + og + ".fa")
    newseqs = {}
    for tip in tips:
        print(og + " " + tip)
        newseqs[str(tip)] = in_seqs[str(tip)]
    out = open(out_seqdir + "/seqs_" + str(index) + ".fa", "w")
    for key, value in newseqs.items():
        out.write(">" + key + "\n" + value + "\n")
    out.close()
    out = open("outgroups", "a")
    out.write(outgroup + "\n")
    out.close()


def main():
    """Using tree reading and writing funciutons, perform the above."""
    sys.setrecursionlimit(20000)
    nodesfile, in_treedir, in_seqdir, out_seqdir = get_args()
    nodeslines = mod.get_file_data(nodesfile)
    i = 0
    for line in nodeslines:
        subtree = tm.extract_treenode(line, in_treedir)
        if subtree:
            write_seqs(subtree, in_seqdir, out_seqdir, i, line)
            print()
        i += 1


if __name__ == "__main__":
    main()
