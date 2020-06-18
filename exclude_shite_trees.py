#!usr/bin/env python3.6
"""Write a file of trees without extremely low branch lengths.

The idea is that these are the trees that orthofinder estimated incorrectly.
"""

import re
import sys
import my_module as mod

def get_args():
    """Get user arguments."""
    if len(sys.argv) == 4:
        return sys.argv[1:]
    else:
        print("\nUSAGE python exclude_shite_trees.py treedir nodes_file outfile\n")

def extract_branch_lengths(tree):
    """Get all the branch lengths on a tree.

    Returns only a list of the lengths - not ordered or anything fancy."""
    lengths = re.findall(r':[0-9\.e\-]+', tree)
    float_lengths = []
    for length in lengths:
        float_lengths.append(float(length[1:]))
    return float_lengths



def main():
    """Do the above."""
    treedir, node_file, outfile = get_args()
    nodelines = mod.get_file_data(node_file)
    nodes = []
    for line in nodelines:
        nodes.append(line.split()[1])
    
    i = 0
    while i < 578:
        flag = 1
        try:
            tree_string = mod.get_file_data(treedir + "/full_"
                                            + str(i) + ".rooted")[0]
        except:
            i += 1
            continue
        branch_lengths = extract_branch_lengths(tree_string)
        for bl in branch_lengths:
            if bl < 0.01:
                flag = 0
                break
        if flag == 1:
            out = open(outfile, "a")
            out.write("full_" + str(i) + ".rooted" + "\t" + nodes[i] + "\n")
            out.close()
        i += 1



if __name__ == "__main__":
    main()
