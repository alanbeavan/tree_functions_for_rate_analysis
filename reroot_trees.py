#!usr/bin/env python3.6
"""
Reroot gene trees according to the root in orthofinder gene tree.

I expect this will look like:
    for each gene tree:
        Extract the tip names from the orthofinder gene tree.
        Modify them according to the iqtree renaming.
        Use newick or ete3 or something to make these taxa the root.
        Write the rooted to tree to a file.

"""

import re
import os
import sys
import newick
import ete3
import dendropy as dp
import tree_module as tm
import my_module as mod

def get_args():
    """Get user arguments or remind user of how to run."""
    if len(sys.argv) == 3:
        return sys.argv[1:]
    print("USAGE: python reroot_trees.py orthofinder_tree_dir iqtree_dir\n\n")
    exit()

def get_og_outgroup(tree):
    """Get the tips from one side of a tree.
    
    RETURNS:
    list of tip names.
    """
    return tree.descendants[0].get_leaf_names()

def translate_tiplabels(labels):
    """Translate the tip labels as iqtree does.

    RETURNS:
    list of translated tip labels.
    """
    translate_dict = {"|": "_"}
    newlabs = []
    for label in labels:
        for key, value in translate_dict.items():
            pattern = re.compile(re.escape(key))
            newlabs.append(re.sub(pattern, value, label))
    return newlabs

def reroot_tree(tree, outgroup):
    """Reroot the tree so that the outgroup is the outgroup.

    Arguments:
    tree - a tree in python newick format
    outgroup - a list of taxa that ought to make up the entirety of the
        outgroup

    Returns:
    rooted_tree - a tree with the correct outgroup in a format.

    Requires:
    Dendropy
    """
    newick_string = newick.dumps(tree)
    tree = ete3.Tree(newick_string, format = 1)
    print(tree)
    try:
        if len(outgroup) == 1:
            tree.set_outgroup(outgroup[0])
        else:
            mrca = tree.get_common_ancestor(outgroup)
            tree.set_outgroup(mrca)
    except:
        taxa = []
        for leaf in tree:
            taxa.append(leaf.name)
        outgroup = list(set(taxa) - set(outgroup))
        #print(outgroup)
        if len(outgroup) == 1:
            tree.set_outgroup(outgroup[0])
        else:
            mrca = tree.get_common_ancestor(outgroup)
            tree.set_outgroup(mrca)
    return tree

def write_tree(tree, outdir, index):
    """Write new rooted tree to a file according to the oudir and index."""
    #tree_string = tree.write(tree, format = 1)
    #out = open(outdir + "/full_" + str(index) + ".rooted", "w")
    #out.write(tree_string)
    #out.close()
    tree.write(format = 1, outfile = outdir + "/full_" + str(index) + ".rooted")

def main():
    """Call functions to achieve the aims of the program."""
    orthofinder_trees, iqtree_trees = get_args()
    n_trees = len(os.listdir(orthofinder_trees))
    i = 0
    while i < n_trees:
        try:
            print(str(i))
            og_tree = tm.read_tree(orthofinder_trees + "/tree_" + str(i) + ".txt")
            iqtree_tree = tm.read_tree(iqtree_trees + "/full_" + str(i) + \
                                   ".treefile")
            outgroup_taxa = get_og_outgroup(og_tree)
            outgroup_taxa = translate_tiplabels(outgroup_taxa)
            
            rooted_tree = reroot_tree(iqtree_tree, outgroup_taxa)

            write_tree(rooted_tree, iqtree_trees, i)

        except:
            print(str(i) + "didn't work for some reason")

        i += 1

if __name__ == "__main__":
    main()
