#!usr/bin/env python3.6
"""Check which trees have the appropriate topology for analysis."""

import glob
import sys
import newick
import tree_module as tm
import my_module as mod

def get_args():
    """Get user arguments."""
    if len(sys.argv) == 4:
        return sys.argv[1:]
    print("\nUSAGE python assess_iqtree_topoliges.py treefiles_dir species_tree output_file\n")
    exit()

def main():
    """
    Write appropriate trees to a file.
    
    For each tree:
    get the nodes of the tree
    look at the duplication nodes
    for each side:
    look at one side to see if it is species-like
    if it is, check if the other one is node with deletion
    if success, add to the file or whatever.
    """

    treedir, species_tree, out_file = get_args()
    species_tree = tm.read_tree(species_tree)
    sp_tree_nodes = tm.get_all_nodes(species_tree, tm.calculate_nnodes(species_tree), [], [])
    treefiles = glob.glob(treedir + "/*rooted")
    good_trees = []
#    for i in range(len(treefiles)):
    for i in [16]:
        tree = tm.read_tree(treefiles[i])
        nodes = tm.get_all_nodes(tree, tm.calculate_nnodes(tree), [], [])    
        backup = newick.dumps(tree)
        for node in tree.descendants:
            tree = newick.loads(backup)
            if tm.is_duplication(node):
                for j in range(2):
                    if tm.is_species_like(node.descendants[j], sp_tree_nodes):
                        print("species_like_side")
                        if tm.node_with_deletion(node.descendants[j], node.descendants[(j + 1) % 2]):
                            good_trees.append(treefiles[i])
    out = open(out_file, "w")
    out.write("\n".join(good_trees))
    out.close()

        
    

if __name__ == "__main__":
    main()
