#!usr/bin/env python3.6
"""Find the families with a motif we are looking for for further analyses."""

import newick
import sys
import my_module as mod
import tree_module as tm

def main():
    """Write all gene trees appropriate to a file."""
    #species tree
    species_tree = tm.read_tree("../proteomes_repeats_removed/OrthoFinder/Results_Jan20/Species_Tree/SpeciesTree_rooted_node_labels.txt")
    sp_nnodes = tm.calculate_nnodes(species_tree)
    species_tree_nodes = tm.get_all_nodes(species_tree, sp_nnodes)
    
    #Gene trees
    #This is the old filterred list
    #candidates = mod.get_file_data("candidates")
    lines = mod.get_file_data("../proteomes_repeats_removed/OrthoFinder/Results_Jan20/Orthogroups/Orthogroups.tsv")
    candidates = []
    for line in lines[1:]:
        candidates.append(line.split("\t")[0])
    

    for family_name in candidates:
        sys.stderr.write(family_name + "\n")
        try:
            tree = tm.read_tree("../proteomes_repeats_removed/OrthoFinder/Results_Jan20/Resolved_Gene_Trees/" + family_name + "_tree.txt")
 
            nnodes = tm.calculate_nnodes(tree)
            nodes = tm.get_all_nodes(tree, nnodes, [], [])
            backup = newick.dumps(tree)
            for node in nodes:
                tree = newick.loads(backup)
                if tm.is_duplication(node):
                    #print(node)
                    # do the old checky
                    sides = node.descendants
                    for i in range(len(sides)):
                        if tm.is_species_like(sides[i], species_tree_nodes):
                            #print("blammo")
                            if tm.node_with_deletion(sides[i], sides[(i + 1) % 2]):
                                print(family_name + "\t" + node.name)
                            break
        except:
            continue
 



if __name__ == "__main__":
    main()
