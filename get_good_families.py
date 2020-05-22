#!usr/bin/env python3.6
"""Find the families with a motif we are looking for for further analyses."""

import newick
import sys
import my_module as mod
import tree_module as tm

def get_args():
    """Get the files from the user arguments."""
    if len(sys.argv) == 6:
        return sys.argv[1:7]
    elif len(sys.argv) == 5:
        args = sys.argv[1:6]
        args.append("_tree.txt")
        return args
    else:
        print("USAGE:  python get_good_families.py species_tree orthogruops.tsv trees_directory outfile gene_tree_extention(default = _tree.txt)")
        exit()



def main():
    """Write all gene trees appropriate to a file."""
    species_treefile, orthogroups_file, gene_tree_dir, outfile, extention = get_args()


    #species tree
    species_tree = tm.read_tree(species_treefile)
    sp_nnodes = tm.calculate_nnodes(species_tree)
    species_tree_nodes = tm.get_all_nodes(species_tree, sp_nnodes)
    
    #Gene trees
    #This is the old filterred list
    #candidates = mod.get_file_data("candidates")
    lines = mod.get_file_data(orthogroups_file)
    candidates = []
    for line in lines[1:]:
        candidates.append(line.split("\t")[0])
    
    out = open(outfile, "w")
    for family_name in candidates:
        sys.stderr.write(family_name + "\n")
        #print(tree_dir + "/" + family_name + extention)
        try:
            tree = tm.read_tree(gene_tree_dir + "/" + family_name + extention)
        except:
            sys.stderr.write(family_name + "not in candidates\n")
            continue
 
        nnodes = tm.calculate_nnodes(tree)
        sys.stderr.write("nodes = " + str(nnodes) + "\n")
        nodes = tm.get_all_nodes(tree, nnodes, [], [])
        backup = newick.dumps(tree)
        for node in nodes:
            sys.stderr.write("node = " + str(node) + "\n")
            tree = newick.loads(backup)
            if tm.is_duplication(node):
                sys.stderr.write("is duplication\n")
                # do the old checky
                sides = node.descendants
                for i in range(len(sides)):
                    if tm.is_species_like(sides[i], species_tree_nodes):
                        #print("blammo")
                        if tm.node_with_deletion(sides[i], sides[(i + 1) % 2]):
                            out.write(family_name + "\t" + node.name + "\n")
                            sys.stderr.write("good node!\n\n")
                        break
        #



if __name__ == "__main__":
    sys.setrecursionlimit(10**6)
    main()
