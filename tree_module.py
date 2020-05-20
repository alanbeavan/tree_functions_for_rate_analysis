#!usr/bin/env python3.6
"""Module to parse and interpret trees..

Requires:
copy
io
itertools
newick
(run from general environment anaconda)
"""

import copy
import io
import itertools
import newick
import my_module as mod

def read_tree(filename):
    """Read tree into memory.
    
    NOTE - The tree file must be in NEWICK format

    Arguemnts:
    filename - the location of the tree file.
    
    Requires:
    newick, io.
    """
    with io.open(filename, encoding = 'utf8') as fp:
        tree = newick.load(fp)
    return tree[0]

def calculate_nnodes(tree):
    """Return the number of nodes (internal + tips) for a tree.
    
    Arguments:
    tree - a tree object in the format from the newick module.
    
    Requires:
    newick.
    """
    leaves = tree.get_leaf_names()
    return len(leaves) * 2 - 1

def get_all_nodes(tree, nnodes, nodes = [], parents = []):
    """Get all the nodes from a tree.
    
    NOTE - nodes need to be labelled uniquely or this will not work.
    
    Arguments:
    tree - the tree to visit the nodes of.
    nnodes - the number of nodes of the tree.
    nodes - an initially empty list that is filled as nodes are visited
        and passed back to get_all_nodes.
    parents - an initially empry list containing all the nodes that are
        ancestors of the current nodes. This is used to traverse back 
        through the tree in subsequent calls of get_all_nodes.

    Requires:
    newick.

    Recursively visit all nodes in a tree, adding them to the list of nodes
    if they are not yet in it. If we reach a tip, we go back to the parent
    node and continue. When all nodes have been visitted we return the list.
    """
    if tree not in nodes:
        nodes.append(tree)
    if len(tree.descendants) > 0:
        for node in tree.descendants:
            if node not in nodes:
                parents.append(tree)
                get_all_nodes(node, nnodes, nodes, parents)
    
    if len(nodes) < nnodes:
        if len(parents) > 0:
            get_all_nodes(parents[-1], nnodes, nodes, parents[:len(parents)-1])
    return nodes

def is_duplication(node):
    """Calculate if a node is a duplication.
    
    Arguments:
    node - the node to be tested.

    Requires:
    newick.

    Look at all the tips on each side of the node, removing all the gene
    information so we just know which species are present on each side. If
    there are species present on both side of the node, return 1, as it is a
    duplication. Otherwise, return 0 (it could strictly be a duplication still
    but for the purpose of this program, we can say it isn't. Further checks 
    are implemented later).
    """
    if len(node.descendants) < 2:
            return 0
    tips = []
    for n in node.descendants:
        tips.append(n.get_leaf_names())
    for set in tips:
        for i in range(len(set)):
            set[i] = set[i].split("_")[0]

    for tip in tips[0]:
        if tip in tips[1]:
            return 1
    return 0

def is_species_like(node, species_tree_nodes):
    """Ask if the node of this gene tree reflects the species tree.
    
    Return 1 if the gene each node in the gene tree describes a monophyletic
    clade in the species tree (ie. by removing the gene names, keeping only
    the species name). I don't think it's possible for this to be true and the
    gene tree to be describing para/polyphyletic clades.

    Otherwise return 0.

    Arguments:
    node - a gene tree node
    species_tree - the underlying species tree
    """
    for n in get_all_nodes(node, calculate_nnodes(node), [], []):
        if len(n.descendants) > 0:
            species = n.get_leaf_names()
            for i in range(len(species)):
                species[i] = species[i].split("_")[0]
            if len(list(set(species))) != len(species):
                return 0 #This means there are duplicate nodes in the tree
            good = 0
            for sp_node in species_tree_nodes:
                if sorted(species) == sorted(sp_node.get_leaf_names()):
                    good = 1
                    break
            if not good:
                return 0
    return 1 #Success

def is_monophyletic(taxa, nodes):
    """Determine if a set of tips are monophyletic.
    
    Arguments:
    taxa - a list of tips
    nodes - a list of nodes in the tree
    
    Requires:
    newick.
    
    Loop through the nodes and ask if it contains all of the tips exclusively.
    Return 1 if any do. Otherwise return 0.
    """
    for node in nodes:
        species = node.get_leaf_names()
        if sorted(species) == sorted(taxa):
            return node
    return 0



def node_with_deletion(old_node, new_node):
    """Ask if a node can be descibed by deleting a node from another node.
    
    I apreaciate that is a confusing sentence.
    
    Arguments:
    old_node - the node that we hope to delete a node from.
    new_node - the node we hope can be described by this.

    Requires:
    Newick

    return 1 is new_node can be descibed by deleting a node from old_node.
    else return 0. First remove all the gene information from the leaf names.
    Then ...  MAYBE make a copy of old, remove a node, remove all node labels
    from both an old copy the new (actually do this before the removing of
    nodes), then see if the nodes are the same. I wonder if the way they are
    presented makes a difference. I hope not.
    """
    for leaf in old_node.get_leaves():
        leaf.name = leaf.name.split("_")[0]
    for leaf in new_node.get_leaves():
        leaf.name = leaf.name.split("_")[0]
    #new_node.remove_internal_names()
    #old_node.remove_internal_names() 

    old_nodes = get_all_nodes(old_node, calculate_nnodes(old_node), [], [])
    new_nodes = get_all_nodes(new_node, calculate_nnodes(new_node), [], [])
    #backup = newick.dumps(old_node)
    
    old_taxa = old_node.get_leaf_names()
    new_taxa = new_node.get_leaf_names()
    old_exclusive = []
    for taxon in old_taxa:
        if taxon not in new_taxa:
            old_exclusive.append(taxon)
    to_remove = is_monophyletic(old_exclusive, old_nodes)
    if to_remove:
        old_node.prune([to_remove])
        old_node.remove_redundant_nodes()
        old_nodes = get_all_nodes(old_node, calculate_nnodes(old_node), [], [])

        for i in range(len(new_taxa)):
            for perm in itertools.permutations(new_taxa, i):
                if is_monophyletic(list(perm), new_nodes) and not is_monophyletic(list(perm), old_nodes):
                    return 0
        return 1
    else:
        return 0
    
    #print(old_copy.ascii_art())
    #for 
    

def main():
    """Just tests."""
    #tree=newick.loads("((Dmel_CG7377:5.71073e-07,Dsim_Dsim\GD12794:0.426781)n1:0.0026795,(Dmel_CG33268:0.022453,(Dsim_Dsim\GD14314:0.015169,Dsec_Dsec\GM25283:0.029888)n3:0.079816)n2:0.0026795)n0;")
    #tree=newick.loads("(((((((((((((((((((((((((A,B)N1,C)N2,D)N3,E)N4,F)N5,G)N6,H)N7,I)N8,J)N9,K)N10,L)N11,M)N12,N)N13,O)N14,P)N15,Q)N16,R)N17,S)N18,T)N19,U)N20,V)N21,W)N22,X)N23,Y)N24,Z)N25;")
    
    #Example tree
    species_tree = read_tree("../proteomes_repeats_removed/OrthoFinder/Results_Jan20/Species_Tree/SpeciesTree_rooted_node_labels.txt")
    sp_nnodes = calculate_nnodes(species_tree)
    species_tree_nodes = get_all_nodes(species_tree, sp_nnodes)
    
    
    tree = read_tree("../proteomes_repeats_removed/OrthoFinder/Results_Jan20/Resolved_Gene_Trees/OG0000028_tree.txt")
    print(tree.ascii_art())
    #tree = read_tree("../proteomes_repeats_removed/OrthoFinder/Results_Jan20/Resolved_Gene_Trees/OG0012151_tree.txt")
 
    nnodes = calculate_nnodes(tree)
    nodes = get_all_nodes(tree, nnodes, [], [])
    backup = newick.dumps(tree)
    for node in nodes:
        tree = newick.loads(backup)
        if is_duplication(node):
            #print(node)
            # do the old checky
            sides = node.descendants
            for i in range(len(sides)):
                if is_species_like(sides[i], species_tree_nodes):
                    #print(sides[i].ascii_art())
                    #print("blammo")
                    if node_with_deletion(sides[i], sides[(i + 1) % 2]):
                        print("it's good")
                    break
                
            
    
if __name__ == "__main__":
    main()
