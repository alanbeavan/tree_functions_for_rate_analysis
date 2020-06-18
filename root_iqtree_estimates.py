#!usr/bin/env python3.6
"""Root gene trees according to the ougroups file."""

import sys
import ete3
import my_module as mod

def get_args():
    """Get user input."""
    if len(sys.argv) == 3:
        return sys.argv[1:]
    else:
        print("USAGE: python root_iqtree_estimates.py treedir outgroups\n")
        exit()


def main():
    """Read in trees and outgroups, write rooted trees."""
    treedir, outgroups_file = get_args()
    outgroups = mod.get_file_data(outgroups_file)
    for i in range(len(outgroups)):
        tree = ete3.Tree(treedir + "/full_" + str(i) + ".treefile", format = 1)
        tree.set_outgroup(outgroups[i])
        out = open(treedir + "/full_" + str(i) + ".rooted", "w")
        out.write(tree.write())
        out.close()



if __name__ == "__main__":
    main()
