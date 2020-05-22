#!usr/bin/env python3.6
"""Find gene families with an internal duplication > 40% support < 100%."""

import sys
import my_module as mod

def get_args():
    """File locations etc."""
    if len(sys.argv) == 3:
        return sys.argv[1:3]
    else:
        print("USAGE python3 duplication_file outfile\n\n")
        exit()

def main():
    """
    Run through the lines printing gene families if they have
    the above specification.
    """
    dups_file, outfile = get_args()

    families = []
    for line in mod.get_file_data(dups_file)[1:]:
        fields = line.split()
        if fields[4] == "Non-Terminal" \
        and float(fields[3]) < 1.0 \
        and float(fields[3]) > 0.4:
            if fields[0] not in families:
                families.append(fields[0])
    out = open(outfile, "w")
    out.write("\n".join(families))
    out.close()

if __name__ == "__main__":
    main()
