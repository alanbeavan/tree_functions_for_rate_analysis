#!usr/bin/env python3.6
"""Find gene families with an internal duplication > 40% support < 100%."""

import my_module as mod

def main():
    """
    Run through the lines printing gene families if they have
    the above specification.
    """
    families = []
    for line in mod.get_file_data("Gene_Duplication_Events/Duplications.tsv")[1:]:
        fields = line.split()
        if fields[4] == "Non-Terminal" \
        and float(fields[3]) < 1.0 \
        and float(fields[3]) > 0.4:
            if fields[0] not in families:
                families.append(fields[0])
    print("\n".join(families))


if __name__ == "__main__":
    main()
