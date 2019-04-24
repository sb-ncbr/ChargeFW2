#!/usr/bin/env python3.7

# Convert bond order information from CCD's data to a more compact format

import sys

with open(sys.argv[1]) as f:
    data = f.read()

order = {'SING': 1, 'DOUB': 2, 'TRIP': 3}


it = iter(data.split('\n'))

try:
    while True:
        line = next(it)
        if line.startswith('#'):
            line = next(it)
            if line.startswith('loop_'):
                line = next(it)
                if line.startswith('_chem_comp_bond'):
                    while line.startswith('_'):
                        line = next(it)
                    print(line.split()[0])
                    while not line.startswith('#'):
                        atom1, atom2, b_order = line.split()[1: 4]
                        b_order = order[b_order]
                        print(f'{atom1:s} {atom2:s} {b_order:d}')
                        line = next(it)
                    print()
            else:
                if line.startswith('_chem_comp_bond'):
                    print(line.split()[1])
                    atom1 = next(it).split()[1]
                    atom2 = next(it).split()[1]
                    b_order = order[next(it).split()[1]]
                    print(f'{atom1:s} {atom2:s} {b_order:d}')
                    print()
except StopIteration:
    pass
