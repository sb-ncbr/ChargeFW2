#!/usr/bin/env python3.7

# Convert bond order information from CCD's data to a more compact format

import sys

with open(sys.argv[1]) as f:
    data = f.read()

it = iter(data.split('\n')) 
for line in it: 
    if line.startswith('_chem_comp_bond.pdbx_ordinal'): 
        line = next(it)
        print(line.split()[0])
        while not line.startswith('#'): 
            atom1, atom2, order = line.split()[1: 4]
            if order == 'SING':
                order = 1
            elif order == 'DOUB':
                order = 2
            elif order == 'TRIP':
                order = 3
            else:
                raise RuntimeError(f'Incorrect bond order {order}')
            
            print(f'{atom1:s} {atom2:s} {order:d}')
            line = next(it) 
        print()
