#!/usr/bin/env python3.7

# Convert bond order information from CCD's data to a more compact format

import sys


AMINO_ACIDS = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
               'SER', 'THR', 'TRP', 'TYR', 'VAL'}

ORDER = {'SING': 1, 'DOUB': 2, 'TRIP': 3}


def set_file(residue_name: str):
    return amino_acids_file if residue_name in AMINO_ACIDS else other_file


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Not enough arguments', file=sys.stderr)
        exit(1)

    amino_acids_filename = sys.argv[2]
    other_filename = sys.argv[3]

    amino_acids_file = open(amino_acids_filename, 'w')
    other_file = open(other_filename, 'w')

    with open(sys.argv[1]) as f:
        data = f.read()

    it = iter(data.split('\n'))

    file = None
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
                        name = line.split()[0]
                        file = set_file(name)
                        print(name, file=file)
                        while not line.startswith('#'):
                            atom1, atom2, b_order = line.split()[1: 4]
                            b_order = ORDER[b_order]
                            print(f'{atom1:s} {atom2:s} {b_order:d}', file=file)
                            line = next(it)
                        print(file=file)
                else:
                    if line.startswith('_chem_comp_bond'):
                        name = line.split()[1]
                        file = set_file(name)
                        print(name, file=file)
                        atom1 = next(it).split()[1]
                        atom2 = next(it).split()[1]
                        b_order = ORDER[next(it).split()[1]]
                        print(f'{atom1:s} {atom2:s} {b_order:d}', file=file)
                        print(file=file)
    except StopIteration:
        pass
    amino_acids_file.close()
    other_file.close()
