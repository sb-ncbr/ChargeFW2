#!/usr/bin/env python3

# Convert bond order information from CCD's data to a more compact format

import argparse
import gemmi

AMINO_ACIDS = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
    'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
    'TYR', 'VAL',
}

# Map gemmi bond types to integer orders
BOND_ORDER = {
    gemmi.BondType.Single: 1,
    gemmi.BondType.Double: 2,
    gemmi.BondType.Triple: 3,
}


def main():
    parser = argparse.ArgumentParser(
        description='Convert CCD bond order information to a compact format.'
    )
    parser.add_argument('ccd_file', help='CCD file (components.cif or components.cif.gz)')
    parser.add_argument('amino_acids_out', help='Output file for amino-acid residues')
    parser.add_argument('other_out', help='Output file for non-amino-acid residues')

    args = parser.parse_args()
    doc = gemmi.cif.read(args.ccd_file)

    with open(args.amino_acids_out, 'w') as aa_out, open(args.other_out, 'w') as other_out:
        for block in doc:
            cc = gemmi.make_chemcomp_from_block(block)
            resname = cc.name

            out_file = aa_out if resname in AMINO_ACIDS else other_out
            print(resname, file=out_file)

            for bond in cc.rt.bonds:
                order = BOND_ORDER.get(bond.type)
                atom1 = bond.id1.atom
                atom2 = bond.id2.atom
                print(f'{atom1} {atom2} {order}', file=out_file)

            print(file=out_file)  # blank line after each residue


if __name__ == '__main__':
    main()
