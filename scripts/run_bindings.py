import chargefw2_python
import argparse

from pathlib import Path



def cli():

    parser = argparse.ArgumentParser()

    parser.add_argument('cif_file', type=str)

    parser.add_argument('-d', '--cif-dir', required=False, type=Path, 
        default=Path('../../../data/proteins/fw2/input'))

    parser.add_argument('-m','--method', required=False, type=str, default='eem')

    parser.add_argument('-p', '--parameter-set', required=False, type=str, 
        default='EEM_10_Cheminf_b3lyp_aim')

    args = parser.parse_args()

    args.cif_file = args.cif_dir / args.cif_file
    args.out_dir = args.cif_dir.parent / "output"

    args.cif_dir.mkdir(0o744, parents=True, exist_ok=True)
    args.out_dir.mkdir(0o744, parents=True, exist_ok=True)

    if args.cif_file.suffix != '.cif':
        raise ValueError('cif_file parameter must end with a .cif extension.')

    return args




if __name__=="__main__":

    args = cli()

    protein = chargefw2_python.Molecules(str(args.cif_file))

    charges = chargefw2_python.calculate_charges(protein, args.method, args.parameter_set)

    chargefw2_python.write_cif(protein, charges, str(args.cif_file), str(args.cif_dir))