import subprocess
import json
import os
import argparse
from pathlib import Path
from collections import Counter, defaultdict

CHARGEFW2_DIR = Path('../chargefw2')
CHARGEFW2_BIN = CHARGEFW2_DIR / 'bin/chargefw2'
CHARGEFW2_PARAMS_DIR = CHARGEFW2_DIR / 'share/parameters'


def cli():

    parser = argparse.ArgumentParser()

    parser.add_argument('cif_file', type=str)

    parser.add_argument('-d', '--cif-dir', required=False, type=Path, 
        default=Path('./input'))

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



def calculate(source, charge_out_dir, method_name='eem', parameters_name='EEM_10_Cheminf_b3lyp_aim.json'):
    logfile = charge_out_dir / 'logs/computations'

    method_file = CHARGEFW2_DIR / 'share/methods.json'

    with method_file.open() as f:
        method_data = json.load(f)

    method_data = method_data["methods"]

    args = [
        str(CHARGEFW2_BIN), '--mode', 'charges', '--method', method_name,
        '--input-file', str(source), '--chg-out-dir', str(charge_out_dir), '--read-hetatm', '--log-file', str(logfile),
        '--permissive-types'
    ]
    if next(m for m in method_data if m['internal_name'] == method_name)['has_parameters']:
        args.extend(['--par-file', str(CHARGEFW2_PARAMS_DIR / parameters_name)])

    calculation = subprocess.run(args, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    print(' '.join(calculation.args))
    return calculation


def calculate_charges(method_name, parameters_name, tmp_dir):
    structures: Dict[str, str] = {}
    charges: Dict[str, str] = {}
    formats: Dict[str, str] = {}
    logs: Dict[str, str] = {}

    for file in os.listdir(os.path.join(tmp_dir, 'input')):
        res = calculate(method_name, parameters_name, os.path.join(tmp_dir, 'input', file),
                        os.path.join(tmp_dir, 'output'))

        stderr = res.stderr.decode('utf-8')

        with open(os.path.join(tmp_dir, 'logs', f'{file}.stdout'), 'w') as f_stdout:
            f_stdout.write(res.stdout.decode('utf-8'))

        with open(os.path.join(tmp_dir, 'logs', f'{file}.stderr'), 'w') as f_stderr:
            f_stderr.write(stderr)

        if stderr.strip():
            logs['stderr'] = stderr

        if res.returncode:
            flash('Computation failed. See logs for details.', 'error')

        _, ext = os.path.splitext(file)
        ext = ext.lower()

        tmp_structures: Dict[str, str] = {}
        with open(os.path.join(tmp_dir, 'input', file)) as f:
            if ext == '.sdf':
                if get_MOL_versions(os.path.join(tmp_dir, 'input', file)) == {'V2000'}:
                    tmp_structures.update(parse_sdf(f))
                    fmt = 'SDF'
                else:
                    tmp_structures.update(convert_to_mmcif(f, 'sdf', file))
                    fmt = 'mmCIF'
            elif ext == '.mol2':
                tmp_structures.update(convert_to_mmcif(f, 'mol2', file))
                fmt = 'mmCIF'
            elif ext == '.pdb':
                tmp_structures.update(parse_pdb(f))
                fmt = 'PDB'
            elif ext == '.cif':
                tmp_structures.update(parse_cif(f))
                fmt = 'mmCIF'
            else:
                raise RuntimeError(f'Not supported format: {ext}')

        for s in tmp_structures:
            formats[s] = fmt

        structures.update(tmp_structures)

        with open(os.path.join(tmp_dir, 'output', f'{file}.txt')) as f:
            charges.update(parse_txt(f))
    return charges, structures, formats, logs


def get_suitable_methods(filepath):

    suitable_methods = Counter()
    args = [
        str(CHARGEFW2_BIN), 
        '--mode', 'suitable-methods', 
        '--read-hetatm',
        '--permissive-types', 
        '--input-file', str(filepath)
    ]

    calculation = subprocess.run(args, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    print(' '.join(calculation.args))
    if calculation.returncode:
        output = calculation.stderr.decode('utf-8').strip()
        print(output)
        error = output.split('\n')[-1]
        raise RuntimeError(error)

    for line in calculation.stdout.decode('utf-8').splitlines():
        if not line.strip():
            continue
        method, *parameters = line.strip().split()
        print(method, parameters)
        if not parameters:
            suitable_methods[(method,)] += 1
        else:
            for p in parameters:
                print('pair: ', method, p)
                suitable_methods[(method, p)] += 1

    # This is if we are looping over a directory and want to make sure we pick a pair that is valid
    # for all proteins in the directory. 
    # all_valid = [
    #   pair for pair in suitable_methods 
    #   if suitable_methods[pair] == len(directory.iterdir())
    # ]

    # remove duplicate methods
    methods = {data[0] for data in suitable_methods}

    parameters = defaultdict(list)
    for pair in suitable_methods:
        if len(pair) == 2:
            parameters[pair[0]].append(pair[1])

    for method, params in parameters.items():
        print(f'method: {method}, params: {params}')

        

# TODO add a cli
if __name__=="__main__":

    method = "eem"

    # parameter = 'EEM65_Ionescu2013_mpa_pcm.json'
    parameter = 'EEM_10_Cheminf_b3lyp_aim.json'

    args = cli()

    if not args.parameter_set.endswith('.json'):
        args.parameter_set += ".json"

    # TODO move to a suitable method cli script
    # get_suitable_methods(source)

    #run a calculation
    print(calculate(str(args.cif_file), args.out_dir), args.method, args.parameter_set)