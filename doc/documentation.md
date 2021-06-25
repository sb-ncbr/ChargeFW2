# ChargeFW2 usage

The following tutorial assumes you have ChargeFW2 installed in `/opt/chargefw2` directory.
See the main [Readme](../README.md) for installation instructions.

The [Python bindings](ChargeFW2%20-%20tutorial.pdf) are documented in a separate Jupyter Notebook.

## Describe a set of molecules 
The application provides several modes of operation. The simplest one is to describe the set of molecules:

```shell
$ /opt/chargefw2/bin/chargefw2 --mode info --input-file molecules.sdf

Number of molecules: 47
Number of atoms: 1531
N  plain  *   : 101
O  plain  *   : 118
C  plain  *   : 545
H  plain  *   : 751
S  plain  *   : 16
```

The output just prints the number of atoms, atomic types and molecules in the particular file.
The phrase `plain *` means that the atomic types are based on element only.
Note that several input formats are supported: SDF with MOL V2000 and V3000, PDB, Mol2 and mmCIF.

## Suitable methods

Not every method or parameter set is suitable for a given set of molecules. To find out what combinations are available,
use the `suitable-methods` mode.

```shell
$ /opt/chargefw2/bin/chargefw2 --mode suitable-methods --input-file molecules.sdf
eem EEM_00_NEEMP_ccd2016_npa.json EEM_00_NEEMP_ccd2016_npa2.json EEM_05_NEEMP_ccd2016_mpa.json EEM_05_NEEMP_ccd2016_mpa2.json EEM_10_Cheminf_b3lyp_aim.json EEM_10_Cheminf_b3lyp_mpa.json EEM_10_Cheminf_b3lyp_npa.json EEM_10_Cheminf_hf_aim.json EEM_10_Cheminf_hf_mpa.json EEM_10_Cheminf_hf_npa.json
qeq QEq_00_original.json
eqeq
sqe SQE_10_Schindler2021_CCD_gen.json
sqeq0 SQEq0_10_Schindler2021_CCD_gen.json
sqeqp SQEqp_10_Schindler2021_CCD_gen.json
mgc
denr DENR_00_from_QEq.json
tsef TSEF_00_from_QEq.json
veem
formal
```

The output lists the methods (always the first word on the line) and the usable parameter sets for each method.
If no parameter set is mentioned, the method does not use parameters.

## Computation of partial atomic charges
The main purpose of ChargeFW2 is the computation of partial atomic charges via a `charges` mode.
This mode requires several additional parameters, mainly the specification of the method (`--method`)
and the output directory (`--chg-out-dir`) where the charges should be stored.


```shell
$ /opt/chargefw2/bin/chargefw2 --mode charges --input-file molecules.sdf --method eem --chg-out-dir /tmp
Method: EEM
Best parameters found: EEM_00_NEEMP_ccd2016_npa.json
Parameters: Racek 2016 (ccd2016_npa)
Common parameters
kappa: 0.512
Atom parameters
C  hbo    1   :  2.722  0.640 
C  hbo    2   :  2.767  0.651 
C  hbo    3   :  2.694  0.678 
O  hbo    1   :  3.150  1.058 
O  hbo    2   :  3.049  0.841 
H  hbo    1   :  2.547  1.164 
N  hbo    1   :  2.975  0.908 
N  hbo    2   :  2.889  0.665 
N  hbo    3   :  3.024  1.424 
F  hbo    1   :  2.998  0.998 
S  hbo    1   :  2.651  0.490 
S  hbo    2   :  2.647  0.451 
Cl hbo    1   :  2.703  1.154 
P  hbo    2   :  2.293  0.576 
Br hbo    1   :  2.626  1.111 

Number of unclassified molecules: 0

Number of molecules: 47
Number of atoms: 1531
C  hbo    1   : 295
C  hbo    2   : 250
O  hbo    1   : 51
O  hbo    2   : 67
H  hbo    1   : 751
N  hbo    1   : 64
N  hbo    2   : 37
S  hbo    1   : 11
S  hbo    2   : 5

```

Optinally, user can specify a parameter set (`--par-file`). If it is not provided, the ChargeFW2 selects
the first one from a particular methods based on the list in the previous section. The parameters are stored in the
`/opt/chargefw2/share/parameters/` directory. See the explicit parameter setting:

```shell
$ /opt/chargefw2/bin/chargefw2 --mode charges --input-file molecules.sdf --method eem --chg-out-dir /tmp --par-file /opt/chargefw2/share/parameters/EEM_20_Bult2002_npa.json
Method: EEM
Parameters: Bultinck 2002 (NPA)
Common parameters
kappa: 0.529
Atom parameters
C  plain  *   :  0.312  0.673 
F  plain  *   :  1.440  3.241 
H  plain  *   :  0.037  0.714 
N  plain  *   :  0.494  0.782 
O  plain  *   :  0.994  1.443 

Number of unclassified molecules: 11

Number of molecules: 36
Number of atoms: 1192
C  plain  *   : 438
H  plain  *   : 593
N  plain  *   : 73
O  plain  *   : 88
```

Note that this parameter set was not able to cover 11 out of the 47 molecules.
The reason is that these parameters does not have values for a sulphur,
so the molecules containing sulphur had to be omitted from the calculation.

Finally, the charges are stored in the `/tmp` directory with a name prefixed by the input-file name:

```shell
$ ls /tmp/molecules.sdf*
/tmp/molecules.sdf.mol2  /tmp/molecules.sdf.txt
```

Generally, two types of files are generated. If the input is recognized to contain a protein structure
(either from PDB or mmCIF file), the PQR file is generated. Otherwise, Mol2 is the standard output format for storing the charges.

In both cases, also a plain text file is available:

```shell
$ head -4 /tmp/molecules.sdf.txt
NSC_100000
-0.47112 -0.64759 -0.68543 -0.41867 -0.33283 -0.34347 -0.58064 -0.53754 0.26420 -0.50466 -0.37041 0.17366 -0.02944 -0.37767 -0.48880 0.20385 0.40290 -1.06418 0.46233 0.44876 0.43404 0.61994 0.48120 0.45123 0.43844 0.45491 0.60701 0.75376 0.65621 
NSC_100013
-0.67515 -0.64035 -0.67342 -0.42289 -0.32961 -0.32024 -0.60965 -0.59315 0.38509 -0.55511 -0.35609 0.22204 -0.03048 -0.38049 0.20857 0.45355 -0.48655 0.48109 0.44285 0.46264 0.73712 0.58084 0.54661 0.45527 0.47393 0.62358
```
