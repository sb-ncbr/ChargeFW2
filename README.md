# ChargeFW2

Application for computing partial atomic charges using selected empirical methods.
ChargeFW2 is the computational core of [Atomic Charge Calculator](https://acc2.ncbr.muni.cz).

See the [short description](https://acc2.ncbr.muni.cz/static/assets/methods.pdf) of implemented methods. 

## Compilation requirements

To compile ChargeFW2, several requirements have to be installed first.
The compilation was tested on Ubuntu 25.10 with these versions of libraries:

- [CMake](https://cmake.org/) 3.31
- [GCC](https://gcc.gnu.org/) 15.2
- [Boost](https://www.boost.org/) 1.83
- [Eigen](http://eigen.tuxfamily.org) 3.4
- [nanoflann](https://github.com/jlblancoc/nanoflann) 1.7
- [JSON for Modern C++](https://github.com/nlohmann/json) 3.11
- [GEMMI](https://github.com/project-gemmi/gemmi) 0.7.4
- [pybind11](https://github.com/pybind/pybind11) 2.11

Other versions of the libraries might work too, but it was not tested. See the `Dockerfile` as the formal description
of the dependencies.

## Installation
After downloading and unpacking the sources, run the following in the ChargeFW2 directory:

```bash
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=<WHERE-TO-INSTALL> -DPYTHON_MODULE=ON -DCMAKE_BUILD_TYPE=Release
make
make install
```

## Docker

Rather than installing all dependencies, you can run ChargeFW2 directly inside a Docker container.

### Using the pre-built image

A ready-to-use image is available on Docker Hub:

```bash
# Explicit pull of the image (optional)
docker pull sbncbr/chargefw2
# Check that it works
docker run --rm sbncbr/chargefw2 --help
```

### Local build
You can also build the Docker image locally from the provided `Dockerfile`.
This version will be optimised for your hardware as opposed to the pre-built general image from Docker Hub.

```bash
docker build -t sbncbr/chargefw2:local .
docker run --rm  sbncbr/chargefw2:local --help
```

### Running

Suppose you want to compute partial atomic charges for a molecule stored in
`molecule.sdf` located at `/home/user/data/molecules` on your host machine,
and you want ChargeFW2 to write the results into `/home/user/data/charges`.

To make these host directories accessible inside the container, bind-mount them
and expose them as `/in` and `/out` inside the container:

```bash
INPUT_DIR="/home/user/data/molecules"
OUTPUT_DIR="/home/user/data/charges"

docker run --rm \
  -v "$INPUT_DIR:/in" \
  -v "$OUTPUT_DIR:/out" \
  sbncbr/chargefw2 \
  --mode charges \
  --input-file /in/molecule.sdf \
  --chg-out-dir /out
```

## Usage

The [documentation](doc/documentation.md) for the application and its [Python bindings](doc/ChargeFW2%20-%20tutorial.pdf) is located in the [doc](doc) folder.

## How to cite
If you found ChargeFW2 or Atomic Charge Calculator helpful, please cite: [Raček, T., Schindler, O., Toušek, D., Horský, V., Berka, K., Koča, J., & Svobodová, R. (2020). Atomic Charge Calculator II: web-based tool for the calculation of partial atomic charges. Nucleic Acids Research](https://doi.org/10.1093/nar/gkaa367).
