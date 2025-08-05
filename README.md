# ChargeFW2

Application for computing partial atomic charges using selected empirical methods.
ChargeFW2 is the computational core of [Atomic Charge Calculator II](https://acc2.ncbr.muni.cz).

See the [short description](https://acc2.ncbr.muni.cz/static/assets/methods.pdf) of implemented methods. 

## Compilation requirements

To compile ChargeFW2, several requirements have to be installed first.
The compilation was tested on Ubuntu 25.04 with these versions of libraries:

- [CMake](https://cmake.org/) 3.31
- [GCC](https://gcc.gnu.org/) 14
- [Boost](https://www.boost.org/) 1.83
- [Eigen](http://eigen.tuxfamily.org) 3.4
- [fmt](https://fmt.dev) 10.1
- [nanoflann](https://github.com/jlblancoc/nanoflann) 1.7
- [JSON for Modern C++](https://github.com/nlohmann/json) 3.11
- [GEMMI](https://github.com/project-gemmi/gemmi) 0.7.3
- [pybind11](https://github.com/pybind/pybind11) 2.11

Other versions of the libraries might work too, but it was not tested. See the `Dockerfile` as the formal description
of dependencies.

## Installation
After downloading and unpacking the sources, run the following in the ChargeFW2 directory:

```bash
$ mkdir build
$ cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=<WHERE-TO-INSTALL>
$ make
$ make install
```

## Docker
Rather than installing all dependencies, you can run ChargeFW2 directly inside a Docker container.

```bash
# Build docker image
$ docker build -t chargefw2 .

# Run docker container
$ docker run -it --rm --entrypoint bash chargefw2
```

### CLI workflow

The Docker container can be used in CLI workflows. However, since containers are run in isolated environments and don't have access to your local host files by default, you will need to create a [volume](https://docs.docker.com/storage/volumes/) or a [bind-mount](https://docs.docker.com/storage/bind-mounts/), which will allow sharing host files with the container.

In this example we bind-mount the current working directory into the container. This allows us to specify input files and output directories using relative paths. This is the most similar to running ChargeFW2 natively.

```bash
$ docker run -it --rm -v $PWD:$PWD chargefw2 --mode charges \
    --input-file $PWD/doc/molecules.sdf --chg-out-dir $PWD/
```

However, it is good practise to only bind-mount necessary directories. In the following example we bind-mount two directories: one directory that contains our input files and the other directory for ChargeFW2 output.

```bash
$ INPUT_DIRECTORY="/path/to/input/directory/"
$ OUTPUT_DIRECTORY="/path/to/output/directory/"

$ docker run -it --rm \
    -v $INPUT_DIRECTORY:$INPUT_DIRECTORY \
    -v $OUTPUT_DIRECTORY:$OUTPUT_DIRECTORY \
    chargefw2 --mode charges \
    --input-file $INPUT_DIRECTORY/path/to/file \
    --chg-out-dir $OUTPUT_DIRECTORY
```

The previous examples used a container which was destroyed once ChargeFW2 finished running. If you don't want to keep recreating the container for each use of ChargeFW2, you can create a detached container that will keep running in the background. You can then use `docker exec` to run commands inside the container.

```bash
# Create detached container
$ CONTAINER_ID=$(docker run -dt --rm --entrypoint bash -v $PWD:$PWD chargefw2)

# Run ChargeFW2 in detached container (can run multiple commands)
$ docker exec -it $CONTAINER_ID chargefw2 --mode charges \
    --input-file $PWD/doc/molecules.sdf --chg-out-dir $PWD/

# Remove container
$ docker rm -f $CONTAINER_ID
```

## Usage

The [documentation](doc/documentation.md) for the application and its [Python bindings](doc/ChargeFW2%20-%20tutorial.pdf) is located in the [doc](doc) folder.

## How to cite
If you found ChargeFW2 or Atomic Charge Calculator II helpful, please cite: [Raček, T., Schindler, O., Toušek, D., Horský, V., Berka, K., Koča, J., & Svobodová, R. (2020). Atomic Charge Calculator II: web-based tool for the calculation of partial atomic charges. Nucleic Acids Research](https://doi.org/10.1093/nar/gkaa367).
