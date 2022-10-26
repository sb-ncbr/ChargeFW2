# ChargeFW2

Application for computing partial atomic charges using selected empirical methods.
ChargeFW2 is the computational core of [Atomic Charge Calculator II](https://acc2.ncbr.muni.cz).

See the [short description](https://acc2.ncbr.muni.cz/static/methods.pdf) of implemented methods. 

## Compilation requirements
- [CMake](https://cmake.org/) 3.17
- [GCC](https://gcc.gnu.org/) 10 or [Clang](https://clang.llvm.org/) 10
- [Boost](https://www.boost.org/) 1.69
- [Eigen](http://eigen.tuxfamily.org) 3.3
- [fmt](https://fmt.dev) 6.2.1
- [nanoflann](https://github.com/jlblancoc/nanoflann) 1.4.3
- [JSON for Modern C++](https://github.com/nlohmann/json) 3.7.3
- [GEMMI](https://github.com/project-gemmi/gemmi) 0.4.7
- [pybind11](https://github.com/pybind/pybind11) 2.5.0

Tested on Fedora 32-36 and Ubuntu 20.04-22.04. Other version of the libraries might work too however this was not tested.

## Installation
After downloading and unpacking the sources, run the following in the ChargeFW2 directory:

```shell script
$ mkdir build
$ cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=<WHERE-TO-INSTALL> -DCMAKE_BUILD_TYPE=Release
$ make
$ make install
```

## Docker
Rather than installing all dependencies, you can run ChargeFW2 directly in a docker container:
```shell script
$ cd docker  
$ docker build -t chargefw2 .
$ docker run -it chargefw2
```

## Usage

The [documentation](doc/documentation.md) for the application and its [Python bindings](doc/ChargeFW2%20-%20tutorial.pdf) is located in the [doc](doc) folder.

## How to cite
If you found ChargeFW2 or Atomic Charge Calculator II helpful, please cite: [Raček, T., Schindler, O., Toušek, D., Horský, V., Berka, K., Koča, J., & Svobodová, R. (2020). Atomic Charge Calculator II: web-based tool for the calculation of partial atomic charges. Nucleic Acids Research](https://doi.org/10.1093/nar/gkaa367).
