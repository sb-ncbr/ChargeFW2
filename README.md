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
- [nanoflann](https://github.com/jlblancoc/nanoflann) 1.3.0
- [JSON for Modern C++](https://github.com/nlohmann/json) 3.7.3
- [GEMMI](https://github.com/project-gemmi/gemmi) 0.3.6
- [pybind11](https://github.com/pybind/pybind11) 2.5.0

Tested on Fedora 32-34 and Ubuntu 20.04. Other version of the libraries might work too however this was not tested.

## Installation
After downloading and unpacking the sources, run the following in the ChargeFW2 directory:

```shell script
$ mkdir build
$ cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=<WHERE-TO-INSTALL> -DCMAKE_BUILD_TYPE=Release
$ make
$ make install
```

## Usage

The [documentation](doc/documentation.md) for the application and its [Python bindings](doc/ChargeFW2%20-%20tutorial.pdf) is located in the [doc](doc) folder.
