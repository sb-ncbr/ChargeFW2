# ChargeFW2

Application for computing partial atomic charges using selected empirical methods.

## Compilation requirements
- [CMake](https://cmake.org/) >= 3.9
- [GCC](https://gcc.gnu.org/) >= 8.2 or [Clang](https://clang.llvm.org/) >= 8
- [Boost](https://www.boost.org/) >= 1.69
- [Eigen](http://eigen.tuxfamily.org) == 3.3
- [fmt](https://fmt.dev) >= 5.3.0
- [nanoflann](https://github.com/jlblancoc/nanoflann) == 1.3.0
- [NLopt](https://github.com/stevengj/nlopt) >= 2.4.2
- [JSON for Modern C++](https://github.com/nlohmann/json) >= 3.6.1
- [GEMMI](https://github.com/project-gemmi/gemmi) >= 0.3.1

## Installation
After downloading and unpacking the sources, run the following in the ChargeFW2 directory:

```shell script
$ mkdir build
$ cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=<WHERE-TO-INSTALL>
$ make
$ make install
```
