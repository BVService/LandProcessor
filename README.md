# Overview

LandProcessor is a GIS tool for representing landcapes into connected homogeneous spatial units.  
It is developped for the BVservice project and is mainly focused on surface hydrology in agricultural landscapes. Resulting connected spatial units are distributed into classes : SU (Surface Units such as plots), LI (Linear Interfaces such as hedges or benches), RS (Reach Segments such as ditches or rivers).  
This code source is mainly a shared library. it also provides simulators for the [OpenFLUID platform](http://www.openfluid-project.org/)


# Building from sources

## Requirements

LandProcessor is written in C++11 and is intented to run on Linux platforms.
The following tools and frameworks are also required:
- GCC 4.9 or higher
- CMake 2.8.12 or higher
- OpenFLUID 2.1.3 or higher
- GDAL 1.11
- GRASS 7.xx
- Qt 4.8.x


## Building and testing

To build from sources, execute the following commands in the source tree using a terminal
```
mkdir _build
cd _build
cmake ..
make
```
To run the provided tests, execute the following commands from the `_build` directory using a terminal:
```
mkdir _build
cd _build
ctest
make
```

# TODO

- [ ] Initialize mp_VectorDriver in constructor once for all
- [ ] Remove absolute paths for temporary files, use openfluid::base::Environment::getTempDir() instead
- [ ] Refactor main source code to be more easily maintainable
- [ ] Introduce correct land use codes in tests
- [ ] Improve performance by moving invariants out of loops
