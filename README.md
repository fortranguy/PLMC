# PLMC

A simple Physics of Liquids Monte Carlo simulation code written in Fortran.

## Getting Started

### Prerequisites

You need:
* a Fortran compiler such as [gfortran 6.1+](https://gcc.gnu.org/wiki/GFortran) or
[ifort 16.0+](https://software.intel.com/en-us/qualify-for-free-software)
* [CMake 3.0+](https://cmake.org/) to build the project
* [json-fortran 5.0+](https://github.com/jacobwilliams/json-fortran) to input simulation parameters.

You may need:
* [FORD 5.0+](https://github.com/cmacmackin/ford) to generate the
[documentation](https://fortranguy.github.io/PLMC/)
* [LuaLaTex 0.95+](https://www.tug.org/texlive/) to draw the
[class diagrams](https://github.com/fortranguy/PLMC/blob/gh-pages/plmc_design.pdf)
* [Julia 0.5+](http://julialang.org/) to create an initial configuration or check
energy consistency
* [Python 3+](https://www.continuum.io/downloads) to analyse dipoles clusters.

### Installation

Make a folder and execute cmake from it:
```bash
mkdir build
cd build
FC=gfortran cmake /path/to/PLMC
```
The default build is debug. You can change it using `make edit_cache`.
From now on `PLMC` will be a shorthand for `/path/to/PLMC`.
Now you can build the project:
```bash
make -j <num_jobs>
```
Binaries should be in `bin/Debug/`.

### Run a simulation

Copy a parameters file to generate a Markov chain:
```bash
cp PLMC/data/hard_spheres.json .
```
Create an initial configuration:
```bash
julia PLMC/scripts/randomCoordinates.jl hard_spheres.json
```
If you want to use `randomCoordinates.jl`, you may need to add this line in your `.juliarc.jl`:
```julia
push!(LOAD_PATH, "PLMC/scripts/")
```
Finally, you can run your simulation:
```bash
plmc_generate hard_spheres.json
```
You can plot translations acceptance ratio to check the progress:
```
gnuplot
plot "generating_box_1/changes_success_1.out"
```

### Analyse the results

Copy another parameters file to explore the generated configurations:
```bash
cp PLMC/data/exploring.json .
```
If you want the
[radial distribution function](https://en.wikipedia.org/wiki/Radial_distribution_function), enter:
```bash
radial hard_spheres.json exploring.json generating_box_1/coordinates/coordinates_*.xyz
```
And you can plot it:
```
gnuplot
plot "radial_1-1.out" w l
```

## Legacy

I hope to defend my thesis soon.
Feel free to fork this project, make it your own and share it.

## Licence

This code is under GPL3 licence.
