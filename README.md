# linequad

Nearly singular quadrature for line integrals in 2D and 3D, using the "singularity swap"
method.

This repository accompanies the paper "Accurate quadrature of nearly singular line
integrals in two and three dimensions" by L. af Klinteberg and A. Barnett. It contains
complete Matlab implementations of the algorithms for 2D and 3D line quadrature. Provided
are also C implementations for all computational steps of the 3D line quadrature. These are linked to Matlab via MEX, but can also be compiled into a library.

The 3D code has the most mature interface, in the form of `line3_near_weights`.

## Getting started

To initialize the required submodules, run `git submodule update --init --recursive`.

In the `matlab` directory, open Matlab and run:
```Matlab
make % To build the MEX code (gcc required)

init % To setup the paths for running things (must be done every session)
```

Then run one of the examples (while still in the `matlab` directory):

* The programs in `matlab/examples` demo the basics of the 3D line quadrature.
* The programs in `matlab/examples` can be used to reproduce all the results (2D and 3D) shown in the paper.
