# Build Instructions

This repository includes performance-critical Fortran extensions that implement the core numerical routines of the HEOM propagation, including sparse matrix operations, friction kernel evaluation, and correlation function evaluation.

These routines are interfaced to Python using f2py and must be compiled before running any simulations.

---

# Compilation Requirements

The code depends on:

- a Fortran compiler (either gfortran or Intel Fortran)
- BLAS/LAPACK support via Intel MKL
- OpenMP support for parallel execution over the nuclear coordinate grid
- a working f2py3 installation (provided via NumPy in most Python environments)

---

# Recommended Environment

The code has been most extensively tested with:

- Anaconda 2022 distribution
- Python 3.x within that environment
- Intel MKL installed via Anaconda or system modules
- gfortran or Intel compiler toolchain

Important: compatibility with newer Anaconda versions is not guaranteed.

---

# Compiler Options

Two compilation paths are supported:

## 1. gfortran (recommended for local machines)

This is the default tested setup for local development. It uses gfortran with MKL linked explicitly.

## 2. Intel compiler (recommended for HPC clusters)

On clusters with Intel CPUs, using the Intel compiler (ifort / ifx toolchain) typically yields better performance.

However, note:

- recent Intel compilers may break compatibility with f2py
- therefore gfortran remains the most reliable option across systems

---

# MKL Dependency

The Fortran extensions require Intel MKL for optimized sparse linear algebra.

Users must ensure that:

- MKL is installed (via Anaconda or Intel oneAPI)
- include and library paths are correctly set in the compilation script

In particular, the paths in compile_f2py.sh may need to be adjusted to match the local system installation.

---

# Compilation Procedure

Before running any simulations, the user must execute the provided build script:

run compile_f2py.sh from the main directory.

This script compiles the following Fortran modules:

- eta_gamma_pade
- sparsity
- sparsity_corrfunc
- sparse_friction
- sparse_corrfunc

These modules implement:

- HEOM bath decomposition routines
- sparse Liouvillian construction
- friction tensor evaluation
- force–force correlation function evaluation
- OpenMP-parallelized coordinate grid operations

---

To run the script, the user should execute it in a terminal using bash (for example: ./compile_f2py.sh).

# Important Notes on the Build Script

The provided compile_f2py.sh script:

- assumes Intel oneAPI environment variables are available
- optionally activates an Anaconda environment
- links explicitly against MKL libraries
- uses OpenMP for parallel sections in selected modules

Users must ensure that:

- the Intel MKL path in the script matches their installation
- the correct compiler flags are used for their system
- only one compiler backend (gfortran OR Intel) is used consistently

The script includes commented Intel compiler versions for HPC systems; these may be enabled manually if supported.

---

# Common Issues

## f2py compatibility

f2py support is known to be fragile across compiler versions. If compilation fails:

- prefer gfortran over Intel oneAPI compilers
- ensure NumPy and Python versions are consistent with Anaconda 2022

---

## MKL linking errors

If linking fails:

- verify MKL installation path
- check that environment variables (LD_LIBRARY_PATH, MKLROOT) are correctly set
- adjust compile_f2py.sh accordingly

---

## OpenMP issues

If parallel execution fails or is unstable:

- ensure -fopenmp (or Intel equivalent) is enabled
- verify thread settings on your system

---

# Summary

The recommended and most stable setup is:

- Anaconda 2022 environment
- gfortran compiler
- Intel MKL (via Anaconda or oneAPI)
- compile via compile_f2py.sh before running any Python scripts
