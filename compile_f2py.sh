#!/bin/bash

source /opt/intel/oneapi/setvars.sh
#module load compiler/intel/
#module load numlib/mkl
source ~/anaconda3/bin/activate

f2py3 --fcompiler=gfortran -m eta_gamma_pade -c eta_gamma_pade.f90 --opt="-Ofast"
f2py3 --fcompiler=gfortran -m sparsity -c sparsity.f90 --f90flags="-c -g -DF2PY_REPORT_ON_ARRAY_COPY -O0 -fcheck=all -fbounds-check -Wall -fbacktrace -ffree-line-length-512" --opt="-Ofast"
f2py3 --fcompiler=gfortran -m sparsity_corrfunc -c sparsity_corrfunc.f90 --f90flags="-ffree-line-length-512" --opt="-Ofast"
f2py3 --fcompiler=gfortran -m sparse_friction -c sparse_friction.f90 --f90flags="-fopenmp -ffree-line-length-512 -I/opt/intel/oneapi/mkl/latest/include" --opt="-O3" -L/opt/intel/oneapi/mkl/latest/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl  
f2py3 --fcompiler=gfortran -m sparse_corrfunc -c sparse_corrfunc.f90 --f90flags="-c -g -DF2PY_REPORT_ON_ARRAY_COPY -O0 -fcheck=all -fbounds-check -Wall -fbacktrace -fopenmp -ffree-line-length-512 -I/opt/intel/oneapi/mkl/latest/include" --opt="-O3" -L/opt/intel/oneapi/mkl/latest/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl  

#f2py3 --fcompiler=intelem -m eta_gamma_pade -c eta_gamma_pade.f90 --opt="-fast"
#f2py3 --fcompiler=intelem -m sparsity -c sparsity.f90 --opt="-fast"
#f2py3 --fcompiler=intelem -m sparsity_corrfunc -c sparsity_corrfunc.f90 --opt="-fast"
#f2py3 --fcompiler=intelem -m sparse_friction -c sparse_friction.f90 --f90flags="-qopenmp -qmkl" --opt="-fast"
#f2py3 --fcompiler=intelem -m sparse_corrfunc -c sparse_corrfunc.f90 --f90flags="-qopenmp -qmkl" --opt="-fast"

