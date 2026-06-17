# heom_electronic_forces
Calculate Markovian and non-Markovian nonadiabatic electronic forces acting on molecules interacting with metal surfaces via the Hieararchical Equations of Motion (HEOM) approach.

This Python file contains the main code used to perform HEOM time-propagation, which is then used to calculate electronic forces. Because implementing HEOM is a complex task, various parts of the implementation are split into Python modules and Fortran subroutines. This code imports all modules and runs them in the correct order. 

USAGE - RUN FROM COMMAND LINE (TERMINAL) WITH ANACONDA:
      First, one must create the Python wrappers from the Fortran subroutines by running ./compile_f2py.sh. This requires an anaconda installation. Then, one runs:
      
      python3 friction_heom_main.py ss
      (for steady state forces)
      python3 friction_heom_main.py markovian
      (for time-dependent Markovian/Nonmarkovian forces)
