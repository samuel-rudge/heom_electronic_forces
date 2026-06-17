# heom_electronic_forces
Calculate Markovian and non-Markovian nonadiabatic electronic forces acting on molecules interacting with metal surfaces via the Hieararchical Equations of Motion (HEOM) approach.

# This Python file contains the main code used to perform HEOM time-propagation. Because implementing
# HEOM is a complex task, various parts of the implementation are split into Python modules and Fortran
# subroutines. This code imports all modules and runs them in the correct order. It also contains the code 
# that determines what we do with the end result; in this case we plot the current and elements of the 
# density matrix.
#
# USAGE - RUN FROM COMMAND LINE (TERMINAL) WITH ANACONDA:
#       Note that one must create the Python wrappers from the Fortran subroutines first by running ./compile_f2py.sh
#       python3 friction_heom_main.py ss
#       (for steady state forces)
#       python3 friction_heom_main.py markovian
#       (for time-dependent Markovian/Nonmarkovian forces)
#

#       Alternatively, one could go into a python environment (type Python into command line) and run each line manually.
