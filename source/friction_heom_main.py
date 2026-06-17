"""
HEOM-Based Electronic Friction Framework
Main driver: friction_heom_main.py

This script is the primary entry point for all simulations in this repository.

It implements a Hierarchical Equations of Motion (HEOM) based framework for
computing electronic forces, currents, and friction kernels in vibronic
quantum transport systems with classical nuclear degrees of freedom and
quantum electronic open-system dynamics.

-----------------------------------------------------------------------
USAGE MODES
-----------------------------------------------------------------------

The script is executed from the terminal in the main directory.

Two main simulation modes are available:

1. Steady-state electronic observables:
   python3 friction_heom_main.py ss

   Computes steady-state quantities on a predefined nuclear coordinate grid,
   including:
   - adiabatic mean forces
   - electronic currents (adiabatic and non-adiabatic contributions)

   Outputs are written as .dat files in the working directory.

2. Markovian friction and correlation functions (Can only be run after ss):
   python3 friction_heom_main.py markovian

   Computes:
   - Markovian electronic friction tensor
   - force–force correlation functions

   Optional non-Markovian integrand evaluation is enabled via the flag
   print_integrand_yn in input_parameters.py.

-----------------------------------------------------------------------
OUTPUT
-----------------------------------------------------------------------

All results are written directly to the working directory. This includes:

- steady-state observables (forces, currents)
- friction tensors and correlation functions
- optional intermediate and diagnostic files
- simulation metadata (simulation_info.dat)

Existing output files may be overwritten upon execution.

-----------------------------------------------------------------------
IMPLEMENTATION DETAILS
-----------------------------------------------------------------------

- Electronic dynamics is treated using HEOM with user-defined hierarchy depth.
- Nuclear degrees of freedom are treated classically on a coordinate grid.
- Heavy numerical routines are implemented in Fortran and interfaced via f2py.
- Sparse linear algebra and HEOM propagation use MKL-accelerated routines.
- OpenMP parallelization is used over nuclear coordinate evaluations in
  Markovian mode.

-----------------------------------------------------------------------
DEPENDENCIES

- Python (Anaconda 2022 recommended)
- NumPy / SciPy
- Fortran compiler (gfortran or Intel)
- Intel MKL (required for optimized sparse routines)
- Precompiled Fortran extensions (see compile_f2py.sh)

-----------------------------------------------------------------------
BUILD STEP (REQUIRED BEFORE RUNNING)

Before executing this script, compile all Fortran extensions by running:
compile_f2py.sh in the main directory using a bash shell (e.g. bash compile_f2py.sh).

The build links against Intel MKL and uses f2py to generate Python interfaces
to the Fortran modules.

-----------------------------------------------------------------------
IMPORTANT

- All output files are written to the working directory and may be overwritten.
- The Markovian mode requires a predefined nuclear coordinate grid.
- For non-Markovian single-point evaluation, the grid must be reduced to one
  point by setting x_max_total = x_min_total + dx_grid in input_parameters.py.
- The code has been primarily tested with Anaconda 2022 environments.
  Compatibility with newer distributions is not guaranteed.
"""

#### IMPORT PYTHON MODULES ####

from re import S
from shutil import make_archive
from tokenize import Double
import numpy as np                                                                                          # Import intrinsic Python modules
import scipy as sc
from time import perf_counter, asctime
import pickle
import sys

#### IMPORT HEOM MODULES ####

from constants import *                                                                                     # Import HEOM modules
import system
from input_parameters import *

#ss_or_markovian_key = "markovian"
ss_or_markovian_key = sys.argv[1]

if ss_or_markovian_key == "ss":
    print("Running code for the steady state parts")
    redo_sparsity = True
    redo_adiabatic_results = True
    redo_ss_spatial_derivative = True
    redo_friction = False
    redo_corrfunc = False
    redo_simulation_info = True
elif ss_or_markovian_key == "markovian":
    print("Running code for the Markovian integration parts")
    redo_sparsity = True
    redo_adiabatic_results = False
    redo_ss_spatial_derivative = False
    redo_friction = True
    redo_corrfunc = True
    redo_simulation_info = False

nondiag_key = False
if small_polaron_yn:
    if ss_or_markovian_key == "markovian":
        nondiag_key = True

sys.stdout.flush()
import Index_pm_filter
import generating_sparsity_class
import gmres_ss_solver
import gmres_ss_spatial_derivative
import markovian_friction
import markovian_corrfunc
import eta_gamma_barycentric
import psutil
from matplotlib import pyplot as plt
import gc
import os
import sys
import scipy.sparse as sparse
import scipy.sparse.linalg as sparse_linalg

gc.collect()

# if (x_min > x_max_total):
#    output_info = open( 'simulation_info.txt', 'w')                                                             # Open txt file called simulation_info.txt
#    output_info.write( 'Simulation unnecessary')
#    sys.exit("Simulation unnecessary")

t_start_program = asctime()                     
if os.path.isfile('simulation_info.dat') and not redo_simulation_info:
    output_info = open('simulation_info.dat','a')
else:
    output_info = open('simulation_info.dat','w') 
    output_info.write('The simulation starts at '+str(t_start_program)+'\n')

# ---------------------------------------------------------------
#            DEFINE OPERATORS IN MOLECULAR HILBERT SPACE
# ---------------------------------------------------------------

system_output = system.system_operators(Single_El_Int,Double_El_Int,Vib_Freq_qu,El_Nuclear_Couplings_qu,
                                El_Nuclear_Couplings_cl,Nel,N_qu_vib_modes,N_cl_vib_modes,max_occ_qu_vib_modes,dim_rho,
                                len_x_vec,x_vec,dx_sssd,small_polaron_yn,nondiag_key)

d_ops,d,ddag,Fock_states,Ham,Ham_x,deriv_Ham = system_output[0:7]
if bool(N_qu_vib_modes):
    b_ops,b,bdag,nvib_qu,nvib_qu_diag,nvib_qu_nondiag = system_output[7:13]
    nvib_qu_log = np.array(nvib_qu,dtype=bool)
    nvib_qu_diag_log = np.array(nvib_qu_diag,dtype=bool)
    nvib_qu_nondiag_log = np.array(nvib_qu_nondiag,dtype=bool)
    d_ops_dressed,FC_Matrix,FC_Matrix_Fock_Space,Kmatrix,Kmatrix_inv = system_output[13:19]

d_ops_comp = d_ops
if bool(small_polaron_yn):
    if ss_or_markovian_key == "ss":
        d_ops_comp = d_ops_dressed

### Define logical versions of all system operators ###

Ham_log = np.logical_or(np.array(Ham_x[:,:,0,0],dtype=bool),np.eye(dim_rho,dtype=bool))
deriv_Ham_log = np.array(deriv_Ham,dtype=bool)
d_ops_log = np.array(d_ops,dtype=bool)
d_ops_comp_log = np.array(d_ops_comp,dtype=bool)
rho_0_log = np.array(rho_0,dtype=bool)
identity_dim_rho = np.eye(dim_rho)

# ---------------------------------------------------------------
#      BATH-CORRELATION EXPANSION - BARYCENTRIC AND PADE
# ---------------------------------------------------------------

EtaGamma = eta_gamma_barycentric.bath_correlation_decomposition(Ncutoff,specwidth,Nsupport_points_barycentric,
                Npoles_pade,symmetrized_fermi_specwidth,Temp,Nleads,Nsign,muvec,tol_Gamma_barycentric,
                tol_fermi_symmetrized_barycentric,wbl_YN,analytic_spectral_function_decomposition,tol_F)
eta_vec_barycentric,gamma_vec_barycentric = EtaGamma.barycentric_bath_correlation_expansion()
eta_vec_pade,gamma_vec_pade = EtaGamma.pade_bath_correlation_expansion()

if pole_choice == "pade":
    eta_vec = eta_vec_pade
    gamma_vec = gamma_vec_pade
elif pole_choice == "barycentric":
    eta_vec = eta_vec_barycentric
    gamma_vec = gamma_vec_barycentric
elif pole_choice == "prony":
    raise ValueError("Prony/MPM decomposition not yet implemented")
    # eta_vec = eta_vec_prony
    # gamma_vec = gamma_vec_prony
else:
    raise ValueError("Choose an appropriate pole decomposition scheme: Options are pade or barycentric")

if wbl_YN == 0:                                                                    # If not under the wide-band limit, they are assumed to have a Lorentzian density of states
    Npoles = len(eta_vec[0,0,:]) - 1
    Nmodes = (Npoles+1)*Nel*Nleads*Nsign                                           # Calculate number of modes outside of wide-band limit
elif wbl_YN == 1:
    Npoles = len(eta_vec[0,0,:])
    Nmodes = Npoles*Nel*Nleads*Nsign         

sys.stdout.flush()
t_start = perf_counter()                                                                               # Return value of performance counter (internal clock with no reference time)

# ---------------------------------------------------------------
#                 INDEX GENERATION OF ADOs IN HEOM
# ---------------------------------------------------------------

Indices = Index_pm_filter.Hierarchy_index(Nmax,Nel,Npoles,Nleads,Nsign,Nmodes,wbl_YN)                              # Define object of Hierarchy_index class with HEOM parameters as input
if wbl_YN == 0:
    if filtering_YN == 1:                                                                                   # Run filtering process if filtering_YN is true
        max_V_km = np.max(el_lead_couplings)                                                                        # Define maximum coupling between leads and electronic levels in the system
        KsigLm_filtered,Un_Ind_filtered,Hier_ind_filtered,Index_minus_filtered,Index_plus_filtered = Indices.Print_Filtered_Ind_Info(tol,eta_vec,gamma_vec,max_V_km)
    else:
        KsigLm,Un_Ind,Hier_ind,Index_Minus,Index_Plus,len_un_ind,len_index_plus,tier_index = Indices.Print_Ind_Info()
                                                                                                            # Return index information from Indices object; see Index_pm_filter for details
elif wbl_YN == 1:
    if filtering_YN == 1:                                                                                   # Run filtering process if filtering_YN is true
        max_V_km = np.max(el_lead_couplings)                                                                        # Define maximum coupling between leads and electronic levels in the system
        KsigLm_filtered,Ksig0m_filtered,Un_Ind_filtered,Hier_ind_filtered,Index_minus_filtered,Index_plus_filtered = Indices.Print_Filtered_Ind_Info(tol,eta_vec,gamma_vec,max_V_km)
    else:
        KsigLm,Ksig0m,Un_Ind,Hier_ind,Index_Minus,Index_Plus,len_un_ind,len_index_plus,tier_index = Indices.Print_Ind_Info()
                                                                                                            # Return index information from Indices object; see Index_pm_filter for details

sys.stdout.flush()
t_end = perf_counter()                                                                                      # Return value of performance counter at end of index generation
output_info.write("Elapsed time of indices generation: " + str(t_end-t_start) +'\n')                        # Write into simulation_info.txt the time taken to perform index generation
output_info.close()
output_info = open('simulation_info.dat','a')

# ---------------------------------------------------------------
#          TRANSFORMING HEOM TO SPARSE REPRESENTATION 
# ---------------------------------------------------------------

sparsity_key = True
if not redo_everything:
    if os.path.isfile("sparsity_ingredients.p"):
        if not bool(redo_sparsity):
            sparsity_key = False

if sparsity_key:
    sparsity_object = generating_sparsity_class.sparsity_heom_liouvillian(ksiglm=KsigLm,tier_index=tier_index,
                            index_minus=Index_Minus,index_plus=Index_Plus,d_ops_comp=d_ops_comp,ham=Ham_x,
                            d_ops_comp_log=d_ops_comp_log,ham_log=Ham_log,rho_0_log=rho_0_log,max_expan_order=max_expan_order,
                            dim_rho=dim_rho,len_index_plus=len_index_plus,len_un_ind=len_un_ind,nmax=Nmax,nel=Nel,
                            wbl_yn=wbl_YN,degenerate_levels=degenerate_levels,atol=atol,rtol=rtol,un_ind=Un_Ind,
                            gamma_vec=gamma_vec,eta_vec=eta_vec,nsign=Nsign,nleads=Nleads,npoles=Npoles,ham_x=Ham_x,
                            len_x_vec=len_x_vec,el_lead_couplings=el_lead_couplings)
    if nondiag_key:
        sparsity_ingredients_nondiag = sparsity_object.return_sparse_heom()
        pickle.dump(sparsity_ingredients_nondiag,open("sparsity_ingredients_nondiag.p","wb"))
        sparsity_ingredients_file = open("sparsity_ingredients.p","rb")
        sparsity_ingredients = pickle.load(sparsity_ingredients_file)
        sparsity_ingredients_file.close()    
    else:
        sparsity_ingredients = sparsity_object.return_sparse_heom()
        pickle.dump(sparsity_ingredients,open("sparsity_ingredients.p","wb"))
    print("Sparsity has been (re)generated")
else:
    sparsity_ingredients_file = open("sparsity_ingredients.p","rb")
    sparsity_ingredients = pickle.load(sparsity_ingredients_file)
    sparsity_ingredients_file.close()
    if nondiag_key:
        sparsity_ingredients_nondiag_file = open("sparsity_ingredients_nondiag.p","rb")
        sparsity_ingredients_nondiag = pickle.load(sparsity_ingredients_nondiag_file)
        sparsity_ingredients_nondiag_file.close()
    print("Sparsity has been loaded and not regenerated")

sys.stdout.flush()

### Generate sparse information on Liouvillian of HEOM ###
pair_info_row_fil,pair_info_col_fil,pair_values_fil,npairs_fil,npairs_uf,\
nnz_elements_sparse_fil,nnz_elements_sparse_zeroth_tier_fil,row_old_indices,\
atol_vec,rtol_vec,rho_nonzeros_sparse,isreal_sparse,complex_coefficients,\
sparse_trace_array,nnz_elements_zeroth_tier,trace_cols,rhs_vector,\
rho_nonzeros,rho_sparsity,nnz_elements,is_connected_array,rho_out = sparsity_ingredients
if nondiag_key:
    pair_info_row_nondiag_fil,pair_info_col_nondiag_fil,pair_values_nondiag_fil,npairs_nondiag_fil,npairs_nondiag_uf,\
    nnz_elements_sparse_nondiag_fil,nnz_elements_sparse_zeroth_tier_nondiag_fil,row_old_indices_nondiag,\
    atol_vec_nondiag,rtol_vec_nondiag,rho_nonzeros_sparse_nondiag,isreal_sparse_nondiag,complex_coefficients_nondiag,\
    sparse_trace_array_nondiag,nnz_elements_zeroth_tier_nondiag,trace_cols_nondiag,rhs_vector_nondiag,\
    rho_nonzeros_nondiag,rho_sparsity_nondiag,nnz_elements_nondiag,is_connected_array_nondiag,rho_out_nondiag = \
        sparsity_ingredients_nondiag

t_start = t_end
t_end = perf_counter()
output_info.write("Elapsed time to transform to sparse representation: " + str(t_end-t_start) +'\n')
output_info.close()
output_info = open('simulation_info.dat','a')

# --------------------------------------------------------------------------
#               STEADY STATE VIA DIRECT SOLVER WITH GMRES  
# --------------------------------------------------------------------------

adiabatic_ss_key = True
if not redo_everything:
    if not bool(redo_adiabatic_results):
        if os.path.isfile("adiabatic_ss_results.p"):
            adiabatic_ss_key = False

if adiabatic_ss_key:
    if bool(N_qu_vib_modes):
        nvib_object = [nvib_qu_diag,nvib_qu_nondiag,nvib_qu_diag_log,nvib_qu_nondiag_log,dim_vib_mode_qu,\
                        Kmatrix,Kmatrix_inv]
    else:
        nvib_object = []
    adiabatic_ss_object = gmres_ss_solver.steady_state_x_grid(pair_info_col=pair_info_col_fil,
        pair_info_row=pair_info_row_fil,pair_values=pair_values_fil,npairs=npairs_fil,
        nnz_elements_sparse=nnz_elements_sparse_fil,sparse_trace_array=sparse_trace_array,
        rhs_vector=rhs_vector,nnz_elements_sparse_zeroth_tier=nnz_elements_sparse_zeroth_tier_fil,
        complex_coefficients=complex_coefficients,d_ops=d_ops_comp,isreal_sparse=isreal_sparse,
        deriv_ham=deriv_Ham,deriv_ham_log=deriv_Ham_log,rho_nonzeros_sparse=rho_nonzeros_sparse,
        tier_index=tier_index,un_ind=Un_Ind,ksiglm=KsigLm,nvib_object=nvib_object,ham_x=Ham_x,
        trace_cols=trace_cols,fock_states=Fock_states)
    adiabatic_ss_results = adiabatic_ss_object.return_ss_all_x_values()
    pickle.dump(adiabatic_ss_results,open("adiabatic_ss_results.p","wb"))
    print("Adiabatic steady state has been (re)generated")
else:
    adiabatic_ss_results_file = open("adiabatic_ss_results.p","rb")
    adiabatic_ss_results = pickle.load(adiabatic_ss_results_file)
    adiabatic_ss_results_file.close()
    if nondiag_key:
        rho_ss_x_arr = adiabatic_ss_results[0]
        ss_nondiag_object = \
            gmres_ss_solver.\
                steady_state_nondiag(rho_ss_x_arr=rho_ss_x_arr,Kmatrix=Kmatrix,Kmatrix_inv=Kmatrix_inv,
                    is_connected_array_nondiag=is_connected_array_nondiag,len_un_ind=len_un_ind,
                    rho_nonzeros_nondiag=rho_nonzeros_nondiag,complex_coefficients=complex_coefficients,
                    nnz_elements_sparse=nnz_elements_sparse_fil,rho_nonzeros_sparse=rho_nonzeros_sparse,
                    nnz_elements_sparse_nondiag_fil=nnz_elements_sparse_nondiag_fil,
                    nnz_elements_nondiag=nnz_elements_nondiag)
        rho_ss_x_arr_nondiag = ss_nondiag_object.return_ss_x_arr_nondiag()
        adiabatic_ss_results.append(rho_ss_x_arr_nondiag)
        pickle.dump(adiabatic_ss_results,open("adiabatic_ss_results.p","wb"))
    print("Adiabatic steady state has been loaded and not regenerated")

sys.stdout.flush()
if bool(N_qu_vib_modes):
    rho_ss_x_arr,rho_system,average_electronic_force_mol_vec,average_electronic_force_molleads_vec,current_vec,\
            average_vib_excitation_vec_diag,average_vib_excitation_vec_nondiag = adiabatic_ss_results[0:7]
    if nondiag_key:
        rho_ss_x_arr_nondiag = adiabatic_ss_results[7]
else:
    rho_ss_x_arr,rho_system,average_electronic_force_mol_vec,\
        average_electronic_force_molleads_vec,current_vec = adiabatic_ss_results

t_start = t_end                                                                                             # Define start time of time-propagation code
t_end = perf_counter()                                                                                      # Define end time of time-propagation code
output_info.write("Elapsed time after calculating adiabatic quantities via GMRES: " + str(t_end-t_start) +'\n')                               # Write the total elapsed computer time for time-propagation code to run
output_info.close()
output_info = open('simulation_info.dat','a')

# --------------------------------------------------------------------------
#                      SPATIAL DERIVATIVE WITH GMRES
# --------------------------------------------------------------------------

ss_spatial_derivative_key = True
if not redo_everything:
    if not bool(redo_ss_spatial_derivative):
        if os.path.isfile("rho_ss_spatial_derivative.npy"):
            ss_spatial_derivative_key = False

if ss_spatial_derivative_key:
    ss_spatial_derivative_object = gmres_ss_spatial_derivative.steady_state_spatial_derivative(pair_info_col=pair_info_col_fil,
        pair_info_row=pair_info_row_fil,pair_values=pair_values_fil,npairs=npairs_fil,
        nnz_elements_sparse=nnz_elements_sparse_fil,
        sparse_trace_array=sparse_trace_array,rhs_vector=rhs_vector,
        nnz_elements_sparse_zeroth_tier=nnz_elements_sparse_zeroth_tier_fil,
        complex_coefficients=complex_coefficients,isreal_sparse=isreal_sparse,trace_cols=trace_cols,
        rho_ss_x_arr=rho_ss_x_arr)
    rho_ss_spatial_derivative = ss_spatial_derivative_object.return_ss_spatial_derivative()
    np.save("rho_ss_spatial_derivative.npy",rho_ss_spatial_derivative)
else:
    rho_ss_spatial_derivative = np.load("rho_ss_spatial_derivative.npy")
    if nondiag_key:
        ss_spatial_derivative_nondiag_object = \
            gmres_ss_spatial_derivative.\
                steady_state_spatial_derivative_nondiag(rho_ss_spatial_derivative=rho_ss_spatial_derivative,\
                    Kmatrix=Kmatrix,Kmatrix_inv=Kmatrix_inv,
                    is_connected_array_nondiag=is_connected_array_nondiag,len_un_ind=len_un_ind,
                    rho_nonzeros_nondiag=rho_nonzeros_nondiag,complex_coefficients=complex_coefficients,
                    nnz_elements_sparse=nnz_elements_sparse_fil,rho_nonzeros_sparse=rho_nonzeros_sparse,
                    nnz_elements_sparse_nondiag_fil=nnz_elements_sparse_nondiag_fil,
                    nnz_elements_nondiag=nnz_elements_nondiag)
        rho_ss_spatial_derivative_nondiag = ss_spatial_derivative_nondiag_object.return_ss_spatial_derivative_nondiag()
        np.save("rho_ss_spatial_derivative_nondiag.npy",rho_ss_spatial_derivative_nondiag)
    print("Spatial derivative of steady state has been loaded and not regenerated")

sys.stdout.flush()
t_start = t_end                                                                                             # Define start time of time-propagation code
t_end = perf_counter()                                                                                      # Define end time of time-propagation code
output_info.write("Elapsed time after generating spatial derivative of steady state: " + str(t_end-t_start) +'\n')                               # Write the total elapsed computer time for time-propagation code to run
output_info.close()
output_info = open('simulation_info.dat','a')

# -----------------------------------------------------------------------
#  CALCULATING MARKOVIAN FRICTION COEFFICIENT WITH HEOM TIME-PROPAGATION  
# -----------------------------------------------------------------------

if nondiag_key:
    pair_info_row_markov = pair_info_row_nondiag_fil
    pair_info_col_markov = pair_info_col_nondiag_fil
    pair_values_markov = pair_values_nondiag_fil
    npairs_markov = npairs_nondiag_fil
    nnz_elements_sparse_markov = nnz_elements_sparse_nondiag_fil
    sparse_trace_array_markov = sparse_trace_array_nondiag
    nnz_elements_sparse_zeroth_tier_markov = nnz_elements_sparse_zeroth_tier_nondiag_fil
    complex_coefficients_markov = complex_coefficients_nondiag
    isreal_sparse_markov = isreal_sparse_nondiag
    trace_cols_markov = trace_cols_nondiag
    rho_nonzeros_sparse_markov = rho_nonzeros_sparse_nondiag
    rho_ss_spatial_derivative_markov = rho_ss_spatial_derivative_nondiag
    rho_ss_x_arr_markov = rho_ss_x_arr_nondiag
    atol_vec_markov = atol_vec_nondiag
    rtol_vec_markov = rtol_vec_nondiag
    is_connected_array_markov = is_connected_array_nondiag
else:
    pair_info_row_markov = pair_info_row_fil
    pair_info_col_markov = pair_info_col_fil
    pair_values_markov = pair_values_fil
    npairs_markov = npairs_fil
    nnz_elements_sparse_markov = nnz_elements_sparse_fil
    sparse_trace_array_markov = sparse_trace_array
    nnz_elements_sparse_zeroth_tier_markov = nnz_elements_sparse_zeroth_tier_fil
    complex_coefficients_markov = complex_coefficients
    isreal_sparse_markov = isreal_sparse
    trace_cols_markov = trace_cols
    rho_nonzeros_sparse_markov = rho_nonzeros_sparse
    rho_ss_spatial_derivative_markov = rho_ss_spatial_derivative
    rho_ss_x_arr_markov = rho_ss_x_arr
    atol_vec_markov = atol_vec
    rtol_vec_markov = rtol_vec
    is_connected_array_markov = is_connected_array

friction_key = True
if not redo_everything:
    if not bool(redo_friction):
        friction_key = False
    else:
        if os.path.isfile("friction.dat"):
            if not bool(redo_friction):
                friction_key = False

if friction_key:
    friction_object = markovian_friction.calculate_markovian_friction(pair_info_col=pair_info_col_markov,
        pair_info_row=pair_info_row_markov,pair_values=pair_values_markov,npairs=npairs_markov,
        nnz_elements_sparse=nnz_elements_sparse_markov,sparse_trace_array=sparse_trace_array_markov,
        nnz_elements_sparse_zeroth_tier=nnz_elements_sparse_zeroth_tier_markov,complex_coefficients=complex_coefficients_markov,
        isreal_sparse=isreal_sparse_markov,trace_cols=trace_cols_markov,rho_nonzeros_sparse=rho_nonzeros_sparse_markov,tier_index=tier_index,
        rho_ss_spatial_derivative=rho_ss_spatial_derivative_markov,deriv_ham=deriv_Ham,deriv_ham_log=deriv_Ham_log,
        len_un_ind=len_un_ind,atol_vec=atol_vec_markov,rtol_vec=rtol_vec_markov,n_prony_terms=n_prony_terms,nmodes=Nmodes,
        un_ind=Un_Ind,ksiglm=KsigLm,d_ops_log=d_ops_comp_log,d_ops=d_ops_comp)
else:
    print("Markovian friction coefficient has not been (re)generated")

sys.stdout.flush()
t_start = t_end                                                                                             # Define start time of time-propagation code
t_end = perf_counter()                                                                                      # Define end time of time-propagation code
output_info.write("Elapsed time after calculating Markovian friction coefficient: " + str(t_end-t_start) +'\n')                               # Write the total elapsed computer time for time-propagation code to run
output_info.close()
output_info = open('simulation_info.dat','a')

# -----------------------------------------------------------------------
#  CALCULATING MARKOVIAN DIFFFUSION COEFFICIENT WITH HEOM TIME-PROPAGATION  
# -----------------------------------------------------------------------

corrfunc_key = True
if not redo_everything:
    if not bool(redo_corrfunc):
        corrfunc_key = False
    else:
        if os.path.isfile("corrfunc.dat"):
            if not bool(redo_corrfunc):
                corrfunc_key = False

if corrfunc_key:
    corrfunc_object = markovian_corrfunc.calculate_markovian_corrfunc(pair_info_col=pair_info_col_markov,
        pair_info_row=pair_info_row_markov,pair_values=pair_values_markov,npairs=npairs_markov,
        nnz_elements_sparse=nnz_elements_sparse_markov,sparse_trace_array=sparse_trace_array_markov,
        nnz_elements_sparse_zeroth_tier=nnz_elements_sparse_zeroth_tier_markov,
        complex_coefficients=complex_coefficients_markov,isreal_sparse=isreal_sparse_markov,
        trace_cols=trace_cols_markov,rho_nonzeros_sparse=rho_nonzeros_sparse_markov,
        tier_index=tier_index,deriv_ham=deriv_Ham,len_un_ind=len_un_ind,atol_vec=atol_vec_markov,
        rtol_vec=rtol_vec_markov,identity_dim_rho=identity_dim_rho,rho_ss_x_arr=rho_ss_x_arr_markov,
        ksiglm=KsigLm,index_minus=Index_Minus,index_plus=Index_Plus,d_ops=d_ops,d_ops_comp=d_ops_comp,
        eta_vec=eta_vec,rho_sparsity=rho_sparsity,rho_nonzeros=rho_nonzeros,nnz_elements=nnz_elements,
        len_index_plus=len_index_plus,len_index_minus=len_un_ind,nmodes=Nmodes,npoles=Npoles,
        ham_log=Ham_log,d_ops_log=d_ops_log,d_ops_comp_log=d_ops_comp_log,
        average_electronic_force_mol_vec=average_electronic_force_mol_vec,
        average_electronic_force_molleads_vec=average_electronic_force_molleads_vec,
        is_connected_array=is_connected_array_markov,deriv_ham_log=deriv_Ham_log,un_ind=Un_Ind)
else:
    print("Markovian correlation function has not been (re)generated")

sys.stdout.flush()
t_start = t_end                                                                                             # Define start time of time-propagation code
t_end = perf_counter()                                                                                      # Define end time of time-propagation code
output_info.write("Elapsed time after calculating Markovian correlation function: " + str(t_end-t_start) +'\n')                               # Write the total elapsed computer time for time-propagation code to run

t_end = asctime()                                                                                           
output_info.write('The simulation ends at '+str(t_end)+'\n')
output_info.close()                                   
