"""
------------------------------------------------------------------------------
STEADY STATE SOLUTION ON X-GRID (GMRES-BASED SPARSE HEOM SOLVER)
------------------------------------------------------------------------------

This module computes the instantaneous electronic steady state of a
position-dependent HEOM (hierarchical equations of motion) problem
for a coupled electron-vibration system on a 1D nuclear coordinate grid.

The steady state at each grid point x is obtained by solving the sparse linear
system:

    (L_HEOM(x) + Trace_Constraint) * rho_ss(x) = rhs

using a preconditioned GMRES solver (scipy.sparse.linalg.gmres).

Core responsibilities:
- Construct and solve the sparse steady-state HEOM equation at each x-grid point
- Compute observables:
  - Reduced density matrix rho(x)
  - Electronic forces (molecular + lead contributions)
  - Lead currents
  - Optional vibrational excitations (if quantum vibrational modes are enabled)
- Use physically informed initial guesses (previous x-point, voltage continuation,
  partition continuation when available)
- Store full x-dependent steady-state solution and derived observables

Key components:
- gmres_counter: iteration monitor for GMRES solver
- steady_state_x_grid:
    - sequential or parallel evaluation over x-grid
    - sparse GMRES solve per grid point
    - extraction of physical observables
    - output to .dat files

Dependencies:
- sparsity-generated HEOM Liouvillian structure
- system operators and Hamiltonian from system.py
- global simulation parameters from input_parameters.py

Main outputs:
- rho_ss_x_arr: sparse steady-state vector over x-grid
- average_electronic_force_*: force contributions
- current_vec: lead currents
- rho_system: reconstructed density matrices
- optional vibrational observables

------------------------------------------------------------------------------
"""

import numpy as np
import scipy
import psutil
import gc
import os
import sys
from time import perf_counter, asctime
import scipy.sparse as sparse
from input_parameters import *

gc.enable()

class gmres_counter(object):
    def __init__(self, disp=True):
        self._disp = disp
        self.niter = 0
    def __call__(self,rk=None):
        self.niter += 1
        if self._disp:
            print('iter %3i\trk = %s' % (self.niter, str(rk)))

class steady_state_x_grid():

    def __init__(self,pair_info_col,pair_info_row,pair_values,npairs,nnz_elements_sparse,
                 sparse_trace_array,rhs_vector,nnz_elements_sparse_zeroth_tier,
                 complex_coefficients,d_ops,
                 isreal_sparse,deriv_ham,deriv_ham_log,rho_nonzeros_sparse,tier_index,
                 un_ind,ksiglm,nvib_object,ham_x,trace_cols,fock_states):

        self.pair_info_col = pair_info_col
        self.pair_info_row = pair_info_row
        self.pair_values = pair_values
        self.npairs = npairs
        self.sparse_trace_array = sparse_trace_array
        self.nnz_elements_sparse = nnz_elements_sparse
        self.rhs_vector = rhs_vector
        self.nnz_elements_sparse_zeroth_tier = nnz_elements_sparse_zeroth_tier
        self.complex_coefficients = complex_coefficients
        self.d_ops = d_ops
        self.isreal_sparse = isreal_sparse
        self.deriv_ham = deriv_ham
        self.deriv_ham_log = deriv_ham_log
        self.rho_nonzeros_sparse = rho_nonzeros_sparse
        self.tier_index = tier_index
        self.un_ind = un_ind
        self.ksiglm = ksiglm
        self.ham_x = ham_x
        self.trace_cols = trace_cols
        self.fock_states = fock_states
        if bool(N_qu_vib_modes):
            self.nvib_qu_diag = nvib_object[0]
            self.nvib_qu_nondiag = nvib_object[1]
            self.nvib_qu_diag_log = nvib_object[2]
            self.nvib_qu_nondiag_log = nvib_object[3]
            self.dim_vib_mode_qu = nvib_object[4]
            if small_polaron_yn:
                self.Kmatrix = nvib_object[5]
                self.Kmatrix_inv = nvib_object[6]
        
        self.voltage_below_complete_yn = False
        self.partition_below_complete_yn = False
        if bool(checking_partitions_and_voltages):
            current_dir = os.getcwd()
            self.voltage_below_dir = f"{current_dir}/../../voltage_{(voltage-dv):.2f}eV/x_partition_{(x_partition_count):.0f}/"
            self.partition_below_dir = f"{current_dir}/../x_partition_{(x_partition_count-1):.0f}/"
            self.partition_below_complete_yn = os.isfile(f"{self.partition_below_dir}/rho_ss_x_arr.npy")
            self.voltage_below_complete_yn = os.isfile(f"{self.voltage_below_dir}/rho_ss_x_arr.npy")

        if parallelize_x_grid:
            self.generate_average_electronic_force_parallel()
        else:
            self.generate_average_electronic_force_sequential()

    def generate_average_electronic_force_sequential(self):

        self.rho_ss_x_arr = np.zeros((len_x_vec,self.nnz_elements_sparse+1),dtype=float)
        self.rho_ss_x_arr[:,0] = x_vec
        self.rho_system = np.zeros((len_x_vec,dim_rho,dim_rho),dtype=complex)
        self.average_electronic_force_mol_vec = np.zeros((len_x_vec,2),dtype=complex)
        self.average_electronic_force_mol_vec[:,0] = x_vec
        self.average_electronic_force_molleads_vec = np.zeros((len_x_vec,2),dtype=complex)
        self.average_electronic_force_molleads_vec[:,0] = x_vec
        self.average_electronic_force_vec = np.zeros((len_x_vec,2),dtype=complex)
        self.average_electronic_force_vec[:,0] = x_vec
        self.current_vec = np.zeros((len_x_vec,Nleads+1),dtype=float)
        self.current_vec[:,0] = x_vec
        if bool(N_qu_vib_modes):
            self.average_vib_excitation_vec_diag = np.zeros((len_x_vec,N_qu_vib_modes+1),dtype=float)
            self.average_vib_excitation_vec_diag[:,0] = x_vec
            self.average_vib_excitation_vec_nondiag = np.zeros((len_x_vec,N_qu_vib_modes+1),dtype=float)
            self.average_vib_excitation_vec_nondiag[:,0] = x_vec
        self.process_x_grid = psutil.Process(os.getpid())
        self.base_memory_usage = self.process_x_grid.memory_info().rss
        for itrx in range(len_x_vec):
            t_x_start = perf_counter()
            self.counter_gmres = gmres_counter(disp=False)
            initial_guess_x = self.generate_initial_guess_x(itrx)
            pair_values_x = self.pair_values[:,itrx,0]
            deriv_ham_x = self.deriv_ham[:,:,0,itrx]
            dV_dx_x = dV_dx[:,:,itrx]
            steady_state_x = self.steady_state_one_x_value(itrx,pair_values_x,initial_guess_x)
            self.rho_ss_x_arr[itrx,1:] = steady_state_x
            steady_state_quantities_x = self.steady_state_quantities_one_x_value(steady_state_x,deriv_ham_x,dV_dx_x)
            self.average_electronic_force_mol_vec[itrx,1:] = steady_state_quantities_x[0]
            self.average_electronic_force_molleads_vec[itrx,1:] = steady_state_quantities_x[1]
            self.average_electronic_force_vec[itrx,1:] = steady_state_quantities_x[0] + steady_state_quantities_x[1]
            self.current_vec[itrx,1:] = steady_state_quantities_x[2]*(6.623618237510/27.211386245988)*(10**3)
            self.rho_system[itrx,:,:] = steady_state_quantities_x[3]
            # print(steady_state_quantities_x[3])
            if bool(N_qu_vib_modes):
                self.average_vib_excitation_vec_diag[itrx,1:] = steady_state_quantities_x[4]
                self.average_vib_excitation_vec_nondiag[itrx,1:] = steady_state_quantities_x[5]
            t_x_end = perf_counter()
            print("Time for x value ",itrx," is ",t_x_end-t_x_start)
            sys.stdout.flush()
        adiabatic_force_mol_file = "adiabatic_force_mol.dat"
        adiabatic_force_molleads_file = "adiabatic_force_molleads.dat"
        adiabatic_force_file = "adiabatic_force.dat"
        adiabatic_force_headings = 'Coord. Val.\tAdiabatic Force'
        current_ad_file = "current_ad.dat"
        current_ad_headings = 'Coord. Val.\tCurrent_0\tCurrent_1'
        if bool(N_qu_vib_modes):
            average_vib_excitation_diag_file = "average_vib_excitation_diag.dat"
            average_vib_excitation_nondiag_file = "average_vib_excitation_nondiag.dat"
            average_vib_excitation_headings = 'Coord. Val.\tVib. Excitation.'
        np.savetxt(adiabatic_force_mol_file,np.real(self.average_electronic_force_mol_vec),fmt='%f',delimiter='\t',header=adiabatic_force_headings,comments='')
        np.savetxt(adiabatic_force_molleads_file,np.real(self.average_electronic_force_molleads_vec),fmt='%f',delimiter='\t',header=adiabatic_force_headings,comments='')
        np.savetxt(adiabatic_force_file,np.real(self.average_electronic_force_vec),fmt='%f',delimiter='\t',header=adiabatic_force_headings,comments='')
        np.savetxt(current_ad_file,self.current_vec,fmt='%f',delimiter='\t',header=current_ad_headings,comments='')
        if bool(N_qu_vib_modes):
            np.savetxt(average_vib_excitation_diag_file,self.average_vib_excitation_vec_diag,fmt='%f',delimiter='\t',header=average_vib_excitation_headings,comments='')
            np.savetxt(average_vib_excitation_nondiag_file,self.average_vib_excitation_vec_nondiag,fmt='%f',delimiter='\t',header=average_vib_excitation_headings,comments='')

    def generate_average_electronic_force_parallel(self):
        
        print("Not yet implemented")
        sys.stdout.flush()

    def steady_state_one_x_value(self,itrx,pair_values_x,initial_guess_x):
            
            sparse_heom_generator_x = sparse.csc_matrix((pair_values_x,(self.pair_info_row,self.pair_info_col)),
                                            shape=(self.nnz_elements_sparse,self.nnz_elements_sparse),dtype=float)
            sparse_heom_w_trace_x = sparse_heom_generator_x + self.sparse_trace_array
            preconditioner_values = sparse_heom_generator_x.diagonal()
            preconditioner_values[0:self.nnz_elements_sparse_zeroth_tier] = 1.0
            preconditioner_values[self.nnz_elements_sparse_zeroth_tier:] = \
                        1.0/preconditioner_values[self.nnz_elements_sparse_zeroth_tier:]
            preconditioner_rows = np.arange(0,self.nnz_elements_sparse)
            preconditioner_cols = preconditioner_rows
            sparse_preconditioner = sparse.csc_matrix((preconditioner_values,(preconditioner_rows,preconditioner_cols)),
                                            shape=(self.nnz_elements_sparse,self.nnz_elements_sparse),dtype=float)
            ### GENERATE STEADY STATE AT 1 GRID POINT ### 
            steady_state_gmres,success_yn = \
                sparse.linalg.gmres(A=sparse_heom_w_trace_x,b=self.rhs_vector,x0=initial_guess_x,
                                    tol=tol_gmres_ss,M=sparse_preconditioner,restart=dim_krylov_space,
                                    maxiter=dim_krylov_space,callback=self.counter_gmres)
            print("Number of GMRES iterations: ",self.counter_gmres.niter)
            print("Found the steady state? ",not bool(success_yn))
            memory_usage = self.process_x_grid.memory_info().rss
            loop_memory_usage = memory_usage - self.base_memory_usage
            print("Memory for x value ",itrx," is ",loop_memory_usage/1024**2,"MB")
            sys.stdout.flush()
            
            return steady_state_gmres
            
    def steady_state_quantities_one_x_value(self,steady_state_x,deriv_ham_x,dV_dx_x):

        rho_system_x = np.zeros((dim_rho,dim_rho),dtype=complex)
        current_x = np.zeros(Nleads,dtype=float)
        average_electronic_force_x_mol = 0.0
        average_electronic_force_x_molleads = 0.0
        if bool(N_qu_vib_modes):
            average_vib_excitation_diag_x = np.zeros(N_qu_vib_modes,dtype=float)
            average_vib_excitation_nondiag_x = np.zeros(N_qu_vib_modes,dtype=float)
        for itrnz in range(self.nnz_elements_sparse):
            indjn = self.rho_nonzeros_sparse[itrnz,0]
            itrn = self.tier_index[indjn]
            nrow = self.rho_nonzeros_sparse[itrnz,1]
            ncol = self.rho_nonzeros_sparse[itrnz,2]
            if (itrn == 0):
                rho_system_x[nrow,ncol] += steady_state_x[itrnz]*self.complex_coefficients[itrnz]
                if (self.deriv_ham_log[ncol,nrow,0,0]):
                    if self.isreal_sparse[itrnz]:
                        average_electronic_force_x_mol -= deriv_ham_x[ncol,nrow]*steady_state_x[itrnz]
                    if bool(wbl_YN):
                        if self.isreal_sparse[itrnz]:
                            for itrlead in range(Nleads):
                                if (degenerate_levels == False):
                                    for itrel1 in range(self.nel):
                                        for itrel2 in range(self.nel):
                                            sign_val = 0.0
                                            for ndash in range(dim_rho):
                                                sign_val = sign_val + self.d_ops[ndash,nrow,itrel1,0]*self.d_ops[ncol,ndash,itrel2,1] \
                                                                    - self.d_ops[ndash,nrow,itrel2,1]*self.d_ops[ncol,ndash,itrel1,0]
                                            current_x[itrlead] += sign_val*np.pi*el_lead_couplings[itrlead,itrel1]*\
                                                                el_lead_couplings[itrlead,itrel2]*steady_state_x[itrnz]
                                else:
                                    for itrel in range(self.nel):
                                        sign_val = 0.0
                                        for ndash in range(dim_rho):
                                            sign_val = sign_val + self.d_ops[ndash,nrow,itrel,0]*self.d_ops[ncol,ndash,itrel,1] \
                                                                - self.d_ops[ndash,nrow,itrel,1]*self.d_ops[ncol,ndash,itrel,0]
                                        current_x[itrlead] += sign_val*np.pi*(el_lead_couplings[itrlead,itrel]**2)*\
                                                                            self.steady_state_x[itrnz]
                        if bool(N_qu_vib_modes):
                            for itr_nvib_qu in range(N_qu_vib_modes):
                                if (self.nvib_qu_diag_log[ncol,nrow,N_qu_vib_modes]):
                                        average_vib_excitation_diag_x[itr_nvib_qu] += \
                                            self.nvib_qu_diag[ncol,nrow,itr_nvib_qu]*steady_state_x[itrnz]
                                if (self.nvib_qu_nondiag_log[ncol,nrow,N_qu_vib_modes]):
                                        average_vib_excitation_nondiag_x[itr_nvib_qu] += \
                                            self.nvib_qu_nondiag[ncol,nrow,itr_nvib_qu]*steady_state_x[itrnz]
            elif (itrn == 1):
                    jn = self.un_ind[indjn,0]
                    leads_n = self.ksiglm[jn,0]                                                         
                    sign_n = 1-self.ksiglm[jn,1]                                                         
                    eldash_n = self.ksiglm[jn,3]
                    if self.isreal_sparse[itrnz]:
                        average_electronic_force_x_molleads += -dV_dx_x[leads_n,eldash_n]*2.0*\
                                                self.d_ops[ncol,nrow,eldash_n,sign_n]*steady_state_x[itrnz]
                    else:
                        if not degenerate_levels:
                            for itrel in range(self.Nel):
                                current_x[leads_n] += ((-1.0)**sign_n)*2.0*el_lead_couplings[leads_n,eldash_n,itrx,0]*self.d_ops[ncol,nrow,itrel,sign_n]*steady_state_x[itrnz]
                        else:
                            current_x[leads_n] += ((-1.0)**sign_n)*2.0*el_lead_couplings[leads_n,eldash_n,itrx,0]*self.d_ops[ncol,nrow,eldash_n,sign_n]*steady_state_x[itrnz]
            else:
                break

        if not bool(N_qu_vib_modes):
            return average_electronic_force_x_mol,average_electronic_force_x_molleads,current_x,rho_system_x
        else:
            return average_electronic_force_x_mol,average_electronic_force_x_molleads,current_x,rho_system_x,\
                    average_vib_excitation_diag_x,average_vib_excitation_nondiag_x

    def generate_initial_guess_x(self,itrx):

        initial_guess_key = False
        if itrx > 0:
            initial_guess_x = self.rho_ss_x_arr[itrx-1,1:]
            initial_guess_key = True
        if self.partition_below_complete_yn and not initial_guess_key:
            ind_x_total_partition_below = np.load(f"{self.partition_below_dir}/ind_x_total.npy")
            matching_ind_x_partion_below = np.array((ind_x_total-1) == ind_x_total_partition_below,dtype=bool)
            if (matching_ind_x_partion_below).any():
                initial_guess_x = \
                    np.load(f"{self.partition_below_dir}/rho_ss_x_arr.npy")[matching_ind_x_partion_below,1:]
                initial_guess_key = True
        if self.voltage_below_complete_yn and not initial_guess_key:
            x_vec_voltage_below = np.load(f"{self.voltage_below_dir}/rho_ss_x_arr.npy")[:,0]
            matching_x_vec_voltage_below = np.isclose(x_vec,x_vec_voltage_below)
            if (matching_x_vec_voltage_below).any():
                initial_guess_x = np.load(f"{self.voltage_below_dir}/rho_ss_x_arr.npy")[matching_x_vec_voltage_below,1:]
                initial_guess_key = True
        if not initial_guess_key:
            initial_guess_x = np.zeros(self.nnz_elements_sparse,dtype=float)
            if (np.abs(self.ham_x[1,1,itrx,0]) <= voltage/2):
                initial_guess_x[self.trace_cols[0:2]] = 1/dim_el
            elif (self.ham_x[1,1,itrx,0] > voltage/2):
                if not bool(N_qu_vib_modes):
                    initial_guess_x[self.trace_cols[self.fock_states[:,0]==0]] = 1
                else:
                    initial_guess_x[0] = 1.0
            elif (self.ham_x[1,1,itrx,0] < -voltage/2):
                if not bool(N_qu_vib_modes):
                    initial_guess_x[self.trace_cols[self.fock_states[:,0]==1]] = 1
                else:
                    initial_guess_x[1] = 1.0

        return initial_guess_x

    def return_ss_all_x_values(self):

        if bool(N_qu_vib_modes):
            steady_state_results = [self.rho_ss_x_arr,self.rho_system,\
                    self.average_electronic_force_mol_vec,
                    self.average_electronic_force_molleads_vec,self.current_vec,\
                    self.average_vib_excitation_vec_diag,self.average_vib_excitation_vec_nondiag]
        else:
            steady_state_results = [self.rho_ss_x_arr,self.rho_system,\
                                    self.average_electronic_force_mol_vec,
                                    self.average_electronic_force_molleads_vec,self.current_vec]
                    
        return steady_state_results

class steady_state_nondiag():
    
    def __init__(self,rho_ss_x_arr,Kmatrix,Kmatrix_inv,
                 is_connected_array_nondiag,len_un_ind,rho_nonzeros_nondiag,
                 complex_coefficients,nnz_elements_sparse,rho_nonzeros_sparse,
                 nnz_elements_sparse_nondiag_fil,nnz_elements_nondiag):
        
        self.rho_ss_x_arr = rho_ss_x_arr
        self.Kmatrix = Kmatrix
        self.Kmatrix_inv = Kmatrix_inv
        self.is_connected_array_nondiag = is_connected_array_nondiag
        self.nnz_elements_sparse = nnz_elements_sparse
        self.len_un_ind = len_un_ind
        self.rho_nonzeros_nondiag = rho_nonzeros_nondiag
        self.complex_coefficients = complex_coefficients
        self.nnz_elements_sparse_nondiag_fil = nnz_elements_sparse_nondiag_fil
        self.nnz_elements_nondiag = nnz_elements_nondiag
        self.rho_nonzeros_sparse = rho_nonzeros_sparse

        self.generate_ss_nondiagonal_basis()
        self.generate_ss_nondiagonal_basis_sparse()

    def generate_ss_nondiagonal_basis(self):

        self.rho_ss_nondiag_hilbert_space = np.zeros((len_x_vec,dim_rho,dim_rho,self.len_un_ind),dtype=complex)
        for itrx in range(len_x_vec):
            rho_ss_one_x = self.rho_ss_x_arr[itrx,1:]
            rho_ss_one_x_hilbert_space = np.zeros((dim_rho,dim_rho,self.len_un_ind),dtype=complex)
            for itrnz in range(self.nnz_elements_sparse):
                ado_ind = self.rho_nonzeros_sparse[itrnz,0]
                nrow = self.rho_nonzeros_sparse[itrnz,1]
                ncol = self.rho_nonzeros_sparse[itrnz,2]
                real_complex_val = self.complex_coefficients[itrnz]
                rho_ss_one_x_hilbert_space[nrow,ncol,ado_ind] += real_complex_val*rho_ss_one_x[itrnz]
            rho_ss_one_x_hilbert_space_nondiag = np.zeros((dim_rho,dim_rho,self.len_un_ind),dtype=complex)
            for itr_ado in range(self.len_un_ind):
                rho_ss_one_x_hilbert_space_nondiag[:,:,itr_ado] = np.matmul(self.Kmatrix_inv,\
                                                            np.matmul(rho_ss_one_x_hilbert_space[:,:,itr_ado],\
                                                            self.Kmatrix))
            self.rho_ss_nondiag_hilbert_space[itrx,:,:,:] = rho_ss_one_x_hilbert_space_nondiag

    def generate_ss_nondiagonal_basis_sparse(self):

        rho_ss_nondiag_uf = np.zeros((len_x_vec,2*self.nnz_elements_nondiag),dtype=float)
        for itrnz in range(self.nnz_elements_nondiag):
            ind_ado = self.rho_nonzeros_nondiag[itrnz,0]
            row = self.rho_nonzeros_nondiag[itrnz,1]
            col = self.rho_nonzeros_nondiag[itrnz,2]
            rho_ss_nondiag_uf[:,2*itrnz] = \
                np.real(self.rho_ss_nondiag_hilbert_space[:,row,col,ind_ado])
            rho_ss_nondiag_uf[:,2*itrnz+1] = \
                np.imag(self.rho_ss_nondiag_hilbert_space[:,row,col,ind_ado])
        
        self.rho_ss_nondiag_fil = np.zeros((len_x_vec,self.nnz_elements_sparse_nondiag_fil+1),dtype=float)
        self.rho_ss_nondiag_fil[:,0] = x_vec
        self.rho_ss_nondiag_fil[:,1:] = rho_ss_nondiag_uf[:,self.is_connected_array_nondiag]

    def return_ss_x_arr_nondiag(self):

        return self.rho_ss_nondiag_fil


