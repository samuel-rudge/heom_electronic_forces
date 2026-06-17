import numpy as np
import scipy
import psutil
import gc
import os
from time import perf_counter, asctime
import scipy.sparse as sparse
from input_parameters import *
import sys

gc.enable()

class gmres_counter(object):
    def __init__(self, disp=True):
        self._disp = disp
        self.niter = 0
    def __call__(self,rk=None):
        self.niter += 1
        if self._disp:
            print('iter %3i\trk = %s' % (self.niter, str(rk)))

class steady_state_spatial_derivative():

    def __init__(self,pair_info_col,pair_info_row,pair_values,npairs,nnz_elements_sparse,
                 sparse_trace_array,rhs_vector,nnz_elements_sparse_zeroth_tier,
                 complex_coefficients,isreal_sparse,trace_cols,
                 rho_ss_x_arr):    
        
        self.pair_info_col = pair_info_col
        self.pair_info_row = pair_info_row
        self.pair_values = pair_values
        self.npairs = npairs
        self.sparse_trace_array = sparse_trace_array
        self.nnz_elements_sparse = nnz_elements_sparse
        self.rhs_vector = rhs_vector
        self.nnz_elements_sparse_zeroth_tier = nnz_elements_sparse_zeroth_tier
        self.complex_coefficients = complex_coefficients
        self.isreal_sparse = isreal_sparse
        self.trace_cols = trace_cols
        self.rho_ss_x_arr = rho_ss_x_arr

        if parallelize_x_grid:
            self.generate_ss_spatial_derivative_parallel()
        else:
            self.generate_ss_spatial_derivative_sequential()

    def generate_ss_spatial_derivative_sequential(self):

        self.rho_ss_spatial_derivative = np.zeros((len_x_vec,self.nnz_elements_sparse+1),dtype=float)
        self.rho_ss_spatial_derivative[:,0] = x_vec
        coeff_vec = np.array([-1.0,8.0,-8.0,1.0],dtype=float)
        self.process_x_grid = psutil.Process(os.getpid())
        for itrx in range(len_x_vec):
            t_x_start = perf_counter()
            initial_guess_dx = self.rho_ss_x_arr[itrx,1:]
            rho_deriv_x = np.zeros(self.nnz_elements_sparse,dtype=float)
            self.base_memory_usage = self.process_x_grid.memory_info().rss
            self.counter_gmres = gmres_counter(disp=False)
            for itrdx in range(4):
                print("Finding steady state of itrdx ",itrdx)
                sys.stdout.flush()
                pair_values_dx = self.pair_values[:,itrx,itrdx+1]
                steady_state_dx = self.steady_state_one_dx_value(itrx,pair_values_dx,initial_guess_dx)
                rho_deriv_x += steady_state_dx*coeff_vec[itrdx]/(12*dx_sssd)
            self.rho_ss_spatial_derivative[itrx,1:] = rho_deriv_x
            t_x_end = perf_counter()
            print("Time for x value ",itrx," is ",t_x_end-t_x_start)
            sys.stdout.flush()

    def steady_state_one_dx_value(self,itrx,pair_values_dx,initial_guess_dx):
            
            sparse_heom_generator_dx = sparse.csc_matrix((pair_values_dx,(self.pair_info_row,self.pair_info_col)),
                                            shape=(self.nnz_elements_sparse,self.nnz_elements_sparse),dtype=float)
            sparse_heom_w_trace_dx = sparse_heom_generator_dx + self.sparse_trace_array
            preconditioner_values = sparse_heom_generator_dx.diagonal()
            preconditioner_values[0:self.nnz_elements_sparse_zeroth_tier] = 1.0
            preconditioner_values[self.nnz_elements_sparse_zeroth_tier:] = \
                        1.0/preconditioner_values[self.nnz_elements_sparse_zeroth_tier:]
            preconditioner_rows = np.arange(0,self.nnz_elements_sparse)
            preconditioner_cols = preconditioner_rows
            sparse_preconditioner = sparse.csc_matrix((preconditioner_values,(preconditioner_rows,preconditioner_cols)),
                                            shape=(self.nnz_elements_sparse,self.nnz_elements_sparse),dtype=float)
            ### GENERATE STEADY STATE AT 1 GRID POINT ### 
            steady_state_gmres,success_yn = \
                sparse.linalg.gmres(sparse_heom_w_trace_dx,self.rhs_vector,x0=initial_guess_dx,
                                    tol=tol_gmres_ss,M=sparse_preconditioner,restart=dim_krylov_space,
                                    maxiter=dim_krylov_space,callback=self.counter_gmres)
            print("Number of GMRES iterations: ",self.counter_gmres.niter)
            print("Found the steady state? ",not bool(success_yn))
            memory_usage = self.process_x_grid.memory_info().rss
            loop_memory_usage = memory_usage - self.base_memory_usage
            print("Memory for x value ",itrx," is ",loop_memory_usage/1024**2,"MB")
            sys.stdout.flush()
            gc.collect()

            return steady_state_gmres            

    def generate_ss_spatial_derivative_parallel(self):

        print("Not implemented yet")
        sys.stdout.flush()

    def return_ss_spatial_derivative(self):

        return self.rho_ss_spatial_derivative

class steady_state_spatial_derivative_nondiag():
    
    def __init__(self,rho_ss_spatial_derivative,Kmatrix,Kmatrix_inv,
                 is_connected_array_nondiag,len_un_ind,rho_nonzeros_nondiag,
                 complex_coefficients,nnz_elements_sparse,rho_nonzeros_sparse,
                 nnz_elements_sparse_nondiag_fil,nnz_elements_nondiag):
        
        self.rho_ss_spatial_derivative = rho_ss_spatial_derivative
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

        self.generate_ss_spatial_derivatve_nondiagonal_basis()
        self.generate_ss_spatial_derivatve_nondiagonal_basis_sparse()

    def generate_ss_spatial_derivatve_nondiagonal_basis(self):

        self.rho_ss_spatial_derivative_nondiag_hilbert_space = np.zeros((len_x_vec,dim_rho,dim_rho,self.len_un_ind),dtype=complex)
        for itrx in range(len_x_vec):
            rho_spatial_derivative_one_x = self.rho_ss_spatial_derivative[itrx,1:]
            rho_spatial_derivative_one_x_hilbert_space = np.zeros((dim_rho,dim_rho,self.len_un_ind),dtype=complex)
            for itrnz in range(self.nnz_elements_sparse):
                ado_ind = self.rho_nonzeros_sparse[itrnz,0]
                nrow = self.rho_nonzeros_sparse[itrnz,1]
                ncol = self.rho_nonzeros_sparse[itrnz,2]
                real_complex_val = self.complex_coefficients[itrnz]
                rho_spatial_derivative_one_x_hilbert_space[nrow,ncol,ado_ind] += real_complex_val*\
                            rho_spatial_derivative_one_x[itrnz]
            rho_ss_spatial_derivative_one_x_hilbert_space_nondiag = np.zeros((dim_rho,dim_rho,self.len_un_ind),dtype=complex)
            for itr_ado in range(self.len_un_ind):
                rho_ss_spatial_derivative_one_x_hilbert_space_nondiag[:,:,itr_ado] = np.matmul(self.Kmatrix_inv,\
                                            np.matmul(rho_spatial_derivative_one_x_hilbert_space[:,:,itr_ado],\
                                            self.Kmatrix))
            self.rho_ss_spatial_derivative_nondiag_hilbert_space[itrx,:,:,:] = rho_ss_spatial_derivative_one_x_hilbert_space_nondiag

    def generate_ss_spatial_derivatve_nondiagonal_basis_sparse(self):

        rho_ss_spatial_derivative_nondiag_uf = np.zeros((len_x_vec,2*self.nnz_elements_nondiag),dtype=float)
        for itrnz in range(self.nnz_elements_nondiag):
            ind_ado = self.rho_nonzeros_nondiag[itrnz,0]
            row = self.rho_nonzeros_nondiag[itrnz,1]
            col = self.rho_nonzeros_nondiag[itrnz,2]            
            rho_ss_spatial_derivative_nondiag_uf[:,2*itrnz] = \
                np.real(self.rho_ss_spatial_derivative_nondiag_hilbert_space[:,row,col,ind_ado])
            rho_ss_spatial_derivative_nondiag_uf[:,2*itrnz+1] = \
                np.imag(self.rho_ss_spatial_derivative_nondiag_hilbert_space[:,row,col,ind_ado])
        
        self.rho_ss_spatial_derivative_nondiag = np.zeros((len_x_vec,self.nnz_elements_sparse_nondiag_fil+1),dtype=float)
        self.rho_ss_spatial_derivative_nondiag[:,0] = x_vec
        self.rho_ss_spatial_derivative_nondiag[:,1:] = rho_ss_spatial_derivative_nondiag_uf[:,self.is_connected_array_nondiag]

    def return_ss_spatial_derivative_nondiag(self):

        return self.rho_ss_spatial_derivative_nondiag



