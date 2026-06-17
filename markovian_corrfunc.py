# markovian_corrfunc.py
#
# This module computes the Markovian corrfunc kernel from a HEOM-based electronic structure
# calculation on a nuclear coordinate grid.
#
# The core idea is to propagate a linear-response–like auxiliary dynamics (via a sparse HEOM
# Liouvillian constructed in Fortran through sparse_corrfunc.f90), extract the time-dependent
# corrfunc integrand, and then evaluate the corrfunc coefficient by fitting the integrand
# to a sum of exponentials using a Prony decomposition.
#
# Main features:
# - Sequential or (future) parallel evaluation over nuclear coordinates (x-grid)
# - Adaptive time propagation of the corrfunc integrand using sparse HEOM propagation
# - Construction of corrfunc contributions from molecular and lead-resolved channels
# - On-the-fly convergence checking of the time integral
# - Prony-based compression of the corrfunc kernel into exponential modes
# - Storage of corrfunc, current, and convergence diagnostics as x-dependent arrays
#
# External dependencies:
# - sparse_corrfunc.f90 (compiled Fortran backend for HEOM propagation)
# - numpy, scipy (sparse linear algebra, interpolation, optimization)
# - input_parameters.py (global physical and numerical parameters)
#
# Output:
# - corrfunc_mol_vec: molecular contribution to corrfunc
# - corrfunc_molleads_vec: lead-induced contribution
# - corrfunc_vec: total Markovian corrfunc
# - current_na_vec: associated non-adiabatic current
# - prop_info_corrfunc: propagation diagnostics (timesteps, time, convergence, etc.)

import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy
from scipy.interpolate import CubicSpline
import scipy.sparse as sparse
import sparse_corrfunc
import sparsity_corrfunc
from input_parameters import *

class calculate_markovian_corrfunc():

    def __init__(self,pair_info_col,pair_info_row,pair_values,npairs,nnz_elements_sparse,
                    sparse_trace_array,nnz_elements_sparse_zeroth_tier,complex_coefficients,
                    isreal_sparse,trace_cols,rho_nonzeros_sparse,tier_index,
                    deriv_ham,len_un_ind,atol_vec,rtol_vec,identity_dim_rho,rho_ss_x_arr,ksiglm,
                    index_minus,index_plus,d_ops,d_ops_comp,eta_vec,
                    rho_sparsity,rho_nonzeros,nnz_elements,len_index_plus,len_index_minus,nmodes,npoles,
                    ham_log,d_ops_log,d_ops_comp_log,average_electronic_force_mol_vec,average_electronic_force_molleads_vec,
                    is_connected_array,deriv_ham_log,un_ind):

        self.pair_info_col = pair_info_col
        self.pair_info_row = pair_info_row
        self.pair_values = pair_values
        self.npairs = npairs
        self.sparse_trace_array = sparse_trace_array
        self.nnz_elements_sparse = nnz_elements_sparse
        self.nnz_elements_sparse_zeroth_tier = nnz_elements_sparse_zeroth_tier
        self.complex_coefficients = complex_coefficients
        self.isreal_sparse = isreal_sparse
        self.trace_cols = trace_cols
        self.rho_nonzeros_sparse = rho_nonzeros_sparse
        self.tier_index = tier_index
        self.nnz_elements = nnz_elements_sparse
        self.len_un_ind = len_un_ind
        self.deriv_ham = deriv_ham
        self.deriv_ham_log = deriv_ham_log
        self.atol_vec = atol_vec
        self.rtol_vec = rtol_vec
        self.identity_dim_rho = identity_dim_rho
        self.rho_ss_x_arr = rho_ss_x_arr
        self.ksiglm = ksiglm
        self.un_ind = un_ind
        self.index_minus = index_minus
        self.index_plus = index_plus
        self.d_ops = d_ops
        self.d_ops_comp = d_ops_comp
        self.eta_vec = eta_vec
        self.rho_sparsity = rho_sparsity
        self.rho_nonzeros = rho_nonzeros
        self.nnz_elements = nnz_elements
        self.len_index_plus = len_index_plus
        self.len_index_minus = len_index_minus
        self.nmodes = nmodes
        self.npoles = npoles
        self.ham_log = ham_log
        self.d_ops_log=d_ops_log
        self.d_ops_comp_log=d_ops_comp_log
        self.average_electronic_force_mol_vec=average_electronic_force_mol_vec
        self.average_electronic_force_molleads_vec=average_electronic_force_molleads_vec
        self.is_connected_array = is_connected_array

        self.generate_ic_sparsity()
        if parallelize_x_grid:
            self.generate_markovian_corrfunc_parallel()
        else:
            if checking_timesteps_yn:
                self.corrfunction_ic()
            self.generate_markovian_corrfunc_sequential()

    def generate_markovian_corrfunc_parallel(self):

        print("Not yet implemented")

    def generate_markovian_corrfunc_sequential(self):
        
        self.corrfunc_vec = np.zeros((len_x_vec,2),dtype=complex)
        self.corrfunc_vec[:,0] = x_vec
        self.corrfunc_mol_vec = np.zeros((len_x_vec,2),dtype=complex)
        self.corrfunc_mol_vec[:,0] = x_vec
        self.corrfunc_molleads_vec = np.zeros((len_x_vec,2),dtype=complex)
        self.corrfunc_molleads_vec[:,0] = x_vec
        if checking_timesteps_yn:
            if min_time > checking_time_interval:
                raise ValueError("Minimum propagation time cannot be smaller than the time interval between \
                                successive convergence checks")
            self.prop_info_corrfunc = np.zeros((len_x_vec,6),dtype=float)
            self.prop_info_corrfunc[:,0] = self.x_vec
            for itrx in range(len_x_vec):
                print("Calculating Markovian corrfunc coefficient of "+str(itrx)+"th coordinate from "+str(len_x_vec))
                pair_values_x = self.pair_values[:,itrx,0]
                rel_diff_corrfunc = 1.0
                ic_prop = self.rho_0_diff_x_arr[itrx,1:]
                dforce_x = self.dforce[:,:,itrx]
                dt_input = dt_init
                facmax_input = facmax_init
                dt_tot_input = 0.0
                rel_diff_input = 1.0
                markovian_corrfunc_coefficient_input = 0.0
                norm_input = 0.0
                time_input = 0.0
                itrt_input = int(0)
                corrfunc_estimate_old = 1.0
                self.convergence_info = []
                time_vec = []
                corrfunc_integrand_data = []
                continue_checking = True
                first_timecheck = True
                while ((rel_diff_corrfunc >= tol_checking_time) and continue_checking):
                    propagation_outputs = \
                        sparse_corrfunc.sparse_markovian_corrfunc_propagation_w_checking(pair_info_row=self.pair_info_row,
                            pair_info_col=self.pair_info_col,pair_values=pair_values_x,
                            rho_nonzeros=self.rho_nonzeros_sparse,
                            tier_index=self.tier_index,rho_0=ic_prop,tol_corrfunc=tol_corrfunc,
                            dforce=dforce_x,dforce_log=self.dforce_log,max_expan_order=max_expan_order,dim_rho=dim_rho,nnz_elements=self.nnz_elements_sparse,
                            npairs=self.npairs,len_un_ind=self.len_un_ind,nthreads_liouvillian=nthreads_liouvillian,
                            rk_coeff=rk_coeff,rk_coeffhat=rk_coeffhat,
                            facmin=facmin,fac=fac,atol_vec=self.atol_vec,rtol_vec=self.rtol_vec,
                            max_time=max_time,dt_min=dt_min,complex_coefficients=self.complex_coefficients,min_time=min_time,
                            time_input=time_input,dt_tot_input=dt_tot_input,norm_input=norm_input,
                            rel_diff_input=rel_diff_input,itrt_input=itrt_input,dt_input=dt_input,facmax_input=facmax_input,
                            markovian_corrfunc_coefficient_input=markovian_corrfunc_coefficient_input,
                            checking_time_interval=checking_time_interval,first_timecheck=first_timecheck)
                    continue_checking,ic_prop,itrt_input,norm_input,rel_diff_input,\
                    markovian_corrfunc_coefficient_input,time_input,dt_tot_input,dt_input = propagation_outputs
                    if (continue_checking):
                        new_integrand_data = np.loadtxt('corrfunc_integrand_new.dat')
                        time_vec.extend(new_integrand_data[:,0].tolist())
                        corrfunc_integrand_data.extend(new_integrand_data[:,1].tolist())
                        integrand_func = CubicSpline(time_vec,corrfunc_integrand_data)
                        time_vec_even = np.linspace(time_vec[0],time_vec[-1],2*len(time_vec))
                        corrfunc_integrand_even = integrand_func(time_vec_even)
                        weights,frequencies = self.prony_function(time_vec_even,corrfunc_integrand_even,self.n_prony_terms)
                        corrfunc_estimate = -0.5*np.real(np.sum(weights/frequencies))
                        rel_diff_corrfunc = np.abs(corrfunc_estimate - corrfunc_estimate_old)/np.abs(corrfunc_estimate_old)
                        self.convergence_info.append(rel_diff_corrfunc)
                        corrfunc_estimate_old = corrfunc_estimate
                    else:
                        corrfunc_estimate = markovian_corrfunc_coefficient_input
                    first_timecheck = False
                np.savetxt("corrfunc_integrand_total.dat",np.concatenate((time_vec,corrfunc_integrand_data)).reshape(len(time_vec),2,order='F'))
                self.corrfunc_vec[itrx,1] = corrfunc_estimate
                self.prop_info_corrfunc[itrx,1:6] = [itrt_input,time_input,norm_input,rel_diff_input,(dt_tot_input/itrt_input)]
                print("Number of timesteps: "+str(itrt_input))
                print("Propagation time: "+str(time_input))
                print("Norm at end of propagation: "+str(norm_input))
                print("Correlation function: "+str(np.real(corrfunc_estimate)))
                print("Conv. of integral: "+str(rel_diff_input))
                print("Conv. of timestep checking: "+str(self.convergence_info[-1]))
                print("Average timestep size of propagation: ",str(dt_tot_input/itrt_input))
        else:
            self.prop_info_corrfunc,corrfunc_mol_vec_wout_x,corrfunc_molleads_vec_wout_x = \
                sparse_corrfunc.spatially_dependent_corrfunc(pair_info_row=self.pair_info_row,
                    pair_info_col=self.pair_info_col,pair_values=self.pair_values,rho_nonzeros=self.rho_nonzeros_sparse,
                    tier_index=self.tier_index,tol_corrfunc=tol_corrfunc,max_expan_order=max_expan_order,
                    dim_rho=dim_rho,nnz_elements=self.nnz_elements_sparse,npairs=self.npairs,len_un_ind=self.len_un_ind,
                    len_x_vec=len_x_vec,deriv_ham=self.deriv_ham[:,:,0,:],deriv_ham_log=self.deriv_ham_log[:,:,0,0],
                    nthreads_x_grid=nthreads_x_grid,nthreads_liouvillian=nthreads_liouvillian,
                    rk_coeff=rk_coeff,rk_coeffhat=rk_coeffhat,dt_init=dt_init,
                    facmax_init=facmax_init,facmin=facmin,fac=fac,atol_vec=self.atol_vec,rtol_vec=self.rtol_vec,
                    printing_timestep=printing_timestep,max_time=max_time,dt_min=dt_min,
                    complex_coefficients=self.complex_coefficients,min_time=min_time,
                    isreal_sparse=self.isreal_sparse,deriv_v=dV_dx,nel=Nel,nleads=Nleads,nmodes=self.nmodes,
                    nmax=Nmax,un_ind=self.un_ind,ksiglm=self.ksiglm,d_ops_log=self.d_ops_log,d_ops=self.d_ops,
                    print_integrand_yn=print_integrand_yn,ic_corrfunc=self.ic_corrfunc)
            self.corrfunc_mol_vec[:,1] = corrfunc_mol_vec_wout_x
            self.corrfunc_molleads_vec[:,1] = corrfunc_molleads_vec_wout_x
            self.corrfunc_vec[:,1] = corrfunc_mol_vec_wout_x + corrfunc_molleads_vec_wout_x
        prop_info_file = "prop_info_corrfunc.dat"
        prop_info_headings = ['Coord. Val.','No. Timesteps','Time','Convergence','Av. Timestep']
        np.savetxt(prop_info_file,self.prop_info_corrfunc,fmt='%d',delimiter='\t',
                   header=', '.join(prop_info_headings),comments='')
        if checking_timesteps_yn:
            corrfunc_file = "corrfunc.dat"
            corrfunc_headings = 'Coord. Val.\tCorr. Func.'
            np.savetxt(corrfunc_file,np.real(self.corrfunc_vec),fmt='%f',delimiter='\t',
                    header=corrfunc_headings,comments='')
        corrfunc_mol_file = "corrfunc_mol.dat"
        corrfunc_molleads_file = "corrfunc_molleads.dat"
        corrfunc_total_file = "corrfunc.dat"
        corrfunc_headings = 'Coord. Val.\tCorr. Func.'
        np.savetxt(corrfunc_mol_file,np.real(self.corrfunc_mol_vec),fmt='%.12f',delimiter='\t',header=corrfunc_headings,comments='')
        np.savetxt(corrfunc_molleads_file,np.real(self.corrfunc_molleads_vec),fmt='%.12f',delimiter='\t',header=corrfunc_headings,comments='')
        np.savetxt(corrfunc_total_file,np.real(self.corrfunc_vec),fmt='%.12f',delimiter='\t',header=corrfunc_headings,comments='')

    def corrfunction_ic(self):

        self.rho_0_diff_x_arr = np.zeros((len_x_vec,self.nnz_elements+1),dtype=float)
        self.rho_0_diff_x_arr[:,0] = x_vec
        self.dforce = np.zeros((dim_rho,dim_rho,len_x_vec),dtype=float)
        for itrx in range(len_x_vec):
            rho_ss_x_hilbert = np.zeros((dim_rho,dim_rho,self.len_un_ind),dtype=complex)
            rho_0_diff_hilbert = np.zeros((dim_rho,dim_rho,self.len_un_ind),dtype=complex)
            self.dforce[:,:,itrx] = self.deriv_ham[:,:,0] + self.average_electronic_force_vec[itrx,1]*self.identity_dim_rho
            for itrnz in range(self.nnz_elements):
                indjn = self.rho_nonzeros_sparse[itrnz,0]
                nrow = self.rho_nonzeros_sparse[itrnz,1]
                ncol = self.rho_nonzeros_sparse[itrnz,2]
                itrn = self.tier_index[indjn]
                rho_ss_x_hilbert[nrow,ncol,indjn] += self.rho_ss_x_arr[itrx,itrnz+1]*self.complex_coefficients[itrnz]
            for itr_ado in range(self.len_un_ind):
                rho_0_diff_hilbert[:,:,itr_ado] = \
                        np.matmul(self.dforce[:,:,itrx],rho_ss_x_hilbert[:,:,itr_ado]) + \
                        np.matmul(rho_ss_x_hilbert[:,:,itr_ado],self.dforce[:,:,itrx])
            for itrnz in range(self.nnz_elements):
                indjn = self.rho_nonzeros_sparse[itrnz,0]                                               
                nrow = self.rho_nonzeros_sparse[itrnz,1]                                                
                ncol = self.rho_nonzeros_sparse[itrnz,2]
                itrn = self.tier_index[indjn]
                if self.isreal_sparse[itrnz]:
                    self.rho_0_diff_x_arr[itrx,itrnz+1] = np.real(rho_0_diff_hilbert[nrow,ncol,indjn])
                else:
                    self.rho_0_diff_x_arr[itrx,itrnz+1] = np.imag(rho_0_diff_hilbert[nrow,ncol,indjn])

    def generate_ic_sparsity(self):

        if not bool(wbl_YN): 
            npairs_corrfunc,n_indnz2_this_indnz1_max_corrfunc,n_indnz2_prev_indnz1_vec_corrfunc = \
                sparsity_corrfunc.sparse_matrix_elements_a(ksiglm=self.ksiglm,tier_index=self.tier_index,
                    index_minus=self.index_minus,index_plus=self.index_plus,d_ops_log=self.d_ops_comp_log,
                    ham_log=self.ham_log,rho_sparsity=self.rho_sparsity,rho_nonzeros=self.rho_nonzeros,
                    nnz_elements=self.nnz_elements,
                    dim_rho=dim_rho,len_index_plus=self.len_index_plus,len_un_ind=self.len_un_ind,
                    len_index_minus=self.len_un_ind,nmax=Nmax,nmodes=self.nmodes,nel=Nel,
                    degenerate_levels=degenerate_levels)
            self.pair_info_row_corrfunc,self.pair_info_col_corrfunc,self.pair_values_corrfunc = \
                sparsity_corrfunc.sparse_matrix_elements_b(ksiglm=self.ksiglm,tier_index=self.tier_index,
                    index_minus=self.index_minus,index_plus=self.index_plus,d_ops=self.d_ops_comp,
                    eta_vec=self.eta_vec,rho_sparsity=self.rho_sparsity,rho_nonzeros=self.rho_nonzeros,
                    nnz_elements=self.nnz_elements,dim_rho=dim_rho,len_index_plus=self.len_index_plus,
                    len_un_ind=self.len_un_ind,len_index_minus=self.len_index_minus,
                    nmax=Nmax,nmodes=self.nmodes,nel=Nel,nsign=Nsign,npairs=npairs_corrfunc,
                    nleads=Nleads,npoles=self.npoles,ham_log=self.ham_log,d_ops_log=self.d_ops_comp_log,
                    degenerate_levels=degenerate_levels,n_indnz2_this_indnz1_max=n_indnz2_this_indnz1_max_corrfunc,
                    deriv_ham=self.deriv_ham,len_x_vec=len_x_vec,n_indnz2_prev_indnz1_vec=n_indnz2_prev_indnz1_vec_corrfunc,
                    dv_dx=dV_dx,average_electronic_force_mol_vec=self.average_electronic_force_mol_vec[:,1],
                    average_electronic_force_molleads_vec=self.average_electronic_force_molleads_vec[:,1])
        else:
            raise ValueError("WBL not yet implemented")

        npairs_corrfunc_uf = self.pair_info_row_corrfunc.shape[0]
        zero_couplings = []
        for itr_pair in range(npairs_corrfunc_uf):
            if (np.abs(self.pair_values_corrfunc[itr_pair,:]) <= 1e-14).all():
                zero_couplings.append(itr_pair)

        zero_couplings = np.array(zero_couplings,dtype=int)
        pair_values_corrfunc_uf = np.delete(self.pair_values_corrfunc,zero_couplings,axis=0)
        pair_info_row_corrfunc_uf = np.delete(self.pair_info_row_corrfunc,zero_couplings,axis=0)
        pair_info_col_corrfunc_uf = np.delete(self.pair_info_col_corrfunc,zero_couplings,axis=0)
        self.npairs_corrfunc_uf = pair_info_row_corrfunc_uf.shape[0]

        nnz_elements_sparse_uf = 2*self.nnz_elements
        self.pair_info_row_corrfunc_fil = np.copy(pair_info_row_corrfunc_uf)
        self.pair_info_col_corrfunc_fil = np.copy(pair_info_col_corrfunc_uf)
        self.pair_values_corrfunc_fil = np.copy(pair_values_corrfunc_uf)
        for itr_nnz in range(nnz_elements_sparse_uf):
            if (self.is_connected_array[itr_nnz] == False):
                pair_indices_row_above = np.where(pair_info_row_corrfunc_uf > itr_nnz)[0]
                pair_indices_col_above = np.where(pair_info_col_corrfunc_uf > itr_nnz)[0]
                self.pair_info_row_corrfunc_fil[pair_indices_row_above] -= 1
                self.pair_info_col_corrfunc_fil[pair_indices_col_above] -= 1

        pair_indices_remove = np.array([],dtype=int)
        for itr_nnz in range(nnz_elements_sparse_uf):
            if not self.is_connected_array[itr_nnz]:
                pair_indices_remove = np.concatenate((np.where(pair_info_row_corrfunc_uf == itr_nnz)[0],\
                                                        pair_indices_remove),axis=0)
                pair_indices_remove = np.concatenate((np.where(pair_info_col_corrfunc_uf == itr_nnz)[0],\
                                                        pair_indices_remove),axis=0)

        self.pair_values_corrfunc_fil = np.delete(self.pair_values_corrfunc_fil,pair_indices_remove,axis=0)
        self.pair_info_row_corrfunc_fil = np.delete(self.pair_info_row_corrfunc_fil,pair_indices_remove,axis=0)
        self.pair_info_col_corrfunc_fil = np.delete(self.pair_info_col_corrfunc_fil,pair_indices_remove,axis=0)
        self.npairs_corrfunc_fil = self.pair_info_row_corrfunc_fil.shape[0]
        if not checking_timesteps_yn:
            self.ic_corrfunc = np.zeros((len_x_vec,self.nnz_elements_sparse+1),dtype=float)
            self.ic_corrfunc[:,0] = x_vec
            for itrx in range(len_x_vec):
                ic_corrfunc_one_x = self.generate_ic_corrfunc_one_x(itrx)
                self.ic_corrfunc[itrx,1:] = ic_corrfunc_one_x

        np.save("ic_corrfunc.npy",self.ic_corrfunc)

    def generate_ic_corrfunc_one_x(self,itrx):
        
        pair_values_corrfunc_one_x = self.pair_values_corrfunc_fil[:,itrx]
        sparse_ic_corrfunc_generator_one_x = sparse.csc_matrix((pair_values_corrfunc_one_x,\
                                            (self.pair_info_row_corrfunc_fil,self.pair_info_col_corrfunc_fil)),
                                            shape=(self.nnz_elements_sparse,self.nnz_elements_sparse),dtype=float)
        ic_corrfunc_one_x = sparse_ic_corrfunc_generator_one_x.dot(self.rho_ss_x_arr[itrx,1:])

        return ic_corrfunc_one_x

    def prony_function(self,t, F, m):
        # Solve LLS problem in step 1
        # Amat is (N-m)*m and bmat is N-m*1
        N = len(t)
        Amat = np.zeros((N-m, m),dtype=complex)
        bmat = F[m:N]
        for jcol in range(m):
            Amat[:, jcol] = F[m-jcol-1:N-1-jcol]
        sol = np.linalg.lstsq(Amat, bmat,rcond=None)
        d = sol[0]
        # Solve the roots of the polynomial in step 2
        # first, form the polynomial coefficients
        c = np.zeros(m+1,dtype=complex)
        c[m] = 1.
        for i in range(1,m+1):
            c[m-i] = -d[i-1]
        u = poly.polyroots(c)
        b_est = np.log(u)/(t[1] - t[0])
        # Set up LLS problem to find the "a"s in step 3
        Amat = np.zeros((N, m),dtype=complex)
        bmat = F
        for irow in range(N):
            Amat[irow, :] = u**irow
        sol = np.linalg.lstsq(Amat, bmat,rcond=None)
        a_est = sol[0]
        return a_est, b_est
