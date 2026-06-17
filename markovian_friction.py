import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy
from scipy.interpolate import CubicSpline
import sparse_friction
from input_parameters import *
import math

class calculate_markovian_friction():

    def __init__(self,pair_info_col,pair_info_row,pair_values,npairs,nnz_elements_sparse,
                    sparse_trace_array,nnz_elements_sparse_zeroth_tier,
                    complex_coefficients,isreal_sparse,trace_cols,
                    rho_nonzeros_sparse,tier_index,rho_ss_spatial_derivative,deriv_ham,deriv_ham_log,len_un_ind,
                    atol_vec,rtol_vec,n_prony_terms,nmodes,un_ind,ksiglm,d_ops_log,d_ops):

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
        self.rho_ss_spatial_derivative = rho_ss_spatial_derivative
        self.nnz_elements = nnz_elements_sparse
        self.len_un_ind = len_un_ind
        self.deriv_ham = deriv_ham
        self.deriv_ham_log = deriv_ham_log
        self.atol_vec = atol_vec
        self.rtol_vec = rtol_vec
        self.n_prony_terms = n_prony_terms
        self.nmodes = nmodes
        self.un_ind = un_ind
        self.ksiglm = ksiglm
        self.d_ops_log = d_ops_log
        self.d_ops = d_ops

        if parallelize_x_grid:
            self.generate_markovian_friction_parallel()
        else:
            self.generate_markovian_friction_sequential()

    def generate_markovian_friction_parallel(self):

        print("Not yet implemented")

    def generate_markovian_friction_sequential(self):

        self.friction_vec = np.zeros((len_x_vec,2),dtype=complex)
        self.friction_vec[:,0] = x_vec
        self.friction_mol_vec = np.zeros((len_x_vec,2),dtype=complex)
        self.friction_mol_vec[:,0] = x_vec
        self.friction_molleads_vec = np.zeros((len_x_vec,2),dtype=complex)
        self.friction_molleads_vec[:,0] = x_vec
        self.current_na_vec = np.zeros((len_x_vec,Nleads+1),dtype=complex)
        self.current_na_vec[:,0] = x_vec
        if checking_timesteps_yn:
            if min_time > checking_time_interval:
                raise ValueError("Minimum propagation time cannot be smaller than the time interval between \
                                successive convergence checks")
            self.prop_info_friction = np.zeros((len_x_vec,6),dtype=float)
            self.prop_info_friction[:,0] = x_vec
            for itrx in range(len_x_vec):
                print("Calculating Markovian friction coefficient of "+str(itrx)+"th coordinate from "+str(len_x_vec))
                pair_values_x = self.pair_values[:,itrx,0]
                deriv_ham_one_x = self.deriv_ham[:,:,0,itrx]
                deriv_V_one_x = dV_dx[:,:,itrx]
                rel_diff_friction = 1.0
                ic_prop = self.rho_ss_spatial_derivative[itrx,1:]
                dt_input = dt_init
                facmax_input = facmax_init
                dt_tot_input = 0.0
                rel_diff_input = 1.0
                markovian_friction_coefficient_mol_input = 0.0
                markovian_friction_coefficient_molleads_input = 0.0
                norm_input = 0.0
                time_input = 0.0
                itrt_input = int(0)
                friction_estimate_old = 1.0
                self.convergence_info = []
                time_vec = []
                friction_integrand_data_mol = []
                friction_integrand_data_molleads = []
                continue_checking = True
                first_timecheck = True
                while ((rel_diff_friction >= tol_checking_time) and continue_checking):
                    propagation_outputs = \
                        sparse_friction.sparse_markovian_friction_propagation_w_checking(pair_info_row=self.pair_info_row,
                            pair_info_col=self.pair_info_col,pair_values=pair_values_x,rho_nonzeros=self.rho_nonzeros_sparse,
                            tier_index=self.tier_index,rho_0=ic_prop,tol_friction=tol_friction,
                            deriv_ham_one_x=deriv_ham_one_x,deriv_ham_log=self.deriv_ham_log[:,:,0,0],
                            max_expan_order=max_expan_order,dim_rho=dim_rho,nnz_elements=self.nnz_elements_sparse,
                            npairs=self.npairs,len_un_ind=self.len_un_ind,nthreads_liouvillian=nthreads_liouvillian,
                            rk_coeff=rk_coeff,rk_coeffhat=rk_coeffhat,
                            facmin=facmin,fac=fac,atol_vec=self.atol_vec,rtol_vec=self.rtol_vec,
                            max_time=max_time,dt_min=dt_min,complex_coefficients=self.complex_coefficients,min_time=min_time,
                            time_input=time_input,dt_tot_input=dt_tot_input,norm_input=norm_input,
                            rel_diff_input=rel_diff_input,itrt_input=itrt_input,dt_input=dt_input,facmax_input=facmax_input,
                            markovian_friction_coefficient_mol_input=markovian_friction_coefficient_mol_input,
                            markovian_friction_coefficient_molleads_input=markovian_friction_coefficient_molleads_input,
                            checking_time_interval=checking_time_interval,first_timecheck=first_timecheck,
                            isreal_sparse=self.isreal_sparse,deriv_v_one_x=deriv_V_one_x,nel=Nel,nleads=Nleads,
                            nmodes=self.nmodes,nmax=Nmax,un_ind=self.un_ind,ksiglm=self.ksiglm,
                            d_ops_log=self.d_ops_log,d_ops=self.d_ops)
                    continue_checking,ic_prop,itrt_input,norm_input,rel_diff_input,\
                    markovian_friction_coefficient_mol_input,markovian_friction_coefficient_molleads_input,\
                    time_input,dt_tot_input,dt_input = propagation_outputs
                    if (continue_checking):
                        new_integrand_data_mol = np.loadtxt('friction_integrand_mol_new.dat')
                        new_integrand_data_molleads = np.loadtxt('friction_integrand_molleads_new.dat')
                        time_vec.extend(new_integrand_data_mol[:,0].tolist())
                        friction_integrand_data_mol.extend(new_integrand_data_mol[:,1].tolist())
                        friction_integrand_data_molleads.extend(new_integrand_data_molleads[:,1].tolist())
                        integrand_func_mol = CubicSpline(time_vec,friction_integrand_data_mol)
                        integrand_func_molleads = CubicSpline(time_vec,friction_integrand_data_molleads)
                        time_vec_even = np.linspace(time_vec[0],time_vec[-1],2*len(time_vec))
                        friction_integrand_even_mol = integrand_func_mol(time_vec_even)
                        friction_integrand_even_molleads = integrand_func_molleads(time_vec_even)
                        weights_mol,frequencies_mol = self.prony_function(time_vec_even,friction_integrand_even_mol,self.n_prony_terms)
                        weights_molleads,frequencies_molleads = self.prony_function(time_vec_even,friction_integrand_even_molleads,self.n_prony_terms)
                        friction_estimate_mol = np.real(np.sum(weights_mol/frequencies_mol))
                        friction_estimate_molleads = np.real(np.sum(weights_molleads/frequencies_molleads))
                        if math.isnan(friction_estimate_molleads):
                            friction_estimate_molleads = 0.0
                        friction_estimate = friction_estimate_mol + friction_estimate_molleads
                        rel_diff_friction = np.abs(friction_estimate - friction_estimate_old)/np.abs(friction_estimate_old)
                        print("Rel. Diff at time",time_vec[-1]," is ",rel_diff_friction)
                        self.convergence_info.append(rel_diff_friction)
                        friction_estimate_old = friction_estimate
                    else:
                        friction_mol_estimate = markovian_friction_coefficient_mol_input
                        friction_molleads_estimate = markovian_friction_coefficient_molleads_input
                        friction_estimate = markovian_friction_coefficient_mol_input + markovian_friction_coefficient_molleads_input
                    first_timecheck = False
                np.savetxt("friction_integrand_mol_total.dat",\
                           np.concatenate((time_vec,friction_integrand_data_mol)).reshape(len(time_vec),2,order='F'))
                np.savetxt("friction_integrand_molleads_total.dat",\
                           np.concatenate((time_vec,friction_integrand_data_molleads)).reshape(len(time_vec),2,order='F'))
                np.savetxt("friction_integrand_total.dat",\
                           np.concatenate((time_vec,[x + y for x, y in zip(friction_integrand_data_mol,\
                                            friction_integrand_data_molleads)])).reshape(len(time_vec),2,order='F'))
                self.friction_mol_vec[itrx,1] = friction_estimate_mol
                self.friction_molleads_vec[itrx,1] = friction_estimate_molleads
                self.friction_vec[itrx,1] = friction_estimate
                self.prop_info_friction[itrx,1:6] = [itrt_input,time_input,norm_input,rel_diff_input,(dt_tot_input/itrt_input)]
                print("Number of timesteps: "+str(itrt_input))
                print("Propagation time: "+str(time_input))
                print("Norm at end of propagation: "+str(norm_input))
                print("Friction: "+str(np.real(friction_estimate)))
                print("Conv. of integral: "+str(rel_diff_input))
                print("Conv. of timestep checking: "+str(self.convergence_info[-1]))
                print("Average timestep size of propagation: ",str(dt_tot_input/itrt_input))
        else:
            self.prop_info_friction,friction_mol_vec_wout_x,friction_molleads_vec_wout_x,current_na_vec_wout_x = \
                        sparse_friction.spatially_dependent_friction(pair_info_row=self.pair_info_row,
                            pair_info_col=self.pair_info_col,pair_values=self.pair_values,rho_nonzeros=self.rho_nonzeros_sparse,
                            tier_index=self.tier_index,x_vec=x_vec,rho_ss_spatial_derivative=self.rho_ss_spatial_derivative,
                            tol_friction=tol_friction,max_expan_order=max_expan_order,
                            dim_rho=dim_rho,nnz_elements=self.nnz_elements_sparse,npairs=self.npairs,len_un_ind=self.len_un_ind,
                            len_x_vec=len_x_vec,deriv_ham=self.deriv_ham[:,:,0,:],deriv_ham_log=self.deriv_ham_log[:,:,0,0],
                            nthreads_x_grid=nthreads_x_grid,nthreads_liouvillian=nthreads_liouvillian,
                            rk_coeff=rk_coeff,rk_coeffhat=rk_coeffhat,dt_init=dt_init,
                            facmax_init=facmax_init,facmin=facmin,fac=fac,atol_vec=self.atol_vec,rtol_vec=self.rtol_vec,
                            printing_timestep=printing_timestep,max_time=max_time,dt_min=dt_min,
                            complex_coefficients=self.complex_coefficients,min_time=min_time,isreal_sparse=self.isreal_sparse,
                            deriv_v=dV_dx,nel=Nel,nleads=Nleads,nmodes=self.nmodes,nmax=Nmax,un_ind=self.un_ind,
                            ksiglm=self.ksiglm,d_ops_log=self.d_ops_log,d_ops=self.d_ops,print_integrand_yn=print_integrand_yn,
                            degenerate_levels=degenerate_levels,el_lead_couplings=el_lead_couplings[:,:,:,0])
            self.friction_mol_vec[:,1] = friction_mol_vec_wout_x
            self.friction_molleads_vec[:,1] = friction_molleads_vec_wout_x
            self.friction_vec[:,1] = friction_mol_vec_wout_x + friction_molleads_vec_wout_x
            self.current_na_vec[:,1:] = current_na_vec_wout_x
        prop_info_file = "prop_info_friction.dat"
        prop_info_headings = ['Coord. Val.','No. Timesteps','Time','Convergence','Av. Timestep']
        np.savetxt(prop_info_file,self.prop_info_friction,fmt='%d',delimiter='\t',
                   header=', '.join(prop_info_headings),comments='')
        friction_mol_file = "friction_mol.dat"
        friction_molleads_file = "friction_molleads.dat"
        friction_total_file = "friction.dat"
        friction_headings = 'Coord. Val.\tFriction'
        current_na_file = "current_na.dat"
        current_na_headings = 'Coord. Val.\tCurrent_0\tCurrent_1'
        np.savetxt(friction_mol_file,np.real(self.friction_mol_vec),fmt='%.12f',delimiter='\t',header=friction_headings,comments='')
        np.savetxt(friction_molleads_file,np.real(self.friction_molleads_vec),fmt='%.12f',delimiter='\t',header=friction_headings,comments='')
        np.savetxt(friction_total_file,np.real(self.friction_vec),fmt='%.12f',delimiter='\t',header=friction_headings,comments='')
        np.savetxt(current_na_file,np.real(self.current_na_vec),fmt='%.12f',delimiter='\t',header=current_na_headings,comments='')

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
