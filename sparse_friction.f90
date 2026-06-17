! ---------------------------------------------------------------------------
! 
!                     PERFORM TIME-PROPAGATION OF ADOs
!
! ---------------------------------------------------------------------------
!
! This Fortran subroutine is meant to be converted to a Python wrapper using f2py, for use 
! in the main Python code. The various subroutines detail the sparse nature of the HQME.
!
! USAGE - RUN FROM TERMINAL (LINUX) TO CREATE PYTHON WRAPPER:
!       Check fortran compilers available in your platform:  f2py -c --help-fcompiler
!       For 'ifort': f2py -c -m sparse_propagation sparse_propagation.f90 --opt='-O3' --fcompiler=intelem --f90flags='-openmp -D__OPENMP' -liomp5
!       For 'gfortran': f2py -c -m sparse_propagation sparse_propagation.f90 --opt='-O3' --fcompiler=gnu95 --f90flags='-fopenmp -D__OPENMP' -lgomp
!
! It uses parallel programming and the previously generated sparse representation of the HQME to 
! propagate an initial density matrix in time, keeping track of both the time-dependent current
! and the density matrix itself.
!
! USAGE - RUN FROM MAIN PYTHON CODE ONCE WRAPPED:
!       rho_system,current = sparse_propagation.sparse_propagation(pair_info=pair_info,pair_values=pair_values,ham_pair_info=ham_pair_info,
!                               rho_nonzeros=rho_nonzeros,gamma_values=gamma_values,ham=Ham,tier_index=tier_index,un_ind=Un_Ind,d_ops=d_ops,
!                               ksiglm=KsigLm,rho_0=rho_0,dt=dt,nsteps=nsteps,dim_rho=dim_rho,max_expan_order=max_expan_order,nnz_elements=nnz_elements,
!                               npairs=npairs,nhampairs=nhampairs,len_un_ind=len_un_ind,nmax=Nmax,nmodes=Nmodes,nel=Nel,nleads=Nleads)
!
! INPUTS:
!       pair_info                                   Array of size [npairs,4] containing information about the pairs of coupled nonzero elements in the HQME.
!                                                   Each row corresponds to a different coupled pair;
!                                                   Column 1 contains the nonzero index of the LHS ADO element, 
!                                                   Column 2 contains the nonzero index of the RHS ADO element,
!                                                   Column 3 contains 1 if the connection between these two ADOs requires a hermiticity relation, and 0 if not
!                                                   Column 4 contains ??? FINISH
!       pair_values                                 Array of size [1,npairs] containing the value of this connection (i.e. the corresponding element of e^{L*dt} in Liouville space)            
!
!       ham_pair_info                               The same as pair_info, but just for the coherent part containing the commutator with the Hamiltonian
!
!       rho_nonzeros                                Array of size [nnz_elements,3] containing information about the nonzero elements in the HQME
!                                                   Each row corresponds to a different nonzero element and the columns contain the index of that ADO
!                                                   to which it belongs, and its row and and column, in that order.
!
!       gamma_values                                Array of size [1,nnz_elements] containing the sum over gamma values for each ADO in the HQME                                
! 
!       ham                                         Array of size [dim_rho,dim_rho] containing system Hamiltonian
!
!       tier_index                                  Array of size [1,len_un_ind] containing tier of each unique ADO
!
!       un_ind                                      Array of size [len_un_ind,nmax] containing the modes of each unique ADO    
!      
!       d_ops                                       Array of size [dim_rho,dim_rho,nel,2] containing annihilation and creation operators for each electronic
!                                                   level in the system
!
!       ksiglm                                      Array of size [1,nmodes] containing all modes j = {K,sigma,l,m}
!
!       rho_0                                       Array of size [dim_rho,dim_rho] containing initial density matrix of system
!
!       dt                                          Scalar that specifies size of each timestep
!       
!       nsteps                                      Scalar that specifies the number of timesteps to propagate
!   
!       dim_rho                                     Scalar that specifies the number of Fock states defining the system
!
!       max_expan_order                             Scalar that determines the maximum order to go to in the aprpoximate expansion of e^{L*dt}
!
!       nnz_elements                                Number of ADO elements that are required to be propagated
!       
!       npairs                                      Number of connected pairs of elements in HQME. Say we have an ADO of index i. If we consider the element of rho_{i} in the 
!                                                   ath row and bth column, rho_{i}(a,b), then its time-evolution could depend on the cth row and dth column of ADO with index j: rho_{j}(c,d).
!                                                   This connection may come from any part of the dissipative part of the HQME. The total number of these connected pairs is npairs.
!
!       nhampairs                                   The same, but just for the coherent part of the HQME: the commutator [H,rho]
!
!       len_un_ind                                  Scalar determining the number of unique ADOs required in the HQME
!
!       nmax                                        Scalar giving the maximum tier of the hierarchy
!
!       nmodes                                      Scalar giving the number of modes, j = {K,sigma,l,m}, included in HQME
!
!       nel                                         Scalar giving the number of electronic levels in the system
!
!       nleads                                      calar giving the number of electrodes
!
! OUTPUTS:
!       rho_all_ADOs_1_nuc_frame
!

include 'mkl_spblas.f90'

subroutine spatially_dependent_friction(pair_info_row,pair_info_col,pair_values,rho_nonzeros,tier_index,&
            x_vec,rho_ss_spatial_derivative,tol_friction,max_expan_order,dim_rho,nnz_elements,npairs,&
            len_un_ind,len_x_vec,deriv_ham,deriv_ham_log,nthreads_x_grid,nthreads_liouvillian,&
            prop_info_friction,rk_coeff,rk_coeffhat,dt_init,facmax_init,facmin,fac,atol_vec,rtol_vec,&
            printing_timestep,max_time,dt_min,complex_coefficients,min_time,isreal_sparse,deriv_v,&
            nel,nleads,nmodes,nmax,un_ind,ksiglm,d_ops_log,d_ops,markovian_friction_coefficient_mol_vec,&
            markovian_friction_coefficient_molleads_vec,print_integrand_yn,current_na_vec,degenerate_levels,&
            el_lead_couplings)
                                        
    implicit none
    
    include 'omp_lib.h'                                                                             ! Import the OpenMP (OMP) library for parallel programming
    integer, intent(in) :: dim_rho,max_expan_order,nnz_elements,npairs,len_un_ind
    integer, intent(in) :: len_x_vec,nthreads_x_grid,nthreads_liouvillian,nel,nleads,nmodes,nmax
    double precision, intent(in), dimension(0:len_x_vec-1) :: x_vec
    double precision, intent(in), dimension(0:npairs-1,0:len_x_vec-1,0:4) :: pair_values
    integer, intent(in), dimension(0:npairs-1) :: pair_info_row,pair_info_col
    integer, intent(in), dimension(0:nnz_elements-1,0:2) :: rho_nonzeros
    double complex, intent(in), dimension(0:nnz_elements-1) :: complex_coefficients
    logical, intent(in), dimension(0:nnz_elements-1) :: isreal_sparse
    integer, intent(in), dimension(0:len_un_ind-1) :: tier_index
    integer, intent(in), dimension(0:len_un_ind-1,0:nmax-1) :: un_ind
    integer, intent(in), dimension(0:nmodes-1,0:3) :: ksiglm
    double complex, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:len_x_vec-1) :: deriv_ham
    double precision, intent(in), dimension(0:nleads-1,0:nel-1,0:len_x_vec-1) :: deriv_v
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: deriv_ham_log         
    double precision, intent(in), dimension(0:len_x_vec-1,0:nnz_elements) :: rho_ss_spatial_derivative
    double precision, intent(in), dimension(0:max_expan_order-1) :: rk_coeff,rk_coeffhat
    double precision, intent(in) :: dt_init,facmax_init,facmin,fac,tol_friction,printing_timestep,max_time,dt_min,min_time
    double precision, intent(in), dimension(0:nnz_elements-1) :: atol_vec,rtol_vec
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops_log
    double precision, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops
    logical, intent(in) :: print_integrand_yn,degenerate_levels
    double precision, dimension(0:nleads-1,0:nel-1,0:len_x_vec-1) :: el_lead_couplings

    double precision, intent(out), dimension(0:len_x_vec-1):: markovian_friction_coefficient_mol_vec
    double precision, intent(out), dimension(0:len_x_vec-1):: markovian_friction_coefficient_molleads_vec
    double precision, intent(out), dimension(0:len_x_vec-1,0:nleads-1) :: current_na_vec
    double precision, intent(out), dimension(0:len_x_vec-1,0:4) :: prop_info_friction

    double precision, dimension(0:nnz_elements-1) :: rho_0
    integer :: itrx,itrt
    double complex :: markovian_friction_coefficient_mol,markovian_friction_coefficient_molleads
    double precision :: norm,time,rel_diff,dt_av
    double precision, dimension(0:npairs-1) :: pair_values_one_x
    double precision, dimension(0:nleads-1,0:nel-1) :: el_lead_couplings_x
    double complex, dimension(0:dim_rho-1,0:dim_rho-1) :: deriv_ham_one_x
    double precision, dimension(0:nleads-1,0:nel-1) :: deriv_v_one_x
    double precision, dimension(0:nleads-1) :: current_na

    markovian_friction_coefficient_mol_vec = 0.d0
    markovian_friction_coefficient_molleads_vec = 0.d0
    current_na_vec = 0.d0
    prop_info_friction = 0.d0
    call omp_set_num_threads(nthreads_x_grid)
    !$OMP PARALLEL SHARED(markovian_friction_coefficient_mol_vec,markovian_friction_coefficient_molleads_vec,prop_info_friction)
    !$OMP DO PRIVATE(itrx,pair_values_one_x,markovian_friction_coefficient_mol,&
    !$OMP            markovian_friction_coefficient_molleads,rho_0,itrt,norm,time,rel_diff,&
    !$OMP            deriv_ham_one_x,deriv_v_one_x,current_na,el_lead_couplings_x)
    do itrx = 0,(len_x_vec-1)
        print *,"Doing ",itrx,"th value of ",len_x_vec," x values"
        rho_0 = rho_ss_spatial_derivative(itrx,1:)
        el_lead_couplings_x = el_lead_couplings(:,:,itrx)
        pair_values_one_x = pair_values(:,itrx,0)
        deriv_ham_one_x = deriv_ham(:,:,itrx)
        deriv_v_one_x = deriv_v(:,:,itrx)
        call sparse_markovian_friction_propagation(pair_info_row,pair_info_col,pair_values_one_x,rho_nonzeros,tier_index,&
                            rho_0,tol_friction,deriv_ham_one_x,deriv_ham_log,markovian_friction_coefficient_mol,&
                            markovian_friction_coefficient_molleads,max_expan_order,dim_rho,nnz_elements,npairs,&
                            len_un_ind,nthreads_liouvillian,itrt,norm,rk_coeff,rk_coeffhat,dt_init,facmax_init,&
                            facmin,fac,atol_vec,rtol_vec,rel_diff,time,dt_av,printing_timestep,max_time,dt_min,&
                            complex_coefficients,min_time,isreal_sparse,deriv_v_one_x,nel,nleads,nmodes,nmax,un_ind,&
                            ksiglm,d_ops_log,d_ops,print_integrand_yn,current_na,degenerate_levels,el_lead_couplings_x)
        prop_info_friction(itrx,0) = float(itrt)
        prop_info_friction(itrx,1) = norm
        prop_info_friction(itrx,2) = time
        prop_info_friction(itrx,3) = rel_diff
        prop_info_friction(itrx,4) = dt_av
        markovian_friction_coefficient_mol_vec(itrx) = dble(markovian_friction_coefficient_mol)
        markovian_friction_coefficient_molleads_vec(itrx) = dble(markovian_friction_coefficient_molleads)
        current_na_vec(itrx,:) = current_na*(6.623618237510/27.211386245988)*(10**3)
    enddo
    !$OMP END DO 
    !$OMP END PARALLEL

    ! open(100,file='friction_mol.dat')                                                                   ! and information about time the program takes to run
    ! open(200,file='friction_molleads.dat')                                                                   ! and information about time the program takes to run
    ! do itrx = 0,(len_x_vec-1)
    !     write(100,'(f15.5,a,f15.10)') x_vec(itrx)," ",dble(markovian_friction_coefficient_mol_vec(itrx))
    !     write(200,'(f15.5,a,f15.10)') x_vec(itrx)," ",dble(markovian_friction_coefficient_molleads_vec(itrx))
    ! enddo
    ! close(100)
    ! close(200)

end subroutine spatially_dependent_friction

subroutine sparse_markovian_friction_propagation(pair_info_row,pair_info_col,pair_values_one_x,rho_nonzeros,tier_index,&
            rho_0,tol_friction,deriv_ham_one_x,deriv_ham_log,markovian_friction_coefficient_mol,&
            markovian_friction_coefficient_molleads,&
            max_expan_order,dim_rho,nnz_elements,npairs,len_un_ind,nthreads_liouvillian,itrt,&
            norm,rk_coeff,rk_coeffhat,dt_init,facmax_init,facmin,fac,atol_vec,&
            rtol_vec,rel_diff,time,dt_av,printing_timestep,max_time,dt_min,&
            complex_coefficients,min_time,isreal_sparse,deriv_v_one_x,nel,nleads,nmodes,nmax,un_ind,&
            ksiglm,d_ops_log,d_ops,print_integrand_yn,current_na,degenerate_levels,el_lead_couplings_x)

    use mkl_spblas

    implicit none

    integer, intent(in) :: dim_rho,max_expan_order,nnz_elements,npairs,nthreads_liouvillian
    integer, intent(in) :: len_un_ind,nel,nleads,nmodes,nmax
    double complex, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: deriv_ham_one_x
    double precision, intent(in), dimension(0:nleads-1,0:nel-1) :: deriv_v_one_x
    double precision, intent(in), dimension(0:nnz_elements-1) :: rho_0
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: deriv_ham_log
    double precision, intent(in), dimension(0:npairs-1) :: pair_values_one_x
    integer, intent(in), dimension(0:npairs-1) :: pair_info_row,pair_info_col
    integer, intent(in), dimension(0:nnz_elements-1,0:2) :: rho_nonzeros
    double complex, intent(in), dimension(0:nnz_elements-1) :: complex_coefficients
    logical, intent(in), dimension(0:nnz_elements-1) :: isreal_sparse
    integer, intent(in), dimension(0:len_un_ind-1) :: tier_index
    integer, intent(in), dimension(0:len_un_ind-1,0:nmax-1) :: un_ind
    integer, intent(in), dimension(0:nmodes-1,0:3) :: ksiglm
    double precision, intent(in), dimension(0:max_expan_order-1) :: rk_coeff,rk_coeffhat
    double precision, intent(in) :: dt_init,facmax_init,facmin,fac,tol_friction,printing_timestep,max_time,dt_min,min_time
    double precision, intent(in), dimension(0:nnz_elements-1) :: atol_vec,rtol_vec
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops_log
    double precision, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops
    logical, intent(in) :: print_integrand_yn,degenerate_levels
    double precision, dimension(0:nleads-1,0:nel-1) :: el_lead_couplings_x

    double complex, intent(out) :: markovian_friction_coefficient_mol,markovian_friction_coefficient_molleads
    double precision, dimension(0:nleads-1), intent(out) :: current_na
    double precision, intent(out) :: norm,rel_diff,time,dt_av
    integer, intent(out) :: itrt

    double precision, dimension(0:nnz_elements-1) :: rho,rho_old,rho_temp,rho_deriv,sc,rho_diff_adaptive
    double precision, dimension(0:nnz_elements-1) :: rho_new,rhohat_new,abs_rho,abs_rho_new,fac_a
    double complex :: trace_friction_mol,trace_friction_molleads,trace_friction_mol_old
    double complex :: trace_friction_molleads_old,markovian_friction_coefficient_total
    double complex :: markovian_friction_coefficient_total_old
    double precision, dimension(0:nleads-1) :: trace_current_na,trace_current_na_old,current_na_old
    integer :: indjn,itrnz,nrow,ncol,itrl,jn,leads_n,sign_n,eldash_n
    integer :: indnz1,indnz2,itrn,count_adaptive_timestep,ik_coo,ik_csr,itrel
    double precision :: err_adaptive_timestep,facmax,dt,dt_tot,dt_new,dnrm2
    double precision :: dt_old,time_of_last_print,time_since_last_print
    double precision, dimension(0:nleads-1) :: rel_diff_current_na
    logical :: already_t_min

    TYPE(SPARSE_MATRIX_T) :: sparse_handle_coo,sparse_handle_csr
    TYPE(MATRIX_DESCR) :: descra
    descra % TYPE = SPARSE_MATRIX_TYPE_GENERAL

    rho = rho_0                                                                                     ! Create array containing initial values of all nonzero ADO elements. 
                                                                                                    ! rho is a temporary array that will be updated each timestep to contain all nonzero 
                                                                                                    ! ADO elements at that timstep
    trace_friction_mol = 0.d0
    trace_friction_molleads = 0.d0
    trace_current_na = 0.d0
    do itrnz = 0,nnz_elements-1                                                                     ! Loop through nonzero elements (elements that must be propagated)
        indjn = rho_nonzeros(itrnz,0)                                                               ! For this nonzero element, find the ADO index and row/column value from rho_nonzero
        nrow = rho_nonzeros(itrnz,1)                                                                ! So find the row and column of this ADO element
        ncol = rho_nonzeros(itrnz,2)
        itrn = tier_index(indjn)
        if ((itrn == 0) .and. (deriv_ham_log(ncol,nrow) .eqv. .true.)) then
            if (isreal_sparse(itrnz) .eqv. .true.) then
                trace_friction_mol = trace_friction_mol + deriv_ham_one_x(ncol,nrow)*rho(itrnz)
            endif
        elseif (itrn == 1) then
            jn = un_ind(indjn,0)                                                                ! Find mode of this 1st tier ADO
            leads_n = ksiglm(jn,0)                                                              ! Find lead index of this mode
            sign_n = 1-ksiglm(jn,1)                                                             ! Find \bar{\sigma} (opposite sign) of this mode
            eldash_n = ksiglm(jn,3)
            if (d_ops_log(ncol,nrow,eldash_n,sign_n) .eqv. .true.) then
                if (isreal_sparse(itrnz) .eqv. .true.) then
                    trace_friction_molleads = trace_friction_molleads + deriv_v_one_x(leads_n,eldash_n)*2.d0*&
                                    dble(d_ops(ncol,nrow,eldash_n,sign_n))*rho(itrnz)
                else
                    if (degenerate_levels .eqv. .false.) then
                        do itrel = 0,nel-1
                            trace_current_na(leads_n) = trace_current_na(leads_n) + &
                                    ((-1.0)**sign_n)*2.d0*d_ops(ncol,nrow,itrel,sign_n)*rho(itrnz)*el_lead_couplings_x(leads_n,itrel)
                        enddo
                    else
                        trace_current_na(leads_n) = trace_current_na(leads_n) + &
                                ((-1.0)**sign_n)*2.d0*d_ops(ncol,nrow,eldash_n,sign_n)*rho(itrnz)*el_lead_couplings_x(leads_n,eldash_n)
                    endif
                endif
            endif
        elseif (itrn > 1) then
            exit
        endif
    enddo

    markovian_friction_coefficient_mol = 0.d0
    markovian_friction_coefficient_molleads = 0.d0
    markovian_friction_coefficient_total = 0.d0
    current_na = 0.d0
    time = 0.d0
    if (print_integrand_yn .eqv. .true.) then
        open(300,file='friction_integrand_heom.dat')
        open(400,file='current_na_integrand_heom.dat')
        write(300,'(f15.5,a,f15.10)') time," ",dble(trace_friction_mol + trace_friction_molleads)
        write(400,'(f15.5,a,f15.10)') time," ",dble(trace_current_na(0))
    endif
    rel_diff=1.d0
    itrt = 0
    facmax = facmax_init
    dt = dt_init
    dt_tot = 0.d0
    norm = 0.d0
    time_of_last_print = 0.d0
    time_since_last_print = 0.d0

    ik_coo = mkl_sparse_d_create_coo(sparse_handle_coo,SPARSE_INDEX_BASE_ZERO,nnz_elements,nnz_elements,&
                                    npairs,pair_info_row,pair_info_col,pair_values_one_x)
    ik_csr = mkl_sparse_convert_csr(sparse_handle_coo,SPARSE_OPERATION_NON_TRANSPOSE,sparse_handle_csr)
    if (ik_csr.ne.0) write(*,*) 'Create handle for sparse matrix is wrong', ik_csr
    ik_csr = MKL_SPARSE_OPTIMIZE( sparse_handle_csr )
    if (ik_csr.ne.0) write(*,*) 'OPTIMIZE wrong', ik_csr
    ik_csr = mkl_sparse_copy(sparse_handle_csr,descra,sparse_handle_csr)
    call mkl_set_num_threads(nthreads_liouvillian)
    call omp_set_dynamic(0)
    call mkl_set_dynamic(0)
    already_t_min = .false.
    do while (((rel_diff >= tol_friction) .and. (time <= max_time)) .or. (time <= min_time))
        err_adaptive_timestep = 2.d0
        count_adaptive_timestep = 0
        do while (err_adaptive_timestep > 1.d0)
            call dcopy(nnz_elements,rho,1,rho_temp,1)
            call dcopy(nnz_elements,rho,1,rho_new,1)
            call dcopy(nnz_elements,rho,1,rhohat_new,1)
            do itrl = 0,max_expan_order-1                                                               ! Loop through the Taylor series expansion of e^{L*dt}
                rho_deriv = 0.d0
                ik_csr = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,dt,sparse_handle_csr,descra,rho_temp,1.0d0,rho_deriv)
                call daxpy(nnz_elements,rk_coeff(itrl),rho_deriv,1,rho_new,1)
                call daxpy(nnz_elements,rk_coeffhat(itrl),rho_deriv,1,rhohat_new,1)
                call dcopy(nnz_elements,rho_deriv,1,rho_temp,1)
            enddo                  
            call vdAbs(nnz_elements,rho,abs_rho)
            call vdAbs(nnz_elements,rho_new,abs_rho_new)
            call dcopy(nnz_elements,atol_vec,1,sc,1)
            call vdFmax (nnz_elements,abs_rho,abs_rho_new,fac_a)
            call vdmul(nnz_elements,fac_a,rtol_vec,fac_a)
            call daxpy(nnz_elements,1.d0,fac_a,1,sc,1)
            call dcopy(nnz_elements,rho_new,1,rho_diff_adaptive,1)
            call daxpy(nnz_elements,-1.d0,rhohat_new,1,rho_diff_adaptive,1)
            call vddiv(nnz_elements,rho_diff_adaptive,sc,rho_diff_adaptive)
            err_adaptive_timestep = dnrm2(nnz_elements,rho_diff_adaptive,1)/sqrt(dble(nnz_elements))
            dt_new = dt*min(facmax,max(facmin,fac*((1/err_adaptive_timestep)**(1.d0/5.d0))))
            if (err_adaptive_timestep > 1.d0) then
                facmax = 1.d0
                if (dt_new < dt_min) then
                    if (already_t_min .neqv. .true.) then
                        dt_old = dt_min
                        dt = dt_min
                        call dcopy(nnz_elements,rho,1,rho_temp,1)
                        call dcopy(nnz_elements,rho,1,rho_new,1)
                        call dcopy(nnz_elements,rho,1,rhohat_new,1)
                        do itrl = 0,max_expan_order-1                                                               ! Loop through the Taylor series expansion of e^{L*dt}
                            rho_deriv = 0.d0
                            ik_csr = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,dt,sparse_handle_csr,descra,rho_temp,1.0d0,rho_deriv)
                            call daxpy(nnz_elements,rk_coeff(itrl),rho_deriv,1,rho_new,1)
                            call daxpy(nnz_elements,rk_coeffhat(itrl),rho_deriv,1,rhohat_new,1)
                            call dcopy(nnz_elements,rho_deriv,1,rho_temp,1)
                        enddo
                        already_t_min = .true.
                        time = time + dt
                        err_adaptive_timestep = 0.d0
                        count_adaptive_timestep = count_adaptive_timestep + 1
                    else
                        dt_old = dt_min
                        dt = dt_min
                        time = time + dt
                        err_adaptive_timestep = 0.d0
                    endif
                else
                    dt_old = dt
                    dt = dt_new
                    already_t_min = .false.
                endif
            else
                dt_old = dt
                time = time + dt
                if (dt_new < dt_min) then
                    dt = dt_min
                    already_t_min = .true.
                else
                    dt = dt_new
                    already_t_min = .false.
                endif
            endif
            count_adaptive_timestep = count_adaptive_timestep + 1
        enddo
        
        rho = rho_new
        if (count_adaptive_timestep == 1) then
            facmax = facmax_init
        endif

        trace_friction_mol_old = trace_friction_mol
        trace_friction_molleads_old = trace_friction_molleads
        trace_current_na_old = trace_current_na
        trace_friction_mol = 0.d0
        trace_friction_molleads = 0.d0
        trace_current_na = 0.d0
        norm = 0.d0
        do itrnz = 0,nnz_elements-1                                                                 ! Now we need to identify the elements of the system density matrix at this timestep, so loop
                                                                                                    ! through all nonzero elements
            indjn = rho_nonzeros(itrnz,0)                                                           ! Find the ADO index of this nonzero element
            itrn = tier_index(indjn)                                                                ! Find the tier of this ADO
            nrow = rho_nonzeros(itrnz,1)                                                            ! Find the row and column of this nonzero ADO element
            ncol = rho_nonzeros(itrnz,2)
            if (itrn == 0) then
                if ((deriv_ham_log(ncol,nrow) .eqv. .true.)) then
                    if (isreal_sparse(itrnz) .eqv. .true.) then
                        trace_friction_mol = trace_friction_mol + deriv_ham_one_x(ncol,nrow)*rho(itrnz)
                    endif
                endif
                if (nrow == ncol) then
                    norm = norm + rho(itrnz)
                endif
            elseif (itrn == 1) then
                jn = un_ind(indjn,0)                                                                ! Find mode of this 1st tier ADO
                leads_n = ksiglm(jn,0)                                                              ! Find lead index of this mode
                sign_n = 1-ksiglm(jn,1)                                                             ! Find \bar{\sigma} (opposite sign) of this mode
                eldash_n = ksiglm(jn,3)
                if (d_ops_log(ncol,nrow,eldash_n,sign_n) .eqv. .true.) then
                    if (isreal_sparse(itrnz) .eqv. .true.) then
                        trace_friction_molleads = trace_friction_molleads + deriv_v_one_x(leads_n,eldash_n)*2.d0*&
                                        dble(d_ops(ncol,nrow,eldash_n,sign_n))*rho(itrnz)
                    else
                        if (degenerate_levels .eqv. .false.) then
                            do itrel = 0,nel-1
                                trace_current_na(leads_n) = trace_current_na(leads_n) + &
                                       ((-1.0)**sign_n)*2.d0*d_ops(ncol,nrow,itrel,sign_n)*rho(itrnz)*el_lead_couplings_x(leads_n,itrel)
                            enddo
                        else
                            trace_current_na(leads_n) = trace_current_na(leads_n) + &
                                    ((-1.0)**sign_n)*2.d0*d_ops(ncol,nrow,eldash_n,sign_n)*rho(itrnz)*el_lead_couplings_x(leads_n,eldash_n)
                        endif
                    endif
                endif
            elseif (itrn > 1) then
                exit
            endif
        enddo

        time_since_last_print = time - time_of_last_print
        if (print_integrand_yn .eqv. .true.) then
            if (time_since_last_print >= printing_timestep) then
                write(300,'(f15.5,a,f15.10)') time," ",dble(trace_friction_mol + trace_friction_molleads)
                time_of_last_print = time
            endif
        endif
        markovian_friction_coefficient_total_old = markovian_friction_coefficient_total
        current_na_old = current_na
        markovian_friction_coefficient_mol = markovian_friction_coefficient_mol - &
                        (dt_old/2.d0)*(trace_friction_mol + trace_friction_mol_old)
        markovian_friction_coefficient_molleads = markovian_friction_coefficient_molleads - &
                        (dt_old/2.d0)*(trace_friction_molleads + trace_friction_molleads_old)
        markovian_friction_coefficient_total = markovian_friction_coefficient_mol + markovian_friction_coefficient_molleads
        current_na = current_na - (dt_old/2.d0)*(trace_current_na + trace_current_na_old)
        rel_diff = abs(markovian_friction_coefficient_total - markovian_friction_coefficient_total_old)&
                        /abs(markovian_friction_coefficient_total)
        rel_diff_current_na = abs(current_na - current_na_old)/abs(current_na)
        itrt = itrt + 1
        dt_tot = dt_tot + dt

    enddo

    dt_av = dt_tot/itrt

    print *,"Number of timesteps: ",itrt
    print *,"Propagation time: ",time
    print *,"Norm at end of propagation: ",norm
    print *,"Friction: ",dble(markovian_friction_coefficient_total)
    print *,"NA Current: ",current_na
    print *,"Conv. of integral: ",rel_diff
    print *,"Conv. of NA current integral: ",rel_diff_current_na
    print *,"Average timestep size of propagation: ",dt_av

    if (print_integrand_yn .eqv. .true.) then
        close(300)
        close(400)
    endif

end subroutine sparse_markovian_friction_propagation

subroutine sparse_markovian_friction_propagation_w_checking(pair_info_row,pair_info_col,&
    pair_values,rho_nonzeros,tier_index,rho_0,tol_friction,deriv_ham_one_x,deriv_ham_log,&
    max_expan_order,dim_rho,nnz_elements,npairs,len_un_ind,nthreads_liouvillian,&
    rk_coeff,rk_coeffhat,facmin,fac,atol_vec,rtol_vec,max_time,dt_min,&
    complex_coefficients,min_time,time_input,dt_tot_input,norm_input,rel_diff_input,&
    itrt_input,dt_input,facmax_input,markovian_friction_coefficient_mol_input,&
    markovian_friction_coefficient_molleads_input,checking_time_interval,&
    first_timecheck,continue_checking,integrand_output,itrt,norm,rel_diff,markovian_friction_coefficient_mol,&
    markovian_friction_coefficient_molleads,time,dt_tot,dt,isreal_sparse,deriv_v_one_x,nel,nleads,&
    nmodes,nmax,un_ind,ksiglm,d_ops_log,d_ops)

    use mkl_spblas

    implicit none

    integer, intent(in) :: dim_rho,max_expan_order,nnz_elements,npairs,nthreads_liouvillian
    integer, intent(in) :: len_un_ind,nel,nleads,nmodes,nmax
    double complex, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: deriv_ham_one_x                         ! Define input variables and arrays    
    double precision, intent(in), dimension(0:nleads-1,0:nel-1) :: deriv_v_one_x
    double precision, intent(in), dimension(0:nnz_elements-1) :: rho_0
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: deriv_ham_log                         ! Define input variables and arrays   
    double precision, intent(in), dimension(0:npairs-1) :: pair_values
    integer, intent(in), dimension(0:npairs-1) :: pair_info_row,pair_info_col
    integer, intent(in), dimension(0:nnz_elements-1,0:2) :: rho_nonzeros
    double complex, intent(in), dimension(0:nnz_elements-1) :: complex_coefficients
    logical, intent(in), dimension(0:nnz_elements-1) :: isreal_sparse
    integer, intent(in), dimension(0:len_un_ind-1) :: tier_index
    integer, intent(in), dimension(0:len_un_ind-1,0:nmax-1) :: un_ind
    integer, intent(in), dimension(0:nmodes-1,0:3) :: ksiglm
    double precision, intent(in), dimension(0:max_expan_order-1) :: rk_coeff,rk_coeffhat
    double precision, intent(in) :: facmin,fac,tol_friction
    double precision, intent(in) :: max_time,dt_min,min_time
    double precision, intent(in), dimension(0:nnz_elements-1) :: atol_vec,rtol_vec
    double precision, intent(in) :: time_input,dt_tot_input,norm_input,rel_diff_input,dt_input
    double complex, intent(in) :: markovian_friction_coefficient_mol_input,markovian_friction_coefficient_molleads_input
    double precision, intent(in) :: facmax_input,checking_time_interval
    logical, intent(in) :: first_timecheck
    integer, intent(in) :: itrt_input
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops_log
    double precision, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops

    double complex, intent(out) :: markovian_friction_coefficient_mol,markovian_friction_coefficient_molleads
    double precision, intent(out) :: norm,rel_diff,time,dt_tot,dt
    integer, intent(out) :: itrt
    logical, intent(out) :: continue_checking
    double precision, intent(out), dimension(0:nnz_elements-1) :: integrand_output

    double precision, dimension(0:nnz_elements-1) :: rho,rho_old,rho_temp,rho_deriv,sc,rho_diff_adaptive
    double precision, dimension(0:nnz_elements-1) :: rho_new,rhohat_new,abs_rho,abs_rho_new,fac_a
    double complex :: trace_friction_mol,trace_friction_molleads,trace_friction_mol_old
    double complex :: trace_friction_molleads_old,markovian_friction_coefficient_total
    double complex :: markovian_friction_coefficient_total_old
    integer :: indjn,itrnz,nrow,ncol,itrl,jn,leads_n,sign_n,eldash_n
    integer :: indnz1,indnz2,itrn,count_adaptive_timestep,ik_coo,ik_csr
    double precision :: err_adaptive_timestep,facmax,dt_new,dnrm2
    double precision :: dt_old,time_of_last_check,time_since_last_check
    logical :: already_t_min

    TYPE(SPARSE_MATRIX_T) :: sparse_handle_coo,sparse_handle_csr
    TYPE(MATRIX_DESCR) :: descra
    descra % TYPE = SPARSE_MATRIX_TYPE_GENERAL

    rho = rho_0                                                                                     ! Create array containing initial values of all nonzero ADO elements. 
                                                                                                    ! rho is a temporary array that will be updated each timestep to contain all nonzero 
                                                                                                    ! ADO elements at that timstep
    trace_friction_mol = 0.d0
    trace_friction_molleads = 0.d0
    do itrnz = 0,nnz_elements-1                                                                     ! Loop through nonzero elements (elements that must be propagated)
        indjn = rho_nonzeros(itrnz,0)                                                               ! For this nonzero element, find the ADO index and row/column value from rho_nonzero
        nrow = rho_nonzeros(itrnz,1)                                                                ! So find the row and column of this ADO element
        ncol = rho_nonzeros(itrnz,2)
        itrn = tier_index(indjn)
        if ((itrn == 0) .and. (deriv_ham_log(ncol,nrow) .eqv. .true.)) then
            trace_friction_mol = trace_friction_mol + deriv_ham_one_x(ncol,nrow)*rho(itrnz)*complex_coefficients(itrnz)
        elseif (itrn == 1) then
            jn = un_ind(indjn,0)                                                                ! Find mode of this 1st tier ADO
            leads_n = ksiglm(jn,0)                                                              ! Find lead index of this mode
            sign_n = 1-ksiglm(jn,1)                                                             ! Find \bar{\sigma} (opposite sign) of this mode
            eldash_n = ksiglm(jn,3)
            if (d_ops_log(ncol,nrow,eldash_n,sign_n) .eqv. .true.) then
                if (isreal_sparse(itrnz) .eqv. .true.) then
                    trace_friction_molleads = trace_friction_molleads + deriv_v_one_x(leads_n,eldash_n)*2.d0*&
                                    dble(d_ops(ncol,nrow,eldash_n,sign_n))*rho(itrnz)
                endif
            endif
    
            exit
        endif
    enddo

    markovian_friction_coefficient_mol = markovian_friction_coefficient_mol_input
    markovian_friction_coefficient_molleads = markovian_friction_coefficient_molleads_input
    markovian_friction_coefficient_total = markovian_friction_coefficient_mol + markovian_friction_coefficient_molleads
    time = time_input
    open(400,file='friction_integrand_mol_new.dat')
    open(500,file='friction_integrand_molleads_new.dat')
    if (first_timecheck .eqv. .true.) then
        write(400,'(f15.5,a,f15.10)') time," ",dble(trace_friction_mol)
        write(500,'(f15.5,a,f15.10)') time," ",dble(trace_friction_molleads)
    endif
    rel_diff = rel_diff_input
    itrt = itrt_input + 1
    facmax = facmax_input
    dt = dt_input
    dt_tot = dt_tot_input
    norm = norm_input
    time_of_last_check = time
    time_since_last_check = 0.d0

    ik_coo = mkl_sparse_d_create_coo(sparse_handle_coo,SPARSE_INDEX_BASE_ZERO,nnz_elements,nnz_elements,&
                                    npairs,pair_info_row,pair_info_col,pair_values)
    ik_csr = mkl_sparse_convert_csr(sparse_handle_coo,SPARSE_OPERATION_NON_TRANSPOSE,sparse_handle_csr)
    if (ik_csr.ne.0) write(*,*) 'Create handle for sparse matrix is wrong', ik_csr
    ik_csr = MKL_SPARSE_OPTIMIZE( sparse_handle_csr )
    if (ik_csr.ne.0) write(*,*) 'OPTIMIZE wrong', ik_csr
    ik_csr = mkl_sparse_copy(sparse_handle_csr,descra,sparse_handle_csr)
    call mkl_set_num_threads(nthreads_liouvillian)
    call omp_set_dynamic(0)
    call mkl_set_dynamic(0)
    already_t_min = .false.
    do while (((rel_diff >= tol_friction) .and. (time <= max_time) .and. &
                (time_since_last_check < checking_time_interval)) .or. (time <= min_time))
        err_adaptive_timestep = 2.d0
        count_adaptive_timestep = 0
        do while (err_adaptive_timestep > 1.d0)
            call dcopy(nnz_elements,rho,1,rho_temp,1)
            call dcopy(nnz_elements,rho,1,rho_new,1)
            call dcopy(nnz_elements,rho,1,rhohat_new,1)
            do itrl = 0,max_expan_order-1                                                               ! Loop through the Taylor series expansion of e^{L*dt}
                rho_deriv = 0.d0
                ik_csr = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,dt,sparse_handle_csr,descra,rho_temp,1.0d0,rho_deriv)
                call daxpy(nnz_elements,rk_coeff(itrl),rho_deriv,1,rho_new,1)
                call daxpy(nnz_elements,rk_coeffhat(itrl),rho_deriv,1,rhohat_new,1)
                call dcopy(nnz_elements,rho_deriv,1,rho_temp,1)
            enddo                  
            call vdAbs(nnz_elements,rho,abs_rho)
            call vdAbs(nnz_elements,rho_new,abs_rho_new)
            call dcopy(nnz_elements,atol_vec,1,sc,1)
            call vdFmax (nnz_elements,abs_rho,abs_rho_new,fac_a)
            call vdmul(nnz_elements,fac_a,rtol_vec,fac_a)
            call daxpy(nnz_elements,1.d0,fac_a,1,sc,1)
            call dcopy(nnz_elements,rho_new,1,rho_diff_adaptive,1)
            call daxpy(nnz_elements,-1.d0,rhohat_new,1,rho_diff_adaptive,1)
            call vddiv(nnz_elements,rho_diff_adaptive,sc,rho_diff_adaptive)
            err_adaptive_timestep = dnrm2(nnz_elements,rho_diff_adaptive,1)/sqrt(dble(nnz_elements))
            dt_new = dt*min(facmax,max(facmin,fac*((1/err_adaptive_timestep)**(1.d0/5.d0))))
            if (err_adaptive_timestep > 1.d0) then
                facmax = 1.d0
                if (dt_new < dt_min) then
                    if (already_t_min .neqv. .true.) then
                        dt_old = dt_min
                        dt = dt_min
                        call dcopy(nnz_elements,rho,1,rho_temp,1)
                        call dcopy(nnz_elements,rho,1,rho_new,1)
                        call dcopy(nnz_elements,rho,1,rhohat_new,1)
                        do itrl = 0,max_expan_order-1                                                               ! Loop through the Taylor series expansion of e^{L*dt}
                            rho_deriv = 0.d0
                            ik_csr = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,dt,sparse_handle_csr,descra,rho_temp,1.0d0,rho_deriv)
                            call daxpy(nnz_elements,rk_coeff(itrl),rho_deriv,1,rho_new,1)
                            call daxpy(nnz_elements,rk_coeffhat(itrl),rho_deriv,1,rhohat_new,1)
                            call dcopy(nnz_elements,rho_deriv,1,rho_temp,1)
                        enddo
                        already_t_min = .true.
                        time = time + dt
                        err_adaptive_timestep = 0.d0
                        count_adaptive_timestep = count_adaptive_timestep + 1
                    else
                        dt_old = dt_min
                        dt = dt_min
                        time = time + dt
                        err_adaptive_timestep = 0.d0
                    endif
                else
                    dt_old = dt
                    dt = dt_new
                    already_t_min = .false.
                endif
            else
                dt_old = dt
                time = time + dt
                if (dt_new < dt_min) then
                    dt = dt_min
                    already_t_min = .true.
                else
                    dt = dt_new
                    already_t_min = .false.
                endif
            endif
            count_adaptive_timestep = count_adaptive_timestep + 1
        enddo
        
        rho = rho_new
        if (count_adaptive_timestep == 1) then
            facmax = facmax_input
        endif

        trace_friction_mol_old = trace_friction_mol
        trace_friction_molleads_old = trace_friction_molleads
        trace_friction_mol = 0.d0
        trace_friction_molleads = 0.d0
        norm = 0.d0
        do itrnz = 0,nnz_elements-1                                                                 ! Now we need to identify the elements of the system density matrix at this timestep, so loop
                                                                                                    ! through all nonzero elements
            indjn = rho_nonzeros(itrnz,0)                                                           ! Find the ADO index of this nonzero element
            itrn = tier_index(indjn)                                                                ! Find the tier of this ADO
            nrow = rho_nonzeros(itrnz,1)                                                            ! Find the row and column of this nonzero ADO element
            ncol = rho_nonzeros(itrnz,2)
            if (itrn == 0) then
                if ((deriv_ham_log(ncol,nrow) .eqv. .true.)) then
                    trace_friction_mol = trace_friction_mol + deriv_ham_one_x(ncol,nrow)*rho(itrnz)*complex_coefficients(itrnz)
                endif
                if (nrow == ncol) then
                    norm = norm + rho(itrnz)
                endif
            elseif (itrn == 1) then
                jn = un_ind(indjn,0)                                                                ! Find mode of this 1st tier ADO
                leads_n = ksiglm(jn,0)                                                              ! Find lead index of this mode
                sign_n = 1-ksiglm(jn,1)                                                             ! Find \bar{\sigma} (opposite sign) of this mode
                eldash_n = ksiglm(jn,3)
                if (d_ops_log(ncol,nrow,eldash_n,sign_n) .eqv. .true.) then
                    if (isreal_sparse(itrnz) .eqv. .true.) then
                        trace_friction_molleads = trace_friction_molleads + deriv_v_one_x(leads_n,eldash_n)*2.d0*&
                                        dble(d_ops(ncol,nrow,eldash_n,sign_n))*rho(itrnz)
                    endif
                endif
            elseif (itrn > 1) then
                exit
            endif
        enddo

        time_since_last_check = time - time_of_last_check
        write(400,'(f15.5,a,f15.10)') time," ",dble(trace_friction_mol)
        write(500,'(f15.5,a,f15.10)') time," ",dble(trace_friction_molleads)
        markovian_friction_coefficient_total_old = markovian_friction_coefficient_total
        markovian_friction_coefficient_mol = markovian_friction_coefficient_mol - &
                        (dt_old/2.d0)*(trace_friction_mol + trace_friction_mol_old)
        markovian_friction_coefficient_molleads = markovian_friction_coefficient_molleads - &
                        (dt_old/2.d0)*(trace_friction_molleads + trace_friction_molleads_old)
        markovian_friction_coefficient_total = markovian_friction_coefficient_mol + markovian_friction_coefficient_molleads
        rel_diff = abs(markovian_friction_coefficient_total - markovian_friction_coefficient_total_old)&
                        /abs(markovian_friction_coefficient_total)
        itrt = itrt + 1
        dt_tot = dt_tot + dt

    enddo

    integrand_output = rho
    continue_checking = .true.
    if ((rel_diff < tol_friction) .or. (time > max_time)) then
        continue_checking = .false.
    endif

    close(400)
    close(500)

end subroutine sparse_markovian_friction_propagation_w_checking

subroutine abs_value(a,n,b)

    implicit none

    integer, intent(in) :: n
    double precision, intent(in), dimension(0:n-1) :: a
    
    double precision, intent(out) :: b

    b = sqrt(sum(a**2))

end subroutine abs_value
