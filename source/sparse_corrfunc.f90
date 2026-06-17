! -----------------------------------------------------------------------------
!
!        SPATIALLY-RESOLVED HEOM CORRELATION FUNCTION EVALUATION
!
! -----------------------------------------------------------------------------
!
! This module provides the Fortran backend for the computation of
! electronic force–force correlation functions within the HEOM framework,
! evaluated on a grid of nuclear (vibrational) coordinates.
!
! It is designed to be wrapped using f2py and called from Python.
!
! -----------------------------------------------------------------------------
!
! OVERVIEW
! -----------------------------------------------------------------------------
!
! The routine computes spatially dependent electronic correlation functions
! entering the construction of the electronic friction tensor.
!
! In particular, it evaluates:
!
!   - Markovian limit correlation contributions
!   - Molecular contribution to the force–force correlation function
!   - Molecule–lead coupling contribution
!
! These quantities are computed as a function of the nuclear coordinate grid
! x_vec, allowing reconstruction of position-dependent friction kernels.
!
! The implementation is based on:
!
!   - Sparse HEOM representation of the Liouvillian
!   - Time propagation of auxiliary density operators (ADOs)
!   - OpenMP-parallel evaluation over nuclear coordinate grid points
!   - MKL-accelerated linear algebra kernels (via linked dependencies)
!
! -----------------------------------------------------------------------------
!
! MAIN SUBROUTINE
! -----------------------------------------------------------------------------
!
! spatially_dependent_corrfunc(...)
!
! This routine evaluates correlation functions for each point on a nuclear
! coordinate grid using a preconstructed sparse HEOM representation.
!
! -----------------------------------------------------------------------------
!
! INPUT DATA STRUCTURE
! -----------------------------------------------------------------------------
!
! HEOM sparse structure:
!
!   pair_info_row / pair_info_col
!       Connectivity structure of sparse HEOM Liouvillian.
!
!   pair_values
!       Coupling coefficients defining transitions between HEOM elements
!       as a function of nuclear coordinate x.
!
!   rho_nonzeros
!       Index mapping of nonzero density matrix / ADO elements.
!
!   tier_index
!   un_ind
!       Hierarchy structure defining ADO organization.
!
! -----------------------------------------------------------------------------
!
! SYSTEM OPERATORS
! -----------------------------------------------------------------------------
!
!   deriv_ham
!       Derivative of system Hamiltonian with respect to nuclear coordinate.
!
!   deriv_v
!       Derivative of molecule–lead coupling with respect to nuclear coordinate.
!
!   d_ops
!   d_ops_log
!       Electronic annihilation/creation operators in system Fock basis.
!
! -----------------------------------------------------------------------------
!
! NUMERICAL PARAMETERS
! -----------------------------------------------------------------------------
!
!   rk_coeff, rk_coeffhat
!       Runge–Kutta / exponential integrator coefficients.
!
!   dt_init, dt_min
!       Initial and minimum timestep for adaptive propagation.
!
!   fac, facmin, facmax_init
!       Adaptive timestep control parameters.
!
!   atol_vec, rtol_vec
!       Absolute and relative tolerances per HEOM element.
!
!   max_expan_order
!       Maximum order of exponential/Liouvillian expansion.
!
! -----------------------------------------------------------------------------
!
! GRID DEPENDENCE
! -----------------------------------------------------------------------------
!
!   len_x_vec
!       Number of nuclear coordinate grid points.
!
! The computation is parallelized over this grid using OpenMP.
!
! -----------------------------------------------------------------------------
!
! OUTPUTS
! -----------------------------------------------------------------------------
!
!   markovian_corrfunc_coefficient_mol_vec
!       Molecular contribution to Markovian correlation function
!       as a function of nuclear coordinate.
!
!   markovian_corrfunc_coefficient_molleads_vec
!       Molecule–lead coupling contribution.
!
!   prop_info_corrfunc
!       Propagation diagnostics (timesteps, convergence, norms, etc.)
!
! -----------------------------------------------------------------------------
!
! OPTIONAL OUTPUT CONTROL
! -----------------------------------------------------------------------------
!
!   print_integrand_yn
!       If enabled, additional intermediate HEOM integrand quantities
!       are computed for diagnostic purposes (single-point evaluation mode
!       in nuclear coordinate space).
!
! -----------------------------------------------------------------------------
!
! IMPLEMENTATION DETAILS
! -----------------------------------------------------------------------------
!
! - Sparse HEOM Liouvillian propagation
! - Coordinate-dependent electronic structure
! - OpenMP parallelization over nuclear grid points
! - Mixed real/complex arithmetic depending on sparsity structure
! - Designed for f2py Python integration
!
! -----------------------------------------------------------------------------

include 'mkl_spblas.f90'

subroutine spatially_dependent_corrfunc(pair_info_row,pair_info_col,pair_values,rho_nonzeros,tier_index,&
            tol_corrfunc,max_expan_order,dim_rho,nnz_elements,npairs,len_un_ind,len_x_vec,&
            deriv_ham,deriv_ham_log,nthreads_x_grid,nthreads_liouvillian,&
            prop_info_corrfunc,rk_coeff,rk_coeffhat,dt_init,facmax_init,facmin,fac,atol_vec,&
            rtol_vec,printing_timestep,max_time,dt_min,complex_coefficients,min_time,&
            isreal_sparse,deriv_v,nel,nleads,nmodes,nmax,un_ind,ksiglm,d_ops_log,d_ops,&
            markovian_corrfunc_coefficient_mol_vec,markovian_corrfunc_coefficient_molleads_vec,&
            print_integrand_yn,ic_corrfunc)

    implicit none
    
    include 'omp_lib.h'                                                                             ! Import the OpenMP (OMP) library for parallel programming
    integer, intent(in) :: dim_rho,max_expan_order,nnz_elements,npairs,len_un_ind
    integer, intent(in) :: len_x_vec,nthreads_x_grid,nthreads_liouvillian,nel,nleads,nmodes,nmax
    double precision, intent(in), dimension(0:npairs-1,0:len_x_vec-1,0:4) :: pair_values
    integer, intent(in), dimension(0:npairs-1) :: pair_info_row,pair_info_col
    integer, intent(in), dimension(0:nnz_elements-1,0:2) :: rho_nonzeros    
    logical, intent(in), dimension(0:nnz_elements-1) :: isreal_sparse
    integer, intent(in), dimension(0:len_un_ind-1) :: tier_index
    integer, intent(in), dimension(0:len_un_ind-1,0:nmax-1) :: un_ind
    integer, intent(in), dimension(0:nmodes-1,0:3) :: ksiglm
    double complex, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:len_x_vec-1) :: deriv_ham
    double precision, intent(in), dimension(0:nleads-1,0:nel-1,0:len_x_vec-1) :: deriv_v
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: deriv_ham_log 
    double complex, intent(in), dimension(0:nnz_elements-1) :: complex_coefficients
    double precision, intent(in), dimension(0:max_expan_order-1) :: rk_coeff,rk_coeffhat
    double precision, intent(in) :: dt_init,facmax_init,facmin,fac,tol_corrfunc
    double precision, intent(in) :: printing_timestep,max_time,dt_min,min_time
    double precision, intent(in), dimension(0:nnz_elements-1) :: atol_vec,rtol_vec
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops_log
    double precision, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops
    logical, intent(in) :: print_integrand_yn
    double precision, intent(in), dimension(0:len_x_vec-1,0:nnz_elements) :: ic_corrfunc

    double precision, intent(out), dimension(0:len_x_vec-1):: markovian_corrfunc_coefficient_mol_vec
    double precision, intent(out), dimension(0:len_x_vec-1):: markovian_corrfunc_coefficient_molleads_vec
    double precision, intent(out), dimension(0:len_x_vec-1,0:4) :: prop_info_corrfunc
    
    double precision, dimension(0:nnz_elements-1) :: ic_corrfunc_one_x
    double precision, dimension(0:npairs-1) :: pair_values_one_x
    double complex, dimension(0:dim_rho-1,0:dim_rho-1) :: deriv_ham_one_x
    double precision, dimension(0:nleads-1,0:nel-1) :: deriv_v_one_x
    integer :: itrx,itrt,itrnz,indjn,nrow,ncol,itrn,itr_ado
    double precision :: rel_diff,time,dt_av,norm
    double complex :: markovian_corrfunc_coefficient_mol,markovian_corrfunc_coefficient_molleads
    
    markovian_corrfunc_coefficient_mol_vec = 0.d0
    markovian_corrfunc_coefficient_molleads_vec = 0.d0
    prop_info_corrfunc = 0.d0
    call omp_set_num_threads(nthreads_x_grid)
    !$OMP PARALLEL SHARED(markovian_corrfunc_coefficient_mol_vec,markovian_corrfunc_coefficient_molleads_vec,prop_info_corrfunc)
    !$OMP DO PRIVATE(itrx,markovian_corrfunc_coefficient_mol,markovian_corrfunc_coefficient_molleads,&
    !$OMP            ic_corrfunc_one_x,itrt,rel_diff,time,norm,dt_av,deriv_ham_one_x,deriv_v_one_x,pair_values_one_x)
    do itrx = 0,(len_x_vec-1)
        print *,itrx,"th value of ",len_x_vec,"x values"
        ic_corrfunc_one_x = ic_corrfunc(itrx,1:)
        pair_values_one_x = pair_values(:,itrx,0)
        deriv_ham_one_x = deriv_ham(:,:,itrx)
        deriv_v_one_x = deriv_v(:,:,itrx)
        call sparse_markovian_corrfunc_propagation(pair_info_row,pair_info_col,&
            pair_values_one_x,rho_nonzeros,tier_index,ic_corrfunc_one_x,tol_corrfunc,&
            markovian_corrfunc_coefficient_mol,markovian_corrfunc_coefficient_molleads,&
            max_expan_order,dim_rho,nnz_elements,npairs,len_un_ind,nthreads_liouvillian,&
            itrt,rel_diff,time,rk_coeff,rk_coeffhat,dt_init,facmax_init,facmin,fac,&
            atol_vec,rtol_vec,dt_av,printing_timestep,max_time,norm,dt_min,&
            complex_coefficients,min_time,isreal_sparse,deriv_v_one_x,nel,nleads,nmodes,nmax,un_ind,&
            ksiglm,d_ops_log,d_ops,print_integrand_yn,deriv_ham_one_x,deriv_ham_log)
        prop_info_corrfunc(itrx,0) = float(itrt)
        prop_info_corrfunc(itrx,1) = time
        prop_info_corrfunc(itrx,2) = rel_diff
        prop_info_corrfunc(itrx,3) = norm
        prop_info_corrfunc(itrx,4) = dt_av
        markovian_corrfunc_coefficient_mol_vec(itrx) = dble(markovian_corrfunc_coefficient_mol)
        markovian_corrfunc_coefficient_molleads_vec(itrx) = dble(markovian_corrfunc_coefficient_molleads)
    enddo
    !$OMP END DO 
    !$OMP END PARALLEL

    ! open(100,file='corrfunc_mol.dat')
    ! open(200,file='corrfunc_molleads.dat')
    ! do itrx = 0,(len_x_vec-1)
    !     write(100,'(f15.5,a,f15.10)') x_vec(itrx)," ",markovian_corrfunc_coefficient_mol_vec(itrx)
    !     write(200,'(f15.5,a,f15.10)') x_vec(itrx)," ",markovian_corrfunc_coefficient_molleads_vec(itrx)
    ! enddo
    ! close(100)
    ! close(200)

end subroutine spatially_dependent_corrfunc

subroutine sparse_markovian_corrfunc_propagation(pair_info_row,pair_info_col,&
        pair_values_one_x,rho_nonzeros,tier_index,ic_corrfunc_one_x,tol_corrfunc,&
        markovian_corrfunc_coefficient_mol,markovian_corrfunc_coefficient_molleads,&
        max_expan_order,dim_rho,nnz_elements,npairs,len_un_ind,nthreads_liouvillian,&
        itrt,rel_diff,time,rk_coeff,rk_coeffhat,dt_init,facmax_init,facmin,fac,&
        atol_vec,rtol_vec,dt_av,printing_timestep,max_time,norm,dt_min,&
        complex_coefficients,min_time,isreal_sparse,deriv_v_one_x,nel,nleads,nmodes,nmax,un_ind,&
        ksiglm,d_ops_log,d_ops,print_integrand_yn,deriv_ham_one_x,deriv_ham_log)

    use mkl_spblas  
                                                
    implicit none

    integer, intent(in) :: dim_rho,max_expan_order,nnz_elements,npairs,nthreads_liouvillian
    integer, intent(in) :: len_un_ind,nel,nleads,nmodes,nmax
    double complex, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: deriv_ham_one_x
    double precision, intent(in), dimension(0:nleads-1,0:nel-1) :: deriv_v_one_x
    double precision, intent(in), dimension(0:nnz_elements-1) :: ic_corrfunc_one_x
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: deriv_ham_log
    double precision, intent(in), dimension(0:npairs-1) :: pair_values_one_x
    integer, intent(in), dimension(0:npairs-1) :: pair_info_row,pair_info_col
    integer, intent(in), dimension(0:nnz_elements-1,0:2) :: rho_nonzeros
    integer, intent(in), dimension(0:len_un_ind-1) :: tier_index
    double precision, intent(in), dimension(0:max_expan_order-1) :: rk_coeff,rk_coeffhat
    double precision, intent(in) :: dt_init,facmax_init,facmin,fac,tol_corrfunc
    double precision, intent(in) :: printing_timestep,max_time,dt_min,min_time
    double precision, intent(in), dimension(0:nnz_elements-1) :: atol_vec,rtol_vec
    double complex, intent(in), dimension(0:nnz_elements-1) :: complex_coefficients
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops_log
    double precision, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops
    logical, intent(in) :: print_integrand_yn
    integer, intent(in), dimension(0:len_un_ind-1,0:nmax-1) :: un_ind
    integer, intent(in), dimension(0:nmodes-1,0:3) :: ksiglm
    logical, intent(in), dimension(0:nnz_elements-1) :: isreal_sparse

    double complex, intent(out) :: markovian_corrfunc_coefficient_mol,markovian_corrfunc_coefficient_molleads
    double precision, intent(out) :: rel_diff,time,dt_av,norm
    integer, intent(out) :: itrt

    double precision, dimension(0:nnz_elements-1) :: rho,rho_old,rho_temp,rho_deriv,sc,rho_diff_adaptive
    double precision, dimension(0:nnz_elements-1) :: rho_new,rhohat_new,abs_rho,abs_rho_new,fac_a
    double complex :: trace_corrfunc_mol,trace_corrfunc_molleads,markovian_corrfunc_coefficient_total
    double complex :: trace_corrfunc_mol_old,trace_corrfunc_molleads_old,markovian_corrfunc_coefficient_mol_old
    double complex :: markovian_corrfunc_coefficient_molleads_old,markovian_corrfunc_coefficient_total_old
    integer :: indjn,itrnz,nrow,ncol,itrl,nnz_elements_sparse,itrn,count_adaptive_timestep
    integer :: jn,leads_n,sign_n,eldash_n
    double precision :: err_adaptive_timestep,facmax,dt,dt_tot,dt_old,ik_coo,ik_csr
    double precision :: time_since_last_print,time_of_last_print,dnrm2,dt_new
    logical :: already_t_min

    TYPE(SPARSE_MATRIX_T) :: sparse_handle_coo,sparse_handle_csr
    TYPE(MATRIX_DESCR) :: descra
    descra % TYPE = SPARSE_MATRIX_TYPE_GENERAL

    rho = ic_corrfunc_one_x
    trace_corrfunc_mol = 0.d0
    trace_corrfunc_molleads = 0.d0
    do itrnz = 0,nnz_elements-1        
        indjn = rho_nonzeros(itrnz,0)  
        nrow = rho_nonzeros(itrnz,1)   
        ncol = rho_nonzeros(itrnz,2)
        itrn = tier_index(indjn)
        if (itrn == 0) then
            if (deriv_ham_log(ncol,nrow) .eqv. .true.) then
                trace_corrfunc_mol = trace_corrfunc_mol + deriv_ham_one_x(ncol,nrow)*rho(itrnz)*complex_coefficients(itrnz)
            endif
        elseif (itrn == 1) then
            jn = un_ind(indjn,0)                  
            leads_n = ksiglm(jn,0)                
            sign_n = 1-ksiglm(jn,1)               
            eldash_n = ksiglm(jn,3)
            trace_corrfunc_molleads = trace_corrfunc_molleads + deriv_v_one_x(leads_n,eldash_n)*2.d0*&     
                                dble(d_ops(ncol,nrow,eldash_n,sign_n))*rho(itrnz)
        elseif (itrn>1)then
            exit
        endif
    enddo


    markovian_corrfunc_coefficient_mol = 0.d0
    markovian_corrfunc_coefficient_molleads = 0.d0
    time = 0.d0
    if (print_integrand_yn .eqv. .true.) then
        open(300,file='corrfunc_integrand_heom.dat')                                                                   ! and information about time the program takes to run
        write(300,'(f15.5,a,f15.10)') time," ",dble(trace_corrfunc_mol + trace_corrfunc_molleads)
    endif
    rel_diff=1.d0
    norm = 0.d0
    itrt = 0
    facmax = facmax_init
    dt = dt_init
    dt_tot = 0.d0
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
    do while (((rel_diff >= tol_corrfunc) .and. (time <= max_time)) .or. (time <= min_time))
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

        trace_corrfunc_mol_old = trace_corrfunc_mol
        trace_corrfunc_molleads_old = trace_corrfunc_molleads
        trace_corrfunc_mol = 0.d0
        trace_corrfunc_molleads = 0.d0
        norm = 0.d0
        do itrnz = 0,nnz_elements-1                                                                 ! Now we need to identify the elements of the system density matrix at this timestep, so loop
                                                                    ! through all nonzero elements
            indjn = rho_nonzeros(itrnz,0)                                                           ! Find the ADO index of this nonzero element
            itrn = tier_index(indjn)                                                                ! Find the tier of this ADO
            nrow = rho_nonzeros(itrnz,1)                                                            ! Find the row and column of this nonzero ADO element
            ncol = rho_nonzeros(itrnz,2)
            if (itrn == 0) then
                if (deriv_ham_log(ncol,nrow) .eqv. .true.) then
                    trace_corrfunc_mol = trace_corrfunc_mol + deriv_ham_one_x(ncol,nrow)*rho(itrnz)*complex_coefficients(itrnz)
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
                        trace_corrfunc_molleads = trace_corrfunc_molleads + deriv_v_one_x(leads_n,eldash_n)*2.d0*&
                                        dble(d_ops(ncol,nrow,eldash_n,sign_n))*rho(itrnz)
                    endif
                endif
            elseif (itrn > 1) then
                exit
            endif
        enddo

        time_since_last_print = time - time_of_last_print
        if (time_since_last_print >= printing_timestep) then
            write(300,'(f15.5,a,f15.10)') time," ",dble(trace_corrfunc_mol + trace_corrfunc_molleads)
            time_of_last_print = time
        endif
        markovian_corrfunc_coefficient_total_old = markovian_corrfunc_coefficient_total
        markovian_corrfunc_coefficient_mol = markovian_corrfunc_coefficient_mol + &
                        (dt_old/2.d0)*(trace_corrfunc_mol + trace_corrfunc_mol_old)
        markovian_corrfunc_coefficient_molleads = markovian_corrfunc_coefficient_molleads + &
                        (dt_old/2.d0)*(trace_corrfunc_molleads + trace_corrfunc_molleads_old)
        markovian_corrfunc_coefficient_total = markovian_corrfunc_coefficient_mol + markovian_corrfunc_coefficient_molleads
        rel_diff = abs(markovian_corrfunc_coefficient_total - markovian_corrfunc_coefficient_total_old)&
                        /abs(markovian_corrfunc_coefficient_total)
        itrt = itrt + 1
        dt_tot = dt_tot + dt

    enddo

    dt_av = dt_tot/itrt

    print *,"Number of timesteps: ",itrt
    print *,"Propagation time: ",time
    print *,"Norm at end of propagation: ",norm
    print *,"Correlation function: ",dble(markovian_corrfunc_coefficient_mol + markovian_corrfunc_coefficient_molleads)
    print *,"Conv. of integral: ",rel_diff
    print *,"Average timestep size of propagation: ",dt_av

    if (print_integrand_yn .eqv. .true.) then
        close(300)
    endif

end subroutine sparse_markovian_corrfunc_propagation

subroutine sparse_markovian_corrfunc_propagation_w_checking(pair_info_row,pair_info_col,&
    pair_values,rho_nonzeros,tier_index,rho_0,tol_corrfunc,dforce,dforce_log,&
    max_expan_order,dim_rho,nnz_elements,npairs,len_un_ind,nthreads_liouvillian,&
    rk_coeff,rk_coeffhat,facmin,fac,atol_vec,rtol_vec,max_time,dt_min,&
    complex_coefficients,min_time,time_input,dt_tot_input,norm_input,rel_diff_input,&
    itrt_input,dt_input,facmax_input,markovian_corrfunc_coefficient_input,checking_time_interval,&
    first_timecheck,continue_checking,integrand_output,itrt,norm,rel_diff,markovian_corrfunc_coefficient,&
    time,dt_tot,dt)

    use mkl_spblas

    implicit none

    integer, intent(in) :: dim_rho,max_expan_order,nnz_elements,npairs,nthreads_liouvillian
    integer, intent(in) :: len_un_ind
    double complex, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: dforce
    double precision, intent(in), dimension(0:nnz_elements-1) :: rho_0
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: dforce_log                         ! Define input variables and arrays   
    double precision, intent(in), dimension(0:npairs-1) :: pair_values
    integer, intent(in), dimension(0:npairs-1) :: pair_info_row,pair_info_col
    integer, intent(in), dimension(0:nnz_elements-1,0:2) :: rho_nonzeros
    double complex, intent(in), dimension(0:nnz_elements-1) :: complex_coefficients
    integer, intent(in), dimension(0:len_un_ind-1) :: tier_index
    double precision, intent(in), dimension(0:max_expan_order-1) :: rk_coeff,rk_coeffhat
    double precision, intent(in) :: facmin,fac,tol_corrfunc
    double precision, intent(in) :: max_time,dt_min,min_time
    double precision, intent(in), dimension(0:nnz_elements-1) :: atol_vec,rtol_vec
    double precision, intent(in) :: time_input,dt_tot_input,norm_input,rel_diff_input,dt_input
    double precision, intent(in) :: markovian_corrfunc_coefficient_input,facmax_input,checking_time_interval
    logical, intent(in) :: first_timecheck
    integer, intent(in) :: itrt_input

    double complex, intent(out) :: markovian_corrfunc_coefficient
    double precision, intent(out) :: norm,rel_diff,time,dt_tot,dt
    integer, intent(out) :: itrt
    logical, intent(out) :: continue_checking
    double precision, intent(out), dimension(0:nnz_elements-1) :: integrand_output

    double precision, dimension(0:nnz_elements-1) :: rho,rho_old,rho_temp,rho_deriv,sc,rho_diff_adaptive
    double precision, dimension(0:nnz_elements-1) :: rho_new,rhohat_new,abs_rho,abs_rho_new,fac_a
    double complex :: trace_corrfunc,markovian_corrfunc_coefficient_old,trace_corrfunc_old
    integer :: indjn,itrnz,nrow,ncol,itrl
    integer :: indnz1,indnz2,itrn,count_adaptive_timestep,ik_coo,ik_csr
    double precision :: err_adaptive_timestep,facmax,dt_new,dnrm2
    double precision :: dt_old,time_of_last_check,time_since_last_check
    logical :: already_t_min

    TYPE(SPARSE_MATRIX_T) :: sparse_handle_coo,sparse_handle_csr
    TYPE(MATRIX_DESCR) :: descra
    descra % TYPE = SPARSE_MATRIX_TYPE_GENERAL

    rho = rho_0
    trace_corrfunc = 0.d0
    do itrnz = 0,nnz_elements-1                                                                     ! Loop through nonzero elements (elements that must be propagated)
        indjn = rho_nonzeros(itrnz,0)                                                               ! For this nonzero element, find the ADO index and row/column value from rho_nonzero
        nrow = rho_nonzeros(itrnz,1)                                                                ! So find the row and column of this ADO element
        ncol = rho_nonzeros(itrnz,2)
        itrn = tier_index(indjn)
        if ((itrn == 0) .and. (dforce_log(ncol,nrow) .eqv. .true.)) then
            trace_corrfunc = trace_corrfunc + dforce(ncol,nrow)*rho(itrnz)*complex_coefficients(itrnz)
        elseif (itrn > 0) then
            exit
        endif
    enddo

    markovian_corrfunc_coefficient = markovian_corrfunc_coefficient_input
    time = time_input
    open(400,file='corrfunc_integrand_new.dat')
    if (first_timecheck .eqv. .true.) then
        write(400,'(f15.5,a,f15.10)') time," ",dble(trace_corrfunc)
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
    do while (((rel_diff >= tol_corrfunc) .and. (time <= max_time) .and. &
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

        trace_corrfunc_old = trace_corrfunc
        trace_corrfunc = 0.d0
        norm = 0.d0
        do itrnz = 0,nnz_elements-1                                                                 ! Now we need to identify the elements of the system density matrix at this timestep, so loop
                                                                                                    ! through all nonzero elements
            indjn = rho_nonzeros(itrnz,0)                                                           ! Find the ADO index of this nonzero element
            itrn = tier_index(indjn)                                                                ! Find the tier of this ADO
            nrow = rho_nonzeros(itrnz,1)                                                            ! Find the row and column of this nonzero ADO element
            ncol = rho_nonzeros(itrnz,2)
            if (itrn == 0) then
                if ((dforce_log(ncol,nrow) .eqv. .true.)) then
                    trace_corrfunc = trace_corrfunc + dforce(ncol,nrow)*rho(itrnz)*complex_coefficients(itrnz)
                endif
                if (nrow == ncol) then
                    norm = norm + rho(itrnz)
                endif
            elseif (itrn > 0) then
                exit
            endif
        enddo

        time_since_last_check = time - time_of_last_check
        write(400,'(f15.5,a,f15.10)') time," ",dble(trace_corrfunc)
        markovian_corrfunc_coefficient_old = markovian_corrfunc_coefficient
        markovian_corrfunc_coefficient = markovian_corrfunc_coefficient + 0.5*(dt_old/2.d0)*&
                                            (trace_corrfunc+trace_corrfunc_old)
        rel_diff = abs(markovian_corrfunc_coefficient - markovian_corrfunc_coefficient_old)/abs(markovian_corrfunc_coefficient)
        itrt = itrt + 1
        dt_tot = dt_tot + dt

    enddo

    integrand_output = rho
    continue_checking = .true.
    if ((rel_diff < tol_corrfunc) .or. (time > max_time)) then
        continue_checking = .false.
    endif
    close(400)

end subroutine sparse_markovian_corrfunc_propagation_w_checking

subroutine abs_value(a,n,b)

    implicit none

    integer, intent(in) :: n
    double precision, intent(in), dimension(0:n-1) :: a
    
    double precision, intent(out) :: b

    b = sqrt(sum(a**2))

end subroutine abs_value
