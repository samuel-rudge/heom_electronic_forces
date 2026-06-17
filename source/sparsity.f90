! ---------------------------------------------------------------------------
! 
!           SERIES OF SUBROUTINES DETAILING SPARSITY OF HEOM
!
! ---------------------------------------------------------------------------
!
! The Fortran subroutines in this program are meant to be converted to a Python wrapper using f2py, for use 
! in the main Python code. The various subroutines detail the sparse nature of the HEOM.
!
! USAGE - RUN FROM TERMINAL (LINUX) TO CREATE PYTHON WRAPPER:
!       Check fortran compilers available in your platform:  f2py -c --help-fcompiler
!       For 'ifort': f2py -c -m sparsity sparsity.f90 --opt='-O3' --fcompiler=intelem --f90flags='-openmp -D__OPENMP' -liomp5
!       For 'gfortran': f2py -c -m sparsity sparsity.f90 --opt='-O3' --fcompiler=gnu95 --f90flags='-fopenmp -D__OPENMP' -lgomp
! ----------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! 
!           CALCULATE THE NUMBER OF NONZERO ELEMENTS IN THE HEOM
!
! ---------------------------------------------------------------------------
!
! This Fortran subroutine calculates the number of coupled (nonzero) ADO elements that need to be propagated
! in time.
!
! USAGE - RUN FROM MAIN PYTHON CODE ONCE WRAPPED:
!       rho_sparsity,nnz_elements,rho_out = sparsity.nnz(ksiglm=KsigLm,tier_index=tier_index,index_minus=Index_Minus,index_plus=Index_Plus,d_ops=d_ops_log,ham=Ham_log,
!                                        rho_0=rho_0_log,max_expan_order=max_expan_order,dim_rho=dim_rho,len_index_plus=len_index_plus,len_un_ind=len_un_ind,
!                                        len_index_minus=len_un_ind,nmax=Nmax,nmodes=Nmodes,nel=Nel) 
!
! INPUTS:
!       ksiglm                                      Array of size [1,nmodes] containing all modes j = {K,sigma,l,m}
!
!       tier_index                                  Array of size [1,len_un_ind] containing tier of each unique ADO
!
!       index_minus                                 Array of size [len_un_ind,nmax,4] allowing fast indexing between tier n and n-1
!
!       index_plus                                  Array of size [len_index_plus,nmodes,3] allowing fast indexing between tier n and n+1
!
!       d_ops                                       Array of size [dim_rho,dim_rho,nel,2] containing annihilation and creation operators for each electronic
!                                                   level in the system, converted to logical representation (filled elements = 1 and unfilled elements = 0)
! 
!       ham                                         Array of size [dim_rho,dim_rho] containing system Hamiltonian, converted to logical representation 
!                                                   (filled elements = 1 and unfilled elements = 0)
!
!       rho_0                                       Array of size [dim_rho,dim_rho] containing initial density matrix of system, converted to logical 
!                                                   representation (filled elements = 1 and unfilled elements = 0)
!
!       max_expan_order                             Scalar that determines the maximum order to go to in the aprpoximate expansion of e^{L*dt}
!      
!       dim_rho                                     Scalar that specifies the number of Fock states defining the system
!
!       len_index_plus                              Scalar giving the size of Index_plus
!
!       len_un_ind                                  Scalar determining the number of unique ADOs required in the HEOM
!
!       len_index_minus                             Scalar giving the size of Index_minus (same as len_un_ind)
!   
!       nmax                                        Scalar giving the maximum tier of the hierarchy
!   
!       nmodes                                      Scalar giving the number of modes, j = {K,sigma,l,m}, included in HEOM
!
!       nel                                         Scalar giving the number of electronic levels in the system
!       
!       wbl_yn                                      Integer with value of 0 or 1 for whether wide-band limit is applied
!
! OUTPUTS:
!       rho_sparsity                                Array of size [dim_rho,dim_rho,len_un_ind] containing all ADOs after one timestep, converted to sparse representation.
!                                                   (unfilled elements = -1 and filled elements = nnz, where nnz is their corresponding nonzero index - counting goes across
!                                                   all columns then rows of each ADO, then across all ADOs.)
!       
!       nnz_elements                                Scalar giving total number of elements in ADOs that need to be propagated (i.e. that are nonzero after 1 timestep)
!
!       rho_out                                     An array of size [dim_rho,dim_rho,len_un_ind] containing all ADOs after one timestep, converted to logical 
!                                                   representation (filled elements = 1 and unfilled elements = 0)
 
subroutine nnz(ksiglm,tier_index,index_minus,index_plus,d_ops,ham,&
                rho_0,rho_sparsity,nnz_elements,rho_out,max_expan_order,dim_rho,&
                len_index_plus,len_un_ind,len_index_minus,nmax,nmodes,nel,wbl_yn,&
                degenerate_levels)

    implicit none                                                                                   ! Prevents Fortran from treating all variables that start with the letters i, j, k, l, m and n
                                                                                                    ! as integers and all other variables as real arguments.

    integer, intent(in) :: len_index_plus,len_index_minus,len_un_ind                              ! Define all input variables with their types as well (integer,logical,etc)
    integer, intent(in) :: nmax,nmodes,nel
    integer, intent(in) :: max_expan_order,dim_rho
    logical, intent(in) :: wbl_yn,degenerate_levels
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops                    ! All array indexing starts at 0 to be consistent with Python
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: ham      
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: rho_0
    integer, intent(in), dimension(0:nmodes-1,0:3) :: ksiglm
    integer, intent(in), dimension(0:len_un_ind-1) :: tier_index
    integer, intent(in), dimension(0:len_index_plus-1,0:nmodes-1,0:2) :: index_plus
    integer, intent(in), dimension(0:len_index_minus-1,0:nmax-1,0:3) :: index_minus

    integer, intent(out), dimension(0:dim_rho-1,0:dim_rho-1,0:len_un_ind-1) :: rho_sparsity         ! Define all output variables with their types as well
    logical, intent(out), dimension(0:dim_rho-1,0:dim_rho-1,0:len_un_ind-1) :: rho_out
    integer, intent(out) :: nnz_elements

    logical, dimension(0:dim_rho-1,0:dim_rho-1,0:len_un_ind-1) :: rho_in,rho_logical                ! Define other necessary variables
    ! logical, dimension(0:dim_rho-1,0:dim_rho-1,0:len_un_ind-1) :: rho_out
    integer :: nrow,ncol,indjn

    rho_in = .false.                                                                                ! rho_in is an array of size [dim_rho,dim_rho,len_un_ind] containing all elements of unique ADOs
    rho_logical = .false.
    rho_in(:,:,0) = rho_0(:,:)
    do while (any(rho_logical .neqv. rho_in))
        rho_logical=rho_in
        call one_step_propagation(ksiglm,tier_index,index_minus,index_plus,&
                                    d_ops,ham,rho_in,rho_out,&        ! Call one_step_propagation function that returns rho_out = e^{L*dt}*rho_in,
                                    max_expan_order,dim_rho,len_index_plus,len_un_ind,&                     ! where rho_out is also a logical array. It tells us which elements of the unique ADOs are 
                                    len_index_minus,nmax,nmodes,nel,wbl_yn,degenerate_levels)
        rho_in = rho_out
    enddo

    nnz_elements = 0                                                                                ! Start number of nonzero elements (nnz) count at 0
    rho_sparsity = -1                                                                               ! Create rho_sparsity and fill all elements with -1
    ! rho_out = rho_in
    do indjn = 0,len_un_ind-1                                                                       ! Loop through all indices of unique ADOs
        do nrow = 0,dim_rho-1                                                                       ! Loop through rows of indjn-th ADO
            do ncol = 0,dim_rho-1                                                                   ! Loop through columns of indjn-th ADO
                if (rho_out(nrow,ncol,indjn)) then                                                  ! Check if this element of rho_out is nonzero (true). If false, do nothing
                    rho_sparsity(nrow,ncol,indjn) = nnz_elements                                    ! If it is true, put the current count of nonzero elements in this place
                    nnz_elements = nnz_elements + 1                                                 ! Increase the count of nonzero elements by 1
                endif
            enddo
        enddo
    enddo
    
end subroutine nnz

! ---------------------------------------------------------------------------
! 
!               PERFORM HEOM PROPAGATION BY ONE TIME STEP
!
! ---------------------------------------------------------------------------
!
! This Fortran subroutine takes an array containing all ADOs, rho_in, and uses the HEOM to calculate 
! the ADOs at the next timestep, rho_out = e^{L*dt}*rho_in. Essentially, it evaluates an approximation
! of e^{L*dt}, the time-propagation operator. The propagation is not performed for the actual numerical values,
! but rather for rho_in and rho_out converted to logical arrays; at this point we only care about whether a 
! particular element needs to be propagated in the time-evolution, not its actual value. Everything is also performed in the Hilbert
! space of the ADOs, although in principle one can easily convert it to Liouville space. The approximation is based
! on the 4th-order Runge-Kutta method, which for a system of coupled 1st-order ODEs translates to expanding 
! e^{L*dt} in a Taylor series and keeping up terms up to n = 4:
! rho_out ~= rho_in + dt*L*rho_in + (dt^2)*(L^2)*rho_in/2 + (dt^3)*(L^3)*rho_in/6 + (dt^4)*(L^4)*rho_in/24
! For an explanation of the algorithm, consult the 
! accompanying text. 
!
! USAGE - RUN IN ABOVE FORTRAN SUBROUTINE nnz(...):
!           call one_step_propagation(ksiglm,tier_index,index_minus,index_plus,d_ops,ham,rho_in,rho_out,&
!                                      max_expan_order,dim_rho,len_index_plus,len_un_ind,&                    
!                                       len_index_minus,nmax,nmodes,nel) 
!
! INPUTS:
!       ksiglm,tier_index,index_minus,              Same input parameters as in nnz(...) above
!       index_plus,d_ops,ham,max_expan_order
!       dim_rho,len_index_plus,len_und_ind,
!       len_index_minus,nmax,nmodes,nel,wbl_yn             
!
!       rho_in                                      rho_in is an array of size [dim_rho,dim_rho,len_un_ind] containing all elements of unique ADOs
!                                                   contained in the hierarchy, in logical representation (0 or false if empty, 
!                                                   1 or true if filled).Array of size [1,len_un_ind] containing tier of each unique ADO
!
! OUTPUTS:
!       rho_out                                     rho_out is an array of size [dim_rho,dim_rho,len_un_ind], defined as rho_out=e^{L*dt}*rho_in

subroutine one_step_propagation(ksiglm,tier_index,index_minus,index_plus,&
                                d_ops,ham,rho_in,rho_out,&                                                          ! Call one_step_propagation function that returns rho_out = e^{L*dt}*rho_in,
                                max_expan_order,dim_rho,len_index_plus,len_un_ind,&                     ! where rho_out is also a logical array. It tells us which elements of the unique ADOs are 
                                len_index_minus,nmax,nmodes,nel,wbl_yn,degenerate_levels)

    implicit none

    integer, intent(in) :: len_index_plus,len_index_minus,len_un_ind
    logical, intent(in) :: wbl_yn,degenerate_levels
    integer, intent(in) :: nmax,nmodes,nel
    integer, intent(in) :: max_expan_order,dim_rho
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops                    ! All array indexing starts at 0 to be consistent with Python
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: ham      
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:len_un_ind-1) :: rho_in
    integer, intent(in), dimension(0:nmodes-1,0:3) :: ksiglm
    integer, intent(in), dimension(0:len_un_ind-1) :: tier_index
    integer, intent(in), dimension(0:len_index_plus-1,0:nmodes-1,0:2) :: index_plus
    integer, intent(in), dimension(0:len_index_minus-1,0:nmax-1,0:3) :: index_minus

    logical, intent(out), dimension(0:dim_rho-1,0:dim_rho-1,0:len_un_ind-1) :: rho_out              ! Define output array
    logical, dimension(0:dim_rho-1,0:dim_rho-1,0:len_un_ind-1) :: rho_deriv,rho_temp                ! Define other necessary variables and arrays. rho_deriv and rho_temp work in tandem
                                                                                                    ! to update rho_out for each iteration of the Taylor series expansion. First, before 
                                                                                                    ! loop itrl starts, rho_temp is set to the updating value of the previous loop.
                                                                                                    ! Then, during the loop, rho_deriv is evaluated as rho_deriv = dt*L*rho_temp/itrl
                                                                                                    ! Then, at the end of the loop, rho_out is evaluated as rho_out = rho_out + rho_deriv.
                                                                                                    ! The final step (which is the same as the first step and thus closes the loop), is to 
                                                                                                    ! set rho_temp = rho_deriv in preparation for loop itrl + 1.
    integer :: itrl,itrn,indjn,itrjnm1,jnm1,indjnm1,sign_nm1,el_nm1
    integer :: jnp1,indjnp1,sign_np1,el_np1,conj_nm1,conj_np1
    integer :: itrel1,itrel2,itrel

    rho_temp = rho_in                                                                               ! rho_temp is an array the same size as rho_in and rho_out that will be updated for
                                                                                                    ! each iteration of the Taylor series expansion to contain 
                                                                                                    ! rho_temp = (dt^{itrl})*(L^{itrl})*rho_in/(itrl!). The 0th tier rho_temp is rho_in
    rho_out = rho_in                                                                                ! Apply 0th-order term of approximation (rho_out = rho_in + dt*L*rho_{n} + ... )
    do itrl = 1,max_expan_order                                                                     ! Loop through all orders of the Taylor series expansion of e^{L*dt}
        rho_deriv = .false.                                                                         ! rho_deriv is initiated as entirely false (unfilled) at the start of the loop
        do indjn = 0,len_un_ind-1
            itrn = tier_index(indjn)                                                             ! Set itrn = tier of the current ADO
            rho_deriv(:,:,indjn) = rho_deriv(:,:,indjn) .or. (matmul(ham,rho_temp(:,:,indjn))&      ! Evaluate coherent part of HEOM with commutator [H,rho_temp]. Evaluation is done with logical
                                                        .or. matmul(rho_temp(:,:,indjn),ham))       ! operations so to add it to rho_deriv we need to take the union of rho_deriv and [H,rho].                   
            if (wbl_yn .eqv. .true.) then 
                if (degenerate_levels .eqv. .true.) then
                    do itrel = 0,nel-1
                        rho_deriv(:,:,indjn) = rho_deriv(:,:,indjn) .or. &
                                                        (matmul(d_ops(:,:,itrel,1),&
                                                        matmul(rho_temp(:,:,indjn),&
                                                        d_ops(:,:,itrel,0))) .or.&
                                                        matmul(d_ops(:,:,itrel,0),&
                                                        matmul(rho_temp(:,:,indjn),&
                                                        d_ops(:,:,itrel,1))))
                    enddo
                else
                    do itrel1 = 0,nel-1
                        do itrel2 = 0,nel-1
                            rho_deriv(:,:,indjn) = rho_deriv(:,:,indjn) .or. &
                                                            (matmul(d_ops(:,:,itrel1,0),&
                                                            matmul(rho_temp(:,:,indjn),&
                                                            d_ops(:,:,itrel2,1))) .or.&
                                                            matmul(d_ops(:,:,itrel1,1),&
                                                            matmul(rho_temp(:,:,indjn),&
                                                            d_ops(:,:,itrel2,0))))
                        enddo
                    enddo
                endif
                rho_deriv(:,:,indjn) = rho_deriv(:,:,indjn) .or. rho_temp(:,:,indjn)                ! Incorporate coupling between nth tier ADO and itself (gamma terms in HEOM)
            endif
            if (itrn > 0) then                                                                      ! Only run this next part for tiers > 0
                rho_deriv(:,:,indjn) = rho_deriv(:,:,indjn) .or. rho_temp(:,:,indjn)                ! Incorporate coupling between nth tier ADO and itself (gamma terms in HEOM)
                do itrjnm1 = 0,itrn-1                                                               ! Loop through mode indices in current ADO of tier itrn
                    jnm1 = index_minus(indjn,itrjnm1,0)                                             ! Calculate mode to be removed corresponding to this mode index
                    sign_nm1 = ksiglm(jnm1,1)                                                       ! Find sigma (sign) of this mode
                    el_nm1 = ksiglm(jnm1,3)                                                         ! Find electronic level of this mode
                    indjnm1 = index_minus(indjn,itrjnm1,1)                                          ! Find index of new ADO created by removing this mode (or index of unique ADO this new
                                                                                                    ! ADO can be expressed as)
                    conj_nm1 = index_minus(indjn,itrjnm1,3)                                         ! Find whether a hermiticity relation is used to connect these two ADOs
                    if (conj_nm1 == 1) then                                                         ! If so, then we need to take the Hermitian conjugate of the new (itrn-1)th tier ADO
                        rho_deriv(:,:,indjn) = rho_deriv(:,:,indjn) .or. &
                                        (matmul(d_ops(:,:,el_nm1,sign_nm1),transpose(rho_temp(:,:,indjnm1))) .or. &
                                        matmul(transpose(rho_temp(:,:,indjnm1)),d_ops(:,:,el_nm1,sign_nm1)))
                                                                                                    ! Apply part of HEOM with C operator, connecting nth and (n-1)th tier ADOs,
                                                                                                    ! performed logically. We have to use the Hermitian conjugate of the connecting ADO
                                                                                                    ! for this case, as a hermiticity relation was used.
                    elseif (conj_nm1 == 0) then
                        rho_deriv(:,:,indjn) = rho_deriv(:,:,indjn) .or. &
                                        (matmul(d_ops(:,:,el_nm1,sign_nm1),rho_temp(:,:,indjnm1)) .or. &
                                        matmul(rho_temp(:,:,indjnm1),d_ops(:,:,el_nm1,sign_nm1)))
                                                                                                    ! Do the same, but without conjugation of the ADO, because there was no hermiticity relation
                                                                                                    ! applied to connect these two ADOs.
                    endif
                enddo
            endif

            if (itrn < nmax) then                                                                ! Connect the current ADOs to ADOs in the tier itrn+1, as long as we are not at nmax already
                do jnp1 = 0,nmodes-1                                                                ! Loop through all possible modes to add to the current ADO
                    indjnp1 = index_plus(indjn,jnp1,0)                                              ! Find index of new (itrn+1)-th tier ADO created by adding mode jnp1 to current ADO
                    if (indjnp1 .ne. -1) then                                                       ! Some ADOs are explicitly included from hierarchy (those with the same Grassmann number/mode
                                                                                                    ! twice.) So if the index returns -1, we can immediately move to the next mode 
                        conj_np1 = index_plus(indjn,jnp1,2)                                         ! Determine whether a hermiticity relation was applied to connect these two ADOs
                        sign_np1 = 1-ksiglm(jnp1,1)                                                 ! Find \bar{sigma} (opposite sign) of new (itrn+1)-th tier ADO
                        if (degenerate_levels .eqv. .true.) then
                            el_np1 = ksiglm(jnp1,3)
                            if (conj_np1 == 1) then
                                rho_deriv(:,:,indjn) = rho_deriv(:,:,indjn) .or. (matmul(d_ops(:,:,el_np1,sign_np1),&
                                                transpose(rho_temp(:,:,indjnp1))).or.&
                                                matmul(transpose(rho_temp(:,:,indjnp1)),&
                                                d_ops(:,:,el_np1,sign_np1)))
                                                                                                        ! Apply part of HEOM with A operator, connecting nth and (n+1)th tier ADOs,                                                                                                    ! Apply part of HEOM with C operator, connecting nth and (n-1)th tier ADOs,
                                                                                                        ! performed logically. We have to use the Hermitian conjugate of the connecting ADO
                                                                                                        ! for this case, as a hermiticity relation was used.
                            elseif (conj_np1 == 0) then
                                rho_deriv(:,:,indjn) = rho_deriv(:,:,indjn) .or. (matmul(d_ops(:,:,el_np1,sign_np1),&
                                                rho_temp(:,:,indjnp1)).or.matmul(rho_temp(:,:,indjnp1),&
                                                d_ops(:,:,el_np1,sign_np1)))
                                                                                                        ! Do the same, but without conjugation of the ADO, because there was no hermiticity relation
                                                                                                        ! applied to connect these two ADOs.
                
                            endif
                        else
                            do el_np1 = 0,nel-1
                                if (conj_np1 == 1) then
                                    rho_deriv(:,:,indjn) = rho_deriv(:,:,indjn) .or. (matmul(d_ops(:,:,el_np1,sign_np1),&
                                                    transpose(rho_temp(:,:,indjnp1))).or.&
                                                    matmul(transpose(rho_temp(:,:,indjnp1)),&
                                                    d_ops(:,:,el_np1,sign_np1)))
                                                                                                            ! Apply part of HEOM with A operator, connecting nth and (n+1)th tier ADOs,                                                                                                    ! Apply part of HEOM with C operator, connecting nth and (n-1)th tier ADOs,
                                                                                                            ! performed logically. We have to use the Hermitian conjugate of the connecting ADO
                                                                                                            ! for this case, as a hermiticity relation was used.
                                elseif (conj_np1 == 0) then
                                    rho_deriv(:,:,indjn) = rho_deriv(:,:,indjn) .or. (matmul(d_ops(:,:,el_np1,sign_np1),&
                                                    rho_temp(:,:,indjnp1)).or.matmul(rho_temp(:,:,indjnp1),&
                                                    d_ops(:,:,el_np1,sign_np1)))
                                                                                                            ! Do the same, but without conjugation of the ADO, because there was no hermiticity relation
                                                                                                            ! applied to connect these two ADOs.
                    
                                endif
                            enddo
                        endif
                    endif
                enddo
            endif
        enddo
                                                                                                    ! rho_deriv now contains all logical connections applied by the HEOM
        rho_out = rho_deriv .or. rho_out                                                            ! Update rho_out with rho_deriv. If an element of rho_out is not filled (i.e. false)
                                                                                                    ! and the corresponding element of rho_deriv is filled (i.e. true), then it updates 
                                                                                                    ! that element of rho_out to be true.
        rho_temp = rho_deriv                                                                        ! Update rho_temp to be equal to rho_deriv for the next loop (next term in the Taylor series
                                                                                                    ! expansion)

    enddo

end subroutine one_step_propagation

! ---------------------------------------------------------------------------
! 
!       CALCULATE NUMBER OF CONNECTIONS BETWEEN ADO ELEMENTS IN HEOM
!
! ---------------------------------------------------------------------------
!
! This Fortran subroutine calculates the number of connections between nonzero ADO elements. If we were in Liouville
! space, this would be all nonzero elements of e^{L*dt}. 
!
! USAGE - RUN IN MAIN PYTHON SCRIPT:
!       rho_nonzeros,npairs,nhampairs = sparsity.sparse_matrix_elements_a(ksiglm=KsigLm,tier_index=tier_index,
!                                un_ind=Un_Ind,index_minus=Index_Minus,index_plus=Index_Plus,d_ops=d_ops,ham=Ham,gamma_vec=gamma_vec,eta_vec=eta_vec,
!                                rho_sparsity=rho_sparsity,nnz_elements=nnz_elements,dim_rho=dim_rho,len_index_plus=len_index_plus,len_un_ind=len_un_ind,
!                                len_index_minus=len_un_ind,nmax=Nmax,nmodes=Nmodes,nel=Nel,nsign=Nsign,max_pairs=max_pairs,
!                                max_ham_pairs=max_ham_pairs,nleads=Nleads,npoles=Npoles)
!
! INPUTS:
!       ksiglm,tier_index,index_minus,              Same input parameters as in nnz(...) above
!       index_plus,ham,max_expan_order
!       dim_rho,len_index_plus,len_un_ind,
!       len_index_minus,nmax,nmodes,nel             
!
!       d_ops                                       Annihilation and creation operators in an array of size [dim_rho,dim_rho,nel,2]. Unlike in nnz, these are the actual operators,
!                                                   not their logical format
!
!       ham                                         System Hamiltonian (array of size [dim_rho,dim_rho]), NOT in logical format
!
!       rho_sparsity                                Array of size [dim_rho,dim_rho,len_un_ind] containing all ADOs after one timestep, converted to sparse representation.
!                                                   (unfilled elements = -1 and filled elements = nnz, where nnz is their corresponding nonzero index - counting goes across
!                                                   all columns then rows of each ADO, then across all ADOs.)
!
!       nnz_elements                                Number of nonzero elements in HEOM
!
! OUTPUTS:
!       rho_nonzeros                                Array of size [nnz_elements,3] containing information about the nonzero elements in the HEOM
!                                                   Each row corresponds to a different nonzero element and the columns contain the index of that ADO
!                                                   to which it belongs, and its row and and column, in that order.
!
!       npairs                                      Number of connected pairs of elements in HEOM. Say we have an ADO of index i. If we consider the element of rho_{i} in the 
!                                                   ath row and bth column, rho_{i}(a,b), then its time-evolution could depend on the cth row and dth column of ADO with index j: rho_{j}(c,d).
!                                                   This connection may come from any part of the dissipative part of the HEOM. The total number of these connected pairs is npairs.
!
!       nhampairs                                   The same, but just for the coherent part of the HEOM: the commutator [H,rho]

subroutine sparse_matrix_elements_a(ksiglm,tier_index,index_minus,index_plus,d_ops_log,ham_log,&
                                    rho_sparsity,rho_nonzeros,npairs,nnz_elements,dim_rho,len_index_plus,&
                                    len_un_ind,len_index_minus,nmax,nmodes,nel,degenerate_levels,&
                                    n_indnz2_this_indnz1_max,n_indnz2_prev_indnz1_vec)

    implicit none

    integer, intent(in) :: len_index_plus,len_index_minus,len_un_ind,nnz_elements,nmax,nmodes,nel,dim_rho
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops_log                    ! Define input variables
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: ham_log
    logical, intent(in) :: degenerate_levels
    integer, intent(in), dimension(0:nmodes-1,0:3) :: ksiglm
    integer, intent(in), dimension(0:len_un_ind-1) :: tier_index
    integer, intent(in), dimension(0:len_index_plus-1,0:nmodes-1,0:2) :: index_plus
    integer, intent(in), dimension(0:len_index_minus-1,0:nmax-1,0:3) :: index_minus
    integer, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:len_un_ind-1) :: rho_sparsity

    integer, intent(out), dimension(0:nnz_elements-1,0:2) :: rho_nonzeros                           ! Define output variables
    integer, intent(out) :: npairs,n_indnz2_this_indnz1_max
    integer, intent(out), dimension(0:nnz_elements):: n_indnz2_prev_indnz1_vec

    integer, dimension(0:nnz_elements-1) :: unique_indnz2
    integer :: ind_nz,nrow,ncol,itrn,indjn,ndash,indnz2,indnz1                      ! Define other necessary variables
    integer :: itrjnm1,jnm1,sign_nm1,el_nm1,indjnm1,conj_nm1,conj_np1
    integer :: jnp1,sign_np1,el_np1,eldash_np1,indjnp1,count_pairs_total
    integer :: n_indnz2_this_indnz1,itr_indnz2,indnz2_compare
    logical :: ham_value,operator_value,already_indnz2

    ind_nz = -1                                                                                     ! Start index of nonzero elements at -1 (we define the first index as 0 to agree with Python)
    do indjn = 0,len_un_ind-1                                                                       ! Loop through all unique ADOs
        do nrow = 0,dim_rho-1                                                                       ! Loop through rows and columns of this ADO
            do ncol = 0,dim_rho-1
                if (rho_sparsity(nrow,ncol,indjn) /= -1) then                                       ! Check if this element is nonzero in the HEOM
                    ind_nz = ind_nz + 1                                                             ! If it is, then add 1 to the nonzero index
                    rho_nonzeros(ind_nz,0) = indjn                                                  ! Fill row ind_nz of rho_nonzeros with the ADO index and row/column values
                    rho_nonzeros(ind_nz,1) = nrow
                    rho_nonzeros(ind_nz,2) = ncol
                endif
            enddo
        enddo
    enddo

    count_pairs_total = -1                                                                          ! Start both pair counts at -1 
    n_indnz2_prev_indnz1_vec = 0

    do indnz1 = 0,nnz_elements-1                                                                    ! Loop through all nonzero elements (i.e. elements of ADOs on LHS of HEOM)
        unique_indnz2 = -2
        n_indnz2_this_indnz1 = -1
        indjn = rho_nonzeros(indnz1,0)                                                              ! For this nonzero element, find the ADO index and row/column value from rho_nonzeros
        nrow = rho_nonzeros(indnz1,1)
        ncol = rho_nonzeros(indnz1,2)
        itrn = tier_index(indjn)                                                                    ! Set itrn = tier of the current ADO
                                                                                                    ! Find the tier of this ADO
        do ndash = 0,dim_rho-1                                                                      ! Loop through all rows/columns of this ADO to see which elements connect (RHS of HEOM)
            indnz2 = rho_sparsity(ndash,ncol,indjn)                                                 ! This part pertains to the H*rho part of the commutator. Essentially evaluating whether 
                                                                                                    ! rho_{indjn}(nrow,ncol) = H(nrow,ndash)*rho_{indjn}(ndash,ncol) connects two nonzero 
                                                                                                    ! ADO elements via a nonzero element of the Hamiltonian
            ham_value = ham_log(nrow,ndash)                                                             ! Find the corresponding element of the Hamiltonian
            if ((indnz2 .ne. -1) .and. (ham_value .eqv. .true.)) then                                     ! If the connecting element is nonzero and the Hamiltonian value is also nonzero,
                already_indnz2 = .false.
                do itr_indnz2 = 0,n_indnz2_this_indnz1
                    indnz2_compare = unique_indnz2(itr_indnz2)
                    if (indnz2 == indnz2_compare) then
                        already_indnz2 = .true.
                        exit
                    endif
                enddo
                if (already_indnz2 .eqv. .false.) then
                    count_pairs_total = count_pairs_total + 1
                    n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                    unique_indnz2(n_indnz2_this_indnz1) = indnz2
                endif
            endif

            indnz2 = rho_sparsity(nrow,ndash,indjn)                                                 ! Do the same as before, but now for the rho*H part of the commutator:
                                                                                                    ! rho_{i}(nrow,ncol) = rho_{i}(nrow,ndash)*H(ndash,ncol)
            ham_value = ham_log(ndash,ncol)
            if ((indnz2 .ne. -1) .and. (ham_value .eqv. .true.)) then
                already_indnz2 = .false.
                do itr_indnz2 = 0,n_indnz2_this_indnz1
                    indnz2_compare = unique_indnz2(itr_indnz2)
                    if (indnz2 == indnz2_compare) then
                        already_indnz2 = .true.
                        exit
                    endif
                enddo
                if (already_indnz2 .eqv. .false.) then
                    count_pairs_total = count_pairs_total + 1
                    n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                    unique_indnz2(n_indnz2_this_indnz1) = indnz2
                endif
            endif
        enddo

        if (itrn > 0) then                                                                          ! Next, we move on to pairs in the other parts of the HEOM, starting with connections between 
                                                                                                    ! tier n and tier n-1. We do not need to do the gamma part as we already know that each element
                                                                                                    ! of nonzero ADOs couples to itself in this part.
            indnz2 = indnz1
            already_indnz2 = .false.
            do itr_indnz2 = 0,n_indnz2_this_indnz1
                indnz2_compare = unique_indnz2(itr_indnz2)
                if (indnz2 == indnz2_compare) then
                    already_indnz2 = .true.
                    exit
                endif
            enddo
            if (already_indnz2 .eqv. .false.) then
                count_pairs_total = count_pairs_total + 1
                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                unique_indnz2(n_indnz2_this_indnz1) = indnz2
            endif
            do itrjnm1 = 0,itrn-1                                                                   ! Same process as in one_step_propagation: loop through indices of modes to remove from 
                                                                                                    ! {j_{n},...j_{1}}
                jnm1 = index_minus(indjn,itrjnm1,0)                                                 ! Find removed mode
                sign_nm1 = ksiglm(jnm1,1)                                                           ! Find sigma associated with removed mode
                el_nm1 = ksiglm(jnm1,3)                                                             ! Find electronic level of removed mode
                indjnm1 = index_minus(indjn,itrjnm1,1)                                              ! Find ADO index of new ADO (after it has been transformed to a unique one, if necessary)
                                                                                                    ! created by removing this mode from the current ADO
                conj_nm1 = index_minus(indjn,itrjnm1,3)                                             ! Determine whether a hermiticity relation was required to transform the new ADO to a unique one
                if (conj_nm1 == 1) then                                                             ! If yes, then we need to work with the conjugate of the new ADO
                    do ndash = 0,dim_rho-1                                                          ! If a hermiticity relationship has been applied, loop through all rows/columns of connecting ADO
                        indnz2 = rho_sparsity(ncol,ndash,indjnm1)                                   ! This part pertains to the d^{\sigma}_{m}*rho part of the n,n-1 part of the HEOM. Essentially, 
                                                                                                    ! it evaluates whether 
                                                                                                    ! rho_{indjn}(nrow,ncol) = d^{sign_nm1}_{el_nm1}(nrow,ndash)*rho_{indjnm1}(ncol,ndash)
                                                                                                    ! connects two nonzero ADO elements via a nonzero element of the ann./cre. operator.
                                                                                                    ! Specifically, this line returns -1 if rho_{indjnm1}(ncol,ndash) is a zero element, and its 
                                                                                                    ! nonzero index if it is nonzero. Note that we have had to transpose rho_{indjnm1} because 
                                                                                                    ! of the hermiticity relation.
                        operator_value = d_ops_log(nrow,ndash,el_nm1,sign_nm1)                      ! Find the corresponding ann./cre. operator value
                        if ((indnz2 .ne. - 1) .and. (operator_value .eqv. .true.)) then             ! If all elements are nonzero,
                            already_indnz2 = .false.
                            do itr_indnz2 = 0,n_indnz2_this_indnz1
                                indnz2_compare = unique_indnz2(itr_indnz2)
                                if (indnz2 == indnz2_compare) then
                                    already_indnz2 = .true.
                                    exit
                                endif
                            enddo
                            if (already_indnz2 .eqv. .false.) then
                                count_pairs_total = count_pairs_total + 1
                                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                unique_indnz2(n_indnz2_this_indnz1) = indnz2
                            endif
                        endif

                        indnz2 = rho_sparsity(ndash,nrow,indjnm1)                                   ! This does the same, except now for 
                                                                                                    ! rho_{indjn}(nrow,ncol) = rho_{indjnm1}(ndash,nrow)*d^{sign_nm1}_{el_nm1}(ndash,ncol)
                        operator_value = d_ops_log(ndash,ncol,el_nm1,sign_nm1)
                        if ((indnz2 .ne. - 1) .and. (operator_value .eqv. .true.)) then
                            already_indnz2 = .false.
                            do itr_indnz2 = 0,n_indnz2_this_indnz1
                                indnz2_compare = unique_indnz2(itr_indnz2)
                                if (indnz2 == indnz2_compare) then
                                    already_indnz2 = .true.
                                    exit
                                endif
                            enddo
                            if (already_indnz2 .eqv. .false.) then
                                count_pairs_total = count_pairs_total + 1
                                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                unique_indnz2(n_indnz2_this_indnz1) = indnz2
                            endif
                        endif
                    enddo
                elseif (conj_nm1 == 0) then                                                         ! If we did not need to apply a hermiticity relation to relate the two ADOs, this section runs
                                                                                                    ! instead. It does the same thing, but without transposing rho_{indjnm1}
                    do ndash = 0,dim_rho-1
                        indnz2 = rho_sparsity(ndash,ncol,indjnm1)
                        operator_value = d_ops_log(nrow,ndash,el_nm1,sign_nm1)
                        if ((indnz2 .ne. - 1) .and. (operator_value .eqv. .true.)) then
                            already_indnz2 = .false.
                            do itr_indnz2 = 0,n_indnz2_this_indnz1
                                indnz2_compare = unique_indnz2(itr_indnz2)
                                if (indnz2 == indnz2_compare) then
                                    already_indnz2 = .true.
                                    exit
                                endif
                            enddo
                            if (already_indnz2 .eqv. .false.) then
                                count_pairs_total = count_pairs_total + 1
                                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                unique_indnz2(n_indnz2_this_indnz1) = indnz2
                            endif
                        endif

                        indnz2 = rho_sparsity(nrow,ndash,indjnm1)
                        operator_value = d_ops_log(ndash,ncol,el_nm1,sign_nm1)
                        if ((indnz2 .ne. - 1) .and. (operator_value .eqv. .true.)) then
                            already_indnz2 = .false.
                            do itr_indnz2 = 0,n_indnz2_this_indnz1
                                indnz2_compare = unique_indnz2(itr_indnz2)
                                if (indnz2 == indnz2_compare) then
                                    already_indnz2 = .true.
                                    exit
                                endif
                            enddo
                            if (already_indnz2 .eqv. .false.) then
                                count_pairs_total = count_pairs_total + 1
                                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                unique_indnz2(n_indnz2_this_indnz1) = indnz2
                            endif
                        endif
                    enddo
                endif
            enddo
        endif

        if (itrn < nmax) then                                                                       ! Now we connect the nth and (n+1)th tier ADOs
            do jnp1 = 0,nmodes-1                                                                    ! Loop through all possible modes to add to the current ADO
                indjnp1 = index_plus(indjn,jnp1,0)                                                  ! Find index of new (itrn+1)-th tier ADO created by adding mode jnp1 to current ADO
                if (indjnp1 .ne. -1) then                                                           ! Some ADOs are explicitly included from hierarchy (those with the same Grassmann number/mode
                                                                                                    ! twice.) So if the index returns -1, we can immediately move to the next mode 
                    conj_np1 = index_plus(indjn,jnp1,2)                                             ! Determine whether a hermiticity relation was applied to connect these two ADOs
                    sign_np1 = 1-ksiglm(jnp1,1)                                                     ! Find \bar{sigma} (opposite sign) of new (itrn+1)-th tier ADO
                    eldash_np1 = ksiglm(jnp1,3)                                                         ! Find electronic level of new (itrn+1)-th tier ADO 
                    if (degenerate_levels .eqv. .true.) then
                        el_np1 = eldash_np1
                        if (conj_np1 == 1) then                                                         ! Run this section if a hermiticity relation was applied to relate the two ADOs
                            do ndash = 0,dim_rho-1                                                      ! Loop through all rows/columns of rho_{indjnp1} that could connect to rho_{indjn}(nrow,ncol)
                                indnz2 = rho_sparsity(ncol,ndash,indjnp1)                               ! This part pertains to the d^{\bar{\sigma}}_{m}*rho part of the n,n+1 part of the HEOM.
                                                                                                        ! It evaluates whether 
                                                                                                        ! rho_{indjn}(nrow,ncol) = d^{1-sign_np1}_{el_np1}(nrow,ndash)*rho_{indjnp1}(ncol,ndash)
                                                                                                        ! connects two nonzero ADO elements via a nonzero element of the ann./cre. operator.
                                                                                                        ! Specifically, this line returns -1 if rho_{indjnp1}(ncol,ndash) is a zero element, and its 
                                                                                                        ! nonzero index if it is nonzero. Note that we have had to transpose rho_{indjnp1} because 
                                                                                                        ! of the hermiticity relation.
                                operator_value = d_ops_log(nrow,ndash,el_np1,sign_np1)
                                if ((indnz2 .ne. -1) .and. (operator_value .eqv. .true.)) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                    endif
                                endif

                                indnz2 = rho_sparsity(ndash,nrow,indjnp1)                               ! Does the same, but for the rho*d^{\bar{\sigma}}_{m} part of the n,n+1 part of the HEOM:
                                                                                                        ! rho_{indjn}(nrow,ncol) = rho_{indjnp1}(ndash,nrow)*d^{1-sign_np1}_{el_np1}(ndash,ncol)
                                operator_value = d_ops_log(ndash,ncol,el_np1,sign_np1)
                                if ((indnz2 .ne. -1) .and. (operator_value .eqv. .true.)) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                    endif
                                endif
                            enddo
                        elseif (conj_np1 == 0) then                                                     ! If we did not need to apply a hermiticity relation to relate the two ADOs, this section runs
                                                                                                        ! instead. It does the same thing, but without transposing rho_{indjnp1}
                            do ndash = 0,dim_rho-1
                                indnz2 = rho_sparsity(ndash,ncol,indjnp1)
                                operator_value = d_ops_log(nrow,ndash,el_np1,sign_np1)
                                if ((indnz2 .ne. - 1) .and. (operator_value .eqv. .true.)) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                    endif
                                endif

                                indnz2 = rho_sparsity(nrow,ndash,indjnp1)
                                operator_value = d_ops_log(ndash,ncol,el_np1,sign_np1)
                                if ((indnz2 .ne. - 1) .and. (operator_value .eqv. .true.)) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                    endif
                                endif
                            enddo
                        endif
                    else
                        do el_np1 = 0,nel-1
                            if (conj_np1 == 1) then                                                         ! Run this section if a hermiticity relation was applied to relate the two ADOs
                                do ndash = 0,dim_rho-1                                                      ! Loop through all rows/columns of rho_{indjnp1} that could connect to rho_{indjn}(nrow,ncol)
                                    indnz2 = rho_sparsity(ncol,ndash,indjnp1)                               ! This part pertains to the d^{\bar{\sigma}}_{m}*rho part of the n,n+1 part of the HEOM.
                                                                                                            ! It evaluates whether 
                                                                                                            ! rho_{indjn}(nrow,ncol) = d^{1-sign_np1}_{el_np1}(nrow,ndash)*rho_{indjnp1}(ncol,ndash)
                                                                                                            ! connects two nonzero ADO elements via a nonzero element of the ann./cre. operator.
                                                                                                            ! Specifically, this line returns -1 if rho_{indjnp1}(ncol,ndash) is a zero element, and its 
                                                                                                            ! nonzero index if it is nonzero. Note that we have had to transpose rho_{indjnp1} because 
                                                                                                            ! of the hermiticity relation.
                                    operator_value = d_ops_log(nrow,ndash,el_np1,sign_np1)
                                    if ((indnz2 .ne. -1) .and. (operator_value .eqv. .true.)) then
                                        already_indnz2 = .false.
                                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                                            indnz2_compare = unique_indnz2(itr_indnz2)
                                            if (indnz2 == indnz2_compare) then
                                                already_indnz2 = .true.
                                                exit
                                            endif
                                        enddo
                                        if (already_indnz2 .eqv. .false.) then
                                            count_pairs_total = count_pairs_total + 1
                                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                            unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                        endif
                                    endif

                                    indnz2 = rho_sparsity(ndash,nrow,indjnp1)                               ! Does the same, but for the rho*d^{\bar{\sigma}}_{m} part of the n,n+1 part of the HEOM:
                                                                                                            ! rho_{indjn}(nrow,ncol) = rho_{indjnp1}(ndash,nrow)*d^{1-sign_np1}_{el_np1}(ndash,ncol)
                                    operator_value = d_ops_log(ndash,ncol,el_np1,sign_np1)
                                    if ((indnz2 .ne. -1) .and. (operator_value .eqv. .true.)) then
                                        already_indnz2 = .false.
                                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                                            indnz2_compare = unique_indnz2(itr_indnz2)
                                            if (indnz2 == indnz2_compare) then
                                                already_indnz2 = .true.
                                                exit
                                            endif
                                        enddo
                                        if (already_indnz2 .eqv. .false.) then
                                            count_pairs_total = count_pairs_total + 1
                                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                            unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                        endif
                                    endif
                                enddo
                            elseif (conj_np1 == 0) then                                                     ! If we did not need to apply a hermiticity relation to relate the two ADOs, this section runs
                                                                                                            ! instead. It does the same thing, but without transposing rho_{indjnp1}
                                do ndash = 0,dim_rho-1
                                    indnz2 = rho_sparsity(ndash,ncol,indjnp1)
                                    operator_value = d_ops_log(nrow,ndash,el_np1,sign_np1)
                                    if ((indnz2 .ne. - 1) .and. (operator_value .eqv. .true.)) then
                                        already_indnz2 = .false.
                                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                                            indnz2_compare = unique_indnz2(itr_indnz2)
                                            if (indnz2 == indnz2_compare) then
                                                already_indnz2 = .true.
                                                exit
                                            endif
                                        enddo
                                        if (already_indnz2 .eqv. .false.) then
                                            count_pairs_total = count_pairs_total + 1
                                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                            unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                        endif
                                    endif

                                    indnz2 = rho_sparsity(nrow,ndash,indjnp1)
                                    operator_value = d_ops_log(ndash,ncol,el_np1,sign_np1)
                                    if ((indnz2 .ne. - 1) .and. (operator_value .eqv. .true.)) then
                                        already_indnz2 = .false.
                                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                                            indnz2_compare = unique_indnz2(itr_indnz2)
                                            if (indnz2 == indnz2_compare) then
                                                already_indnz2 = .true.
                                                exit
                                            endif
                                        enddo
                                        if (already_indnz2 .eqv. .false.) then
                                            count_pairs_total = count_pairs_total + 1
                                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                            unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                        endif
                                    endif
                                enddo
                            endif
                        enddo
                    endif
                endif
            enddo
        endif
        if (n_indnz2_this_indnz1 > n_indnz2_this_indnz1_max) then
            n_indnz2_this_indnz1_max = n_indnz2_this_indnz1
        endif
        n_indnz2_prev_indnz1_vec(indnz1+1) = n_indnz2_prev_indnz1_vec(indnz1) + n_indnz2_this_indnz1 + 1
    enddo

    npairs = count_pairs_total + 1 

end subroutine sparse_matrix_elements_a

! ---------------------------------------------------------------------------
! 
!    GENERATE VALUES AND INFORMATION OF EACH COUPLED PAIR OF ADO ELEMENTS 
!
! ---------------------------------------------------------------------------
!
! This Fortran subroutine is similar to the previous sparse_matrix_elements_a subroutine, which calculated
! the number of coupled pairs of ADO elements in the HEOM time-propagation. This subroutine goes further and 
! records information about each coupled pair: the value connecting them and their position in the hierarchy.
! For example, if we know that rho_{i}(a,b) = e^{L*dt}(a,c)*rho_{j}(c,b), then the ADO elements 
! rho_{i}(a,b) (ath row and bth column of ADO with index i) and rho_{j}(c,b) (cth row and bth column of ADO with
! index j) are connected via the ath row and cth column of e^{L*dt}. We need to know this value and these positions
! in order to propagate the HEOM. Note that it is necessary to run sparse_matrix_elements_a before 
! sparse_matrix_elements_b in order to work out the size of the arrays containing the information. Fortran cannot
! dynamically change the shape of arrays (it kind of can but it is slow) and it is much too slow to allocate an 
! array size larger than necessary and then cut it down at the end.
!
! USAGE - RUN IN ABOVE FORTRAN SUBROUTINE nnz(...):
!               pair_info,pair_values,ham_pair_info,gamma_values = sparsity.sparse_matrix_elements_b(ksiglm=KsigLm,tier_index=tier_index,
!                                un_ind=Un_Ind,index_minus=Index_Minus,index_plus=Index_Plus,d_ops=d_ops,ham=Ham,gamma_vec=gamma_vec,eta_vec=eta_vec,
!                                rho_sparsity=rho_sparsity,nnz_elements=nnz_elements,dim_rho=dim_rho,len_index_plus=len_index_plus,len_un_ind=len_un_ind,
!                                len_index_minus=len_un_ind,nmax=Nmax,nmodes=Nmodes,nel=Nel,nsign=Nsign,nleads=Nleads,npoles=Npoles)
!
! INPUTS:
!       ksiglm,tier_index,index_minus,              Same input parameters as in sparse_matrix_elements_a above
!       index_plus,ham,max_expan_order
!       dim_rho,len_index_plus,len_un_ind,
!       len_index_minus,nmax,nmodes,nel,d_ops
!       ham,rho_nonzeros,nnz_elements,
!       rho_sparsity
!
!       gamma_vec,eta_vec                           Arrays containing the exponents and coefficients of the bath-correlation function expansion
!
!       nsign,nleads,npoles                         nsign = 2 (+=0, -=1), nleads =  number of electronic leads, npoles = number of Pade poles
!
! OUTPUTS:
!       pair_info                                   Array of size [npairs,4] containing information about the pairs of coupled nonzero elements in the HEOM.
!                                                   Each row corresponds to a different coupled pair;
!                                                   Column 1 contains the nonzero index of the LHS ADO element, 
!                                                   Column 2 contains the nonzero index of the RHS ADO element,
!                                                   Column 3 contains 1 if the connection between these two ADOs requires a hermiticity relation, and 0 if not
!                                                   Column 4 contains ??? FINISH
!
!       pair_values                                 Array of size [1,npairs] containing the value of this connection (i.e. the corresponding element of e^{L*dt} in Liouville space)            
!
!       ham_pair_info                               The same as pair_info, but just for the coherent part containing the commutator with the Hamiltonian
!
!       gamma_values                                Array of size [1,nnz_elements] containing the sum over gamma values for each ADO in the HEOM


subroutine sparse_matrix_elements_b(ksiglm,tier_index,un_ind,index_minus,index_plus,d_ops,gamma_vec,eta_vec,&
                                rho_sparsity,rho_nonzeros,pair_info_row,pair_info_col,pair_values,&
                                nnz_elements,dim_rho,len_index_plus,len_un_ind,&
                                len_index_minus,nmax,nmodes,nel,nsign,npairs,nleads,npoles,&
                                ham_log,d_ops_log,degenerate_levels,n_indnz2_this_indnz1_max,ham_x,len_x_vec,&
                                n_indnz2_prev_indnz1_vec,el_lead_couplings)

    implicit none

    integer, intent(in) :: len_index_plus,len_index_minus,nsign,len_un_ind,nnz_elements,npairs      ! Define input variables and arrays
    integer, intent(in) :: nmax,nmodes,nel,dim_rho,npoles,nleads,n_indnz2_this_indnz1_max,len_x_vec
    logical, intent(in) :: degenerate_levels
    double precision, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops
    double precision, intent(in), dimension(0:nleads-1,0:nel-1,0:len_x_vec-1,0:4) :: el_lead_couplings
    complex*16, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:len_x_vec-1,0:4) :: ham_x
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops_log
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: ham_log
    integer, intent(in), dimension(0:nmodes-1,0:3) :: ksiglm
    integer, intent(in), dimension(0:len_un_ind-1) :: tier_index
    integer, intent(in), dimension(0:len_un_ind-1,0:nmax-1) :: un_ind
    integer, intent(in), dimension(0:len_index_plus-1,0:nmodes-1,0:2) :: index_plus
    integer, intent(in), dimension(0:len_index_minus-1,0:nmax-1,0:3) :: index_minus
    complex*16, intent(in),dimension(0:nleads-1,0:nsign-1,0:npoles) :: gamma_vec       ! Define complex array output variables
    complex*16, intent(in),dimension(0:nleads-1,0:nsign-1,0:npoles) :: eta_vec
    integer, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:len_un_ind-1) :: rho_sparsity
    integer, intent(in), dimension(0:nnz_elements-1,0:2) :: rho_nonzeros
    integer, intent(in), dimension(0:nnz_elements) :: n_indnz2_prev_indnz1_vec

    integer, intent(out), dimension(0:4*npairs-1) :: pair_info_row,pair_info_col
    double precision, intent(out), dimension(0:4*npairs-1,0:len_x_vec-1,0:4) :: pair_values
    
    integer :: nrow,ncol,itrn,indjn,ndash,indnz2,indnz1,leads_np1                                   ! Define other necessary variables
    integer :: itrgamma,leads_gamma,sign_gamma,j_gamma,poles_gamma,poles_np1
    integer :: itrjnm1,jnm1,sign_nm1,el_nm1,indjnm1,leads_nm1,poles_nm1,conj_np1
    integer :: jnp1,sign_np1,el_np1,eldash_np1,indjnp1,count_pairs_total,conj_nm1
    integer :: n_indnz2_this_indnz1,itr_indnz2,indnz2_compare,pair_index,itrx
    integer :: si_real,si_imag,itrdx
    integer, dimension(0:n_indnz2_this_indnz1_max,0:1) :: unique_indnz2
    complex*16 :: gamma_sum,ham_value
    double precision :: perm_np1,perm_nm1,conj_button_nm1,conj_button_np1,operator_value
    logical :: ham_value_log,operator_value_log,already_indnz2
    complex*16, parameter :: ci=(0.d0,1.d0)                                                      ! Define imaginary number (=sqrt(-1))

    count_pairs_total = -1                                                                       ! Start the pair count at -1
    pair_info_row = -1                                                                           ! Set the indices columns in pair_info = -1 to start
    pair_info_col = -1                                                                           ! Set the indices columns in pair_info = -1 to start
    pair_values = 0.d0                                                                           ! Set pair and gamma values initially to 0 (double precision)

    do indnz1 = 0,nnz_elements-1                                                                 ! Loop through all nonzero elements (i.e. elements of ADOs on LHS of HEOM)
        unique_indnz2 = -2
        n_indnz2_this_indnz1 = -1
        si_real = 4*n_indnz2_prev_indnz1_vec(indnz1)
        si_imag = si_real + 2*(n_indnz2_prev_indnz1_vec(indnz1+1) - n_indnz2_prev_indnz1_vec(indnz1))
        indjn = rho_nonzeros(indnz1,0)                                                           ! For this nonzero element, find the ADO index and row/column value from rho_nonzeros
        nrow = rho_nonzeros(indnz1,1)
        ncol = rho_nonzeros(indnz1,2)
        itrn = tier_index(indjn)                                                                ! Set itrn = tier of the current ADO
                                                                                                ! Find the tier of this ADO
        do ndash = 0,dim_rho-1                                                                      
            indnz2 = rho_sparsity(ndash,ncol,indjn)
            ham_value_log = ham_log(nrow,ndash)
            if ((indnz2 .ne. -1) .and. (ham_value_log .eqv. .true.)) then                                     ! Testing whether rho_{indjn}(nrow,ncol) = H(nrow,ndash)*rho_{indjn}(ndash,ncol)
                already_indnz2 = .false.
                do itr_indnz2 = 0,n_indnz2_this_indnz1
                    indnz2_compare = unique_indnz2(itr_indnz2,0)
                    if (indnz2 == indnz2_compare) then
                        already_indnz2 = .true.
                        exit
                    endif
                enddo
                if (already_indnz2 .eqv. .false.) then
                    count_pairs_total = count_pairs_total + 1
                    n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                    unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                    unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                    pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                    pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                    pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                    pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                    pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                    pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                    pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                    pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                    do itrx = 0,(len_x_vec-1)
                        do itrdx = 0,4
                            ham_value = ham_x(nrow,ndash,itrx,itrdx)
                            pair_values(si_real+2*n_indnz2_this_indnz1,itrx,itrdx) = dble(-ci*ham_value)
                            pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx,itrdx) = -aimag(-ci*ham_value)
                            pair_values(si_imag+2*n_indnz2_this_indnz1,itrx,itrdx) = aimag(-ci*ham_value)
                            pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx,itrdx) =  dble(-ci*ham_value)
                        enddo
                    enddo
                else
                    pair_index = unique_indnz2(itr_indnz2,1)
                    do itrx = 0,(len_x_vec-1)
                        do itrdx = 0,4
                            ham_value = ham_x(nrow,ndash,itrx,itrdx)
                            pair_values(si_real+2*pair_index,itrx,itrdx) = pair_values(si_real+2*pair_index,itrx,itrdx) + &
                                                                             dble(-ci*ham_value)
                            pair_values(si_real+2*pair_index+1,itrx,itrdx) = pair_values(si_real+2*pair_index+1,itrx,itrdx) -&
                                                                             aimag(-ci*ham_value)
                            pair_values(si_imag+2*pair_index,itrx,itrdx) = pair_values(si_imag+2*pair_index,itrx,itrdx) + &
                                                                             aimag(-ci*ham_value)
                            pair_values(si_imag+2*pair_index+1,itrx,itrdx) = pair_values(si_imag+2*pair_index+1,itrx,itrdx) + &
                                                                             dble(-ci*ham_value)
                        enddo
                    enddo
                endif
            endif
    
            indnz2 = rho_sparsity(nrow,ndash,indjn)
            ham_value_log = ham_log(ndash,ncol)
            if ((indnz2 .ne. -1) .and. (ham_value_log .eqv. .true.)) then
                already_indnz2 = .false.
                do itr_indnz2 = 0,n_indnz2_this_indnz1
                    indnz2_compare = unique_indnz2(itr_indnz2,0)
                    if (indnz2 == indnz2_compare) then
                        already_indnz2 = .true.
                        exit
                    endif
                enddo
                if (already_indnz2 .eqv. .false.) then
                    count_pairs_total = count_pairs_total + 1
                    n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                    unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                    unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                    pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                    pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                    pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                    pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                    pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                    pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                    pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                    pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                    do itrx = 0,(len_x_vec-1)
                        do itrdx = 0,4
                            ham_value = ham_x(ndash,ncol,itrx,itrdx)
                            pair_values(si_real+2*n_indnz2_this_indnz1,itrx,itrdx) = dble(ci*ham_value)
                            pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx,itrdx) = -aimag(ci*ham_value)
                            pair_values(si_imag+2*n_indnz2_this_indnz1,itrx,itrdx) = aimag(ci*ham_value)
                            pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx,itrdx) = dble(ci*ham_value)
                        enddo
                    enddo
                else
                    pair_index = unique_indnz2(itr_indnz2,1)
                    do itrx = 0,(len_x_vec-1)
                        do itrdx = 0,4
                            ham_value = ham_x(ndash,ncol,itrx,itrdx)
                            pair_values(si_real+2*pair_index,itrx,itrdx) = pair_values(si_real+2*pair_index,itrx,itrdx) +&
                                                                             dble(ci*ham_value)
                            pair_values(si_real+2*pair_index+1,itrx,itrdx) = pair_values(si_real+2*pair_index+1,itrx,itrdx) -&
                                                                             aimag(ci*ham_value)
                            pair_values(si_imag+2*pair_index,itrx,itrdx) = pair_values(si_imag+2*pair_index,itrx,itrdx) +&
                                                                             aimag(ci*ham_value)
                            pair_values(si_imag+2*pair_index+1,itrx,itrdx) = pair_values(si_imag+2*pair_index+1,itrx,itrdx) +&
                                                                             dble(ci*ham_value)
                        enddo
                    enddo
                endif
            endif
        enddo
        
        if (itrn > 0) then
                                                                                                    ! Each nonzero ADO element directly couples to itself via the gamm term in the HEOM
            gamma_sum = 0.d0                                                                        ! Set sum of gammes to double precision 0 to start
            do itrgamma = 0,itrn-1                                                                  ! Loop through all mode indices of the ADO of this ADO element
                j_gamma = un_ind(indjn,itrgamma)                                                    ! Find the mode corresponding to this mode index
                leads_gamma = ksiglm(j_gamma,0)                                                     ! Find the lead corresponding to this mode
                sign_gamma = ksiglm(j_gamma,1)                                                      ! Find the sigma (sign) corresponding to this mode
                poles_gamma = ksiglm(j_gamma,2)                                                     ! Find the Pade pole corresponding to this mode
                gamma_sum = gamma_sum + gamma_vec(leads_gamma,sign_gamma,poles_gamma)               ! Update the sum over gammas with the gamma corresponding to this mode
            enddo
            indnz2 = indnz1
            already_indnz2 = .false.
            do itr_indnz2 = 0,n_indnz2_this_indnz1
                indnz2_compare = unique_indnz2(itr_indnz2,0)
                if (indnz2 == indnz2_compare) then
                    already_indnz2 = .true.
                    exit
                endif
            enddo
            if (already_indnz2 .eqv. .false.) then
                count_pairs_total = count_pairs_total + 1
                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = dble(-gamma_sum)
                pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = -aimag(-gamma_sum)
                pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = aimag(-gamma_sum)
                pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = dble(-gamma_sum)
            else
                pair_index = unique_indnz2(itr_indnz2,1)
                pair_values(si_real+2*pair_index,:,:) = pair_values(si_real+2*pair_index,:,:) + dble(-gamma_sum)
                pair_values(si_real+2*pair_index+1,:,:) = pair_values(si_real+2*pair_index+1,:,:) - aimag(-gamma_sum)
                pair_values(si_imag+2*pair_index,:,:) = pair_values(si_imag+2*pair_index,:,:) + aimag(-gamma_sum)
                pair_values(si_imag+2*pair_index+1,:,:) = pair_values(si_imag+2*pair_index+1,:,:) + dble(-gamma_sum)
            endif
            
            do itrjnm1 = 0,itrn-1                                                                   ! This section assesses couplings between ADO elements of tier n and tier n-1. First, loop
                                                                                                    ! through mode indices to remove from ADO of current element
                jnm1 = index_minus(indjn,itrjnm1,0)                                                 ! Find mode being removed in this loop
                leads_nm1 = ksiglm(jnm1,0)                                                          ! Find lead index of mode being removed
                sign_nm1 = ksiglm(jnm1,1)                                                           ! Find sigma index (sign) of mode being removed
                poles_nm1 = ksiglm(jnm1,2)                                                          ! Find Pade pole of mode being removed
                el_nm1 = ksiglm(jnm1,3)                                                             ! Find electronic level index of mode being removed
                indjnm1 = index_minus(indjn,itrjnm1,1)                                              ! Find index of new (n-1)-th tier ADO created after this mode is removed
                perm_nm1 = (-1.d0)**(index_minus(indjn,itrjnm1,2))                                  ! Calculate permutation prefactor required to connect new ADO to current ADO
                conj_nm1 = index_minus(indjn,itrjnm1,3)                                             ! Determine whether a hermiticity relation was applied to connect new ADO to current ADO
                if (conj_nm1 == 1) then                                                             ! If a hermiticity relaiton was applied, we need to use the transpose of the new ADO
                    conj_button_nm1 = (-1.d0)**(floor((dble(itrn)-1.d0)/2.d0))                      ! Calculate hermiticity prefactor
                    do ndash = 0,dim_rho-1                                                          ! Loop through all rows/columns that could possibly connect the two ADOs:
                                                                                                    ! rho_{indjn}(nrow,ncol) = eta_{l}*d^{\sigma}_{m}(nrow,ndash)*rho_{indjnm1}(ncol,ndash)
                        indnz2 = rho_sparsity(ncol,ndash,indjnm1)                                   ! Find sparsity index of RHS ADO element: rho_{indjnm1}(ncol,ndash)
                        operator_value = d_ops(nrow,ndash,el_nm1,sign_nm1)                          ! Find corresponding element of ann./cre. operator
                        operator_value_log = d_ops_log(nrow,ndash,el_nm1,sign_nm1)
                        if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then          ! If this ADO element and ann./cre. element are both nonzero,
                            already_indnz2 = .false.
                            do itr_indnz2 = 0,n_indnz2_this_indnz1
                                indnz2_compare = unique_indnz2(itr_indnz2,0)
                                if (indnz2 == indnz2_compare) then
                                    already_indnz2 = .true.
                                    exit
                                endif
                            enddo
                            if (already_indnz2 .eqv. .false.) then
                                count_pairs_total = count_pairs_total + 1
                                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                do itrx = 0,len_x_vec-1
                                    do itrdx = 0,4
                                        pair_values(si_real+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                                    dble(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                                    conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                                    el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                                    aimag(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                                    conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                                    el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                                    aimag(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                                    conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                                    el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                                    -dble(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                                    conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                                    el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                    enddo
                                enddo
                            else
                                pair_index = unique_indnz2(itr_indnz2,1)
                                do itrx = 0,len_x_vec-1
                                    do itrdx = 0,4
                                        pair_values(si_real+2*pair_index,itrx,itrdx) = &
                                                                    pair_values(si_real+2*pair_index,itrx,itrdx) + &
                                                                    dble(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                                    conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                                    el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_real+2*pair_index+1,itrx,itrdx) = &
                                                                    pair_values(si_real+2*pair_index+1,itrx,itrdx) + &
                                                                    aimag(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                                    conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                                    el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_imag+2*pair_index,itrx,itrdx) = &
                                                                    pair_values(si_imag+2*pair_index,itrx,itrdx) + &
                                                                    aimag(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                                    conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                                    el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_imag+2*pair_index+1,itrx,itrdx) = &
                                                                    pair_values(si_imag+2*pair_index+1,itrx,itrdx) - & 
                                                                    dble(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                                    conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                                    el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                    enddo
                                enddo
                            endif
                        endif

                        indnz2 = rho_sparsity(ndash,nrow,indjnm1)                                   ! Do the same, but for 
                                                                                                    ! rho_{indjn}(nrow,ncol) = eta^{*}_{l}*rho_{indjnm1}(ndash,nrow)*d^{\sigma}_{m}(ndash,ncol)
                        operator_value = d_ops(ndash,ncol,el_nm1,sign_nm1)
                        operator_value_log = d_ops_log(ndash,ncol,el_nm1,sign_nm1)
                        if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                            already_indnz2 = .false.
                            do itr_indnz2 = 0,n_indnz2_this_indnz1
                                indnz2_compare = unique_indnz2(itr_indnz2,0)
                                if (indnz2 == indnz2_compare) then
                                    already_indnz2 = .true.
                                    exit
                                endif
                            enddo
                            if (already_indnz2 .eqv. .false.) then
                                count_pairs_total = count_pairs_total + 1
                                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                do itrx = 0,len_x_vec-1
                                    do itrdx = 0,4
                                        pair_values(si_real+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                            dble(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conj_button_nm1*conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                            aimag(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conj_button_nm1*conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                            aimag(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conj_button_nm1*conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                            -dble(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conj_button_nm1*conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                    enddo
                                enddo
                            else
                                pair_index = unique_indnz2(itr_indnz2,1)
                                do itrx = 0,len_x_vec-1
                                    do itrdx = 0,4
                                        pair_values(si_real+2*pair_index,itrx,itrdx) = &
                                            pair_values(si_real+2*pair_index,itrx,itrdx) + &
                                            dble(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conj_button_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_real+2*pair_index+1,itrx,itrdx) = &
                                            pair_values(si_real+2*pair_index+1,itrx,itrdx) + &
                                            aimag(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conj_button_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_imag+2*pair_index,itrx,itrdx) = &
                                            pair_values(si_imag+2*pair_index,itrx,itrdx) + &
                                            aimag(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conj_button_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_imag+2*pair_index+1,itrx,itrdx) = &
                                            pair_values(si_imag+2*pair_index+1,itrx,itrdx) - & 
                                            dble(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conj_button_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                    enddo
                                enddo               
                            endif
                        endif
                    enddo
                elseif (conj_nm1 == 0) then                                                         ! Do the same process as above, but run only if no hermiticity relation is applied and we do 
                                                                                                    ! not need to work with the conjugate-transposed ADO
                    do ndash = 0,dim_rho-1
                        indnz2 = rho_sparsity(ndash,ncol,indjnm1)
                        operator_value = d_ops(nrow,ndash,el_nm1,sign_nm1)
                        operator_value_log = d_ops_log(nrow,ndash,el_nm1,sign_nm1)
                        if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                            already_indnz2 = .false.
                            do itr_indnz2 = 0,n_indnz2_this_indnz1
                                indnz2_compare = unique_indnz2(itr_indnz2,0)
                                if (indnz2 == indnz2_compare) then
                                    already_indnz2 = .true.
                                    exit
                                endif
                            enddo
                            if (already_indnz2 .eqv. .false.) then
                                count_pairs_total = count_pairs_total + 1
                                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                do itrx = 0,len_x_vec-1
                                    do itrdx = 0,4
                                        pair_values(si_real+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                            dble(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)                                                
                                        pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                            -aimag(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                            aimag(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                            dble(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                    enddo
                                enddo
                            else
                                pair_index = unique_indnz2(itr_indnz2,1)
                                do itrx = 0,len_x_vec-1
                                    do itrdx = 0,4
                                        pair_values(si_real+2*pair_index,itrx,itrdx) = &
                                            pair_values(si_real+2*pair_index,itrx,itrdx) + &
                                            dble(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_real+2*pair_index+1,itrx,itrdx) = &
                                            pair_values(si_real+2*pair_index+1,itrx,itrdx) - &
                                            aimag(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_imag+2*pair_index,itrx,itrdx) = &
                                            pair_values(si_imag+2*pair_index,itrx,itrdx) + &
                                            aimag(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_imag+2*pair_index+1,itrx,itrdx) = &
                                            pair_values(si_imag+2*pair_index+1,itrx,itrdx) + & 
                                            dble(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                    enddo                        
                                enddo
                            endif
                        endif
                        indnz2 = rho_sparsity(nrow,ndash,indjnm1)
                        operator_value = d_ops(ndash,ncol,el_nm1,sign_nm1)
                        operator_value_log = d_ops_log(ndash,ncol,el_nm1,sign_nm1)
                        if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                            already_indnz2 = .false.
                            do itr_indnz2 = 0,n_indnz2_this_indnz1
                                indnz2_compare = unique_indnz2(itr_indnz2,0)
                                if (indnz2 == indnz2_compare) then
                                    already_indnz2 = .true.
                                    exit
                                endif
                            enddo
                            if (already_indnz2 .eqv. .false.) then
                                count_pairs_total = count_pairs_total + 1
                                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                do itrx = 0,len_x_vec-1
                                    do itrdx = 0,4
                                        pair_values(si_real+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                            dble(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                            -aimag(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                            aimag(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                            dble(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                    enddo
                                enddo
                            else
                                pair_index = unique_indnz2(itr_indnz2,1)
                                do itrx = 0,len_x_vec-1
                                    do itrdx = 0,4
                                        pair_values(si_real+2*pair_index,itrx,itrdx) = &
                                            pair_values(si_real+2*pair_index,itrx,itrdx) + &
                                            dble(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_real+2*pair_index+1,itrx,itrdx) = &
                                            pair_values(si_real+2*pair_index+1,itrx,itrdx) - &
                                            aimag(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_imag+2*pair_index,itrx,itrdx) = &
                                            pair_values(si_imag+2*pair_index,itrx,itrdx) + &
                                            aimag(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                        pair_values(si_imag+2*pair_index+1,itrx,itrdx) = &
                                            pair_values(si_imag+2*pair_index+1,itrx,itrdx) + & 
                                            dble(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            el_lead_couplings(leads_nm1,el_nm1,itrx,itrdx)*operator_value)
                                    enddo
                                enddo
                            endif
                        endif
                    enddo
                endif
            enddo
        endif

        if (itrn < nmax) then                                                                       ! Now we connect the nth and (n+1)th tier ADOs
            do jnp1 = 0,nmodes-1                                                                    ! Loop through all possible modes to add to the current ADO
                indjnp1 = index_plus(indjn,jnp1,0)                                                  ! Find index of new (itrn+1)-th tier ADO created by adding mode jnp1 to current ADO
                if (indjnp1 .ne. -1) then
                    perm_np1 = (-1.d0)**(index_plus(indjn,jnp1,1))                                  ! Calculate the permutation prefactor generated when connecting these two ADOs
                    conj_np1 = index_plus(indjn,jnp1,2)                                             ! Determine whether a hermiticity relation was applied to connect these two ADOs
                    leads_np1 = ksiglm(jnp1,0)                                                      ! Find lead index of mode being added
                    sign_np1 = 1-ksiglm(jnp1,1)                                                     ! Find \bar{\sigma} (1-sign) index of mode being added
                    poles_np1 = ksiglm(jnp1,2)                                                      ! Find Pade pole of mode being added
                    eldash_np1 = ksiglm(jnp1,3)                                                         ! Find electronic level of mode being added
                    if (degenerate_levels .eqv. .true.) then
                        el_np1 = eldash_np1
                        if (conj_np1 == 1) then                                                         ! Run this code if a hermiticity relation was applied
                            conj_button_np1 = (-1.d0)**(floor((dble(itrn)+1.d0)/2.d0))                         ! Calculate the hermiticity prefactor generated when connecting these two ADOs
                            do ndash = 0,dim_rho-1                                                      ! Loop through all rows/columns of rho_{indjnp1} that could connect to rho_{indjn}(nrow,ncol):
                                                                                                        ! rho_{indjn}(nrow,ncol) = d^{\bar{\sigma}}_{m}(nrow,ndash)*rho_{indjnp1}(ncol,ndash)
                                indnz2 = rho_sparsity(ncol,ndash,indjnp1)
                                operator_value = d_ops(nrow,ndash,el_np1,sign_np1)
                                operator_value_log = d_ops_log(nrow,ndash,el_np1,sign_np1)
                                if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2,0)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                        unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        do itrx = 0,len_x_vec-1
                                            do itrdx = 0,4
                                                pair_values(si_real+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                    dble(-ci*perm_np1*conj_button_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                    aimag(-ci*perm_np1*conj_button_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_imag+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                    aimag(-ci*perm_np1*conj_button_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                    -dble(-ci*perm_np1*conj_button_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                            enddo
                                        enddo
                                    else
                                        pair_index = unique_indnz2(itr_indnz2,1)
                                        do itrx = 0,len_x_vec-1
                                            do itrdx = 0,4
                                                pair_values(si_real+2*pair_index,itrx,itrdx) = &
                                                    pair_values(si_real+2*pair_index,itrx,itrdx) + &
                                                    dble(-ci*perm_np1*conj_button_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_real+2*pair_index+1,itrx,itrdx) = &
                                                    pair_values(si_real+2*pair_index+1,itrx,itrdx) + &
                                                    aimag(-ci*perm_np1*conj_button_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_imag+2*pair_index,itrx,itrdx) = &
                                                    pair_values(si_imag+2*pair_index,itrx,itrdx) + &
                                                    aimag(-ci*perm_np1*conj_button_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_imag+2*pair_index+1,itrx,itrdx) = &
                                                    pair_values(si_imag+2*pair_index+1,itrx,itrdx) - & 
                                                    dble(-ci*perm_np1*conj_button_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                            enddo
                                        enddo
                                    endif
                                endif
                                indnz2 = rho_sparsity(ndash,nrow,indjnp1)                               ! Do the same, but for 
                                                                                                        ! rho_{indjn}(nrow,ncol) = rho_{indjnp1}(ndash,nrow)*d^{\sigma}_{m}(ndash,ncol)
                                operator_value = d_ops(ndash,ncol,el_np1,sign_np1)
                                operator_value_log = d_ops_log(ndash,ncol,el_np1,sign_np1)
                                if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2,0)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                        unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        do itrx = 0,len_x_vec-1
                                            do itrdx = 0,4
                                                pair_values(si_real+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                    dble(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                    aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_imag+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                    aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                    -dble(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                            enddo
                                        enddo
                                    else
                                        pair_index = unique_indnz2(itr_indnz2,1)
                                        do itrx = 0,len_x_vec-1
                                            do itrdx = 0,4
                                                pair_values(si_real+2*pair_index,itrx,itrdx) = &
                                                    pair_values(si_real+2*pair_index,itrx,itrdx) + &
                                                    dble(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_real+2*pair_index+1,itrx,itrdx) = &
                                                    pair_values(si_real+2*pair_index+1,itrx,itrdx) + &
                                                    aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_imag+2*pair_index,itrx,itrdx) = &
                                                    pair_values(si_imag+2*pair_index,itrx,itrdx) + &
                                                    aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_imag+2*pair_index+1,itrx,itrdx) = &
                                                    pair_values(si_imag+2*pair_index+1,itrx,itrdx) - & 
                                                    dble(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                            enddo
                                        enddo
                                    endif
                                endif
                            enddo
                        elseif (conj_np1 == 0) then                                                     ! Run this part instead if hermiticity relation is not applied to connect ADOs 
                                                                                                        ! and we do not need to take the conjugate transpose 
                            do ndash = 0,dim_rho-1
                                indnz2 = rho_sparsity(ndash,ncol,indjnp1)
                                operator_value = d_ops(nrow,ndash,el_np1,sign_np1)
                                operator_value_log = d_ops_log(nrow,ndash,el_np1,sign_np1)
                                if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2,0)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                        unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        do itrx = 0,len_x_vec-1
                                            do itrdx = 0,4
                                                pair_values(si_real+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                    dble(-ci*perm_np1*el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                    -aimag(-ci*perm_np1*el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_imag+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                    aimag(-ci*perm_np1*el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                    dble(-ci*perm_np1*el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                            enddo
                                        enddo
                                    else
                                        pair_index = unique_indnz2(itr_indnz2,1)
                                        do itrx = 0,len_x_vec-1
                                            do itrdx = 0,4
                                                pair_values(si_real+2*pair_index,itrx,itrdx) = &
                                                        pair_values(si_real+2*pair_index,itrx,itrdx) + &
                                                        dble(-ci*perm_np1*el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_real+2*pair_index+1,itrx,itrdx) = &
                                                        pair_values(si_real+2*pair_index+1,itrx,itrdx) - &
                                                        aimag(-ci*perm_np1*el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_imag+2*pair_index,itrx,itrdx) = &
                                                        pair_values(si_imag+2*pair_index,itrx,itrdx) + &
                                                        aimag(-ci*perm_np1*el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_imag+2*pair_index+1,itrx,itrdx) = &
                                                        pair_values(si_imag+2*pair_index+1,itrx,itrdx) + & 
                                                        dble(-ci*perm_np1*el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                            enddo
                                        enddo
                                    endif
                                endif
                    
                                indnz2 = rho_sparsity(nrow,ndash,indjnp1)
                                operator_value = d_ops(ndash,ncol,el_np1,sign_np1)
                                operator_value_log = d_ops_log(ndash,ncol,el_np1,sign_np1)
                                if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2,0)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                        unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        do itrx = 0,len_x_vec-1
                                            do itrdx = 0,4
                                                pair_values(si_real+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                    dble(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                    -aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_imag+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                    aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                    dble(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                            enddo
                                        enddo
                                    else
                                        pair_index = unique_indnz2(itr_indnz2,1)
                                        do itrx = 0,len_x_vec-1
                                            do itrdx = 0,4
                                                pair_values(si_real+2*pair_index,itrx,itrdx) = &
                                                    pair_values(si_real+2*pair_index,itrx,itrdx) + &
                                                    dble(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_real+2*pair_index+1,itrx,itrdx) = &
                                                    pair_values(si_real+2*pair_index+1,itrx,itrdx) - &
                                                    aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_imag+2*pair_index,itrx,itrdx) = &
                                                    pair_values(si_imag+2*pair_index,itrx,itrdx) + &
                                                    aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                pair_values(si_imag+2*pair_index+1,itrx,itrdx) = &
                                                    pair_values(si_imag+2*pair_index+1,itrx,itrdx) + & 
                                                    dble(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                    el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                            enddo
                                        enddo
                                    endif
                                endif
                            enddo
                        endif
                    else
                        do el_np1 = 0,nel-1
                            if (conj_np1 == 1) then                                                         ! Run this code if a hermiticity relation was applied
                                conj_button_np1 = (-1.d0)**(floor((dble(itrn)+1.d0)/2.d0))                         ! Calculate the hermiticity prefactor generated when connecting these two ADOs
                                do ndash = 0,dim_rho-1                                                      ! Loop through all rows/columns of rho_{indjnp1} that could connect to rho_{indjn}(nrow,ncol):
                                                                                                            ! rho_{indjn}(nrow,ncol) = d^{\bar{\sigma}}_{m}(nrow,ndash)*rho_{indjnp1}(ncol,ndash)
                                    indnz2 = rho_sparsity(ncol,ndash,indjnp1)
                                    operator_value = d_ops(nrow,ndash,el_np1,sign_np1)
                                    operator_value_log = d_ops_log(nrow,ndash,el_np1,sign_np1)
                                    if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                                        already_indnz2 = .false.
                                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                                            indnz2_compare = unique_indnz2(itr_indnz2,0)
                                            if (indnz2 == indnz2_compare) then
                                                already_indnz2 = .true.
                                                exit
                                            endif
                                        enddo
                                        if (already_indnz2 .eqv. .false.) then
                                            count_pairs_total = count_pairs_total + 1
                                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                            unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                            unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                            pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                            pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                            pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                            pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                            pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                            pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                            pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                            pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                            do itrx = 0,len_x_vec-1
                                                do itrdx = 0,4
                                                    pair_values(si_real+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                            dble(-ci*perm_np1*conj_button_np1*&
                                                            el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                            aimag(-ci*perm_np1*conj_button_np1*&
                                                            el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_imag+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                            aimag(-ci*perm_np1*conj_button_np1*&
                                                            el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                            -dble(-ci*perm_np1*conj_button_np1*&
                                                            el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                enddo
                                            enddo
                                        else
                                            pair_index = unique_indnz2(itr_indnz2,1)
                                            do itrx = 0,len_x_vec-1
                                                do itrdx = 0,4
                                                    pair_values(si_real+2*pair_index,itrx,itrdx) = &
                                                            pair_values(si_real+2*pair_index,itrx,itrdx) + &
                                                            dble(-ci*perm_np1*conj_button_np1*&
                                                            el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_real+2*pair_index+1,itrx,itrdx) = &
                                                            pair_values(si_real+2*pair_index+1,itrx,itrdx) + &
                                                            aimag(-ci*perm_np1*conj_button_np1*&
                                                            el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_imag+2*pair_index,itrx,itrdx) = &
                                                            pair_values(si_imag+2*pair_index,itrx,itrdx) + &
                                                            aimag(-ci*perm_np1*conj_button_np1*&
                                                            el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_imag+2*pair_index+1,itrx,itrdx) = &
                                                            pair_values(si_imag+2*pair_index+1,itrx,itrdx) - & 
                                                            dble(-ci*perm_np1*conj_button_np1*&
                                                            el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                enddo
                                            enddo
                                        endif
                                    endif
                                    indnz2 = rho_sparsity(ndash,nrow,indjnp1)                               ! Do the same, but for 
                                                                                                            ! rho_{indjn}(nrow,ncol) = rho_{indjnp1}(ndash,nrow)*d^{\sigma}_{m}(ndash,ncol)
                                    operator_value = d_ops(ndash,ncol,el_np1,sign_np1)
                                    operator_value_log = d_ops_log(ndash,ncol,el_np1,sign_np1)
                                    if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                                        already_indnz2 = .false.
                                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                                            indnz2_compare = unique_indnz2(itr_indnz2,0)
                                            if (indnz2 == indnz2_compare) then
                                                already_indnz2 = .true.
                                                exit
                                            endif
                                        enddo
                                        if (already_indnz2 .eqv. .false.) then
                                            count_pairs_total = count_pairs_total + 1
                                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                            unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                            unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                            pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                            pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                            pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                            pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                            pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                            pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                            pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                            pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                            do itrx = 0,len_x_vec-1
                                                do itrdx = 0,4
                                                    pair_values(si_real+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                            dble(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                            el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                            aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                            el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_imag+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                            aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                            el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                            -dble(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                            el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                enddo
                                            enddo
                                        else
                                            pair_index = unique_indnz2(itr_indnz2,1)
                                            do itrx = 0,len_x_vec-1
                                                do itrdx = 0,4
                                                    pair_values(si_real+2*pair_index,itrx,itrdx) = &
                                                            pair_values(si_real+2*pair_index,itrx,itrdx) + &
                                                            dble(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                            el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_real+2*pair_index+1,itrx,itrdx) = &
                                                            pair_values(si_real+2*pair_index+1,itrx,itrdx) + &
                                                            aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                            el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_imag+2*pair_index,itrx,itrdx) = &
                                                            pair_values(si_imag+2*pair_index,itrx,itrdx) + &
                                                            aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                            el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_imag+2*pair_index+1,itrx,itrdx) = &
                                                            pair_values(si_imag+2*pair_index+1,itrx,itrdx) - & 
                                                            dble(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                            el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                enddo
                                            enddo
                                        endif
                                    endif
                                enddo
                            elseif (conj_np1 == 0) then                                                     ! Run this part instead if hermiticity relation is not applied to connect ADOs 
                                                                                                            ! and we do not need to take the conjugate transpose 
                                do ndash = 0,dim_rho-1
                                    indnz2 = rho_sparsity(ndash,ncol,indjnp1)
                                    operator_value = d_ops(nrow,ndash,el_np1,sign_np1)
                                    operator_value_log = d_ops_log(nrow,ndash,el_np1,sign_np1)
                                    if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                                        already_indnz2 = .false.
                                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                                            indnz2_compare = unique_indnz2(itr_indnz2,0)
                                            if (indnz2 == indnz2_compare) then
                                                already_indnz2 = .true.
                                                exit
                                            endif
                                        enddo
                                        if (already_indnz2 .eqv. .false.) then
                                            count_pairs_total = count_pairs_total + 1
                                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                            unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                            unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                            pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                            pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                            pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                            pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                            pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                            pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                            pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                            pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                            do itrx = 0,len_x_vec-1
                                                do itrdx = 0,4
                                                    pair_values(si_real+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                        dble(-ci*perm_np1*el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                        -aimag(-ci*perm_np1*el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_imag+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                        aimag(-ci*perm_np1*el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                        dble(-ci*perm_np1*el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                enddo
                                            enddo
                                        else
                                            pair_index = unique_indnz2(itr_indnz2,1)
                                            do itrx = 0,len_x_vec-1
                                                do itrdx = 0,4
                                                    pair_values(si_real+2*pair_index,itrx,itrdx) = &
                                                        pair_values(si_real+2*pair_index,itrx,itrdx) + &
                                                        dble(-ci*perm_np1*el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_real+2*pair_index+1,itrx,itrdx) = &
                                                        pair_values(si_real+2*pair_index+1,itrx,itrdx) - &
                                                        aimag(-ci*perm_np1*el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_imag+2*pair_index,itrx,itrdx) = &
                                                        pair_values(si_imag+2*pair_index,itrx,itrdx) + &
                                                        aimag(-ci*perm_np1*el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_imag+2*pair_index+1,itrx,itrdx) = &
                                                        pair_values(si_imag+2*pair_index+1,itrx,itrdx) + & 
                                                        dble(-ci*perm_np1*el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                enddo
                                            enddo
                                        endif
                                    endif
                        
                                    indnz2 = rho_sparsity(nrow,ndash,indjnp1)
                                    operator_value = d_ops(ndash,ncol,el_np1,sign_np1)
                                    operator_value_log = d_ops_log(ndash,ncol,el_np1,sign_np1)
                                    if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                                        already_indnz2 = .false.
                                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                                            indnz2_compare = unique_indnz2(itr_indnz2,0)
                                            if (indnz2 == indnz2_compare) then
                                                already_indnz2 = .true.
                                                exit
                                            endif
                                        enddo
                                        if (already_indnz2 .eqv. .false.) then
                                            count_pairs_total = count_pairs_total + 1
                                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                            unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                            unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                            pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                            pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                            pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                            pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                            pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                            pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                            pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                            pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                            do itrx = 0,len_x_vec-1
                                                do itrdx = 0,4
                                                    pair_values(si_real+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                        dble(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                        el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                        -aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                        el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_imag+2*n_indnz2_this_indnz1,itrx,itrdx) = &
                                                        aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                        el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx,itrdx) = &
                                                        dble(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                        el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                enddo
                                            enddo
                                        else
                                            pair_index = unique_indnz2(itr_indnz2,1)
                                            do itrx = 0,len_x_vec-1
                                                do itrdx = 0,4
                                                    pair_values(si_real+2*pair_index,itrx,itrdx) = &
                                                        pair_values(si_real+2*pair_index,itrx,itrdx) + &
                                                        dble(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                        el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_real+2*pair_index+1,itrx,itrdx) = &
                                                        pair_values(si_real+2*pair_index+1,itrx,itrdx) - &
                                                        aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                        el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_imag+2*pair_index,itrx,itrdx) = &
                                                        pair_values(si_imag+2*pair_index,itrx,itrdx) + &
                                                        aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                        el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                    pair_values(si_imag+2*pair_index+1,itrx,itrdx) = &
                                                        pair_values(si_imag+2*pair_index+1,itrx,itrdx) + & 
                                                        dble(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                        el_lead_couplings(leads_np1,el_np1,itrx,itrdx)*operator_value)
                                                enddo
                                            enddo
                                        endif
                                    endif
                                enddo
                            endif
                        enddo
                    endif
                endif
            enddo
        endif
    enddo

end subroutine sparse_matrix_elements_b
! ---------------------------------------------------------------------------
! 
!       CALCULATE NUMBER OF CONNECTIONS BETWEEN ADO ELEMENTS IN HEOM
!                      UNDER THE WIDE-BAND LIMIT
!
! ---------------------------------------------------------------------------
!
! Same as sparse_matrix_elements_a, except under the WBL


subroutine sparse_matrix_elements_a_wbl(ksiglm,tier_index,index_minus,index_plus,ham_log,d_ops_log,&
                                        rho_sparsity,rho_nonzeros,npairs,nnz_elements,&
                                        dim_rho,len_index_plus,len_un_ind,len_index_minus,nmax,nmodes,nel,&
                                        degenerate_levels,n_indnz2_this_indnz1_max,n_indnz2_prev_indnz1_vec)

    implicit none

    integer, intent(in) :: len_index_plus,len_index_minus,len_un_ind,nnz_elements,nmax,nmodes,nel,dim_rho
    logical, intent(in) :: degenerate_levels
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops_log                    ! Define input variables
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: ham_log
    integer, intent(in), dimension(0:nmodes-1,0:3) :: ksiglm
    integer, intent(in), dimension(0:len_un_ind-1) :: tier_index
    integer, intent(in), dimension(0:len_index_plus-1,0:nmodes-1,0:2) :: index_plus
    integer, intent(in), dimension(0:len_index_minus-1,0:nmax-1,0:3) :: index_minus
    integer, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:len_un_ind-1) :: rho_sparsity

    integer, intent(out), dimension(0:nnz_elements-1,0:2) :: rho_nonzeros                           ! Define output variables
    integer, intent(out) :: npairs,n_indnz2_this_indnz1_max
    integer, intent(out), dimension(0:nnz_elements):: n_indnz2_prev_indnz1_vec

    integer, dimension(0:nnz_elements-1) :: unique_indnz2
    integer :: ind_nz,nrow,ncol,itrn,indjn,ndash,indnz2,indnz1,itrel2                      ! Define other necessary variables
    integer :: itrjnm1,jnm1,sign_nm1,el_nm1,indjnm1,conj_nm1,conj_np1,itrel
    integer :: jnp1,sign_np1,el_np1,eldash_np1,indjnp1,count_pairs_total,itrel1
    integer :: n_indnz2_this_indnz1,nddash,itr_indnz2,indnz2_compare
    logical :: ham_value,operator_value,already_indnz2,wbl_value

    ind_nz = -1                                                                                     ! Start index of nonzero elements at -1 (we define the first index as 0 to agree with Python)
    do indjn = 0,len_un_ind-1                                                                       ! Loop through all unique ADOs
        do nrow = 0,dim_rho-1                                                                       ! Loop through rows and columns of this ADO
            do ncol = 0,dim_rho-1
                if (rho_sparsity(nrow,ncol,indjn) /= -1) then                                       ! Check if this element is nonzero in the HEOM
                    ind_nz = ind_nz + 1                                                             ! If it is, then add 1 to the nonzero index
                    rho_nonzeros(ind_nz,0) = indjn                                                  ! Fill row ind_nz of rho_nonzeros with the ADO index and row/column values
                    rho_nonzeros(ind_nz,1) = nrow
                    rho_nonzeros(ind_nz,2) = ncol
                endif
            enddo
        enddo
    enddo

    count_pairs_total = -1    
    n_indnz2_prev_indnz1_vec = 0
    
    do indnz1 = 0,nnz_elements-1                                                                    ! Loop through all nonzero elements (i.e. elements of ADOs on LHS of HEOM)
        unique_indnz2 = -2
        n_indnz2_this_indnz1 = -1
        indjn = rho_nonzeros(indnz1,0)                                                              ! For this nonzero element, find the ADO index and row/column value from rho_nonzeros
        nrow = rho_nonzeros(indnz1,1)
        ncol = rho_nonzeros(indnz1,2)
        itrn = tier_index(indjn)                                                                    ! Set itrn = tier of the current ADO
                                                                                                    ! Find the tier of this ADO
        do ndash = 0,dim_rho-1                                                                      ! Loop through all rows/columns of this ADO to see which elements connect (RHS of HEOM)
            indnz2 = rho_sparsity(ndash,ncol,indjn)                                                 ! This part pertains to the H*rho part of the commutator. Essentially evaluating whether 
                                                                                                    ! rho_{indjn}(nrow,ncol) = H(nrow,ndash)*rho_{indjn}(ndash,ncol) connects two nonzero 
                                                                                                    ! ADO elements via a nonzero element of the Hamiltonian
            ham_value = ham_log(nrow,ndash)                                                             ! Find the corresponding element of the Hamiltonian
            if ((indnz2 .ne. -1) .and. (ham_value .eqv. .true.)) then                                     ! If the connecting element is nonzero and the Hamiltonian value is also nonzero,
                already_indnz2 = .false.
                do itr_indnz2 = 0,n_indnz2_this_indnz1
                    indnz2_compare = unique_indnz2(itr_indnz2)
                    if (indnz2 == indnz2_compare) then
                        already_indnz2 = .true.
                        exit
                    endif
                enddo
                if (already_indnz2 .eqv. .false.) then
                    count_pairs_total = count_pairs_total + 1
                    n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                    unique_indnz2(n_indnz2_this_indnz1) = indnz2
                endif
            endif

            indnz2 = rho_sparsity(nrow,ndash,indjn)                                                 ! Do the same as before, but now for the rho*H part of the commutator:
                                                                                                    ! rho_{i}(nrow,ncol) = rho_{i}(nrow,ndash)*H(ndash,ncol)
            ham_value = ham_log(ndash,ncol)
            if ((indnz2 .ne. -1) .and. (ham_value .eqv. .true.)) then
                already_indnz2 = .false.
                do itr_indnz2 = 0,n_indnz2_this_indnz1
                    indnz2_compare = unique_indnz2(itr_indnz2)
                    if (indnz2 == indnz2_compare) then
                        already_indnz2 = .true.
                        exit
                    endif
                enddo
                if (already_indnz2 .eqv. .false.) then
                    count_pairs_total = count_pairs_total + 1
                    n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                    unique_indnz2(n_indnz2_this_indnz1) = indnz2
                endif
            endif
        enddo

        if (degenerate_levels .eqv. .true.) then
            do ndash = 0,dim_rho-1
                do nddash = 0,dim_rho-1
                    indnz2 = rho_sparsity(ndash,nddash,indjn)
                    if (indnz2 .ne. -1) then
                        do itrel = 0,nel-1
                            wbl_value = d_ops_log(nrow,ndash,itrel,1) .and. d_ops_log(nddash,ncol,itrel,0)
                            if (wbl_value .eqv. .true.) then
                                already_indnz2 = .false.
                                do itr_indnz2 = 0,n_indnz2_this_indnz1
                                    indnz2_compare = unique_indnz2(itr_indnz2)
                                    if (indnz2 == indnz2_compare) then
                                        already_indnz2 = .true.
                                        exit
                                    endif
                                enddo
                                if (already_indnz2 .eqv. .false.) then
                                    count_pairs_total = count_pairs_total + 1
                                    n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                    unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                endif
                            endif
                            wbl_value = d_ops_log(nrow,ndash,itrel,0) .and. d_ops_log(nddash,ncol,itrel,1)
                            if (wbl_value .eqv. .true.) then
                                already_indnz2 = .false.
                                do itr_indnz2 = 0,n_indnz2_this_indnz1
                                    indnz2_compare = unique_indnz2(itr_indnz2)
                                    if (indnz2 == indnz2_compare) then
                                        already_indnz2 = .true.
                                        exit
                                    endif
                                enddo
                                if (already_indnz2 .eqv. .false.) then
                                    count_pairs_total = count_pairs_total + 1
                                    n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                    unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                endif
                            endif
                        enddo
                    endif
                enddo
            enddo
        else
            do ndash = 0,dim_rho-1
                do nddash = 0,dim_rho-1
                    indnz2 = rho_sparsity(ndash,nddash,indjn)
                    if (indnz2 .ne. -1) then
                        do itrel1 = 0,nel-1
                            do itrel2 = 0,nel-1
                                wbl_value = d_ops_log(nrow,ndash,itrel1,0) .and. d_ops_log(nddash,ncol,itrel2,1)
                                if (wbl_value .eqv. .true.) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                    endif
                                endif
                                wbl_value = d_ops_log(nrow,ndash,itrel1,1) .and. d_ops_log(nddash,ncol,itrel2,0)
                                if (wbl_value .eqv. .true.) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                    endif
                                endif
                            enddo
                        enddo
                    endif
                enddo
            enddo
        endif

        indnz2 = indnz1
        already_indnz2 = .false.
        do itr_indnz2 = 0,n_indnz2_this_indnz1
            indnz2_compare = unique_indnz2(itr_indnz2)
            if (indnz2 == indnz2_compare) then
                already_indnz2 = .true.
                exit
            endif
        enddo
        if (already_indnz2 .eqv. .false.) then
            count_pairs_total = count_pairs_total + 1
            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
            unique_indnz2(n_indnz2_this_indnz1) = indnz2
        endif

        if (itrn > 0) then                                                                          ! Next, we move on to pairs in the other parts of the HEOM, starting with connections between 
                                                                                                    ! tier n and tier n-1. We do not need to do the gamma part as we already know that each element
                                                                                                    ! of nonzero ADOs couples to itself in this part.
            do itrjnm1 = 0,itrn-1                                                                   ! Same process as in one_step_propagation: loop through indices of modes to remove from 
                                                                                                    ! {j_{n},...j_{1}}
                jnm1 = index_minus(indjn,itrjnm1,0)                                                 ! Find removed mode
                sign_nm1 = ksiglm(jnm1,1)                                                           ! Find sigma associated with removed mode
                el_nm1 = ksiglm(jnm1,3)                                                             ! Find electronic level of removed mode
                indjnm1 = index_minus(indjn,itrjnm1,1)                                              ! Find ADO index of new ADO (after it has been transformed to a unique one, if necessary)
                                                                                                    ! created by removing this mode from the current ADO
                conj_nm1 = index_minus(indjn,itrjnm1,3)                                             ! Determine whether a hermiticity relation was required to transform the new ADO to a unique one
                if (conj_nm1 == 1) then                                                             ! If yes, then we need to work with the conjugate of the new ADO
                    do ndash = 0,dim_rho-1                                                          ! If a hermiticity relationship has been applied, loop through all rows/columns of connecting ADO
                        indnz2 = rho_sparsity(ncol,ndash,indjnm1)                                   ! This part pertains to the d^{\sigma}_{m}*rho part of the n,n-1 part of the HEOM. Essentially, 
                                                                                                    ! it evaluates whether 
                                                                                                    ! rho_{indjn}(nrow,ncol) = d^{sign_nm1}_{el_nm1}(nrow,ndash)*rho_{indjnm1}(ncol,ndash)
                                                                                                    ! connects two nonzero ADO elements via a nonzero element of the ann./cre. operator.
                                                                                                    ! Specifically, this line returns -1 if rho_{indjnm1}(ncol,ndash) is a zero element, and its 
                                                                                                    ! nonzero index if it is nonzero. Note that we have had to transpose rho_{indjnm1} because 
                                                                                                    ! of the hermiticity relation.
                        operator_value = d_ops_log(nrow,ndash,el_nm1,sign_nm1)                      ! Find the corresponding ann./cre. operator value
                        if ((indnz2 .ne. - 1) .and. (operator_value .eqv. .true.)) then             ! If all elements are nonzero,
                            already_indnz2 = .false.
                            do itr_indnz2 = 0,n_indnz2_this_indnz1
                                indnz2_compare = unique_indnz2(itr_indnz2)
                                if (indnz2 == indnz2_compare) then
                                    already_indnz2 = .true.
                                    exit
                                endif
                            enddo
                            if (already_indnz2 .eqv. .false.) then
                                count_pairs_total = count_pairs_total + 1
                                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                unique_indnz2(n_indnz2_this_indnz1) = indnz2
                            endif
                        endif

                        indnz2 = rho_sparsity(ndash,nrow,indjnm1)                                   ! This does the same, except now for 
                                                                                                    ! rho_{indjn}(nrow,ncol) = rho_{indjnm1}(ndash,nrow)*d^{sign_nm1}_{el_nm1}(ndash,ncol)
                        operator_value = d_ops_log(ndash,ncol,el_nm1,sign_nm1)
                        if ((indnz2 .ne. - 1) .and. (operator_value .eqv. .true.)) then
                            already_indnz2 = .false.
                            do itr_indnz2 = 0,n_indnz2_this_indnz1
                                indnz2_compare = unique_indnz2(itr_indnz2)
                                if (indnz2 == indnz2_compare) then
                                    already_indnz2 = .true.
                                    exit
                                endif
                            enddo
                            if (already_indnz2 .eqv. .false.) then
                                count_pairs_total = count_pairs_total + 1
                                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                unique_indnz2(n_indnz2_this_indnz1) = indnz2
                            endif
                        endif
                    enddo
                elseif (conj_nm1 == 0) then                                                         ! If we did not need to apply a hermiticity relation to relate the two ADOs, this section runs
                                                                                                    ! instead. It does the same thing, but without transposing rho_{indjnm1}
                    do ndash = 0,dim_rho-1
                        indnz2 = rho_sparsity(ndash,ncol,indjnm1)
                        operator_value = d_ops_log(nrow,ndash,el_nm1,sign_nm1)
                        if ((indnz2 .ne. - 1) .and. (operator_value .eqv. .true.)) then
                            already_indnz2 = .false.
                            do itr_indnz2 = 0,n_indnz2_this_indnz1
                                indnz2_compare = unique_indnz2(itr_indnz2)
                                if (indnz2 == indnz2_compare) then
                                    already_indnz2 = .true.
                                    exit
                                endif
                            enddo
                            if (already_indnz2 .eqv. .false.) then
                                count_pairs_total = count_pairs_total + 1
                                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                unique_indnz2(n_indnz2_this_indnz1) = indnz2
                            endif
                        endif

                        indnz2 = rho_sparsity(nrow,ndash,indjnm1)
                        operator_value = d_ops_log(ndash,ncol,el_nm1,sign_nm1)
                        if ((indnz2 .ne. - 1) .and. (operator_value .eqv. .true.)) then
                            already_indnz2 = .false.
                            do itr_indnz2 = 0,n_indnz2_this_indnz1
                                indnz2_compare = unique_indnz2(itr_indnz2)
                                if (indnz2 == indnz2_compare) then
                                    already_indnz2 = .true.
                                    exit
                                endif
                            enddo
                            if (already_indnz2 .eqv. .false.) then
                                count_pairs_total = count_pairs_total + 1
                                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                unique_indnz2(n_indnz2_this_indnz1) = indnz2
                            endif
                        endif
                    enddo
                endif
            enddo
        endif

        if (itrn < nmax) then                                                                       ! Now we connect the nth and (n+1)th tier ADOs
            do jnp1 = 0,nmodes-1                                                                    ! Loop through all possible modes to add to the current ADO
                indjnp1 = index_plus(indjn,jnp1,0)                                                  ! Find index of new (itrn+1)-th tier ADO created by adding mode jnp1 to current ADO
                if (indjnp1 .ne. -1) then                                                           ! Some ADOs are explicitly included from hierarchy (those with the same Grassmann number/mode
                                                                                                    ! twice.) So if the index returns -1, we can immediately move to the next mode 
                    conj_np1 = index_plus(indjn,jnp1,2)                                             ! Determine whether a hermiticity relation was applied to connect these two ADOs
                    sign_np1 = 1-ksiglm(jnp1,1)                                                     ! Find \bar{sigma} (opposite sign) of new (itrn+1)-th tier ADO
                    eldash_np1 = ksiglm(jnp1,3)                                                         ! Find electronic level of new (itrn+1)-th tier ADO 
                    if (degenerate_levels .eqv. .true.) then
                        el_np1 = eldash_np1
                        if (conj_np1 == 1) then                                                         ! Run this section if a hermiticity relation was applied to relate the two ADOs
                            do ndash = 0,dim_rho-1                                                      ! Loop through all rows/columns of rho_{indjnp1} that could connect to rho_{indjn}(nrow,ncol)
                                indnz2 = rho_sparsity(ncol,ndash,indjnp1)                               ! This part pertains to the d^{\bar{\sigma}}_{m}*rho part of the n,n+1 part of the HEOM.
                                                                                                        ! It evaluates whether 
                                                                                                        ! rho_{indjn}(nrow,ncol) = d^{1-sign_np1}_{el_np1}(nrow,ndash)*rho_{indjnp1}(ncol,ndash)
                                                                                                        ! connects two nonzero ADO elements via a nonzero element of the ann./cre. operator.
                                                                                                        ! Specifically, this line returns -1 if rho_{indjnp1}(ncol,ndash) is a zero element, and its 
                                                                                                        ! nonzero index if it is nonzero. Note that we have had to transpose rho_{indjnp1} because 
                                                                                                        ! of the hermiticity relation.
                                operator_value = d_ops_log(nrow,ndash,el_np1,sign_np1)
                                if ((indnz2 .ne. -1) .and. (operator_value .eqv. .true.)) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                    endif
                                endif

                                indnz2 = rho_sparsity(ndash,nrow,indjnp1)                               ! Does the same, but for the rho*d^{\bar{\sigma}}_{m} part of the n,n+1 part of the HEOM:
                                                                                                        ! rho_{indjn}(nrow,ncol) = rho_{indjnp1}(ndash,nrow)*d^{1-sign_np1}_{el_np1}(ndash,ncol)
                                operator_value = d_ops_log(ndash,ncol,el_np1,sign_np1)
                                if ((indnz2 .ne. -1) .and. (operator_value .eqv. .true.)) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                    endif
                                endif
                            enddo
                        elseif (conj_np1 == 0) then                                                     ! If we did not need to apply a hermiticity relation to relate the two ADOs, this section runs
                                                                                                        ! instead. It does the same thing, but without transposing rho_{indjnp1}
                            do ndash = 0,dim_rho-1
                                indnz2 = rho_sparsity(ndash,ncol,indjnp1)
                                operator_value = d_ops_log(nrow,ndash,el_np1,sign_np1)
                                if ((indnz2 .ne. - 1) .and. (operator_value .eqv. .true.)) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                    endif
                                endif

                                indnz2 = rho_sparsity(nrow,ndash,indjnp1)
                                operator_value = d_ops_log(ndash,ncol,el_np1,sign_np1)
                                if ((indnz2 .ne. - 1) .and. (operator_value .eqv. .true.)) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                    endif
                                endif
                            enddo
                        endif
                    else
                        do el_np1 = 0,nel-1
                            if (conj_np1 == 1) then                                                         ! Run this section if a hermiticity relation was applied to relate the two ADOs
                                do ndash = 0,dim_rho-1                                                      ! Loop through all rows/columns of rho_{indjnp1} that could connect to rho_{indjn}(nrow,ncol)
                                    indnz2 = rho_sparsity(ncol,ndash,indjnp1)                               ! This part pertains to the d^{\bar{\sigma}}_{m}*rho part of the n,n+1 part of the HEOM.
                                                                                                            ! It evaluates whether 
                                                                                                            ! rho_{indjn}(nrow,ncol) = d^{1-sign_np1}_{el_np1}(nrow,ndash)*rho_{indjnp1}(ncol,ndash)
                                                                                                            ! connects two nonzero ADO elements via a nonzero element of the ann./cre. operator.
                                                                                                            ! Specifically, this line returns -1 if rho_{indjnp1}(ncol,ndash) is a zero element, and its 
                                                                                                            ! nonzero index if it is nonzero. Note that we have had to transpose rho_{indjnp1} because 
                                                                                                            ! of the hermiticity relation.
                                    operator_value = d_ops_log(nrow,ndash,el_np1,sign_np1)
                                    if ((indnz2 .ne. -1) .and. (operator_value .eqv. .true.)) then
                                        already_indnz2 = .false.
                                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                                            indnz2_compare = unique_indnz2(itr_indnz2)
                                            if (indnz2 == indnz2_compare) then
                                                already_indnz2 = .true.
                                                exit
                                            endif
                                        enddo
                                        if (already_indnz2 .eqv. .false.) then
                                            count_pairs_total = count_pairs_total + 1
                                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                            unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                        endif
                                    endif

                                    indnz2 = rho_sparsity(ndash,nrow,indjnp1)                               ! Does the same, but for the rho*d^{\bar{\sigma}}_{m} part of the n,n+1 part of the HEOM:
                                                                                                            ! rho_{indjn}(nrow,ncol) = rho_{indjnp1}(ndash,nrow)*d^{1-sign_np1}_{el_np1}(ndash,ncol)
                                    operator_value = d_ops_log(ndash,ncol,el_np1,sign_np1)
                                    if ((indnz2 .ne. -1) .and. (operator_value .eqv. .true.)) then
                                        already_indnz2 = .false.
                                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                                            indnz2_compare = unique_indnz2(itr_indnz2)
                                            if (indnz2 == indnz2_compare) then
                                                already_indnz2 = .true.
                                                exit
                                            endif
                                        enddo
                                        if (already_indnz2 .eqv. .false.) then
                                            count_pairs_total = count_pairs_total + 1
                                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                            unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                        endif
                                    endif
                                enddo
                            elseif (conj_np1 == 0) then                                                     ! If we did not need to apply a hermiticity relation to relate the two ADOs, this section runs
                                                                                                            ! instead. It does the same thing, but without transposing rho_{indjnp1}
                                do ndash = 0,dim_rho-1
                                    indnz2 = rho_sparsity(ndash,ncol,indjnp1)
                                    operator_value = d_ops_log(nrow,ndash,el_np1,sign_np1)
                                    if ((indnz2 .ne. - 1) .and. (operator_value .eqv. .true.)) then
                                        already_indnz2 = .false.
                                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                                            indnz2_compare = unique_indnz2(itr_indnz2)
                                            if (indnz2 == indnz2_compare) then
                                                already_indnz2 = .true.
                                                exit
                                            endif
                                        enddo
                                        if (already_indnz2 .eqv. .false.) then
                                            count_pairs_total = count_pairs_total + 1
                                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                            unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                        endif
                                    endif

                                    indnz2 = rho_sparsity(nrow,ndash,indjnp1)
                                    operator_value = d_ops_log(ndash,ncol,el_np1,sign_np1)
                                    if ((indnz2 .ne. - 1) .and. (operator_value .eqv. .true.)) then
                                        already_indnz2 = .false.
                                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                                            indnz2_compare = unique_indnz2(itr_indnz2)
                                            if (indnz2 == indnz2_compare) then
                                                already_indnz2 = .true.
                                                exit
                                            endif
                                        enddo
                                        if (already_indnz2 .eqv. .false.) then
                                            count_pairs_total = count_pairs_total + 1
                                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                            unique_indnz2(n_indnz2_this_indnz1) = indnz2
                                        endif
                                    endif
                                enddo
                            endif
                        enddo
                    endif
                endif
            enddo
        endif
        if (n_indnz2_this_indnz1 > n_indnz2_this_indnz1_max) then
            n_indnz2_this_indnz1_max = n_indnz2_this_indnz1
        endif
        n_indnz2_prev_indnz1_vec(indnz1+1) = n_indnz2_prev_indnz1_vec(indnz1) + n_indnz2_this_indnz1 + 1
    enddo

    npairs = count_pairs_total + 1 

end subroutine sparse_matrix_elements_a_wbl

! ---------------------------------------------------------------------------
! 
!    GENERATE VALUES AND INFORMATION OF EACH COUPLED PAIR OF ADO ELEMENTS
!                   UNDER THE WIDE-BAND LIMIT 
!
! ---------------------------------------------------------------------------
!
! Same as a sparse_matrix_elements_b, except under the WBL

subroutine sparse_matrix_elements_b_wbl(ksiglm,tier_index,un_ind,index_minus,index_plus,d_ops,gamma_vec,eta_vec,&
                                        rho_sparsity,rho_nonzeros,nnz_elements,dim_rho,len_index_plus,&
                                        len_un_ind,len_index_minus,nmax,nmodes,nel,nsign,npairs,nleads,npoles,&
                                        d_ops_log,ham_log,degenerate_levels,n_indnz2_this_indnz1_max,ham_x,&
                                        len_x_vec,n_indnz2_prev_indnz1_vec,pair_info_row,pair_info_col,pair_values,&
                                        el_lead_couplings)

    implicit none

    integer, intent(in) :: len_index_plus,len_index_minus,nsign,len_un_ind,nnz_elements,npairs      ! Define input variables and arrays
    integer, intent(in) :: nmax,nmodes,nel,dim_rho,npoles,nleads,n_indnz2_this_indnz1_max,len_x_vec
    logical, intent(in) :: degenerate_levels
    double precision, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops
    double precision, intent(in), dimension(0:nleads-1,0:nel-1) :: el_lead_couplings
    complex*16, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:len_x_vec-1,0:4) :: ham_x
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops_log
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: ham_log
    integer, intent(in), dimension(0:nmodes-1,0:3) :: ksiglm
    integer, intent(in), dimension(0:len_un_ind-1) :: tier_index
    integer, intent(in), dimension(0:len_un_ind-1,0:nmax-1) :: un_ind
    integer, intent(in), dimension(0:len_index_plus-1,0:nmodes-1,0:2) :: index_plus
    integer, intent(in), dimension(0:len_index_minus-1,0:nmax-1,0:3) :: index_minus
    complex*16, intent(in),dimension(0:nleads-1,0:nsign-1,0:npoles-1) :: gamma_vec       ! Define complex array output variables
    complex*16, intent(in),dimension(0:nleads-1,0:nsign-1,0:npoles-1) :: eta_vec
    integer, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:len_un_ind-1) :: rho_sparsity
    integer, intent(in), dimension(0:nnz_elements-1,0:2) :: rho_nonzeros
    integer, intent(in), dimension(0:nnz_elements) :: n_indnz2_prev_indnz1_vec

    integer, intent(out), dimension(0:4*npairs-1) :: pair_info_row,pair_info_col
    double precision, intent(out), dimension(0:4*npairs-1,0:len_x_vec-1,0:4) :: pair_values
    
    integer :: nrow,ncol,itrn,indjn,ndash,indnz2,indnz1,leads_np1                                   ! Define other necessary variables
    integer :: itrgamma,leads_gamma,sign_gamma,j_gamma,poles_gamma,poles_np1
    integer :: itrjnm1,jnm1,sign_nm1,el_nm1,indjnm1,leads_nm1,poles_nm1,conj_np1
    integer :: jnp1,sign_np1,el_np1,eldash_np1,indjnp1,count_pairs_total,conj_nm1
    integer :: n_indnz2_this_indnz1,itr_indnz2,indnz2_compare,pair_index,itrx,itrdx
    integer :: si_real,si_imag,nddash,itrel,itrel1,itrel2,itrlead
    integer, dimension(0:n_indnz2_this_indnz1_max-1,0:1) :: unique_indnz2
    complex*16 :: gamma_sum,ham_value
    double precision :: perm_np1,perm_nm1,conj_button_nm1,conj_button_np1,operator_value,wbl_value
    logical :: ham_value_log,operator_value_log,already_indnz2,wbl_value_log
    complex*16, parameter :: ci=(0.d0,1.d0)    
    double precision, parameter :: pi=dacos(-1.d0)

    count_pairs_total = -1
    pair_info_row = -1
    pair_info_col = -1
    pair_values = 0.d0

    do indnz1 = 0,nnz_elements-1                                                                    ! Loop through all nonzero elements (i.e. elements of ADOs on LHS of HEOM)
        unique_indnz2 = -2
        n_indnz2_this_indnz1 = -1
        si_real = 4*n_indnz2_prev_indnz1_vec(indnz1)
        si_imag = si_real + 2*(n_indnz2_prev_indnz1_vec(indnz1+1) - n_indnz2_prev_indnz1_vec(indnz1))
        indjn = rho_nonzeros(indnz1,0)                                                           ! For this nonzero element, find the ADO index and row/column value from rho_nonzeros
        nrow = rho_nonzeros(indnz1,1)
        ncol = rho_nonzeros(indnz1,2)
        itrn = tier_index(indjn)                                                                ! Set itrn = tier of the current ADO
                                                                                                ! Find the tier of this ADO
        do ndash = 0,dim_rho-1                                                                      
            indnz2 = rho_sparsity(ndash,ncol,indjn)
            ham_value_log = ham_log(nrow,ndash)
            if ((indnz2 .ne. -1) .and. (ham_value_log .eqv. .true.)) then                                     ! Testing whether rho_{indjn}(nrow,ncol) = H(nrow,ndash)*rho_{indjn}(ndash,ncol)
                do itrx = 0,(len_x_vec-1)
                    do itrdx = 0,4
                        ham_value = ham_x(nrow,ndash,itrx,itrdx)
                        already_indnz2 = .false.
                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                            indnz2_compare = unique_indnz2(itr_indnz2,0)
                            if (indnz2 == indnz2_compare) then
                                already_indnz2 = .true.
                                exit
                            endif
                        enddo
                        if (already_indnz2 .eqv. .false.) then
                            count_pairs_total = count_pairs_total + 1
                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                            unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                            unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                            pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                            pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                            pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                            pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                            pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                            pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                            pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                            pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                            pair_values(si_real+2*n_indnz2_this_indnz1,itrx,itrdx) = dble(-ci*ham_value)
                            pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx,itrdx) = -aimag(-ci*ham_value)
                            pair_values(si_imag+2*n_indnz2_this_indnz1,itrx,itrdx) = aimag(-ci*ham_value)
                            pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx,itrdx) = dble(-ci*ham_value)
                        else
                            pair_index = unique_indnz2(itr_indnz2,1)
                            pair_values(si_real+2*pair_index,itrx,itrdx) = pair_values(si_real+2*pair_index,itrx,itrdx) + dble(-ci*ham_value)
                            pair_values(si_real+2*pair_index+1,itrx,itrdx) = pair_values(si_real+2*pair_index+1,itrx,itrdx) - aimag(-ci*ham_value)
                            pair_values(si_imag+2*pair_index,itrx,itrdx) = pair_values(si_imag+2*pair_index,itrx,itrdx) + aimag(-ci*ham_value)
                            pair_values(si_imag+2*pair_index+1,itrx,itrdx) = pair_values(si_imag+2*pair_index+1,itrx,itrdx) + dble(-ci*ham_value)
                        endif
                    enddo
                enddo
            endif
    
            indnz2 = rho_sparsity(nrow,ndash,indjn)
            ham_value_log = ham_log(ndash,ncol)
            if ((indnz2 .ne. -1) .and. (ham_value_log .eqv. .true.)) then
                do itrx = 0,(len_x_vec-1)
                    do itrdx = 0,4
                        ham_value = ham_x(ndash,ncol,itrx,itrdx)
                        already_indnz2 = .false.
                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                            indnz2_compare = unique_indnz2(itr_indnz2,0)
                            if (indnz2 == indnz2_compare) then
                                already_indnz2 = .true.
                                exit
                            endif
                        enddo
                        if (already_indnz2 .eqv. .false.) then
                            count_pairs_total = count_pairs_total + 1
                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                            unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                            unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                            pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                            pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                            pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                            pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                            pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                            pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                            pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                            pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                            pair_values(si_real+2*n_indnz2_this_indnz1,itrx,itrdx) = dble(ci*ham_value)
                            pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx,itrdx) = -aimag(ci*ham_value)
                            pair_values(si_imag+2*n_indnz2_this_indnz1,itrx,itrdx) = aimag(ci*ham_value)
                            pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx,itrdx) = dble(ci*ham_value)
                        else
                            pair_index = unique_indnz2(itr_indnz2,1)
                            pair_values(si_real+2*pair_index,itrx,itrdx) = pair_values(si_real+2*pair_index,itrx,itrdx) + dble(ci*ham_value)
                            pair_values(si_real+2*pair_index+1,itrx,itrdx) = pair_values(si_real+2*pair_index+1,itrx,itrdx) - aimag(ci*ham_value)
                            pair_values(si_imag+2*pair_index,itrx,itrdx) = pair_values(si_imag+2*pair_index,itrx,itrdx) + aimag(ci*ham_value)
                            pair_values(si_imag+2*pair_index+1,itrx,itrdx) = pair_values(si_imag+2*pair_index+1,itrx,itrdx) + dble(ci*ham_value)
                        endif
                    enddo
                enddo
            endif
        enddo

        if (degenerate_levels .eqv. .true.) then 
            do ndash = 0,dim_rho-1
                do nddash = 0,dim_rho-1
                    indnz2 = rho_sparsity(ndash,nddash,indjn)
                    if (indnz2 .ne. -1) then
                        do itrel = 0,nel-1
                            wbl_value = d_ops(nrow,ndash,itrel,0)*d_ops(nddash,ncol,itrel,1)
                            wbl_value_log = d_ops_log(nrow,ndash,itrel,0) .and. d_ops_log(nddash,ncol,itrel,1)
                            if (wbl_value_log .eqv. .true.) then
                                already_indnz2 = .false.
                                do itr_indnz2 = 0,n_indnz2_this_indnz1
                                    indnz2_compare = unique_indnz2(itr_indnz2,0)
                                    if (indnz2 == indnz2_compare) then
                                        already_indnz2 = .true.
                                        exit
                                    endif
                                enddo
                                if (already_indnz2 .eqv. .false.) then
                                    count_pairs_total = count_pairs_total + 1
                                    n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                    unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                    unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                    pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                    pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                    pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                    pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                    pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                    pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                    pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                    pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                    pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = 0.d0
                                    pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = 0.d0
                                    pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = 0.d0
                                    pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = 0.d0
                                    do itrlead = 0,nleads-1
                                        pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = &
                                            pair_values(si_real+2*n_indnz2_this_indnz1,:,:) &
                                            - pi*((-1.d0)**(itrn+1))*wbl_value*(el_lead_couplings(itrlead,itrel)**2)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = &
                                        pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) &
                                            - pi*((-1.d0)**(itrn+1))*wbl_value*(el_lead_couplings(itrlead,itrel)**2)
                                    enddo
                                else
                                    pair_index = unique_indnz2(itr_indnz2,1)
                                    do itrlead = 0,nleads-1
                                        pair_values(si_real+2*pair_index,:,:) = &
                                            pair_values(si_real+2*pair_index,:,:) &
                                            - pi*((-1.d0)**(itrn+1))*wbl_value*(el_lead_couplings(itrlead,itrel)**2)
                                        pair_values(si_imag+2*pair_index+1,:,:) = &
                                            pair_values(si_imag+2*pair_index+1,:,:)&
                                            - pi*((-1.d0)**(itrn+1))*wbl_value*(el_lead_couplings(itrlead,itrel)**2)
                                    enddo
                                endif
                            endif
                            wbl_value = d_ops(nrow,ndash,itrel,1)*d_ops(nddash,ncol,itrel,0)
                            wbl_value_log = d_ops_log(nrow,ndash,itrel,1) .and. d_ops_log(nddash,ncol,itrel,0)
                            if (wbl_value_log .eqv. .true.) then
                                already_indnz2 = .false.
                                do itr_indnz2 = 0,n_indnz2_this_indnz1
                                    indnz2_compare = unique_indnz2(itr_indnz2,0)
                                    if (indnz2 == indnz2_compare) then
                                        already_indnz2 = .true.
                                        exit
                                    endif
                                enddo
                                if (already_indnz2 .eqv. .false.) then
                                    count_pairs_total = count_pairs_total + 1
                                    n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                    unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                    unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                    pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                    pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                    pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                    pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                    pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                    pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                    pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                    pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                    pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = 0.d0
                                    pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = 0.d0
                                    pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = 0.d0
                                    pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = 0.d0
                                    do itrlead = 0,nleads-1
                                        pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = &
                                            pair_values(si_real+2*n_indnz2_this_indnz1,:,:) &
                                            - pi*((-1.d0)**(itrn+1))*wbl_value*(el_lead_couplings(itrlead,itrel)**2)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = &
                                        pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) &
                                            - pi*((-1.d0)**(itrn+1))*wbl_value*(el_lead_couplings(itrlead,itrel)**2)
                                    enddo
                                else
                                    pair_index = unique_indnz2(itr_indnz2,1)
                                    do itrlead = 0,nleads-1
                                        pair_values(si_real+2*pair_index,:,:) = &
                                                pair_values(si_real+2*pair_index,:,:) &
                                                - pi*((-1.d0)**(itrn+1))*wbl_value*(el_lead_couplings(itrlead,itrel)**2)
                                        pair_values(si_imag+2*pair_index+1,:,:) = pair_values(si_imag+2*pair_index+1,:,:)&
                                                - pi*((-1.d0)**(itrn+1))*wbl_value*(el_lead_couplings(itrlead,itrel)**2)
                                    enddo
                                endif
                            endif
                        enddo
                    endif
                enddo
            enddo
        else
            do ndash = 0,dim_rho-1
                do nddash = 0,dim_rho-1
                    indnz2 = rho_sparsity(ndash,nddash,indjn)
                    if (indnz2 .ne. -1) then
                        do itrel1 = 0,nel-1
                            do itrel2 = 0,nel-1
                                wbl_value = d_ops(nrow,ndash,itrel1,0)*d_ops(nddash,ncol,itrel2,1)
                                wbl_value_log = d_ops_log(nrow,ndash,itrel1,0) .and. &
                                                d_ops_log(nddash,ncol,itrel2,1)
                                if (wbl_value_log .eqv. .true.) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2,0)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                        unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = 0.d0
                                        pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = 0.d0
                                        pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = 0.d0
                                        pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = 0.d0
                                        do itrlead = 0,nleads-1
                                            pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = &
                                                pair_values(si_real+2*n_indnz2_this_indnz1,:,:) &
                                                - pi*((-1.d0)**(itrn+1))*wbl_value*&
                                                el_lead_couplings(itrlead,itrel1)*el_lead_couplings(itrlead,itrel2)
                                            pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = &
                                            pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) &
                                                - pi*((-1.d0)**(itrn+1))*wbl_value*&
                                                el_lead_couplings(itrlead,itrel1)*el_lead_couplings(itrlead,itrel2)
                                        enddo
                                    else
                                        pair_index = unique_indnz2(itr_indnz2,1)
                                        do itrlead = 0,nleads-1
                                            pair_values(si_real+2*pair_index,:,:) = &
                                                    pair_values(si_real+2*pair_index,:,:) &
                                                    - pi*((-1.d0)**(itrn+1))*wbl_value*&
                                                    el_lead_couplings(itrlead,itrel1)*el_lead_couplings(itrlead,itrel2)
                                            pair_values(si_imag+2*pair_index+1,:,:) = pair_values(si_imag+2*pair_index+1,:,:)&
                                                    - pi*((-1.d0)**(itrn+1))*wbl_value*&
                                                    el_lead_couplings(itrlead,itrel1)*el_lead_couplings(itrlead,itrel2)
                                        enddo
                                    endif
                                endif
                                wbl_value = d_ops(nrow,ndash,itrel1,1)*d_ops(nddash,ncol,itrel2,0)
                                wbl_value_log = d_ops_log(nrow,ndash,itrel1,1) .and. &
                                                d_ops_log(nddash,ncol,itrel2,0)
                                if (wbl_value_log .eqv. .true.) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2,0)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                        unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = 0.d0
                                        pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = 0.d0
                                        pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = 0.d0
                                        pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = 0.d0
                                        do itrlead = 0,nleads-1
                                            pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = &
                                                pair_values(si_real+2*n_indnz2_this_indnz1,:,:) &
                                                - pi*((-1.d0)**(itrn+1))*wbl_value*&
                                                el_lead_couplings(itrlead,itrel1)*el_lead_couplings(itrlead,itrel2)
                                            pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = &
                                            pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) &
                                                - pi*((-1.d0)**(itrn+1))*wbl_value*&
                                                el_lead_couplings(itrlead,itrel1)*el_lead_couplings(itrlead,itrel2)
                                        enddo
                                    else
                                        pair_index = unique_indnz2(itr_indnz2,1)
                                        do itrlead = 0,nleads-1
                                            pair_values(si_real+2*pair_index,:,:) = &
                                                    pair_values(si_real+2*pair_index,:,:) &
                                                    - pi*((-1.d0)**(itrn+1))*wbl_value*&
                                                    el_lead_couplings(itrlead,itrel1)*el_lead_couplings(itrlead,itrel2)
                                            pair_values(si_imag+2*pair_index+1,:,:) = pair_values(si_imag+2*pair_index+1,:,:)&
                                                    - pi*((-1.d0)**(itrn+1))*wbl_value*&
                                                    el_lead_couplings(itrlead,itrel1)*el_lead_couplings(itrlead,itrel2)
                                        enddo
                                    endif
                                endif
                            enddo
                        enddo
                    endif
                enddo
            enddo
        endif

        indnz2 = indnz1
        already_indnz2 = .false.
        do itr_indnz2 = 0,n_indnz2_this_indnz1
            indnz2_compare = unique_indnz2(itr_indnz2,0)
            if (indnz2 == indnz2_compare) then
                already_indnz2 = .true.
                exit
            endif
        enddo
        if (already_indnz2 .eqv. .false.) then
            count_pairs_total = count_pairs_total + 1
            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
            unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
            unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
            pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
            pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
            pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
            pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
            pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
            pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
            pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
            pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
            gamma_sum = 0.d0
            if (itrn > 0) then
                do itrgamma = 0,itrn-1                                                                  ! Loop through all mode indices of the ADO of this ADO element
                    j_gamma = un_ind(indjn,itrgamma)                                                    ! Find the mode corresponding to this mode index
                    leads_gamma = ksiglm(j_gamma,0)                                                     ! Find the lead corresponding to this mode
                    sign_gamma = ksiglm(j_gamma,1)                                                      ! Find the sigma (sign) corresponding to this mode
                    poles_gamma = ksiglm(j_gamma,2)                                                     ! Find the Pade pole corresponding to this mode
                    gamma_sum = gamma_sum + gamma_vec(leads_gamma,sign_gamma,poles_gamma)               ! Update the sum over gammas with the gamma corresponding to this mode
                enddo
            endif
            pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = dble(-gamma_sum)
            pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = -aimag(-gamma_sum)
            pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = aimag(-gamma_sum)
            pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = dble(-gamma_sum)
            do itrlead = 0,nleads-1
                do itrel = 0,nel-1
                    pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = &
                            pair_values(si_real+2*n_indnz2_this_indnz1,:,:) &
                            - pi*((el_lead_couplings(itrlead,itrel)**2))
                    pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = &
                            pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) &
                            - pi*((el_lead_couplings(itrlead,itrel)**2))
                enddo
            enddo
        else
            pair_index = unique_indnz2(itr_indnz2,1)
            gamma_sum = 0.d0
            if (itrn > 0) then
                do itrgamma = 0,itrn-1                                                                  ! Loop through all mode indices of the ADO of this ADO element
                    j_gamma = un_ind(indjn,itrgamma)                                                    ! Find the mode corresponding to this mode index
                    leads_gamma = ksiglm(j_gamma,0)                                                     ! Find the lead corresponding to this mode
                    sign_gamma = ksiglm(j_gamma,1)                                                      ! Find the sigma (sign) corresponding to this mode
                    poles_gamma = ksiglm(j_gamma,2)                                                     ! Find the Pade pole corresponding to this mode
                    gamma_sum = gamma_sum + gamma_vec(leads_gamma,sign_gamma,poles_gamma)               ! Update the sum over gammas with the gamma corresponding to this mode
                enddo
            endif
            pair_values(si_real+2*pair_index,:,:) = pair_values(si_real+2*pair_index,:,:) + dble(-gamma_sum)
            pair_values(si_real+2*pair_index+1,:,:) = pair_values(si_real+2*pair_index+1,:,:) - aimag(-gamma_sum)
            pair_values(si_imag+2*pair_index,:,:) = pair_values(si_imag+2*pair_index,:,:) + aimag(-gamma_sum)
            pair_values(si_imag+2*pair_index+1,:,:) = pair_values(si_imag+2*pair_index+1,:,:) + dble(-gamma_sum)
            do itrlead = 0,nleads-1
                do itrel = 0,nel-1
                    pair_values(si_real+2*pair_index,:,:) = pair_values(si_real+2*pair_index,:,:) &
                            - pi*((el_lead_couplings(itrlead,itrel)**2))
                    pair_values(si_imag+2*pair_index+1,:,:) = pair_values(si_imag+2*pair_index+1,:,:) + &
                            - pi*((el_lead_couplings(itrlead,itrel)**2))
                enddo
            enddo
        endif

        if (itrn > 0) then
            do itrjnm1 = 0,itrn-1                                                                   ! This section assesses couplings between ADO elements of tier n and tier n-1. First, loop
                                                                                                    ! through mode indices to remove from ADO of current element
                jnm1 = index_minus(indjn,itrjnm1,0)                                                 ! Find mode being removed in this loop
                leads_nm1 = ksiglm(jnm1,0)                                                          ! Find lead index of mode being removed
                sign_nm1 = ksiglm(jnm1,1)                                                           ! Find sigma index (sign) of mode being removed
                poles_nm1 = ksiglm(jnm1,2)                                                          ! Find Pade pole of mode being removed
                el_nm1 = ksiglm(jnm1,3)                                                             ! Find electronic level index of mode being removed
                indjnm1 = index_minus(indjn,itrjnm1,1)                                              ! Find index of new (n-1)-th tier ADO created after this mode is removed
                perm_nm1 = (-1.d0)**(index_minus(indjn,itrjnm1,2))                                  ! Calculate permutation prefactor required to connect new ADO to current ADO
                conj_nm1 = index_minus(indjn,itrjnm1,3)                                             ! Determine whether a hermiticity relation was applied to connect new ADO to current ADO
                if (conj_nm1 == 1) then                                                             ! If a hermiticity relaiton was applied, we need to use the transpose of the new ADO
                    conj_button_nm1 = (-1.d0)**(floor((dble(itrn)-1.d0)/2.d0))                      ! Calculate hermiticity prefactor
                    do ndash = 0,dim_rho-1                                                          ! Loop through all rows/columns that could possibly connect the two ADOs:
                                                                                                    ! rho_{indjn}(nrow,ncol) = eta_{l}*d^{\sigma}_{m}(nrow,ndash)*rho_{indjnm1}(ncol,ndash)
                        indnz2 = rho_sparsity(ncol,ndash,indjnm1)                                   ! Find sparsity index of RHS ADO element: rho_{indjnm1}(ncol,ndash)
                        operator_value = d_ops(nrow,ndash,el_nm1,sign_nm1)                          ! Find corresponding element of ann./cre. operator
                        operator_value_log = d_ops_log(nrow,ndash,el_nm1,sign_nm1)
                        if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then          ! If this ADO element and ann./cre. element are both nonzero,
                            already_indnz2 = .false.
                            do itr_indnz2 = 0,n_indnz2_this_indnz1
                                indnz2_compare = unique_indnz2(itr_indnz2,0)
                                if (indnz2 == indnz2_compare) then
                                    already_indnz2 = .true.
                                    exit
                                endif
                            enddo
                            if (already_indnz2 .eqv. .false.) then
                                count_pairs_total = count_pairs_total + 1
                                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = &
                                                            dble(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                            operator_value)
                                pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = &
                                                            aimag(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                            operator_value)
                                pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = &
                                                            aimag(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                            operator_value)
                                pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = &
                                                            -dble(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                            operator_value)
                            else
                                pair_index = unique_indnz2(itr_indnz2,1)
                                pair_values(si_real+2*pair_index,:,:) = &
                                                            pair_values(si_real+2*pair_index,:,:) + &
                                                            dble(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                            operator_value)
                                pair_values(si_real+2*pair_index+1,:,:) = &
                                                            pair_values(si_real+2*pair_index+1,:,:) + &
                                                            aimag(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                            operator_value)
                                pair_values(si_imag+2*pair_index,:,:) = &
                                                            pair_values(si_imag+2*pair_index,:,:) + &
                                                            aimag(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                            operator_value)
                                pair_values(si_imag+2*pair_index+1,:,:) = &
                                                            pair_values(si_imag+2*pair_index+1,:,:) - & 
                                                            dble(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                            operator_value)
                            endif
                        endif

                        indnz2 = rho_sparsity(ndash,nrow,indjnm1)                                   ! Do the same, but for 
                                                                                                    ! rho_{indjn}(nrow,ncol) = eta^{*}_{l}*rho_{indjnm1}(ndash,nrow)*d^{\sigma}_{m}(ndash,ncol)
                        operator_value = d_ops(ndash,ncol,el_nm1,sign_nm1)
                        operator_value_log = d_ops_log(ndash,ncol,el_nm1,sign_nm1)
                        if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                            already_indnz2 = .false.
                            do itr_indnz2 = 0,n_indnz2_this_indnz1
                                indnz2_compare = unique_indnz2(itr_indnz2,0)
                                if (indnz2 == indnz2_compare) then
                                    already_indnz2 = .true.
                                    exit
                                endif
                            enddo
                            if (already_indnz2 .eqv. .false.) then
                                count_pairs_total = count_pairs_total + 1
                                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = &
                                                            dble(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conj_button_nm1*conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                            operator_value)
                                pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = &
                                                            aimag(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conj_button_nm1*conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                            operator_value)
                                pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = &
                                                            aimag(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conj_button_nm1*conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                            operator_value)
                                pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = &
                                                            -dble(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conj_button_nm1*conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                            operator_value)
                            else
                                pair_index = unique_indnz2(itr_indnz2,1)
                                pair_values(si_real+2*pair_index,:,:) = pair_values(si_real+2*pair_index,:,:) + &
                                                            dble(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conj_button_nm1*&
                                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                            operator_value)
                                pair_values(si_real+2*pair_index+1,:,:) = pair_values(si_real+2*pair_index+1,:,:) + &
                                                            aimag(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conj_button_nm1*&
                                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                            operator_value)
                                pair_values(si_imag+2*pair_index,:,:) = pair_values(si_imag+2*pair_index,:,:) + &
                                                            aimag(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conj_button_nm1*&
                                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                            operator_value)
                                pair_values(si_imag+2*pair_index+1,:,:) = pair_values(si_imag+2*pair_index+1,:,:) - & 
                                                            dble(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conj_button_nm1*&
                                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                            operator_value)
                            endif
                        endif
                    enddo
                elseif (conj_nm1 == 0) then                                                         ! Do the same process as above, but run only if no hermiticity relation is applied and we do 
                                                                                                    ! not need to work with the conjugate-transposed ADO
                    do ndash = 0,dim_rho-1
                        indnz2 = rho_sparsity(ndash,ncol,indjnm1)
                        operator_value = d_ops(nrow,ndash,el_nm1,sign_nm1)
                        operator_value_log = d_ops_log(nrow,ndash,el_nm1,sign_nm1)
                        if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                            already_indnz2 = .false.
                            do itr_indnz2 = 0,n_indnz2_this_indnz1
                                indnz2_compare = unique_indnz2(itr_indnz2,0)
                                if (indnz2 == indnz2_compare) then
                                    already_indnz2 = .true.
                                    exit
                                endif
                            enddo
                            if (already_indnz2 .eqv. .false.) then
                                count_pairs_total = count_pairs_total + 1
                                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = &
                                                        dble(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                        eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                        operator_value)                                                
                                pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = &
                                                        -aimag(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                        eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                        operator_value)
                                pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = &
                                                        aimag(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                        eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                        operator_value)
                                pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = &
                                                        dble(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                        eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                        operator_value)
                            else
                                pair_index = unique_indnz2(itr_indnz2,1)
                                pair_values(si_real+2*pair_index,:,:) = &
                                                        pair_values(si_real+2*pair_index,:,:) + &
                                                        dble(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                        eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                        operator_value)
                                pair_values(si_real+2*pair_index+1,:,:) = &
                                                        pair_values(si_real+2*pair_index+1,:,:) - &
                                                        aimag(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                        eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                        operator_value)
                                pair_values(si_imag+2*pair_index,:,:) = &
                                                        pair_values(si_imag+2*pair_index,:,:) + &
                                                        aimag(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                        eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                        operator_value)
                                pair_values(si_imag+2*pair_index+1,:,:) = &
                                                        pair_values(si_imag+2*pair_index+1,:,:) + & 
                                                        dble(-ci*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                        eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                        operator_value)
                            endif
                        endif
                        indnz2 = rho_sparsity(nrow,ndash,indjnm1)
                        operator_value = d_ops(ndash,ncol,el_nm1,sign_nm1)
                        operator_value_log = d_ops_log(ndash,ncol,el_nm1,sign_nm1)
                        if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                            already_indnz2 = .false.
                            do itr_indnz2 = 0,n_indnz2_this_indnz1
                                indnz2_compare = unique_indnz2(itr_indnz2,0)
                                if (indnz2 == indnz2_compare) then
                                    already_indnz2 = .true.
                                    exit
                                endif
                            enddo
                            if (already_indnz2 .eqv. .false.) then
                                count_pairs_total = count_pairs_total + 1
                                n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = &
                                                            dble(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                            operator_value)
                                pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = &
                                                            -aimag(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                            operator_value)
                                pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = &
                                                            aimag(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                            operator_value)
                                pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = &
                                                            dble(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                            operator_value)
                            else
                                pair_index = unique_indnz2(itr_indnz2,1)
                                pair_values(si_real+2*pair_index,:,:) = pair_values(si_real+2*pair_index,:,:) + &
                                                            dble(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                            operator_value)
                                pair_values(si_real+2*pair_index+1,:,:) = pair_values(si_real+2*pair_index+1,:,:) - &
                                                            aimag(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                            operator_value)
                                pair_values(si_imag+2*pair_index,:,:) = pair_values(si_imag+2*pair_index,:,:) + &
                                                            aimag(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                            operator_value)
                                pair_values(si_imag+2*pair_index+1,:,:) = pair_values(si_imag+2*pair_index+1,:,:) + & 
                                                            dble(ci*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                            operator_value)
                            endif
                        endif
                    enddo
                endif
            enddo
        endif

        if (itrn < nmax) then                                                                       ! Now we connect the nth and (n+1)th tier ADOs
            do jnp1 = 0,nmodes-1                                                                    ! Loop through all possible modes to add to the current ADO
                indjnp1 = index_plus(indjn,jnp1,0)                                                  ! Find index of new (itrn+1)-th tier ADO created by adding mode jnp1 to current ADO
                if (indjnp1 .ne. -1) then
                    perm_np1 = (-1.d0)**(index_plus(indjn,jnp1,1))                                  ! Calculate the permutation prefactor generated when connecting these two ADOs
                    conj_np1 = index_plus(indjn,jnp1,2)                                             ! Determine whether a hermiticity relation was applied to connect these two ADOs
                    leads_np1 = ksiglm(jnp1,0)                                                      ! Find lead index of mode being added
                    sign_np1 = 1-ksiglm(jnp1,1)                                                     ! Find \bar{\sigma} (1-sign) index of mode being added
                    poles_np1 = ksiglm(jnp1,2)                                                      ! Find Pade pole of mode being added
                    eldash_np1 = ksiglm(jnp1,3)                                                         ! Find electronic level of mode being added
                    if (degenerate_levels .eqv. .true.) then
                        el_np1 = eldash_np1
                        if (conj_np1 == 1) then                                                         ! Run this code if a hermiticity relation was applied
                            conj_button_np1 = (-1.d0)**(floor((dble(itrn)+1.d0)/2.d0))                         ! Calculate the hermiticity prefactor generated when connecting these two ADOs
                            do ndash = 0,dim_rho-1                                                      ! Loop through all rows/columns of rho_{indjnp1} that could connect to rho_{indjn}(nrow,ncol):
                                                                                                        ! rho_{indjn}(nrow,ncol) = d^{\bar{\sigma}}_{m}(nrow,ndash)*rho_{indjnp1}(ncol,ndash)
                                indnz2 = rho_sparsity(ncol,ndash,indjnp1)
                                operator_value = d_ops(nrow,ndash,el_np1,sign_np1)
                                operator_value_log = d_ops_log(nrow,ndash,el_np1,sign_np1)
                                if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2,0)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                        unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = &
                                                                dble(-ci*perm_np1*conj_button_np1*&
                                                                operator_value)
                                        pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = &
                                                                aimag(-ci*perm_np1*conj_button_np1*&
                                                                operator_value)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = &
                                                                aimag(-ci*perm_np1*conj_button_np1*&
                                                                operator_value)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = &
                                                                -dble(-ci*perm_np1*conj_button_np1*&
                                                                operator_value)
                                    else
                                        pair_index = unique_indnz2(itr_indnz2,1)
                                        pair_values(si_real+2*pair_index,:,:) = &
                                                                pair_values(si_real+2*pair_index,:,:) + &
                                                                dble(-ci*perm_np1*conj_button_np1*&
                                                                operator_value)
                                        pair_values(si_real+2*pair_index+1,:,:) = &
                                                                pair_values(si_real+2*pair_index+1,:,:) + &
                                                                aimag(-ci*perm_np1*conj_button_np1*&
                                                                operator_value)
                                        pair_values(si_imag+2*pair_index,:,:) = &
                                                                pair_values(si_imag+2*pair_index,:,:) + &
                                                                aimag(-ci*perm_np1*conj_button_np1*&
                                                                operator_value)
                                        pair_values(si_imag+2*pair_index+1,:,:) = &
                                                                pair_values(si_imag+2*pair_index+1,:,:) - & 
                                                                dble(-ci*perm_np1*conj_button_np1*&
                                                                operator_value)
                                    endif
                                endif
                                indnz2 = rho_sparsity(ndash,nrow,indjnp1)                               ! Do the same, but for 
                                                                                                        ! rho_{indjn}(nrow,ncol) = rho_{indjnp1}(ndash,nrow)*d^{\sigma}_{m}(ndash,ncol)
                                operator_value = d_ops(ndash,ncol,el_np1,sign_np1)
                                operator_value_log = d_ops_log(ndash,ncol,el_np1,sign_np1)
                                if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2,0)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                        unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = &
                                                                dble(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                operator_value)
                                        pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = &
                                                                aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                operator_value)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = &
                                                                aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                operator_value)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = &
                                                                -dble(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                operator_value)
                                    else
                                        pair_index = unique_indnz2(itr_indnz2,1)
                                        pair_values(si_real+2*pair_index,:,:) = &
                                                                pair_values(si_real+2*pair_index,:,:) + &
                                                                dble(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                operator_value)
                                        pair_values(si_real+2*pair_index+1,:,:) = &
                                                                pair_values(si_real+2*pair_index+1,:,:) + &
                                                                aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                operator_value)
                                        pair_values(si_imag+2*pair_index,:,:) = &
                                                                pair_values(si_imag+2*pair_index,:,:) + &
                                                                aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                operator_value)
                                        pair_values(si_imag+2*pair_index+1,:,:) = &
                                                                pair_values(si_imag+2*pair_index+1,:,:) - & 
                                                                dble(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                operator_value)
                                    endif
                                endif
                            enddo
                        elseif (conj_np1 == 0) then                                                     ! Run this part instead if hermiticity relation is not applied to connect ADOs 
                                                                                                        ! and we do not need to take the conjugate transpose 
                            do ndash = 0,dim_rho-1
                                indnz2 = rho_sparsity(ndash,ncol,indjnp1)
                                operator_value = d_ops(nrow,ndash,el_np1,sign_np1)
                                operator_value_log = d_ops_log(nrow,ndash,el_np1,sign_np1)
                                if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2,0)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                        unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = &
                                                                dble(-ci*perm_np1*operator_value)
                                        pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = &
                                                                -aimag(-ci*perm_np1*operator_value)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = &
                                                                aimag(-ci*perm_np1*operator_value)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = &
                                                                dble(-ci*perm_np1*operator_value)
                                    else
                                        pair_index = unique_indnz2(itr_indnz2,1)
                                        pair_values(si_real+2*pair_index,:,:) = &
                                                                pair_values(si_real+2*pair_index,:,:) + &
                                                                dble(-ci*perm_np1*operator_value)
                                        pair_values(si_real+2*pair_index+1,:,:) = &
                                                                pair_values(si_real+2*pair_index+1,:,:) - &
                                                                aimag(-ci*perm_np1*operator_value)
                                        pair_values(si_imag+2*pair_index,:,:) = &
                                                                pair_values(si_imag+2*pair_index,:,:) + &
                                                                aimag(-ci*perm_np1*operator_value)
                                        pair_values(si_imag+2*pair_index+1,:,:) = &
                                                                pair_values(si_imag+2*pair_index+1,:,:) + & 
                                                                dble(-ci*perm_np1*operator_value)
                                    endif
                                endif
                    
                                indnz2 = rho_sparsity(nrow,ndash,indjnp1)
                                operator_value = d_ops(ndash,ncol,el_np1,sign_np1)
                                operator_value_log = d_ops_log(ndash,ncol,el_np1,sign_np1)
                                if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                                    already_indnz2 = .false.
                                    do itr_indnz2 = 0,n_indnz2_this_indnz1
                                        indnz2_compare = unique_indnz2(itr_indnz2,0)
                                        if (indnz2 == indnz2_compare) then
                                            already_indnz2 = .true.
                                            exit
                                        endif
                                    enddo
                                    if (already_indnz2 .eqv. .false.) then
                                        count_pairs_total = count_pairs_total + 1
                                        n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                        unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                        unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                        pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                        pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                        pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                        pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = &
                                                                dble(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                                operator_value)
                                        pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = &
                                                                -aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                                operator_value)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = &
                                                                aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                                operator_value)
                                        pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = &
                                                                dble(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                                operator_value)
                                    else
                                        pair_index = unique_indnz2(itr_indnz2,1)
                                        pair_values(si_real+2*pair_index,:,:) = &
                                                                pair_values(si_real+2*pair_index,:,:) + &
                                                                dble(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                                operator_value)
                                        pair_values(si_real+2*pair_index+1,:,:) = &
                                                                pair_values(si_real+2*pair_index+1,:,:) - &
                                                                aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                                operator_value)
                                        pair_values(si_imag+2*pair_index,:,:) = &
                                                                pair_values(si_imag+2*pair_index,:,:) + &
                                                                aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                                operator_value)
                                        pair_values(si_imag+2*pair_index+1,:,:) = &
                                                                pair_values(si_imag+2*pair_index+1,:,:) + & 
                                                                dble(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                                operator_value)
                                    endif
                                endif
                            enddo
                        endif
                    else
                        do el_np1 = 0,nel-1
                            if (conj_np1 == 1) then                                                         ! Run this code if a hermiticity relation was applied
                                conj_button_np1 = (-1.d0)**(floor((dble(itrn)+1.d0)/2.d0))                         ! Calculate the hermiticity prefactor generated when connecting these two ADOs
                                do ndash = 0,dim_rho-1                                                      ! Loop through all rows/columns of rho_{indjnp1} that could connect to rho_{indjn}(nrow,ncol):
                                                                                                            ! rho_{indjn}(nrow,ncol) = d^{\bar{\sigma}}_{m}(nrow,ndash)*rho_{indjnp1}(ncol,ndash)
                                    indnz2 = rho_sparsity(ncol,ndash,indjnp1)
                                    operator_value = d_ops(nrow,ndash,el_np1,sign_np1)
                                    operator_value_log = d_ops_log(nrow,ndash,el_np1,sign_np1)
                                    if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                                        already_indnz2 = .false.
                                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                                            indnz2_compare = unique_indnz2(itr_indnz2,0)
                                            if (indnz2 == indnz2_compare) then
                                                already_indnz2 = .true.
                                                exit
                                            endif
                                        enddo
                                        if (already_indnz2 .eqv. .false.) then
                                            count_pairs_total = count_pairs_total + 1
                                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                            unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                            unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                            pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                            pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                            pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                            pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                            pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                            pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                            pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                            pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                            pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = &
                                                                    dble(-ci*perm_np1*conj_button_np1*&
                                                                    operator_value)
                                            pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = &
                                                                    aimag(-ci*perm_np1*conj_button_np1*&
                                                                    operator_value)
                                            pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = &
                                                                    aimag(-ci*perm_np1*conj_button_np1*&
                                                                    operator_value)
                                            pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = &
                                                                    -dble(-ci*perm_np1*conj_button_np1*&
                                                                    operator_value)
                                        else
                                            pair_index = unique_indnz2(itr_indnz2,1)
                                            pair_values(si_real+2*pair_index,:,:) = &
                                                                    pair_values(si_real+2*pair_index,:,:) + &
                                                                    dble(-ci*perm_np1*conj_button_np1*&
                                                                    operator_value)
                                            pair_values(si_real+2*pair_index+1,:,:) = &
                                                                    pair_values(si_real+2*pair_index+1,:,:) + &
                                                                    aimag(-ci*perm_np1*conj_button_np1*&
                                                                    operator_value)
                                            pair_values(si_imag+2*pair_index,:,:) = &
                                                                    pair_values(si_imag+2*pair_index,:,:) + &
                                                                    aimag(-ci*perm_np1*conj_button_np1*&
                                                                    operator_value)
                                            pair_values(si_imag+2*pair_index+1,:,:) = &
                                                                    pair_values(si_imag+2*pair_index+1,:,:) - & 
                                                                    dble(-ci*perm_np1*conj_button_np1*&
                                                                    operator_value)
                                        endif
                                    endif
                                    indnz2 = rho_sparsity(ndash,nrow,indjnp1)                               ! Do the same, but for 
                                                                                                            ! rho_{indjn}(nrow,ncol) = rho_{indjnp1}(ndash,nrow)*d^{\sigma}_{m}(ndash,ncol)
                                    operator_value = d_ops(ndash,ncol,el_np1,sign_np1)
                                    operator_value_log = d_ops_log(ndash,ncol,el_np1,sign_np1)
                                    if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                                        already_indnz2 = .false.
                                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                                            indnz2_compare = unique_indnz2(itr_indnz2,0)
                                            if (indnz2 == indnz2_compare) then
                                                already_indnz2 = .true.
                                                exit
                                            endif
                                        enddo
                                        if (already_indnz2 .eqv. .false.) then
                                            count_pairs_total = count_pairs_total + 1
                                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                            unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                            unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                            pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                            pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                            pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                            pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                            pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                            pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                            pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                            pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                            pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = &
                                                                    dble(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                    operator_value)
                                            pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = &
                                                                    aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                    operator_value)
                                            pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = &
                                                                    aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                    operator_value)
                                            pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = &
                                                                    -dble(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                    operator_value)
                                        else
                                            pair_index = unique_indnz2(itr_indnz2,1)
                                            pair_values(si_real+2*pair_index,:,:) = &
                                                                    pair_values(si_real+2*pair_index,:,:) + &
                                                                    dble(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                    operator_value)
                                            pair_values(si_real+2*pair_index+1,:,:) = &
                                                                    pair_values(si_real+2*pair_index+1,:,:) + &
                                                                    aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                    operator_value)
                                            pair_values(si_imag+2*pair_index,:,:) = &
                                                                    pair_values(si_imag+2*pair_index,:,:) + &
                                                                    aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                    operator_value)
                                            pair_values(si_imag+2*pair_index+1,:,:) = &
                                                                    pair_values(si_imag+2*pair_index+1,:,:) - & 
                                                                    dble(-ci*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                    operator_value)
                                        endif
                                    endif
                                enddo
                            elseif (conj_np1 == 0) then                                                     ! Run this part instead if hermiticity relation is not applied to connect ADOs 
                                                                                                            ! and we do not need to take the conjugate transpose 
                                do ndash = 0,dim_rho-1
                                    indnz2 = rho_sparsity(ndash,ncol,indjnp1)
                                    operator_value = d_ops(nrow,ndash,el_np1,sign_np1)
                                    operator_value_log = d_ops_log(nrow,ndash,el_np1,sign_np1)
                                    if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                                        already_indnz2 = .false.
                                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                                            indnz2_compare = unique_indnz2(itr_indnz2,0)
                                            if (indnz2 == indnz2_compare) then
                                                already_indnz2 = .true.
                                                exit
                                            endif
                                        enddo
                                        if (already_indnz2 .eqv. .false.) then
                                            count_pairs_total = count_pairs_total + 1
                                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                            unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                            unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                            pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                            pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                            pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                            pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                            pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                            pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                            pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                            pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                            pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = &
                                                                    dble(-ci*perm_np1*operator_value)
                                            pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = &
                                                                    -aimag(-ci*perm_np1*operator_value)
                                            pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = &
                                                                    aimag(-ci*perm_np1*operator_value)
                                            pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = &
                                                                    dble(-ci*perm_np1*operator_value)
                                        else
                                            pair_index = unique_indnz2(itr_indnz2,1)
                                            pair_values(si_real+2*pair_index,:,:) = &
                                                                    pair_values(si_real+2*pair_index,:,:) + &
                                                                    dble(-ci*perm_np1*operator_value)
                                            pair_values(si_real+2*pair_index+1,:,:) = &
                                                                    pair_values(si_real+2*pair_index+1,:,:) - &
                                                                    aimag(-ci*perm_np1*operator_value)
                                            pair_values(si_imag+2*pair_index,:,:) = &
                                                                    pair_values(si_imag+2*pair_index,:,:) + &
                                                                    aimag(-ci*perm_np1*operator_value)
                                            pair_values(si_imag+2*pair_index+1,:,:) = &
                                                                    pair_values(si_imag+2*pair_index+1,:,:) + & 
                                                                    dble(-ci*perm_np1*operator_value)
                                        endif
                                    endif

                                    indnz2 = rho_sparsity(nrow,ndash,indjnp1)
                                    operator_value = d_ops(ndash,ncol,el_np1,sign_np1)
                                    operator_value_log = d_ops_log(ndash,ncol,el_np1,sign_np1)
                                    if ((indnz2 .ne. -1) .and. (operator_value_log .eqv. .true.)) then
                                        already_indnz2 = .false.
                                        do itr_indnz2 = 0,n_indnz2_this_indnz1
                                            indnz2_compare = unique_indnz2(itr_indnz2,0)
                                            if (indnz2 == indnz2_compare) then
                                                already_indnz2 = .true.
                                                exit
                                            endif
                                        enddo
                                        if (already_indnz2 .eqv. .false.) then
                                            count_pairs_total = count_pairs_total + 1
                                            n_indnz2_this_indnz1 = n_indnz2_this_indnz1 + 1
                                            unique_indnz2(n_indnz2_this_indnz1,0) = indnz2
                                            unique_indnz2(n_indnz2_this_indnz1,1) = n_indnz2_this_indnz1
                                            pair_info_row(si_real+2*n_indnz2_this_indnz1) = 2*indnz1
                                            pair_info_row(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz1
                                            pair_info_row(si_imag+2*n_indnz2_this_indnz1) = 2*indnz1+1
                                            pair_info_row(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz1+1
                                            pair_info_col(si_real+2*n_indnz2_this_indnz1) = 2*indnz2
                                            pair_info_col(si_real+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                            pair_info_col(si_imag+2*n_indnz2_this_indnz1) = 2*indnz2
                                            pair_info_col(si_imag+2*n_indnz2_this_indnz1+1) = 2*indnz2+1
                                            pair_values(si_real+2*n_indnz2_this_indnz1,:,:) = &
                                                                    dble(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                                    operator_value)
                                            pair_values(si_real+2*n_indnz2_this_indnz1+1,:,:) = &
                                                                    -aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                                    operator_value)
                                            pair_values(si_imag+2*n_indnz2_this_indnz1,:,:) = &
                                                                    aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                                    operator_value)
                                            pair_values(si_imag+2*n_indnz2_this_indnz1+1,:,:) = &
                                                                    dble(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                                    operator_value)
                                        else
                                            pair_index = unique_indnz2(itr_indnz2,1)
                                            pair_values(si_real+2*pair_index,:,:) = &
                                                                    pair_values(si_real+2*pair_index,:,:) + &
                                                                    dble(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                                    operator_value)
                                            pair_values(si_real+2*pair_index+1,:,:) = &
                                                                    pair_values(si_real+2*pair_index+1,:,:) - &
                                                                    aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                                    operator_value)
                                            pair_values(si_imag+2*pair_index,:,:) = &
                                                                    pair_values(si_imag+2*pair_index,:,:) + &
                                                                    aimag(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                                    operator_value)
                                            pair_values(si_imag+2*pair_index+1,:,:) = &
                                                                    pair_values(si_imag+2*pair_index+1,:,:) + & 
                                                                    dble(-ci*((-1.d0)**(itrn+1))*perm_np1*&
                                                                    operator_value)
                                        endif
                                    endif
                                enddo
                            endif
                        enddo
                    endif
                endif
            enddo
        endif
    enddo

end subroutine sparse_matrix_elements_b_wbl
