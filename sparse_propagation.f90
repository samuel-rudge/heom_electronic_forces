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
!       rho_system                                  Array of size [nsteps,dim_rho,dim_rho] containing the system density matrix (zeroth tier ADO) at each timestep
!       
!       current                                     Array of size [nsteps,nleads] containing the current through each lead at each timestep

subroutine sparse_propagation(pair_info,pair_values,ham_pair_info,rho_nonzeros,gamma_values,&
                                ham,tier_index,rho_0,dt,nsteps,&
                                max_expan_order,dim_rho,nnz_elements,npairs,nhampairs,&
                                len_un_ind,nthreads)

    implicit none

    include 'omp_lib.h'                                                                             ! Import the OpenMP (OMP) library for parallel programming
    integer, intent(in) :: nsteps,dim_rho,max_expan_order,nnz_elements,npairs,nhampairs,len_un_ind
    integer, intent(in) :: nthreads
    complex*16, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: ham,rho_0                         ! Define input variables and arrays    
    complex*16, intent(in), dimension(0:nnz_elements-1) :: gamma_values
    complex*16, intent(in), dimension(0:npairs-1) :: pair_values
    integer, intent(in), dimension(0:npairs-1,0:3) :: pair_info
    integer, intent(in), dimension(0:nhampairs-1,0:4) :: ham_pair_info
    integer, intent(in), dimension(0:nnz_elements-1,0:2) :: rho_nonzeros
    integer, intent(in), dimension(0:len_un_ind-1) :: tier_index
    real*8, intent(in) :: dt

    complex*16, dimension(0:dim_rho-1,0:dim_rho-1) :: rho_system                ! Define output variables and arrays
    
    complex*16, dimension(0:nnz_elements-1) :: rho                                                    ! Define other necessary arrays
    complex*16, dimension(0:nnz_elements-1) :: rho_temp
    complex*16, dimension(0:nnz_elements-1) :: rho_deriv
    complex*16, parameter :: ci=(0.d0,1.d0)                                                         ! Define imaginary number (=sqrt(-1))   
    integer :: indjn,itrnz,nrow,ncol,itrt,itrl,itrpair,itrhpair,conj                               ! Define other necessary varaibles
    integer :: leftright,indnz1,indnz2,itrn,itr_rho
    real*8 :: time,norm
    character(len=32) :: output_fmt_rho

    write(output_fmt_rho,'(a,i0,a)') '(f15.5,',dim_rho,'(f20.15),f10.6)' 

    rho = 0.d0                                                                                      ! Create array containing initial values of all nonzero ADO elements. 
                                                                                ! rho is a temporary array that will be updated each timestep to contain all nonzero 
                                                                                ! ADO elements at that timstep
    do itrnz = 0,nnz_elements-1                                                                     ! Loop through nonzero elements (elements that must be propagated)
        indjn = rho_nonzeros(itrnz,0)                                                           ! For this nonzero element, find the ADO index and row/column value from rho_nonzero
        itrn = tier_index(indjn)
        if (itrn == 0) then                                                                         ! If it is a 0th tier ADO element, it is part of the initial system density matrix
            nrow = rho_nonzeros(itrnz,1)                                                            ! So find the row and column of this ADO element
            ncol = rho_nonzeros(itrnz,2)
            rho(itrnz) = rho_0(nrow,ncol)                                                           ! Fill the corresponding element of rho with the initial condition
        else
            exit
        endif
    enddo

    open(100, file = 'rho.txt')                                                                       ! Open 3 text files to write into information about the system density matrix, the current, 
    time = 0.d0
    rho_system = rho_0
    write(100,output_fmt_rho) time,(dreal(rho_system(itr_rho,itr_rho)),itr_rho=0,dim_rho-1)

    call omp_set_num_threads(nthreads)
    do itrt = 1,nsteps                                                                              ! Loop through timesteps
        print *,itrt
        time = dt*dble(itrt)                                                                        ! Calculate 'in-simulation' time 
        rho_temp = rho                                                                              ! At this point we need to perform the one-step-propagation algorithm, which is outlined in 
                                                                                    ! great detail in sparsity.f90. 
        do itrl = 1,max_expan_order                                                                 ! Loop through the Taylor series expansion of e^{L*dt}
            rho_deriv = 0.d0                                                                        ! rho_deriv initialized as 0
            !$OMP PARALLEL SHARED(rho_deriv,rho_temp)                                                           
                                                                                        ! Construct parallel region of code with variables rho_deriv and rho_temp to be shared between
                                                                                        ! parallel-running threads
            !$OMP DO PRIVATE(indnz1,indnz2,nrow,ncol,leftright,itrhpair)  
                                                                                        ! Let OMP know that the next line will start a for (do) loop that is to be run in parallel,
                                                                                        ! with each thread having its own indnz1,indnz2,nrow,ncol,lefright,itrhpair
            do itrhpair = 0,nhampairs-1         
                indnz1 = ham_pair_info(itrhpair,0)                                                  ! Return nonzero index of LHS ADO element
                indnz2 = ham_pair_info(itrhpair,1)                                                  ! Return nonzero index of RHS ADO element
                nrow = ham_pair_info(itrhpair,2)                                                    ! Return row and column of Hamiltonian coupling these two ADO elements    
                ncol = ham_pair_info(itrhpair,3)
                leftright = ham_pair_info(itrhpair,4)                                               ! Return 0 if Hamiltonian acts on left and 1 if Hamiltonian acts on right of ADO

                rho_deriv(indnz1) = rho_deriv(indnz1) - ci*((-1)**(leftright))*ham(nrow,ncol)*rho_temp(indnz2)
                                                                                            ! Add contribution to rho_deriv. (-1)**(leftright) accounts for the different sign 
                                                                                            ! of the two commutator terms
            enddo
            !$OMP END DO

            !$OMP DO PRIVATE(itrnz)
            do itrnz = 0,nnz_elements-1                                                             ! This parallelized for loop runs through all nonzero ADO elements and incorporates 
                                                                                        ! their contribution from the gamma term in the HQME
                rho_deriv(itrnz) = rho_deriv(itrnz) + gamma_values(itrnz)*rho_temp(itrnz)
            enddo
            !$OMP END DO

            !$OMP DO PRIVATE(itrpair,indnz1,indnz2,conj)
            do itrpair = 0,npairs-1                                                                 ! This parallelized for loop runs through all other connections in the HQME
                indnz1 = pair_info(itrpair,0)                                                       ! Generate the nonzero index of the ADO element from the LHS of the HQME
                indnz2 = pair_info(itrpair,1)                                                       ! Generate the nonzero index of the ADO element from the RHS of the HQME
                conj = pair_info(itrpair,2)                                                         ! Determine whether a hermiticity relation was applied to ADO element with nonzero index
                                                                                                ! indnz2 
                if (conj == 1) then                                                                 ! If a hermiticity relation was applied, we need to take the conjugate of this ADO element
                    rho_deriv(indnz1) = rho_deriv(indnz1) + pair_values(itrpair)*conjg(rho_temp(indnz2))
                                                                                            ! Update the ADO element on the LHS with this contribution
                elseif (conj == 0) then                                                             ! Run this in the case of no hermiticity relation, so no conjugation necessary
                    rho_deriv(indnz1) = rho_deriv(indnz1) + pair_values(itrpair)*rho_temp(indnz2)
                endif
            enddo
            !$OMP END DO      
            !$OMP END PARALLEL  

            rho_deriv = rho_deriv*dt/dble(itrl)                                                     ! Multiply by timestep and divide by l
            rho = rho + rho_deriv                                                                   ! Update rho at this timestep
            rho_temp = rho_deriv                                                                    ! Set rho_temp = rho_deriv to prepare next loop in Taylor series expansion

        enddo                                                                                       ! At the end of this process, rho contains all nonzero ADO elements for this timestep

        norm = 0.d0
        do itrnz = 0,nnz_elements-1                                                                 ! Now we need to identify the elements of the system density matrix at this timestep, so loop
                                                                                                    ! through all nonzero elements
            indjn = rho_nonzeros(itrnz,0)                                                           ! For this nonzero element, find the ADO index and row/column value from rho_nonzero
            nrow = rho_nonzeros(itrnz,1)
            ncol = rho_nonzeros(itrnz,2)
            itrn = tier_index(indjn)                                                                ! Set itrn = tier of the current ADO
            if (itrn == 0) then
                rho_system(nrow,ncol) = rho(itrnz)                                             ! If it corresponds to the 0th tier, then put this element in the correct place in the system
                if (nrow == ncol) then
                    norm = norm + dreal(rho(itrnz))
                endif
            else
                exit
            endif
        enddo
        
        write(100,output_fmt_rho) time,(dreal(rho_system(itr_rho,itr_rho)),itr_rho=0,dim_rho-1)

        if (abs(norm - 1.d0) > 0.01) then
            exit
        endif
    enddo

    close(100)                                                                                        ! Close text files after all timesteps are calculated

end subroutine

! ---------------------------------------------------------------------------
! 
!                     PERFORM TIME-PROPAGATION OF ADOs
!                         IN THE WIDE-BAND LIMIT
!
! ---------------------------------------------------------------------------
!
! Same as sparse_propagation, but now under the WBL.

subroutine sparse_propagation_wbl(pair_info,pair_values,ham_pair_info,wbl_pair_info,wbl_pair_values,&
                                rho_nonzeros,gamma_values,ham,tier_index,&
                                el_lead_couplings,rho_0,dt,nsteps,max_expan_order,nleads,dim_rho,&
                                nnz_elements,npairs,nhampairs,nwblpairs,len_un_ind,nel,nthreads)
                                
    implicit none

    include 'omp_lib.h'                                                                             ! Import the OpenMP (OMP) library for parallel programming
    integer, intent(in) :: nsteps,dim_rho,max_expan_order,nnz_elements,npairs,nhampairs,nwblpairs
    integer, intent(in) :: len_un_ind,nel,nleads,nthreads
    complex*16, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: ham,rho_0                         ! Define input variables and arrays    
    complex*16, intent(in), dimension(0:nnz_elements-1) :: gamma_values
    complex*16, intent(in), dimension(0:npairs-1) :: pair_values
    real*8, intent(in), dimension(0:nwblpairs-1) :: wbl_pair_values
    integer, intent(in), dimension(0:npairs-1,0:3) :: pair_info
    integer, intent(in), dimension(0:nwblpairs-1,0:3) :: wbl_pair_info
    integer, intent(in), dimension(0:nhampairs-1,0:4) :: ham_pair_info
    integer, intent(in), dimension(0:nnz_elements-1,0:2) :: rho_nonzeros
    integer, intent(in), dimension(0:len_un_ind-1) :: tier_index
    real*8, intent(in), dimension(0:nleads-1,0:nel-1) :: el_lead_couplings
    real*8, intent(in) :: dt

    complex*16, dimension(0:dim_rho-1,0:dim_rho-1) :: rho_system                ! Define output variables and arrays

    complex*16, dimension(0:nnz_elements-1) :: rho                                                    ! Define other necessary arrays
    complex*16, dimension(0:nnz_elements-1) :: rho_temp
    complex*16, dimension(0:nnz_elements-1) :: rho_deriv

    complex*16, parameter :: ci=(0.d0,1.d0)                                                         ! Define imaginary number (=sqrt(-1))   
    integer :: indjn,itrnz,nrow,ncol,itrt,itrl,itrpair,itrhpair,conj                               ! Define other necessary varaibles
    integer :: leftright,indnz1,indnz2,itrn
    integer :: itrwblpair,itrel,itrlead,itr_rho
    real*8 :: time,wbl_value,norm        
    character(len=32) :: output_fmt_rho

    write(output_fmt_rho,'(a,i0,a)') '(f15.5,',dim_rho,'(f20.15),f10.6)' 

    rho = 0.d0                                                                                      ! Create array containing initial values of all nonzero ADO elements. 
                                                                                ! rho is a temporary array that will be updated each timestep to contain all nonzero 
                                                                                ! ADO elements at that timstep
    do itrnz = 0,nnz_elements-1                                                                     ! Loop through nonzero elements (elements that must be propagated)
        indjn = rho_nonzeros(itrnz,0)                                                           ! For this nonzero element, find the ADO index and row/column value from rho_nonzero
        itrn = tier_index(indjn)
        if (itrn == 0) then                                                                         ! If it is a 0th tier ADO element, it is part of the initial system density matrix
            nrow = rho_nonzeros(itrnz,1)                                                            ! So find the row and column of this ADO element
            ncol = rho_nonzeros(itrnz,2)
            rho(itrnz) = rho_0(nrow,ncol)                                                           ! Fill the corresponding element of rho with the initial condition
        else                                                                      ! Since all elements are ordered, once we get to the 1st tier ADOs there are no more nonzero
                                                                                    ! elements from the zeroth tier, so can exit the for loop
            exit
        endif
    enddo

    open(100, file = 'rho.txt')                                                                       ! Open 3 text files to write into information about the system density matrix, the current, 
    time = 0.d0
    rho_system = rho_0
    write(100,output_fmt_rho) time,(dreal(rho_system(itr_rho,itr_rho)),itr_rho=0,dim_rho-1)

    call omp_set_num_threads(nthreads)
    do itrt = 1,nsteps                                                                              ! Loop through timesteps
        print *,itrt
        time = dt*dble(itrt)                                                                        ! Calculate 'in-simulation' time 
        rho_temp = rho                                                                              ! At this point we need to perform the one-step-propagation algorithm, which is outlined in 
                                                                                    ! great detail in sparsity.f90. 
        do itrl = 1,max_expan_order                                                                 ! Loop through the Taylor series expansion of e^{L*dt}
            rho_deriv = 0.d0                                                                        ! rho_deriv initialized as 0
            !$OMP PARALLEL SHARED(rho_deriv,rho_temp)                                                           
                                                                                        ! Construct parallel region of code with variables rho_deriv and rho_temp to be shared between
                                                                                        ! parallel-running threads
            !$OMP DO PRIVATE(indnz1,indnz2,nrow,ncol,leftright,itrhpair)  
                                                                                        ! Let OMP know that the next line will start a for (do) loop that is to be run in parallel,
                                                                                        ! with each thread having its own indnz1,indnz2,nrow,ncol,lefright,itrhpair
            do itrhpair = 0,nhampairs-1                                                 
                indnz1 = ham_pair_info(itrhpair,0)                                                  ! Return nonzero index of LHS ADO element
                indnz2 = ham_pair_info(itrhpair,1)                                                  ! Return nonzero index of RHS ADO element
                nrow = ham_pair_info(itrhpair,2)                                                    ! Return row and column of Hamiltonian coupling these two ADO elements    
                ncol = ham_pair_info(itrhpair,3)
                leftright = ham_pair_info(itrhpair,4)                                               ! Return 0 if Hamiltonian acts on left and 1 if Hamiltonian acts on right of ADO

                rho_deriv(indnz1) = rho_deriv(indnz1) - ci*((-1)**(leftright))*ham(nrow,ncol)*rho_temp(indnz2)
                                                                                            ! Add contribution to rho_deriv. (-1)**(leftright) accounts for the different sign 
                                                                                            ! of the two commutator terms
            enddo
            !$OMP END DO

            !$OMP DO PRIVATE(itrnz,itrel,itrlead)
            do itrnz = 0,nnz_elements-1                                                             ! This parallelized for loop runs through all nonzero ADO elements and incorporates 
                                                                                        ! their contribution from the gamma term in the HQME
                rho_deriv(itrnz) = rho_deriv(itrnz) + gamma_values(itrnz)*rho_temp(itrnz)
                do itrlead = 0,nleads-1
                    do itrel = 0,nel-1
                        rho_deriv(itrnz) = rho_deriv(itrnz) - ((el_lead_couplings(itrlead,itrel)**2)/2.d0)*&
                                            rho_temp(itrnz)
                    enddo
                enddo
            enddo
            !$OMP END DO

            !$OMP DO PRIVATE(itrpair,indnz1,indnz2,conj)
            do itrpair = 0,npairs-1                                                                 ! This parallelized for loop runs through all other connections in the HQME
                indnz1 = pair_info(itrpair,0)                                                       ! Generate the nonzero index of the ADO element from the LHS of the HQME
                indnz2 = pair_info(itrpair,1)                                                       ! Generate the nonzero index of the ADO element from the RHS of the HQME
                conj = pair_info(itrpair,2)                                                         ! Determine whether a hermiticity relation was applied to ADO element with nonzero index
                                                                                                ! indnz2 
                if (conj == 1) then                                                                 ! If a hermiticity relation was applied, we need to take the conjugate of this ADO element
                    rho_deriv(indnz1) = rho_deriv(indnz1) + pair_values(itrpair)*conjg(rho_temp(indnz2))
                                                                                            ! Update the ADO element on the LHS with this contribution
                elseif (conj == 0) then                                                             ! Run this in the case of no hermiticity relation, so no conjugation necessary
                    rho_deriv(indnz1) = rho_deriv(indnz1) + pair_values(itrpair)*rho_temp(indnz2)
                endif
            enddo
            !$OMP END DO 

            !$OMP DO PRIVATE(indnz1,indnz2,itrel,itrwblpair,itrlead,wbl_value)
            do itrwblpair = 0,nwblpairs-1
                indnz1 = wbl_pair_info(itrwblpair,0)
                indnz2 = wbl_pair_info(itrwblpair,1)
                itrel = wbl_pair_info(itrwblpair,2)
                wbl_value = wbl_pair_values(itrwblpair)
                do itrlead = 0,nleads-1
                    rho_deriv(indnz1) = rho_deriv(indnz1) - ((el_lead_couplings(itrlead,itrel)**2)/2.d0)&
                                        *wbl_value*rho_temp(indnz2)
                enddo
            enddo
            !$OMP END DO     

            !$OMP END PARALLEL

            rho_deriv = rho_deriv*dt/dble(itrl)                                                     ! Multiply by timestep and divide by l
            rho = rho + rho_deriv                                                                   ! Update rho at this timestep
            rho_temp = rho_deriv                                                                    ! Set rho_temp = rho_deriv to prepare next loop in Taylor series expansion

        enddo                                                                                       ! At the end of this process, rho contains all nonzero ADO elements for this timestep
                                                                                     ! At the end of this process, rho contains all nonzero ADO elements for this timestep
        norm = 0.d0
        do itrnz = 0,nnz_elements-1                                                                 ! Now we need to identify the elements of the system density matrix at this timestep, so loop
                                                                                                    ! through all nonzero elements
            indjn = rho_nonzeros(itrnz,0)                                                           ! For this nonzero element, find the ADO index and row/column value from rho_nonzero
            nrow = rho_nonzeros(itrnz,1)
            ncol = rho_nonzeros(itrnz,2)
            itrn = tier_index(indjn)                                                                ! Set itrn = tier of the current ADO
            if (itrn == 0) then
                rho_system(nrow,ncol) = rho(itrnz)                                             ! If it corresponds to the 0th tier, then put this element in the correct place in the system
                if (nrow == ncol) then
                    norm = norm + dreal(rho(itrnz))
                endif
            else
                exit
            endif
        enddo
        
        write(100,output_fmt_rho) time,(dreal(rho_system(itr_rho,itr_rho)),itr_rho=0,dim_rho-1)

        if (abs(norm - 1.d0) > 0.01) then
            exit
        endif
    enddo

    close(100)                                                                                        ! Close text files after all timesteps are calculated
    
end subroutine