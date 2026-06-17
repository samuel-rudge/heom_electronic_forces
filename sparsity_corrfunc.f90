

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

    integer, intent(in), dimension(0:nnz_elements-1,0:2) :: rho_nonzeros                           ! Define output variables
    integer, intent(out) :: npairs,n_indnz2_this_indnz1_max
    integer, intent(out), dimension(0:nnz_elements):: n_indnz2_prev_indnz1_vec

    integer, dimension(0:nnz_elements-1) :: unique_indnz2
    integer :: nrow,ncol,itrn,indjn,ndash,indnz2,indnz1                      ! Define other necessary variables
    integer :: itrjnm1,jnm1,sign_nm1,el_nm1,indjnm1,conj_nm1,conj_np1
    integer :: jnp1,sign_np1,el_np1,eldash_np1,indjnp1,count_pairs_total
    integer :: n_indnz2_this_indnz1,itr_indnz2,indnz2_compare
    logical :: ham_value,operator_value,already_indnz2

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

subroutine sparse_matrix_elements_b(ksiglm,tier_index,index_minus,index_plus,d_ops,eta_vec,&
                                    rho_sparsity,rho_nonzeros,pair_info_row,pair_info_col,pair_values,&
                                    nnz_elements,dim_rho,len_index_plus,len_un_ind,len_index_minus,&
                                    nmax,nmodes,nel,nsign,npairs,nleads,npoles,ham_log,d_ops_log,&
                                    degenerate_levels,n_indnz2_this_indnz1_max,deriv_ham,len_x_vec,&
                                    n_indnz2_prev_indnz1_vec,dv_dx,average_electronic_force_mol_vec,&
                                    average_electronic_force_molleads_vec)

    implicit none

    integer, intent(in) :: len_index_plus,len_index_minus,nsign,len_un_ind,nnz_elements,npairs      ! Define input variables and arrays
    integer, intent(in) :: nmax,nmodes,nel,dim_rho,npoles,nleads,n_indnz2_this_indnz1_max,len_x_vec
    logical, intent(in) :: degenerate_levels
    real*8, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops
    complex*16, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:len_x_vec-1) :: deriv_ham
    double precision, intent(in), dimension(0:len_x_vec-1) :: average_electronic_force_mol_vec
    double precision, intent(in), dimension(0:len_x_vec-1) :: average_electronic_force_molleads_vec
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:nel-1,0:1) :: d_ops_log
    logical, intent(in), dimension(0:dim_rho-1,0:dim_rho-1) :: ham_log
    integer, intent(in), dimension(0:nmodes-1,0:3) :: ksiglm
    integer, intent(in), dimension(0:len_un_ind-1) :: tier_index
    integer, intent(in), dimension(0:len_index_plus-1,0:nmodes-1,0:2) :: index_plus
    integer, intent(in), dimension(0:len_index_minus-1,0:nmax-1,0:3) :: index_minus
    complex*16, intent(in),dimension(0:nleads-1,0:nsign-1,0:npoles) :: eta_vec
    integer, intent(in), dimension(0:dim_rho-1,0:dim_rho-1,0:len_un_ind-1) :: rho_sparsity
    integer, intent(in), dimension(0:nnz_elements-1,0:2) :: rho_nonzeros
    integer, intent(in), dimension(0:nnz_elements) :: n_indnz2_prev_indnz1_vec
    double precision, intent(in), dimension(0:nleads-1,0:nel-1,0:len_x_vec-1) :: dv_dx

    integer, intent(out), dimension(0:4*npairs-1) :: pair_info_row,pair_info_col
    real*8, intent(out), dimension(0:4*npairs-1,0:len_x_vec-1) :: pair_values

    integer :: nrow,ncol,itrn,indjn,ndash,indnz2,indnz1,leads_np1,poles_np1
    integer :: itrjnm1,jnm1,sign_nm1,el_nm1,indjnm1,leads_nm1,poles_nm1,conj_np1
    integer :: jnp1,sign_np1,el_np1,eldash_np1,indjnp1,count_pairs_total,conj_nm1
    integer :: n_indnz2_this_indnz1,itr_indnz2,indnz2_compare,pair_index,itrx
    integer :: si_real,si_imag
    integer, dimension(0:n_indnz2_this_indnz1_max,0:1) :: unique_indnz2
    complex*16 :: ham_value,average_electronic_force_value
    real*8 :: perm_np1,perm_nm1,conj_button_nm1,conj_button_np1,operator_value
    logical :: ham_value_log,operator_value_log,already_indnz2
    complex*16, parameter :: ci=(0.d0,1.d0)                                                      ! Define imaginary number (=sqrt(-1))

    count_pairs_total = -1                                                                       ! Start the pair count at -1
    pair_info_row = -1                                                                           ! Set the indices columns in pair_info = -1 to start
    pair_info_col = -1                                                                           ! Set the indices columns in pair_info = -1 to start
    pair_values = 0.d0                       
                                                        ! Set pair and gamma values initially to 0 (double precision)

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
                        ham_value = deriv_ham(nrow,ndash,itrx)
                        pair_values(si_real+2*n_indnz2_this_indnz1,itrx) = dble(0.5*ham_value)
                        pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx) = -aimag(0.5*ham_value)
                        pair_values(si_imag+2*n_indnz2_this_indnz1,itrx) = aimag(0.5*ham_value)
                        pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx) =  dble(0.5*ham_value)
                    enddo
                else
                    pair_index = unique_indnz2(itr_indnz2,1)
                    do itrx = 0,(len_x_vec-1)
                        ham_value = deriv_ham(nrow,ndash,itrx)
                        pair_values(si_real+2*pair_index,itrx) = pair_values(si_real+2*pair_index,itrx) + dble(0.5*ham_value)
                        pair_values(si_real+2*pair_index+1,itrx) = pair_values(si_real+2*pair_index+1,itrx) - aimag(0.5*ham_value)
                        pair_values(si_imag+2*pair_index,itrx) = pair_values(si_imag+2*pair_index,itrx) + aimag(0.5*ham_value)
                        pair_values(si_imag+2*pair_index+1,itrx) = pair_values(si_imag+2*pair_index+1,itrx) + dble(0.5*ham_value)
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
                        ham_value = deriv_ham(ndash,ncol,itrx)
                        pair_values(si_real+2*n_indnz2_this_indnz1,itrx) = dble(0.5*ham_value)
                        pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx) = -aimag(0.5*ham_value)
                        pair_values(si_imag+2*n_indnz2_this_indnz1,itrx) = aimag(0.5*ham_value)
                        pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx) = dble(0.5*ham_value)
                    enddo
                else
                    pair_index = unique_indnz2(itr_indnz2,1)
                    do itrx = 0,(len_x_vec-1)
                        ham_value = deriv_ham(ndash,ncol,itrx)
                        pair_values(si_real+2*pair_index,itrx) = &
                                    pair_values(si_real+2*pair_index,itrx) + dble(0.5*ham_value)
                        pair_values(si_real+2*pair_index+1,itrx) = &
                                    pair_values(si_real+2*pair_index+1,itrx) - aimag(0.5*ham_value)
                        pair_values(si_imag+2*pair_index,itrx) = &
                                    pair_values(si_imag+2*pair_index,itrx) + aimag(0.5*ham_value)
                        pair_values(si_imag+2*pair_index+1,itrx) = &
                                    pair_values(si_imag+2*pair_index+1,itrx) + dble(0.5*ham_value)
                    enddo
                endif
            endif
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
            do itrx = 0,len_x_vec-1
                average_electronic_force_value = (average_electronic_force_mol_vec(itrx) + &
                                                    average_electronic_force_molleads_vec(itrx))
                pair_values(si_real+2*n_indnz2_this_indnz1,itrx) = dble(average_electronic_force_value)
                pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx) = -aimag(average_electronic_force_value)
                pair_values(si_imag+2*n_indnz2_this_indnz1,itrx) = aimag(average_electronic_force_value)
                pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx) = dble(average_electronic_force_value)
            enddo
        else
            do itrx = 0,len_x_vec-1
                pair_index = unique_indnz2(itr_indnz2,1)
                average_electronic_force_value = (average_electronic_force_mol_vec(itrx) + &
                                                    average_electronic_force_molleads_vec(itrx))
                pair_values(si_real+2*pair_index,itrx) = &
                    pair_values(si_real+2*pair_index,itrx) + dble(average_electronic_force_value)
                pair_values(si_real+2*pair_index+1,itrx) = &
                    pair_values(si_real+2*pair_index+1,itrx) - aimag(average_electronic_force_value)
                pair_values(si_imag+2*pair_index,itrx) = &
                    pair_values(si_imag+2*pair_index,itrx) + aimag(average_electronic_force_value)
                pair_values(si_imag+2*pair_index+1,itrx) = &
                    pair_values(si_imag+2*pair_index+1,itrx) + dble(average_electronic_force_value)
            enddo
        endif

        if (itrn > 0) then
                                                                                                    ! Each nonzero ADO element directly couples to itself via the gamm term in the HEOM
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
                                    pair_values(si_real+2*n_indnz2_this_indnz1,itrx) = &
                                            dble(0.5*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx) = &
                                            aimag(0.5*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_imag+2*n_indnz2_this_indnz1,itrx) = &
                                            aimag(0.5*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx) = &
                                            -dble(0.5*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                enddo
                            else
                                do itrx = 0,len_x_vec-1
                                    pair_index = unique_indnz2(itr_indnz2,1)
                                    pair_values(si_real+2*pair_index,itrx) = &
                                                pair_values(si_real+2*pair_index,itrx) + &
                                                dble(0.5*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_real+2*pair_index+1,itrx) = &
                                                pair_values(si_real+2*pair_index+1,itrx) + &
                                                aimag(0.5*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_imag+2*pair_index,itrx) = &
                                                pair_values(si_imag+2*pair_index,itrx) + &
                                                aimag(0.5*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_imag+2*pair_index+1,itrx) = &
                                                pair_values(si_imag+2*pair_index+1,itrx) - & 
                                                dble(0.5*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                conj_button_nm1*eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                                dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
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
                                    pair_values(si_real+2*n_indnz2_this_indnz1,itrx) = &
                                                dble(0.5*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                conj_button_nm1*conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx) = &
                                                aimag(0.5*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                conj_button_nm1*conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_imag+2*n_indnz2_this_indnz1,itrx) = &
                                                aimag(0.5*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                conj_button_nm1*conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx) = &
                                                -dble(0.5*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                conj_button_nm1*conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                enddo
                            else
                                pair_index = unique_indnz2(itr_indnz2,1)
                                do itrx = 0,len_x_vec-1   
                                    pair_values(si_real+2*pair_index,itrx) = pair_values(si_real+2*pair_index,itrx) + &
                                                dble(0.5*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                conj_button_nm1*&
                                                conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_real+2*pair_index+1,itrx) = pair_values(si_real+2*pair_index+1,itrx) + &
                                                aimag(0.5*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                conj_button_nm1*&
                                                conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_imag+2*pair_index,itrx) = pair_values(si_imag+2*pair_index,itrx) + &
                                                aimag(0.5*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                conj_button_nm1*&
                                                conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_imag+2*pair_index+1,itrx) = pair_values(si_imag+2*pair_index+1,itrx) - & 
                                                dble(0.5*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                                conj_button_nm1*&
                                                conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                                dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
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
                                    pair_values(si_real+2*n_indnz2_this_indnz1,itrx) = &
                                            dble(0.5*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)                                                
                                    pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx) = &
                                            -aimag(0.5*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_imag+2*n_indnz2_this_indnz1,itrx) = &
                                            aimag(0.5*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx) = &
                                            dble(0.5*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                enddo
                            else
                                pair_index = unique_indnz2(itr_indnz2,1)
                                do itrx = 0,len_x_vec-1   
                                    pair_values(si_real+2*pair_index,itrx) = &
                                            pair_values(si_real+2*pair_index,itrx) + &
                                            dble(0.5*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_real+2*pair_index+1,itrx) = &
                                            pair_values(si_real+2*pair_index+1,itrx) - &
                                            aimag(0.5*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_imag+2*pair_index,itrx) = &
                                            pair_values(si_imag+2*pair_index,itrx) + &
                                            aimag(0.5*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_imag+2*pair_index+1,itrx) = &
                                            pair_values(si_imag+2*pair_index+1,itrx) + & 
                                            dble(0.5*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            eta_vec(leads_nm1,sign_nm1,poles_nm1)*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
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
                                    pair_values(si_real+2*n_indnz2_this_indnz1,itrx) = &
                                            dble(0.5*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx) = &
                                            -aimag(0.5*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_imag+2*n_indnz2_this_indnz1,itrx) = &
                                            aimag(0.5*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx) = &
                                            dble(0.5*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                enddo
                            else
                                do itrx = 0,len_x_vec-1   
                                    pair_index = unique_indnz2(itr_indnz2,1)
                                    pair_values(si_real+2*pair_index,itrx) = pair_values(si_real+2*pair_index,itrx) + &
                                            dble(0.5*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_real+2*pair_index+1,itrx) = pair_values(si_real+2*pair_index+1,itrx) - &
                                            aimag(0.5*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_imag+2*pair_index,itrx) = pair_values(si_imag+2*pair_index,itrx) + &
                                            aimag(0.5*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
                                    pair_values(si_imag+2*pair_index+1,itrx) = pair_values(si_imag+2*pair_index+1,itrx) + & 
                                            dble(0.5*((-1.d0)**(itrn-1))*((-1.d0)**(itrn-itrjnm1-1))*perm_nm1*&
                                            conjg(eta_vec(leads_nm1,1-sign_nm1,poles_nm1))*&
                                            dv_dx(leads_nm1,el_nm1,itrx)*operator_value)
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
                                            pair_values(si_real+2*n_indnz2_this_indnz1,itrx) = &
                                                                dble(0.5*perm_np1*conj_button_np1*&
                                                                dv_dx(leads_np1,el_np1,itrx)*operator_value)
                                            ! pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx) = &
                                            !                     aimag(0.5*perm_np1*conj_button_np1*&
                                            !                     dv_dx(leads_np1,el_np1,itrx)*operator_value)
                                            ! pair_values(si_imag+2*n_indnz2_this_indnz1,itrx) = &
                                            !                     aimag(0.5*perm_np1*conj_button_np1*&
                                            !                     dv_dx(leads_np1,el_np1,itrx)*operator_value)
                                            pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx) = &
                                                                -dble(0.5*perm_np1*conj_button_np1*&
                                                                dv_dx(leads_np1,el_np1,itrx)*operator_value)
                                        enddo
                                    else
                                        do itrx = 0,len_x_vec-1 
                                            pair_index = unique_indnz2(itr_indnz2,1)
                                            pair_values(si_real+2*pair_index,itrx) = &
                                                                    pair_values(si_real+2*pair_index,itrx) + &
                                                                    dble(0.5*perm_np1*conj_button_np1*&
                                                                    dv_dx(leads_np1,el_np1,itrx)*operator_value)
                                            ! pair_values(si_real+2*pair_index+1,itrx) = &
                                            !                         pair_values(si_real+2*pair_index+1,itrx) + &
                                            !                         aimag(0.5*perm_np1*conj_button_np1*&
                                            !                         dv_dx(leads_np1,el_np1,itrx)*operator_value)
                                            ! pair_values(si_imag+2*pair_index,itrx) = &
                                            !                         pair_values(si_imag+2*pair_index,itrx) + &
                                            !                         aimag(0.5*perm_np1*conj_button_np1*&
                                            !                         dv_dx(leads_np1,el_np1,itrx)*operator_value)
                                            pair_values(si_imag+2*pair_index+1,itrx) = &
                                                                    pair_values(si_imag+2*pair_index+1,itrx) - & 
                                                                    dble(0.5*perm_np1*conj_button_np1*&
                                                                    dv_dx(leads_np1,el_np1,itrx)*operator_value)
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
                                            pair_values(si_real+2*n_indnz2_this_indnz1,itrx) = &
                                                                    dble(-0.5*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                    dv_dx(leads_np1,el_np1,itrx)*operator_value)
                                            ! pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx) = &
                                            !                         aimag(0.5*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                            !                         dv_dx(leads_np1,el_np1,itrx)*operator_value)
                                            ! pair_values(si_imag+2*n_indnz2_this_indnz1,itrx) = &
                                            !                         aimag(0.5*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                            !                         dv_dx(leads_np1,el_np1,itrx)*operator_value)
                                            pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx) = &
                                                                    -dble(-0.5*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                    dv_dx(leads_np1,el_np1,itrx)*operator_value)
                                        enddo
                                    else
                                        pair_index = unique_indnz2(itr_indnz2,1)
                                        do itrx = 0,len_x_vec-1   
                                            pair_values(si_real+2*pair_index,itrx) = &
                                                                    pair_values(si_real+2*pair_index,itrx) + &
                                                                    dble(-0.5*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                    dv_dx(leads_np1,el_np1,itrx)*operator_value)
                                            ! pair_values(si_real+2*pair_index+1,itrx) = &
                                            !                         pair_values(si_real+2*pair_index+1,itrx) + &
                                            !                         aimag(0.5*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                            !                         dv_dx(leads_np1,el_np1,itrx)*operator_value)
                                            ! pair_values(si_imag+2*pair_index,itrx) = &
                                            !                         pair_values(si_imag+2*pair_index,itrx) + &
                                            !                         aimag(0.5*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                            !                         dv_dx(leads_np1,el_np1,itrx)*operator_value)
                                            pair_values(si_imag+2*pair_index+1,itrx) = &
                                                                    pair_values(si_imag+2*pair_index+1,itrx) - & 
                                                                    dble(-0.5*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                    dv_dx(leads_np1,el_np1,itrx)*operator_value)
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
                                            pair_values(si_real+2*n_indnz2_this_indnz1,itrx) = &
                                                                dble(0.5*perm_np1*operator_value*dv_dx(leads_np1,el_np1,itrx))
                                            ! pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx) = &
                                            !                         -aimag(0.5*perm_np1*operator_value*dv_dx(leads_np1,el_np1,itrx))
                                            ! pair_values(si_imag+2*n_indnz2_this_indnz1,itrx) = &
                                            !                         aimag(0.5*perm_np1*operator_value*dv_dx(leads_np1,el_np1,itrx))
                                            pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx) = &
                                                                    dble(0.5*perm_np1*operator_value*dv_dx(leads_np1,el_np1,itrx))
                                        enddo
                                    else
                                        do itrx = 0,len_x_vec-1   
                                            pair_index = unique_indnz2(itr_indnz2,1)
                                            pair_values(si_real+2*pair_index,itrx) = &
                                                                    pair_values(si_real+2*pair_index,itrx) + &
                                                                    dble(0.5*perm_np1*operator_value*dv_dx(leads_np1,el_np1,itrx))
                                            ! pair_values(si_real+2*pair_index+1,itrx) = &
                                            !                         pair_values(si_real+2*pair_index+1,itrx) - &
                                            !                         aimag(0.5*perm_np1*operator_value*dv_dx(leads_np1,el_np1,itrx))
                                            ! pair_values(si_imag+2*pair_index,itrx) = &
                                            !                         pair_values(si_imag+2*pair_index,itrx) + &
                                            !                         aimag(0.5*perm_np1*operator_value*dv_dx(leads_np1,el_np1,itrx))
                                            pair_values(si_imag+2*pair_index+1,itrx) = &
                                                                    pair_values(si_imag+2*pair_index+1,itrx) + & 
                                                                    dble(0.5*perm_np1*operator_value*dv_dx(leads_np1,el_np1,itrx))
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
                                            pair_values(si_real+2*n_indnz2_this_indnz1,itrx) = &
                                                                dble(-0.5*((-1.d0)**(itrn+1))*perm_np1*&
                                                                operator_value*dv_dx(leads_np1,el_np1,itrx))
                                            ! pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx) = &
                                            !                         -aimag(0.5*((-1.d0)**(itrn+1))*perm_np1*&
                                            !                         operator_value*dv_dx(leads_np1,el_np1,itrx))
                                            ! pair_values(si_imag+2*n_indnz2_this_indnz1,itrx) = &
                                            !                         aimag(0.5*((-1.d0)**(itrn+1))*perm_np1*&
                                            !                         operator_value*dv_dx(leads_np1,el_np1,itrx))
                                            pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx) = &
                                                                    dble(-0.5*((-1.d0)**(itrn+1))*perm_np1*&
                                                                    operator_value*dv_dx(leads_np1,el_np1,itrx))
                                        enddo
                                    else
                                        pair_index = unique_indnz2(itr_indnz2,1)
                                        do itrx = 0,len_x_vec-1           
                                            pair_values(si_real+2*pair_index,itrx) = &
                                                                    pair_values(si_real+2*pair_index,itrx) + &
                                                                    dble(-0.5*((-1.d0)**(itrn+1))*perm_np1*&
                                                                    operator_value*dv_dx(leads_np1,el_np1,itrx))
                                            ! pair_values(si_real+2*pair_index+1,itrx) = &
                                            !                         pair_values(si_real+2*pair_index+1,itrx) - &
                                            !                         aimag(0.5*((-1.d0)**(itrn+1))*perm_np1*&
                                            !                         operator_value*dv_dx(leads_np1,el_np1,itrx))
                                            ! pair_values(si_imag+2*pair_index,itrx) = &
                                            !                         pair_values(si_imag+2*pair_index,itrx) + &
                                            !                         aimag(0.5*((-1.d0)**(itrn+1))*perm_np1*&
                                            !                         operator_value*dv_dx(leads_np1,el_np1,itrx))
                                            pair_values(si_imag+2*pair_index+1,itrx) = &
                                                                    pair_values(si_imag+2*pair_index+1,itrx) + & 
                                                                    dble(-0.5*((-1.d0)**(itrn+1))*perm_np1*&
                                                                    operator_value*dv_dx(leads_np1,el_np1,itrx))
                                        enddo
                                    endif
                                endif
                            enddo
                        endif
                    else
                        do el_np1 = 0,nel-1
                            if (conj_np1 == 1) then                                                         ! Run this code if a hermiticity relation was applied
                                conj_button_np1 = (-1.d0)**(floor((dble(itrn)+1.d0)/2.d0))                  ! Calculate the hermiticity prefactor generated when connecting these two ADOs
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
                                                pair_values(si_real+2*n_indnz2_this_indnz1,itrx) = &
                                                                        dble(0.5*perm_np1*conj_button_np1*&
                                                                        operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                ! pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx) = &
                                                !                         aimag(0.5*perm_np1*conj_button_np1*&
                                                !                         operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                ! pair_values(si_imag+2*n_indnz2_this_indnz1,itrx) = &
                                                !                         aimag(0.5*perm_np1*conj_button_np1*&
                                                !                         operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx) = &
                                                                        -dble(0.5*perm_np1*conj_button_np1*&
                                                                        operator_value*dv_dx(leads_np1,el_np1,itrx))
                                            enddo
                                        else
                                            pair_index = unique_indnz2(itr_indnz2,1)
                                            do itrx = 0,len_x_vec-1   
                                                pair_values(si_real+2*pair_index,itrx) = &
                                                                        pair_values(si_real+2*pair_index,itrx) + &
                                                                        dble(0.5*perm_np1*conj_button_np1*&
                                                                        operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                ! pair_values(si_real+2*pair_index+1,itrx) = &
                                                !                         pair_values(si_real+2*pair_index+1,itrx) + &
                                                !                         aimag(0.5*perm_np1*conj_button_np1*&
                                                !                         operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                ! pair_values(si_imag+2*pair_index,itrx) = &
                                                !                         pair_values(si_imag+2*pair_index,itrx) + &
                                                !                         aimag(0.5*perm_np1*conj_button_np1*&
                                                !                         operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                pair_values(si_imag+2*pair_index+1,itrx) = &
                                                                        pair_values(si_imag+2*pair_index+1,itrx) - & 
                                                                        dble(0.5*perm_np1*conj_button_np1*&
                                                                        operator_value*dv_dx(leads_np1,el_np1,itrx))
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
                                                pair_values(si_real+2*n_indnz2_this_indnz1,itrx) = &
                                                        dble(-0.5*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                        operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                ! pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx) = &
                                                !                         aimag(0.5*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                !                         operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                ! pair_values(si_imag+2*n_indnz2_this_indnz1,itrx) = &
                                                !                         aimag(0.5*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                !                         operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx) = &
                                                        -dble(-0.5*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                        operator_value*dv_dx(leads_np1,el_np1,itrx))
                                            enddo
                                        else
                                            do itrx = 0,len_x_vec-1   
                                                pair_index = unique_indnz2(itr_indnz2,1)
                                                pair_values(si_real+2*pair_index,itrx) = &
                                                                        pair_values(si_real+2*pair_index,itrx) + &
                                                                        dble(-0.5*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                        operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                ! pair_values(si_real+2*pair_index+1,itrx) = &
                                                !                         pair_values(si_real+2*pair_index+1,itrx) + &
                                                !                         aimag(0.5*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                !                         operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                ! pair_values(si_imag+2*pair_index,itrx) = &
                                                !                         pair_values(si_imag+2*pair_index,itrx) + &
                                                !                         aimag(0.5*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                !                         operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                pair_values(si_imag+2*pair_index+1,itrx) = &
                                                                        pair_values(si_imag+2*pair_index+1,itrx) - & 
                                                                        dble(-0.5*((-1.d0)**(itrn+1))*perm_np1*conj_button_np1*&
                                                                        operator_value*dv_dx(leads_np1,el_np1,itrx))
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
                                                pair_values(si_real+2*n_indnz2_this_indnz1,itrx) = &
                                                                    dble(0.5*perm_np1*operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                ! pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx) = &
                                                !                         -aimag(0.5*perm_np1*operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                ! pair_values(si_imag+2*n_indnz2_this_indnz1,itrx) = &
                                                !                         aimag(0.5*perm_np1*operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx) = &
                                                                        dble(0.5*perm_np1*operator_value*dv_dx(leads_np1,el_np1,itrx))
                                            enddo
                                        else
                                            do itrx = 0,len_x_vec-1   
                                                pair_index = unique_indnz2(itr_indnz2,1)
                                                pair_values(si_real+2*pair_index,itrx) = &
                                                                        pair_values(si_real+2*pair_index,itrx) + &
                                                                        dble(0.5*perm_np1*operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                ! pair_values(si_real+2*pair_index+1,itrx) = &
                                                !                         pair_values(si_real+2*pair_index+1,itrx) - &
                                                !                         aimag(0.5*perm_np1*operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                ! pair_values(si_imag+2*pair_index,itrx) = &
                                                !                         pair_values(si_imag+2*pair_index,itrx) + &
                                                !                         aimag(0.5*perm_np1*operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                pair_values(si_imag+2*pair_index+1,itrx) = &
                                                                        pair_values(si_imag+2*pair_index+1,itrx) + & 
                                                                        dble(0.5*perm_np1*operator_value*dv_dx(leads_np1,el_np1,itrx))
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
                                                pair_values(si_real+2*n_indnz2_this_indnz1,itrx) = &
                                                                        dble(-0.5*((-1.d0)**(itrn+1))*perm_np1*&
                                                                        operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                ! pair_values(si_real+2*n_indnz2_this_indnz1+1,itrx) = &
                                                !                         -aimag(0.5*((-1.d0)**(itrn+1))*perm_np1*&
                                                !                         operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                ! pair_values(si_imag+2*n_indnz2_this_indnz1,itrx) = &
                                                !                         aimag(0.5*((-1.d0)**(itrn+1))*perm_np1*&
                                                !                         operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                pair_values(si_imag+2*n_indnz2_this_indnz1+1,itrx) = &
                                                                        dble(-0.5*((-1.d0)**(itrn+1))*perm_np1*&
                                                                        operator_value*dv_dx(leads_np1,el_np1,itrx))
                                            enddo
                                        else
                                            pair_index = unique_indnz2(itr_indnz2,1)
                                            do itrx = 0,len_x_vec-1       
                                                pair_values(si_real+2*pair_index,itrx) = &
                                                                        pair_values(si_real+2*pair_index,itrx) + &
                                                                        dble(-0.5*((-1.d0)**(itrn+1))*perm_np1*&
                                                                        operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                ! pair_values(si_real+2*pair_index+1,itrx) = &
                                                !                         pair_values(si_real+2*pair_index+1,itrx) - &
                                                !                         aimag(0.5*((-1.d0)**(itrn+1))*perm_np1*&
                                                !                         operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                ! pair_values(si_imag+2*pair_index,itrx) = &
                                                !                         pair_values(si_imag+2*pair_index,itrx) + &
                                                !                         aimag(0.5*((-1.d0)**(itrn+1))*perm_np1*&
                                                !                         operator_value*dv_dx(leads_np1,el_np1,itrx))
                                                pair_values(si_imag+2*pair_index+1,itrx) = &
                                                                        pair_values(si_imag+2*pair_index+1,itrx) + & 
                                                                        dble(-0.5*((-1.d0)**(itrn+1))*perm_np1*&
                                                                        operator_value*dv_dx(leads_np1,el_np1,itrx))
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

