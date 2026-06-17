! ---------------------------------------------------------------------------
! 
!     GENERATE COEFFICIENTS AND EXPONENTS OF BATH-CORRELATION EXPANSION
!
! ---------------------------------------------------------------------------
!
! This FORTRAN subroutine is meant to be converted to a Python wrapper using f2py, for use 
! in the main Python code. It takes the system parameters and results of the Pade approximation
! and outputs the exponents and coefficients of the power series expansion of the bath-correlation functions,
! which are calculated from residue theory.
!
! USAGE - RUN FROM TERMINAL (LINUX) TO CREATE PYTHON WRAPPER:
!       Check fortran compilers available in your platform:  f2py -c --help-fcompiler
!       For 'ifort': f2py -c -m eta_gamma eta_gamma.f90 --opt='-O3' --fcompiler=intelem --f90flags='-openmp -D__OPENMP' -liomp5
!       For 'gfortran': f2py -c -m eta_gamma eta_gamma.f90 --opt='-O3' --fcompiler=gnu95 --f90flags='-fopenmp -D__OPENMP' -lgomp
!
! USAGE - RUN FROM MAIN PYTHON CODE ONCE WRAPPED:
!       gamma_vec,eta_vec = eta_gamma.genbathcoeff(kappa,zeta,specwidth,muvec,el_lead_couplings,temp,nleads,npoles,nel,nsign)
!
! INPUTS:
!       kappa                               Array of size [1,npoles] containing frequencies from Pade expansion
!
!       zeta                                Array of size [1,npoles] containing roots from Pade expansion
!
!       specwidth                           Array of size [1,nleads] that specifies width of Lorentzian spectral density function
!                                           for the leads
!
!       muvec                               Array of size [1,nleads] that specifies chemical potential of the leads
!
!       el_lead_couplings                   Array of size [nleads,nel] that specifies energy independent coupling (V_{K,m}) between each lead and electronic level
! 
!       temp                                Scalar that specifies temperature of leads (assumed to be the same for all leads)
!
!       nleads                              Scalar that specifies the number of leads
!
!       npoles                              Scalar that specifies the number of Pade poles used in the expansion
!       
!       nel                                 Scalar that specifies the number of electronic levels in the system
!
!       nsign                               Scalar that is always nsign = 2, representing the creation (+=0) and annihilation (-=1) operators
!
! OUTPUTS:
!       eta_vec                             Array of size [nleads,nel,npoles+1] containing the coefficients for each bath-correlation function connecting
!                                           the leads and electronic levels. Although technically it only depends on the Pade pole index, here is where we put
!                                           the system-lead coupling, V_{K,m}, so it needs these two indices as well
!       
!       gamma_vec                           Array of size [nleads,npoles,nsign] containing the exponents of the exponential functions in the bath-correlation
!                                           expansion.

subroutine genbathcoeff(kappa,zeta,specwidth,muvec,temp,eta_vec,gamma_vec,nleads,npoles,nsign)

    implicit none                                                                       ! Prevents Fortran from treating all variables that start with the letters i, j, k, l, m and n
                                                                                        ! as integers and all other variables as real arguments.
    external genfermiapp                                                                ! Import genfermiapp subroutine from below

    integer, intent(in) :: nleads,npoles,nsign                                      ! Define integer input variables
    real*8, intent(in) :: temp                                                          ! Define real scalar input variables    
    real*8, intent(in),dimension(0:npoles-1) :: kappa,zeta                              ! Define real array input variables
    real*8, intent(in),dimension(0:nleads-1) :: specwidth,muvec
    ! real*8, intent(in),dimension(0:nleads-1,0:nel-1) :: el_lead_couplings
    complex*16, intent(out),dimension(0:nleads-1,0:nsign-1,0:npoles) :: gamma_vec       ! Define complex array output variables
    complex*16, intent(out),dimension(0:nleads-1,0:nsign-1,0:npoles) :: eta_vec
    complex*16, parameter :: ci=(0.d0,1.d0)                                             ! Define the imaginary number == sqrt(-1)
    real*8 :: pi=4.D0*DATAN(1.D0)
    integer :: itri,itrj,itrk                                                              ! Define further necessary integer variables
    complex*16 :: F_lower,F_upper                                                                   ! Define further necessary complex variable
    
    gamma_vec = 0.d0                                                                    ! Initialize gamma_vec and eta_vec as 0 with double precision
    eta_vec = 0.d0

    do itri = 0,(nleads-1)                                                              ! Loop through the leads index
        gamma_vec(itri,0,0) = specwidth(itri) - ci*muvec(itri)                          ! Generate the gamma exponent for this lead, the '0th' Pade pole, and sigma = -
        gamma_vec(itri,1,0) = specwidth(itri) + ci*muvec(itri)                          ! Generate the gamma exponent for this lead, the '0th' Pade pole, and sigma = -
        do itrj = 1,npoles                                                              ! Loop through the remaining Pade poles
            gamma_vec(itri,0,itrj) = sqrt(-zeta(itrj-1))*temp - ci*muvec(itri)          ! Generate the gamma exponent for this lead, the itrj-th Pade pole, and sigma = -
            gamma_vec(itri,1,itrj) = sqrt(-zeta(itrj-1))*temp + ci*muvec(itri)          ! Generate the gamma exponent for this lead, the itrj-th Pade pole, and sigma = +
        enddo
        
        call genfermiapp(ci*specwidth(itri)/temp,kappa,zeta,F_upper,npoles)                   ! Calculate f_{approx}(iW)
        call genfermiapp(-ci*specwidth(itri)/temp,kappa,zeta,F_lower,npoles)                   ! Calculate f_{approx}(iW)
        ! do itrj = 0,(nsign-1)                                                             ! Loop through the electron levels
        eta_vec(itri,0,0) = pi*specwidth(itri)*(0.5 + F_upper)
        eta_vec(itri,1,0) = pi*specwidth(itri)*(0.5 - F_lower)
                                                                                    ! For each electron level, generate the eta coefficient of the '0th' Pade pole
        do itrk = 1,npoles                                                          ! Loop through the Pade poles
            eta_vec(itri,0:1,itrk) = -ci*2*pi*temp*kappa(itrk-1) & 
                                    *(specwidth(itri)**2)/(zeta(itrk-1)*(temp**2) + specwidth(itri)**2)
                                                                                    ! Generate eta_{l}*V_{K,m} (K = itri,m = itrj,l = itrk)
        enddo
        ! enddo
    enddo
    
end subroutine genbathcoeff

! ---------------------------------------------------------------------------
! 
!     GENERATE COEFFICIENTS AND EXPONENTS OF BATH-CORRELATION EXPANSION
!                   UNDER THE WIDE-BAND APPROXIMATION
!
! ---------------------------------------------------------------------------
!
! Same as genbathcoeff, except under the wide-band approximation

subroutine genbathcoeff_wbl(kappa,zeta,muvec,temp,eta_vec,gamma_vec,nleads,npoles,nsign)

    implicit none                                                                       ! Prevents Fortran from treating all variables that start with the letters i, j, k, l, m and n
                                                                                        ! as integers and all other variables as real arguments.
    external genfermiapp                                                                ! Import genfermiapp subroutine from below

    integer, intent(in) :: nleads,npoles,nsign                                          ! Define integer input variables
    real*8, intent(in) :: temp                                                          ! Define real scalar input variables    
    real*8, intent(in),dimension(0:npoles-1) :: kappa,zeta                              ! Define real array input variables
    real*8, intent(in),dimension(0:nleads-1) :: muvec
    complex*16, intent(out),dimension(0:nleads-1,0:nsign-1,0:npoles-1) :: gamma_vec       ! Define complex array output variables
    complex*16, intent(out),dimension(0:nleads-1,0:nsign-1,0:npoles-1) :: eta_vec
    complex*16, parameter :: ci=(0.d0,1.d0)                                             ! Define the imaginary number == sqrt(-1)
    real*8, parameter :: pi=dacos(-1.d0)
    integer itri,itrj,itrk                                                              ! Define further necessary integer variables
        
    gamma_vec = 0.d0                                                                    ! Initialize gamma_vec and eta_vec as 0 with double precision
    eta_vec = 0.d0

    do itri = 0,nleads-1                                                              ! Loop through the leads index
        do itrj = 0,npoles-1                                                            ! Loop through the remaining Pade poles
            gamma_vec(itri,0,itrj) = sqrt(-zeta(itrj))*temp - ci*muvec(itri)            ! Generate the gamma exponent for this lead, the itrj-th Pade pole, and sigma = -
            gamma_vec(itri,1,itrj) = sqrt(-zeta(itrj))*temp + ci*muvec(itri)            ! Generate the gamma exponent for this lead, the itrj-th Pade pole, and sigma = +
        enddo
        
        do itrj = 0,(nsign-1)                                                           ! Loop through the electron levels
            do itrk = 0,npoles-1                                                        ! Loop through the Pade poles
                eta_vec(itri,itrj,itrk) = -ci*2*pi*temp*kappa(itrk)                     ! Generate eta_{l}*V_{K,m} (K = itri,m = itrj,l = itrk)
            enddo
        enddo
    enddo
    
end subroutine genbathcoeff_wbl

! ---------------------------------------------------------------------------
!
!     APPROXIMATE FERMI-DIRAC FUNCTION BASED ON PADE DECOMPOSITION
!
! ---------------------------------------------------------------------------
!
! This FORTRAN subroutine is meant to be converted to a Python wrapper using f2py, for use 
! in the main Python code. It takes the system parameters and results of the Pade approximation
! and outputs the exponents and coefficients of the power series expansion of the bath-correlation functions,
! which are calculated from residue theory.
!
! USAGE (IN FORTRAN):
!       call genfermiapp(w,kappa,zeta,F,npoles)


subroutine genfermiapp(w,kappa,zeta,F,npoles)

    implicit none

    complex*16, intent(in) :: w                                                         ! Energy at which to evaluate Fermi-Dirac function
    integer, intent(in) :: npoles                                                       ! Define input: number of Pade poles used in approximation
    real*8, intent(in), dimension(0:npoles-1) :: kappa,zeta                             ! Define input: parameters used in Pade expansion
    ! integer, intent(in) :: len_F
    complex*16, intent(out) :: F                                                        ! Define output: result of Fermi-Dirac approximation
    ! complex*16, intent(out), dimension(len_F) :: F
    ! complex*16, dimension(len_F) :: x
    integer :: itri

    ! x = w/temp
    F = 0.d0
    do itri = 0,npoles-1                                                                ! Implement Pade approximation
        F = F - (2.d0*kappa(itri)*w)/(-zeta(itri) + w**2)
    enddo

end subroutine genfermiapp
