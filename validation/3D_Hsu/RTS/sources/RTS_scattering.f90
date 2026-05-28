module scatteringdata

    !---------------------------------------------------------------------------------------!
    !                                                                                       !
    !              /\\\\\\\\\      /\\\\\\\\\\\\\\\     /\\\\\\\\\\\                        !
    !             /\\\///////\\\   \///////\\\/////    /\\\/////////\\\                     !
    !             \/\\\     \/\\\         \/\\\        \//\\\      \///                     !
    !              \/\\\\\\\\\\\/          \/\\\         \////\\\                           !
    !               \/\\\//////\\\          \/\\\            \////\\\                       !
    !                \/\\\    \//\\\         \/\\\               \////\\\                   !
    !                 \/\\\     \//\\\        \/\\\        /\\\      \//\\\                 !
    !                  \/\\\      \//\\\       \/\\\       \///\\\\\\\\\\\/                 !
    !                   \///        \///        \///          \///////////                  !
    !                              Radiative Transfer Simulator                             !
    !---------------------------------------------------------------------------------------!
    !                      This code was developed with gfortran 10.1                       !
    !                           2020-2023 Gustavo Silva Rodrigues                           !
    !---------------------------------------------------------------------------------------!
    !  Permission is hereby granted, free of charge, to any person obtaining a copy of this !
    !     software and associated documentation files, to deal in  the Software without     !
    !    restriction, including without limitation the rights to use, copy, modify, merge,  !
    !  publish, distribute, sublicense, and/or sell  copies of the Software, and to permit  !
    !  persons to whom the Software is furnished to do so, subject to the following         !
    !  conditions:                                                                          !
    !        The above copyright notice and this permission notice shall be included        !
    !                 in all copies or substantial portions of the Software.                !
    !                                                                                       !
    ! The software is provided "as is" in the hope that it will be useful, WITHOUT WARRANTY !
    !    OF ANY KIND, express or implied, including but not limited to the warranties of    !
    !  MERCHANTABILITY, FITNESS FITNESS and noninfringement. In no event shall the authors  !
    !   or copyright holders be liable for any claim, damages or other liability, whether   !
    !   in an action of contract, tort or otherwise, arising from, out of or in connection  !
    !            with the software or the use or other dealings in the software.            !
    !                                                                                       !
    !   Report errors or contributions via GitHub in https://github.com/gusirosx/RTS        !
    !---------------------------------------------------------------------------------------!
    
    !========================================================================================
    use global
    implicit none
    contains
    !========================================================================================
    !                              Phase function subroutines
    !========================================================================================
    !----------------------------------------------------------
    !           Selects the scattering phase function 
    !---------------------- Description -----------------------
    ! <INPUT>
    ! angle      - Is the angle measured between S and S'
    ! C_rad      - Linear anisotropic phase function coefficient
    ! <OUTPUT>
    ! phase_func - Selected phase function
    ! <EXTERNAL ROUTINES>
    ! legendre_phase - Computes the polynomial Legendre phase function
    !----------------------------------------------------------
    ! ID 1 - Mie Scattering
    ! ID 2 - Rayleigh Scattering
    ! ID 3 - Linear Scattering
    ! ID 4 - Henyey-Greenstein Scattering
    !----------------------------------------------------------
    subroutine phase_functions(phase_func,angle)
        double precision :: phase_func,angle
        if(scat_flag == 1)then
            phase_func = legendre_phase(angle)
        else if(scat_flag == 2)then
            phase_func = 0.75d0*(1.0d0 + angle**2.0d0)
        else if(scat_flag == 3)then
            phase_func = 1.0d0 + C_rad*angle
        else if(scat_flag == 4)then
            phase_func = (1.0d0 - C_rad**2.0d0)/((1.0d0 + C_rad**2.0d0 - 2.0d0*C_rad*angle)**1.5d0)
        end if
    end subroutine phase_functions
    !========================================================================================
    !----------------------------------------------------------
    !      Computes the polynomial Legendre phase function
    !---------------------- Description -----------------------
    ! <INPUT>
    ! angle          - Is the angle measured between S and S'
    ! <OUTPUT>
    ! legendre_phase - Polynomial Legendre phase function
    ! <INTERNAL>
    ! A_phase        - Expansion coefficients array
    ! nl             - Size of A_phase
    ! LSUM           - Legendre summation auxiliary variable
    ! <EXTERNAL ROUTINES>
    ! legendre_coeff - Select expansion coefficients 
    ! Pncos          - Computes the Legendre function
    !----------------------------------------------------------
    function legendre_phase(angle)
        integer :: nl,k
        double precision :: LSUM,angle,legendre_phase
        call legendre_coeff(nl)
        LSUM = 1.0d0
        do k=1,nl
            LSUM = LSUM + A_phase(k)*Pncos(k,angle)
        end do
        legendre_phase = LSUM
    end function legendre_phase
    !========================================================================================
    !----------------------------------------------------------
    ! Select the expansion coefficients for the phase function
    !---------------------- Description -----------------------
    !  For more details see "Effect of anisotropic scattering 
    !       on radiative heat transfer in two-dimensional 
    !        rectangular enclosures" by Kim and Lee (1987)
    ! More expansion coefficients for the phase function can be
    ! found in the sources below:
    ! doi: https://doi.org/10.1016/0017-9310(88)90283-9
    ! doi: https://doi.org/10.1364/JOSA.47.000081
    ! doi: https://doi.org/10.1115/1.3245070
    !----------------------------------------------------------
    ! <OUTPUT>
    ! A_phase - Expansion coefficients array
    ! nl      - Size of A_phase
    !----------------------------------------------------------
    ! ID 1  - F0 Forward linear anisotropic
    ! ID 2  - F1 Forward function 1
    ! ID 3  - F2 Forward function 2
    ! ID 4  - F3 Highly forward 
    ! ID 5  - F4 Highly forward 
    ! ID 6  - F5 Forward 
    ! ID 7  - F6 Forward with small backward lobe
    ! ID 8  - F7 Forward 
    ! ID 9  - B0 Backward linear
    ! ID 10 - B1 Backward
    ! ID 11 - B2 Backward
    ! ID 12 - B3 Backward with forward lobe
    !----------------------------------------------------------
    subroutine legendre_coeff(nl)
        integer :: nl
        !Foward Scattering Phase Functions
        if(SPF_flag == 1)then
            A_phase =(/1.0/)!F0
        else if(SPF_flag == 2)then
            A_phase =(/2.5360217,3.5654900,3.9797626,4.0029206,3.6640084,3.0160117,2.2330437,&
                       1.3025078,0.5346286,0.2013563,0.0547964,0.0109929/) !F1
        else if(SPF_flag == 3)then
            A_phase =(/2.0091653,1.5633900,0.6740690,0.2221484,0.0472529,0.0067132,0.0006743,&
                       0.0000494/)!F2
        else if(SPF_flag == 4)then
            A_phase =(/2.45677,3.84181,5.00553,5.96331,6.72161,7.28629,7.66244,7.85767,7.88060,&
                       7.73731,7.44517,7.00558,6.44908,5.77136,5.01092,4.19407,3.30067,2.43442,&
                       1.67637,1.03522,0.546797,0.241290,0.0887443,0.0275242,0.00729357/) !F3
        else if(SPF_flag == 5)then
            A_phase =(/2.78197,4.25856,5.38653,6.19015,6.74492,7.06711,7.20999,7.20063,7.03629,&
                       6.76587,6.35881,5.83351,5.22997,4.47918,3.69000,2.81577,1.92305,1.11502,&
                       0.50766,0.20927,0.07138,0.02090,0.00535,0.00120,0.00024,0.00004/)  !F4
        else if(SPF_flag == 6)then
            A_phase =(/2.30790,3.23995,3.79936,4.03230,4.01175,3.81199,3.49134,3.10532,2.69153,     &
                       2.28376,1.89851,1.55180,1.24622,0.986923,0.769385,0.592926,0.450530,0.340380,&
                       0.254852,0.189456,0.140870,0.104504,0.0777386,0.0575404,0.0427188,0.0313215, &
                       0.0229587,0.0165931,0.0119847,0.00851234,0.00603327,0.00418646,0.00288973,   &
                       0.00193843,0.00126717/)!F5
        else if(SPF_flag == 7)then
            A_phase =(/1.45687,2.15059,1.85048,2.35120,1.85811,2.16895,1.58175,1.75408,1.21806,     &
                       1.33574,0.890579,0.932577,0.585578,0.554820,0.338866,0.340838,0.185501,      &
                       0.182498,0.111208,0.100875,0.0608943,0.0587712,0.0362599,0.0314965,0.0221769,&
                       0.0177937,0.00953461,0.0105159,0.00504125,0.00468271,0.00244863,0.00211119,  &
                       0.00108605/)!F6
        else if(SPF_flag == 8)then
            A_phase =(/1.5,0.5/)!F7
        !Backward Scattering Phase Functions
        else if(SPF_flag == 9)then
            A_phase =(/-1.0/)!B0
        else if(SPF_flag == 10)then
            A_phase =(/-0.56524,0.29783,0.08571,0.01003,0.00063/)!B1
        else if(SPF_flag == 11)then
            A_phase =(/-1.2,0.5/)!B2
        else if(SPF_flag == 12)then
            A_phase =(/1.0,1.0,-0.75/)!B3
        end if
        nl = size(A_phase)
    end subroutine legendre_coeff
    !========================================================================================
    !----------------------------------------------------------
    !           Computes the renormalized associated  
    !                polynomial function Legendre
    !---------------------- Description -----------------------
    !    For more details see equation (6.7.8) in the book
    !    Numerical Recipes: The Art of Scientific Computing 
    !               3ed by W.H.PRESS et al.(2007)
    ! NOTE: angle must be is in the range from -1 to +1
    !----------------------------------------------------------
    function Pncos(n,angle)
        integer :: n,i
        double precision :: Pncos,angle,pll,pmm,pmmp
        pmm = 1.0
        pmmp = angle
        if(n == 1) then
            Pncos = pmmp
        else
            do i=2,n
                pll = (angle*(2.0d0*i - 1.0d0)*pmmp - (i - 1.0d0)*pmm)/i
                pmm = pmmp
                pmmp = pll
            end do
            Pncos = pll
        end if
    end function Pncos
    !========================================================================================
    !-------------------------------------------------------------
    ! Calculates the anisotropic scattering phase function for FAM
    !------------------------ Description ------------------------
    !    For more details see "Advances in Numerical Heat 
    !          Transfer": Vol 2 Chapter 4 pag. 125
    !          by J.C Chai and S.V. Patankar(2000)
    !-------------------------------------------------------------
    ! <INPUT>
    ! dco        - Integral over the solid angle
    ! nsub       - number of angular subdivisions 
    ! <OUTPUT>
    ! phase_f    - Scattering phase function for FAM
    ! <INTERNAL>
    ! PSUM       - Summation auxiliary variable
    ! angle      - Is the angle measured between S and S'
    ! thetam     - theta value at the subcontrol-angle center
    ! phim       - phi value at the subcontrol-angle center
    ! dphi       - Control-angle length 
    ! d_om       - Delta omega
    ! thetas     - subdivided theta-direction grid 
    ! phis       - subdivided phi-direction grid
    ! mu         - x-direction cossines
    ! eta        - y-direction cossines
    ! xi         - z-direction cossines
    ! phase_func - Phase function value
    ! <EXTERNAL ROUTINES>
    ! sub_angular_grid - Performs subdivision at the control-angles
    ! phase_norm_fam   - Normalizes the scattering phase function
    ! phase_functions  - Calculates the selected phase function
    !-------------------------------------------------------------
    subroutine anisotropic_fam
        integer :: l,m,ll,mm,ls,ms,lls,mms
        double precision :: d_om,d_oms,mu,mus,eta,etas,xi,xis,phim,phims,dphi,dphis,thetam,&
                            angle,PSUM,phase_func,phis(np+1,nsub+1),thetas(nt+1,nsub+1)
        call sub_angular_grid(phis,thetas)
        do m=1,np
            do l=1,nt
                do mm=1,np
                    do ll=1,nt
                        PSUM = 0.0d0
                        do ms=1,nsub
                            phim = 0.5d0*(phis(m,ms) + phis(m,ms+1))
                            dphi = phis(m,ms+1) - phis(m,ms)
                            do ls=1,nsub
                                thetam = 0.5d0*(thetas(l,ls) + thetas(l,ls+1))
                                d_om=-dphi*(cos(thetas(l,ls+1)) - cos(thetas(l,ls)))
                                mu  = sin(thetam)*cos(phim)
                                eta = sin(thetam)*sin(phim)
                                xi  = cos(thetam)
                                do mms=1,nsub
                                    phims = 0.5d0*(phis(mm,mms) + phis(mm,mms+1))
                                    dphis = phis(mm,mms+1) - phis(mm,mms)
                                    do lls=1,nsub
                                        thetam= 0.5d0*(thetas(ll,lls) + thetas(ll,lls+1))
                                        d_oms=-dphis*(cos(thetas(ll,lls+1)) - cos(thetas(ll,lls)))
                                        mus  = sin(thetam)*cos(phims)
                                        etas = sin(thetam)*sin(phims)
                                        xis  = cos(thetam)
                                        angle= mu*mus + eta*etas + xi*xis
                                        call phase_functions(phase_func,angle)
                                        PSUM = PSUM + phase_func*d_oms*d_om
                                    end do
                                end do
                            end do
                        end do
                        phase_f(ll,mm,l,m) = PSUM/(dco(ll,mm)*dco(l,m))
                    end do
                end do
            end do
        end do
        call phase_norm_fam
    end subroutine anisotropic_fam
    !========================================================================================
    !----------------------------------------------------------
    !         Performs subdivision at the control-angles
    !---------------------- Description -----------------------
    !    For more details see "Advances in Numerical Heat 
    !          Transfer": Vol 2 Chapter 4 pag. 125
    !          by J.C Chai and S.V. Patankar(2000)
    !----------------------------------------------------------
    ! <INPUT>
    ! theta   - theta-direction grid 
    ! phi     - phi-direction grid
    ! nsub    - number of angular subdivisions for the phase function
    ! <INTERNAL>
    ! dphi    - Control-angle length 
    ! dtheta  - Control-angle height
    ! <OUTPUT>
    ! thetas  - subdivided theta-direction grid 
    ! phis    - subdivided phi-direction grid
    !----------------------------------------------------------
    subroutine sub_angular_grid(phis,thetas)
        integer :: i,l,m
        double precision :: dphi,dtheta
        double precision :: phis(np+1,nsub+1),thetas(nt+1,nsub+1)
        !----------- phi-direction subdivisions -----------!
        do m=1,np
            dphi = (phi(m+1) - phi(m))/nsub
            do i=1,nsub+1
                phis(m,i) = phi(m) + (i-1.0d0)*dphi
            end do
        end do
        !---------- theta-direction subdivisions ----------!
        do l=1,nt
            dtheta = (theta(l+1) - theta(l))/nsub
            do i=1,nsub+1
                thetas(l,i) = theta(l) + (i-1.0d0)*dtheta
            end do
        end do
        !--------------------------------------------------!
    end subroutine sub_angular_grid
    !========================================================================================
    subroutine phase_norm_fam
    !----------------------------------------------------------
    !     Normalizes the scattering phase function for FAM
    !---------------------- Description -----------------------
    ! The inaccuracy of anisotropic scattering solutions is 
    ! mainly due to the error in the integration of the phase
    ! function. The normalization of the phase function is done 
    ! to minimize such errors.
    !----------------------------------------------------------
    ! <INPUT>
    ! phase   - Scattering phase function
    ! dco     - Integral over the solid angle
    ! <OUTPUT>
    ! phase   - Normalized scattering phase function
    ! <INTERNAL>
    ! PSUM   - Summation auxiliary variable
    !----------------------------------------------------------
        integer :: l,m,ll,mm
        double precision :: PSUM,factor
        do l=1,nt
            do m=1,np
                PSUM = 0.0
                do ll=1,nt
                    do mm=1,np
                        PSUM = PSUM + phase_f(ll,mm,l,m)*dco(ll,mm)
                    end do
                end do
                factor = (PSUM/PI4)
                do ll=1,nt
                    do mm=1,np
                        phase_f(ll,mm,l,m) = phase_f(ll,mm,l,m)/(factor + small)
                    end do
                end do
            end do
        end do
    end subroutine phase_norm_fam
    !========================================================================================
    !--------------------------------------------------------------
    ! Calculates the anisotropic scattering phase function for DOM
    !------------------------ Description -------------------------
    ! <INPUT>
    ! mux        - x ordinate cossine direction quadrature
    ! etay       - y ordinate cossine direction quadrature
    ! xiz        - z ordinate cossine direction quadrature
    ! <OUTPUT>
    ! phase_d    - Normalized scattering phase function
    ! <INTERNAL>
    ! angle      - Is the angle measured between S and S'
    ! phase_func - Phase function value
    ! <EXTERNAL ROUTINES>
    ! phase_functions - Calculates the selected phase function
    ! phase_norm_dom  - Normalizes the scattering phase function
    !--------------------------------------------------------------
    subroutine anisotropic_dom
        integer :: l,ls
        double precision :: angle,phase_func
        do l=1,nq
            do ls=1,nq
                angle = mux(l)*mux(ls) + etay(l)*etay(ls) + xiz(l)*xiz(ls)
                call phase_functions(phase_func,angle)
                phase_d(ls,l) = phase_func
            end do
        end do
        call phase_norm_dom
    end subroutine anisotropic_dom
    !========================================================================================
    !----------------------------------------------------------
    !     Normalizes the scattering phase function for DOM
    !---------------------- Description -----------------------
    ! The inaccuracy of anisotropic scattering solutions is 
    ! mainly due to the error in the integration of the phase
    ! function. The normalization of the phase function is done 
    ! to minimize such errors.
    !----------------------------------------------------------
    ! <INPUT>
    ! phase_d - Scattering phase function
    ! Wq      - Quadrature weights
    ! <OUTPUT>
    ! phase_d - Normalized scattering phase function for DOM
    ! <INTERNAL>
    ! PSUM    - Summation auxiliary variable
    !----------------------------------------------------------
    subroutine phase_norm_dom
        integer :: l,ls
        double precision :: PSUM,factor
        do l=1,nq
            PSUM = 0.0d0
            do ls=1,nq
                PSUM = PSUM + phase_d(ls,l)*Wq(ls)
            end do
            factor = (PSUM/PI4)
            do ls=1,nq
                phase_d(ls,l) = phase_d(ls,l)/(factor + small)
            end do
        end do
    end subroutine phase_norm_dom
    !========================================================================================
end module scatteringdata
