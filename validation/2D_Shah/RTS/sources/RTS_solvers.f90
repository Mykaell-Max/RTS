module solvers
    
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
    !----------------------------------------------------------
    !            Calculates the internal coefficients 
    !                  of the elliptic equation
    !---------------------- Description -----------------------
    ! <INPUT>
    !           Coefficient originally present in the RHS
    ! A_RHS   - that must be added to Ap if it depends on
    !           the variable transported, like:
    !           Ae*PHIe + Aw*PHIw - Ap*PHIp = RHS + A_RHS*PHIp
    !
    ! A_dt    - LHS coefficient of the transient elliptic equation
    ! D       - Diffusion coefficient of the elliptic equation
    ! dxp     - Control-volume length
    ! dyp     - Control-volume height
    ! dzp     - Control-volume width 
    ! <OUTPUT>
    ! Aw      - West   based coefficient for matrix A
    ! Ae      - East   based coefficient for matrix A
    ! As      - South  based coefficient for matrix A
    ! An      - North  based coefficient for matrix A
    ! Ab      - Bottom based coefficient for matrix A
    ! At      - Top    based coefficient for matrix A
    ! Ap      - Center based coefficient for matrix A
    ! RHS     - Right hand side for matrix B
    ! <INTERNAL>
    ! dxx    - Auxiliary variable
    ! dyy    - Auxiliary variable
    ! dzz    - Auxiliary variable
    !----------------------------------------------------------
    subroutine elliptic_coefficients(D,A_dt,A_RHS,Ab,At,Aw,Ae,As,An,Ap)
        integer          :: i,j,k,IM,IP,JM,JP,KM,KP
        double precision :: Aw(nxi,nyi,nzi),As(nxi,nyi,nzi),Ab(nxi,nyi,nzi),  &
                            Ae(nxi,nyi,nzi),An(nxi,nyi,nzi),At(nxi,nyi,nzi),  &
                            Ap(nxi,nyi,nzi),D(nxi,nyi,nzi),A_RHS(nxi,nyi,nzi),&
                            A_dt,dxx,dyy,dzz
        do k = 2,nzp
            KM = k - 1 + int(2/k)          !k-1 with a minimum of  2  at k = 2
            KP = k + 1 - int(k/nzp)        !k+1 with a maximum of nzp at k = nzp
            dzz = dzp(k)**2.0d0
            do j = 2,nyp
                JM= j - 1 + int(2/j)       !j-1 with a minimum of  2  at j = 2 
                JP= j + 1 - int(j/nyp)     !j+1 with a maximum of nyp at j = nyp
                dyy = dyp(j)**2.0d0
                do i = 2,nxp
                    IM = i - 1 + int(2/i)  !i-1 with a minimum of  2  at i = 2 
                    IP = i + 1 - int(i/nxp)!i+1 with a maximum of nxp at i = nxp
                    dxx = dxp(i)**2.0d0
                    Aw(i,j,k) = A_dt*(D(IM,j,k) + D(i,j,k))/(2.0d0*dxx)!coefficient for phi(i-1,j,k)
                    Ae(i,j,k) = A_dt*(D(i,j,k) + D(IP,j,k))/(2.0d0*dxx)!coefficient for phi(i+1,j,k)
                    As(i,j,k) = A_dt*(D(i,JM,k) + D(i,j,k))/(2.0d0*dyy)!coefficient for phi(i,j-1,k)
                    An(i,j,k) = A_dt*(D(i,j,k) + D(i,JP,k))/(2.0d0*dyy)!coefficient for phi(i,j+1,k)
                    Ab(i,j,k) = A_dt*(D(i,j,KM) + D(i,j,k))/(2.0d0*dzz)!coefficient for phi(i,j,k-1)
                    At(i,j,k) = A_dt*(D(i,j,k) + D(i,j,KP))/(2.0d0*dzz)!coefficient for phi(i,j,k+1)
                    Ap(i,j,k) = Aw(i,j,k) + Ae(i,j,k) + As(i,j,k) + &   !coefficient for phi(i,j,k)
                                An(i,j,k) + Ab(i,j,k) + At(i,j,k) + A_RHS(i,j,k)
                end do
            end do
        end do
    end subroutine elliptic_coefficients
    !========================================================================================
    !----------------------------------------------------------
    !            Successive Over Relaxation (SOR)
    !   Solve a system of linear equations of the type A*x=B
    !           The system is written in the form: 
    !             Ae*PHIe + Aw*PHIw - Ap*PHIp = RHS
    !---------------------- Description -----------------------
    ! <INPUT>
    ! phi     - Initial guess or data from previous call
    ! Aw      - West   based coefficient for matrix A
    ! Ae      - East   based coefficient for matrix A
    ! As      - South  based coefficient for matrix A
    ! An      - North  based coefficient for matrix A
    ! Ab      - Bottom based coefficient for matrix A
    ! At      - Top    based coefficient for matrix A
    ! Ap      - Center based coefficient for matrix A
    ! RHS     - Right hand side for matrix B
    ! ITMAX   - Max number of allowed iterations
    ! omega   - The over-ralaxation factor
    ! sc      - Stopping criteria
    ! <OUTPUT>
    ! phi     - Solution
    ! eps     - Convergence tolerance residue 
    ! nitac   - Number of iterations to achieve the convergence
    ! <INTERNAL>
    ! phi_old - data from previous iteration
    ! phip    - Auxiliary variable
    ! <EXTERNAL ROUTINES>
    ! ETANORM - Residue assessment
    !----------------------------------------------------------
    subroutine SOR(phi,Ab,At,Aw,Ae,As,An,Ap,RHS,nitac,eps,sc,omega)
        integer          :: i,j,k,IM,IP,JM,JP,KM,KP,nitac
        double precision :: Aw(nxi,nyi,nzi),As(nxi,nyi,nzi),Ab(nxi,nyi,nzi),  &
                            Ae(nxi,nyi,nzi),An(nxi,nyi,nzi),At(nxi,nyi,nzi),  &
                            Ap(nxi,nyi,nzi),phi(nxi,nyi,nzi),RHS(nxi,nyi,nzi),&
                            phi_old(nxi,nyi,nzi),phip,eps,sc,omega
        do nitac = 1,ITMAX 
            phi_old = phi
            do k = 2,nzp
                KM = k - 1 + int(2/k)
                KP = k + 1 - int(k/nzp)
                do j = 2,nyp
                    JM = j - 1 + int(2/j)
                    JP = j + 1 - int(j/nyp)
                    do i = 2,nxp
                        IM = I - 1 + int(2/i)
                        IP = I + 1 - int(i/nxp)
                        phip = As(i,j,k)*phi(i,JM,k) + Aw(i,j,k)*phi(IM,j,k) + Ab(i,j,k)*phi(i,j,KM) + &
                        An(i,j,k)*phi(i,JP,k) + Ae(i,j,k)*phi(IP,j,k) + At(i,j,k)*phi(i,j,KP) - RHS(i,j,k)
                        phi(i,j,k) = phip*(omega/Ap(i,j,k)) + (1.0d0 - omega)*phi(i,j,k)
                    end do
                end do
            end do
            eps = ETANORM(phi_old,phi)
            if (eps < sc) exit
        end do
        if(nitac > ITMAX) write(*,'(a)') 'Convergence has not been achieved!'
    end subroutine SOR
    !========================================================================================
    !----------------------------------------------------------
    !                     Jacobi’s Method
    !   Solve a system of linear equations of the type A*x=B
    !           The system is written in the form: 
    !             AePHIe + AwPHIw - ApPHIp = RHS
    !---------------------- Description -----------------------
    ! <INPUT>
    ! phi     - Initial guess or data from previous call
    ! Aw      - West   based coefficient for matrix A
    ! Ae      - East   based coefficient for matrix A
    ! As      - South  based coefficient for matrix A
    ! An      - North  based coefficient for matrix A
    ! Ab      - Bottom based coefficient for matrix A
    ! At      - Top    based coefficient for matrix A
    ! Ap      - Center based coefficient for matrix A
    ! RHS     - Right hand side for matrix B
    ! ITMAX   - Max number of allowed iterations
    ! sc      - Stopping criteria
    ! <OUTPUT>
    ! phi     - Solution
    ! eps     - Convergence tolerance 
    ! nitac   - Number of iterations to achieve the convergence
    ! <INTERNAL>
    ! phi_old - data from previous iteration
    ! phip    - Auxiliary variable
    ! <EXTERNAL ROUTINES>
    ! ETANORM - Residue assessment
    !----------------------------------------------------------
    subroutine JACOBI(phi,Ab,At,Aw,Ae,As,An,Ap,RHS,nitac,eps,sc)
        integer          :: i,j,k,IM,IP,JM,JP,KM,KP,nitac
        double precision :: Aw(nxi,nyi,nzi),As(nxi,nyi,nzi),Ab(nxi,nyi,nzi),  &
                            Ae(nxi,nyi,nzi),An(nxi,nyi,nzi),At(nxi,nyi,nzi),  &
                            Ap(nxi,nyi,nzi),phi(nxi,nyi,nzi),RHS(nxi,nyi,nzi),&
                            phi_old(nxi,nyi,nzi),phip,eps,sc
        do nitac = 1,ITMAX 
            phi_old = phi
            do k = 2,nzp
                KM = k - 1 + int(2/k)
                KP = k + 1 - int(k/nzp)
                do j = 2,nyp
                    JM = j - 1 + int(2/j)
                    JP = j + 1 - int(j/nyp)
                    do i = 2,nxp
                        IM = I - 1 + int(2/i)
                        IP = I + 1 - int(i/nxp)
                        phip = As(i,j,k)*phi_old(i,JM,k) + Aw(i,j,k)*phi_old(IM,j,k) + Ab(i,j,k)*phi_old(i,j,KM) + &
                        An(i,j,k)*phi_old(i,JP,k) + Ae(i,j,k)*phi_old(IP,j,k) + At(i,j,k)*phi_old(i,j,KP) - RHS(i,j,k)
                        phi(i,j,k) = (phip/Ap(i,j,k))
                    end do
                end do
            end do
            eps = ETANORM(phi_old,phi)
            if (eps <= sc) exit
        end do
        if(nitac > ITMAX) write(*,'(a)') 'Convergence has not been achieved!'
    end subroutine JACOBI
    !========================================================================================
    !----------------------------------------------------------
    !                         Norm L2
    !           Find the residual between the current 
    !           and the past iteration data by norm L2
    !---------------------- Description -----------------------
    ! <INPUT>
    ! old     - data from previous iteration
    ! ruling  - data from  atual   iteration
    ! <OUTPUT>
    ! l2norm  - Norm L2
    ! <INTERNAL>
    ! Vsum    - Auxiliary variable
    !----------------------------------------------------------
    function L2NORM(ruling,old)
        implicit none
        integer :: i,j,k
        double precision :: ruling(nxi,nyi,nzi),old(nxi,nyi,nzi),l2norm,Vsum
        Vsum = 0.0d0
        do k = 2,nzp
            do j = 2,nyp
                do i = 2,nxp 
                    Vsum =  Vsum + (ruling(i,j,k) - old(i,j,k))**2.0d0
                end do
            end do
        end do
        l2norm = sqrt(Vsum/(nx*ny*nz))
    end function L2NORM
    !========================================================================================
    !----------------------------------------------------------
    !                         Norm L1
    !           Find the residual between the current 
    !          and the past iteration data by norm L1
    !---------------------- Description -----------------------
    ! <INPUT>
    ! old     - data from previous iteration
    ! ruling  - data from  atual   iteration
    ! <OUTPUT>
    ! l1norm  - Norm L1
    ! <INTERNAL>
    ! Vsum    - Auxiliary variable
    !----------------------------------------------------------
    function L1NORM(ruling,old)
        implicit none
        integer :: i,j,k
        double precision :: ruling(nxi,nyi,nzi),old(nxi,nyi,nzi),l1norm,Vsum
        Vsum = 0.0d0
        do k = 2,nzp
            do j = 2,nyp
                do i = 2,nxp
                    Vsum = Vsum + abs(ruling(i,j,k) - old(i,j,k))
                end do
            end do
        end do
        l1norm = Vsum/(nx*ny*nz)
    end function L1NORM
    !========================================================================================
    !----------------------------------------------------------
    !                         Eta Norm
    !           Find the residual between the current 
    !           and the past iteration data by eta norm
    !---------------------- Description -----------------------
    ! <INPUT>
    ! old     - data from previous iteration
    ! ruling  - data from  atual   iteration
    ! <OUTPUT>
    ! etanorm  - Eta Norm
    ! <INTERNAL>
    ! Vsum    - Auxiliary variable
    !----------------------------------------------------------
    function ETANORM(ruling,old)
        implicit none
        integer :: i,j,k
        double precision :: ruling(nxi,nyi,nzi),old(nxi,nyi,nzi),etanorm,Vsum
        Vsum = 0.0d0
        do k = 2,nzp
            do j = 2,nyp
                do i = 2,nxp
                    Vsum = Vsum + (ruling(i,j,k) - old(i,j,k))**2.0d0
                end do
            end do
        end do
        etanorm = Vsum/MAXVAL(ruling)
    end function ETANORM
    !========================================================================================
    !----------------------------------------------------------
    !                         Max Norm
    !           Find the residual between the current 
    !           and the past iteration data by eta norm
    !---------------------- Description -----------------------
    ! <INPUT>
    ! old     - data from previous iteration
    ! ruling  - data from  atual   iteration
    ! <OUTPUT>
    ! maxnorm - Max Norm
    ! <INTERNAL>
    ! Vsum    - Auxiliary variable
    !----------------------------------------------------------
    function MAXNORM(ruling,old)
        implicit none
        double precision :: ruling(nxi,nyi,nzi),old(nxi,nyi,nzi),maxnorm
        maxnorm = maxval(abs(ruling - old))
    end function MAXNORM
    !========================================================================================
    !----------------------------------------------------------
    !                      Relative Norm
    !           Find the residual between the current 
    !        and the past iteration data by relative norm
    !---------------------- Description -----------------------
    ! <INPUT>
    ! old     - Data from previous iteration
    ! ruling  - Data from  atual   iteration
    ! <OUTPUT>
    ! difnorm - Relative Norm
    ! <INTERNAL>
    ! small   - A small Number
    ! dif     - Auxiliary variable
    !----------------------------------------------------------
    function DIFNORM(ruling,old)
        implicit none
        integer :: i,j,k
        double precision :: ruling(nxi,nyi,nzi),old(nxi,nyi,nzi),difnorm,dif
        difnorm = 0.0d0
        do k = 2,nzp
            do j = 2,nyp
                do i = 2,nxp
                    dif = abs(ruling(i,j,k) - old(i,j,k) )/(ruling(i,j,k) + small)
                    difnorm = max(difnorm,DIF)
                end do
            end do
        end do
    end function DIFNORM
    !========================================================================================
end module solvers
