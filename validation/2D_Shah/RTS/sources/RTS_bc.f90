module boundary_conditions

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
    use functions
    implicit none
    contains
    !========================================================================================
    !--------------------------------------------------------------------
    !            Boundary conditions coefficients for P1 method 
    !--------------------------- Description ----------------------------
    ! <INPUT>
    ! GAMA      - Gamma coefficient for P1
    ! T_energy  - Temperature
    ! boltz     - Stefan boltzmann constant
    ! BC_flags_G- Boundary flags for incident radiation
    ! <OUTPUT>
    ! a_bxG     - a boundary coefficient x-direction walls
    ! a_byG     - a boundary coefficient y-direction walls
    ! a_bzG     - a boundary coefficient z-direction walls
    ! b_bxG     - b boundary coefficient x-direction walls
    ! b_byG     - b boundary coefficient y-direction walls
    ! b_bzG     - b boundary coefficient z-direction walls
    ! g_bxG     - g boundary coefficient x-direction walls
    ! g_byG     - g boundary coefficient y-direction walls
    ! g_bzG     - g boundary coefficient z-direction walls
    ! <EXTERNAL ROUTINES>
    ! zetaPONE  - P1 method boundary auxiliary function
    !--------------------------------------------------------------------
    subroutine bound_P1
        ! Declare variables
        integer :: i,j,k
        ! Assign values to boundary condition flags
        BC_flags_G = 3
        !------------------ X Boundaries ------------------!
        do k=2,nzp
            do j=2,nyp
                !West Boundary (x = 0)
                a_bxG(1,j,k) = zetaPONE(1,j,k)
                b_bxG(1,j,k) = - GAMA(2,j,k)
                g_bxG(1,j,k) = 4.0d0*zetaPONE(1,j,k)*boltz*T_energy(1,j,k)**4.0d0
                !East Boundary (x = L)
                a_bxG(2,j,k) = zetaPONE(nxi,j,k)
                b_bxG(2,j,k) = GAMA(nxp,j,k)
                g_bxG(2,j,k) = 4.0d0*zetaPONE(nxi,j,k)*boltz*T_energy(nxi,j,k)**4.0d0
            end do
        end do
        !------------------ Y Boundaries ------------------!
        do k=2,nzp
            do i=2,nxp
                !South Boundary (y = 0)
                a_byG(1,i,k) = zetaPONE(i,1,k)
                b_byG(1,i,k) = - GAMA(i,2,k)
                g_byG(1,i,k) = 4.0d0*zetaPONE(i,1,k)*boltz*T_energy(i,1,k)**4.0d0
                !North Boundary (y = L)
                a_byG(2,i,k) = zetaPONE(i,nyi,k)
                b_byG(2,i,k) = GAMA(i,nyp,k)
                g_byG(2,i,k) = 4.0d0*zetaPONE(i,nyi,k)*boltz*T_energy(i,nyi,k)**4.0d0
            end do
        end do
        !------------------ Z Boundaries ------------------!
        do i=2,nxp
            do j=2,nyp
                !Bottom Boundary (z = 0)
                a_bzG(1,i,j) = zetaPONE(i,j,1)
                b_bzG(1,i,j) = - GAMA(i,j,2)
                g_bzG(1,i,j) = 4.0d0*zetaPONE(i,j,1)*boltz*T_energy(i,j,1)**4.0d0
                !Top Boundary (z = L)
                a_bzG(2,i,j) = zetaPONE(i,j,nzi)
                b_bzG(2,i,j) = GAMA(i,j,nzp)
                g_bzG(2,i,j) = 4.0d0*zetaPONE(i,j,nzi)*boltz*T_energy(i,j,nzi)**4.0d0
            end do
        end do
    end subroutine bound_P1
    !========================================================================================
    !--------------------------------------------------------------------
    !           Boundary conditions coefficients for temperature
    !--------------------------- Description ----------------------------
    ! <INPUT>
    ! BC_flags_t- Boundary flags for temperature
    ! wall_temp - Wall temperature
    ! K_term    - Thermal conductivity
    ! h_conv    - Convective heat transfer coefficient
    ! T_inf     - Room temperature
    ! <OUTPUT>
    ! a_bxT     - a boundary coefficient x-direction walls
    ! a_byT     - a boundary coefficient y-direction walls
    ! a_bzT     - a boundary coefficient z-direction walls
    ! b_bxT     - b boundary coefficient x-direction walls
    ! b_byT     - b boundary coefficient y-direction walls
    ! b_bzT     - b boundary coefficient z-direction walls
    ! g_bxT     - g boundary coefficient x-direction walls
    ! g_byT     - g boundary coefficient y-direction walls
    ! g_bzT     - g boundary coefficient z-direction walls
    !--------------------------------------------------------------------
    subroutine bound_temp
        integer :: i,j,k
        !------------------ X Boundaries ------------------!
        !West Boundary (x = 0)
        if(BC_flags_t(1) == 1)then
            do k=2,nzp
                do j=2,nyp
                    a_bxT(1,j,k) = 1.0d0
                    b_bxT(1,j,k) = 0.0d0
                    g_bxT(1,j,k) = wall_temp(1)
                end do
            end do
        elseif(BC_flags_t(1) == 2)then
            do k=2,nzp
                do j=2,nyp
                    a_bxT(1,j,k) = 0.0d0
                    b_bxT(1,j,k) = 1.0d0
                    g_bxT(1,j,k) = 0.0d0
                end do
            end do
        elseif(BC_flags_t(1) == 3)then
            do k=2,nzp
                do j=2,nyp
                    a_bxT(1,j,k) = -h_conv(1)
                    b_bxT(1,j,k) = K_term(2,j,k)
                    g_bxT(1,j,k) = -h_conv(1)*T_inf
                end do
            end do
        end if
        !East Boundary (x = L)
        if(BC_flags_t(2) == 1)then
            do k=2,nzp
                do j=2,nyp
                    a_bxT(2,j,k) = 1.0d0
                    b_bxT(2,j,k) = 0.0d0
                    g_bxT(2,j,k) = wall_temp(2)
                end do
            end do
        elseif(BC_flags_t(2) == 2)then
            do k=2,nzp
                do j=2,nyp
                    a_bxT(2,j,k) = 0.0d0
                    b_bxT(2,j,k) = 1.0d0
                    g_bxT(2,j,k) = 0.0d0
                end do
            end do
        elseif(BC_flags_t(2) == 3)then
            do k=2,nzp
                do j=2,nyp
                    a_bxT(2,j,k) = h_conv(2)
                    b_bxT(2,j,k) = K_term(nxp,j,k)
                    g_bxT(2,j,k) = h_conv(2)*T_inf
                end do
            end do
        end if
        !------------------ Y Boundaries ------------------!
        !South Boundary (y = 0)
        if(BC_flags_t(3) == 1)then
            do k=2,nzp
                do i=2,nxp
                    a_byT(1,i,k) = 1.0d0
                    b_byT(1,i,k) = 0.0d0
                    g_byT(1,i,k) = wall_temp(3)
                end do
            end do
        elseif(BC_flags_t(3) == 2)then
            do k=2,nzp
                do i=2,nxp
                    a_byT(1,i,k) = 0.0d0
                    b_byT(1,i,k) = 1.0d0
                    g_byT(1,i,k) = 0.0d0
                end do
            end do
        elseif(BC_flags_t(3) == 3)then
            do k=2,nzp
                do i=2,nxp
                    a_byT(1,i,k) = -h_conv(3)
                    b_byT(1,i,k) = K_term(i,2,k)
                    g_byT(1,i,k) = -h_conv(3)*T_inf
                end do
            end do
        end if
        !North Boundary (y = L)
        if(BC_flags_t(4) == 1)then
            do k=2,nzp
                do i=2,nxp
                    a_byT(2,i,k) = 1.0d0
                    b_byT(2,i,k) = 0.0d0
                    g_byT(2,i,k) = wall_temp(4)
                end do
            end do
        elseif(BC_flags_t(4) == 2)then
            do k=2,nzp
                do i=2,nxp
                    a_byT(2,i,k) = 0.0d0
                    b_byT(2,i,k) = 1.0d0
                    g_byT(2,i,k) = 0.0d0
                end do
            end do
        elseif(BC_flags_t(4) == 3)then
            do k=2,nzp
                do i=2,nxp
                    a_byT(2,i,k) = h_conv(4)
                    b_byT(2,i,k) = K_term(i,nyp,k)
                    g_byT(2,i,k) = h_conv(4)*T_inf
                end do
            end do
        end if
        !------------------ Z Boundaries ------------------!
        !Bottom Boundary (z = 0)
        if(BC_flags_t(5) == 1) then
            do i=2,nxp
                do j=2,nyp
                    a_bzT(1,i,j) = 1.0d0
                    b_bzT(1,i,j) = 0.0d0
                    g_bzT(1,i,j) = wall_temp(5)
                end do
            end do
        elseif(BC_flags_t(5) == 2) then
            do i=2,nxp
                do j=2,nyp
                    a_bzT(1,i,j) = 0.0d0
                    b_bzT(1,i,j) = 1.0d0
                    g_bzT(1,i,j) = 0.0d0
                end do
            end do
        elseif(BC_flags_t(5) == 3) then
            do i=2,nxp
                do j=2,nyp
                    a_bzT(1,i,j) = -h_conv(5)
                    b_bzT(1,i,j) = K_term(i,j,2)
                    g_bzT(1,i,j) = -h_conv(5)*T_inf
                end do
            end do
        end if
        !Top Boundary (z = L)
        if(BC_flags_t(6) == 1) then
            do i=2,nxp
                do j=2,nyp
                    a_bzT(2,i,j) = 1.0d0
                    b_bzT(2,i,j) = 0.0d0
                    g_bzT(2,i,j) = wall_temp(6)
                end do
            end do
        elseif(BC_flags_t(6) == 2) then
            do i=2,nxp
                do j=2,nyp
                    a_bzT(2,i,j) = 0.0d0
                    b_bzT(2,i,j) = 1.0d0
                    g_bzT(2,i,j) = 0.0d0
                end do
            end do
        elseif(BC_flags_t(6) == 3) then
            do i=2,nxp
                do j=2,nyp
                    a_bzT(2,i,j) = h_conv(6)
                    b_bzT(2,i,j) = K_term(i,j,nzp)
                    g_bzT(2,i,j) = h_conv(6)*T_inf
                end do
            end do
        end if
    end subroutine bound_temp
    !========================================================================================
    !----------------------------------------------------------
    !      Applies boundary conditions to the coefficients
    !---------------------- Description -----------------------
    ! <INPUT>
    ! D       - Diffusion coefficient of the elliptic equation
    ! dxp     - Control-volume length
    ! dyp     - Control-volume height
    ! dzp     - Control-volume width
    ! a_bx    - a boundary coefficient x-direction walls
    ! a_by    - a boundary coefficient y-direction walls
    ! a_bz    - a boundary coefficient z-direction walls
    ! b_bx    - b boundary coefficient x-direction walls
    ! b_by    - b boundary coefficient y-direction walls
    ! b_bz    - b boundary coefficient z-direction walls
    ! g_bx    - g boundary coefficient x-direction walls
    ! g_by    - g boundary coefficient y-direction walls
    ! g_bz    - g boundary coefficient z-direction walls
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
    ! D_wall  - scalar coefficient at the wall
    ! DEN     - Auxiliary denominator
    ! dxx     - Auxiliary variable
    ! dyy     - Auxiliary variable
    ! dzz     - Auxiliary variable
    !----------------------------------------------------------
    subroutine Boundary(Ab,At,Aw,Ae,As,An,Ap,RHS,a_bx,b_bx,g_bx,a_by,b_by,g_by,a_bz,b_bz,g_bz,D)
        integer :: i,j,k
        double precision :: a_bx(2,nyi,nzi),b_bx(2,nyi,nzi),g_bx(2,nyi,nzi),&
                            a_by(2,nxi,nzi),b_by(2,nxi,nzi),g_by(2,nxi,nzi),&
                            a_bz(2,nxi,nyi),b_bz(2,nxi,nyi),g_bz(2,nxi,nyi),&
                            Aw(nxi,nyi,nzi),As(nxi,nyi,nzi),Ab(nxi,nyi,nzi),&
                            Ae(nxi,nyi,nzi),An(nxi,nyi,nzi),At(nxi,nyi,nzi),&
                            Ap(nxi,nyi,nzi),D(nxi,nyi,nzi),RHS(nxi,nyi,nzi),&
                            DEN,D_wall,dxx,dyy,dzz
        !------------------------------- X Boundaries -------------------------------!
        !West Boundary (x = 0)
        dxx = dxp(2)**2.0d0
        do k=2,nzp
            do j=2,nyp
                DEN = (3.0d0*a_bx(1,j,k)*dxx - 8.0d0*b_bx(1,j,k)*dxp(2))
                D_wall = (3.0d0*D(2,j,k) - D(3,j,k))/2.0d0
                Ap(2,j,k) = Ap(2,j,k) + 9.0d0*(a_bx(1,j,k)/DEN)*D_wall - Aw(2,j,k)
                Ae(2,j,k) = Ae(2,j,k) + (a_bx(1,j,k)/DEN)*D_wall
                RHS(2,j,k)= RHS(2,j,k)- 8.0d0*(g_bx(1,j,k)/DEN)*D_wall
                Aw(2,j,k) = 0.0d0
            end do
        end do
        !East Boundary (x = L)
        dxx = dxp(nxp)**2.0d0
        do k=2,nzp
            do j=2,nyp
                DEN = (3.0d0*a_bx(2,j,k)*dxx + 8.0d0*b_bx(2,j,k)*dxp(nxp))
                D_wall = (3.0d0*D(nxp,j,k) - D(nx,j,k))/2.0d0
                Ap(nxp,j,k) = Ap(nxp,j,k) + 9.0d0*(a_bx(2,j,k)/DEN)*D_wall - Ae(nxp,j,k)
                Aw(nxp,j,k) = Aw(nxp,j,k) + (a_bx(2,j,k)/DEN)*D_wall
                RHS(nxp,j,k)= RHS(nxp,j,k)- 8.0d0*(g_bx(2,j,k)/DEN)*D_wall
                Ae(nxp,j,k) = 0.0d0
            end do
        end do
        !------------------------------- Y Boundaries -------------------------------!
        if(DIMEN == 2 .or. DIMEN == 3)then
            !South Boundary (y = 0)
            dyy = dyp(2)**2.0d0
            do k=2,nzp
                do i=2,nxp
                    DEN = (3.0d0*a_by(1,i,k)*dyy - 8.0d0*b_by(1,i,k)*dyp(2))
                    D_wall = (3.0d0*D(i,2,k) - D(i,3,k))/2.0d0
                    Ap(i,2,k) = Ap(i,2,k) + 9.0d0*(a_by(1,i,k)/DEN)*D_wall - As(i,2,k)
                    An(i,2,k) = An(i,2,k) + (a_by(1,i,k)/DEN)*D_wall
                    RHS(i,2,k)= RHS(i,2,k)- 8.0d0*(g_by(1,i,k)/DEN)*D_wall
                    As(i,2,k) = 0.0d0
                end do
            end do
            !North Boundary (y = L)
            dyy = dyp(nyp)**2.0d0
            do k=2,nzp
                do i=2,nxp
                    DEN = (3.0d0*a_by(2,i,k)*dyy + 8.0d0*b_by(2,i,k)*dyp(nyp))
                    D_wall = (3.0d0*D(i,nyp,k) - D(i,ny,k))/2.0d0
                    Ap(i,nyp,k) = Ap(i,nyp,k) + 9.0d0*(a_by(2,i,k)/DEN)*D_wall - An(i,nyp,k)
                    As(i,nyp,k) = As(i,nyp,k) + (a_by(2,i,k)/DEN)*D_wall
                    RHS(i,nyp,k)= RHS(i,nyp,k)- 8.0d0*(g_by(2,i,k)/DEN)*D_wall
                    An(i,nyp,k) = 0.0d0
                end do
            end do
        end if
        !------------------------------- Z Boundaries -------------------------------!
        if(DIMEN == 3)then
            !Bottom Boundary (z = 0)
            dzz = dzp(2)**2.0d0
            do j=2,nyp
                do i=2,nxp
                    DEN = (3.0d0*a_bz(1,i,j)*dzz - 8.0d0*b_bz(1,i,j)*dzp(2))
                    D_wall = (3.0d0*D(i,j,2) - D(i,j,3))/2.0d0
                    Ap(i,j,2) = Ap(i,j,2) + 9.0d0*(a_bz(1,i,j)/DEN)*D_wall - Ab(i,j,2)
                    At(i,j,2) = At(i,j,2) + (a_bz(1,i,j)/DEN)*D_wall
                    RHS(i,j,2)= RHS(i,j,2)- 8.0d0*(g_bz(1,i,j)/DEN)*D_wall
                    Ab(i,j,2) = 0.0d0
                end do
            end do
            !Top Boundary (z = L)
            dzz = dzp(nzp)**2.0d0
            do j=2,nyp
                do i=2,nxp
                    DEN = (3.0d0*a_bz(2,i,j)*dzz + 8.0d0*b_bz(2,i,j)*dzp(nzp))
                    D_wall = (3.0d0*D(i,j,nzp) - D(i,j,nz))/2.0d0
                    Ap(i,j,nzp) = Ap(i,j,nzp) + 9.0d0*(a_bz(2,i,j)/DEN)*D_wall - At(i,j,nzp)
                    Ab(i,j,nzp) = At(i,j,nzp) + (a_bz(2,i,j)/DEN)*D_wall
                    RHS(i,j,nzp)= RHS(i,j,nzp)- 8.0d0*(g_bz(2,i,j)/DEN)*D_wall
                    At(i,j,nzp) = 0.0d0
                end do
            end do
        end if
    end subroutine Boundary
    !========================================================================================
    !--------------------------------------------------------------
    ! Calculates the value of the transported variables at the walls
    !------------------------ Description -------------------------
    ! <INPUT>
    ! dxp     - Control-volume length
    ! dyp     - Control-volume height
    ! dzp     - Control-volume width
    ! a_bx    - a boundary coefficient x-direction walls
    ! a_by    - a boundary coefficient y-direction walls
    ! a_bz    - a boundary coefficient z-direction walls
    ! b_bx    - b boundary coefficient x-direction walls
    ! b_by    - b boundary coefficient y-direction walls
    ! b_bz    - b boundary coefficient z-direction walls
    ! g_bx    - g boundary coefficient x-direction walls
    ! g_by    - g boundary coefficient y-direction walls
    ! g_bz    - g boundary coefficient z-direction walls
    ! <OUTPUT>
    ! phi     - Transported scalar
    ! <INTERNAL>
    ! DEN     - Auxiliary denominator
    ! NUM     - Auxiliary numerator
    ! aux     - Auxiliary variable
    !--------------------------------------------------------------
    subroutine wall_properties(phi,a_bx,b_bx,g_bx,a_by,b_by,g_by,a_bz,b_bz,g_bz)
        integer          :: i,j,k
        double precision :: a_bx(2,nyi,nzi),b_bx(2,nyi,nzi),g_bx(2,nyi,nzi),&
                            a_by(2,nxi,nzi),b_by(2,nxi,nzi),g_by(2,nxi,nzi),&
                            a_bz(2,nxi,nyi),b_bz(2,nxi,nyi),g_bz(2,nxi,nyi),&
                            phi(nxi,nyi,nzi),DEN,NUN,aux
        !------------------------------- X Boundaries -------------------------------!
        !West Boundary (x = 0)
        do k=2,nzp
            do j=2,nyp                
                DEN = a_bx(1,j,k) - (8.0d0*b_bx(1,j,k))/(3.0d0*dxp(2))
                aux = b_bx(1,j,k)/(3.0d0*dxp(2))
                NUN = g_bx(1,j,k) - aux*(9.0d0*phi(2,j,k) - phi(3,j,k))
                phi(1,j,k) = NUN/DEN
            end do
        end do
        !East Boundary (x = L)
        do k=2,nzp
            do j=2,nyp 
                DEN = a_bx(2,j,k) + (8.0d0*b_bx(2,j,k))/(3.0d0*dxp(nxp))
                aux = b_bx(2,j,k)/(3.0d0*dxp(nxp))
                NUN = g_bx(2,j,k) + aux*(9.0d0*phi(nxp,j,k) - phi(nx,j,k))
                phi(nxi,j,k) = NUN/DEN
            end do
        end do
        !------------------------------- Y Boundaries -------------------------------!
        !South Boundary (y = 0)
        do k=2,nzp
            do i=2,nxp
                DEN = a_by(1,i,k) - (8.0d0*b_by(1,i,k))/(3.0d0*dyp(2))
                aux = b_by(1,i,k)/(3.0d0*dyp(2))
                NUN = g_by(1,i,k) - aux*(9.0d0*phi(i,2,k) - phi(i,3,k))
                phi(i,1,k) = NUN/DEN
            end do
        end do
        !North Boundary (y = L)
        do k=2,nzp
            do i=2,nxp
                DEN = a_by(2,i,k) + (8.0d0*b_by(2,i,k))/(3.0d0*dyp(nyp))
                aux = b_by(2,i,k)/(3.0d0*dyp(nyp))
                NUN = g_by(2,i,k) + aux*(9.0d0*phi(i,nyp,k) - phi(i,ny,k))
                phi(i,nyi,k) = NUN/DEN
            end do
        end do
        !------------------------------- Z Boundaries -------------------------------!
        !Bottom Boundary (z = 0)
        do j=2,nyp
            do i=2,nxp
                DEN = a_bz(1,i,j) - (8.0d0*b_bz(1,i,j))/(3.0d0*dzp(2))
                aux = b_bz(1,i,j)/(3.0d0*dzp(2))
                NUN = g_bz(1,i,j) - aux*(9.0d0*phi(i,j,2) - phi(i,j,3))
                phi(i,j,1) = NUN/DEN
            end do
        end do
        !Top Boundary (z = L)
        do j=2,nyp
            do i=2,nxp
                DEN = a_bz(2,i,j) + (8.0d0*b_bz(2,i,j))/(3.0d0*dzp(nzp))
                aux = b_bz(2,i,j)/(3.0d0*dzp(nzp))
                NUN = g_bz(2,i,j) + aux*(9.0d0*phi(i,j,nzp) - phi(i,j,nz))
                phi(i,j,nzi) = NUN/DEN
            end do
        end do
    end subroutine wall_properties
    !========================================================================================
    !--------------------------------------------------------------
    !       Calculates the scalars variables at the walls
    !------------------------ Description -------------------------
    ! <INPUT>
    ! D - Scalar coefficient
    ! <OUTPUT>
    ! D - Scalar coefficient with updated boundaries
    !--------------------------------------------------------------
    subroutine wall_scalars(D)
        integer          :: i,j,k
        double precision :: D(nxi,nyi,nzi)
        !------------------------------- X Boundaries -------------------------------!
        !West Boundary (x = 0)
        do k=2,nzp
            do j=2,nyp
                D(1,j,k) = (3.0d0*D(2,j,k) - D(3,j,k))*0.5d0
            end do
        end do
        !East Boundary (x = L)
        do k=2,nzp
            do j=2,nyp
                D(nxi,j,k) = (3.0d0*D(nxp,j,k) - D(nx,j,k))*0.5d0
            end do
        end do
        !------------------------------- Y Boundaries -------------------------------!
        !South Boundary (y = 0)
        do k=2,nzp
            do i=2,nxp
                D(i,1,k) = (3.0d0*D(i,2,k) - D(i,3,k))*0.5d0
            end do
        end do
        !North Boundary (y = L)
        do k=2,nzp
            do i=2,nxp
                D(i,nyi,k) = (3.0d0*D(i,nyp,k) - D(i,ny,k))*0.5d0
            end do
        end do
        !------------------------------- Z Boundaries -------------------------------!
        !Bottom Boundary (z = 0)
        do j=2,nyp
            do i=2,nxp
                D(i,j,1) = (3.0d0*D(i,j,2) - D(i,j,3))*0.5d0
            end do
        end do
        !Top Boundary (z = L)
        do j=2,nyp
            do i=2,nxp
                D(i,j,nzi) = (3.0d0*D(i,j,nzp) - D(i,j,nz))*0.5d0
            end do
        end do
    end subroutine wall_scalars
    !========================================================================================
    !--------------------------------------------------------------------
    ! Applies the wall boundary values to the temperature field variable
    !--------------------------- Description ----------------------------
    ! <INPUT>
    ! wall_temp - Wall temperature
    ! nuffwall  - Non-uniform field flag on the wall
    ! <OUTPUT>
    ! T_energy  - Temperature field variable
    ! <EXTERNAL ROUTINES>
    ! wallTemp_west   - Set the non-uniform temperature at the west wall
    ! wallTemp_east   - Set the non-uniform temperature at the east wall
    ! wallTemp_south  - Set the non-uniform temperature at the south wall
    ! wallTemp_north  - Set the non-uniform temperature at the north wall
    ! wallTemp_bottom - Set the non-uniform temperature at the bottom wall
    ! wallTemp_top    - Set the non-uniform temperature at the top wall
    !--------------------------------------------------------------------
    subroutine temp_walls()
        integer :: i,j,k
        !------------------ X Boundaries ------------------!
        !West Boundary (x = 0)
        if (nuffwall(1))then
            T_energy(1,2:nyp,2:nzp) = wall_temp(1)  !constant throughout
        else
            do k=2,nzp
                do j=2,nyp
                    T_energy(1,j,k) = wallTemp_west(j,k)
                end do
            end do
        end if
        !East Boundary (x = L)
        if (nuffwall(2))then
            T_energy(nxi,2:nyp,2:nzp) = wall_temp(2)!constant throughout
        else
            do k=2,nzp
                do j=2,nyp
                    T_energy(nxi,j,k) = wallTemp_east(j,k)
                end do
            end do
        end if
        !------------------ Y Boundaries ------------------!
        !South Boundary (y = 0)
        if (nuffwall(3))then
            T_energy(2:nxp,1,2:nzp) = wall_temp(3)  !constant throughout
        else
            do k=2,nzp
                do i=2,nxp
                    T_energy(i,1,k) = wallTemp_south(i,k)
                end do
            end do
        end if
        !North Boundary (y = L)
        if (nuffwall(4))then
            T_energy(2:nxp,nyi,2:nzp) = wall_temp(4)!constant throughout
        else
            do k=2,nzp
                do i=2,nxp
                    T_energy(i,nyi,k) = wallTemp_north(i,k)
                end do
            end do
        end if
        !------------------ Z Boundaries ------------------!
        !Bottom Boundary (z = 0)
        if (nuffwall(5))then
            T_energy(2:nxp,2:nyp,1) = wall_temp(5)  !constant throughout
        else
            do i=2,nxp
                do j=2,nyp
                    T_energy(i,j,1) = wallTemp_bottom(i,j)
                end do
            end do
        end if
        !Top Boundary (z = L)
        if (nuffwall(6))then
            T_energy(2:nxp,2:nyp,nzi) = wall_temp(6)!constant throughout
        else
            do i=2,nxp
                do j=2,nyp
                    T_energy(i,j,nzi) = wallTemp_top(i,j)
                end do
            end do
        end if
    end subroutine temp_walls
    !========================================================================================
    !--------------------------------------------------------------------
    ! Applies the wall boundary values to the emissivity field variable
    !--------------------------- Description ----------------------------
    ! <INPUT>
    ! epsilon_w   - Constant radiation walls emissivities
    ! nuffwall    - Non-uniform field flag on the wall
    ! <OUTPUT>
    ! epsilon_rad - emissivities at the walls and domain(gas)
    ! <EXTERNAL ROUTINES>
    ! wallTemp_west   - Set the non-uniform emissivity at the west wall
    ! wallTemp_east   - Set the non-uniform emissivity at the east wall
    ! wallTemp_south  - Set the non-uniform emissivity at the south wall
    ! wallTemp_north  - Set the non-uniform emissivity at the north wall
    ! wallTemp_bottom - Set the non-uniform emissivity at the bottom wall
    ! wallTemp_top    - Set the non-uniform emissivity at the top wall
    !--------------------------------------------------------------------
    subroutine emissivity_walls()
        integer :: i,j,k
        !------------------ X Boundaries ------------------!
        !West Boundary (x = 0)
        if (nuffwall(1))then
            epsilon_rad(1,2:nyp,2:nzp) = epsilon_w(1)  !constant throughout
        else
            do k=2,nzp
                do j=2,nyp
                    epsilon_rad(1,j,k) = wallemissivity_west(j,k)
                end do
            end do
        end if
        !East Boundary (x = L)
        if (nuffwall(2))then
            epsilon_rad(nxi,2:nyp,2:nzp) = epsilon_w(2)!constant throughout
        else
            do k=2,nzp
                do j=2,nyp
                    epsilon_rad(nxi,j,k) = wallemissivity_east(j,k)
                end do
            end do
        end if
        !------------------ Y Boundaries ------------------!
        !South Boundary (y = 0)
        if (nuffwall(3))then
            epsilon_rad(2:nxp,1,2:nzp) = epsilon_w(3)  !constant throughout
        else
            do k=2,nzp
                do i=2,nxp
                    epsilon_rad(i,1,k) = wallemissivity_south(i,k)
                end do
            end do
        end if
        !North Boundary (y = L)
        if (nuffwall(4))then
            epsilon_rad(2:nxp,nyi,2:nzp) = epsilon_w(4)!constant throughout
        else
            do k=2,nzp
                do i=2,nxp
                    epsilon_rad(i,nyi,k) = wallemissivity_north(i,k)
                end do
            end do
        end if
        !------------------ Z Boundaries ------------------!
        !Bottom Boundary (z = 0)
        if (nuffwall(5))then
            epsilon_rad(2:nxp,2:nyp,1) = epsilon_w(5)  !constant throughout
        else
            do i=2,nxp
                do j=2,nyp
                    epsilon_rad(i,j,1) = wallemissivity_bottom(i,j)
                end do
            end do
        end if
        !Top Boundary (z = L)
        if (nuffwall(6))then
            epsilon_rad(2:nxp,2:nyp,nzi) = epsilon_w(6)!constant throughout
        else
            do i=2,nxp
                do j=2,nyp
                    epsilon_rad(i,j,nzi) = wallemissivity_top(i,j)
                end do
            end do
        end if
    end subroutine emissivity_walls
    !========================================================================================
    !----------------------------------------------------------
    !           P1 method boundary auxiliary function
    !---------------------- Description -----------------------
    ! <INPUT>
    ! i   - x-direction index
    ! j   - y-direction index
    ! k   - z-direction index
    ! epsilon_rad - emissivities at the walls and domain(gas)
    ! <OUTPUT>
    ! zetaPONE - P1 method boundary auxiliary variable
    !----------------------------------------------------------
    function zetaPONE(i,j,k)
        implicit none
        integer :: i,j,k
        double precision :: zetaPONE
        
        zetaPONE = (0.5d0*epsilon_rad(i,j,k))/(2.0d0-epsilon_rad(i,j,k))
        
    end function zetaPONE
    !========================================================================================
end module boundary_conditions
