module radiation
    
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
    use solvers
    use boundary_conditions
    use absorptiondata
    implicit none
    contains
    !========================================================================================
    !-----------------------------------------------------------------------
    !               Solves the Radiative Transfer Equation
    !------------------------------ Description ----------------------------
    ! <EXTERNAL ROUTINES>
    ! radiative_properties - Computes the radiative properties of the medium
    ! P1                   - Compute radiation via P1 Method
    ! DOM                  - Compute radiation via DOM Method
    ! FAM                  - Compute radiation via FAM Method
    !-----------------------------------------------------------------------
    subroutine rtesolve
        if(gas_prop) then
            call radiative_properties
        end if        
        if(model_stc == 1) then
            call P1
        else if(model_stc == 2) then
            call DOM
        else if(model_stc == 3) then
            call FAM
        end if        
    end subroutine rtesolve
    !========================================================================================
    !                                       P1 Method
    !========================================================================================
    !-----------------------------------------------------------------
    !                            P1 Method
    !-------------------------- Description --------------------------
    ! <EXTERNAL ROUTINES>
    ! P1_nongray   - P1 Method using the non-gray approach
    ! P1_gray      - P1 Method using the gray approach
    !-----------------------------------------------------------------
    subroutine P1()
        if(nongray_flag)then
            call P1_nongray
        else
            call P1_gray
        end if
    end subroutine P1
    !========================================================================================
    !-----------------------------------------------------------------------
    !                 P1 Method using the non-gray approach
    !------------------------------ Description ----------------------------
    ! <INPUT>
    ! G       - Initial guess or data from previous call
    ! Aw      - West   radiative based coefficient for matix A
    ! Ae      - East   radiative based coefficient for matix A
    ! As      - South  radiative based coefficient for matix A
    ! An      - North  radiative based coefficient for matix A
    ! Ab      - Bottom radiative based coefficient for matix A
    ! At      - Top    radiative based coefficient for matix A
    ! Ap      - Center radiative based coefficient for matix A
    ! A_RHS   - cappa term thad must be added to Ap
    ! RHS     - Right hand side for matix B
    ! a_bxG   - a boundary coefficient x-direction walls
    ! a_byG   - a boundary coefficient y-direction walls
    ! a_bzG   - a boundary coefficient z-direction walls
    ! b_bxG   - b boundary coefficient x-direction walls
    ! b_byG   - b boundary coefficient y-direction walls
    ! b_bzG   - b boundary coefficient z-direction walls
    ! g_bxG   - g boundary coefficient x-direction walls
    ! g_byG   - g boundary coefficient y-direction walls
    ! g_bzG   - g boundary coefficient z-direction walls
    ! cappa   - Absorption coefficient for a given band
    ! cappaBND- cappa for all bands
    ! sigma   - Scattering coefficient
    ! T_energy- Temperature for all nodes
    ! rad_tol - Radiation convergence criteria
    ! boltz   - Stefan boltzmann constant
    ! nsb     - number of spectral bands (nsb = 1 Gray , nsb = 5 WSGG)
    ! BBF     - blackbody radiation fraction
    ! <OUTPUT>
    ! G       - Radiant intensity (i,j,k)
    ! S_rad   - Divergent radiative flux
    ! res_rad - Residue for radiative equation
    ! ITP_G   - Number of iterations required in P1 equation
    ! <INTERNAL>
    ! Gband   - Radiant intensity for a given band
    ! Sband   - Divergent radiative flux for a given band
    ! GBND    - Radiant intensity value of the previous iteration for all bands
    ! GAMA    - Gamma coefficient for P1
    ! beta    - Extinction coefficient
    ! omega   - The over-ralaxation factor
    ! <EXTERNAL ROUTINES>
    ! elliptic_coefficients - Calculates the elliptic equation coefficients 
    ! bound_P1         - Boundary conditions for P1 
    ! boundary         - Applies boundary conditions to the coefficients
    ! SOR              - Solves the linear system via SOR
    ! wall_properties  - Calculates the value of the variables at the walls
    ! wall_fluxes_pone - Calculates all wall fluxes according to P1
    !-----------------------------------------------------------------------
    subroutine P1_nongray()
        integer :: i,j,k,IBND
        double precision :: Aw(nxi,nyi,nzi),As(nxi,nyi,nzi),Ab(nxi,nyi,nzi),&
                            Ae(nxi,nyi,nzi),An(nxi,nyi,nzi),At(nxi,nyi,nzi),&
                            Ap(nxi,nyi,nzi),RHS(nxi,nyi,nzi),A_RHS(nxi,nyi,nzi),&
                            Gband(nxi,nyi,nzi),Sband(nxi,nyi,nzi),A_dt,omega
        A_dt = 1.0d0; omega = 1.9d0        
        G = 0.0d0; S_rad = 0.0d0 !clears the values so as not to interfere with the summation
        write(*,'(a)') '=========================================='
        write(*,'(a)') '            P1 Method Residuals           '
        write(*,'(a)') '=========================================='
        BAND_LOOP: do IBND = 1,nsb
            Gband(:,:,:) = GBND(:,:,:,IBND)!Initial guess or the value of the previous iteration
            cappa = cappaBND(:,:,:,IBND)   !cappa for each band
            beta  = cappa + sigma          !beta for each band
            GAMA  = 1.0d0/(3.0d0*beta - C_rad*sigma + small)!gamma for each band
            do k = 2,nzp
                do j = 2,nyp
                    do i = 2,nxp
                        A_RHS(i,j,k) = cappa(i,j,k)
                        RHS(i,j,k)= -4.0d0*boltz*BBF(i,j,k,IBND)*cappa(i,j,k)*T_energy(i,j,k)**4.0d0 !coefficient for RHS(i,j,k)
                    end do
                end do
            end do
            call elliptic_coefficients(GAMA,A_dt,A_RHS,Ab,At,Aw,Ae,As,An,Ap)
            call bound_P1
            call boundary(Ab,At,Aw,Ae,As,An,Ap,RHS,a_bxG,b_bxG,&
                        g_bxG,a_byG,b_byG,g_byG,a_bzG,b_bzG,g_bzG,GAMA)
            !-------------------------------------------------------------------
            !Converting from meters to centimeters
            Ab = Ab*1.0d-2;At = At*1.0d-2;Ae = Ae*1.0d-2;Ap = Ap*1.0d-2
            Aw = Aw*1.0d-2;An = An*1.0d-2;As = As*1.0d-2;RHS = RHS*1.0d-6
            Gband = Gband*1.0d-4;
            call SOR(Gband,Ab,At,Aw,Ae,As,An,Ap,RHS,ITP_G,res_rad,rad_tol,omega)
            !Converting from centimeters to meters
            Gband = Gband*1.0d4;
            !-------------------------------------------------------------------
            Sband = cappa*(4.0d0*BBF(:,:,:,IBND)*boltz*T_energy**4.0d0 - Gband)
            G     = G + Gband        !Calculates the sum of all bands for G
            S_rad = S_rad + Sband    !Calculates the sum of all bands for S_rad
            write(*,'(a5,I2,a11,I5,a9,ES10.3)')'Band:',IBND,'Iteration:',ITP_G,'Residue:',res_rad
            GBND(:,:,:,IBND) = Gband !Updates the value of the current iteration
        end do BAND_LOOP
        write(*,'(a)') '------------------------------------------'
        call wall_properties(G,a_bxG,b_bxG,g_bxG,a_byG,b_byG,g_byG,a_bzG,b_bzG,g_bzG)
        call wall_fluxes_pone
    end subroutine P1_nongray
    !========================================================================================
    !-----------------------------------------------------------------------
    !                   P1 Method using the gray approach
    !------------------------------ Description ----------------------------
    ! <INPUT>
    ! G       - Initial guess or data from previous call
    ! Aw      - West   radiative based coefficient for matix A
    ! Ae      - East   radiative based coefficient for matix A
    ! As      - South  radiative based coefficient for matix A
    ! An      - North  radiative based coefficient for matix A
    ! Ab      - Bottom radiative based coefficient for matix A
    ! At      - Top    radiative based coefficient for matix A
    ! Ap      - Center radiative based coefficient for matix A
    ! A_RHS   - cappa term thad must be added to Ap
    ! RHS     - Right hand side for matix B
    ! a_bxG   - a boundary coefficient x-direction walls
    ! a_byG   - a boundary coefficient y-direction walls
    ! a_bzG   - a boundary coefficient z-direction walls
    ! b_bxG   - b boundary coefficient x-direction walls
    ! b_byG   - b boundary coefficient y-direction walls
    ! b_bzG   - b boundary coefficient z-direction walls
    ! g_bxG   - g boundary coefficient x-direction walls
    ! g_byG   - g boundary coefficient y-direction walls
    ! g_bzG   - g boundary coefficient z-direction walls
    ! cappa   - Absorption coefficient for a given band
    ! cappaBND- cappa for all bands
    ! sigma   - Scattering coefficient
    ! T_energy- Temperature for all nodes
    ! rad_tol - Radiation convergence criteria
    ! boltz   - Stefan boltzmann constant
    ! <OUTPUT>
    ! G       - Radiant intensity (i,j,k)
    ! S_rad   - Divergent radiative flux
    ! res_rad - Residue for radiative equation
    ! ITP_G   - Number of iterations required in P1 equation
    ! <INTERNAL>
    ! GAMA    - Gamma coefficient for P1
    ! beta    - Extinction coefficient
    ! omega   - The over-ralaxation factor
    ! <EXTERNAL ROUTINES>
    ! elliptic_coefficients - Calculates the elliptic equation coefficients 
    ! bound_P1         - Boundary conditions for P1 
    ! boundary         - Applies boundary conditions to the coefficients
    ! SOR              - Solves the linear system via SOR
    ! wall_properties  - Calculates the value of the variables at the walls
    ! wall_fluxes_pone - Calculates all wall fluxes according to P1
    !-----------------------------------------------------------------------
    subroutine P1_gray()
        integer :: i,j,k
        double precision :: Aw(nxi,nyi,nzi),As(nxi,nyi,nzi),Ab(nxi,nyi,nzi),&
                            Ae(nxi,nyi,nzi),An(nxi,nyi,nzi),At(nxi,nyi,nzi),&
                            Ap(nxi,nyi,nzi),RHS(nxi,nyi,nzi),A_RHS(nxi,nyi,nzi),&
                            A_dt,omega
        A_dt = 1.0d0; omega = 1.9d0
        beta  = cappa + sigma !extinction coefficient
        do k = 2,nzp
            do j = 2,nyp
                do i = 2,nxp
                    A_RHS(i,j,k) = cappa(i,j,k)
                    RHS(i,j,k)= -4.0d0*boltz*cappa(i,j,k)*T_energy(i,j,k)**4.0d0 !coefficient for RHS(i,j,k)
                end do
            end do
        end do
        GAMA = 1.0d0/(3.0d0*beta - C_rad*sigma + small)
        call elliptic_coefficients(GAMA,A_dt,A_RHS,Ab,At,Aw,Ae,As,An,Ap)
        call bound_P1
        call Boundary(Ab,At,Aw,Ae,As,An,Ap,RHS,a_bxG,b_bxG,&
                      g_bxG,a_byG,b_byG,g_byG,a_bzG,b_bzG,g_bzG,GAMA)
        !-------------------------------------------------------------------
        !Converting from meters to centimeters
        Ab = Ab*1.0d-2;At = At*1.0d-2;Ae = Ae*1.0d-2;Ap = Ap*1.0d-2
        Aw = Aw*1.0d-2;An = An*1.0d-2;As = As*1.0d-2;RHS = RHS*1.0d-6
        G = G*1.0d-4;
        call SOR(G,Ab,At,Aw,Ae,As,An,Ap,RHS,ITP_G,res_rad,rad_tol,omega)
        !Converting from centimeters to meters
        G = G*1.0d4;
        !-------------------------------------------------------------------
        S_rad = cappa*(4.0d0*boltz*T_energy**4.0d0 - G)
        call wall_properties(G,a_bxG,b_bxG,g_bxG,a_byG,b_byG,g_byG,a_bzG,b_bzG,g_bzG)
        call wall_fluxes_pone
        write(*,'(a)') '================================'
        write(*,'(a)') '       P1 Method Residuals      '
        write(*,'(a)') '================================'
        write(*,'(a10,I3,a9,ES10.3)')'Iteration:',ITP_G,'Residue:',res_rad
        write(*,'(a)') '================================'  
    end subroutine P1_gray
    !========================================================================================
    !----------------------------------------------------------
    !            Calculates all wall fluxes using P1          
    !---------------------- Description -----------------------
    ! <INPUT>
    ! T_energy- Temperature for all nodes
    ! boltz   - Stefan boltzmann constant
    ! G       - Radiant intensity (i,j,k)
    ! <OUTPUT>
    ! Q_radw  - Walls radiative fluxes
    ! <EXTERNAL ROUTINES>
    ! zetaPONE  - P1 method boundary auxiliary function
    !----------------------------------------------------------
    subroutine wall_fluxes_pone()
        integer :: i,j,k
        !---------------- x-direction walls ----------------!
        do k=2,nzp
            do j=2,nyp
                Q_radw(1,j,k)   = (4.0d0*boltz*T_energy(1,j,k)**4.0d0 - G(1,j,k))*zetaPONE(1,j,k)
                Q_radw(nxi,j,k) = (G(nxi,j,k)- 4.0d0*boltz*T_energy(nxi,j,k)**4.0d0)*zetaPONE(nxi,j,k)
            end do
        end do
        !---------------- y-direction walls ----------------!
        if(DIMEN == 2 .or. DIMEN == 3)then
            do k=2,nzp
                do i=2,nxp
                    Q_radw(i,1,k)   = (4.0d0*boltz*T_energy(i,1,k)**4.0d0 - G(i,1,k))*zetaPONE(i,1,k)
                    Q_radw(i,nyi,k) = (G(i,nyi,k) - 4.0d0*boltz*T_energy(i,nyi,k)**4.0d0)*zetaPONE(i,nyi,k)
                end do
            end do
        end if
        !---------------- y-direction walls ----------------!
        if(DIMEN == 3)then
            do j=2,nyp
                do i=2,nxp
                    Q_radw(i,j,1)   = (4.0d0*boltz*T_energy(i,j,1)**4.0d0 - G(i,j,1))*zetaPONE(i,j,1)
                    Q_radw(i,j,nzi) = (G(i,j,nzi) - 4.0d0*boltz*T_energy(i,j,nzi)**4.0d0)*zetaPONE(i,j,nzi)
                end do
            end do
        end if
    end subroutine wall_fluxes_pone
    !========================================================================================
    !                                  Finite Angle Method
    !========================================================================================
    !-----------------------------------------------------------------
    !                       Finite Angle Method
    !-------------------------- Description --------------------------
    ! <EXTERNAL ROUTINES>
    ! FAM_nongray  - Finite Angle Method using the non-gray approach
    ! FAM_gray     - Finite Angle Method using the gray approach
    !-----------------------------------------------------------------
    subroutine FAM()
        if(nongray_flag)then
            call FAM_nongray
        else
            call FAM_gray
        end if
    end subroutine FAM
    !========================================================================================
    !          Finite Angle Method using the non-gray approach
    !-----------------------------------------------------------------
    !-------------------------- Description --------------------------
    ! <INPUT>
    ! T_energy- Temperature for all nodes
    ! rad_tol - Radiation convergence criteria
    ! boltz   - Stefan boltzmann constant
    ! cappa   - Absorption coefficient for a given band
    ! cappaBND- cappa for all bands
    ! sigma   - Scattering coefficient
    ! nsb     - number of spectral bands (nsb = 1 Gray , nsb = 5 WSGG)
    ! BBF     - blackbody radiation fraction
    ! <OUTPUT>
    ! G       - Radiant intensity (i,j,k)
    ! S_rad   - Divergent radiative flux
    ! res_rad - Residue for radiative equation
    ! ITP_G   - Number of iterations required to solve the RTE
    ! <INTERNAL>
    ! IBlack  - Blackbody radiation intensity
    ! IBND    - Radiant intensity value of the previous iteration for all bands
    ! beta    - Extinction coefficient
    ! Sm      - RHS of the RTE angular term
    ! IBFF    - Blackbody radiation fraction for each band
    ! IGband  - Radiant intensity for a given band (i,j,k,l,m)
    ! Gband   - Radiant intensity for a given band
    ! Sband   - Divergent radiative flux for a given band
    ! <EXTERNAL ROUTINES>
    ! RHS_SM_FAM      - Calculates the RHS of the discretized equation
    ! agular_loop     - Performs the sweep on each octant
    ! FAMbound_in     - Inflow boundary conditions for FAM 
    ! FAMbound_out    - Outflow boundary conditions for FAM 
    ! G_FAM           - Calculates the intensity for the spatial mesh
    ! wall_fluxes_FAM - Calculates all wall fluxes according to FAM
    !-----------------------------------------------------------------
    subroutine FAM_nongray
        integer :: IBND
        double precision :: IBFF(nxi,nyi,nzi),Gband(nxi,nyi,nzi),Sband(nxi,nyi,nzi),&
                            Sm(nxi,nyi,nzi,nt,np),IGband(nxi,nyi,nzi,nt,np)
        IBlack = (boltz/PI)*T_energy**4.0d0!Computes the radiation emissions
        G = 0.0d0; S_rad = 0.0d0 !clears the values so as not to interfere with the summation
        write(*,'(a)') '=========================================='
        write(*,'(a)') '      Finite Angle Method Residuals       '
        write(*,'(a)') '=========================================='
        BAND_LOOP: do IBND = 1,nsb
            IBFF = IBlack*BBF(:,:,:,IBND) !blackbody radiation fraction for each band
            IGband = IGBND(:,:,:,:,:,IBND)!Initial guess or the value of the previous iteration
            cappa = cappaBND(:,:,:,IBND)  !cappa for each band
            beta  = cappa + sigma         !beta for each band
            do ITP_G = 1,ITMAX
                call FAMbound_in(IGband,IBFF) !Calculates the inflow boundary intensity
                call RHS_SM_FAM(Sm,IGband,IBFF)
                call agular_loop(Sm,IGband,res_rad)
                call FAMbound_out(IGband) !Computes the outflow boundary intensity
                write(*,'(a5,I2,a11,I5,a9,ES10.3)')'Band:',IBND,'Iteration:',ITP_G,'Residue:',res_rad
                if(res_rad < rad_tol) exit
            end do
            IGBND(:,:,:,:,:,IBND) = IGband!Updates the value of the current iteration
            call G_FAM(Gband,IGband)
            Sband = cappa*(4.0d0*BBF(:,:,:,IBND)*boltz*T_energy**4.0d0 - Gband)
            G     = G + Gband             !Calculates the sum of all bands for G
            S_rad = S_rad + Sband         !Calculates the sum of all bands for S_rad
            write(*,'(a)') '------------------------------------------'
        end do BAND_LOOP
        !----------------------------------------------
        !Computes the wall fluxes considering all bands
        !----------------------------------------------
        IGband = 0.0d0
        do IBND = 1,nsb
            IGband = IGband + IGBND(:,:,:,:,:,IBND)
        end do
        call wall_fluxes_FAM(IGband)
        !----------------------------------------------
    end subroutine FAM_nongray
    !========================================================================================
    !-----------------------------------------------------------------
    !           Finite Angle Method using the gray approach
    !-------------------------- Description --------------------------
    ! <INPUT>
    ! T_energy- Temperature for all nodes
    ! rad_tol - Radiation convergence criteria
    ! boltz   - Stefan boltzmann constant
    ! cappa   - Absorption coefficient for a given band
    ! sigma   - Scattering coefficient
    ! <OUTPUT>
    ! G       - Radiant intensity (i,j,k)
    ! S_rad   - Divergent radiative flux
    ! res_rad - Residue for radiative equation
    ! ITP_G   - Number of iterations required to solve the RTE
    ! <INTERNAL>
    ! IG      - Radiant intensity (i,j,k,l,m)
    ! IBlack  - Blackbody radiation intensity
    ! beta    - Extinction coefficient
    ! Sm      - RHS of the RTE angular term
    ! <EXTERNAL ROUTINES>
    ! RHS_SM_FAM      - Calculates the RHS of the discretized equation
    ! agular_loop     - Performs the sweep on each octant
    ! FAMbound_in     - Inflow boundary conditions for FAM 
    ! FAMbound_out    - Outflow boundary conditions for FAM 
    ! G_FAM           - Calculates the intensity for the spatial mesh
    ! wall_fluxes_FAM - Calculates all wall fluxes according to FAM
    !-----------------------------------------------------------------
    subroutine FAM_gray
        double precision :: Sm(nxi,nyi,nzi,nt,np)
        IBlack = (boltz/PI)*T_energy**4.0d0 !Blackbody radiation emissions
        beta  = cappa + sigma       !extinction coefficient
        write(*,'(a)') '================================'
        write(*,'(a)') '  Finite Angle Method Residuals '
        write(*,'(a)') '================================'
        do ITP_G = 1,ITMAX
            call FAMbound_in(IG,IBlack) !Computes the inflow boundary intensity
            call RHS_SM_FAM(Sm,IG,IBlack)
            call agular_loop(Sm,IG,res_rad)
            call FAMbound_out(IG)   !Computes the outflow boundary intensity
            write(*,'(a10,I3,a9,ES10.3)')'Iteration:',ITP_G,'Residue:',res_rad
            if(res_rad < rad_tol) exit
        end do
        write(*,'(a)') '================================'        
        call G_FAM(G,IG)
        S_rad = cappa*(4.0d0*boltz*T_energy**4.0d0 - G)
        call wall_fluxes_FAM(IG)        
    end subroutine FAM_gray
    !========================================================================================
    !----------------------------------------------------------
    !          Calculates the RHS of the discretized 
    !           equation by the Finite Angle Method
    !---------------------- Description -----------------------
    ! <INPUT>
    ! IG      - Radiant intensity
    ! dco     - Integral over the solid angle
    ! phase_f - Scattering phase function for FAM
    ! cappa   - Absorption coefficient 
    ! sigma   - Scattering coefficient
    ! Iblack  - Blackbody radiation intensity
    ! <OUTPUT>
    ! Sm      - RHS of the RTE angular term
    ! <INTERNAL>
    ! SMSUM   - Summation auxiliary variable
    !----------------------------------------------------------
    subroutine RHS_SM_FAM(Sm,IG,Iblack)

        double precision :: Iblack(nxi,nyi,nzi),Sm(nxi,nyi,nzi,nt,np),IG(nxi,nyi,nzi,nt,np),SMSUM
        integer :: i,j,k,l,m,ll,mm
        do k=2,nzp
            do j=2,nyp
                do i=2,nxp
                    do l=1,nt
                        do m=1,np
                            SMSUM = 0.0d0
                            if(sigma(i,j,k) /= 0.0d0)then
                                do ll=1,nt
                                    do mm=1,np
                                        SMSUM = SMSUM + IG(i,j,k,ll,mm)*phase_f(ll,mm,l,m)*dco(ll,mm)
                                    end do
                                end do
                            end if
                            Sm(i,j,k,l,m) = cappa(i,j,k)*Iblack(i,j,k) + sigma(i,j,k)*(SMSUM/PI4)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine RHS_SM_FAM
    !========================================================================================
    !----------------------------------------------------------
    !             Agular Loop of Finite Angle Method
    !---------------------- Description -----------------------
    ! <INPUT>
    ! IG      - Radiant intensity (i,j,k,l,m) initial 
    !           guess or data from previous call
    ! Sm      - RHS of the RTE angular term
    ! beta    - Extinction coefficient for all nodes
    ! volom   - Volume of the spatial-angular cell (i,j,k,l,m)
    ! Ax      - x-direction coefficients
    ! Ay      - y-direction coefficients
    ! Az      - z-direction coefficients
    ! P2      - First  azimuthal region index (0 < phi < PIBY2)
    ! P3      - Second azimuthal region index (PIBY2 < phi < PI)
    ! P4      - Third  azimuthal region index (PI < phi < PI32)
    ! T2      - First polar region index (0 < theta < PIBY2)
    ! <OUTPUT>
    ! IG      - Radiant intensity (i,j,k,l,m) solution
    ! eps     - Convergence tolerance residue 
    ! <INTERNAL>
    ! num     - Auxiliary variable
    ! deno    - Auxiliary variable
    ! dif     - Auxiliary variable
    ! Ip      - data from previous iteration
    ! <EXTERNAL ROUTINES>
    ! rad_scheme - Calculates the spatial discretizing scheme
    !----------------------------------------------------------
    subroutine agular_loop(Sm,IG,eps)
        integer :: i,j,k,l,m
        double precision :: Sm(nxi,nyi,nzi,nt,np),IG(nxi,nyi,nzi,nt,np),num,deno,dif,Ip,Iw,Is,Ib,eps
        eps = -1.0d0
        !1st octant (i=2,nxp;j=2,nyp;k=2,nzp;l=1,T2;m=1,P2)
        do k=2,nzp
            do j=2,nyp
                do i=2,nxp
                    do l=1,T2
                        do m=1,P2
                            Ip = IG(i,j,k,l,m)
                            call rad_scheme(IG(i-1,j,k,l,m),IG(i,j-1,k,l,m),IG(i,j,k-1,l,m),Ip,Iw,Is,Ib)
                            num  = Ax(j,k,l,m)*Iw + Ay(i,k,l,m)*Is + Az(i,j,l,m)*Ib + alp_r*volom(i,j,k,l,m)*Sm(i,j,k,l,m) 
                            deno = Ax(j,k,l,m) + Ay(i,k,l,m) + Az(i,j,l,m) + alp_r*beta(i,j,k)*volom(i,j,k,l,m)
                            IG(i,j,k,l,m) = num/(deno + small)
                            dif = abs(IG(i,j,k,l,m) - Ip)/(IG(i,j,k,l,m) + small)
                            eps = max(eps,dif)
                        end do
                    end do
                end do
            end do
        end do
        !2nd octant (i=nxp,2;j=2,nyp;k=2,nzp;l=1,T2;m=P2+1,P3)
        do k=2,nzp
            do j=2,nyp
                do i=nxp,2,-1
                    do l=1,T2
                        do m=P2+1,P3
                            Ip = IG(i,j,k,l,m)
                            call rad_scheme(IG(i+1,j,k,l,m),IG(i,j-1,k,l,m),IG(i,j,k-1,l,m),Ip,Iw,Is,Ib)
                            num  = Ax(j,k,l,m)*Iw + Ay(i,k,l,m)*Is + Az(i,j,l,m)*Ib + alp_r*volom(i,j,k,l,m)*Sm(i,j,k,l,m) 
                            deno = Ax(j,k,l,m) + Ay(i,k,l,m) + Az(i,j,l,m) + alp_r*beta(i,j,k)*volom(i,j,k,l,m)
                            IG(i,j,k,l,m) = num/(deno + small)
                            dif = abs(IG(i,j,k,l,m) - Ip)/(IG(i,j,k,l,m) + small)
                            eps = max(eps,dif)
                        end do
                    end do
                end do
            end do
        end do
        !3rd octant (i=nxp,2;j=nyp,2;k=2,nzp;l=1,T2;m=P3+1,P4)
        do k=2,nzp
            do j=nyp,2,-1
                do i=nxp,2,-1
                    do l=1,T2
                        do m=P3+1,P4
                            Ip = IG(i,j,k,l,m)
                            call rad_scheme(IG(i+1,j,k,l,m),IG(i,j+1,k,l,m),IG(i,j,k-1,l,m),Ip,Iw,Is,Ib)
                            num  = Ax(j,k,l,m)*Iw + Ay(i,k,l,m)*Is + Az(i,j,l,m)*Ib + alp_r*volom(i,j,k,l,m)*Sm(i,j,k,l,m) 
                            deno = Ax(j,k,l,m) + Ay(i,k,l,m) + Az(i,j,l,m) + alp_r*beta(i,j,k)*volom(i,j,k,l,m)
                            IG(i,j,k,l,m) = num/(deno + small)
                            dif = abs(IG(i,j,k,l,m) - Ip)/(IG(i,j,k,l,m) + small)
                            eps = max(eps,dif)
                        end do
                    end do
                end do
            end do
        end do
        !4th octant (i=2,nxp;j=nyp,2;k=2,nzp;l=1,T2;m=P4+1,np)
        do k=2,nzp
            do j=nyp,2,-1
                do i=2,nxp
                    do l=1,T2
                        do m=P4+1,np
                            Ip = IG(i,j,k,l,m)
                            call rad_scheme(IG(i-1,j,k,l,m),IG(i,j+1,k,l,m),IG(i,j,k-1,l,m),Ip,Iw,Is,Ib)
                            num  = Ax(j,k,l,m)*Iw + Ay(i,k,l,m)*Is + Az(i,j,l,m)*Ib + alp_r*volom(i,j,k,l,m)*Sm(i,j,k,l,m) 
                            deno = Ax(j,k,l,m) + Ay(i,k,l,m) + Az(i,j,l,m) + alp_r*beta(i,j,k)*volom(i,j,k,l,m)
                            IG(i,j,k,l,m) = num/(deno + small)
                            dif = abs(IG(i,j,k,l,m) - Ip)/(IG(i,j,k,l,m) + small)
                            eps = max(eps,dif)
                        end do
                    end do
                end do
            end do
        end do
        !5th octant (i=2,nxp;j=2,nyp;k=nzp,2;l=T2+1,nt;m=1,P2)
        do k=nzp,2,-1
            do j=2,nyp
                do i=2,nxp
                    do l=T2+1,nt
                        do m=1,P2
                            Ip = IG(i,j,k,l,m)
                            call rad_scheme(IG(i-1,j,k,l,m),IG(i,j-1,k,l,m),IG(i,j,k+1,l,m),Ip,Iw,Is,Ib)
                            num  = Ax(j,k,l,m)*Iw + Ay(i,k,l,m)*Is + Az(i,j,l,m)*Ib + alp_r*volom(i,j,k,l,m)*Sm(i,j,k,l,m) 
                            deno = Ax(j,k,l,m) + Ay(i,k,l,m) + Az(i,j,l,m) + alp_r*beta(i,j,k)*volom(i,j,k,l,m)
                            IG(i,j,k,l,m) = num/(deno + small)
                            dif = abs(IG(i,j,k,l,m) - Ip)/(IG(i,j,k,l,m) + small)
                            eps = max(eps,dif)
                        end do
                    end do
                end do
            end do
        end do
        !6th octant (i=nxp,2;j=2,nyp;k=nzp,2;l=T2+1,nt;m=P2+1,P3)
        do k=nzp,2,-1
            do j=2,nyp
                do i=nxp,2,-1
                    do l=T2+1,nt
                        do m=P2+1,P3
                            Ip = IG(i,j,k,l,m)
                            call rad_scheme(IG(i+1,j,k,l,m),IG(i,j-1,k,l,m),IG(i,j,k+1,l,m),Ip,Iw,Is,Ib)
                            num  = Ax(j,k,l,m)*Iw + Ay(i,k,l,m)*Is + Az(i,j,l,m)*Ib + alp_r*volom(i,j,k,l,m)*Sm(i,j,k,l,m) 
                            deno = Ax(j,k,l,m) + Ay(i,k,l,m) + Az(i,j,l,m) + alp_r*beta(i,j,k)*volom(i,j,k,l,m)
                            IG(i,j,k,l,m) = num/(deno + small)
                            dif = abs(IG(i,j,k,l,m) - Ip)/(IG(i,j,k,l,m) + small)
                            eps = max(eps,dif)
                        end do
                    end do
                end do
            end do
        end do
        !7th octant (i=nxp,2;j=nyp,2;k=nzp,2;l=T2+1,nt;m=P3+1,P4)
        do k=nzp,2,-1
            do j=nyp,2,-1
                do i=nxp,2,-1
                    do l=T2+1,nt
                        do m=P3+1,P4
                            Ip = IG(i,j,k,l,m)
                            call rad_scheme(IG(i+1,j,k,l,m),IG(i,j+1,k,l,m),IG(i,j,k+1,l,m),Ip,Iw,Is,Ib)
                            num  = Ax(j,k,l,m)*Iw + Ay(i,k,l,m)*Is + Az(i,j,l,m)*Ib + alp_r*volom(i,j,k,l,m)*Sm(i,j,k,l,m) 
                            deno = Ax(j,k,l,m) + Ay(i,k,l,m) + Az(i,j,l,m) + alp_r*beta(i,j,k)*volom(i,j,k,l,m)
                            IG(i,j,k,l,m) = num/(deno + small)
                            dif = abs(IG(i,j,k,l,m) - Ip)/(IG(i,j,k,l,m) + small)
                            eps = max(eps,dif)
                        end do
                    end do
                end do
            end do
        end do
        !8th octant (i=2,nxp;j=nyp,2;k=nzp,2;l=T2+1,nt;m=P4+1,np)
        do k=nzp,2,-1
            do j=nyp,2,-1
                do i=2,nxp
                    do l=T2+1,nt
                        do m=P4+1,np
                            Ip = IG(i,j,k,l,m)
                            call rad_scheme(IG(i-1,j,k,l,m),IG(i,j+1,k,l,m),IG(i,j,k+1,l,m),Ip,Iw,Is,Ib)
                            num  = Ax(j,k,l,m)*Iw + Ay(i,k,l,m)*Is + Az(i,j,l,m)*Ib + alp_r*volom(i,j,k,l,m)*Sm(i,j,k,l,m) 
                            deno = Ax(j,k,l,m) + Ay(i,k,l,m) + Az(i,j,l,m) + alp_r*beta(i,j,k)*volom(i,j,k,l,m)
                            IG(i,j,k,l,m) = num/(deno + small)
                            dif = abs(IG(i,j,k,l,m) - Ip)/(IG(i,j,k,l,m) + small)
                            eps = max(eps,dif)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine agular_loop
    !========================================================================================
    !----------------------------------------------------------
    !    Calculates the information at the center of the 
    !       control volumes according to the selected 
    !              spatial discretizing scheme
    !---------------------- Description -----------------------
    ! <INPUT>
    ! Ix - Radiant intensity in a given i position
    ! Iy - Radiant intensity in a given j position
    ! Iz - Radiant intensity in a given k position
    ! Ip - Radiant intensity in a given (i,j,k) position 
    ! <OUTPUT>
    ! Iw - West   face intensity value
    ! Is - South  face intensity value
    ! Ib - Bottom face intensity value
    !----------------------------------------------------------
    subroutine rad_scheme(Ix,Iy,Iz,Ip,Iw,Is,Ib)
        double precision :: Iw,Is,Ib,Ix,Iy,Iz,Ip
        Iw = alp_r*Ix + (1.0d0 - alp_r)*Ip
        Is = alp_r*Iy + (1.0d0 - alp_r)*Ip
        Ib = alp_r*Iz + (1.0d0 - alp_r)*Ip
    end subroutine rad_scheme
    !========================================================================================
    !----------------------------------------------------------
    !   Calculates the inflow boundary intensity keeping the 
    !      energy balance conserved between energy inflow 
    !        and outflow for the overlapped boundaries
    !---------------------- Description -----------------------
    ! <INPUT>
    ! IG          - Radiant intensity (i,j,k,l,m)
    ! Iblack      - Blackbody radiation intensity
    ! epsilon_rad - emissivities at the walls and domain(gas)
    ! P2          - First  azimuthal region index (0 < phi < PIBY2)
    ! P3          - Second azimuthal region index (PIBY2 < phi < PI)
    ! P4          - Third  azimuthal region index (PI < phi < PI32)
    ! T2          - First polar region index (0 < theta < PIBY2)
    ! SBCwall     - Symmetry boundary conditions flags
    ! <INTERNAL>
    ! GSUM        - Summation auxiliary variable
    ! bound       - Auxiliary boundary variable
    ! ia          - Auxiliary index for symmetry boundary conditions
    ! <OUTPUT>
    ! IG          - Updated boundary radiant intensity (i,j,k,l,m)
    ! <EXTERNAL ROUTINES>
    ! IGSFAM      - Calculates the summation of radiant intensity 
    !----------------------------------------------------------
    subroutine FAMbound_in(IG,Iblack)
        double precision :: Iblack(nxi,nyi,nzi),IG(nxi,nyi,nzi,nt,np),GSUM,bound
        integer :: i,j,k,l,m,ia
        !---------------------- X Boundaries ----------------------!
        do k=2,nzp
            do j=2,nyp
                !West Boundary
                if(SBCwall(1)) then
                    call IGSFAM(IG,1,j,k,-1,dcx,GSUM)
                    bound = epsilon_rad(1,j,k)*Iblack(1,j,k) + (1.0d0 - epsilon_rad(1,j,k))*(GSUM/PI)
                    do l=1,nt
                        !Interval (0 < phi < PIBY2)
                        do m=1,P2 
                            IG(1,j,k,l,m) = bound
                        end do
                        !Interval (PI32 < phi < 2PI)
                        do m=P4+1,np
                            IG(1,j,k,l,m) = bound
                        end do
                    end do
                else
                    do l=1,nt
                        !Interval (0 < phi < PIBY2)
                        do m=1,P2 
                            ia = P3 - (m - 1)
                            IG(1,j,k,l,m) = IG(1,j,k,l,ia)
                        end do
                        !Interval (PI32 < phi < 2PI)
                        do m=P4+1,np
                            ia = P4 - (m - (P4+1))
                            IG(1,j,k,l,m) = IG(1,j,k,l,ia)
                        end do
                    end do
                end if
                !East Boundary (PIBY2 < phi < PI32)
                if(SBCwall(2)) then
                    call IGSFAM(IG,nxi,j,k,1,dcx,GSUM)
                    bound = epsilon_rad(nxi,j,k)*Iblack(nxi,j,k) + (1.0d0 - epsilon_rad(nxi,j,k))*(GSUM/PI)
                    do l=1,nt
                        do m=P2+1,P4
                            IG(nxi,j,k,l,m) = bound
                        end do
                    end do
                else
                    do l=1,nt
                        do m=P2+1,P3
                            ia = P2 - (m - (P2+1))
                            IG(nxi,j,k,l,m) = IG(nxi,j,k,l,ia)
                        end do
                        do m=P3+1,P4
                            ia = np - (m - (P3+1))
                            IG(nxi,j,k,l,m) = IG(nxi,j,k,l,ia)
                        end do
                    end do
                end if
            end do
        end do
        !---------------------- Y Boundaries ----------------------!
        if(DIMEN == 2 .or. DIMEN == 3)then
            do k=2,nzp
                do i=2,nxp
                    !South Boundary (0 < phi < PI)
                    if(SBCwall(3)) then
                        call IGSFAM(IG,i,1,k,-1,dcy,GSUM)
                        bound = epsilon_rad(i,1,k)*Iblack(i,1,k) + (1.0d0 - epsilon_rad(i,1,k))*(GSUM/PI)
                        do l=1,nt
                            do m=1,P3
                                IG(i,1,k,l,m) = bound
                            end do
                        end do
                    else
                        do l=1,nt
                            do m=1,P3
                                ia = np - (m - 1)
                                IG(i,1,k,l,m) = IG(i,1,k,l,ia)
                            end do
                        end do
                    end if
                    !North Boundary (PI < phi < 2PI)
                    if(SBCwall(4)) then
                        call IGSFAM(IG,i,nyi,k,1,dcy,GSUM)
                        bound = epsilon_rad(i,nyi,k)*Iblack(i,nyi,k) + (1.0d0 - epsilon_rad(i,nyi,k))*(GSUM/PI)
                        do l=1,nt
                            do m=P3+1,np
                                IG(i,nyi,k,l,m) = bound
                            end do
                        end do
                    else
                        do l=1,nt
                            do m=P3+1,np
                                ia = np - (m - 1)
                                IG(i,nyi,k,l,m) = IG(i,nyi,k,l,ia)
                            end do
                        end do
                    end if
                end do
            end do
        end if
        !---------------------- Z Boundaries ----------------------!
        if(DIMEN == 3)then
            do j=2,nyp
                do i=2,nxp
                    !Bottom Boundary (0 < theta < PIBY2)
                    if(SBCwall(5)) then
                        call IGSFAM(IG,i,j,1,-1,dcz,GSUM)
                        bound = epsilon_rad(i,j,1)*Iblack(i,j,1) + (1.0d0 - epsilon_rad(i,j,1))*(GSUM/PI)
                        do m=1,np
                            do l=1,T2
                                IG(i,j,1,l,m) = bound
                            end do
                        end do
                    else
                        do m=1,np
                            do l=1,T2
                                ia = nt - (l - 1)
                                IG(i,j,1,l,m) = IG(i,j,1,ia,m)
                            end do                            
                        end do
                    end if
                    !Top Boundary (PIBY2 < theta < PI)
                    if(SBCwall(6)) then
                        call IGSFAM(IG,i,j,nzi,1,dcz,GSUM)
                        bound = epsilon_rad(i,j,nzi)*Iblack(i,j,nzi) + (1.0d0 - epsilon_rad(i,j,nzi))*(GSUM/PI)
                        do m=1,np
                            do l=T2+1,nt
                                IG(i,j,nzi,l,m) = bound
                            end do
                        end do
                    else
                        do m=1,np
                            do l=T2+1,nt
                                ia = nt - (l - 1)
                                IG(i,j,nzi,l,m) = IG(i,j,nzi,ia,m)
                            end do
                        end do
                    end if
                end do
            end do
        end if
    end subroutine FAMbound_in
    !========================================================================================
    !----------------------------------------------------------
    !       Calculates the summation of radiant intensity 
    !       at the boundaries, keeping the energy balance 
    !        conserved between energy inflow and outflow
    !---------------------- Description -----------------------
    ! <INPUT>
    ! i       - x-direction index
    ! j       - y-direction index
    ! k       - z-direction index
    ! IG      - Radiant intensity (i,j,k,l,m)
    ! D       - Integral over the solid angle in a given direction
    ! indx    - Inflow/outflow boundary intensity index
    ! <OUTPUT>
    ! GSUM    - Summation auxiliary variable
    !----------------------------------------------------------
    subroutine IGSFAM(IG,i,j,k,indx,D,GSUM)
        double precision :: IG(nxi,nyi,nzi,nt,np),D(nt,np),GSUM
        integer :: i,j,k,l,m,indx
        GSUM = 0.0d0
        if(indx > 0) then
            do l=1,nt
                do m=1,np
                    GSUM = GSUM + IG(i,j,k,l,m)*max(D(l,m),0.0d0) 
                end do
            end do
        else
            do l=1,nt
                do m=1,np
                    GSUM = GSUM + IG(i,j,k,l,m)*max(-D(l,m),0.0d0)
                end do
            end do
        end if
    end subroutine IGSFAM
    !========================================================================================
    !----------------------------------------------------------
    ! Update the outflow boundary intensity keeping the energy
    !    balance conserved between energy inflow and outflow
    !               for the overlapped boundaries
    !---------------------- Description -----------------------
    ! <INPUT>
    ! P2      - First  azimuthal region index (0 < phi < PIBY2)
    ! P3      - Second azimuthal region index (PIBY2 < phi < PI)
    ! P4      - Third  azimuthal region index (PI < phi < PI32)
    ! T2      - First polar region index (0 < theta < PIBY2)
    ! IG      - Radiant intensity (i,j,k,l,m)
    ! <OUTPUT>
    ! IG      - Update the radiant intensity (i,j,k,l,m)
    !----------------------------------------------------------
    subroutine FAMbound_out(IG)
        integer :: i,j,k,l,m
        double precision :: IG(nxi,nyi,nzi,nt,np)
        !---------------------- X Boundaries ----------------------!
        do k=2,nzp
            do j=2,nyp
                do l=1,nt
                    !West Boundary (PIBY2 < phi < PI32)
                    do m=P2+1,P4
                        IG(1,j,k,l,m) = IG(2,j,k,l,m)
                    end do
                    !East Boundary (0 < phi < PIBY2)
                    do m=1,P2
                        IG(nxi,j,k,l,m) = IG(nxp,j,k,l,m)
                    end do
                    !East Boundary (PI32 < phi < 2PI)
                    do m=P4+1,np
                        IG(nxi,j,k,l,m) = IG(nxp,j,k,l,m)
                    end do
                end do
            end do
        end do
        !---------------------- Y Boundaries ----------------------!
        if(DIMEN == 2 .or. DIMEN == 3)then
            do k=2,nzp
                do i=2,nxp
                    do l=1,nt
                        !South Boundary (PI < phi < 2PI)
                        do m=P3+1,np
                            IG(i,1,k,l,m) = IG(i,2,k,l,m)
                        end do
                        !North Boundary (0 < phi < PI)
                        do m=1,P3
                            IG(i,nyi,k,l,m) = IG(i,nyp,k,l,m)
                        end do
                    end do
                end do
            end do
        end if
        !---------------------- Z Boundaries ----------------------!
        if(DIMEN == 3)then
            do j=2,nyp
                do i=2,nxp
                    do m=1,np
                        !Bottom Boundary (PIBY2 < theta < PI)
                        do l=T2+1,nt
                            IG(i,j,1,l,m) = IG(i,j,2,l,m)
                        end do
                        !Top Boundary (0 < theta < PIBY2)
                        do l=1,T2
                            IG(i,j,nzi,l,m) = IG(i,j,nzp,l,m)
                        end do
                    end do
                end do
            end do
        end if
    end subroutine FAMbound_out
    !========================================================================================
    !----------------------------------------------------------
    !          Calculates the Radiant intensity 
    !        at the spatial mesh (i,j,k) for FAM
    !---------------------- Description -----------------------
    ! <INPUT>
    ! IG      - Radiant intensity (i,j,k,l,m)
    ! dco     - Integral over the solid angle
    ! <OUTPUT>
    ! Grad    - Radiant intensity (i,j,k)
    ! <INTERNAL>
    ! GSUM    - Summation auxiliary variable
    !----------------------------------------------------------
    subroutine G_FAM(Grad,IG)
        integer :: i,j,k,l,m
        double precision :: GSUM,Grad(nxi,nyi,nzi),IG(nxi,nyi,nzi,nt,np)
        do k=1,nzi
            do i=1,nxi
                do j=1,nyi
                    GSUM = 0.0d0
                    do l=1,nt
                        do m=1,np
                            GSUM = GSUM + IG(i,j,k,l,m)*dco(l,m)
                        end do
                    end do
                    Grad(i,j,k) = GSUM
                end do
            end do
        end do
    end subroutine G_FAM
    !========================================================================================
    !----------------------------------------------------------
    !             Calculates all wall fluxes using 
    !                   Finite Angle Method 
    !---------------------- Description -----------------------
    ! <INPUT>
    ! IG      - Radiant intensity (i,j,k,l,m)
    ! dcx     - Integral over the solid angle in the x-direction
    ! dcy     - Integral over the solid angle in the y-direction
    ! dcz     - Integral over the solid angle in the z-direction
    ! <OUTPUT>
    ! Q_radw  - Walls radiative fluxes
    ! <EXTERNAL ROUTINES>
    ! IGFLUX  - Calculates the summation of radiant intensity 
    !----------------------------------------------------------
    subroutine wall_fluxes_FAM(IG)
        integer :: i,j,k
        double precision :: IG(nxi,nyi,nzi,nt,np)
        !---------------- x-direction walls ----------------!
        do k=2,nzp
            do j=2,nyp
                Q_radw(1,j,k)   = IGFLUX(1,j,k,dcx,IG)   !West Wall
                Q_radw(nxi,j,k) = IGFLUX(nxi,j,k,dcx,IG) !East Wall
            end do
        end do
        !---------------- y-direction walls ----------------!
        if(DIMEN == 2 .or. DIMEN == 3)then
            do k=2,nzp
                do i=2,nxp
                    Q_radw(i,1,k)   = IGFLUX(i,1,k,dcy,IG)   !South Wall
                    Q_radw(i,nyi,k) = IGFLUX(i,nyi,k,dcy,IG) !North Wall
                end do 
            end do
        end if
        !---------------- z-direction walls ----------------!
        if(DIMEN == 3)then
            do j=2,nyp
                do i=2,nxp
                    Q_radw(i,j,1)   = IGFLUX(i,j,1,dcz,IG)   !Bottom Wall
                    Q_radw(i,j,nzi) = IGFLUX(i,j,nzi,dcz,IG) !Top Wall
                end do
            end do
        end if
    end subroutine wall_fluxes_FAM
    !========================================================================================
    !----------------------------------------------------------
    !       Calculates the summation of radiant intensity 
    !                 balance at the boundaries
    !------------------------ Description -------------------------
    ! <INPUT>
    ! i       - x-direction index
    ! j       - y-direction index
    ! k       - z-direction index
    ! IG      - Radiant intensity (i,j,k,l,m)
    ! D       - Integral over the solid angle in a given direction
    ! IGIN    - Inflow  intensity summation 
    ! IGOUT   - Outflow intensity summation 
    ! <OUTPUT>
    ! IGFLUX  - Intensity balance between energy inflow and outflow
    !--------------------------------------------------------------
    function IGFLUX(i,j,k,D,IG)
        implicit none
        integer :: i,j,k,l,m
        double precision :: IGFLUX,IGOUT,IGIN,D(nt,np),IG(nxi,nyi,nzi,nt,np)
        IGIN=0.0d0;IGOUT=0.0d0
        do l=1,nt
            do m=1,np
                if(D(l,m) > 0.0d0)then
                    IGIN = IGIN  + abs(D(l,m))*IG(i,j,k,l,m)
                else
                    IGOUT= IGOUT + abs(D(l,m))*IG(i,j,k,l,m)
                end if
            end do
        end do
        IGFLUX = IGIN - IGOUT
    end function IGFLUX
    !========================================================================================
    !                               Discrete Ordinates Method                               
    !========================================================================================
    !-----------------------------------------------------------------
    !                    Discrete Ordinates Method
    !-------------------------- Description --------------------------
    ! <EXTERNAL ROUTINES>
    ! DOM_nongray  - Discrete Ordinates Method using the non-gray approach
    ! DOM_gray     - Discrete Ordinates Method using the gray approach
    !-----------------------------------------------------------------
    subroutine DOM()
        if(nongray_flag)then
            call DOM_nongray
        else
            call DOM_gray
        end if
    end subroutine DOM
    !========================================================================================
    !-----------------------------------------------------------------
    !      Discrete Ordinates Method using the non-gray approach
    !-------------------------- Description --------------------------
    ! <INPUT>
    ! T_energy- Temperature for all nodes
    ! rad_tol - Radiation convergence criteria
    ! boltz   - Stefan boltzmann constant
    ! cappa   - Absorption coefficient for a given band
    ! cappaBND- cappa for all bands
    ! sigma   - Scattering coefficient
    ! nsb     - number of spectral bands (nsb = 1 Gray , nsb = 5 WSGG)
    ! BBF     - blackbody radiation fraction
    ! <OUTPUT>
    ! G       - Radiant intensity (i,j,k)
    ! S_rad   - Divergent radiative flux
    ! res_rad - Residue for radiative equation
    ! ITP_G   - Number of iterations required to solve the RTE
    ! <INTERNAL>
    ! IBlack  - Blackbody radiation intensity
    ! IBND    - Radiant intensity value of the previous iteration for all bands
    ! beta    - Extinction coefficient
    ! Sm      - RHS of the RTE angular term
    ! IBFF    - Blackbody radiation fraction for each band
    ! IGband  - Radiant intensity for a given band (i,j,k,l,m)
    ! Gband   - Radiant intensity for a given band
    ! Sband   - Divergent radiative flux for a given band
    ! <EXTERNAL ROUTINES>
    ! RHS_SM_DOM      - Calculates the RHS of the discretized equation
    ! orthogonal_loop - Performs the sweep on each octant
    ! DOMbound_in     - Inflow boundary conditions for DOM
    ! DOMbound_out    - Outflow boundary conditions for DOM
    ! G_DOM           - Calculates the intensity for the spatial mesh
    ! wall_fluxes_DOM - Calculates all wall fluxes according to DOM
    !-----------------------------------------------------------------
    subroutine DOM_nongray
        integer :: IBND
        double precision :: IBFF(nxi,nyi,nzi),Gband(nxi,nyi,nzi),Sband(nxi,nyi,nzi),&
                            Sm(nxi,nyi,nzi,nq,8),IGband(nxi,nyi,nzi,nq,8)
        IBlack = (boltz/PI)*T_energy**4.0d0!Computes the radiation emissions
        G = 0.0d0; S_rad = 0.0d0 !clears the values so as not to interfere with the summation
        write(*,'(a)') '=========================================='
        write(*,'(a)') '    Discrete Ordinates Method Residuals   '
        write(*,'(a)') '=========================================='
        BAND_LOOP: do IBND = 1,nsb
            IBFF = IBlack*BBF(:,:,:,IBND) !blackbody radiation fraction for each band
            IGband = IGBND(:,:,:,:,:,IBND)!Initial guess or the value of the previous iteration
            cappa = cappaBND(:,:,:,IBND)  !cappa for each band
            beta  = cappa + sigma         !beta for each band
            do ITP_G = 1,ITMAX
                call DOMbound_in(IGband,IBFF) !Calculates the inflow boundary intensity
                call RHS_SM_DOM(Sm,IGband,IBFF)
                call orthogonal_loop(Sm,IGband,res_rad)
                call DOMbound_out(IGband)
                write(*,'(a5,I2,a11,I5,a9,ES10.3)')'Band:',IBND,'Iteration:',ITP_G,'Residue:',res_rad
                if(res_rad < rad_tol) exit
            end do
            IGBND(:,:,:,:,:,IBND) = IGband!Updates the value of the current iteration
            call G_DOM(Gband,IGband)
            Sband = cappa*(4.0d0*BBF(:,:,:,IBND)*boltz*T_energy**4.0d0 - Gband)
            G     = G + Gband             !Calculates the sum of all bands for G
            S_rad = S_rad + Sband         !Calculates the sum of all bands for S_rad
            write(*,'(a)') '------------------------------------------'
        end do BAND_LOOP
        !----------------------------------------------
        !Computes the wall fluxes considering all bands
        !----------------------------------------------
        IGband = 0.0d0
        do IBND = 1,nsb
            IGband = IGband + IGBND(:,:,:,:,:,IBND)
        end do
        call wall_fluxes_DOM(IGband)
        !----------------------------------------------
    end subroutine DOM_nongray
    !========================================================================================
    !-----------------------------------------------------------------
    !        Discrete Ordinates Method using the gray approach
    !-------------------------- Description --------------------------
    ! <INPUT>
    ! T_energy- Temperature for all nodes
    ! rad_tol - Radiation convergence criteria
    ! boltz   - Stefan boltzmann constant
    ! cappa   - Absorption coefficient for a given band
    ! sigma   - Scattering coefficient
    ! <OUTPUT>
    ! G       - Radiant intensity (i,j,k)
    ! S_rad   - Divergent radiative flux
    ! res_rad - Residue for radiative equation
    ! ITP_G   - Number of iterations required to solve the RTE
    ! <INTERNAL>
    ! IG      - Radiant intensity (i,j,k,l,m)
    ! IBlack  - Blackbody radiation intensity
    ! beta    - Extinction coefficient
    ! Sm      - RHS of the RTE angular term
    ! <EXTERNAL ROUTINES>
    ! RHS_SM_DOM      - Calculates the RHS of the discretized equation
    ! orthogonal_loop - Performs the sweep on each octant
    ! DOMbound_in     - Inflow boundary conditions for DOM
    ! DOMbound_out    - Outflow boundary conditions for DOM
    ! G_DOM           - Calculates the intensity for the spatial mesh
    ! wall_fluxes_DOM - Calculates all wall fluxes according to DOM
    !-----------------------------------------------------------------
    subroutine DOM_gray
        double precision :: Sm(nxi,nyi,nzi,nq,8)
        IBlack = (boltz/PI)*T_energy**4.0d0 !Blackbody radiation emissions
        beta  = cappa + sigma               !extinction coefficient
        write(*,'(a)') '====================================='
        write(*,'(a)') ' Discrete Ordinates Method Residuals '
        write(*,'(a)') '====================================='
        do ITP_G = 1,ITMAX
            call DOMbound_in(IG,IBlack)
            call RHS_SM_DOM(Sm,IG,IBlack)
            call orthogonal_loop(Sm,IG,res_rad)
            call DOMbound_out(IG)
            write(*,'(a10,I3,a9,ES10.3)')'Iteration:',ITP_G,'Residue:',res_rad
            if(res_rad < rad_tol) exit
        end do
        write(*,'(a)') '====================================='    
        call G_DOM(G,IG)
        S_rad = cappa*(4.0d0*boltz*T_energy**4.0d0 - G)
        call wall_fluxes_DOM(IG)
    end subroutine DOM_gray
    !========================================================================================
    !----------------------------------------------------------
    !          Calculates the RHS of the discretized 
    !        equation by the Discrete Ordinates Method
    !---------------------- Description -----------------------
    ! <INPUT>
    ! IG      - Radiant intensity
    ! Wq      - Sn quadrature weigth quadrature
    ! phase_d - Scattering phase function for DOM
    ! cappa   - Absorption coefficient 
    ! sigma   - Scattering coefficient
    ! Iblack  - Blackbody radiation intensity
    ! aux     - Auxiliary variable for emission term
    ! <OUTPUT>
    ! Sm      - RHS of the RTE angular term
    ! <INTERNAL>
    ! SMSUM   - Summation auxiliary variable
    !----------------------------------------------------------
    subroutine RHS_SM_DOM(Sm,IG,Iblack)
        double precision :: Iblack(nxi,nyi,nzi),Sm(nxi,nyi,nzi,nq,8),IG(nxi,nyi,nzi,nq,8),SMSUM,aux
        integer :: i,j,k,l,ls,m
        do k=2,nzp
            do j=2,nyp
                do i=2,nxp
                    aux = cappa(i,j,k)*Iblack(i,j,k)
                    do m=1,8
                        do l=1,nq
                            SMSUM = 0.0d0
                            if(sigma(i,j,k) /= 0.0d0)then
                                do ls=1,nq
                                    SMSUM = SMSUM + IG(i,j,k,ls,m)*phase_d(ls,l)*Wq(ls)
                                end do
                            end if
                            Sm(i,j,k,l,m) = aux + sigma(i,j,k)*(SMSUM/PI4)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine RHS_SM_DOM
    !========================================================================================
    !----------------------------------------------------------
    !      Orthogonal Loop of Discrete Ordinates Method
    !---------------------- Description -----------------------
    ! <INPUT>
    ! IG      - Radiant intensity (i,j,k,l,m) initial 
    !           guess or data from previous call
    ! Sm      - RHS of the RTE angular term
    ! beta    - Extinction coefficient for all nodes
    ! vol     - Volume of the spatial cell (i,j,k)
    ! Axd     - x-direction coefficients
    ! Ayd     - y-direction coefficients
    ! Azd     - z-direction coefficients
    ! <OUTPUT>
    ! IG      - Radiant intensity (i,j,k,l,m) solution
    ! eps     - Convergence tolerance residue 
    ! <INTERNAL>
    ! num     - Auxiliary variable
    ! deno    - Auxiliary variable
    ! dif     - Auxiliary variable
    ! Ip      - data from previous iteration
    !----------------------------------------------------------
    subroutine orthogonal_loop(Sm,IG,eps)
        integer :: i,j,k,l,m
        double precision :: Sm(nxi,nyi,nzi,nq,8),IG(nxi,nyi,nzi,nq,8),num,deno,dif,Ip,Iw,Is,Ib,eps
        eps = -1.0d0
        !1st octant (i=2,nxp;j=2,nyp;k=2,nzp;l=1,nq;m=1)
        m = 1
        do k=2,nzp
            do j=2,nyp
                do i=2,nxp
                    do l=1,nq
                        Ip = IG(i,j,k,l,m)
                        call rad_scheme(IG(i-1,j,k,l,m),IG(i,j-1,k,l,m),IG(i,j,k-1,l,m),Ip,Iw,Is,Ib)
                        num  = Axd(j,k,l)*Iw + Ayd(i,k,l)*Is + Azd(i,j,l)*Ib + alp_r*vol(i,j,k)*Sm(i,j,k,l,m)
                        deno = Axd(j,k,l) + Ayd(i,k,l) + Azd(i,j,l) + alp_r*beta(i,j,k)*vol(i,j,k)
                        IG(i,j,k,l,m) = num/(deno + small)
                        dif = abs(IG(i,j,k,l,m) - Ip)/(IG(i,j,k,l,m) + small)
                        eps = max(eps,dif)
                    end do
                end do
            end do
        end do
        !2nd octant (i=nxp,2;j=2,nyp;k=2,nzp;l=1,nq;m=2)
        m = 2
        do k=2,nzp
            do j=2,nyp
                do i=nxp,2,-1
                    do l=1,nq
                        Ip = IG(i,j,k,l,m)
                        call rad_scheme(IG(i+1,j,k,l,m),IG(i,j-1,k,l,m),IG(i,j,k-1,l,m),Ip,Iw,Is,Ib)
                        num  = Axd(j,k,l)*Iw + Ayd(i,k,l)*Is + Azd(i,j,l)*Ib + alp_r*vol(i,j,k)*Sm(i,j,k,l,m)
                        deno = Axd(j,k,l) + Ayd(i,k,l) + Azd(i,j,l) + alp_r*beta(i,j,k)*vol(i,j,k)
                        IG(i,j,k,l,m) = num/(deno + small)
                        dif = abs(IG(i,j,k,l,m) - Ip)/(IG(i,j,k,l,m) + small)
                        eps = max(eps,dif)
                    end do
                end do
            end do
        end do
        if(DIMEN == 2 .or. DIMEN == 3)then
            !3rd octant (i=nxp,2;j=nyp,2;k=2,nzp;l=1,nq;m=3)
            m = 3
            do k=2,nzp
                do j=nyp,2,-1
                    do i=nxp,2,-1
                        do l=1,nq
                            Ip = IG(i,j,k,l,m)
                            call rad_scheme(IG(i+1,j,k,l,m),IG(i,j+1,k,l,m),IG(i,j,k-1,l,m),Ip,Iw,Is,Ib)
                            num  = Axd(j,k,l)*Iw + Ayd(i,k,l)*Is + Azd(i,j,l)*Ib + alp_r*vol(i,j,k)*Sm(i,j,k,l,m)
                            deno = Axd(j,k,l) + Ayd(i,k,l) + Azd(i,j,l) + alp_r*beta(i,j,k)*vol(i,j,k)
                            IG(i,j,k,l,m) = num/(deno + small)
                            dif = abs(IG(i,j,k,l,m) - Ip)/(IG(i,j,k,l,m) + small)
                            eps = max(eps,dif)
                        end do
                    end do
                end do
            end do
            !4th octant (i=2,nxp;j=nyp,2;k=2,nzp;l=1,nq;m=4)
            m = 4
            do k=2,nzp
                do j=nyp,2,-1
                    do i=2,nxp
                        do l=1,nq
                            Ip = IG(i,j,k,l,m)
                            call rad_scheme(IG(i-1,j,k,l,m),IG(i,j+1,k,l,m),IG(i,j,k-1,l,m),Ip,Iw,Is,Ib)
                            num  = Axd(j,k,l)*Iw + Ayd(i,k,l)*Is + Azd(i,j,l)*Ib + alp_r*vol(i,j,k)*Sm(i,j,k,l,m)
                            deno = Axd(j,k,l) + Ayd(i,k,l) + Azd(i,j,l) + alp_r*beta(i,j,k)*vol(i,j,k)
                            IG(i,j,k,l,m) = num/(deno + small)
                            dif = abs(IG(i,j,k,l,m) - Ip)/(IG(i,j,k,l,m) + small)
                            eps = max(eps,dif)
                        end do
                    end do
                end do
            end do
        end if
        if(DIMEN == 3)then
            !5th octant (i=2,nxp;j=2,nyp;k=nzp,2;l=1,nq;m=5)
            m = 5
            do k=nzp,2,-1
                do j=2,nyp
                    do i=2,nxp
                        do l=1,nq
                            Ip = IG(i,j,k,l,m)
                            call rad_scheme(IG(i-1,j,k,l,m),IG(i,j-1,k,l,m),IG(i,j,k+1,l,m),Ip,Iw,Is,Ib)
                            num  = Axd(j,k,l)*Iw + Ayd(i,k,l)*Is + Azd(i,j,l)*Ib + alp_r*vol(i,j,k)*Sm(i,j,k,l,m)
                            deno = Axd(j,k,l) + Ayd(i,k,l) + Azd(i,j,l) + alp_r*beta(i,j,k)*vol(i,j,k)
                            IG(i,j,k,l,m) = num/(deno + small)
                            dif = abs(IG(i,j,k,l,m) - Ip)/(IG(i,j,k,l,m) + small)
                            eps = max(eps,dif)
                        end do
                    end do
                end do
            end do
            !6th octant (i=nxp,2;j=2,nyp;k=nzp,2;l=1,nq;m=6)
            m = 6
            do k=nzp,2,-1
                do j=2,nyp
                    do i=nxp,2,-1
                        do l=1,nq
                            Ip = IG(i,j,k,l,m)
                            call rad_scheme(IG(i+1,j,k,l,m),IG(i,j-1,k,l,m),IG(i,j,k+1,l,m),Ip,Iw,Is,Ib)
                            num  = Axd(j,k,l)*Iw + Ayd(i,k,l)*Is + Azd(i,j,l)*Ib + alp_r*vol(i,j,k)*Sm(i,j,k,l,m)
                            deno = Axd(j,k,l) + Ayd(i,k,l) + Azd(i,j,l) + alp_r*beta(i,j,k)*vol(i,j,k)
                            IG(i,j,k,l,m) = num/(deno + small)
                            dif = abs(IG(i,j,k,l,m) - Ip)/(IG(i,j,k,l,m) + small)
                            eps = max(eps,dif)
                        end do
                    end do
                end do
            end do
            !7th octant (i=nxp,2;j=nyp,2;k=nzp,2;l=1,nq;m=7)
            m = 7
            do k=nzp,2,-1
                do j=nyp,2,-1
                    do i=nxp,2,-1
                        do l=1,nq
                            Ip = IG(i,j,k,l,m)
                            call rad_scheme(IG(i+1,j,k,l,m),IG(i,j+1,k,l,m),IG(i,j,k+1,l,m),Ip,Iw,Is,Ib)
                            num  = Axd(j,k,l)*Iw + Ayd(i,k,l)*Is + Azd(i,j,l)*Ib + alp_r*vol(i,j,k)*Sm(i,j,k,l,m)
                            deno = Axd(j,k,l) + Ayd(i,k,l) + Azd(i,j,l) + alp_r*beta(i,j,k)*vol(i,j,k)
                            IG(i,j,k,l,m) = num/(deno + small)
                            dif = abs(IG(i,j,k,l,m) - Ip)/(IG(i,j,k,l,m) + small)
                            eps = max(eps,dif)
                        end do
                    end do
                end do
            end do
            !8th octant (i=2,nxp;j=nyp,2;k=nzp,2;l=1,nq;m=8)
            m = 8
            do k=nzp,2,-1
                do j=nyp,2,-1
                    do i=2,nxp
                        do l=1,nq
                            Ip = IG(i,j,k,l,m)
                            call rad_scheme(IG(i-1,j,k,l,m),IG(i,j+1,k,l,m),IG(i,j,k+1,l,m),Ip,Iw,Is,Ib)
                            num  = Axd(j,k,l)*Iw + Ayd(i,k,l)*Is + Azd(i,j,l)*Ib + alp_r*vol(i,j,k)*Sm(i,j,k,l,m)
                            deno = Axd(j,k,l) + Ayd(i,k,l) + Azd(i,j,l) + alp_r*beta(i,j,k)*vol(i,j,k)
                            IG(i,j,k,l,m) = num/(deno + small)
                            dif = abs(IG(i,j,k,l,m) - Ip)/(IG(i,j,k,l,m) + small)
                            eps = max(eps,dif)
                        end do
                    end do
                end do
            end do
        end if
    end subroutine orthogonal_loop
    !========================================================================================
    !----------------------------------------------------------
    !   Calculates the inflow boundary intensity keeping the 
    !      energy balance conserved between energy inflow 
    !             and outflow accounting reflection
    !---------------------- Description -----------------------
    ! <INPUT>
    ! IG          - Radiant intensity (i,j,k,l,m)
    ! Iblack      - Blackbody radiation intensity
    ! epsilon_rad - emissivities at the walls and domain(gas)
    ! octants_in  - inflow  boundary octants index for account reflection
    ! octants_out - outflow boundary octants index for account reflection
    ! SBCwall     - Symmetry boundary conditions flags
    ! <INTERNAL>
    ! GSUM        - Summation auxiliary variable
    ! bound       - Auxiliary boundary variable
    ! ia          - Auxiliary index for symmetry boundary conditions
    ! <OUTPUT>
    ! IG          - Updated boundary radiant intensity (i,j,k,l,m)
    ! <EXTERNAL ROUTINES>
    ! IGSDOM      - Calculates the summation of radiant intensity 
    !----------------------------------------------------------
    subroutine DOMbound_in(IG,Iblack)
        double precision :: Iblack(nxi,nyi,nzi),IG(nxi,nyi,nzi,nq,8),GSUM,bound
        integer, dimension(4)  :: octants_in,octants_out
        integer :: i,j,k,l,m,o,ia
        !---------------------- X Boundaries ----------------------!
        !West Boundary
        octants_in = (/1,4,5,8/); octants_out = (/2,3,6,7/)
        if(SBCwall(1)) then
            do k=2,nzp
                do j=2,nyp
                    call IGSDOM(IG,1,j,k,octants_out,mux,GSUM)
                    bound = epsilon_rad(1,j,k)*Iblack(1,j,k) + (1.0d0 - epsilon_rad(1,j,k))*(GSUM/PI)
                    do l=1,nq
                        do o=1,4
                            m = octants_in(o)
                            IG(1,j,k,l,m) = bound
                        end do
                    end do
                end do
            end do
        else
            do k=2,nzp
                do j=2,nyp
                    do l=1,nq
                        do o=1,4
                            m = octants_in(o)
                            ia=octants_out(o)
                            IG(1,j,k,l,m) = IG(1,j,k,l,ia)
                        end do
                    end do
                end do
            end do
        end if
        !East Boundary
        octants_in = (/2,3,6,7/); octants_out = (/1,4,5,8/)
        if(SBCwall(2)) then
            do k=2,nzp
                do j=2,nyp
                    call IGSDOM(IG,nxi,j,k,octants_out,mux,GSUM)
                    bound = epsilon_rad(nxi,j,k)*Iblack(nxi,j,k) + (1.0d0 - epsilon_rad(nxi,j,k))*(GSUM/PI)
                    do l=1,nq
                        do o=1,4
                            m = octants_in(o)
                            IG(nxi,j,k,l,m) = bound
                        end do
                    end do
                end do
            end do
        else
            do k=2,nzp
                do j=2,nyp
                    do l=1,nq
                        do o=1,4
                            m = octants_in(o)
                            ia=octants_out(o)
                            IG(nxi,j,k,l,m) = IG(nxi,j,k,l,ia)
                        end do
                    end do
                end do
            end do
        end if
        !---------------------- Y Boundaries ----------------------!
        if(DIMEN == 2 .or. DIMEN == 3)then
            !South Boundary
            octants_in = (/1,2,5,6/); octants_out = (/3,4,7,8/)
            if(SBCwall(3)) then
                do k=2,nzp
                    do i=2,nxp
                        call IGSDOM(IG,i,1,k,octants_out,etay,GSUM)
                        bound = epsilon_rad(i,1,k)*Iblack(i,1,k) + (1.0d0 - epsilon_rad(i,1,k))*(GSUM/PI)
                        do l=1,nq
                            do o=1,4
                                m = octants_in(o)
                                IG(i,1,k,l,m) = bound
                            end do
                        end do
                    end do
                end do
            else
                do k=2,nzp
                    do i=2,nxp
                        do l=1,nq
                            do o=1,4
                                m = octants_in(o)
                                ia=octants_out(o)
                                IG(i,1,k,l,m) = IG(i,1,k,l,ia)
                            end do
                        end do
                    end do
                end do
            end if
            !North Boundary
            octants_in = (/3,4,7,8/); octants_out = (/1,2,5,6/)
            if(SBCwall(4)) then
                do k=2,nzp
                    do i=2,nxp
                        call IGSDOM(IG,i,nyi,k,octants_out,etay,GSUM)
                        bound = epsilon_rad(i,nyi,k)*Iblack(i,nyi,k) + (1.0d0 - epsilon_rad(i,nyi,k))*(GSUM/PI)
                        do l=1,nq
                            do o=1,4
                                m = octants_in(o)
                                IG(i,nyi,k,l,m) = bound
                            end do
                        end do
                    end do
                end do
            else
                do k=2,nzp
                    do i=2,nxp
                        do l=1,nq
                            do o=1,4
                                m = octants_in(o)
                                ia=octants_out(o)
                                IG(i,nyi,k,l,m) = IG(i,nyi,k,l,ia)
                            end do
                        end do
                    end do
                end do
            end if
        end if
        !---------------------- Z Boundaries ----------------------!
        if(DIMEN == 3)then
            !Bottom Boundary
            octants_in = (/1,2,3,4/); octants_out = (/5,6,7,8/)
            if(SBCwall(5)) then
                do j=2,nyp
                    do i=2,nxp
                        call IGSDOM(IG,i,j,1,octants_out,xiz,GSUM)
                        bound = epsilon_rad(i,j,1)*Iblack(i,j,1) + (1.0d0 - epsilon_rad(i,j,1))*(GSUM/PI)
                        do l=1,nq
                            do o=1,4
                                m = octants_in(o)
                                IG(i,j,1,l,m) = bound
                            end do
                        end do
                    end do
                end do
            else
                do j=2,nyp
                    do i=2,nxp
                        do l=1,nq
                            do o=1,4
                                m = octants_in(o)
                                ia=octants_out(o)
                                IG(i,j,1,l,m) = IG(i,j,1,l,ia)
                            end do
                        end do
                    end do
                end do
            end if
            !Top Boundary
            octants_in = (/5,6,7,8/); octants_out = (/1,2,3,4/)
            if(SBCwall(6)) then
                do j=2,nyp
                    do i=2,nxp
                        call IGSDOM(IG,i,j,nzi,octants_out,xiz,GSUM)
                        bound = epsilon_rad(i,j,nzi)*Iblack(i,j,nzi) + (1.0d0 - epsilon_rad(i,j,nzi))*(GSUM/PI)
                        do l=1,nq
                            do o=1,4
                                m = octants_in(o)
                                IG(i,j,nzi,l,m) = bound
                            end do
                        end do
                    end do
                end do
            else
                do j=2,nyp
                    do i=2,nxp
                        do l=1,nq
                            do o=1,4
                                m = octants_in(o)
                                ia=octants_out(o)
                                IG(i,j,nzi,l,m) = IG(i,j,nzi,l,ia)
                            end do
                        end do
                    end do
                end do
            end if
        end if
    end subroutine DOMbound_in
    !========================================================================================
    !----------------------------------------------------------
    !       Calculates the summation of radiant intensity 
    !       at the boundaries, accounting for reflection
    !---------------------- Description -----------------------
    ! <INPUT>
    ! i       - x-direction index
    ! j       - y-direction index
    ! k       - z-direction index
    ! IG      - Radiant intensity (i,j,k,l,m)
    ! D       - Solid angle ordinate in a given direction
    ! indx    - outflow boundary octants index for account reflection
    ! <OUTPUT>
    ! GSUM    - Summation auxiliary variable
    !----------------------------------------------------------      
    subroutine IGSDOM(IG,i,j,k,indx,D,GSUM)  
        double precision :: IG(nxi,nyi,nzi,nq,8),D(nq),GSUM
        integer :: i,j,k,l,m,o,indx(4)
        GSUM = 0.0d0
        do l=1,nq
            do o =1,4
                m = indx(o)
                GSUM = GSUM + Wq(l)*abs(D(l))*IG(i,j,k,l,m)
            end do
        end do        
    end subroutine IGSDOM
    !========================================================================================
    !----------------------------------------------------------
    ! Update the outflow boundary intensity keeping the energy
    !    balance conserved between energy inflow and outflow
    !               for accounting reflection
    !---------------------- Description -----------------------
    ! <INPUT>
    ! IG          - Radiant intensity (i,j,k,l,m)
    ! octants_out - outflow boundary octants index 
    ! <OUTPUT>
    ! IG          - Update the radiant intensity (i,j,k,l,m)
    !----------------------------------------------------------
    subroutine DOMbound_out(IG)
        double precision :: IG(nxi,nyi,nzi,nq,8)
        integer :: i,j,k,l,m,o,octants_out(4)
        !---------------------- X Boundaries ----------------------!
        !West Boundary
        octants_out = (/2,3,6,7/)
        do k=2,nzp
            do j=2,nyp
                do l=1,nq
                    do o=1,4
                        m = octants_out(o)
                        IG(1,j,k,l,m) = IG(2,j,k,l,m)
                    end do
                end do
            end do
        end do
        !East Boundary
        octants_out = (/1,4,5,8/)
        do k=2,nzp
            do j=2,nyp
                do l=1,nq
                    do o=1,4
                        m = octants_out(o)
                        IG(nxi,j,k,l,m) = IG(nxp,j,k,l,m)
                    end do
                end do
            end do
        end do
        !---------------------- Y Boundaries ----------------------!
        if(DIMEN == 2 .or. DIMEN == 3)then
            !South Boundary
            octants_out = (/3,4,7,8/)
            do k=2,nzp
                do i=2,nxp
                    do l=1,nq
                        do o=1,4
                            m = octants_out(o)
                            IG(i,1,k,l,m) = IG(i,2,k,l,m)
                        end do
                    end do
                end do
            end do
            !North Boundary
            octants_out = (/1,2,5,6/)
            do k=2,nzp
                do i=2,nxp
                    do l=1,nq
                        do o=1,4
                            m = octants_out(o)
                            IG(i,nyi,k,l,m) = IG(i,nyp,k,l,m)
                        end do
                    end do
                end do
            end do
        end if
        !---------------------- Z Boundaries ----------------------!
        if(DIMEN == 3)then
            !Bottom Boundary
            octants_out = (/5,6,7,8/)
            do j=2,nyp
                do i=2,nxp
                    do l=1,nq
                        do o=1,4
                            m = octants_out(o)
                            IG(i,j,1,l,m) = IG(i,j,2,l,m)
                        end do
                    end do
                end do
            end do
            !Top Boundary
            octants_out = (/1,2,3,4/)
            do j=2,nyp
                do i=2,nxp
                    do l=1,nq
                        do o=1,4
                            m = octants_out(o)
                            IG(i,j,nzi,l,m) = IG(i,j,nzp,l,m)
                        end do
                    end do
                end do
            end do
        end if
    end subroutine DOMbound_out
    !========================================================================================
    !----------------------------------------------------------
    !          Calculates the Radiant intensity 
    !        at the spatial mesh (i,j,k) for DOM
    !---------------------- Description -----------------------
    ! <INPUT>
    ! IG    - Radiant intensity (i,j,k,l,m)
    ! Wq    - Sn quadrature weigth quadrature
    ! nq    - Number of quadratures per octant
    ! <OUTPUT>
    ! Grad  - Radiant intensity (i,j,k)
    ! <INTERNAL>
    ! GSUM  - Summation auxiliary variable
    !----------------------------------------------------------
    subroutine G_DOM(Grad,IG)
        integer :: i,j,k,l,m
        double precision :: GSUM,Grad(nxi,nyi,nzi),IG(nxi,nyi,nzi,nq,8)
        do k=1,nzi
            do i=1,nxi
                do j=1,nyi
                    GSUM = 0.0d0
                    do l=1,nq
                        do m=1,8
                            GSUM = GSUM + IG(i,j,k,l,m)*Wq(l)
                        end do
                    end do
                    Grad(i,j,k) = GSUM
                end do
            end do
        end do
    end subroutine G_DOM
    !========================================================================================
    !-----------------------------------------------------------
    ! Calculates all wall fluxes using Discrete Ordinates Method 
    !---------------------- Description ------------------------
    ! <INPUT>
    ! IG          - Radiant intensity (i,j,k,l,m)
    ! IBlack      - Blackbody radiation intensity
    ! mux         - x ordinate direction quadrature
    ! etay        - y ordinate direction quadrature
    ! xiz         - z ordinate direction quadrature
    ! epsilon_rad - emissivities at the walls and domain(gas)
    ! <INTERNAL>
    ! oct_m       - Downstream outflow intensity index
    ! oct_p       - Upstream outflow intensity index
    ! MSUM        - Downstream summation auxiliary variable
    ! PSUM        - Upstream summation auxiliary variable
    ! <OUTPUT>
    ! Q_radw      - Walls radiative fluxes
    ! <EXTERNAL ROUTINES>
    ! IGSDOM      - Calculates the summation of radiant intensity 
    !-----------------------------------------------------------
    subroutine wall_fluxes_DOM(IG)
        integer :: i,j,k,oct_p(4),oct_m(4)
        double precision :: IG(nxi,nyi,nzi,nt,np),MSUM,PSUM
        !---------------- x-direction walls ----------------!
        oct_m = (/2,3,6,7/); oct_p = (/1,4,5,8/)
        do k=2,nzp
            do j=2,nyp
                call IGSDOM(IG,1,j,k,oct_m,mux,MSUM)
                call IGSDOM(IG,1,j,k,oct_p,mux,PSUM)
                Q_radw(1,j,k)   =  PSUM - MSUM!West Wall
                call IGSDOM(IG,nxi,j,k,oct_m,mux,MSUM)
                call IGSDOM(IG,nxi,j,k,oct_p,mux,PSUM)
                Q_radw(nxi,j,k) =  PSUM - MSUM!East Wall
            end do
        end do
        !---------------- y-direction walls ----------------!
        if(DIMEN == 2 .or. DIMEN == 3)then
            oct_m = (/3,4,7,8/); oct_p = (/1,2,5,6/)
            do k=2,nzp
                do i=2,nxp
                    call IGSDOM(IG,i,1,k,oct_m,etay,MSUM)
                    call IGSDOM(IG,i,1,k,oct_p,etay,PSUM)
                    Q_radw(i,1,k)   =  PSUM - MSUM!South Wall
                    call IGSDOM(IG,i,nyi,k,oct_m,etay,MSUM)
                    call IGSDOM(IG,i,nyi,k,oct_p,etay,PSUM)
                    Q_radw(i,nyi,k) =  PSUM - MSUM!North Wall                
                end do 
            end do 
        end if
        !---------------- z-direction walls ----------------!
        if(DIMEN == 3)then
            oct_m = (/5,6,7,8/); oct_p = (/1,2,3,4/)
            do j=2,nyp
                do i=2,nxp
                    call IGSDOM(IG,i,j,1,oct_m,xiz,MSUM)
                    call IGSDOM(IG,i,j,1,oct_p,xiz,PSUM)
                    Q_radw(i,j,1)   =  PSUM - MSUM!Bottom Wall
                    call IGSDOM(IG,i,j,nzi,oct_m,xiz,MSUM)
                    call IGSDOM(IG,i,j,nzi,oct_p,xiz,PSUM)
                    Q_radw(i,j,nzi) =  PSUM - MSUM!Top Wall
                end do
            end do
        end if
    end subroutine wall_fluxes_DOM
    !========================================================================================
end module radiation
