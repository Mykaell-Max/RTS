module start_variables

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
    use input
    use functions
    use boundary_conditions
    use scatteringdata
    use absorptiondata
    implicit none
    
    procedure(Temp_field), pointer :: sub_function => NULL()
    
    contains
    !========================================================================================
    !----------------------------------------------------------
    !              Start and setup all variables
    !---------------------- Description -----------------------
    ! <INPUT>
    ! a_bxT   - a boundary coefficient x-direction walls
    ! a_byT   - a boundary coefficient y-direction walls
    ! a_bzT   - a boundary coefficient z-direction walls
    ! b_bxT   - b boundary coefficient x-direction walls
    ! b_byT   - b boundary coefficient y-direction walls
    ! b_bzT   - b boundary coefficient z-direction walls
    ! g_bxT   - g boundary coefficient x-direction walls
    ! g_byT   - g boundary coefficient y-direction walls
    ! g_bzT   - g boundary coefficient z-direction walls
    ! K_term  - Thermal conductivity     
    ! T_energy- Temperature
    ! <EXTERNAL ROUTINES>
    ! input_data     - Loads the input parameters from input.rts
    ! dim_sizes      - Setup dimension and sizes
    ! constants      - Setup the global constants
    ! mesh           - Builds the spatial mesh
    ! initialization - Set data for all internal nodes
    ! rad_start      - Start and setup the radiative variables
    !----------------------------------------------------------
    subroutine start()
        call input_data !Read the input file
        call dim_sizes  !Setup of dimension and sizes
        call constants  !Setup all constants
        call mesh       !Grid constructor
        !Allocation of all thermal variables
        allocate(a_bxT(2,nyi,nzi),b_bxT(2,nyi,nzi),g_bxT(2,nyi,nzi),&
                 a_byT(2,nxi,nzi),b_byT(2,nxi,nzi),g_byT(2,nxi,nzi),&
                 a_bzT(2,nxi,nyi),b_bzT(2,nxi,nyi),g_bzT(2,nxi,nyi))
        allocate(K_term(nxi,nyi,nzi),T_energy(nxi,nyi,nzi))
        call rad_start      !Allocation of all radiative variables
        call initialization !Initializes the key variables
    end subroutine start
    !========================================================================================
    !-----------------------------------------------------------------------------------
    !                      Start and setup all radiative variables
    !----------------------------------- Description -----------------------------------
    ! <INPUT>
    ! G         - Radiant intensity (i,j,k)
    ! S_rad     - Divergent radiative flux
    ! IBlack    - Blackbody radiation intensity
    ! sigma     - Scattering coefficient
    ! beta      - Extinction coefficient
    ! cappa     - Absorption coefficient 
    ! Q_radw    - Walls radiative fluxes (only the wall positions are filled)
    ! a_bxG     - a boundary coefficient x-direction walls
    ! a_byG     - a boundary coefficient y-direction walls
    ! a_bzG     - a boundary coefficient z-direction walls
    ! b_bxG     - b boundary coefficient x-direction walls
    ! b_byG     - b boundary coefficient y-direction walls
    ! b_bzG     - b boundary coefficient z-direction walls
    ! g_bxG     - g boundary coefficient x-direction walls
    ! g_byG     - g boundary coefficient y-direction walls
    ! g_bzG     - g boundary coefficient z-direction walls
    ! Ax        - x-direction coefficients for FAM
    ! Ay        - y-direction coefficients for FAM
    ! Az        - z-direction coefficients for FAM
    ! Axd       - x-direction coefficients for FAM
    ! Ayd       - y-direction coefficients for FAM
    ! Azd       - z-direction coefficients for FAM
    ! dcx       - Integral over the solid angle in the x-direction
    ! dcy       - Integral over the solid angle in the y-direction
    ! dcz       - Integral over the solid angle in the z-direction
    ! dco       - Integral over the solid angle
    ! volom     - Volume of the spatial-angular cell (i,j,k,l,m)
    ! IG        - Radiant intensity (i,j,k,l,m)
    ! phase_f   - Scattering phase function for FAM
    ! phase_d   - Scattering phase function for DOM
    ! BBF       - blackbody radiation fraction array
    ! <INTERNAL>
    ! nsb       - number of spectral bands (nsb = 1 Gray , nsb = 5 WSGG)
    ! <EXTERNAL ROUTINES>
    ! emissivity_walls        - Construct the emissivity field variable
    ! angular_grid            - Construct the angular grid for FAM
    ! angular_coefficients    - Calculates the angular coefficients for FAM
    ! anisotropic_dom         - Calculates the anisotropic scattering phase function for DOM
    ! anisotropic_fam         - Calculates the anisotropic scattering phase function for FAM
    ! quadrature_sets         - Select the quadrature sets and weigths for the DOM
    ! orthogonal_coefficients - Calculates the orthogonal coefficients for DOM
    ! WSGG_polynomials        - Loads the WSGG coefficients for a CO2-H2O mixture
    !-----------------------------------------------------------------------------------
    subroutine rad_start()
        ! Check if gas properties are being used
        if(gas_prop)then 
            ! If using WSGG model, call WSGG_polynomials subroutine
            if(wsgg_model) call WSGG_polynomials !Weighted Sum of Gray Gases Module Polynomials
            ! Allocate species variables
            allocate(Y_species(nxi,nyi,nzi,2), P_species(nxi,nyi,nzi)) !Species Variables
        end if
        !-----------------------------------------------------------------------------------
        ! Allocate main variables
        !-----------------------------------------------------------------------------------
        allocate(sigma(nxi,nyi,nzi),beta(nxi,nyi,nzi),cappa(nxi,nyi,nzi),Q_radw(nxi,nyi,nzi))
        allocate(G(nxi,nyi,nzi),S_rad(nxi,nyi,nzi),IBlack(nxi,nyi,nzi),epsilon_rad(nxi,nyi,nzi))
        ! Allocate variables for non-gray radiation
        if(nongray_flag) allocate(BBF(nxi,nyi,nzi,nsb),cappaBND(nxi,nyi,nzi,nsb),GBND(nxi,nyi,nzi,nsb))
        ! Call emissivity_walls subroutine to construct the emissivity field variable
        call emissivity_walls !Construct the emissivity field variable
        !-----------------------------------------------------------------------------------
        !P1 Method
        !-----------------------------------------------------------------------------------
        if(model_stc == 1) then
            allocate(a_bxG(2,nyi,nzi),b_bxG(2,nyi,nzi),g_bxG(2,nyi,nzi),&
                     a_byG(2,nxi,nzi),b_byG(2,nxi,nzi),g_byG(2,nxi,nzi),&
                     a_bzG(2,nxi,nyi),b_bzG(2,nxi,nyi),g_bzG(2,nxi,nyi))
            allocate(GAMA(nxi,nyi,nzi))
            if(aniso_flag .eqv. .false.) C_rad = 0.0d0
        !-----------------------------------------------------------------------------------
        !Discrete Ordinates Method
        !-----------------------------------------------------------------------------------
        else if(model_stc == 2) then
            call quadrature_sets
            allocate(IG(nxi,nyi,nzi,nq,8),phase_d(nq,nq))
            if(nongray_flag) allocate(IGBND(nxi,nyi,nzi,nq,8,nsb))
            allocate(Axd(nyi,nzi,nq),Ayd(nxi,nzi,nq),Azd(nxi,nyi,nq))
            call orthogonal_coefficients
            if(aniso_flag) then
                call anisotropic_dom
            else
                phase_d = 1.0d0
            end if
        !-----------------------------------------------------------------------------------
        !Finite Angle Method
        !-----------------------------------------------------------------------------------
        else if(model_stc == 3) then
            allocate(Ax(nyi,nzi,nt,np),Ay(nxi,nzi,nt,np),Az(nxi,nyi,nt,np),phase_f(nt,np,nt,np))
            allocate(IG(nxi,nyi,nzi,nt,np),volom(nxi,nyi,nzi,nt,np))
            if(nongray_flag) allocate(IGBND(nxi,nyi,nzi,nt,np,nsb))
            allocate(dcx(nt,np),dcy(nt,np),dcz(nt,np),dco(nt,np))
            call angular_grid
            call angular_coefficients
            if(aniso_flag) then
                call anisotropic_fam
            else
                phase_f = 1.0d0
            end if
        end if
        !-----------------------------------------------------------------------------------
        if(scheme_flag == 1)then
            alp_r = 1.0d0
        else if(scheme_flag == 2)then
            alp_r = 0.5d0
        end if
        !-----------------------------------------------------------------------------------
    end subroutine rad_start    
    !========================================================================================
    !----------------------------------------------------------
    !              Initializes the key variables
    !---------------------- Description -----------------------
    ! <EXTERNAL ROUTINES>
    ! initi_func - Initializes a given variable for all internal nodes
    !----------------------------------------------------------
    ! ID 1  - (T_energy) temperature for all nodes (K)
    ! ID 2  - (G) incident radiation for all nodes (w/m^2)
    ! ID 3  - (S_rad) divergent radiative flux (w/m^3)
    ! ID 4  - (IBlack) blackbody radiation intensity
    ! ID 5  - (cappa) absorption coefficient for all nodes (m^-1)
    ! ID 6  - (sigma) scattering coefficient for all nodes (m^-1)
    ! ID 7  - (beta) extinction coefficient for all nodes (m^-1)
    ! ID 8  - (K_term) thermal conductivity (W/m*K)
    ! ID 11 - (P_species) Pressure of the mixture (atm)
    ! ID 12 - (Y_species) CO2 mole fraction
    ! ID 13 - (Y_species) H2O mole fraction
    !----------------------------------------------------------
    subroutine initialization()
        integer:: IBND
        ! Initialize T_energy and G
        call initi_func(T_energyc,T_energy,1) !T_energy initialization
        call initi_func(Gc,G,2)               !G initialization
         ! If nongray_flag is set, initialize GBND
        if(nongray_flag)then
            do IBND = 1,nsb
                GBND(:,:,:,IBND) = G
            end do
        end if
        ! If gas_prop is false, initialize cappa
        if(gas_prop.eqv..false.) call initi_func(k_rad,cappa,5)
        ! Initialize sigma and K_term
        call initi_func(sig_rad,sigma,6)                        !sigma initialization
        call initi_func(K_termc,K_term,8)                       !K_term initialization
        ! If gas_prop is true, initialize P_g, XCO2_g, and XH2O_g
        if(gas_prop)then
            call initi_func(P_g,P_species,11)                   !Pressure of the mixture (atm)
            call initi_func(XCO2_g,Y_species(:,:,:,1),12)       !CO2 mole fraction
            call initi_func(XH2O_g,Y_species(:,:,:,2),13)       !H2O mole fraction
            ! If nongray_flag is false and wsgg_model is true, calculate Lm
            if(nongray_flag.eqv..false..and. wsgg_model)then
                Lm = 1.8d0*(lx*ly*lz)/(lx*ly + lx*lz + ly*lz)
            end if
        end if
        ! Set Q_radw to 0.0
        Q_radw = 0.0d0
        ! If model_stc is 2 or 3, initialize IG and IGBND if nongray_flag is true
        if(model_stc == 2 .or. model_stc == 3)then
            IG = 0.0d0 
            if(nongray_flag) IGBND = 0.0d0 
        end if
        !----------------------------------------------------------
    end subroutine initialization
    !========================================================================================
    !----------------------------------------------------------
    !          Select initialization function pointer
    !---------------------- Description -----------------------
    ! <INPUT>
    ! ID          - Variable ID
    ! <EXTERNAL ROUTINES>
    ! Temp_field  - Temperature space-function
    ! G_field     - Incident radiation space-function
    ! cappa_field - Absorption coefficient space-function
    ! sigma_field - Scattering coefficient space-function
    ! Kterm_field - Thermal conductivity space-function
    ! Pk_field    - Species total pressure space-function
    ! Y_CO2field  - CO2 mole fraction space-function
    ! Y_H2Ofield  - H2O mole fraction space-function
    !----------------------------------------------------------
    subroutine select_exact_pointer(ID)
        integer :: ID
        select case (ID)
            case (1)
            sub_function => Temp_field
            case (2)
            sub_function => G_field
            case (5)
            sub_function => cappa_field
            case (6)
            sub_function => sigma_field
            case (8)
            sub_function => Kterm_field
            case (11)
            sub_function => Pk_field
            case (12)
            sub_function => Y_CO2field
            case (13)
            sub_function => Y_H2Ofield
            case default
            print *, "Invalid ID value"
            stop
        end select
    end subroutine select_exact_pointer
    !========================================================================================
    !----------------------------------------------------------
    !   Initializes a given variable for all internal nodes
    !---------------------- Description -----------------------
    ! <INPUT>
    ! ID         - Variable ID
    ! VAR        - constant value to be used at the initialization
    ! variable   - variable field to be initialized
    ! field_flag - flag to indicate a uniform field or not
    ! <EXTERNAL ROUTINES>
    ! select_exact_pointer - Select initialization function pointer
    !----------------------------------------------------------
    subroutine initi_func(VAR,variable,ID)
        double precision, intent (in out) :: VAR,variable(:,:,:)
        integer :: i,j,k,ID
        call select_exact_pointer(ID)
        if (field_flag(ID))then
            variable = VAR     !constant throughout 
        else
            do k=2,nzp
                do j=2,nyp
                    do i=2,nxp
                        variable(i,j,k) = sub_function(i,j,k)
                    end do
                end do
            end do
        end if
        nullify(sub_function)
    end subroutine initi_func
    !========================================================================================
    !----------------------------------------------------------
    !                Setup dimension and sizes
    !
    !     |---|---|---|---|---|---|---|---|---|---|---|---|
    !     1   2                                      nxp nxi
    ! NOTE: 1 and nxi are information at the boundaries
    ! NOTE: the range between [2-nxp] define the points inside 
    ! NOTE: the domain, all loops must occur in that range
    !---------------------- Description -----------------------
    ! <INPUT>
    ! nx      - Number of control volumes in x-direction
    ! ny      - Number of control volumes in y-direction
    ! nz      - Number of control volumes in z-direction
    ! <OUTPUT>
    ! nxp     - Auxiliary index
    ! nyp     - Auxiliary index
    ! nzp     - Auxiliary index
    ! nxi     - Number of control volumes in 
    !           x-direction considering the walls
    ! nyi     - Number of control volumes in 
    !           y-direction considering the walls
    ! nzi     - Number of control volumes in 
    !           z-direction considering the walls
    !----------------------------------------------------------
    subroutine dim_sizes()
        if(DIMEN == 1)then
            ny = 1; nz = 1
        elseif(DIMEN == 2)then
            nz = 1
        end if
        !-----------
        nxp = nx + 1
        nyp = ny + 1
        nzp = nz + 1
        !-----------
        nxi = nx + 2
        nyi = ny + 2
        nzi = nz + 2
        !-----------
    end subroutine dim_sizes
    !========================================================================================
    !----------------------------------------------------------
    !              Setup the global constants
    ! Made by:     G.S.Rodrigues (October 2020)
    ! Modified by: G.S.Rodrigues (May 2021)
    !---------------------- Description -----------------------
    ! <OUTPUT>
    ! small   - A small number
    ! big     - A big number
    ! boltz   - Stefan boltzmann constant
    ! PI      - pi number
    ! PIBY2   - pi/2
    ! PI32    - 3pi/2
    ! PI4     - 4pi
    !----------------------------------------------------------
    subroutine constants
        small= 1.0d-10
        big  = 1.0d10
        boltz= 5.670374419d-8
        PI = 4.0d0*atan(1.0d0)
        PIBY2= PI*0.5d0
        PI32 = 1.5d0*PI
        PI4  = 4.0d0*PI
    end subroutine constants
    !========================================================================================
    !----------------------------------------------------------
    !                Build the spatial mesh
    !---------------------- Description -----------------------
    ! <OUTPUT>
    ! x       - x-direction grid
    ! y       - y-direction grid
    ! z       - z-direction grid
    ! xc      - x value at the control-volume center
    ! yc      - y value at the control-volume center
    ! zc      - z value at the control-volume center
    ! dxp     - Control-volume length
    ! dyp     - Control-volume height
    ! dzp     - Control-volume width 
    ! <INTERNAL>
    ! dx      - Control-volume length
    ! dy      - Control-volume height
    ! dz      - Control-volume width 
    ! <EXTERNAL ROUTINES>
    ! pathfinder - Build the non-uniform mesh following the paths
    !----------------------------------------------------------
    subroutine mesh()
        integer :: i,j,k
        double precision :: dx,dy,dz
        allocate(x(nxp),y(nyp),z(nzp))
        x = 0.0d0; y = 0.0d0; z = 0.0d0
        !------------ x-direction uniform grid ------------!
        dx = lx/nx
        do i = 1,nxp
            x(i) = (i-1.0d0)*dx
        end do
        !------------ y-direction uniform grid ------------!
        if(DIMEN == 2 .or. DIMEN == 3)then
            dy = ly/ny
            do j = 1,nyp
                y(j) = (j-1.0d0)*dy
            end do
        else
            y =(/0.0d0, lx*0.05d0/)
        end if
        !------------ z-direction uniform grid ------------!
        if(DIMEN == 3)then
            dz = lz/nz
            do k = 1, nzp
                z(k) = (k-1.0d0)*dz
            end do
        else
            z =(/0.0d0, 1.0d-5/)
        end if
        !Path Additions
        if(path_flag) call pathfinder(dx,dy,dz)
        !Espatial position at the center of the control volumes
        allocate(xc(nxi),yc(nyi),zc(nzi))
        xc  = 0.0d0; yc  = 0.0d0; zc  = 0.0d0
        !----- x-direction grid at the center of C.V. -----!
        do i=2,nxp
            xc(i) = (x(i-1) + x(i))*0.5d0
        end do
        !----- y-direction grid at the center of C.V. -----!
        do j=2,nyp
            yc(j) = (y(j-1) + y(j))*0.5d0
        end do
        !----- z-direction grid at the center of C.V. -----!
        do k=2,nzp
            zc(k) = (z(k-1) + z(k))*0.5d0
        end do
        !--------------------------------------------------!
        xc(1)   = x(1);   yc(1)   =  y(1);   zc(1)   = z(1)
        xc(nxi) = x(nxp); yc(nyi) =  y(nyp); zc(nzi) = z(nzp)
        !NOTE: xc(1),yc(1),Zc(1),xc(nxi),yc(nyi) and Zc(nzi)
        !      are the phisical positions on the walls
        !--------------------------------------------------!
        !Direfences vectors (dx,dy,dz) 
        allocate(dxp(nxi),dyp(nyi),dzp(nzi),vol(nxi,nyi,nzi))
        dxp = 0.0d0; dyp = 0.0d0; dzp = 0.0d0; vol = 0.0d0
        !--------- x-direction direfence vectors ----------!
        do i=2,nxp
            dxp(i) = x(i) - x(i-1)
        end do
        !--------- y-direction direfence vectors ----------!
        if(DIMEN == 2 .or. DIMEN == 3) then
            do j=2,nyp
                dyp(j) = y(j) - y(j-1)
            end do
        else
            dyp = 1.0d0                   !for 1D domain
            if(model_stc.eq.1) dyp = big  !only for P1
        end if
        !--------- z-direction direfence vectors ----------!
        if(DIMEN == 3) then
            do k=2,nzp
                dzp(k) = z(k) - z(k-1)
            end do
        else
            dzp = 1.0d0                   !for 2D/1D domain
            if(model_stc.eq.1) dzp = big  !only for P1
        end if
        !--------------------------------------------------!
        !WARNING: dxp(1)  ,dyp(1)  ,dzp(1)
        !         dxp(nxi),dyp(nyi),dzp(nzi) are not used!
        !------------------- Cell Volume ------------------!
        do i=2,nxp
            do j=2,nyp
                do k=2,nzp
                    vol(i,j,k) = dxp(i)*dyp(j)*dzp(k)
                end do
            end do
        end do
    end subroutine mesh
    !========================================================================================
    !----------------------------------------------------------
    !     Build the non-uniform mesh following the paths
    !---------------------- Description -----------------------
    ! <INPUT>
    ! dx      - Control-volume length
    ! dy      - Control-volume height
    ! dz      - Control-volume width 
    ! Lpath   - The paths coordinates
    ! <OUTPUT>
    ! x       - x-direction grid
    ! y       - y-direction grid
    ! z       - z-direction grid
    ! <INTERNAL>
    ! x_path  - x-direction grid for path i
    ! y_path  - y-direction grid for path i
    ! z_path  - z-direction grid for path i
    ! dx_path - Path control-volume length
    ! dy_path - Path control-volume height
    ! dz_path - Path control-volume width
    ! nx_path - Number of control volumes in x-direction for path i
    ! ny_path - Number of control volumes in y-direction for path i
    ! nz_path - Number of control volumes in z-direction for path i        
    ! <EXTERNAL ROUTINES>
    ! rm_dups    - Remove array duplicate values
    ! shell_Sort - Sorts values in ascending order 
    ! dim_sizes  - Setup dimension and sizes
    !----------------------------------------------------------
    subroutine pathfinder(dx,dy,dz)
        implicit none
        integer :: i,j,k,p,nx_path,ny_path,nz_path
        double precision :: dx,dy,dz,dx_path,dy_path,dz_path
        double precision, allocatable,dimension(:) :: x_path,y_path,z_path
        !---------- x-direction non-uniform grid ----------!
        do p = 1,npath
            dx_path = dx/Lpath(7,p)
            Nx_path = int((Lpath(2,p) - Lpath(1,p))/dx_path)
            allocate(x_path(nx_path + 1))
            do i = 1, nx_path + 1
                x_path(i) = Lpath(1,p) + (i-1.0d0)*dx_path
            end do
            x = [x,x_path]      !Concatenation
            deallocate(x_path)
        end do
        x = nint(x*100)/100.0d0 !To avoid duplicate values
        call rm_dups(x)         !Removing duplicate values
        call shell_Sort(x)      !Sorting
        nx = size(x) - 1
        !---------- y-direction non-uniform grid ----------!
        if(DIMEN == 2 .or. DIMEN == 3)then
            do p = 1,npath
                dy_path = dy/Lpath(7,p)
                ny_path = int((Lpath(4,p) - Lpath(3,p))/dy_path)
                allocate(y_path(ny_path + 1))
                do j = 1, ny_path + 1
                    y_path(j) = Lpath(3,p) + (j-1.0d0)*dy_path
                end do
                y = [y,y_path]  !Concatenation
                deallocate(y_path)            
            end do
            y = nint(y*100)/100.0d0!To avoid duplicate values
            call rm_dups(y)     !Removing duplicate values
            call shell_Sort(y)  !Sorting
            ny = size(y) - 1
        end if
        !---------- z-direction non-uniform grid ----------!
        if(DIMEN == 3)then
            do p = 1,npath
                dz_path = dz/Lpath(7,p)
                nz_path = int((Lpath(6,p) - Lpath(5,p))/dz_path)
                allocate(z_path(nz_path + 1))
                do k = 1, nz_path + 1
                    z_path(k) = Lpath(5,p) + (k-1.0d0)*dz_path
                end do
                z = [z,z_path]  !Concatenation
                deallocate(z_path)
            end do
            z = nint(z*100)/100.0d0 !To avoid duplicate values
            call rm_dups(z)     !Removing duplicate values
            call shell_Sort(z)  !Sorting
            nz = size(z) - 1
        end if
        !--------------------------------------------------!
        call dim_sizes  !Corrects the fundamental sizes
    end subroutine pathfinder
    !========================================================================================
    !----------------------------------------------------------
    !              Construct the angular grid for
    !                  Finite Angle Method
    !---------------------- Description -----------------------
    ! <INPUT>
    !
    ! <OUTPUT>
    ! dcx     - Integral over the solid angle in the x-direction
    ! dcy     - Integral over the solid angle in the y-direction
    ! dcz     - Integral over the solid angle in the z-direction
    ! dco     - Integral over the solid angle
    ! volom   - Volume of the spatial-angular cell (i,j,k,l,m)
    ! <INTERNAL>
    ! lp      - phi-direction length  2pi
    ! lt      - theta-direction length pi
    ! theta   - theta-direction grid - value at the control-angle face
    ! phi     - phi-direction   grid - value at the control-angle face
    ! dphi    - Control-angle length 
    ! dtheta  - Control-angle height
    ! dxp     - Control-volume length
    ! dyp     - Control-volume height
    ! dzp     - Control-volume width
    ! D_theta - Auxiliary variable
    !----------------------------------------------------------
    subroutine angular_grid
        double precision :: lp,lt,dphi,dtheta,D_theta
        integer :: i,j,k,l,m
        allocate(theta(nt+1),phi(np+1))
        lp = 2.0d0*PI
        lt = PI
        dphi   = lp/np
        dtheta = lt/nt
        !--------------- phi-direction grid ---------------!
        do m=1,np+1
            phi(m) = (m-1.0d0)*dphi
        end do
        !-------------- theta-direction grid --------------!
        do l=1,nt+1
            theta(l) = (l-1.0d0)*dtheta
        end do
        !--------------------------------------------------!
        !Dividing indexes of the angular mesh
        call angular_indexes
        !Integral over the solid angles
        do l=1,nt
            D_theta = 0.5d0*(dtheta - cos(theta(l+1))*sin(theta(l+1)) + cos(theta(l))*sin(theta(l)))
            do m=1,np
                dco(l,m) = (cos(theta(l)) - cos(theta(l+1)))*dphi
                dcx(l,m) = (sin(phi(m+1)) - sin(phi(m)))*D_theta
                dcy(l,m) = (cos(phi(m)) - cos(phi(m+1)))*D_theta
                dcz(l,m) = 0.5d0*dphi*((sin(theta(l+1)))**2.0d0 - (sin(theta(l)))**2.0d0)
            end do
        end do
        !volom = vol*dco
        do i=2,nxp
            do j=2,nyp
                do k=2,nzp
                    do l=1,nt
                        do m=1,np
                            volom(i,j,k,l,m) = vol(i,j,k)*dco(l,m)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine angular_grid
    !========================================================================================
    !----------------------------------------------------------
    !              Determines the dividing indexes
    !                   of the angular mesh
    !---------------------- Description -----------------------
    ! <INTERNAL>
    ! thetam  - theta value at the control-angle center
    ! phim    - phi value at the control-angle center
    ! <OUTPUT>
    ! P2      - First  azimuthal region index (0 < phi < PIBY2)
    ! P3      - Second azimuthal region index (PIBY2 < phi < PI)
    ! P4      - Third  azimuthal region index (PI < phi < PI32)
    ! T2      - First polar region index (0 < theta < PIBY2)
    !----------------------------------------------------------
    subroutine angular_indexes
        integer :: l,m
        double precision :: phim, thetam
        !------------- phi-direction indexes --------------!
        do m=1,np
            phim = 0.5d0*(PHI(m) + PHI(m+1))
            if(phim < PIBY2)then 
                P2 = m
            elseif(phim < PI)then 
                P3 = m
            elseif(phim < PI32)then 
                P4 = m
            end if
        end do
        !------------- theta-direction indexes ------------!
        do l=1,nt
            thetam = 0.5d0*(theta(l) + theta(l+1))
            if(thetam < PIBY2)then 
                T2 = l
            end if
        end do
    end subroutine angular_indexes
    !========================================================================================
    !----------------------------------------------------------
    !             Calculates the internal angular 
    !          coefficients for Finite Angle Method
    !---------------------- Description -----------------------
    ! <INPUT>
    ! dcx     - Integral over the solid angle in the x-direction
    ! dcy     - Integral over the solid angle in the y-direction
    ! dcz     - Integral over the solid angle in the z-direction
    ! <OUTPUT>
    ! Ax      - x-direction coefficients for FAM
    ! Ay      - y-direction coefficients for FAM
    ! Az      - z-direction coefficients for FAM
    ! <INTERNAL>
    ! dxp     - Control-volume length
    ! dyp     - Control-volume height
    ! dzp     - Control-volume width
    !----------------------------------------------------------
    subroutine angular_coefficients
        integer :: i,j,k,l,m
        Ax = 0.0d0; Ay = 0.0d0; Az = 0.0d0
        !------------ x-direction coefficients ------------!
        do k=2,nzp
            do j=2,nyp
                do l=1,nt
                    do m=1,np
                        Ax(j,k,l,m) = abs(dcx(l,m))*dyp(j)*dzp(k)
                    end do
                end do
            end do
        end do
        !------------ y-direction coefficients ------------!
        if(DIMEN == 2 .or. DIMEN == 3)then
            do k=2,nzp
                do i=2,nxp
                    do l=1,nt
                        do m=1,np
                            Ay(i,k,l,m) = abs(dcy(l,m))*dxp(i)*dzp(k)
                        end do
                    end do
                end do
            end do
        end if
        !------------ z-direction coefficients ------------!
        if(DIMEN == 3)then
            do j=2,nyp
                do i=2,nxp
                    do l=1,nt
                        do m=1,np
                            Az(i,j,l,m) = abs(dcz(l,m))*dxp(i)*dyp(j)
                        end do
                    end do
                end do
            end do
        end if
    end subroutine angular_coefficients
    !========================================================================================
    !----------------------------------------------------------
    !    Select the quadrature sets and weigths for the DOM
    !---------------------- Description -----------------------
    ! <INPUT>
    ! QM_stc  - Selected the model for the quadrature
    ! N_quad  - Quadrature order
    ! <OUTPUT>
    ! Wq      - Quadrature weights
    ! <EXTERNAL ROUTINES>
    ! Sn_sets          - Select the Sn quadrature sets
    ! Tn_sets          - Calculates the Tn quadrature sets
    ! Qn_sets          - Calculates the Qn quadrature sets
    ! Oned_quadratures - Adapts the quadrature sets for 1D domain
    !----------------------------------------------------------
    subroutine quadrature_sets
        if(QM_stc == 1)then
            call SN_sets(N_quad)
        else if(QM_stc == 2)then
            call TN_sets(N_quad,nq)
        else if(QM_stc == 3)then
            call QN_sets(N_quad)
        end if
        !Adapts the quadrature sets for 1D or 2D domain
        if(DIMEN == 1)then
            call Oned_quadratures
        else if(DIMEN == 2)then
            Wq = 2.0d0*Wq
        end if
    end subroutine quadrature_sets
    !========================================================================================
    !----------------------------------------------------------
    !          Select the Sn quadrature sets for DOM
    !---------------------- Description -----------------------
    !  For more details see "Three-dimensional radiative heat
    !   transfer solutions by the discrete-ordinates method"
    !                    by Fiveland (1988)
    !            doi: https://doi.org/10.2514/3.105
    !----------------------------------------------------------
    ! <INPUT>
    ! SN_flag - Discretization order
    ! <OUTPUT>
    ! mux   - x ordinate direction quadrature
    ! etay  - y ordinate direction quadrature
    ! xiz   - z ordinate direction quadrature
    ! Wq    - Quadrature weights
    ! nq    - Number of quadratures per octant
    !----------------------------------------------------------
    subroutine Sn_sets(SN_flag)
        integer :: SN_flag
        select case (SN_flag)
        case (2)   ! S2
            mux = (/0.5773503/)
            etay= (/0.5773503/)
            xiz = (/0.5773503/)
            Wq  = (/1.5707963/)
            nq  = 1
        case (4)   ! S4
            mux = (/0.2958759,0.9082483,0.2958759/)
            etay= (/0.2958759,0.2958759,0.9082483/)
            xiz = (/0.9082483,0.2958759,0.2958759/)
            Wq  = (/0.5235987,0.5235987,0.5235987/)
            nq  = 3
        case (6)   ! S6
            mux = (/0.1838670,0.6950514,0.9656013,0.183867,0.69505140,0.1838670/)
            etay= (/0.1838670,0.1838670,0.1838670,0.6950514,0.6950514,0.9656013/)
            xiz = (/0.9656013,0.6950514,0.1838670,0.6950514,0.1838670,0.1838670/)
            Wq  = (/0.1609517,0.3626469,0.1609517,0.3626469,0.3626469,0.1609517/)
            nq  = 6
        case (8)   ! S8
            mux = (/0.1422555,0.5773503,0.8040087,0.9795543,0.1422555,0.5773503,0.8040087,0.1422555,0.5773503,0.1422555/)
            etay= (/0.1422555,0.1422555,0.1422555,0.1422555,0.5773503,0.5773503,0.5773503,0.8040087,0.8040087,0.9795543/)
            xiz = (/0.9795543,0.8040087,0.5773503,0.1422555,0.8040087,0.5773503,0.1422555,0.5773503,0.1422555,0.1422555/)
            Wq  = (/0.1712359,0.0992284,0.0992284,0.1712359,0.0992284,0.4617179,0.0992284,0.0992284,0.0992284,0.1712359/)
            nq  = 10
        case (10)   ! S10
            mux = (/0.9809754,0.8523177,0.8523177,0.7004129,0.7004129,0.7004129,0.5046889,& 
                    0.5046889,0.5046889,0.5046889,0.1372719,0.1372715,0.1372719,0.1372719,0.1372719/)
            etay= (/0.1372719,0.1372719,0.5046889,0.1372719,0.5046889,0.7004129,0.1372719,&
                    0.5046889,0.7004129,0.8523177,0.1372719,0.5046889,0.7004129,0.8523177,0.9809754/)
            xiz = (/0.1372719,0.5046889,0.1372719,0.7004129,0.5046889,0.1372719,0.8523177,&
                    0.7004129,0.5046889,0.1372719,0.9809754,0.8523177,0.7004129,0.5046889,0.1372719/)
            Wq  = (/0.0944411,0.1483950,0.1483950,0.0173701,0.1149972,0.0173701,0.1483950,&
                    0.1149972,0.1149972,0.1483950,0.0944411,0.1483950,0.0173701,0.1483950,0.0944411/)
            nq  = 15
        case default
            print *, "Invalid SN_flag flag"
            stop
        end select
    end subroutine Sn_sets
    !========================================================================================
    !----------------------------------------------------------
    !     Calculation of Thurgood quadrature sets for DOM
    !---------------------- Description -----------------------
    !  For more details see "The TN Quadrature Set for the 
    !     Discrete Ordinates Method" by Thurgood et al.(1995)     
    !          doi: https://doi.org/10.1115/1.2836285
    !----------------------------------------------------------
    ! <INPUT>
    ! N     - Discretization order
    ! <INTERNAL>
    ! nv    - Number of vertices
    ! nst   - Number of sub-triangles
    ! Lt    - Axis distance between two flat vertices
    ! ct    - sub-triangles counter variable
    ! ver   - Vertices of plane triangles
    ! verp  - Projected vertices in the unit sphere
    ! Cm    - Centers of gravity of plane triangles
    ! subT  - Sub-triangles vertices indicator
    ! Cm_N2 - Norm 2 of Cm
    ! ver_N2- Norm 2 of verp
    ! alpha - Solid angle verticies
    ! P     - Half perimeter
    ! T_flag- Upward/Downward triangle indicator
    ! r     - Auxiliary counter variable
    ! Tu    - Upward triangles counter variable
    ! Td    - Downward triangles counter variable
    ! Tr    - Downward triangles auxiliary counter variable
    ! aux   - Auxiliary variable
    ! <OUTPUT>
    ! mux   - x ordinate direction quadrature
    ! etay  - y ordinate direction quadrature
    ! xiz   - z ordinate direction quadrature
    ! Wq    - Quadrature weights
    !----------------------------------------------------------
    subroutine Tn_sets(N,nst)
        integer :: i,j,k,N,nst,nv,ct,r,Tu,Td,Tr,T_flag
        double precision :: Lt,Cm_N2,ver_N2,P,aux,alpha(3)
        double precision,allocatable,dimension(:,:) :: ver,verp,Cm
        integer,allocatable,dimension(:,:) :: subT
        nv = int((N+2)*(N+1)/2) !Number of vertices
        nst= N*N                !Number of sub-triangles
        Lt = 1.0d0/N !Axis distance between two flat vertices
        allocate(ver(3,nv),verp(3,nv),subT(3,nst),Cm(3,nst))
        allocate(mux(nst),etay(nst),xiz(nst),Wq(nst))
        !Determination of the basal triangle verticies
        k = 1
        do i = 1,N+1
            do j = 1,i
                ver(1,k) = 1 - (i-1)*Lt
                ver(2,k) = (j - 1)*Lt
                ver(3,k) = (i - j)*Lt
                k = k + 1
            end do
        end do
        !First Triangle Loop
        do ct = 1,3
            subT(ct,1) = ct
        end do
        !WARNING:For your own sake,DO NOT modify any of these values
        Tu = 1;Td = 0;Tr = 2;ct = 2;r = 0
        !Remaining Triangles Loop
        do i = 2,N
            T_flag = 0;Tu = Tu - 1
            do j = 1,i + r + 1
                if(T_flag == 0)then!Upward Triangles
                    subT(1,ct) =   i + Tu 
                    subT(2,ct) = 2*i + Tu 
                    subT(3,ct) = subT(2,ct) + 1
                    Tu = Tu + 1 
                    T_flag = 1
                else if(T_flag == 1)then!Downward Triangles
                    subT(1,ct) = i + Td
                    subT(2,ct) = subT(1,ct) + 1
                    subT(3,ct) = subT(2,ct) + Tr
                    Td = Td + 1 
                    T_flag = 0
                end if
                ct = ct + 1
            end do
            Tr = Tr + 1;r = r + 1
        end do
        !Determination of the centers of mass of the flat triangles
        do ct = 1,nst !Triangles Loop
            do i = 1,3!Coordinates Loop
                Cm(i,ct) = 0.0d0
                do k = 1,3 !Vertices Loop
                    Cm(i,ct) = Cm(i,ct) + ver(i,subT(k,ct))
                end do
                Cm(i,ct) = Cm(i,ct)/3.0d0
            end do
        end do
        !Determination of the discrete directions cosines
        do ct = 1,nst !Triangles Loop
            Cm_N2 = 0.0d0
            do i = 1,3 !Coordinates Loop
                Cm_N2 = Cm_N2 + Cm(i,ct)**2.0d0
            end do
            Cm_N2 = sqrt(Cm_N2)
            mux(ct) = Cm(1,ct)/Cm_N2
            etay(ct)= Cm(2,ct)/Cm_N2
            xiz(ct) = Cm(3,ct)/Cm_N2
        end do
        !Projection of verticies in the unit sphere
        do k = 1,nv
            ver_N2 = 0.0d0
            do i = 1,3 !Coordinates Loop
                ver_N2 = ver_N2 + ver(i,k)**2.0d0
            end do
            ver_N2 = sqrt(ver_N2)
            do i = 1,3 !Coordinates Loop
                verp(i,k) = ver(i,k)/ver_N2
            end do
        end do
        !Determination of the discrete weights
        do ct = 1,nst !Triangles Loop
            alpha = 0.0d0
            do i = 1,3 !Coordinates Loop
                alpha(1) = alpha(1) + verp(i,subT(1,ct))*verp(i,subT(2,ct)) !alpha_12
                alpha(2) = alpha(2) + verp(i,subT(1,ct))*verp(i,subT(3,ct)) !alpha_13
                alpha(3) = alpha(3) + verp(i,subT(2,ct))*verp(i,subT(3,ct)) !alpha_23
            end do
            alpha = acos(alpha)
            P = 0.5d0*sum(alpha)
            aux = tan(0.5d0*P)
            do i = 1,3 !Angle Loop
                aux = aux*tan(0.5d0*(P-alpha(i)))
            end do
            Wq(ct) = 4.0d0*atan(sqrt(aux))
        end do
    end subroutine Tn_sets
    !========================================================================================
    !----------------------------------------------------------
    !        Calculation of Qn quadrature sets for DOM
    !---------------------- Description -----------------------
    !  For more details see "Three-dimensional Radiation in 
    !     Absorbing-Emitting-Scattering Medium Using the 
    !  Discrete-Ordinates Approximation" by Wei et al.(1998)     
    !      doi: https://doi.org/10.1007/s11630-998-0035-8
    !----------------------------------------------------------
    ! <INPUT>
    ! N      - Discretization order
    ! <INTERNAL>
    ! alphaq - The unit sphere longitude
    ! betaq  - The unit sphere latitude
    ! thetam - Theta value in a given position
    ! phim   - phi value in a given position
    ! <OUTPUT>
    ! mux    - x ordinate direction quadrature
    ! etay   - y ordinate direction quadrature
    ! xiz    - z ordinate direction quadrature
    ! Wq     - Quadrature weights
    !----------------------------------------------------------
    subroutine Qn_sets(N)
        integer :: k,l,m,N
        double precision :: phim,thetam,alphaq(N+1),betaq(N+1)        
        nq = N*N
        allocate(mux(nq),etay(nq),xiz(nq),Wq(nq))
        phim = (PIBY2/N)
        do l=1,N+1
            alphaq(l) = (l-1)*phim
            betaq(l)  = (l-1)*phim
        end do
        k = 1
        do l=1,N
            phim = 0.5d0*(betaq(l) + betaq(l+1))
            do m=1,N
                thetam = 0.5d0*(alphaq(m) + alphaq(m+1))
                mux(k) = cos(phim)*cos(thetam)
                etay(k)= sin(phim)
                xiz(k) = cos(phim)*sin(thetam)
                Wq(k)  = (PIBY2/N)*(sin(alphaq(l+1)) - sin(alphaq(l)))
                k = k + 1
            end do
        end do
    end subroutine Qn_sets
    !========================================================================================
    !----------------------------------------------------------
    !         Adapts the quadrature sets for 1D domain
    !---------------------- Description -----------------------
    ! <OUTPUT>
    ! mux   - x ordinate direction quadrature
    ! etay  - y ordinate direction quadrature
    ! xiz   - z ordinate direction quadrature
    ! Wq    - Quadrature weights
    ! nq    - Number of quadratures per octant
    ! <INTERNAL>
    ! Wq_n  - New quadrature weights
    ! nq_n  - New number of quadratures per octant
    ! <EXTERNAL ROUTINES>
    ! rm_dups    - Remove array duplicate values
    ! shell_Sort - Sorts values in ascending order 
    !----------------------------------------------------------
    subroutine Oned_quadratures
        integer :: i,j,nq_n
        double precision,allocatable, dimension(:) :: Wq_n
        !Rounds the arrays so that functions rm_dups and shell_Sort work properly
        do i=1,nq
            etay(i) = real(etay(i), 4)
        end do
        mux = etay !Copies the ordinates values of the vertical axis to mux
        call rm_dups(mux)     !Removing duplicate values
        call shell_Sort(mux)  !Sorting
        nq_n = size(mux)
        !Adapts the quadrature sets
        allocate(Wq_n(nq_n)); Wq_n = 0.0d0
        do j=1,nq_n
            do i=1,nq
                if(etay(i) == mux(j))then 
                    Wq_n(j) = Wq_n(j) + Wq(i) !Found a match
                end if
            end do
        end do
        !Updates the quadrature values
        nq = nq_n
        deallocate(Wq); allocate(Wq(nq))
        Wq = 4.0d0*Wq_n
        etay = 0.0d0; xiz = 0.0d0 !to avoid numerical errors
    end subroutine Oned_quadratures
    !========================================================================================
    !----------------------------------------------------------
    !             Calculates the internal orthogonal 
    !        coefficients for Discrete Ordinates Method
    !---------------------- Description -----------------------
    ! <INPUT>
    ! mux     - x ordinate direction quadrature
    ! etay    - y ordinate direction quadrature
    ! xiz     - z ordinate direction quadrature
    ! nq      - Number of quadratures per octant
    ! <OUTPUT>
    ! Axd     - x-direction coefficients for DOM
    ! Ayd     - y-direction coefficients for DOM
    ! Azd     - z-direction coefficients for DOM
    ! <INTERNAL>
    ! dxp     - Control-volume length
    ! dyp     - Control-volume height
    ! dzp     - Control-volume width
    !----------------------------------------------------------
    subroutine orthogonal_coefficients
        integer :: i,j,k,l
        Axd = 0.0d0; Ayd = 0.0d0; Azd = 0.0d0
        !------------ x-direction coefficients ------------!
        do k=2,nzp
            do j=2,nyp
                do l=1,nq
                    Axd(j,k,l) = abs(mux(l))*dyp(j)*dzp(k)
                end do
            end do
        end do
        !------------ y-direction coefficients ------------!
        if(DIMEN == 2 .or. DIMEN == 3)then
            do k=2,nzp
                do i=2,nxp
                    do l=1,nq
                        Ayd(i,k,l) = abs(etay(l))*dxp(i)*dzp(k)
                    end do
                end do
            end do
        end if
        !------------ z-direction coefficients ------------!
        if(DIMEN == 3)then
            do j=2,nyp
                do i=2,nxp
                    do l=1,nq
                        Azd(i,j,l) = abs(xiz(l))*dxp(i)*dyp(j)
                    end do
                end do
            end do
        end if
    end subroutine orthogonal_coefficients
    !========================================================================================
    !----------------------------------------------------------
    !             Removes duplicates from an array
    !---------------------- Description -----------------------
    ! <INPUT>
    ! a       - array with duplicates to be removed
    ! <OUTPUT>
    ! a       - array without duplicates
    !----------------------------------------------------------
    subroutine rm_dups(a)
        implicit none        
        double precision , allocatable :: a(:)
        double precision :: res(size(a))
        integer :: i,j,k
        k = 1; res = 9.0d8
        res(1) = a(1)
        outer: do i=2,size(a)
            do j=1,k
                if (res(j) == a(i)) then
                    cycle outer !Found a match so start looking again
                end if
            end do
            k = k + 1
            res(k) = a(i)       !No match found so add it to the res
        end do outer
        deallocate(a); allocate(a(k))
        a(1:k) = res(1:k)
    end subroutine rm_dups
    !========================================================================================
    !----------------------------------------------------------
    !             Sort a array in ascending order             
    !---------------------- Description -----------------------
    ! <INPUT>
    ! a       - array to be sort in ascending order
    ! <OUTPUT>
    ! a       - ordered array
    !----------------------------------------------------------
    subroutine shell_Sort(a)
        implicit none
        integer :: i, j, increment,aux
        double precision :: temp
        double precision, intent (in out) :: a(:)
        increment = size(a) / 2
        do while (increment > 0)
            do i = increment+1, size(a)
                j = i
                temp = a(i)
                aux = j-increment
                do while (j >= increment+1 .and. a(aux) > temp)
                    a(j) = a(j-increment)
                    j = j - increment
                    !To avoid negative or null indexes
                    if(j-increment <= 0) then
                        aux = 1
                    else
                        aux = j-increment
                    end if
                end do
                a(j) = temp
            end do
            if(increment == 2) then
                increment = 1
            else
                increment = increment * 5 / 11
            end if     
        end do
    end subroutine shell_Sort
    !========================================================================================
end module start_variables
