module input
    
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
    !       Loads the input parameters from input.rts
    !          for more details see RTS_global.f90
    !----------------------------------------------------------
    subroutine input_data
        implicit none
        Character(len = 60) :: aux
        integer:: i,error
        logical :: aux_flag
        open(unit=60,file='input/input.rts',form='formatted',status='old',&
             access = 'sequential',action='read',iostat=error)
        if(error /= 0) then
            print*, "error during opening of input.rts"
            stop
        end if
        !=================== RTS Logo ===================
        do i=1,11
            read(60,*)aux
        end do
        !========== Mesh Generation Parameters ==========
        read(60,*)DIMEN
        read(60,*)nx
        read(60,*)ny
        read(60,*)nz
        read(60,*)lx
        read(60,*)ly
        read(60,*)lz
        read(60,*)path_flag
        read(60,*)aux
        !===== Radiation Method for Solving the RTE =====
        read(60,*)rad_model
        read(60,*)schm_model
        read(60,*)aux
        !---------- Scattering Phase Function -----------
        read(60,*)aniso_flag
        read(60,*)scat_model
        read(60,*)SPF_func
        read(60,*)C_rad
        read(60,*)aux
        !============= Finite Angle Method  =============
        read(60,*)np
        read(60,*)nt
        read(60,*)nsub
        read(60,*)aux
        !========== Discrete Ordinates Method  ==========
        read(60,*)quadrature_DOM
        read(60,*)N_quad
        read(60,*)aux
        !======= Radiative Gas Properties Module ========
        read(60,*)gas_prop
        read(60,*)absrp_model
        read(60,*)nongray_flag
        read(60,*)aux; read(60,*)aux; read(60,*)aux; read(60,*)aux
        !=============== Input Parameters ===============!
        !---------- Initialization Parameters -----------!
        read(60,*)field_flag(1), T_energyc  !Temperature
        read(60,*)field_flag(2), Gc         !Incident radiation
        read(60,*)field_flag(5), k_rad      !Absorption coefficient
        read(60,*)field_flag(6), sig_rad    !Scattering coefficient
        read(60,*)field_flag(12),XCO2_g     !CO2 mole fraction
        read(60,*)field_flag(13),XH2O_g     !H2O mole fraction
        read(60,*)field_flag(11),P_g        !Pressure of the mixture
        read(60,*)aux; read(60,*)aux; read(60,*)aux
        !=========== Information at the Walls ===========
        SBCwall = .true. !by default all walls are solid walls
        do i=1,6
            read(60,*)wall_temp(i),epsilon_w(i),nuffwall(i),aux_flag
            if(aux_flag) SBCwall(i) = .false.
        end do
        read(60,*)aux; read(60,*)aux
        !============= Convergence Criteria =============
        read(60,*)rad_tol   
        read(60,*)ITMAX
        read(60,*)aux
        close(60)        
        if(path_flag) call input_paths  !load the file paths.rts
        call input_energy !load the file energy.rts
        call recorded_variables
        call internal_input
        call the_sanity_checks
    end subroutine input_data
    !========================================================================================
    !----------------------------------------------------------
    !       Loads the input paths parameters from paths.rts
    !          for more details see RTS_global.f90
    !----------------------------------------------------------
    subroutine input_paths
        !---------------------- Description -----------------------
        implicit none
        Character(len = 60) :: aux
        integer:: i,j,error
        open(unit=60,file='input/paths.rts',form='formatted',status='old',&
             access = 'sequential',action='read',iostat=error)
        if(error /= 0) then
            print*, "error during opening of paths.rts"
            stop
        end if
        !=================== RTS Logo ===================
        do i=1,11
            read(60,*)aux
        end do
        !================== Pathfinder ==================
        read(60,*)npath; allocate(Lpath(7,npath))
        read(60,*)aux
        read(60,*)aux
        do j=1,npath
            read(60,*)(Lpath(i,j),i=1,7) !Path i
        end do
        close(60)
    end subroutine input_paths
    !========================================================================================
    !----------------------------------------------------------
    !       Loads the input parameters from energy.rts
    !          for more details see RTS_global.f90
    !----------------------------------------------------------
    subroutine input_energy
        implicit none
        Character(len = 60) :: aux
        integer:: i,error
        open(unit=60,file='input/energy.rts',form='formatted',status='old',&
             access = 'sequential',action='read',iostat=error)
        if(error /= 0) then
            print*, "error during opening of energy.rts"
            stop
        end if
        !=================== RTS Logo ===================
        do i=1,11
            read(60,*)aux
        end do
        !=========== Thermal Input Parameters ===========
        read(60,*)energy_flag
        read(60,*)rte_flag
        read(60,*)trans_flag
        read(60,*)T_inf
        read(60,*)aux; read(60,*)aux; read(60,*)aux; read(60,*)aux
        !=============== Input Parameters ===============!
        read(60,*)field_flag(11), rho
        read(60,*)field_flag(12), Cp
        read(60,*)field_flag(8) , K_termc  !Thermal conductivity (W/m*K)
        if(energy_flag) then
            read(60,*)aux; read(60,*)aux; read(60,*)aux
            !=========== Information at the Walls ===========
            do i=1,6
                read(60,*)wall_temp(i), h_conv(i)
            end do
            read(60,*)aux; read(60,*)aux; read(60,*)aux
            !============= Boundary Conditions  =============
            do i=1,6
                read(60,*)BC_flags_t(i)
            end do
            read(60,*)aux
            !============= Convergence Criteria =============
            read(60,*)ene_tol
            read(60,*)aux
            !=============== Saving Results =================
            read(60,*)f_time 
            read(60,*)CTS
            read(60,*)aux
        end if
        close(60)
    end subroutine input_energy
    !========================================================================================
    !----------------------------------------------------------
    !       Loads the input parameters from recorded_var.rts
    !          for more details see RTS_global.f90
    !---------------------- Description -----------------------
    ! ID 1  - (T_energy) temperature for all nodes (K)
    ! ID 2  - (G) incident radiation for all nodes (w/m^2)
    ! ID 3  - (S_rad) divergent radiative flux (w/m^3)
    ! ID 4  - (IBlack) blackbody radiation intensity
    ! ID 5  - (cappa) absorption coefficient for all nodes (m^-1)
    ! ID 6  - (sigma) scattering coefficient for all nodes (m^-1)
    ! ID 7  - (beta) extinction coefficient for all nodes (m^-1)
    ! ID 8  - (K_term) thermal conductivity (W/m*K)
    ! ID 9  - (epsilon_rad) radiation emissivities
    !----------------------------------------------------------
    subroutine recorded_variables
        implicit none
        Character(len = 100) :: aux
        integer:: i,j,error
        open(unit=60,file='input/output.rts',form='formatted',status='old',&
             access = 'sequential',action='read',iostat=error)
        if(error /= 0) then
            print*, "error during opening of recorded_var.rts"
            stop
        end if
        !=================== RTS Logo ===================
        do i=1,13
            read(60,*)aux
        end do
        !=============== Saving Results =================
        read(60,*)vtk_flag
        read(60,*)dat_flag
        read(60,*)xyz_save
        read(60,*)cstm_save
        read(60,*)nuss_flag
        read(60,*)aux;read(60,*)aux;read(60,*)aux
        !============ Recorded Variables List ===========
        do i=1,9
            read(60,*)Rec_flag(i) !ID i
        end do
        read(60,*)aux;read(60,*)aux;read(60,*)aux;read(60,*)aux
        !========== Walls Custom Variable Save ==========
        ! j  - Variable (see output.rts)
        ! 1  - (T_energy) temperature for all nodes (K)
        ! 2  - (G) incident radiation for all nodes (w/m^2)
        ! 3  - (IBlack) blackbody radiation intensity
        ! 4  - Radiative walls fluxes
        ! 5  - Emissivities at the walls
        !------------------------------------------------
        ! i  - faces 1-6 (West,South,North,Botto,Top)
        !------------------------------------------------
        do j=1,5
            read(60,*)(crec_flag(i,j),i=1,6)
        end do
        read(60,*)aux;read(60,*)aux;read(60,*)aux
        !================== Slices Save ===================
        read(60,*)slice_flag
        if(slice_flag)then
            read(60,*)OPslice
            read(60,*)nslices
            read(60,*)aux;read(60,*)aux;
            allocate(slc_point(nslices),slc_perc(nslices),&
            slc_node(nslices),slc_ID(nslices),slc_axis(nslices))
            do j=1,nslices
                read(60,*)slc_axis(j),slc_perc(j),slc_point(j),slc_node(j),slc_ID(j)
            end do
        end if
        !==================================================
        do j=1,100
            read(60,*)aux
            if(aux =='!<>') exit
        end do
        read(60,*)aux;read(60,*)aux;read(60,*)aux
        !================== Lineup Save ===================
        read(60,*)lnp_flag
        if(lnp_flag)then
            read(60,*)OPlnp
            read(60,*)nlnp
            read(60,*)aux;read(60,*)aux;read(60,*)aux
            read(60,*)aux;read(60,*)aux;read(60,*)aux
            allocate(lnp_axis(nlnp),lnp_point(nlnp,2),&
                     lnp_node(nlnp,2),lnp_ID(nlnp))
            do j=1,nlnp
                read(60,*)lnp_axis(j),lnp_point(j,1),lnp_point(j,2),lnp_node(j,1),lnp_node(j,2),lnp_ID(j)
            end do
        end if
        close(60)
    end subroutine recorded_variables
    !========================================================================================
    !----------------------------------------------------------
    !               Handles the internal inputs
    !----------------------------------------------------------
    subroutine internal_input
        integer :: i
        !---------------- Radiation module ----------------!
        if(rad_model == 'FAM' .or. rad_model == 'fam')then
            model_stc = 3
        else if(rad_model == 'DOM' .or. rad_model == 'dom')then
            model_stc = 2
        else if(rad_model == 'P1' .or. rad_model == 'p1')then
            model_stc = 1
        end if
        !----------- Spatial Scheme for the RTE -----------!
        if(schm_model == 'Upwind' .or. schm_model == 'upwind' .or. schm_model == 'UPWIND' .or. schm_model == 'FOU')then
            scheme_flag = 1
        else if(schm_model == 'CDS' .or. schm_model == 'cds')then
            scheme_flag = 2
        end if
        !---------------- Scattering Model ----------------!
        if(scat_model == 'Mie' .or. scat_model == 'MIE' .or. scat_model == 'mie')then
            scat_flag = 1
        else if(scat_model == 'Rayleigh' .or. scat_model == 'RAYLEIGH' .or. scat_model == 'rayleigh')then
            scat_flag = 2
        else if(scat_model == 'Linear' .or. scat_model == 'LINEAR' .or. scat_model == 'linear')then
            scat_flag = 3
        else if(scat_model == 'HG' .or. scat_model == 'hg')then
            scat_flag = 4
        end if
        !--------------- Scattering Function --------------!
        if(SPF_func == 'F0' .or. SPF_func == 'f0')then
            SPF_flag = 1
        else if(SPF_func == 'F1' .or. SPF_func == 'f1')then
            SPF_flag = 2
        else if(SPF_func == 'F2' .or. SPF_func == 'f2')then
            SPF_flag = 3
        else if(SPF_func == 'F3' .or. SPF_func == 'f3')then
            SPF_flag = 4
        else if(SPF_func == 'F4' .or. SPF_func == 'f4')then
            SPF_flag = 5
        else if(SPF_func == 'F5' .or. SPF_func == 'f5')then
            SPF_flag = 6
        else if(SPF_func == 'F6' .or. SPF_func == 'f6')then
            SPF_flag = 7
        else if(SPF_func == 'F7' .or. SPF_func == 'f7')then
            SPF_flag = 8
        else if(SPF_func == 'B0' .or. SPF_func == 'b0')then
            SPF_flag = 9
        else if(SPF_func == 'B1' .or. SPF_func == 'b1')then
            SPF_flag = 10
        else if(SPF_func == 'B2' .or. SPF_func == 'b2')then
            SPF_flag = 11
        else if(SPF_func == 'B3' .or. SPF_func == 'b3')then
            SPF_flag = 12
        end if
        !---------------- Quadrature Mode -----------------!
        if(quadrature_DOM == 'Sn' .or. quadrature_DOM == 'SN' .or. quadrature_DOM == 'sn')then
            QM_stc = 1
        else if(quadrature_DOM == 'Tn' .or. quadrature_DOM == 'TN' .or. quadrature_DOM == 'tn')then
            QM_stc = 2
        else if(quadrature_DOM == 'Qn' .or. quadrature_DOM == 'QN' .or. quadrature_DOM == 'qn')then
            QM_stc = 3
        end if
        !---------------- Absorption Model ----------------!
        if(absrp_model == 'WSGG' .or. absrp_model == 'wsgg')then
            wsgg_model = .true.
            absrp_stc = 2
        else if(absrp_model == 'GRAY' .or. absrp_model == 'Gray' .or. absrp_model == 'gray')then
            gray_model = .true.
            absrp_stc = 1
        end if
        !---- Radiative Properties of Combustion Gases ----!
        const_conc = .true.
        do i = 11,13
            if(field_flag(i) .eqv. .false.) const_conc = .false.
        end do
        
    end subroutine internal_input
    !========================================================================================
    !----------------------------------------------------------
    !  Checks if the loaded variables are physically coherent
    !----------------------------------------------------------
    subroutine the_sanity_checks
        Character(len = 100) :: errormsg
        logical :: aux_flag
        integer :: i,j,fl = 108
        open(unit=fl,file='output/RTS.log')
        write(fl,'(a)')'---------------------------------------------------------------------------------------'
        write(fl,'(a)')'                                                                                       '
        write(fl,'(a)')'              /\\\\\\\\\      /\\\\\\\\\\\\\\\     /\\\\\\\\\\\                        '
        write(fl,'(a)')'             /\\\///////\\\   \///////\\\/////    /\\\/////////\\\                     '
        write(fl,'(a)')'             \/\\\     \/\\\         \/\\\        \//\\\      \///                     '
        write(fl,'(a)')'              \/\\\\\\\\\\\/          \/\\\         \////\\\                           '
        write(fl,'(a)')'               \/\\\//////\\\          \/\\\            \////\\\                       '
        write(fl,'(a)')'                \/\\\    \//\\\         \/\\\               \////\\\                   '
        write(fl,'(a)')'                 \/\\\     \//\\\        \/\\\        /\\\      \//\\\                 '
        write(fl,'(a)')'                  \/\\\      \//\\\       \/\\\       \///\\\\\\\\\\\/                 '
        write(fl,'(a)')'                   \///        \///        \///          \///////////                  '
        write(fl,'(a)')'                         Radiative Transfer Simulator Log File                         '
        write(fl,'(a)')'---------------------------------------------------------------------------------------'
        !================= Mesh Generation Parameters =================
        if(DIMEN <= 0 .or. DIMEN > 3) then
            errormsg = 'ERROR: Illegal Number of Dimensions!'
            write(*,'(a)') errormsg; write(fl,'(a)')errormsg
            stop_flag = .true.
        else
            if(nx <= 0) then
                errormsg = 'ERROR: Control volumes in x-direction Cannot be Null!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            if(ny <= 0 .and. DIMEN >= 2) then
                errormsg = 'ERROR: Control volumes in y-direction Cannot be Null!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            if(nz <= 0 .and. DIMEN == 3) then
                errormsg = 'ERROR: Control volumes in z-direction Cannot be Null!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            if(lx <= 0.0) then
                errormsg = 'ERROR: Domain Length Cannot be Null!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            if(ly <= 0.0 .and. DIMEN >= 2) then
                errormsg = 'ERROR: Domain Height Cannot be Null!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            if(lz <= 0.0 .and. DIMEN == 3) then
                errormsg = 'ERROR: Domain Width Cannot be Null!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
        end if
        !================== Thermal Input Parameters ==================
        if(energy_flag) then
            if(field_flag(11) .and. rho <= 0.0) then
                errormsg = 'ERROR: Fluid Density Cannot be Null!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            if(field_flag(12) .and. Cp <= 0.0) then
                errormsg = 'ERROR: Specific Heat Cannot be Null!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            if(field_flag(8) .and. K_termc <= 0.0) then               
                errormsg = 'ERROR: Thermal Conductivity Cannot be Null!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            if(T_inf < 0.0) then
                errormsg = 'ERROR: Room Temperature Cannot be Negative!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            !------------------ Boundary Conditions -------------------
            do i=1,6
                if(BC_flags_t(i) <= 0 .or. BC_flags_t(i) > 3 ) then
                    errormsg = 'ERROR: Illegal Type of Boundary Condition Selected'
                    write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                    stop_flag = .true.
                end if
            end do
        end if
        !================= Radiation Input Parameters =================
        if(field_flag(1) .and. T_energyc < 0.0) then
            errormsg = 'ERROR: Temperature Field Cannot be Negative!'
            write(*,'(a)') errormsg; write(fl,'(a)')errormsg
            stop_flag = .true.
        end if
        if(field_flag(2) .and. Gc < 0.0) then
            errormsg = 'ERROR: Incident Radiation Cannot be Negative!'
            write(*,'(a)') errormsg; write(fl,'(a)')errormsg
            stop_flag = .true.
        end if
        if(field_flag(5) .and. k_rad < 0.0) then
            errormsg = 'ERROR: Absorption Coefficient Cannot be Negative!'
            write(*,'(a)') errormsg; write(fl,'(a)')errormsg
            stop_flag = .true.
        end if
        if(field_flag(6) .and. sig_rad < 0.0) then
            errormsg = 'ERROR: Scattering Coefficient Cannot be Negative!'
            write(*,'(a)') errormsg; write(fl,'(a)')errormsg
            stop_flag = .true.
        end if
        if(rad_tol <= 0.0) then
            errormsg = 'ERROR: Radiation Convergence Criteria Cannot be Null!'
            write(*,'(a)') errormsg; write(fl,'(a)')errormsg
            stop_flag = .true.
        end if
        if(aniso_flag) then
            if(scat_flag <= 0 .or. scat_flag > 4) then
                errormsg = 'ERROR: Invalid Scattering Module!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            if(scat_flag == 1) then
                if(SPF_flag <= 0 .or. SPF_flag > 12) then
                    errormsg = 'ERROR: Invalid Mie Scattering Function!'
                    write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                    stop_flag = .true.
                end if
            end if
            if(scat_flag == 3 .or. model_stc == 1) then
                if(C_rad < -1.0d0 .or. C_rad > 1.0d0) then
                    errormsg = 'ERROR: Anisotropic Phase Function Coefficient Cannot Less than -1.0 and Greater than 1.0!'
                    write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                    stop_flag = .true.
                end if
            end if
        end if
        !================== Information at the Walls ==================
        do i=1,6
            if(nuffwall(i)) then
                if(wall_temp(i) < 0.0) then
                    errormsg = 'ERROR: Wall Temperature Cannot be Negative!'
                    write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                    stop_flag = .true.
                end if
                if(epsilon_w(i) < 0.0) then
                    errormsg = 'ERROR: Emissivity Cannot be Negative!'
                    write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                    stop_flag = .true.
                end if
                if(epsilon_w(i) > 1.0) then
                    errormsg = 'ERROR: Emissivity Cannot be Greater Than 1.0!'
                    write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                    stop_flag = .true.
                end if
            end if
        end do
        !======================= Saving Results =======================
        if(vtk_flag.eqv..false.) then
            errormsg = 'WARNING: .vtk Results Will not be Saved!'
            write(*,'(a)') errormsg; write(fl,'(a)')errormsg
        end if
        if(dat_flag.eqv..false.) then
            errormsg = 'WARNING: .dat Results Will not be Saved!'
            write(*,'(a)') errormsg; write(fl,'(a)')errormsg
        end if
        if(trans_flag) then
            if(f_time <= 0) then
                errormsg = 'ERROR: Final Simulation Time Cannot be Null!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            if(CTS <= 0) then
                errormsg = 'ERROR: Illegal save Data Interval!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
        end if
        !============ Radiation Method for Solving the RTE ============
        if(model_stc <= 0 .or. model_stc > 3) then
            errormsg = 'ERROR: Invalid Radiation Module!'
            write(*,'(a)') errormsg; write(fl,'(a)')errormsg
            stop_flag = .true.
        end if
        !Discrete Ordinates Method
        if(model_stc == 2) then
            if(QM_stc <= 0 .or. QM_stc > 3) then
                errormsg = 'ERROR: Invalid Quadrature Model Selected!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            if(QM_stc == 1)then
                if(mod(N_quad,2)/=0) then
                    errormsg = 'ERROR: The Quadrature Order in The Sn Method Must be a Multiple of Two!' 
                    write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                    stop_flag = .true.
                end if
                if(N_quad > 10) then                       
                    errormsg = 'ERROR: The Quadrature Order in The Sn Method at the Moment it Cannot be Greater than 10!'
                    write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                    stop_flag = .true.
                end if
            end if
            if(N_quad == 0) then
                errormsg = 'ERROR: The Quadrature Order Cannot be Null!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
        !Finite Angle Method
        else if(model_stc == 3) then
            if(nt <= 0) then
                errormsg = 'ERROR: Control Volumes in theta-direction Cannot be Null!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            if(np <= 0) then
                errormsg = 'ERROR: Control Volumes in phi-direction Cannot be Null!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            if(mod(nt,2)/=0) then
                errormsg = 'ERROR: The Number of Control Volumes in theta-direction Must be a Multiple of Two!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            if(mod(np,4)/=0) then
                errormsg = 'ERROR: The Number of Control Volumes in phi-direction Must be a Multiple of Four!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
        end if
        if(scheme_flag <= 0 .or. scheme_flag > 2) then
            errormsg = 'ERROR: Invalid Spatial Scheme!'
            write(*,'(a)') errormsg; write(fl,'(a)')errormsg
            stop_flag = .true.
        end if
        !========================= Pathfinder =========================
        if(path_flag) then
            if(npath == 0) then
                errormsg = 'WARNING: The Number of Paths is Null!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
            else if(npath < 0) then
                errormsg = 'ERROR: Illegal Number of Paths!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            else if(npath > 10) then
                errormsg = 'ERROR: Too Many Paths!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            !Path Coordinates
            do i=1,npath
                if(Lpath(1,i) < 0.0 .or. Lpath(3,i) < 0.0 .or. Lpath(5,i) < 0.0) then
                    write(*,'(a50,i1)') 'ERROR: Paths Cannot Begin Before the Domain! Path ', i
                    write(fl,'(a50,i1)')'ERROR: Paths Cannot Begin Before the Domain! Path ', i
                    stop_flag = .true.
                end if
                if(Lpath(2,i) > lx .or. Lpath(4,i) > ly .or. Lpath(6,i) > lz) then
                    write(*,'(a52,i1)') 'ERROR: Paths Cannot be Bigger Than the Domain! Path ', i
                    write(fl,'(a52,i1)')'ERROR: Paths Cannot be Bigger Than the Domain! Path ', i
                    stop_flag = .true.
                end if
                if(mod(Lpath(7,i),2.0)/=0) then
                    write(*,'(a61,i1)') 'ERROR: The Refinement Level Must be a Multiple of Two! Path ', i
                    write(fl,'(a61,i1)')'ERROR: The Refinement Level Must be a Multiple of Two! Path ', i
                    stop_flag = .true.
                end if
                if(Lpath(7,i)>16) then
                    write(*,'(a42,i1)') 'ERROR: The Path Level Its too High! Path ', i
                    write(fl,'(a42,i1)')'ERROR: The Path Level Its too High! Path ', i
                    stop_flag = .true.
                end if
            end do
        end if
        !========= Radiative Properties of Combustion Gases ===========
        if(gas_prop) then
            if(absrp_stc <= 0 .or. absrp_stc > 2) then
                errormsg = 'ERROR: Invalid Gas Properties Module!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            if(const_conc) then
                if(XCO2_g <= 0.0d0 .and. XH2O_g <= 0.0d0) then
                    errormsg = 'ERROR: The CO2 and H2O Molar Fractions Cannot be Zero at The Same Time!'
                    write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                    stop_flag = .true.
                end if
                if(P_g < 0.0d0) then
                    errormsg = 'ERROR: The Total pressure of the mixture Cannot be Negative!'
                    write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                    stop_flag = .true.
                end if
            end if
        else
            if(nongray_flag) then
                errormsg = 'The non-gray formulation cannot be used with the gas properties module deactivated!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
        end if
        !=================== Recorded Variables List ==================
        do i=1,8
            if(Rec_flag(i))then
                aux_flag = .false.
                exit
            else
                aux_flag = .true.
            end if
        end do
        if(aux_flag)then
            errormsg = 'WARNING: No variables will be saved! Go to input\recorded_var.rts'
            write(*,'(a)') errormsg; write(fl,'(a)')errormsg
        end if
        !================= Walls Custom Variable Save =================
        if(cstm_save)then
            do j=1,5
                do i=1,2
                    if(crec_flag(i,j))then !ZY - West,East
                        plane_flag(1) = .true.
                        exit
                    end if
                end do
                do i=3,4
                    if(crec_flag(i,j))then !XZ - South,North
                        plane_flag(2) = .true.
                        exit
                    end if
                end do 
                do i=5,6
                    if(crec_flag(i,j))then !XY - Bottom,Top
                        plane_flag(3) = .true.
                        exit
                    end if
                end do 
            end do
        end if
        !========================= Slices Save ========================
        if(slice_flag)then
            if(OPslice <= 0 .or. OPslice > 3)then
                errormsg = 'ERROR: Invalid Slicing Method!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            if(nslices <= 0)then
                errormsg = 'ERROR: The Number of Slices Cannot be Negative or Null!'
                write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                stop_flag = .true.
            end if
            do j=1,nslices
                if(slc_perc(j) < 0 .and. slc_perc(j) >100)then
                    errormsg = 'ERROR: The percentage must be between 0 and 100!'
                    write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                    stop_flag = .true.
                end if
                if(slc_ID(j) < 1 .and. slc_ID(j) > 9)then
                    errormsg = 'ERROR: The unrecognized variable ID!'
                    write(*,'(a)') errormsg; write(fl,'(a)')errormsg
                    stop_flag = .true.
                end if
                !read(60,*)slc_axis(j),slc_point(j),slc_node(j)
            end do
        end if   
        !==============================================================
        if(stop_flag) then
            errormsg = '------- EXECUTION TERMINATED! -------'
            write(*,'(a)') errormsg; write(fl,'(a)')errormsg
            stop
        else
            write(fl,'(a)')'No potential errors found :)'
        end if
        close(fl)
    end subroutine the_sanity_checks
    !========================================================================================
end module input
