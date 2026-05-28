module output
    
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
    
    interface
        !-----------------------------------------------------------------
        subroutine event(variable,name,face,FID)
            implicit none
            character(len=*), intent (in):: name
            integer, intent (in) ::face,FID
            double precision, intent (in) :: variable(:,:,:)
        end subroutine event
        !-----------------------------------------------------------------
    end interface
    procedure(event), pointer :: sub_pointer
    
    contains
    !========================================================================================
    !----------------------------------------------------------
    !             Main routine for saving results
    !---------------------- Description -----------------------
    ! <EXTERNAL ROUTINES>
    ! VTK_save     - Saves the selected variables in a vtk file
    ! dat_save     - Saves the selected variables in a dat file
    ! main_slice   - Main slices save routine
    ! Nusselt_main - Main routine for calculating the Nusselt number
    ! walls_custom_save - Saves information on walls on different faces
    ! main_lineup  - Main lineup save routine
    !----------------------------------------------------------
    subroutine main_save()
        if(cstm_save .and. DIMEN == 3) call walls_custom_save
        if(csv_flag) call csv_save
        if(vtk_flag) call VTK_save
        if(dat_flag) call dat_save
        if(slice_flag .and. DIMEN == 3) call main_slice
        if(lnp_flag .and. DIMEN > 1) call main_lineup
        if(nuss_flag) call Nusselt_main
    end subroutine main_save

    !========================================================================================
    !----------------------------------------------------------
    !    Saves the selected variables in a CSV file
    !---------------------- Description -----------------------
    !----------------------------------------------------------
    subroutine csv_save
        implicit none
        integer :: i, j, k
        character(len = 1024)::name

        if(trans_flag) then
            write(name,'(A,I6.6,A)') 'output/RTSresults_',ITP,'.csv'
        else
            write(name,'(A,I6.6,A)') 'output/RTSresult.csv'
        end if

        open(unit=778, file = trim(name), action="write", status="replace", access = "sequential")
            write(778,'(a)') adjustl(trim(get_header()))

            do k = 1,nzi
                do j = 1,nyi
                    do i = 1,nxi
                        write(778,'(a)') adjustl(trim(make_row()))
                    end do
                end do
            end do
        close(778)
    contains

    function get_header() result(header)
        implicit none
        character(len = 2048) :: header

        header = "i,j,k,x,y,z"
        if(Rec_flag(1)) header = trim(adjustl(header)) // ',T_energy'
        if(Rec_flag(2)) header = trim(adjustl(header)) // ',G'
        if(Rec_flag(3)) header = trim(adjustl(header)) // ',S_rad'
        if(Rec_flag(4)) header = trim(adjustl(header)) // ',Ib'
        if(Rec_flag(5)) header = trim(adjustl(header)) // ',cappa'
        if(Rec_flag(6)) header = trim(adjustl(header)) // ',sigma_rad'
        if(Rec_flag(7)) header = trim(adjustl(header)) // ',beta'
        if(Rec_flag(8)) header = trim(adjustl(header)) // ',K_term'
        if(Rec_flag(9)) header = trim(adjustl(header)) // ',Emissivity'
        if(Rec_flag(10)) header = trim(adjustl(header)) // ',q_radw'
        if(Rec_flag(11)) header = trim(adjustl(header)) // ',P'
        if(Rec_flag(12)) header = trim(adjustl(header)) // ',XCO2'
        if(Rec_flag(13)) header = trim(adjustl(header)) // ',XH2O'
        header = trim(adjustl(header))
    end function get_header

    function make_row() result(row)
        implicit none
        character(len = 2048) :: row
        character(len = 32) :: value_buffer
        integer :: rf

        write(row,'(*(G30.20,:,","))') i, j, k, xc(i), yc(j), zc(k)

        do rf = 1,13
            if(Rec_flag(rf)) then
                write(value_buffer, "(G30.20)") get_value(rf)
                row = trim(adjustl(row)) // "," // trim(adjustl(value_buffer))
            end if
        end do
    end function make_row

    function get_value(rf) result(value)
        implicit none
        integer, intent(in) :: rf
        double precision :: value

        select case (rf)
            case (1) 
                value = T_energy(i,j,k)
            case (2)
                value = G(i,j,k)
            case (3) 
                value = S_rad(i,j,k)
            case (4) 
                value = IBlack(i,j,k)
            case (5) 
                value = cappa(i,j,k)
            case (6) 
                value = sigma(i,j,k)
            case (7) 
                value = beta(i,j,k)
            case (8) 
                value = K_term(i,j,k)
            case (9) 
                value = epsilon_rad(i,j,k)
            case (10)
                value = Q_radw(i,j,k)
            case (11)
                value = P_species(i,j,k)
            case (12)
                value = Y_species(i,j,k,1)
            case (13)
                value = Y_species(i,j,k,2)
            case default
                value = -1.0d0
        end select
    end function get_value
            
    end subroutine csv_save

    !========================================================================================
    !----------------------------------------------------------
    !       Saves the selected variables in a vtk file
    !---------------------- Description -----------------------
    ! <INPUT>
    ! x       - x-direction grid
    ! y       - y-direction grid
    ! z       - z-direction grid
    !----------------------------------------------------------
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
    subroutine VTK_save()
        integer :: i,j,k,fl = 15
        character(len=60):: name
        if(trans_flag) then
            write(name,'(A,I6.6,A)') 'output/RTSresults_',ITP,'.vtk'
        else
            write(name,'(A,I6.6,A)') 'output/RTSresult.vtk'
        end if
        open(unit=fl, file = name, action = "write", status = "replace")
        write(fl,'(a)')'# vtk DataFile Version 3.0'
        write(fl,'(a)') 'Non-uniform Rectilinear - Rectilinear Grid'
        write(fl,'(a)') 'ASCII'
        write(fl,'(a)') 'DATASET RECTILINEAR_GRID'
        if(trans_flag) then
            write(fl,'(a)') 'FIELD FieldData 2'
            write(fl,'(a)') 'TIME 1 1 double'
            write(fl,'(f14.6)') time
            write(fl,'(a)') 'CYCLE 1 1 int'
            write(fl,'(3i4)') ITP
        end if
        write(fl,'(a,3i4)')'DIMENSIONS', nxp,nyp,nzp
        write(fl,'(a13,1i4,a6)') 'X_COORDINATES',nxp,' float'
        write(fl,'(f14.6)') (x(i),i=1,nxp)
        write(fl,'(a13,1i4,a6)') 'Y_COORDINATES',nyp,' float'
        write(fl,'(f14.6)') (y(j),j=1,nyp)
        write(fl,'(a13,1i4,a6)') 'Z_COORDINATES',nzp,' float'
        write(fl,'(f14.6)') (z(k),k=1,nzp)
        write(fl, '(a, I16)') 'CELL_DATA ', nx*ny*nz
        !-------------------------------------------------------------
        if(Rec_flag(1)) call varsave3D(T_energy,'T_energy',fl)
        if(Rec_flag(2)) call varsave3D(G,'G',fl)
        if(Rec_flag(3)) call varsave3D(S_rad,'S_rad',fl)
        if(Rec_flag(4)) call varsave3D(IBlack,'Ib',fl)
        if(Rec_flag(5)) call varsave3D(cappa,'cappa',fl)
        if(Rec_flag(6)) call varsave3D(sigma,'sigma',fl)
        if(Rec_flag(7)) call varsave3D(beta,'beta',fl)
        if(Rec_flag(8)) call varsave3D(K_term,'K_term',fl)
        if(Rec_flag(9)) call varsave3D(epsilon_rad,'emissivities',fl)
        !-------------------------------------------------------------
        close(fl)
        if(trans_flag .eqv. .false.) write(*,'(a)')'RTSresult.vtk Have Been Successfully Saved in ./output'
    end subroutine VTK_save
    !========================================================================================
    !----------------------------------------------------------
    !       Saves the selected variables in a vtk file
    !---------------------- Description -----------------------
    ! <INPUT>
    ! variable - Variable to be saved
    ! name     - variable name to be displayed
    ! ID       - File ID
    !----------------------------------------------------------
    subroutine varsave3D(variable,name,ID)
        character(len=*):: name
        integer :: i,j,k,ID
        double precision, intent (in) :: variable(:,:,:)
        write(ID,'(a,a,a)')'SCALARS ',trim(name),' float'
        write(ID,'(a)') 'LOOKUP_TABLE default'      
        write(ID,'(f30.14)') (((variable(i,j,k),i=2,nxp),j=2,nyp),k=2,nzp)
    end subroutine varsave3D
    !=======================================================================================
    !----------------------------------------------------------
    !       Saves the selected variables in a dat file
    !---------------------- Description -----------------------
    ! <INPUT>
    ! x       - x-direction grid
    ! y       - y-direction grid
    ! z       - z-direction grid
    ! xc      - x value at the control-volume center
    ! yc      - y value at the control-volume center
    ! zc      - z value at the control-volume center
    ! dxp     - Control-volume length
    ! dyp     - Control-volume height
    ! dzp     - Control-volume width 
    ! <EXTERNAL ROUTINES>
    ! mesh_save   - Saves the mesh position vectors
    ! mesh_dat    - Print the mesh position vectors in a dat file
    ! varsave_DAT - Print the selected variable in a dat file
    ! walls_custom_saveDAT - Saves information on walls on different faces
    !----------------------------------------------------------
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
    subroutine dat_save()
        integer :: k,fl = 4
        if(DIMEN == 2) k = 2           !Central Plane
        if(DIMEN == 3) k = int(nz/2)+1 !Central Plane
        !=============================== Header ===============================!
        open(unit=fl, file = 'output/RTSresult.dat', action = "write", status = "replace")
        write(fl,'(a)')'======================================================================================='
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
        write(fl,'(a)')'                         Radiative Transfer Simulator Results                          '
        write(fl,'(a)')'======================================================================================='
        write(fl,'(a)') '   Mesh Data'
        write(fl,'(a)')'======================================================================================='
        call mesh_dat(x,xc,dxp,nxi,'  i     X           Xc         dx',fl)
        if(DIMEN > 1)  call mesh_dat(y,yc,dyp,nyi,'  j     Y           Yc         dy',fl)
        if(DIMEN == 3) call mesh_dat(z,zc,dzp,nzi,'  k     Z           Zc         dz',fl)
        !========================= Information Fields =========================!
        if(DIMEN == 3) write(fl,'(a60,I2)') 'All recorded information below refer to the central plan k = ',k
        !-------------------------------------------------------------------------------------
        if(Rec_flag(1)) call varsave_DAT(T_energy,'   Medium temperature T_energy (K)',fl,k)
        if(Rec_flag(2)) call varsave_DAT(G,'   Incident Radiation G (W/m2)',fl,k)
        if(Rec_flag(3)) call varsave_DAT(S_rad,'   Divergent Radiative Flux (W/M3)',fl,k)
        if(Rec_flag(4)) call varsave_DAT(IBlack,'   Blackbody Radiation Intensity IBlack (W/m2)',fl,k)
        if(Rec_flag(5)) call varsave_DAT(cappa,'   Absorption Coefficient (1/m)',fl,k)
        if(Rec_flag(6)) call varsave_DAT(sigma,'   Scattering Coefficient (1/m)',fl,k)
        if(Rec_flag(7)) call varsave_DAT(beta,'   Extinction Coefficient (1/m)',fl,k)
        if(Rec_flag(8)) call varsave_DAT(K_term,'   Thermal conductivity (W/m*K)',fl,k)
        if(Rec_flag(9)) call varsave_DAT(epsilon_rad,'   Gas Emissivity ( )',fl,k)
        call walls_custom_saveDAT(fl) !Walls informations
        close(fl)
        write(*,'(a)')'RTSresult.dat Have Been Successfully Saved in ./output'
        if(xyz_save) call mesh_save
    end subroutine dat_save
    !========================================================================================
    !----------------------------------------------------------
    !        Print the selected variable in a dat file
    !---------------------- Description -----------------------
    ! <INPUT>
    ! variable - Variable to be saved
    ! name     - variable name to be displayed
    ! ID       - File ID
    ! k        - z direction plane
    ! cc       - X value at the control-volume center (for 1D)
    !----------------------------------------------------------
    subroutine varsave_DAT(variable,name,ID,k)

        character(len=*):: name
        integer :: i,j,k,ID
        double precision, intent (in) :: variable(:,:,:)
        call header_var_writeDAT(name,ID)
        if(DIMEN == 2 .or. DIMEN == 3)then
            do j=nyi,1,-1
                write(ID,*) (variable(i,j,k),i=1,nxi)
            end do
        else if(DIMEN == 1)then
            write(ID,'(a)')'========================='
            write(ID,'(a)')'     Xc       VARIABLE   '
            write(ID,'(a)')'========================='
            do i=1,nxi
                write(ID,'(ES10.3,F15.4)') xc(i),variable(i,2,2)
            end do
        end if
    end subroutine varsave_DAT
    !========================================================================================
    !----------------------------------------------------------
    !       Print the mesh position vectors in a dat file
    !---------------------- Description -----------------------
    ! <INPUT>
    ! S     - S-direction grid
    ! Sc    - S value at the control-volume center
    ! ds    - Control-volume length|height|width 
    ! ns    - S-direction grid size
    ! FID   - File ID
    ! Title - Title to be displayed
    !----------------------------------------------------------
    subroutine mesh_dat(S,Sc,ds,ns,Title,FID)

        integer :: si,ns,FID
        character(len=*):: Title
        double precision, intent (in), dimension(:) :: S,Sc,ds
        double precision :: Sn(ns)
        Sn(1:ns-1) = S(1:ns-1)
        Sn(ns) = 0.0d0
        write(FID,'(a)') '=====================================' 
        write(FID,'(a)')trim(Title)
        write(FID,'(a)') '====================================='  
        do si=1,ns
            write(FID,'(I3,ES10.3,ES12.3,ES12.3)') si,Sn(si),Sc(si),ds(si)
        end do
    end subroutine mesh_dat
    !=======================================================================================
    !----------------------------------------------------------
    !             Saves the mesh position vectors
    !---------------------- Description -----------------------
    ! <INPUT>
    ! x       - x-direction grid
    ! y       - y-direction grid
    ! z       - z-direction grid
    ! xc      - x value at the control-volume center
    ! yc      - y value at the control-volume center
    ! zc      - z value at the control-volume center
    !----------------------------------------------------------
    subroutine mesh_save()

        integer :: i,j,k,fl=8
        open(UNIT=fl,FILE='output/fmesh.dat')        
        !----------------- direction grids ----------------!
        write(fl,*) (x(i),i=1,nxp)     !x-direction
        if(DIMEN == 2 .or. DIMEN == 3)then
            write(fl,*) (y(j),j=1,nyp) !y-direction
        end if
        if(DIMEN == 3)then
            write(fl,*) (z(k),k=1,nzp) !z-direction
        end if
        close(fl)
        open(UNIT=fl,FILE='output/cmesh.dat') 
        !----------- grids at the center of C.V. ----------!
        write(fl,*) (xc(i),i=1,nxi)     !x-direction
        if(DIMEN == 2 .or. DIMEN == 3)then
            write(fl,*) (yc(j),j=1,nyi) !y-direction
        end if
        if(DIMEN == 3)then
            write(fl,*) (zc(k),k=1,nzi) !z-direction
        end if
        close(fl)
    end subroutine mesh_save
    !========================================================================================
    !----------------------------------------------------------
    !       Saves information on walls on different faces
    !---------------------- Description -----------------------
    ! <INPUT>
    ! x            - x-direction grid
    ! y            - y-direction grid
    ! z            - z-direction grid
    ! plane_flag   - Plane ID flags 1-3 (ZY, XZ, XY)
    ! <INTERNAL>
    ! FID          - File ID
    ! <EXTERNAL ROUTINES>
    ! file_open2D  - Open the file in which the plans will be saved
    ! file_write2D - Main routine to write the variables in the selected plans
    !----------------------------------------------------------
    subroutine walls_custom_save()
        integer :: FID = 8
        sub_pointer => varsave2D_VTK
        if(plane_flag(1))then !ZY - West,East
            call file_open2D(z,y,'zyplane',FID)
            call file_write2D(1,FID)
            write(*,'(a)')'RTS_zyplane.vtk Have Been Successfully Saved in ./output'
        end if
        if(plane_flag(2))then !XZ - South,North
            call file_open2D(x,z,'xzplane',FID)
            call file_write2D(2,FID)
            write(*,'(a)')'RTS_xzplane.vtk Have Been Successfully Saved in ./output'
        end if
        if(plane_flag(3))then !XY - Bottom,Top
            call file_open2D(x,y,'xyplane',FID)
            call file_write2D(3,FID)
            write(*,'(a)')'RTS_xyplane.vtk Have Been Successfully Saved in ./output'
        end if
        nullify(sub_pointer)
    end subroutine walls_custom_save
    !=======================================================================================
    !----------------------------------------------------------
    !      Open the file in which the plans will be saved
    !---------------------- Description -----------------------
    ! <INPUT>
    ! Sx     - Sx-direction grid
    ! Sy     - Sy-direction grid
    ! snxp   - Sx grid size
    ! snyp   - Sy grid size
    ! Fname  - Output file name
    ! FID    - File ID
    !----------------------------------------------------------
    subroutine file_open2D(Sx,Sy,Fname,FID)
        double precision, intent(in), dimension(:) :: Sx,Sy
        integer :: i,j,snxp,snyp,FID
        character(len=*):: Fname
        character(len=60):: name
        snxp = size(Sx); snyp = size(Sy)
        name = trim('output/RTS_') // trim(Fname)//trim('.vtk')
        open(unit=FID, file = name, action="write", status="replace")
        write(FID,'(a)')'# vtk DataFile Version 3.0'
        write(FID,'(a)') 'Non-uniform Rectilinear - Rectilinear Grid'
        write(FID,'(a)') 'ASCII'
        write(FID,'(a)') 'DATASET RECTILINEAR_GRID'
        write(FID,'(a,3i4)')'DIMENSIONS', snxp,snyp,1
        write(FID,'(a13,1i4,a6)') 'X_COORDINATES',snxp,' float'
        write(FID,'(f14.6)') (Sx(i),i=1,snxp)
        write(FID,'(a13,1i4,a6)') 'Y_COORDINATES',snyp,' float'
        write(FID,'(f14.6)') (Sy(j),j=1,snyp)
        write(FID,'(a13,1i4,a6)') 'Z_COORDINATES',1,' float'
        write(FID,'(f14.6)' ) 1.0d0
        write(FID,'(a, I16)') 'CELL_DATA ', (snxp - 1)*(snyp - 1)
    end subroutine file_open2D
    !=======================================================================================
    !----------------------------------------------------------
    ! Main routine to write the variables in the selected plans
    !---------------------- Description -----------------------
    ! <INPUT>
    ! wallID      - Selected face ID
    ! nameID      - Selected variable ID
    ! planeID     - Selected plan ID
    ! FID         - File ID
    ! crec_flag   - Walls custom variables save flags
    ! <EXTERNAL ROUTINES>
    ! var_print2D - Select the variable to be saved
    !----------------------------------------------------------
    subroutine file_write2D(planeID,FID)
        integer :: wallID,nameID,planeID,FID
        do nameID=1,5
            if(planeID == 1)then
                do wallID=1,2
                    if(crec_flag(wallID,nameID)) call var_print2D(nameID,wallID,FID)
                end do
            else if(planeID == 2)then
                do wallID=3,4
                    if(crec_flag(wallID,nameID)) call var_print2D(nameID,wallID,FID)
                end do
            else if(planeID == 3)then
                do wallID=5,6
                    if(crec_flag(wallID,nameID)) call var_print2D(nameID,wallID,FID)
                end do
            end if
        end do
        close(FID) !The file is closed here
    end subroutine file_write2D
    !=======================================================================================
    !----------------------------------------------------------
    !              Select the variable to be saved
    !---------------------- Description -----------------------
    ! <INPUT>
    ! wallID    - Selected face ID
    ! nameID    - Selected variable ID
    ! FID       - File ID
    ! <EXTERNAL ROUTINES>
    ! varsave2D_VTK - Write the selected variables in the selected plans
    !----------------------------------------------------------
    subroutine var_print2D(nameID,wallID,FID)
        integer :: wallID,nameID,FID
        if(nameID == 1)then
            call sub_pointer(T_energy,'T_energy',wallID,FID)
        else if(nameID == 2)then
            call sub_pointer(G,'G',wallID,FID)
        else if(nameID == 3)then
            call sub_pointer(IBlack,'Ib',wallID,FID)
        else if(nameID == 4)then
            call sub_pointer(Q_radw,'q_wr',wallID,FID)
        else if(nameID == 5)then
            call sub_pointer(epsilon_rad,'e_wr',wallID,FID)
        end if
    end subroutine var_print2D
    !=======================================================================================
    !----------------------------------------------------------
    !    Write the selected variables in the selected plans
    !---------------------- Description -----------------------
    ! <INPUT>
    ! face     - Selected face ID
    ! name     - Name of the variable to be displayed
    ! variable - Variable to be saved
    ! FID      - File ID
    ! <EXTERNAL ROUTINES>
    ! var_name_save - Modifies the name to be displayed according to the selected face
    !----------------------------------------------------------
    subroutine varsave2D_VTK(variable,name,face,FID)
        double precision, intent (in) :: variable(:,:,:)
        character(len=*), intent (in):: name
        integer, intent (in):: face,FID
        integer :: i,j,k
        write(FID,'(a,a,a)')'SCALARS ',var_name_save(name,face),' float'
        write(FID,'(a)') 'LOOKUP_TABLE default'      
        if(face == 1)then !WEST
            write(FID,'(f30.14)') ((variable(1,j,k),k=2,nzp),j=2,nyp)
        else if(face == 2)then !EAST
            write(FID,'(f30.14)') ((variable(nxi,j,k),k=2,nzp),j=2,nyp)
        else if(face == 3)then !SOUTH
            write(FID,'(f30.14)') ((variable(i,1,k),i=2,nxp),k=2,nzp)
        else if(face == 4)then !NORTH
            write(FID,'(f30.14)') ((variable(i,nyi,k),i=2,nxp),k=2,nzp)
        else if(face == 5)then !BOTTOM
            write(FID,'(f30.14)') ((variable(i,j,1),i=2,nxp),j=2,nyp)
        else if(face == 6)then !TOP
            write(FID,'(f30.14)') ((variable(i,j,nzi),i=2,nxp),j=2,nyp)
        end if
    end subroutine varsave2D_VTK
    !=======================================================================================
    !----------------------------------------------------------
    ! Modifies the name to be displayed according to the selected face
    !---------------------- Description -----------------------
    ! <INPUT>
    ! name     - Name of the variable to be displayed
    ! var_name_save - Modified name
    !----------------------------------------------------------
    function var_name_save(name,face)
        character(len=60):: var_name_save
        character(len=*):: name
        integer :: face
        if(face == 1)then !WEST
            var_name_save = trim(name)//trim('_WEST')
        else if(face == 2)then !EAST
            var_name_save = trim(name)//trim('_EAST')
        else if(face == 3)then !SOUTH
            var_name_save = trim(name)//trim('_SOUTH')
        else if(face == 4)then !NORTH
            var_name_save = trim(name)//trim('_NORTH')
        else if(face == 5)then !BOTTOM
            var_name_save = trim(name)//trim('_BOTTOM')
        else if(face == 6)then !TOP
            var_name_save = trim(name)//trim('_TOP')
        end if
    end function var_name_save
    !=======================================================================================
    !----------------------------------------------------------
    !       Print an identification bar in the dat file
    !---------------------- Description -----------------------
    ! <INPUT>
    ! FID   - File ID
    ! Title - Title to be displayed
    !----------------------------------------------------------
    subroutine header_var_writeDAT(Title,FID)
        integer :: FID
        character(len=*):: Title
        write(FID,'(a)')'======================================================================================='
        write(FID,'(a)')trim(Title)
        write(FID,'(a)')'======================================================================================='
    end subroutine header_var_writeDAT
    !=======================================================================================
    !----------------------------------------------------------
    !       Saves information on walls on different faces
    !---------------------- Description -----------------------
    ! <INPUT>
    ! wallID        - Selected face ID
    ! nameID        - Selected variable ID
    ! FID           - File ID
    ! crec_flag     - Walls custom variables save flags
    ! <EXTERNAL ROUTINES>
    ! var_print2D   - Select the variable to be saved
    ! varsave3D_DAT - Write the selected variables for 3D problems
    ! varsave2D_DAT - Write the selected variables for 2D problems
    ! varsave1D_DAT - Write the selected variables for 1D problems
    !----------------------------------------------------------
    subroutine walls_custom_saveDAT(FID)
        integer :: wallID,nameID,FID
        if(DIMEN == 3)then
            sub_pointer => varsave3D_DAT
            do nameID=1,5
                do wallID=1,6
                    if(crec_flag(wallID,nameID)) call var_print2D(nameID,wallID,FID)
                end do
            end do
        else if(DIMEN == 2)then
            sub_pointer => varsave2D_DAT
            do nameID=1,5
                do wallID=1,4
                    if(crec_flag(wallID,nameID)) call var_print2D(nameID,wallID,FID)
                end do
            end do
        else if(DIMEN == 1)then
            sub_pointer => varsave1D_DAT
            do nameID=1,5
                do wallID=1,2
                    if(crec_flag(wallID,nameID)) call var_print2D(nameID,wallID,FID)
                end do
            end do
        end if
        nullify(sub_pointer)
    end subroutine walls_custom_saveDAT
    !=======================================================================================
    !----------------------------------------------------------
    !       Write the selected variables for 3D problems
    !---------------------- Description -----------------------
    ! <INPUT>
    ! face     - Selected face ID
    ! name     - Name of the variable to be displayed
    ! variable - Variable to be saved
    ! fl       - File ID
    ! <EXTERNAL ROUTINES>
    ! header_var_writeDAT - Print an identification bar in the dat file
    !----------------------------------------------------------
    subroutine varsave3D_DAT(variable,name,face,fl)
        double precision, intent (in) :: variable(:,:,:)
        character(len=*), intent (in):: name
        integer, intent (in):: face,fl
        integer :: i,j,k
        call header_var_writeDAT(var_name_saveDAT(name,face),fl)
        if(face == 1)then !WEST
            do k=nzp,2,-1
                write(fl,*) (variable(1,j,k),j=2,nyp)
            end do
        else if(face == 2)then !EAST
            do k=nzp,2,-1
                write(fl,*) (variable(nxi,j,k),j=2,nyp)
            end do
        else if(face == 3)then !SOUTH
            do k=nzp,2,-1
                write(fl,*) (variable(i,1,k),i=2,nxp)
            end do
        else if(face == 4)then !NORTH
            do k=nzp,2,-1
                write(fl,*) (variable(i,nyi,k),i=2,nxp)
            end do
        else if(face == 5)then !BOTTOM
            do j=nyp,2,-1
                write(fl,*) (variable(i,j,1),i=2,nxp)
            end do
        else if(face == 6)then !TOP
            do j=nyp,2,-1
                write(fl,*) (variable(i,j,nzi),i=2,nxp)
            end do  
        end if
    end subroutine varsave3D_DAT
    !=======================================================================================
    !----------------------------------------------------------
    !       Write the selected variables for 2D problems
    !---------------------- Description -----------------------
    ! <INPUT>
    ! face     - Selected face ID
    ! name     - Name of the variable to be displayed
    ! variable - Variable to be saved
    ! fl       - File ID
    ! <EXTERNAL ROUTINES>
    !----------------------------------------------------------
    subroutine varsave2D_DAT(variable,name,face,fl)
        double precision, intent (in) :: variable(:,:,:)
        character(len=*), intent (in):: name
        integer, intent (in):: face,fl
        integer :: i,j
        character(len=60) :: display_name    
        if(face == 1 .or. face == 2)then !WEST-EAST
            display_name = trim('   x-direction Wall :')//trim(name)
            write(fl,'(a)') trim(display_name)
            write(fl,'(a)')'==================================' 
            write(fl,'(a)')'     Y          West        East  '
            write(fl,'(a)')'=================================='  
            do j=2,nyp
                write(fl,'(ES10.3,ES12.3,ES12.3)') yc(j),variable(1,j,2),variable(nxi,j,2)
            end do
            write(fl,'(a)')'=================================='
        end if
        if(face == 3 .or. face == 4)then !SOUTH-NORTH
            display_name = trim('   y-direction Wall :')//trim(name)
            write(fl,'(a)') trim(display_name)
            write(fl,'(a)')'==================================' 
            write(fl,'(a)')'     X         South       North  '
            write(fl,'(a)')'=================================='  
            do i=2,nxp
                write(fl,'(ES10.3,ES12.3,ES12.3)') xc(i),variable(i,1,2),variable(i,nyi,2)
            end do
            write(fl,'(a)') '=================================='
        end if
    end subroutine varsave2D_DAT
    !=======================================================================================
    !----------------------------------------------------------
    !       Write the selected variables for 1D problems
    !---------------------- Description -----------------------
    ! <INPUT>
    ! face     - Selected face ID
    ! name     - Name of the variable to be displayed
    ! variable - Variable to be saved
    ! fl       - File ID
    ! <EXTERNAL ROUTINES>
    !----------------------------------------------------------
    subroutine varsave1D_DAT(variable,name,face,fl)
        double precision, intent (in) :: variable(:,:,:)
        character(len=*), intent (in):: name
        integer, intent (in):: face,fl
        if(face == 1 .or. face == 2)then !WEST-EAST
            write(fl,'(a)') trim(name)
            write(fl,'(a)')'========================' 
            write(fl,'(a)')'     West        East   '
            write(fl,'(a)')'========================'  
            write(fl,'(ES12.3,ES12.3)') variable(1,2,2),variable(nxi,2,2)
            write(fl,'(a)')'========================'
        end if
    end subroutine varsave1D_DAT
    !=======================================================================================
    !----------------------------------------------------------
    ! Modifies the name to be displayed according to the wall
    !---------------------- Description -----------------------
    ! <INPUT>
    ! name             - Name of the variable to be displayed
    ! var_name_saveDAT - Modified name
    !----------------------------------------------------------
    function var_name_saveDAT(name,face)
        character(len=60):: var_name_saveDAT
        character(len=*):: name
        integer :: face
        if(face == 1)then !WEST
            var_name_saveDAT = trim('   West Wall :')//trim(name)
        else if(face == 2)then !EAST
            var_name_saveDAT = trim('   East Wall :')//trim(name)
        else if(face == 3)then !SOUTH
            var_name_saveDAT = trim('   South Wall :')//trim(name)
        else if(face == 4)then !NORTH
            var_name_saveDAT = trim('   North Wall :')//trim(name)
        else if(face == 5)then !BOTTOM
            var_name_saveDAT = trim('   Bottom Wall :')//trim(name)
        else if(face == 6)then !TOP
            var_name_saveDAT = trim('   Top Wall :')//trim(name)
        end if
    end function var_name_saveDAT
    !=======================================================================================
    !----------------------------------------------------------
    !                 Main slices save routine
    !---------------------- Description -----------------------
    ! <INPUT>
    ! slc_axis  - Slices axis array
    ! slc_perc  - Slices percent array
    ! slc_point - Slices point array
    ! slc_node  - Slices node array
    ! slc_ID    - Slices variable ID array
    ! FID       - File ID
    ! slc       - Slice counter
    ! nslices   - Number of slices 
    ! <EXTERNAL ROUTINES>
    ! slc_var_select - Select the variable slice to be saved
    !----------------------------------------------------------
    subroutine main_slice()
        integer :: slc,FID=5
        open(unit=FID, file = 'output/RTS_Slices.dat', action="write", status="replace")
        do slc=1,nslices
            call slc_var_select(slc_axis(slc),slc_perc(slc),slc_point(slc),slc_node(slc),slc_ID(slc),FID)
        end do
        close(FID)
        write(*,'(a)')'RTS_Slices.dat Have Been Successfully Saved in ./output'
    end subroutine main_slice
    !=======================================================================================
    !----------------------------------------------------------
    !           Select the variable slice to be saved
    !---------------------- Description -----------------------
    ! <INPUT>
    ! CAID    - Character variable
    ! VID     - Variable ID
    ! FID     - File ID
    ! AID     - Axis ID
    ! PID     - Plane ID
    ! <EXTERNAL ROUTINES>
    ! AxID        - Find the axis ID
    ! get_PID     - Find the plane ID
    ! get_name    - Get the title to be displayed
    ! slice_print - Write the selected slice in the file
    !----------------------------------------------------------
    subroutine slc_var_select(CAID,Percent,Point,Node,VID,FID)
        integer :: AID,PID,Node,Percent,VID,FID
        double precision :: Point
        character(len=*):: CAID
        character(len=70)::name
        AID = AxID(CAID)
        PID = get_PID(AID,Percent,Point,Node)
        name = trim(get_name(PID,AID,CAID))
        if(VID == 1)then
            call slice_print(T_energy,'Temperature :',name,PID,FID,AID)
        else if(VID == 2)then
            call slice_print(G,'Incident Radiation :',name,PID,FID,AID)
        else if(VID == 3)then
            call slice_print(S_rad,'Divergent Radiative Flux:',name,PID,FID,AID)
        else if(VID == 4)then
            call slice_print(IBlack,'Blackbody Radiation :',name,PID,FID,AID)
        else if(VID == 5)then
            call slice_print(cappa,'Absorption Coefficient :',name,PID,FID,AID)
        else if(VID == 6)then
            call slice_print(sigma,'Scattering Coefficient :',name,PID,FID,AID)
        else if(VID == 7)then
            call slice_print(beta,'Extinction Coefficient :',name,PID,FID,AID)
        else if(VID == 8)then
            call slice_print(K_term,'Thermal Conductivity :',name,PID,FID,AID)
        else if(VID == 9)then
            call slice_print(epsilon_rad,'Radiation Emissivities :',name,PID,FID,AID)
        end if
    end subroutine slc_var_select
    !=======================================================================================
    !----------------------------------------------------------
    !              Get the title to be displayed
    !---------------------- Description -----------------------
    ! <INPUT>
    ! CAID     - Character variable
    ! AID      - Axis ID
    ! PID      - Plane ID
    ! <OUTPUT>
    ! get_name - Title to be displayed
    !----------------------------------------------------------
    function get_name(PID,AID,CAID)
        integer :: AID,PID
        character(len=*) :: CAID
        character(len=90):: get_name
        if(AID == 1)then
            write(get_name,'(a,ES10.3,a)') 'x-direction slice (x =',xc(PID),' m)-zy plane'
        else if(AID == 2)then
            write(get_name,'(a,ES10.3,a)') 'y-direction slice (y =',yc(PID),' m)-xz plane'
        else if(AID == 3)then
            write(get_name,'(a,ES10.3,a)') 'z-direction slice (z =',zc(PID),' m)-xy plane'
        end if
    end function get_name
    !=======================================================================================
    !----------------------------------------------------------
    !           Write the selected slice in the file
    !---------------------- Description -----------------------
    ! <INPUT>
    ! variable - Variable to be saved
    ! varname  - Variable name
    ! Title    - Variable ID
    ! FID      - File ID
    ! AID      - Axis ID
    ! PID      - Plane ID
    !----------------------------------------------------------
    subroutine slice_print(variable,varname,Title,PID,FID,AID)
        character(len=*):: Title,varname
        double precision, intent (in) :: variable(:,:,:)
        integer :: i,j,k,PID,FID,AID
        Title = trim(varname)//trim(Title)
        call header_var_writeDAT(Title,FID)
        if(AID == 1)then !x-direction slice
            do j=nyi,1,-1
                write(FID,*) (variable(PID,j,k),k=1,nzi)
            end do
        else if(AID == 2)then !y-direction slice
            do k=nzi,1,-1
                write(FID,*) (variable(i,PID,k),i=1,nxi)
            end do
        else if(AID == 3)then !z-direction slice
            do j=nyi,1,-1
                write(FID,*) (variable(i,j,PID),i=1,nxi)
            end do
        end if
    end subroutine slice_print
    !=======================================================================================
    !----------------------------------------------------------
    !                     Find the axis ID
    !---------------------- Description -----------------------
    ! <INPUT>
    ! CAID    - Character variable
    ! <OUTPUT>
    ! AxID    - Axis ID
    !----------------------------------------------------------
    function AxID(CAID)
        character(len=*):: CAID
        integer :: AxID
        if(CAID == 'x'.or. CAID == 'X')then      !x-direction slice
            AxID = 1
        else if(CAID == 'y'.or. CAID == 'Y')then !y-direction slice
            AxID = 2
        else if(CAID == 'z'.or. CAID == 'Z')then !z-direction slice
            AxID = 3
        end if
    end function AxID
    !=======================================================================================
    !----------------------------------------------------------
    !                    Find the plane ID
    !---------------------- Description -----------------------
    ! <INPUT>
    ! AID     - Axis ID
    ! OPslice - Slicing method
    ! <OUTPUT>
    ! get_PID - Plane ID
    ! <EXTERNAL ROUTINES>
    ! get_PID_percent - Find the plane ID using the percent mode
    ! get_PID_point   - Find the plane ID using the point mode
    !----------------------------------------------------------
    function get_PID(AID,Percent,Point,Node)
        integer :: AID,get_PID,Node,Percent
        double precision :: Point
        if(OPslice == 1)then      !Percent
            get_PID = get_PID_percent(Percent,AID)
        else if(OPslice == 2)then !Point
            get_PID = get_PID_point(Point,AID)
        else if(OPslice == 3)then !Control Volume Node
            get_PID = Node
        end if
    end function get_PID
    !=======================================================================================
    !----------------------------------------------------------
    !         Find the plane ID using the percent mode
    !---------------------- Description -----------------------
    ! <INPUT>
    ! AID             - Axis ID
    ! <OUTPUT>
    ! get_PID_percent - Plane ID
    ! <EXTERNAL ROUTINES>
    ! index_finder - Finds the index corresponding to a given position
    !----------------------------------------------------------
    function get_PID_percent(Percent,AID)

        integer :: AID,get_PID_percent,Percent
        if(AID == 1)then
            get_PID_percent = index_finder((Percent/100.0d0)*lx,xc)
        else if(AID == 2)then
            get_PID_percent = index_finder((Percent/100.0d0)*ly,yc)
        else if(AID == 3)then
            get_PID_percent = index_finder((Percent/100.0d0)*lz,zc)
        end if
    end function get_PID_percent
    !=======================================================================================
    function get_PID_point(Point,AID)
    !----------------------------------------------------------
    !          Find the plane ID using the point mode
    !---------------------- Description -----------------------
    ! <INPUT>
    ! AID           - Axis ID
    ! <OUTPUT>
    ! get_PID_point - Plane ID
    ! <EXTERNAL ROUTINES>
    ! index_finder - Finds the index corresponding to a given position
    !----------------------------------------------------------
        integer :: AID,get_PID_point
        double precision :: Point
        if(AID == 1)then
            get_PID_point = index_finder(Point,xc)
        else if(AID == 2)then
            get_PID_point = index_finder(Point,yc)
        else if(AID == 3)then
            get_PID_point = index_finder(Point,zc)
        end if
    end function get_PID_point
    !=======================================================================================
    !----------------------------------------------------------
    !    Finds the index corresponding to a given position
    !---------------------- Description -----------------------
    ! <INPUT>
    ! sc  - s value at the control-volume center
    ! lp  - Given spatial position
    ! <OUTPUT>
    ! index_finder - Index corresponding to a given position
    !----------------------------------------------------------
    function index_finder(lp,sc)
        implicit none
        integer :: l,index_finder
        double precision :: lp,sc(:),aux
        aux = (sc(2) - sc(1))*0.5d0
        do l = 1,size(sc)
            if(sc(l)>= lp - aux)then
                index_finder = l
                exit
            end if
            if(sc(l)>= lp + aux)then
                index_finder = l
                exit
            end if
        end do
    end function index_finder
    !=======================================================================================
    !----------------------------------------------------------
    !      Main routine for calculating the Nusselt number
    !---------------------- Description -----------------------
    ! <INPUT>
    ! T_energy - Energy temperature
    ! Q_radw   - Walls radiative fluxes
    ! K_termc  - Constant thermal conductivity
    ! <INTERNAL>
    ! NuC      - Convective Nusselt distribution
    ! NuR      - Radiative Nusselt distribution
    ! NuCm     - Averaged convective Nusselt number
    ! NuRm     - Averaged radiative Nusselt number
    ! <EXTERNAL ROUTINES>
    ! q_convec - Calculate the convective flux in the walls
    ! nusselt_CR    - Compute the Nusselt number (Convective/Radiative)
    ! nusselt_AVG   - Compute the average Nusselt number
    ! file_open2D   - Open the file in which the plans will be saved
    ! varsave2D_VTK - Write the selected variables in the selected plans
    ! header_var_writeDAT - Print an identification bar in the dat file
    !----------------------------------------------------------
    subroutine Nusselt_main()
        double precision :: NuC(nxi,nyi,nzi),NuR(nxi,nyi,nzi),NuCm,NuRm
        integer :: j,k,FID=16
        NuC = 0.0d0; NuR = 0.0d0
        !-------------- Nusselt Distribution --------------!
        do k=2,nzp
            do j=2,nyp
                DT = dabs(T_energy(1,j,k) - T_energy(nxi,j,k))
                !West Wall
                NuC(1,j,k) = nusselt_CR(K_termc,lx,DT,q_convec(j,k,K_termc,xc,T_energy,1))
                NuR(1,j,k) = nusselt_CR(K_termc,lx,DT,Q_radw(1,j,k))
                !East Wall
                NuC(nxi,j,k) = nusselt_CR(K_termc,lx,DT,q_convec(j,k,K_termc,xc,T_energy,2))
                NuR(nxi,j,k) = nusselt_CR(K_termc,lx,DT,Q_radw(nxi,j,k))
            end do
        end do
        !---------------- Nusselt vtk file ----------------!
        call file_open2D(z,y,'Nusselt',FID)
        call varsave2D_VTK(NuC,'NuC',1,FID)
        call varsave2D_VTK(NuC,'NuC',2,FID)
        call varsave2D_VTK(NuR,'NuR',1,FID)
        call varsave2D_VTK(NuR,'NuR',2,FID)
        close(FID)
        !-------------- Local Nusselt Number --------------!
        open(unit=FID, file = 'output/RTS_Nusselt.dat', action="write", status="replace")
        call header_var_writeDAT('Convective Nusselt West Wall',FID)
        do j=nyi,1,-1
            write(FID,*) (NuC(1,j,k),k=1,nzi)
        end do
        call header_var_writeDAT('Convective Nusselt East Wall',FID)
        do j=nyi,1,-1
            write(FID,*) (NuC(nxi,j,k),k=1,nzi)
        end do
        call header_var_writeDAT('Radiative Nusselt West Wall',FID)
        do j=nyi,1,-1
            write(FID,*) (NuR(1,j,k),k=1,nzi)
        end do
        call header_var_writeDAT('Radiative Nusselt East Wall',FID)
        do j=nyi,1,-1
            write(FID,*) (NuR(nxi,j,k),k=1,nzi)
        end do
        call header_var_writeDAT('Total Nusselt West Wall',FID)
        do j=nyi,1,-1
            write(FID,*) (NuC(1,j,k)+NuR(1,j,k),k=1,nzi)
        end do
        call header_var_writeDAT('Total Nusselt East Wall',FID)
        do j=nyi,1,-1
            write(FID,*) (NuC(nxi,j,k)+NuR(nxi,j,k),k=1,nzi)
        end do
        close(FID)
        !------------- Average Nusselt Number -------------!
        !Average Nusselt Number
        open(unit=FID, file = 'output/RTS_Nusselt_AVG.dat', action="write", status="replace")
        write(FID,*) 'NuConv,NuRad,Nu'
        write(FID,*) 'West Wall'
        NuCm = nusselt_AVG(1,NuC)
        NuRm = nusselt_AVG(1,NuR)
        write(FID,*) NuCm,NuRm,NuCm+NuRm
        write(FID,*) 'East Wall'
        NuCm = nusselt_AVG(nxi,NuC)
        NuRm = nusselt_AVG(nxi,NuR)
        write(FID,*) NuCm,NuRm,NuCm+NuRm
        close(FID)
        write(*,'(a)')'Nusselt Results Have Been Successfully Saved in ./output'
    end subroutine Nusselt_main
    !=======================================================================================
    !----------------------------------------------------------
    !        Calculate the convective flux in the walls
    !---------------------- Description -----------------------
    ! <INPUT>
    ! Temp     - Energy temperature
    ! Sc       - S value at the control-volume center
    ! Kt       - Constant thermal conductivity
    ! WID      - Wall ID
    ! i,j      - Field indexes
    ! <OUTPUT>
    ! q_convec - Convective flux
    !----------------------------------------------------------
    function q_convec(j,k,Kt,Sc,Temp,WID)
        integer :: j,k,WID
        double precision :: q_convec,Kt,Sc(:),Temp(:,:,:)
        if(WID == 1)then      !WEST
            q_convec = - Kt*(Temp(2,j,k) - Temp(1,j,k))/(Sc(2) - Sc(1))
        else if(WID == 2)then !EAST
            q_convec = - Kt*(Temp(nxi,j,k) - Temp(nxp,j,k))/(Sc(nxi) - Sc(nxp))
        end if
    end function q_convec
    !=======================================================================================
    !----------------------------------------------------------
    !     Compute the Nusselt number (Convective/Radiative)
    !---------------------- Description -----------------------
    ! <INPUT>
    ! DT     - Energy temperature difference
    ! q_flux - Walls (Convective/Radiative) fluxes
    ! Kt     - Constant thermal conductivity
    ! Ls     - Domain length
    ! <OUTPUT>
    ! nusselt_CR - (Convective/Radiative) Nusselt number
    !----------------------------------------------------------
    function nusselt_CR(Kt,Ls,DT,q_flux)
        double precision :: Kt,Ls,DT,q_flux,nusselt_CR
        nusselt_CR = Ls*dabs(q_flux)/(Kt*DT)
    end function nusselt_CR
    !=======================================================================================
    !----------------------------------------------------------
    !            Compute the average Nusselt number
    !---------------------- Description -----------------------
    ! <INPUT>
    ! NUS - Nusselt distribution
    ! i   - Wall index
    ! K_termc  - Constant thermal conductivity
    ! <OUTPUT>
    ! nusselt_AVG - Averaged Nusselt number
    !----------------------------------------------------------
    function nusselt_AVG(i,NUS)
        integer :: i,j,k
        double precision :: Nu,nusselt_AVG,NUS(:,:,:)
        Nu = 0.0d0
        do k=2,nzp
            do j=2,nyp
                Nu = Nu + NUS(i,j,k)*dyp(j)
            end do
        end do
        nusselt_AVG = Nu/ly
    end function nusselt_AVG
    !=======================================================================================
    !----------------------------------------------------------
    !                 Main lineup save routine
    !---------------------- Description -----------------------
    ! <INPUT>
    ! lnp_axis  - Lineup axis array
    ! lnp_point - Lineup point array
    ! lnp_node  - Lineup node array
    ! lnp_ID    - Lineup variable ID array
    ! FID       - File ID
    ! lnp       - Lineup method
    ! nlnp      - Number of lineups
    ! <EXTERNAL ROUTINES>
    ! lineup_select - Select the variable lineup to be saved
    !----------------------------------------------------------
    subroutine main_lineup()
        integer :: lnp,FID=23
        open(unit=FID, file = 'output/RTS_lineup.dat', action="write", status="replace")
        do lnp=1,nlnp
            call lineup_select(lnp_axis(lnp),lnp_point(lnp,1),lnp_point(lnp,2),&
            lnp_node(lnp,1),lnp_node(lnp,2),lnp_ID(lnp),FID)
        end do
        close(FID)
        write(*,'(a)')'RTS_lineup.dat Have Been Successfully Saved in ./output'
    end subroutine main_lineup
    !=======================================================================================
    !----------------------------------------------------------
    !           Select the variable lineup to be saved
    !---------------------- Description -----------------------
    ! <INPUT>
    ! CAID    - Character variable
    ! VID     - Variable ID
    ! FID     - File ID
    ! AID     - Axis ID
    ! PID     - Plane ID
    ! <EXTERNAL ROUTINES>
    ! AxID         - Find the axis ID
    ! ijkp_finder  - Find the lineup coordinates
    ! lineup_write - Write the selected lineup in the file
    !----------------------------------------------------------
    subroutine lineup_select(CAID,a,b,ap,bp,VID,FID)
        integer :: ip,jp,kp,ap,bp,AID,VID,FID
        double precision :: a,b
        character(len=*):: CAID
        ! The variable AID is set equal to the result of the function AxID(CAID)
        AID = AxID(CAID)
        ! Call function ijkp_finder to find ip, jp, kp, ap, and bp
        call ijkp_finder(AID,ip,jp,kp,a,b,ap,bp)
        ! Select the variable to write to the output file based on the input variable ID
        select case (VID)
            case (1) ! If VID is 1, write T_energy
                call lineup_write(T_energy, AID, ip, jp, kp, CAID, 'Temperature', FID)
            case (2) ! If VID is 2, write G
                call lineup_write(G, AID, ip, jp, kp, CAID, 'G', FID)
            case (3) ! If VID is 3, write S_rad
                call lineup_write(S_rad, AID, ip, jp, kp, CAID, 'S_rad', FID)
            case (4) ! If VID is 4, write IBlack
                call lineup_write(IBlack, AID, ip, jp, kp, CAID, 'Ib', FID)
            case (5) ! If VID is 5, write cappa
                call lineup_write(cappa, AID, ip, jp, kp, CAID, 'cappa', FID)
            case (6) ! If VID is 6, write sigma
                call lineup_write(sigma, AID, ip, jp, kp, CAID, 'sigma', FID)
            case (7) ! If VID is 7, write beta
                call lineup_write(beta, AID, ip, jp, kp, CAID, 'beta', FID)
            case (8) ! If VID is 8, write K_term
                call lineup_write(K_term, AID, ip, jp, kp, CAID, 'Kcond', FID)
            case (9) ! If VID is 9, write epsilon_rad
                call lineup_write(epsilon_rad, AID, ip, jp, kp, CAID, 'Emissivity', FID)
            case (10) ! If VID is 10, write Q_radw
                call lineup_write(Q_radw, AID, ip, jp, kp, CAID, 'Q_radw', FID)
            case default ! If VID is not between 1 and 10, write an error message to standard output
                write(*, *) "Error: Invalid variable ID."
        end select
    end subroutine lineup_select
    !=======================================================================================
    !----------------------------------------------------------
    !                Find the lineup coordinates
    !---------------------- Description -----------------------
    ! <INPUT>
    ! AID      - Axis ID
    ! OPlnp    - Lineup method
    ! a,b      - Lineup points
    ! ap,bp    - Lineup nodes
    ! <OUTPUT>
    ! ip,jp,kp - lineup coordinates
    ! <EXTERNAL ROUTINES>
    ! index_finder   - Finds the index corresponding to a given position
    !----------------------------------------------------------
    subroutine ijkp_finder(AID,ip,jp,kp,a,b,ap,bp)
        double precision, intent(in) :: a,b
        integer, intent(out) :: ip,jp,kp
        integer :: AID,ap,bp
        if(OPlnp == 1)then
            if(AID == 1)then      !X:(x,a,b)
                jp = index_finder(a,yc)
                kp = index_finder(b,zc)
            else if(AID == 2)then !Y:(a,y,b) 
                ip = index_finder(a,xc)
                kp = index_finder(b,zc)
            else if(AID == 3)then !Z:(a,b,z) 
                ip = index_finder(a,xc)
                jp = index_finder(b,yc)
            end if
        else if(OPlnp == 2)then
            if(AID == 1)then      !X:(xp,ap,bp)
                jp = ap
                kp = bp
            else if(AID == 2)then !Y:(ap,yp,bp) 
                ip = ap
                kp = bp
            else if(AID == 3)then !Z:(ap,bp,z) 
                ip = ap
                jp = bp
            end if
        end if
        if(DIMEN == 2) kp = 2
    end subroutine ijkp_finder
    !=======================================================================================
    !----------------------------------------------------------
    !           Write the selected lineup in the file
    !---------------------- Description -----------------------
    ! <INPUT>
    ! var      - Variable to be saved
    ! vname    - Variable name
    ! Title    - Variable ID
    ! CAID     - Character variable
    ! FID      - File ID
    ! AID      - Axis ID
    ! ip,jp,kp - lineup coordinates
    ! <EXTERNAL ROUTINES>
    ! get_name_lineup - Get the title to be displayed in lineup
    !----------------------------------------------------------
    subroutine lineup_write(var,AID,ip,jp,kp,CAID,vname,FID)
        character(len=*):: CAID,vname
        integer :: AID,FID,ip,jp,kp,l
        double precision, intent(in):: var(:,:,:)
        write(FID,'(a,a)') trim(vname),get_name_lineup(AID,ip,jp,kp)
        write(FID,'(a)') '=====================================' 
        write(FID,'(a,a,a,a)')'     ',trim(CAID),'            ',trim(vname)
        write(FID,'(a)') '=====================================' 
        if(AID == 1)then
            do l=1,nxi
                write(FID,'(ES10.3,F15.5)') xc(l), var(l,jp,kp)
            end do
        else if(AID == 2)then
            do l=1,nyi
                write(FID,'(ES10.3,F15.5)') yc(l), var(ip,l,kp)
            end do
        else if(AID == 3)then
            do l=1,nzi
                write(FID,'(ES10.3,F15.5)') zc(l), var(ip,jp,l)
            end do
        end if
    end subroutine lineup_write
    !=======================================================================================
    !----------------------------------------------------------
    !          Get the title to be displayed in lineup
    !---------------------- Description -----------------------
    ! <INPUT>
    ! AID      - Axis ID
    ! ip,jp,kp - lineup coordinates
    ! <OUTPUT>
    ! get_name_lineup - Title to be displayed
    !----------------------------------------------------------
    function get_name_lineup(AID,ip,jp,kp)
        integer :: AID,ip,jp,kp
        character(len=70):: get_name_lineup
        if(AID == 1)then
            write(get_name_lineup,'(a,F5.3,a,F5.3,a)') ':(x,',yc(jp),',',zc(kp),')'
        else if(AID == 2)then
            write(get_name_lineup,'(a,F5.3,a,F5.3,a)') ':(',xc(ip),',y,',zc(kp),')'
        else if(AID == 3)then
            write(get_name_lineup,'(a,F5.3,a,F5.3,a)') ':(',xc(ip),',',yc(jp),',z)'
        end if
    end function get_name_lineup
    !=======================================================================================
end module output
