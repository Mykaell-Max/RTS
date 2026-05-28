module functions

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
    !                            Variable initialization functions
    !========================================================================================
    !----------------------------------------------------------
    !             Temperature space-function
    !---------------------- Description -----------------------
    ! <INPUT>
    ! xc  - x value at the control-volume center
    ! yc  - y value at the control-volume center
    ! zc  - z value at the control-volume center
    ! i   - x-direction index
    ! j   - y-direction index
    ! k   - z-direction index
    ! <OUTPUT>
    ! Temp_field - Temperature data
    !NOTE: use xc(i), yc(j) and zc(k) for spatial positions
    !----------------------------------------------------------
    function Temp_field(i,j,k)
        implicit none
        integer, intent(in) :: i,j,k
        double precision :: Temp_field,Tc,r,func

        if(yc(j)>= 0.0d0 .and. yc(j)<= 0.375d0)then
            Tc = 400.0d0 + 1400.0d0*(yc(j)/0.375d0)
        else if(yc(j)> 0.375d0 .and. yc(j)<= 4.0d0)then
            Tc = 1800.0d0 - 250.0d0*yc(j)
        end if   
        
        r = sqrt((xc(i) - 1.0d0)**2.0d0 + (zc(k) - 1.0d0)**2.0d0)
        
        if(r <= 1.000001d0)then
            func = 1.0d0 - 3.0d0*(r)**2.0d0 + 2.0d0*(r)**3.0d0
            Temp_field = (Tc - 800.0d0)*func + 800.0d0
        else
            Temp_field = 800.0d0
        end if

    end function Temp_field
    !========================================================================================
    !----------------------------------------------------------
    !           Incident radiation space-function
    !---------------------- Description -----------------------
    ! <INPUT>
    ! xc  - x value at the control-volume center
    ! yc  - y value at the control-volume center
    ! zc  - z value at the control-volume center
    ! i   - x-direction index
    ! j   - y-direction index
    ! k   - z-direction index
    ! <OUTPUT>
    ! G_field - Incident radiation data
    !NOTE: use xc(i), yc(j) and zc(k) for spatial positions
    !----------------------------------------------------------
    function G_field(i,j,k)
        implicit none
        integer, intent(in) :: i,j,k
        double precision :: G_field
        
        G_field = 0.0d0

    end function G_field
    function cappa_field(i,j,k)
        !----------------------------------------------------------
        !          Absorption coefficient space-function
        ! Made by:     G.S.Rodrigues (November 2020)
        ! Modified by: G.S.Rodrigues (July 2021)
        !---------------------- Description -----------------------
        ! <INPUT>
        ! xc  - x value at the control-volume center
        ! yc  - y value at the control-volume center
        ! zc  - z value at the control-volume center
        ! i   - x-direction index
        ! j   - y-direction index
        ! k   - z-direction index
        ! <OUTPUT>
        ! cappa_field - Absorption coefficient data
        !NOTE: use xc(i), yc(j) and zc(k) for spatial positions
        !----------------------------------------------------------
        implicit none
        integer, intent(in) :: i,j,k
        double precision :: cappa_field,betaux
        
        betaux = 0.9d0*(1.0d0 - 2.0d0*abs(xc(i) - 0.5d0))*(1.0d0 - 2.0d0*abs(yc(j) - 0.5d0))
        betaux = betaux*(1.0d0 - 2.0d0*abs(zc(k) - 0.5d0)) + 0.1d0
        cappa_field = betaux
        !cappa_field = 0.1d0*betaux

    end function cappa_field
    !========================================================================================
    function sigma_field(i,j,k)
        !----------------------------------------------------------
        !           Scattering coefficient space-function
        ! Made by:     G.S.Rodrigues (November 2020)
        ! Modified by: G.S.Rodrigues (July 2021)
        !---------------------- Description -----------------------
        ! <INPUT>
        ! xc  - x value at the control-volume center
        ! yc  - y value at the control-volume center
        ! zc  - z value at the control-volume center
        ! i   - x-direction index
        ! j   - y-direction index
        ! k   - z-direction index
        ! <OUTPUT>
        ! sigma_field - Scattering coefficient data
        !NOTE: use xc(i), yc(j) and zc(k) for spatial positions
        !----------------------------------------------------------
        implicit none
        integer, intent(in) :: i,j,k
        double precision :: sigma_field,betaux
                
        betaux = 0.9d0*(1.0d0 - 2.0d0*abs(xc(i) - 0.5d0))*(1.0d0 - 2.0d0*abs(yc(j) - 0.5d0))
        betaux = betaux*(1.0d0 - 2.0d0*abs(zc(k) - 0.5d0)) + 0.1d0
        sigma_field = 0.9d0*betaux
        
    end function sigma_field
    !========================================================================================
    !========================================================================================
    !----------------------------------------------------------
    !           Thermal conductivity space-function
    !---------------------- Description -----------------------
    ! <INPUT>
    ! xc  - x value at the control-volume center
    ! yc  - y value at the control-volume center
    ! zc  - z value at the control-volume center
    ! i   - x-direction index
    ! j   - y-direction index
    ! k   - z-direction index
    ! <OUTPUT>
    ! Kterm_field - Thermal conductivity data
    !NOTE: use xc(i), yc(j) and zc(k) for spatial positions
    !----------------------------------------------------------
    function Kterm_field(i,j,k)
        implicit none
        integer, intent(in) :: i,j,k
        double precision :: Kterm_field
        
        Kterm_field = 0.0d0

    end function Kterm_field
    !========================================================================================
    !----------------------------------------------------------
    !             CO2 mole fraction space-function
    !---------------------- Description -----------------------
    ! <INPUT>
    ! xc   - x value at the control-volume center
    ! yc   - y value at the control-volume center
    ! zc   - z value at the control-volume center
    ! i    - x-direction index
    ! j    - y-direction index
    ! k    - z-direction index
    ! <OUTPUT>
    ! Y_CO2field - CO2 mole fraction concentration
    !NOTE: use xc(i), yc(j) and zc(k) for spatial positions
    !----------------------------------------------------------
    function Y_CO2field(i,j,k)
        implicit none
        integer, intent(in) :: i,j,k
        double precision :: Y_CO2field
        Y_CO2field = XCO2_g
        
    end function Y_CO2field
    !========================================================================================
    !----------------------------------------------------------
    !             H2O mole fraction space-function
    !---------------------- Description -----------------------
    ! <INPUT>
    ! xc   - x value at the control-volume center
    ! yc   - y value at the control-volume center
    ! zc   - z value at the control-volume center
    ! i    - x-direction index
    ! j    - y-direction index
    ! k    - z-direction index
    ! <OUTPUT>
    ! Y_H2Ofield - H2O mole fraction concentration
    !NOTE: use xc(i), yc(j) and zc(k) for spatial positions
    !----------------------------------------------------------
    function Y_H2Ofield(i,j,k)
        implicit none
        integer, intent(in) :: i,j,k
        double precision :: Y_H2Ofield
        Y_H2Ofield = XH2O_g
        
    end function Y_H2Ofield
    !========================================================================================
    !----------------------------------------------------------
    !           Species total pressure space-function
    !---------------------- Description -----------------------
    ! <INPUT>
    ! xc   - x value at the control-volume center
    ! yc   - y value at the control-volume center
    ! zc   - z value at the control-volume center
    ! i    - x-direction index
    ! j    - y-direction index
    ! k    - z-direction index
    ! <OUTPUT>
    ! Pk_field - Species partial pressure PCO2 + PH2O (atm)
    !NOTE: use xc(i), yc(j) and zc(k) for spatial positions
    !----------------------------------------------------------
    function Pk_field(i,j,k)
        implicit none
        integer, intent(in) :: i,j,k
        double precision :: Pk_field
        Pk_field = P_g
        
    end function Pk_field
    !========================================================================================
    !----------------------------------------------------------
    !     Set the non-uniform temperature at the west wall
    !---------------------- Description -----------------------
    ! <INPUT>
    ! yc  - y value at the control-volume center
    ! zc  - z value at the control-volume center
    ! j   - y-direction index
    ! k   - z-direction index
    ! <OUTPUT>
    ! wallTemp_west - West wall temperature data
    !NOTE: use yc(j) and zc(k) for spatial positions
    !----------------------------------------------------------
    function wallTemp_west(j,k)
        implicit none
        integer :: j,k
        double precision :: wallTemp_west
        
        wallTemp_west = wall_temp(1)

    end function wallTemp_west
    !========================================================================================
    !----------------------------------------------------------
    !     Set the non-uniform temperature at the east wall
    !---------------------- Description -----------------------
    ! <INPUT>
    ! yc  - y value at the control-volume center
    ! zc  - z value at the control-volume center
    ! j   - y-direction index
    ! k   - z-direction index
    ! <OUTPUT>
    ! wallTemp_east - West wall temperature data
    !NOTE: use yc(j) and zc(k) for spatial positions
    !----------------------------------------------------------
    function wallTemp_east(j,k)
        implicit none
        integer :: j,k
        double precision :: wallTemp_east
        
        wallTemp_east = wall_temp(2)

    end function wallTemp_east
    !========================================================================================
    !----------------------------------------------------------
    !    Set the non-uniform temperature at the south wall
    !---------------------- Description -----------------------
    ! <INPUT>
    ! xc  - x value at the control-volume center
    ! zc  - z value at the control-volume center
    ! i   - x-direction index
    ! k   - z-direction index
    ! <OUTPUT>
    ! wallTemp_south - West wall temperature data
    !NOTE: use xc(i) and zc(k) for spatial positions
    !----------------------------------------------------------
    function wallTemp_south(i,k)
        implicit none
        integer :: i,k
        double precision :: wallTemp_south
        
        wallTemp_south = wall_temp(3)

    end function wallTemp_south
    !========================================================================================
    !----------------------------------------------------------
    !    Set the non-uniform temperature at the north wall
    !---------------------- Description -----------------------
    ! <INPUT>
    ! xc  - x value at the control-volume center
    ! zc  - z value at the control-volume center
    ! i   - x-direction index
    ! k   - z-direction index
    ! <OUTPUT>
    ! wallTemp_north - West wall temperature data
    !NOTE: use xc(i) and zc(k) for spatial positions
    !----------------------------------------------------------
    function wallTemp_north(i,k)

        implicit none
        integer :: i,k
        double precision :: wallTemp_north
        
        wallTemp_north = wall_temp(4)

    end function wallTemp_north
    !========================================================================================
    !----------------------------------------------------------
    !    Set the non-uniform temperature at the bottom wall
    !---------------------- Description -----------------------
    ! <INPUT>
    ! xc  - x value at the control-volume center
    ! yc  - y value at the control-volume center
    ! i   - x-direction index
    ! j   - y-direction index
    ! <OUTPUT>
    ! wallTemp_bottom - West wall temperature data
    !NOTE: use xc(i) and yc(j) for spatial positions
    !----------------------------------------------------------
    function wallTemp_bottom(i,j)
        implicit none
        integer :: i,j
        double precision :: wallTemp_bottom
        
        wallTemp_bottom = wall_temp(5)

    end function wallTemp_bottom
    !========================================================================================
    !----------------------------------------------------------
    !     Set the non-uniform temperature at the top wall
    !---------------------- Description -----------------------
    ! <INPUT>
    ! xc  - x value at the control-volume center
    ! yc  - y value at the control-volume center
    ! i   - x-direction index
    ! j   - y-direction index
    ! <OUTPUT>
    ! wallTemp_top - West wall temperature data
    !NOTE: use xc(i) and yc(j) for spatial positions
    !----------------------------------------------------------
    function wallTemp_top(i,j)
        implicit none
        integer :: i,j
        double precision :: wallTemp_top
        
        wallTemp_top = wall_temp(6)

    end function wallTemp_top
    !========================================================================================
    !----------------------------------------------------------
    !     Set the non-uniform emissivity at the west wall
    !---------------------- Description -----------------------
    ! <INPUT>
    ! yc  - y value at the control-volume center
    ! zc  - z value at the control-volume center
    ! j   - y-direction index
    ! k   - z-direction index
    ! <OUTPUT>
    ! wallemissivity_west - West wall emissivity data
    !NOTE: use yc(j) and zc(k) for spatial positions
    !----------------------------------------------------------
    function wallemissivity_west(j,k)
        implicit none
        integer :: j,k
        double precision :: wallemissivity_west
        
        wallemissivity_west = epsilon_w(1)

    end function wallemissivity_west
    !========================================================================================
    !----------------------------------------------------------
    !     Set the non-uniform emissivity at the east wall
    !---------------------- Description -----------------------
    ! <INPUT>
    ! yc  - y value at the control-volume center
    ! zc  - z value at the control-volume center
    ! j   - y-direction index
    ! k   - z-direction index
    ! <OUTPUT>
    ! wallemissivity_east - East wall emissivity data
    !NOTE: use yc(j) and zc(k) for spatial positions
    !----------------------------------------------------------
    function wallemissivity_east(j,k)
        implicit none
        integer :: j,k
        double precision :: wallemissivity_east
        
        wallemissivity_east = epsilon_w(2)

    end function wallemissivity_east
    !========================================================================================
    !----------------------------------------------------------
    !     Set the non-uniform emissivity at the south wall
    !---------------------- Description -----------------------
    ! <INPUT>
    ! xc  - x value at the control-volume center
    ! zc  - z value at the control-volume center
    ! i   - x-direction index
    ! k   - z-direction index
    ! <OUTPUT>
    ! wallemissivity_south - South wall emissivity data
    !NOTE: use xc(i) and zc(k) for spatial positions
    !----------------------------------------------------------
    function wallemissivity_south(i,k)
        implicit none
        integer :: i,k
        double precision :: wallemissivity_south
        
        wallemissivity_south = epsilon_w(3)

    end function wallemissivity_south
    !========================================================================================
    !----------------------------------------------------------
    !     Set the non-uniform emissivity at the south wall
    !---------------------- Description -----------------------
    ! <INPUT>
    ! xc  - x value at the control-volume center
    ! zc  - z value at the control-volume center
    ! i   - x-direction index
    ! k   - z-direction index
    ! <OUTPUT>
    ! wallemissivity_north - South wall emissivity data
    !NOTE: use xc(i) and zc(k) for spatial positions
    !----------------------------------------------------------
    function wallemissivity_north(i,k)
        implicit none
        integer :: i,k
        double precision :: wallemissivity_north
        
        wallemissivity_north = epsilon_w(4)

    end function wallemissivity_north
    !========================================================================================
    !----------------------------------------------------------
    !    Set the non-uniform emissivity at the bottom wall
    !---------------------- Description -----------------------
    ! <INPUT>
    ! xc  - x value at the control-volume center
    ! yc  - y value at the control-volume center
    ! i   - x-direction index
    ! j   - y-direction index
    ! <OUTPUT>
    ! wallemissivity_bottom - South wall emissivity data
    !NOTE: use xc(i) and yc(j) for spatial positions
    !----------------------------------------------------------
    function wallemissivity_bottom(i,j)
        implicit none
        integer :: i,j
        double precision :: wallemissivity_bottom
        
        wallemissivity_bottom = epsilon_w(5)

    end function wallemissivity_bottom
    !========================================================================================
    !----------------------------------------------------------
    !      Set the non-uniform emissivity at the top wall
    !---------------------- Description -----------------------
    ! <INPUT>
    ! xc  - x value at the control-volume center
    ! yc  - y value at the control-volume center
    ! i   - x-direction index
    ! j   - y-direction index
    ! <OUTPUT>
    ! wallemissivity_top - South wall emissivity data
    !NOTE: use xc(i) and yc(j) for spatial positions
    !----------------------------------------------------------
    function wallemissivity_top(i,j)
        implicit none
        integer :: i,j
        double precision :: wallemissivity_top
        
        wallemissivity_top = epsilon_w(6)

    end function wallemissivity_top
    !========================================================================================
end module functions
