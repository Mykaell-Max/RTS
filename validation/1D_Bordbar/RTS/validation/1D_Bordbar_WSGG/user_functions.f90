    !========================================================================================
    function Temp_field(i,j,k)
        !----------------------------------------------------------
        !             Temperature space-function
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
        ! Temp_field - Temperature data
        !NOTE: use xc(i), yc(j) and zc(k) for spatial positions
        !----------------------------------------------------------
        implicit none
        integer, intent(in) :: i,j,k
        double precision :: Temp_field

        Temp_field = 400.0d0 + 1400.0d0*(sin(PI*(xc(i)/lx)))**2.0d0

    end function Temp_field
    !========================================================================================
    function Y_CO2field(i,j,k)
        !----------------------------------------------------------
        !             CO2 mole fraction space-function
        ! Made by:     G.S.Rodrigues (July 2021)
        ! Modified by: 
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
        implicit none
        integer, intent(in) :: i,j,k
        double precision :: Y_CO2field
        Y_CO2field = 1.0d0 - 1.0d-4 + (1 - 1.0d-4)*(sin(PI*(xc(i)/lx)))**2.0d0
        
    end function Y_CO2field
    !========================================================================================
    function Y_H2Ofield(i,j,k)
        !----------------------------------------------------------
        !             H2O mole fraction space-function
        ! Made by:     G.S.Rodrigues (July 2021)
        ! Modified by: 
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
        implicit none
        integer, intent(in) :: i,j,k
        double precision :: Y_H2Ofield
        Y_H2Ofield = 1.0d-4 + (1 - 1.0d-4)*(sin(PI*(xc(i)/lx)))**2.0d0
        
    end function Y_H2Ofield
    !========================================================================================
