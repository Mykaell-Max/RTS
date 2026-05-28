    !========================================================================================
    function Temp_field(i,j,k)
        !----------------------------------------------------------
        !             Temperature space-function
        ! Made by:     G.S.Rodrigues (November 2020)
        ! Modified by: G.S.Rodrigues (May 2021)
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
        double precision :: Temp_field,Tc
        
        Tc = (xc(i) - 0.25d0)**2.0d0 + (yc(j) - 0.25d0)**2.0d0 + (zc(k) - 0.25d0)**2.0d0
        Temp_field = exp(-Tc)*10.0d0 + 300.0d0

    end function Temp_field
    !========================================================================================
