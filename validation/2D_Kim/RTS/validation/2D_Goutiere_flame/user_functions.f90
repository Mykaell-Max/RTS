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
        double precision :: Temp_field,Yo,aux
        
        Yo = abs(1.0d0 - 4.0d0*yc(j))
        aux = (1.0d0 - 3.0d0*Yo**2.0d0 + 2.0d0*Yo**3.0d0)
        if(xc(i) < 0.1d0)then
            Temp_field = (14000.0d0*xc(i) - 400.0d0)*aux + 800.0d0
        elseif(xc(i) >= 0.1d0)then
            Temp_field = -(1.0d4/9.0d0)*(xc(i) - 1.0d0)*aux + 800.0d0
        end if

    end function Temp_field
    !========================================================================================
