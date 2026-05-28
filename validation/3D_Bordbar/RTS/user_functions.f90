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
