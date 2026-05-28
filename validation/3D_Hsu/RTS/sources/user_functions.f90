    !========================================================================================
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
