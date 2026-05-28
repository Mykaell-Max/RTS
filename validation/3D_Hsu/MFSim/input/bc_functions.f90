
! Implement in this file any of the used boundary condition
! functions (use one of the following function names):
!
! - pres property:
!    exact_p
!    dpdx
!    dpdy
!    dpdz
!
! - velu, vu_1, vu_0 properties:
!    exact_u
!    dudx
!    dudy
!    dudz
!    dudt
!
! - velv, vv_1, vv_0 properties:
!    exact_v
!    dvdx
!    dvdy
!    dvdz
!    dvdt
!
! - velw, vw_1, vw_0 properties:
!    exact_w
!    dwdx
!    dwdy
!    dwdz
!    dwdt
!
! - vof property:
!    exact_vof
!
! - sca property:
!    exact_sca
!    dscadx
!    dscady
!    dscadz
!    dscadt
!
! - vari property:
!    exact_var
!    dvardx
!    dvardy
!    dvardz
!    dvardt
!
! - temp property:
!    exact_temp
!    dtempdx
!    dtempdy
!    dtempdz
!    dtempdt
!
! - cond property:
!    exact_cond
!
! - cepe property:
!    exact_cepe
!
! - dens, dens_0, dens_1 properties:
!    exact_rho (if not defined, default is taken from input file)
!
! - visc, visc_0, visc_1 properties:
!    exact_mu (if not defined, default is taken from input file)
!
! - tur1 property:
!    exact_tur1
!    dtur1dx
!    dtur1dy
!    dtur1dz
!    dtur1dt
!
! - tur2 property:
!    exact_tur2
!    dtur2dx
!    dtur2dy
!    dtur2dz
!    dtur2dt
!
!         !!!!!!!!!!!!!!!!!!!!!!
!         !!! Your Code Here !!!
!         !!!!!!!!!!!!!!!!!!!!!!
! end function
  !************************************************************** 
  !> \brief Define a exact value of temperature
  !<

   function cappa_field(t,x,y,z)
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
    double precision :: x,y,z,t,cappa_field,betaux
    
    betaux = 0.9d0*(1.0d0 - 2.0d0*abs(x - 0.5d0))*(1.0d0 - 2.0d0*abs(y - 0.5d0))
    betaux = betaux*(1.0d0 - 2.0d0*abs(z - 0.5d0)) + 0.1d0
    cappa_field = betaux
    !cappa_field = 0.1d0*betaux

end function cappa_field
!========================================================================================
function sigma_rad_field(t,x,y,z)
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
    double precision :: t,x,y,z,sigma_rad_field,betaux
            
    betaux = 0.9d0*(1.0d0 - 2.0d0*abs(x - 0.5d0))*(1.0d0 - 2.0d0*abs(y - 0.5d0))
    betaux = betaux*(1.0d0 - 2.0d0*abs(z - 0.5d0)) + 0.1d0
    sigma_rad_field = 0.9d0*betaux
    
end function sigma_rad_field

  function exact_temp(t,x,y,z)
    implicit none
    double precision :: t, x, y, z, delta 
    double precision :: exact_temp
    logical :: is_boundary
    
    delta = 1.0d-3

    is_boundary = (x < delta .or. x > 1.0d0 - delta) .or. &
                  (y < delta .or. y > 1.0d0 - delta) .or. &
                  (z < delta .or. z > 1.0d0 - delta)
    
    if (is_boundary) then
      exact_temp = 0.0d0
    else 
      exact_temp = 64.8032d0
    end if
  
  end function exact_temp
  

