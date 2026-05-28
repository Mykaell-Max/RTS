
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

  function exact_temp(t,x,y,z)
    implicit none
    double precision :: t, x, y, z 
    double precision :: exact_temp
    double precision :: Tc, r, func
    logical :: is_boundary
    real(kind=8), parameter :: PI = 4.0d0 * atan(1.0d0)
  
    is_boundary = (x < 1.0d-3 .or. x > 1.0d0 - 1.0d-3) !.or. &
                  ! (y < 1.0d-8 .or. y > 1.0d-5 - 1.0d-8) .or. &
                  ! (z < 1.0d-8 .or. z > 1.0d-5 - 1.0d-8)
    
    if (is_boundary) then
      exact_temp = 400.0d0
    else 
        exact_temp = 400.0d0 + 1400.0d0*(sin(PI*(x/1.0d0)))**2.0d0
    end if

    print *, x, y, z, exact_temp

  
  end function exact_temp

function exact_p(t,x,y,z)
  implicit none
  double precision :: t, x, y, z
  double precision :: exact_p

  exact_p = 101325
end function exact_p


function exact_yk(t,x,y,z)
  use data_grid
  use data_concentrations
  use eul_cantera
  implicit none
  double precision :: t, x, y, z 
  double precision :: exact_yk
  
  exact_yk=EnvironmentVec(cont_esc_transp)
end function exact_yk
  

