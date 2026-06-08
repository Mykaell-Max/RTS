module global

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
    implicit none

    !<><><><><><><><><><><><><><><><><><> Constant Parameters  <><><><><><><><><><><><><><><><><>
    double precision :: &
    small,     & !A small number
    big,       & !A big number
    boltz,     & !stefan boltzmann constant
    PI,        & !PI constant
    PIBY2,     & !PI/2
    PI32,      & !3PI/2
    PI4          !4PI
    
    !<><><><><><><><><><><><><><><><><><> Solver Variables <><><><><><><><><><><><><><><><><><><>
    logical :: &
    stop_flag    !Abort code execution
    
    integer :: & 
    ITMAX,     & !maximum number of iterations
    ITP,       & !number of iterations
    CTS,       & !save the data each CTS time step
    ITP_TM,    & !number of iterations required to solve the energy equation
    ITP_G        !number of iterations required to solve the radiative equation
    
    double precision :: &
    dt,        & !time step (s)
    rad_tol,   & !radiation convergence criteria
    ene_tol,   & !energy convergence criteria
    res_tm,    & !residue for energy equation
    res_rad,   & !residue for radiative equation
    time,      & !physical time (s)    
    f_time,    & !final physical time (s)
    cpu_time_s,& !CPU simulation initial time
    cpu_time_f,& !CPU simulation final time
    cpu_time_d   !CPU duration time
    
    !<><><><><><><><><><><><><><><><><><> Output Variables <><><><><><><><><><><><><><><><><><><>
    logical :: &
    csv_flag,  & !save results to csv
    vtk_flag,  & !save results to vtk
    dat_flag,  & !save results to dat
    xyz_save,  & !save mesh positions
    cstm_save, & !save wall results to vtk
    nuss_flag    !Nusselt module flag
    
    logical, dimension(13)::&
    plane_flag,& !Plane ID flags 1-3 (ZY, XZ, XY)
    field_flag,& !It indicates whether the input field is uniform or not
    Rec_flag,  & !It indicates whether the variable will be Recorded or not
    nuffwall,  & !Non-uniform field flag on the wall
    SBCwall      !Symmetry boundary conditions flags
    
    logical, dimension(6,5)::&
    crec_flag    !Walls custom variables save flags
    
    !<><><><><><><><><><><><><><><><><><><> Mesh Variables <><><><><><><><><><><><><><><><><><><>
    logical :: &
    path_flag    !use Paths
    
    integer :: & 
    DIMEN,     & !number of dimensions
    nx,        & !number of control volumes in x-direction
    ny,        & !number of control volumes in y-direction
    nz,        & !number of control volumes in z-direction
    nxp,       & !auxiliary index (nx + 1)
    nyp,       & !auxiliary index (ny + 1)
    nzp,       & !auxiliary index (nz + 1)
    nxi,       & !number of information stored in the x-direction (nx + 2)
    nyi,       & !number of information stored in the y-direction (ny + 2)
    nzi,       & !number of information stored in the z-direction (nz + 2)
    npath        !total number of paths
    
    double precision :: &
    lx,        & !domain length (m)
    ly,        & !domain height (m)
    lz           !domain width  (m)
    
    double precision, allocatable, dimension(:):: &
    x,         & !x-direction grid (m) - value at the control-volume face
    y,         & !y-direction grid (m) - value at the control-volume face
    z,         & !z-direction grid (m) - value at the control-volume face
    xc,        & !x-direction grid (m) - value at the control-volume center
    yc,        & !y-direction grid (m) - value at the control-volume center
    zc,        & !z-direction grid (m) - value at the control-volume center
    dxp,       & !x-direction length of control-volume
    dyp,       & !y-direction height of control-volume
    dzp          !z-direction width  of control-volume
    
    double precision, allocatable :: &
    Lpath(:,:),& !Stores the paths coordinates
    vol(:,:,:)   !cell volume (m^3)
    
    !<><><><><><><><><><><><><><><><><><><> Core Variables <><><><><><><><><><><><><><><><><><><>
    double precision :: &
    T_energyc, & !constant temperature for all nodes (K)
    Gc,        & !constant incident radiation for all nodes (w/m^2)
    K_termc,   & !constant thermal conductivity (W/m*K)
    k_rad,     & !constant for absorption coefficient
    sig_rad      !constant for scattering coefficient
    
    double precision, allocatable, dimension(:,:,:):: &
    T_energy,  & !temperature for all nodes (K)               ID 1
    G,         & !incident radiation for all nodes (w/m^2)    ID 2
    S_rad,     & !divergent radiative flux (w/m^3)            ID 3
    IBlack,    & !blackbody radiation intensity               ID 4
    cappa,     & !absorption coefficient for all nodes (m^-1) ID 5
    sigma,     & !scattering coefficient for all nodes (m^-1) ID 6
    beta,      & !extinction coefficient for all nodes (m^-1) ID 7
    K_term,    & !thermal conductivity (W/m*K)                ID 8
    Q_radw,    & !walls radiative fluxes (w/m^2)              ID 10
    epsilon_rad  !emissivities at the walls and domain(gas)   ID 9
    
    !Boundary Variables
    double precision, dimension(6)::&
    wall_temp, & !wall temperature
    epsilon_w    !radiation emissivities
    
    !<><><><><><><><><><><><><><><><><><> Energy Variables <><><><><><><><><><><><><><><><><><><>
    logical :: &
    rte_flag,  & !use radiation module
    trans_flag,& !transient simulation
    energy_flag  !solve energy equation
    
    integer, dimension(6)::&
    BC_flags_t   !Boundary flags temperature
    
    double precision :: &
    rho,       & !constant fluid density (kg/m^3)
    Cp,        & !constant specific heat (J/kg*K)
    T_inf        !room temperature (K)
    
    !Boundary Variables
    double precision, dimension(6)::&
    h_conv       !convective heat transfer coefficient 
    
    double precision, allocatable, dimension(:,:,:):: &
    !Temperature
    a_bxT,     & !a boundary coefficient x-direction walls
    a_byT,     & !a boundary coefficient y-direction walls
    a_bzT,     & !a boundary coefficient z-direction walls
    b_bxT,     & !b boundary coefficient x-direction walls
    b_byT,     & !b boundary coefficient y-direction walls
    b_bzT,     & !b boundary coefficient z-direction walls
    g_bxT,     & !g boundary coefficient x-direction walls
    g_byT,     & !g boundary coefficient y-direction walls
    g_bzT        !g boundary coefficient z-direction walls

    !<><><><><><><><><><><><><><><><> Common Radiation Variables <><><><><><><><><><><><><><><><>
    logical :: &
    aniso_flag   !anisotropic phase function
    
    character(len=70) :: &
    schm_model,& !Select the spatial scheme used to solve the RTE
    SPF_func,  & !Scattering phase function
    rad_model, & !Selected the model for quadrature sets and weigths
    scat_model   !Select the phase function model
    
    integer :: & 
    model_stc, & !Selected method to resolve radiation
    SPF_flag,  & !Selected scattering phase function
    scat_flag, & !Select the phase function
    scheme_flag  !Select the spatial scheme used to solve the RTE
    
    double precision :: &
    C_rad,     & !Linear anisotropic phase function coefficient
    alp_r        !Radiation scheme variable
    
    !<><><><><><><><><><><><><><><><><><> P1 Method Variables <><><><><><><><><><><><><><><><><><>
    integer, dimension(6)::&
    BC_flags_G   !Boundary flags incident radiation
    
    double precision, allocatable, dimension(:,:,:):: &
    GAMA,      & !gamma coefficient (m) (P1 Method)
    a_bxG,     & !a boundary coefficient x-direction walls
    a_byG,     & !a boundary coefficient y-direction walls
    a_bzG,     & !a boundary coefficient z-direction walls
    b_bxG,     & !b boundary coefficient x-direction walls
    b_byG,     & !b boundary coefficient y-direction walls
    b_bzG,     & !b boundary coefficient z-direction walls
    g_bxG,     & !g boundary coefficient x-direction walls
    g_byG,     & !g boundary coefficient y-direction walls
    g_bzG        !g boundary coefficient z-direction walls
    
    !<><><><><><><><><><><><><><><><> Finite Angle Method Variables  <><><><><><><><><><><><><><><>
    integer :: &
    nt,        & !control volumes in theta-direction
    np,        & !control volumes in phi-direction
    nsub,      & !number of angular subdivisions for the phase function
    P2,        & !P2 First  azimuthal region index (0 < phi < PIBY2)
    P3,        & !P3 Second azimuthal region index (PIBY2 < phi < PI)
    P4,        & !P4 Third  azimuthal region index (PI < phi < PI32)
    T2           !T2 First polar region index (0 < theta < PIBY2)
    
    double precision, allocatable, dimension(:):: &
    theta,     & !theta-direction grid - value at the control-angle face
    phi,       & !phi-direction   grid - value at the control-angle face
    A_phase      !expansion coefficients array 

    double precision, allocatable, dimension(:,:):: &
    dcx,       & !integral over the solid angle in the x-direction
    dcy,       & !integral over the solid angle in the y-direction
    dcz,       & !integral over the solid angle in the z-direction
    dco          !integral over the solid angle
    
    double precision, allocatable, dimension(:,:,:,:):: &
    phase_f,   & !scattering phase function for FAM
    Ax,        & !x-direction coefficients for FAM
    Ay,        & !y-direction coefficients for FAM
    Az           !z-direction coefficients for FAM
    
    double precision, allocatable, dimension(:,:,:,:,:):: &
    IG,        & !radiant intensity (i,j,k,l,m)
    volom        !volume of the spatial-angular cell (i,j,k,l,m)
    
    !<><><><><><><><><><><><><><>  Discrete Ordinates Method Variables <><><><><><><><><><><><><><>
    character(len=70) :: &
    quadrature_DOM !selected the model for quadrature sets and weigths
    
    integer :: &
    QM_stc,    & !selected the model for quadrature sets and weigths
    N_quad,    & !quadrature order
    nq           !number of quadratures per octant
    
    double precision, allocatable, dimension(:):: &
    mux,       & !x ordinate direction
    etay,      & !y ordinate direction
    xiz,       & !z ordinate direction
    Wq           !quadrature weigth
    
    double precision, allocatable, dimension(:,:):: &
    phase_d      !scattering phase function for DOM
    
    double precision, allocatable, dimension(:,:,:):: &
    Axd,       & !x-direction coefficients for DOM
    Ayd,       & !y-direction coefficients for DOM
    Azd          !z-direction coefficients for DOM
    
    !<><><><><><><><><><><><><> Radiative Properties of Combustion Gases <><><><><><><><><><><><><>
    logical :: &
    wsgg_model,& !use WSGG module
    gray_model,& !use gray module
    gas_prop,  & !gas properties flag
    const_conc,& !constant gas concentration flag
    nongray_flag !non-gray radiation flag
    
    character(len=70) :: &
    absrp_model !absorption model selected
    
    integer :: &
    nsb,      & !number of spectral bands (nsb = 1 Gray , nsb = 5 WSGG)
    absrp_stc   !absorption model selected
    
    double precision :: &
    P_g,       & !Global total pressure of the mixture
    XCO2_g,    & !Global CO2 mole fraction
    XH2O_g       !Global H2O mole fraction
    
    !Spectral Bands Variables
    double precision,allocatable,dimension(:,:,:) :: &
    P_species    !Field species partial pressure (i,j,k,specie) ID 11
    
    double precision,allocatable,dimension(:,:,:,:) :: &
    Y_species, & !Field species mole fraction    (i,j,k,specie) ID 12/13 
    cappaBND,  & !cappaBND- cappa for all bands    (i,j,k,IBND)
    BBF,       & !blackbody radiation fraction     (i,j,k,IBND)
    GBND         !incident radiation for all bands (i,j,k,IBND)
    
    double precision, allocatable, dimension(:,:,:,:,:,:):: &
    IGBND        !radiant intensity for all bands (i,j,k,l,m,IBND)
    
    !<><><><><><><><><><><><><><> Weighted Sum of Gray Gases Variables <><><><><><><><><><><><><><>
    integer :: &
    nGi,       & !WSGG polynomial order
    nGj          !WSGG polynomial order of the weight coefficient
    
    double precision :: &
    Lm           !Averaged mean beam length
    
    double precision,allocatable,dimension(:,:)   :: &
    d_wgh        !emissivity polynomial coefficients of the WSGG
    !NOTE: b_wgh (WSGG b coefficient) was moved to a local array inside the
    !      WSGG subroutine (RTS_absorption.f90). It is pure scratch space passed
    !      between WSGG_weights and a_weights within a single WSGG call. Keeping
    !      it as a shared module variable created a race condition when WSGG is
    !      called concurrently by multiple threads (radiative_properties loop).
    
    double precision,allocatable,dimension(:,:,:) :: &
    c_wgh        !asorptivity polynomial coefficients of the WSGG
    
    !<><><><><><><><><><><><><><><><><>  Slices module variables <><><><><><><><><><><><><><><><><>
    logical :: &
    slice_flag   !use slices save module
    
    character,allocatable,dimension(:) :: &
    slc_axis     !Slices axis array
    
    integer :: &
    OPslice,   & !Slicing method
    nslices      !Number of slices 
    
    integer,allocatable,dimension(:) :: &
    slc_perc,  & !Slices percent array
    slc_node,  & !Slices node array
    slc_ID       !Slices variable ID array
    
    double precision,allocatable,dimension(:) :: &
    slc_point    !Slices point array
    
    !<><><><><><><><><><><><><><><><><>  Lineup module variables <><><><><><><><><><><><><><><><><>
    logical :: &
    lnp_flag     !use lineup save module
    
    character,allocatable,dimension(:) :: &
    lnp_axis     !Slices axis array
    
    integer :: &
    OPlnp,     & !Lineup method
    nlnp         !Number of lineups
    
    integer,allocatable :: &
    lnp_node(:,:),  & !Lineup node array
    lnp_ID(:)         !Lineup variable ID array
    
    double precision ,allocatable, dimension(:,:) :: &
    lnp_point    !Lineup point array
    
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
end module global
