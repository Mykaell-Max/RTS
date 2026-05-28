module energy

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
    use output
    use solvers
    use radiation
    use boundary_conditions
    implicit none
    contains
    !========================================================================================
    !-----------------------------------------------------------------------
    !               Solves the conduction energy equation
    !------------------------------ Description ----------------------------
    ! <INPUT>
    ! T_energy- Initial guess or data from previous call
    ! G       - Initial guess or data from previous call
    ! a_bxT   - a boundary coefficient x-direction walls
    ! a_byT   - a boundary coefficient y-direction walls
    ! a_bzT   - a boundary coefficient z-direction walls
    ! b_bxT   - b boundary coefficient x-direction walls
    ! b_byT   - b boundary coefficient y-direction walls
    ! b_bzT   - b boundary coefficient z-direction walls
    ! g_bxT   - g boundary coefficient x-direction walls
    ! g_byT   - g boundary coefficient y-direction walls
    ! g_bzT   - g boundary coefficient z-direction walls
    ! rho     - Fluid density 
    ! K_term  - Thermal conductivity
    ! Cp      - Specific heat
    ! CTS     - Save the data each CTS time step
    ! <INTERNAL>
    ! Aw      - West   difusion based coefficient for matix A
    ! Ae      - East   difusion based coefficient for matix A
    ! As      - South  difusion based coefficient for matix A
    ! An      - North  difusion based coefficient for matix A
    ! Ab      - Bottom difusion based coefficient for matix A
    ! At      - Top    difusion based coefficient for matix A
    ! Ap      - Center difusion based coefficient for matix A
    ! A_RHS   - Auxiliary term thad must be added to Ap
    ! RHS     - Right hand side for matix B
    ! A_dt    - Auxiliary transient coefficient
    ! S_rad   - Divergent radiative flux
    ! res_tm  - Residue for energy equation
    ! res_rad - Residue for radiative equation
    ! omega   - The over-ralaxation factor
    ! dt      - time step
    ! time    - Simulated time
    ! f_time  - Final time
    ! <OUTPUT>
    ! T_energy- Temperature for all nodes
    ! G       - Radiant intensity (i,j,k)
    ! ITP_TM  - Number of iterations required in energy equation
    ! ITP_G   - Number of iterations required in P1 equation
    ! <EXTERNAL ROUTINES>
    ! elliptic_coefficients - Calculates the elliptic equation coefficients 
    ! bound_temp       - Termal boundary conditions
    ! Boundary         - Applies the boundary conditions to the coefficients
    ! SOR              - Solves the linear system via SOR
    ! wall_properties  - Calculates the value of the variables at the walls
    ! rtesolve        - Solves the Radiative Transfer Equation
    ! VTK_save         - Saves the selected variables in a vtk file
    ! dat_save         - Saves the selected variables in a dat file
    !-----------------------------------------------------------------------
    subroutine energy_equation
        double precision :: Aw(nxi,nyi,nzi),As(nxi,nyi,nzi),Ab(nxi,nyi,nzi),    &
                            Ae(nxi,nyi,nzi),An(nxi,nyi,nzi),At(nxi,nyi,nzi),    &
                            Ap(nxi,nyi,nzi),RHS(nxi,nyi,nzi),A_RHS(nxi,nyi,nzi),&
                            A_dt,omega,dx,dy,dz
        !--------------- the over-ralaxation factor --------------
        !omega = 2.0d0/(1.0d0 + sin(3.14d0*(1.0d0/(sqrt((nx*nx+ny*ny+nz*nz)/3.0d0)+1.0d0))))
        omega = 1.9d0
        !----------------- time step calculation -----------------
        dx = maxval(dxp);dy = maxval(dyp);dz = maxval(dzp)
        dt = 2.0d0*(1.0d0/(dx*dx) + 1.0d0/(dy*dy) + 1.0d0/(dz*dz))
        dt = 1.0d0/(dt*(maxval(K_term)/(rho*Cp)))
        !---------------------------------------------------------
        if(trans_flag) then
            A_RHS = 1.0d0
            A_dt = dt/(rho*Cp)
        else
            A_dt = 1.0d0
            A_RHS = 0.0d0
        end if
        call elliptic_coefficients(K_term,A_dt,A_RHS,Ab,At,Aw,Ae,As,An,Ap)
        call bound_temp
         write(*,'(a)') '======== Energy Equation ======='
        do ITP = 1,ITMAX
            if(trans_flag) then
                RHS = -(T_energy + S_rad*(dt/(rho*Cp)))
            else
                RHS = 0.0d0
            end if
            time = dt*ITP
            call Boundary(Ab,At,Aw,Ae,As,An,Ap,RHS,a_bxT,b_bxT,g_bxT,a_byT,b_byT,g_byT,a_bzT,b_bzT,g_bzT,K_term)
            call SOR(T_energy,Ab,At,Aw,Ae,As,An,Ap,RHS,ITP_TM,res_tm,ene_tol,omega)
            call wall_properties(T_energy,a_bxT,b_bxT,g_bxT,a_byT,b_byT,g_byT,a_bzT,b_bzT,g_bzT)
            write(*,'(a25,I6)') 'Residuals for iteration: ',ITP
            write(*,'(a3,E10.3,a7,I5)') 'T:',res_tm,' N_ITE:',ITP_TM
            if(rte_flag .eqv. .false.) then
                call rtesolve
            end if
            if(trans_flag) write(*,'(a14,f9.3,a8)')'Physical time:',time, ' seconds'
            write(*,'(a)') '================================'
            if(mod(ITP,CTS) == 0 .or. ITP == 1) then
                call VTK_save
                if(trans_flag .eqv. .false.) exit
            end if
            if(f_time <= time) then
                call VTK_save
                exit
            end if
        end do
        call dat_save
    end subroutine energy_equation
    !========================================================================================
end module energy
