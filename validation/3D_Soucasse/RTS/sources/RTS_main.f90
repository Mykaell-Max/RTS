program RTS

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
    use start_variables
    use output
    use radiation
    use energy
    !========================================================================================
    !-----------------------------------------------------------------------
    !                            Main Aplication
    !------------------------------ Description ----------------------------
    ! <EXTERNAL ROUTINES>
    ! start      - Start and setup all variables
    ! cpu_time   - Returns the elapsed CPU time in seconds
    ! temp_walls - Apply the termal boundary conditions
    ! rtesolve   - Solves the Radiative Transfer Equation
    ! main_save  - Main routine for saving results
    !-----------------------------------------------------------------------
    call start 
    call cpu_time(cpu_time_s)
    if(energy_flag) then
        call energy_equation
    else
        call temp_walls
        call rtesolve
    end if
    call cpu_time(cpu_time_f)
    cpu_time_d =  (cpu_time_f - cpu_time_s)
    if(cpu_time_d < 60.0d0)then
        write(*,'(a24,f5.2,a8)') 'Simulation Completed in ',cpu_time_d,' Seconds'
    else if(cpu_time_d >= 60.0d0)then
        write(*,'(a24,f5.2,a8)') 'Simulation Completed in ',cpu_time_d/60.0d0,' Minutes'
    else if(cpu_time_d >= 3600.0d0)then
        write(*,'(a24,f5.2,a8)') 'Simulation Completed in ',cpu_time_d/3600.0d0,' Hours'
    end if
    call main_save
    !========================================================================================
end program RTS
