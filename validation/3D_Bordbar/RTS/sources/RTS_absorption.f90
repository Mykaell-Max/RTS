module absorptiondata

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
    implicit none
    contains
    !========================================================================================
    !                               Gas properties subroutines
    !========================================================================================
    !----------------------------------------------------------
    !      Computes the radiative properties of the medium
    !---------------------- Description -----------------------
    ! <INPUT>
    ! T_energy - Temperature (K)
    ! nsb      - number of gray gases plus 1 for the clear gas
    ! <INTERNAL>
    ! K_wgh    - local absorption coefficients of the gray gases
    ! a_wgh    - local array of weights (sum a_wgh = 1)
    ! <OUTPUT>
    ! cappa    - absorption coefficients of the gray gases (1/m)
    ! BBF      - blackbody radiation fraction array
    ! <EXTERNAL ROUTINES>
    ! WSGG - Compute the weighted sum of gray gases model coefficients
    !----------------------------------------------------------
    subroutine radiative_properties()
        integer :: i,j,k,IBND
        double precision :: a_wgh(nsb),K_wgh(nsb)
        !-------------------------------------------------------------------
        !                   Non-gray radiation properties
        !-------------------------------------------------------------------
        if(nongray_flag)then 
            if(wsgg_model)then
                write(*,'(a)') 'Absorption Model: Non-Gray WSGG '
                do k=2,nzp
                    do j=2,nyp
                        do i=2,nxp
                            K_wgh = 0.0d0; a_wgh = 0.0d0
                            call WSGG(T_energy(i,j,k),P_species(i,j,k),Y_species(i,j,k,1),Y_species(i,j,k,2),a_wgh,K_wgh)
                            do IBND=1,nsb
                                cappaBND(i,j,k,IBND) = K_wgh(IBND)
                                BBF(i,j,k,IBND) = a_wgh(IBND)
                            end do
                        end do
                    end do
                end do
            end if
        !-------------------------------------------------------------------
        !                     Gray radiation properties
        !-------------------------------------------------------------------
        else
            if(wsgg_model)then
                write(*,'(a)') ' Absorption Model: Gray WSGG '
                do k=2,nzp
                    do j=2,nyp
                        do i=2,nxp
                            K_wgh = 0.0d0; a_wgh = 0.0d0
                            call WSGG(T_energy(i,j,k),P_species(i,j,k),Y_species(i,j,k,1),Y_species(i,j,k,2),a_wgh,K_wgh)
                            cappa(i,j,k) = K_avg(P_species(i,j,k),Y_species(i,j,k,1),Y_species(i,j,k,2),K_wgh,a_wgh)
                        end do
                    end do
                end do
            else if(gray_model)then
                write(*,'(a)') ' Absorption Model: Gray Gas'
                do k=2,nzp
                    do j=2,nyp
                        do i=2,nxp
                            cappa(i,j,k) = gray_gas(T_energy(i,j,k),P_species(i,j,k),Y_species(i,j,k,1),Y_species(i,j,k,2))
                        end do
                    end do
                end do
            end if
        end if
        !-------------------------------------------------------------------
    end subroutine radiative_properties
    !========================================================================================
    !----------------------------------------------------------
    ! Compute the weighted sum of gray gases model coefficients
    !---------------------- Description -----------------------
    !  For more details see "A line by line based weighted sum 
    !   of gray gases model for inhomogeneous CO2–H2O mixture  
    !     in oxy-fired combustion" by Bordbar et al.(2014)
    ! doi: https://doi.org/10.1016/j.combustflame.2014.03.013
    !----------------------------------------------------------
    ! <INPUT>
    ! T      - Temperature (K)
    ! P      - Total pressure of the mixture: PCO2 + PH2O (atm)
    ! XCO2   - CO2 mole fraction
    ! XH2O   - H2O mole fraction
    ! <INTERNAL>
    ! Mr     - Molar fraction ratio
    ! nsb    - number of gray gases plus 1 for the clear gas
    ! nGi    - polynomial order
    ! nGj    - polynomial order of the weight coefficient
    ! c_wgh  - asorptivity polynomial coefficients of the WSGG
    ! b_wgh  - b asorptivity polynomial coefficient of the WSGG
    ! d_wgh  - emissivity polynomial coefficients of the WSGG
    ! <OUTPUT>
    ! K_wgh  - absorption coefficients of the gray gases (1/m)
    ! a_wgh  - local array of weights (sum a_wgh = 1)
    ! <EXTERNAL ROUTINES>
    ! WSGG_MolRatio - Compute the molar ratio of a CO2-H2O mixture
    ! WSGG_weights  - Compute the asorptivity polynomial coefficient
    ! WSGG_kappa    - Compute the gray gas absorption coefficients
    ! a_weights     - Compute the weighting coefficients for the WSGG        
    !----------------------------------------------------------
    subroutine WSGG(T,P,XCO2,XH2O,a_wgh,K_wgh)
        implicit none
        double precision :: T,P,XCO2,XH2O,Mr,a_wgh(nsb),K_wgh(nsb)
        call WSGG_MolRatio(Mr,XCO2,XH2O)
        call WSGG_weights(Mr,b_wgh)
        call a_weights(T,a_wgh,b_wgh)
        call WSGG_kappa(P,XCO2,XH2O,Mr,K_wgh)
    end subroutine WSGG
    !========================================================================================
    !----------------------------------------------------------
    ! Compute the molar ratio of CO2-H2O mixture
    !---------------------- Description -----------------------
    ! <INPUT>
    ! XCO2   - CO2 mole fraction
    ! XH2O   - H2O mole fraction
    ! <OUTPUT>
    ! Mr     - Molar fraction ratio
    !----------------------------------------------------------
    subroutine WSGG_MolRatio(Mr,XCO2,XH2O)
        double precision :: XCO2,XH2O,Mr
        if(XCO2 == 0.0d0)then
            Mr = 0.0d0
        else
            Mr = XH2O/XCO2      !CO2-H2O Mixture
        end if
        !Sanity checks
        if(Mr < 0.01d0)then     !Only CO2
            Mr = 0.01d0
        else if(Mr > 4.0d0)then !Only H2O
            Mr = 4.0d0
        end if
    end subroutine WSGG_MolRatio
    !========================================================================================
    !----------------------------------------------------------
    ! Compute the asorptivity polynomial coefficient of the WSGG
    !---------------------- Description -----------------------
    !  For more details see: "An extended weighted-sum-of-gray
    !  gases model to account for all CO2 ‐ H2O molar fraction 
    !   ratios in thermal radiation " by Bordbar et al.(2020)
    ! doi: https://doi.org/10.1016/j.icheatmasstransfer.2019.104400
    !----------------------------------------------------------
    ! <INPUT>
    ! nsb    - number of gray gases plus 1 for the clear gas
    ! nGi    - polynomial order
    ! nGj    - polynomial order of the weight coefficient
    ! c_wgh  - asorptivity polynomial coefficients of the WSGG
    ! Mr     - Molar fraction ratio
    ! <INTERNAL>
    ! sum_b  - Summation auxiliary variable
    ! <OUTPUT>
    ! b      - b asorptivity polynomial coefficient of the WSGG
    !----------------------------------------------------------
    subroutine WSGG_weights(Mr,b)
        implicit none
        integer :: i,j,k
        double precision :: Mr,sum_b,b(nGi,nGj)
        if(Mr <= 0.01d0)then               !Only CO2
            b(1,1:5) = (/ 8.495135d-1, -1.496812d+0,  1.361406d+0, -5.551699d-1,  8.076589d-2 /)
            b(2,1:5) = (/-1.103102d-1,  9.363958d-1, -1.250799d+0,  6.527827d-1, -1.206959d-1 /)
            b(3,1:5) = (/ 1.731716d-1, -5.174223d-1,  8.256840d-1, -4.998864d-1,  1.008743d-1 /)
            b(4,1:5) = (/ 3.995426d-2,  1.423006d-1, -1.649481d-1,  5.140768d-2, -3.497246d-3 /)
        else if(Mr >= 4.0d0)then           !Only H2O
            b(1,1:5) = (/ 6.670204d-1, -1.228413d+0,  1.428908d+0, -6.267906d-1,  9.628539d-2 /)
            b(2,1:5) = (/ 2.343433d-1, -3.192256d-1,  8.867348d-1, -5.927787d-1,  1.185824d-1 /)
            b(3,1:5) = (/-1.793041d-1,  1.683454d+0, -2.136989d+0,  1.020422d+0, -1.723960d-1 /)
            b(4,1:5) = (/ 3.455969d-1, -7.510442d-1,  6.313180d-1, -2.416500d-1,  3.530972d-2 /)
        else                               !CO2-H2O Mixture
            do i=1,nGi
                do j=1,nGj
                    sum_b = 0.0d0
                    do k=1,nsb
                        sum_b = sum_b + c_wgh(i,j,k)*Mr**(k-1)
                    end do
                    b(i,j) = sum_b 
                end do
            end do
        end if
    end subroutine WSGG_weights
    !========================================================================================
    !----------------------------------------------------------
    !    Compute the temperature weighting factors for the 
    !                  gray gases coefficients
    !---------------------- Description -----------------------
    ! <INPUT>
    ! T      - Temperature (K)
    ! nsb    - number of gray gases plus 1 for the clear gas
    ! nGi    - polynomial order
    ! nGj    - polynomial order of the weight coefficient
    ! <INTERNAL>
    ! b      - b polynomial coefficient of the WSGG
    ! sum_a  - Summation auxiliary variable
    ! <OUTPUT>
    ! a      - array of weights (sum a = 1)
    !----------------------------------------------------------
    subroutine a_weights(T,a,b)
        integer :: i,j
        double precision :: T,sum_a,a(nsb),b(nGi,nGj)
        !Computing the weighting factors for the gray gases
        a(1) = 1.0d0
        do i=2,nsb
            sum_a = 0.0d0
            do j=1,nGj
                sum_a = sum_a + b(i-1,j)*(T/1200.0d0)**(j-1)
            end do
            a(i) = sum_a
            a(1) = a(1) - sum_a
        end do
    end subroutine a_weights
    !========================================================================================
    !----------------------------------------------------------
    ! Compute the gray gas absorption coefficients
    !---------------------- Description -----------------------
    !  For more details see: "An extended weighted-sum-of-gray
    !  gases model to account for all CO2 ‐ H2O molar fraction 
    !   ratios in thermal radiation " by Bordbar et al.(2020)
    ! doi: https://doi.org/10.1016/j.icheatmasstransfer.2019.104400
    !----------------------------------------------------------
    ! <INPUT>
    ! P      - Total pressure of the mixture: PCO2 + PH2O (atm)
    ! XCO2   - CO2 mole fraction
    ! XH2O   - H2O mole fraction
    ! nsb    - number of gray gases plus 1 for the clear gas
    ! d_wgh  - emissivity polynomial coefficients of the WSGG
    ! Mr     - Molar fraction ratio
    ! <INTERNAL>
    ! sum_k  - Summation auxiliary variable
    ! aux    - Auxiliary variable
    ! <OUTPUT>
    ! K_wgh  - absorption coefficients of the gray gases (1/m)
    !----------------------------------------------------------
    subroutine WSGG_kappa(P,XCO2,XH2O,Mr,K_wgh)
        implicit none
        integer :: i,k
        double precision :: P,XCO2,XH2O,Mr,sum_k,K_wgh(nsb),aux
        aux = (P/101325)*(XH2O + XCO2)
        if(Mr <= 0.01d0)then
            K_wgh = (/0.0d0, 3.272772d-2, 4.229655d-1, 4.905367d0, 1.08544d2/) !Only CO2
        else if(Mr >= 4.0d0)then
            K_wgh = (/0.0d0, 8.047859d-2, 9.557208d-1, 8.005283d0, 7.613186d1/) !Only H2O
        else
            do i=2,nsb
                sum_k = 0.0d0
                do k=1,nsb
                    sum_k = sum_k + d_wgh(i-1,k)*Mr**(k-1)
                end do
                K_wgh(i) = sum_k*aux !CO2-H2O mixture
            end do
        end if
    end subroutine WSGG_kappa
    !========================================================================================
    !----------------------------------------------------------
    !  Loads the coefficients C and D used in the WSGG for a
    !    CO2-H2O mixture according to Bordbar et al.(2014)
    !---------------------- Description -----------------------
    !  For more details see: "A line by line based weighted sum 
    !   of gray gases model for inhomogeneous CO2–H2O mixture  
    !     in oxy-fired combustion" by Bordbar et al.(2014)
    ! doi: https://doi.org/10.1016/j.combustflame.2014.03.013
    !----------------------------------------------------------
    ! <INTERNAL>
    ! nsb    - number of gray gases plus 1 for the clear gas
    ! nGi    - polynomial order
    ! nGj    - polynomial order of the weight coefficient
    ! <OUTPUT>
    ! c_wgh  - asorptivity polynomial coefficients of the WSGG
    ! d_wgh  - emissivity polynomial coefficients of the WSGG
    ! b_wgh  - b asorptivity polynomial coefficient of the WSGG
    !----------------------------------------------------------
    subroutine WSGG_polynomials()
        nsb = 5; nGi = 4; nGj = 5
        allocate(d_wgh(nGi,nsb),c_wgh(nGi,nGj,nsb),b_wgh(nGi,nGj))
        !Mounting the C's array
        c_wgh(1,1,1:5) = (/ 7.412956d-1, -5.244441d-1,  5.822860d-1, -2.096994d-1,  2.420312d-2 /)
        c_wgh(1,2,1:5) = (/-9.412652d-1,  2.799577d-1, -7.672319d-1,  3.204027d-1, -3.910174d-2 /)
        c_wgh(1,3,1:5) = (/ 8.531866d-1,  8.230754d-2,  5.289430d-1, -2.468463d-1,  3.109396d-2 /)
        c_wgh(1,4,1:5) = (/-3.342806d-1,  1.474987d-1, -4.160689d-1,  1.697627d-1, -2.040660d-2 /)
        c_wgh(1,5,1:5) = (/ 4.314362d-2, -6.886217d-2,  1.109773d-1, -4.208608d-2,  4.918817d-3 /)
        c_wgh(2,1,1:5) = (/ 1.552073d-1, -4.862117d-1,  3.668088d-1, -1.055508d-1,  1.058568d-2 /)
        c_wgh(2,2,1:5) = (/ 6.755648d-1,  1.409271d+0, -1.383449d+0,  4.575210d-1, -5.019760d-2 /)
        c_wgh(2,3,1:5) = (/-1.125394d+0, -5.913199d-1,  9.085441d-1, -3.334201d-1,  3.842361d-2 /)
        c_wgh(2,4,1:5) = (/ 6.040543d-1, -5.533854d-2, -1.733014d-1,  7.916083d-2, -9.893357d-3 /)
        c_wgh(2,5,1:5) = (/-1.105453d-1,  4.646634d-2, -1.612982d-3, -3.539835d-3,  6.121277d-4 /)
        c_wgh(3,1,1:5) = (/ 2.550242d-1,  3.805403d-1, -4.249709d-1,  1.429446d-1, -1.574075d-2 /)
        c_wgh(3,2,1:5) = (/-6.065428d-1,  3.494024d-1,  1.853509d-1, -1.013694d-1,  1.302441d-2 /)
        c_wgh(3,3,1:5) = (/ 8.123855d-1, -1.102009d+0,  4.046178d-1, -8.118223d-2,  6.298101d-3 /)
        c_wgh(3,4,1:5) = (/-4.532290d-1,  6.784475d-1, -3.432603d-1,  8.830883d-2, -8.415221d-3 /)
        c_wgh(3,5,1:5) = (/ 8.693093d-2, -1.306996d-1,  7.414464d-2, -2.029294d-2,  2.010969d-3 /)
        c_wgh(4,1,1:5) = (/-3.451994d-2,  2.656726d-1, -1.225365d-1,  3.001508d-2, -2.820525d-3 /)
        c_wgh(4,2,1:5) = (/ 4.112046d-1, -5.728350d-1,  2.924490d-1, -7.980766d-2,  7.996603d-3 /)
        c_wgh(4,3,1:5) = (/-5.055995d-1,  4.579559d-1, -2.616436d-1,  7.648413d-2, -7.908356d-3 /)
        c_wgh(4,4,1:5) = (/ 2.317509d-1, -1.656759d-1,  1.052608d-1, -3.219347d-2,  3.386965d-3 /)
        c_wgh(4,5,1:5) = (/-3.754908d-2,  2.295193d-2, -1.600472d-2,  5.046318d-3, -5.364326d-4 /)
        !Mounting the D's array
        d_wgh(1,1:5) = (/ 3.404288d-2,  6.523048d-2, -4.636852d-2,  1.386835d-2, -1.444993d-3 /)
        d_wgh(2,1:5) = (/ 3.509457d-1,  7.465138d-1, -5.293090d-1,  1.594423d-1, -1.663261d-2 /)
        d_wgh(3,1:5) = (/ 4.570740d+0,  2.168067d+0, -1.498901d+0,  4.917165d-1, -5.429990d-2 /)
        d_wgh(4,1:5) = (/ 1.098169d+2, -5.092359d+1,  2.343236d+1, -5.163892d+0,  4.393889d-1 /)
    end subroutine WSGG_polynomials
    !========================================================================================
    !----------------------------------------------------------
    !       Computes the average absorption coefficient
    !---------------------- Description -----------------------
    ! <INPUT>
    ! T_energy - Temperature (K)
    ! P        - Pressure (Pa)
    ! XCO2     - CO2 mole fraction
    ! XH2O     - H2O mole fraction
    ! K_wgh    - Absorption coefficients of the gray gases (1/m)
    ! a_wgh    - Local array of weights (sum a_wgh = 1)
    ! Lm       - Averaged mean beam length
    ! <INTERNAL>
    ! e_tot    - Total emissivity
    ! <OUTPUT>
    ! K_avg    - Average absorption coefficient
    !----------------------------------------------------------
    function K_avg(P,XCO2,XH2O,K_wgh,a_wgh)
        integer :: i
        double precision :: P,XCO2,XH2O,K_wgh(nsb),a_wgh(nsb),aux,e_tot,K_avg
        if(Lm > 1.0d-4)then
            aux = (P/101325)*(XH2O + XCO2)*Lm
            e_tot = 0.0d0
            do i=1,nsb
                e_tot = e_tot + a_wgh(i)*(1.0d0 - exp(-K_wgh(i)*aux))
            end do
            K_avg = - log(1 - e_tot)/Lm
        elseif(Lm <= 1.0d-4)then
            K_avg = 0.0d0
            do i=1,nsb
                K_avg = K_avg + a_wgh(i)*K_wgh(i)*P
            end do
        end if
    end function K_avg
    !========================================================================================
    !----------------------------------------------------------
    !      Computes the radiative properties of the medium 
    !                  using the Gray Gas Model
    !---------------------- Description -----------------------
    ! For more details see: 
    !        http://www.sandia.gov/TNF/radiation.html
    !----------------------------------------------------------
    ! <INPUT>
    ! Temp     - Temperature (K)
    ! P        - Pressure (Pa)
    ! XCO2     - CO2 mole fraction
    ! XH2O     - H2O mole fraction
    ! <INTERNAL>
    ! C_H2O    - H2O coefficients
    ! C_CO2    - CO2 coefficients
    ! <OUTPUT>
    ! gray_gas - absorption coefficients of the mixture (1/m)
    ! <EXTERNAL ROUTINES>
    ! PMA      - Compute the Planck mean absorption function
    !----------------------------------------------------------
    function gray_gas(Temp,P,XCO2,XH2O)
        double precision :: Temp,P,XCO2,XH2O,C_H2O(6),C_CO2(6),gray_gas
        C_H2O = (/-0.23093d0,-1.12390d0,9.4153d0,-2.99880d0,0.51582d0,-1.8684d-5/) !H2O
        C_CO2 = (/ 18.7410d0,-121.310d0,273.50d0,-194.050d0,56.3100d0,-5.8169d0/)  !CO2
        
        gray_gas = (XH2O*PMA(Temp,C_H2O) + XCO2*PMA(Temp,C_H2O))*(P/101325.0d0)

    end function gray_gas
    !========================================================================================
    !----------------------------------------------------------
    !       Compute the Planck mean absorption function
    !---------------------- Description -----------------------
    ! <INPUT>
    ! T     - Temperature (K)
    ! coef  - Gas coefficients
    ! <INTERNAL>
    ! Tref  - Reference Temperature (K)
    ! <OUTPUT>
    ! PMA   - Planck mean absorption coefficient (m-1atm-1)
    !----------------------------------------------------------
    function PMA(T,Coef)
        integer :: i
        double precision :: PMA,T,Tref,coef(6)
        Tref = (1000.0d0/T)
        PMA = coef(1)
        do i=2,6
            PMA = PMA + coef(i)*Tref**(i-1)
        end do
    end function PMA
    !========================================================================================
end module absorptiondata
