!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c f c _ f l u x _ m o d
!
!  Purpose : air-sea CFCs flux
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit
!
! This file is part of CLIMBER-X.
!
! CLIMBER-X is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! CLIMBER-X is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with CLIMBER-X.  If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module cfc_flux_mod

  use precision, only : wp

  use ocn_grid, only : maxi, maxj, mask_ocn

  implicit none

  private
  public :: cfc_flux

contains

! Compute surface CFC 11 and CFC 12 fluxes 
  subroutine cfc_flux(f_sic,wind,slp,t,s,pcfc11_atm,pcfc12_atm,cfc11,cfc12, &
                      cfc11_flx, cfc12_flx)

    implicit none

    real(wp), intent(in) :: f_sic(:,:)  ! sea ice fraction
    real(wp), intent(in) :: wind(:,:)   ! wind speed [m/s]
    real(wp), intent(in) :: slp(:,:)    ! 
    real(wp), intent(in) :: t(:,:)  ! sea surface temperature [degreeC]
    real(wp), intent(in) :: s(:,:)  ! sea surface salinity [psu]
    real(wp), intent(in) :: pcfc11_atm, pcfc12_atm ! atmospheric concentrations of CFC11 and CFC12 [pptv]
    real(wp), intent(in) :: cfc11(:,:), cfc12(:,:)  ! surface ocean CFC concentrations [mol/m3]

    real(wp), intent(out) :: cfc11_flx(:,:), cfc12_flx(:,:)

    integer :: i, j
    real(wp) :: kw11, kw12
    real(wp) :: alpha11, alpha12
    real(wp) :: Csat11, Csat12
    real(wp) :: Sc11, Sc12
    real(wp), parameter :: Xconvxa = 6.9722e-07_wp
    real(wp), parameter :: Pa0 = 101300._wp ! Pa

      
    do  i=1,maxi
      do  j=1,maxj
        if (mask_ocn(i,j).eq.1) then

          ! Solubility, mol/m3/pptv
          alpha11 = sol_cfc(t(i,j),s(i,j),11)
          alpha12 = sol_cfc(t(i,j),s(i,j),12)

          ! Atmosphere concentration, mol/m3
          Csat11 = alpha11 * pcfc11_atm * slp(i,j)/Pa0
          Csat12 = alpha12 * pcfc12_atm * slp(i,j)/Pa0

          ! Schmidt number
          Sc11 = sc_cfc(t(i,j),11)
          Sc12 = sc_cfc(t(i,j),12)

          ! Piston velocity, m/s
          kw11 = (1._wp - f_sic(i,j)) * Xconvxa * wind(i,j)**2 * sqrt(660._wp/Sc11)
          kw12 = (1._wp - f_sic(i,j)) * Xconvxa * wind(i,j)**2 * sqrt(660._wp/Sc12)

          !               Kw11=(1-FICE2(i,n))*XKW2(i,n)*(Sc11/660)**(-0.5)
          !               Kw12=(1-FICE2(i,n))*XKW2(i,n)*(Sc12/660)**(-0.5)

          ! Surface CFC flux, mol/m2/s, positive out of the ocean
          cfc11_flx(i,j) = -kw11*(Csat11 - cfc11(i,j))
          cfc12_flx(i,j) = -kw12*(Csat12 - cfc12(i,j))

          !              if (i.eq.1 .and. j.eq.10) then
          !                 print *, 'alpha11 ', alpha11
          !                 print *, 'Sc11 ',Sc11
          !                 print *, 'wind ', wind(i,j)
          !                 print *, 'FICE ', f_sic(i,j)
          !                 print *, 'Kw11 ', Kw11
          !                 print *, '1/Kw11 ', 1/Kw11
          !                 print *, 'pCFC11_atm', pCFC11_atm
          !                 print *, 'P sea level ', slp(i,j)
          !                 print *, 'Csat11 ',Csat11
          !                 print *, 'CFC11 ', CFC11(i,j)
          !                 print *, 'CFC11_flx ', CFC11_flx(i,j)
          !              endif

        endif
      enddo
    enddo         

    return

  end subroutine cfc_flux

!********************************************************************

!--------------------------------------------------

      real(wp) function sc_cfc(t,kn)
!---------------------------------------------------
!     CFC 11 and 12 schmidt number as a fonction of temperature. 
!
!     ref: Zheng et al (1998), JGR, vol 103,No C1 
!
!     t: temperature (degree Celcius)
!     kn: = 11 for CFC-11,  12 for CFC-12
!
!     J-C Dutay - LSCE
!---------------------------------------------------

      implicit none

      real(wp), intent(in) :: t
      integer, intent(in) :: kn

      real(wp) :: a1 ( 11: 12), a2 ( 11: 12), a3 ( 11: 12), a4 ( 11: 12)

!   coefficients with t in degree Celcius
!   ------------------------------------
      a1(11) = 3501.8
      a2(11) = -210.31
      a3(11) =    6.1851
      a4(11) =   -0.07513

      a1(12) = 3845.4
      a2(12) = -228.95
      a3(12) =    6.1908
      a4(12) =   -0.067430


      sc_cfc = a1(kn) + a2(kn) * t + a3(kn) *t*t + a4(kn) *t*t*t
  
      return
      end


!********************************************************************


      real(wp) function sol_cfc(pt,ps,kn)
!-------------------------------------------------------------------
!
!     CFC 11 and 12 Solubilities in seawater
!     ref: Warner & Weiss (1985) , Deep Sea Research, vol32
!
!     pt:       temperature (degre Celcius)
!     ps:       salinity    (o/oo)
!     kn:       11 = CFC-11, 12 = CFC-12
!     sol_cfc:  in mol/m3/pptv
!               1 pptv = 1 part per trillion = 10^-12 atm = 1 picoatm
!
!
!     J-C Dutay - LSCE
!-------------------------------------------------------------------

      real(wp), intent(in) :: pt, ps
      integer, intent(in) :: kn

      real(wp) :: ta, d
      real(wp) :: a1 ( 11: 12), a2 ( 11: 12), a3 ( 11: 12), a4 ( 11: 12)
      real(wp) :: b1 ( 11: 12), b2 ( 11: 12), b3 ( 11: 12)


! coefficient for solubility in  mol/l/atm
! ----------------------------------------

!     for CFC 11
!     ----------
      a1 ( 11) = -229.9261
      a2 ( 11) =  319.6552
      a3 ( 11) =  119.4471
      a4 ( 11) =   -1.39165
      b1 ( 11) =   -0.142382
      b2 ( 11) =    0.091459
      b3 ( 11) =   -0.0157274
    
!     for CFC/12
!     ----------
      a1 ( 12) = -218.0971
      a2 ( 12) =  298.9702
      a3 ( 12) =  113.8049
      a4 ( 12) =   -1.39165
      b1 ( 12) =   -0.143566
      b2 ( 12) =    0.091015
      b3 ( 12) =   -0.0153924


      ta   = ( pt + 273.16_wp)* 0.01_wp
      d    = ( b3 ( kn)* ta + b2 ( kn))* ta + b1 ( kn)


      sol_cfc = &
                 exp ( a1 ( kn) &
         +       a2 ( kn)/ ta &
         +       a3 ( kn)* log ( ta ) &
         +       a4 ( kn)* ta * ta  + ps* d )

!     conversion from mol/(l * atm) to mol/(m^3 * atm) 
!     ------------------------------------------------
      sol_cfc = 1000._wp * sol_cfc

!     conversion from mol/(m^3 * atm) to mol/(m3 * pptv) 
!     --------------------------------------------------
      sol_cfc = 1.0e-12_wp * sol_cfc

      end

end module cfc_flux_mod
