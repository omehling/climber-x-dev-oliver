!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  i n i t _ t e m p _ w a t e r _ a g e _ m
!
!> @file
!!
!! Initial temperature, water content and age.
!!
!! @section Copyright
!!
!! Copyright 2009-2017 Ralf Greve, Thorben Dunse
!!
!! @section License
!!
!! This file is part of SICOPOLIS.
!!
!! SICOPOLIS is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! SICOPOLIS is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with SICOPOLIS.  If not, see <http://www.gnu.org/licenses/>.
!<
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------------------------------------------------------------
!> Initial temperature, water content and age.
!<------------------------------------------------------------------------------
module init_temp_water_age_m

  use sico_types_m
  use sico_state
  use sico_grid_mod
  use sico_params, only : sico_par_class, epsi, DELTA_TM_SW 
  use timer, only : sec_year

  implicit none

  private
  public :: init_temp_water_age_1_1, init_temp_water_age_1_2
  public :: init_temp_water_age_1_3, init_temp_water_age_1_4
  public :: init_temp_water_age_2

contains

!-------------------------------------------------------------------------------
!> Initial temperature, water content and age
!! (case ANF_DAT==1, i_TEMP_INIT==1:
!! present-day initial topography, isothermal conditions).
!<------------------------------------------------------------------------------
  subroutine init_temp_water_age_1_1(st,grd,par)

  implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_par_class), intent(in) :: par
  integer :: i, j, kc

!-------- Initial ice temperature --------

  do i=0, grd%IMAX
  do j=0, grd%JMAX

     do kc=0, grd%KCMAX
        st%temp_c(kc,j,i) =  par%temp_init
     end do

  end do
  end do

!-------- Initial lithosphere temperature, water content and age --------

  call init_temp_r(st,grd)
  call init_water(st)
  call init_age(st)

  end subroutine init_temp_water_age_1_1

!-------------------------------------------------------------------------------
!> Initial temperature, water content and age
!! (case ANF_DAT==1, i_TEMP_INIT==2:
!! present-day initial topography,
!! ice temperature equal to local surface temperature).
!<------------------------------------------------------------------------------
  subroutine init_temp_water_age_1_2(st,grd)

  implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
  integer :: i, j, kc

!-------- Initial ice temperature --------

  do i=0, grd%IMAX
  do j=0, grd%JMAX

     do kc=0, grd%KCMAX
        st%temp_c(kc,j,i) = st%temp_s(j,i)
     end do

  end do
  end do

!-------- Initial lithosphere temperature, water content and age --------

  call init_temp_r(st,grd)
  call init_water(st)
  call init_age(st)

  end subroutine init_temp_water_age_1_2

!-------------------------------------------------------------------------------
!> Initial temperature, water content and age
!! (case ANF_DAT==1, i_TEMP_INIT==3:
!! present-day initial topography,
!! ice temperature linearly increasing with depth).
!<------------------------------------------------------------------------------
  subroutine init_temp_water_age_1_3(st,grd)

  use ice_material_properties_m, only : kappa_val

  implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
  integer :: i, j, kc
  real(wp)     :: kappa_const_val
  real(wp)     :: temp_ice_base

!-------- Initial ice temperature --------

  kappa_const_val = kappa_val(-10.0_wp)

  do i=0, grd%IMAX
  do j=0, grd%JMAX

     if (st%maske(j,i)<=2) then

        do kc=0, grd%KCMAX

           st%temp_c(kc,j,i) = st%temp_s(j,i) &
                            + (st%q_geo(j,i)/kappa_const_val) &
                              *st%H_c(j,i)*(1.0_wp-grd%eaz_c_quotient(kc))
                            ! linear temperature distribution according to the
                            ! geothermal heat flux
        end do

        if (st%temp_c(0,j,i) >= -BETA*st%H_c(j,i)) then

           temp_ice_base = -BETA*st%H_c(j,i)

           do kc=0, grd%KCMAX
              st%temp_c(kc,j,i) = st%temp_s(j,i) &
                               + (temp_ice_base-st%temp_s(j,i)) &
                                 *(1.0_wp-grd%eaz_c_quotient(kc))
           end do

        end if

     else   ! maske(j,i)==3, floating ice

        temp_ice_base = -BETA*st%H_c(j,i) - DELTA_TM_SW

        do kc=0, grd%KCMAX
           st%temp_c(kc,j,i) = st%temp_s(j,i) &
                            + (temp_ice_base-st%temp_s(j,i)) &
                              *(1.0_wp-grd%eaz_c_quotient(kc))
        end do

     end if

  end do
  end do

!-------- Initial lithosphere temperature, water content and age --------

  call init_temp_r(st,grd)
  call init_water(st)
  call init_age(st)

  end subroutine init_temp_water_age_1_3

!-------------------------------------------------------------------------------
!> Initial temperature, water content and age
!! (case ANF_DAT==1, i_TEMP_INIT==4:
!! present-day initial topography, ice temperature from Robin (1955) solution).
!<------------------------------------------------------------------------------
  subroutine init_temp_water_age_1_4(st,grd)

  use ice_material_properties_m, only : kappa_val, c_val

  implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
  integer :: i, j, kc
  real(wp)     :: kappa_const_val, c_const_val
  real(wp)     :: as_val, H_val, qgeo_val, K, z_above_base
  real(wp)     :: erf_val_1, erf_val_2
  real(wp)     :: temp_ice_base, temp_scale_factor

!-------- Initial ice temperature --------

  kappa_const_val = kappa_val(-10.0_wp)
      c_const_val =     c_val(-10.0_wp)

  do i=0, grd%IMAX
  do j=0, grd%JMAX

     if (st%maske(j,i)<=2) then
        as_val = max(st%as_perp(j,i), epsi)
     else   ! maske(j,i)==3, floating ice
        as_val = epsi   ! this will produce an almost linear temperature profile
     end if

     H_val    = max(st%H_c(j,i)  , eps)
     qgeo_val = max(st%q_geo(j,i), eps)

     K = sqrt( (kappa_const_val/(RHO*c_const_val)) * (H_val/as_val) )

     erf_val_1 = erf(st%H_c(j,i)/(sqrt(2.0_wp)*K))

     do kc=0, grd%KCMAX
        z_above_base   = st%H_c(j,i)*grd%eaz_c_quotient(kc)
        erf_val_2      = erf(z_above_base/(sqrt(2.0_wp)*K))
        st%temp_c(kc,j,i) = st%temp_s(j,i) &
                          + (qgeo_val/kappa_const_val) &
                            * sqrt(0.5_wp*pi)*K*(erf_val_1-erf_val_2)
     end do

     if ( (st%maske(j,i) <= 2).and.(st%temp_c(0,j,i) >= -BETA*st%H_c(j,i)) ) then
        temp_ice_base     = -BETA*st%H_c(j,i)
        temp_scale_factor = (temp_ice_base-st%temp_s(j,i)) &
                                  /(st%temp_c(0,j,i)-st%temp_s(j,i))
     else if (st%maske(j,i) == 3) then
        temp_ice_base     = -BETA*st%H_c(j,i)-DELTA_TM_SW
        temp_scale_factor = (temp_ice_base-st%temp_s(j,i)) &
                                /(st%temp_c(0,j,i)-st%temp_s(j,i))
     else
        temp_scale_factor = 1.0_wp
     end if

     do kc=0, grd%KCMAX
        st%temp_c(kc,j,i) = st%temp_s(j,i) &
                         + temp_scale_factor*(st%temp_c(kc,j,i)-st%temp_s(j,i))
     end do

  end do
  end do

!-------- Initial lithosphere temperature, water content and age --------

  call init_temp_r(st,grd)
  call init_water(st)
  call init_age(st)

  end subroutine init_temp_water_age_1_4

!-------------------------------------------------------------------------------
!> Initial temperature, water content and age
!! (case ANF_DAT==2: ice-free conditions with relaxed bedrock).
!<------------------------------------------------------------------------------
  subroutine init_temp_water_age_2(st,grd)

  implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
  integer :: i, j, kc

!-------- Initial ice temperature --------

  do i=0, grd%IMAX
  do j=0, grd%JMAX

     do kc=0, grd%KCMAX
        st%temp_c(kc,j,i) = st%temp_s(j,i)
     end do

  end do
  end do

!-------- Initial lithosphere temperature, water content and age --------

  call init_temp_r(st,grd)
  call init_water(st)
  call init_age(st)

  end subroutine init_temp_water_age_2

!-------------------------------------------------------------------------------
!> Initial lithosphere temperature.
!<------------------------------------------------------------------------------
  subroutine init_temp_r(st,grd)

  implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
  integer :: i, j, kr

  do i=0, grd%IMAX
  do j=0, grd%JMAX

     do kr=0, grd%KRMAX
        st%temp_r(kr,j,i) = st%temp_c(0,j,i) &
                         + (st%q_geo(j,i)/KAPPA_R) &
                           *grd%H_R*(1.0_wp-grd%zeta_r(kr))
              ! linear temperature distribution according to the
              ! geothermal heat flux
     end do

  end do
  end do

  end subroutine init_temp_r

!-------------------------------------------------------------------------------
!> Initial water content.
!<------------------------------------------------------------------------------
  subroutine init_water(st)

  implicit none
type(sico_state_class), intent(inout) :: st

  st%omega_c = 0.0_wp   ! only required for the enthalpy method
  st%omega_t = 0.0_wp

  end subroutine init_water

!-------------------------------------------------------------------------------
!> Initial age.
!<------------------------------------------------------------------------------
  subroutine init_age(st)

  implicit none
type(sico_state_class), intent(inout) :: st

  st%age_c = 15000.0_wp*sec_year
  st%age_t = 15000.0_wp*sec_year


  end subroutine init_age

!-------------------------------------------------------------------------------

end module init_temp_water_age_m
!
