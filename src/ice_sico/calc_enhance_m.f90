!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ e n h a n c e _ m
!
!> @file
!!
!! Computation of the flow enhancement factor.
!!
!! @section Copyright
!!
!! Copyright 2009-2017 Ralf Greve
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
!> Computation of the flow enhancement factor.
!<------------------------------------------------------------------------------
module calc_enhance_m

  use sico_types_m
  use sico_state
  use sico_grid_mod
  use sico_timer
  use sico_params, only : sico_par_class
  use timer, only : sec_year

  implicit none

  private
  public :: calc_enhance

contains

!-------------------------------------------------------------------------------
!> Computation of the flow enhancement factor.
!<------------------------------------------------------------------------------
  subroutine calc_enhance(st,grd,tmr,par)

  implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class), intent(in)     :: grd
  type(sico_timer_class), intent(in)    :: tmr
  type(sico_par_class), intent(in)    :: par

  integer :: i, j, kc, kt
  real(wp)     :: age_trans
  real(wp)     :: date_trans1, date_trans2, date_trans3


if (par%enhmod==1) then
!! constant for grounded ice, constant for floating ice.

  st%enh_t = par%ENH_FACT
  st%enh_c = par%ENH_FACT
  if (par%margin==3) then ! /* floating ice */
    call calc_enhance_floating_const(st,grd,par)
  endif

else if (par%enhmod==2) then
!! two different values depending on age for grounded ice,
!! constant for floating ice.

  age_trans = par%AGE_TRANS_0*sec_year

  do i=0, grd%IMAX
  do j=0, grd%JMAX

     do kt=0, grd%KTMAX
        if (st%age_t(kt,j,i) < age_trans) then
           st%enh_t(kt,j,i) = par%ENH_INTG   ! Holocene ice
        else
           st%enh_t(kt,j,i) = par%ENH_FACT   ! Pleistocene ice
        end if
     end do

     do kc=0, grd%KCMAX
        if (st%age_c(kc,j,i) < age_trans) then
           st%enh_c(kc,j,i) = par%ENH_INTG   ! Holocene ice
        else
           st%enh_c(kc,j,i) = par%ENH_FACT   ! Pleistocene ice
        end if
     end do

  end do
  end do

  if (par%margin==3) then ! /* floating ice */
    call calc_enhance_floating_const(st,grd,par)
  endif

else if (par%enhmod==3) then
!! two different values depending on time of deposition for grounded ice,
!! constant for floating ice.

  date_trans1 = par%DATE_TRANS1_0*sec_year
  date_trans2 = par%DATE_TRANS2_0*sec_year
  date_trans3 = par%DATE_TRANS3_0*sec_year

  do i=0, grd%IMAX
  do j=0, grd%JMAX

     do kt=0, grd%KTMAX
        if ( (tmr%time-st%age_t(kt,j,i)) < date_trans1 ) then
           st%enh_t(kt,j,i) = par%ENH_FACT   ! pre-Eemian ice
        else if ( ((tmr%time-st%age_t(kt,j,i)) >= date_trans1).and. &
                  ((tmr%time-st%age_t(kt,j,i)) <  date_trans2) ) then
           st%enh_t(kt,j,i) = par%ENH_INTG   ! Eemian ice
        else if ( ((tmr%time-st%age_t(kt,j,i)) >= date_trans2).and. &
                  ((tmr%time-st%age_t(kt,j,i)) <  date_trans3) ) then
           st%enh_t(kt,j,i) = par%ENH_FACT   ! Weichselian ice
        else
           st%enh_t(kt,j,i) = par%ENH_INTG   ! Holocene ice
        end if
     end do

     do kc=0, grd%KCMAX
        if ( (tmr%time-st%age_c(kc,j,i)) < date_trans1 ) then
           st%enh_c(kc,j,i) = par%ENH_FACT   ! pre-Eemian ice
        else if ( ((tmr%time-st%age_c(kc,j,i)) >= date_trans1).and. &
                  ((tmr%time-st%age_c(kc,j,i)) <  date_trans2) ) then
           st%enh_c(kc,j,i) = par%ENH_INTG   ! Eemian ice
        else if ( ((tmr%time-st%age_c(kc,j,i)) >= date_trans2).and. &
                  ((tmr%time-st%age_c(kc,j,i)) <  date_trans3) ) then
           st%enh_c(kc,j,i) = par%ENH_FACT   ! Weichselian ice
        else
           st%enh_c(kc,j,i) = par%ENH_INTG   ! Holocene ice
        end if
     end do

  end do
  end do

  if (par%margin==3) then   ! /* floating ice */
    call calc_enhance_floating_const(st,grd,par)
  endif

else if (par%enhmod==4) then
!! minimal anisotropic enhancement factor for grounded ice,
!! constant for floating ice.

  call calc_enhance_aniso(st,grd,par)

  if (par%margin==3) then ! /* floating ice */
    call calc_enhance_floating_const(st,grd,par)
  endif

else if (par%enhmod==5) then
!! minimal anisotropic enhancement factor for grounded and floating ice.

  call calc_enhance_aniso(st,grd,par)

endif

  end subroutine calc_enhance


!-------------------------------------------------------------------------------
!> Minimal anisotropic flow enhancement factor.
!<------------------------------------------------------------------------------
  subroutine calc_enhance_aniso(st,grd,par)

  implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class), intent(in)     :: grd
  type(sico_par_class), intent(in)    :: par

  integer :: i, j, kc, kt
  real(wp)     :: enh_shear, enh_compr
  real(wp)     :: enh_shear_compr_diff

  enh_shear = par%ENH_SHEAR

  enh_compr = par%ENH_COMPR

  enh_shear_compr_diff = enh_shear-enh_compr

  do i=0, grd%IMAX
  do j=0, grd%JMAX

     do kt=0, grd%KTMAX
        st%enh_t(kt,j,i) = enh_compr &
                        + enh_shear_compr_diff*st%lambda_shear_t(kt,j,i)**2
     end do

     do kc=0, grd%KCMAX
        st%enh_c(kc,j,i) = enh_compr &
                        + enh_shear_compr_diff*st%lambda_shear_c(kc,j,i)**2
     end do

  end do
  end do

  end subroutine calc_enhance_aniso

!-------------------------------------------------------------------------------
!> Constant, prescribed flow enhancement factor for floating ice.
!<------------------------------------------------------------------------------
  subroutine calc_enhance_floating_const(st,grd,par)

  implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class), intent(in)     :: grd
  type(sico_par_class), intent(in)    :: par

  integer :: i, j, kc, kt
  real(wp)     :: enh_shelf

  enh_shelf = par%ENH_SHELF

  do i=0, grd%IMAX
  do j=0, grd%JMAX

     if ( st%maske(j,i)==3 ) then   ! floating ice

        do kt=0, grd%KTMAX
           st%enh_t(kt,j,i) = enh_shelf
        end do

        do kc=0, grd%KCMAX
           st%enh_c(kc,j,i) = enh_shelf
        end do

     end if

  end do
  end do

  end subroutine calc_enhance_floating_const


end module calc_enhance_m
!
