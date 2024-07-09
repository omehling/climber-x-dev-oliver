!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ t e m p _ e n t h _ m
!
!> @file
!!
!! Computation of temperature, water content and age with the enthalpy method.
!!
!! @section Copyright
!!
!! Copyright 2013-2017 Ralf Greve, Heinz Blatter
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
!> Computation of temperature, water content and age with the enthalpy method.
!<------------------------------------------------------------------------------
module calc_temp_enth_m

  use sico_types_m
  use sico_state
  use sico_grid_mod
  use sico_timer
  use sico_params, only : sico_par_class, delta_tm_sw, omega_max
  use timer, only : sec_year

  implicit none

  private
  public :: calc_temp_enth

contains

!-------------------------------------------------------------------------------
!> Main subroutine of calc_temp_enth_m:
!! Computation of temperature, water content and age with the enthalpy method.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth(st,grd,tmr,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class),  intent(in)    :: grd
type(sico_timer_class), intent(in)    :: tmr
type(sico_par_class), intent(in)      :: par

integer :: i, j, kc, kt, kr, ii, jj


!-------- Computation loop --------

!$omp parallel do &
!$omp private (i,j,kc) 
do i=1, grd%IMAX-1   ! skipping domain margins
do j=1, grd%JMAX-1   ! skipping domain margins

   if (st%maske(j,i)==0) then   ! glaciated land

!  ------ Old vertical column cold

      if (st%n_cts(j,i) == -1) then

         st%n_cts_neu(j,i)  = -1
         st%kc_cts_neu(j,i) = 0
         st%zm_neu(j,i)     = st%zb(j,i)
         st%H_c_neu(j,i)    = st%H_c(j,i)
         st%H_t_neu(j,i)    = 0.0_wp

         call calc_temp_enth_1(st,grd,tmr,par,i, j)

!    ---- Check whether base has become temperate

         if (st%temp_c_neu(0,j,i) > st%temp_c_m(0,j,i)-eps) then

            st%n_cts_neu(j,i)  = 0
            st%kc_cts_neu(j,i) = 0

            call calc_temp_enth_2(st,grd,tmr,par,i, j)

         end if

!  ------ Old vertical column with temperate base

      else if (st%n_cts(j,i) == 0) then

         st%n_cts_neu(j,i)  = 0
         st%kc_cts_neu(j,i) = st%kc_cts(j,i)
         st%zm_neu(j,i)     = st%zb(j,i)
         st%H_c_neu(j,i)    = st%H_c(j,i)
         st%H_t_neu(j,i)    = st%H_t(j,i)

         call calc_temp_enth_2(st,grd,tmr,par,i, j)

!    ---- Check whether temperate base becomes cold

         if ( (st%temp_c_neu(1,j,i)-st%temp_c_neu(0,j,i)) < (grd%am1*st%H_c(j,i)) ) then

            st%n_cts_neu(j,i)  = -1
            st%kc_cts_neu(j,i) = 0

            call calc_temp_enth_1(st,grd,tmr,par,i, j)

            if (st%temp_c_neu(0,j,i) > st%temp_c_m(0,j,i)-eps) then

               st%n_cts_neu(j,i)  = 0
               st%kc_cts_neu(j,i) = 0

               call calc_temp_enth_2(st,grd,tmr,par,i, j)

            end if

         end if

      end if

   else if (st%maske(j,i)==3 .and. par%margin==3) then   ! floating ice

      st%n_cts_neu(j,i)  = -1
      st%kc_cts_neu(j,i) = 0
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = 0.0_wp

      call calc_temp_enth_ssa(st,grd,tmr,par,i, j)

!  ------ Reset temperatures above melting to the melting point
!         and water contents above zero to zero
!         (should not occur, but just in case)

      do kc=0, grd%KCMAX
         if (st%temp_c_neu(kc,j,i) > st%temp_c_m(kc,j,i)) &
                    st%temp_c_neu(kc,j,i) = st%temp_c_m(kc,j,i)
         if (st%omega_c_neu(kc,j,i) > 0.0_wp) &
                    st%omega_c_neu(kc,j,i) = 0.0_wp
      end do

   else   ! st%maske(j,i) == 1,2 (ice-free land or sea point)

      st%n_cts_neu(j,i)  = -1
      st%kc_cts_neu(j,i) = 0
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = 0.0_wp

      call calc_temp_enth_r(st,grd,par,tmr,i, j)

   end if

end do
end do 
!$omp end parallel do

!-------- Extrapolate values on margins --------

!  ------ Lower left corner

i=0
j=0

if ( (st%maske(j,i) == 0).or.(st%maske(j,i) == 3) ) then
                              ! glaciated land or floating ice
   ii=i+1
   jj=j+1

   do kc=0, grd%KCMAX
      st%enth_c_neu(kc,j,i)  = st%enth_c_neu(kc,jj,ii)
      st%temp_c_neu(kc,j,i)  = st%temp_c_neu(kc,jj,ii)
      st%omega_c_neu(kc,j,i) = st%omega_c_neu(kc,jj,ii)
      st%age_c_neu(kc,j,i)   = st%age_c_neu(kc,jj,ii)
   end do

   do kt=0, grd%KTMAX
            ! redundant, lower (kt) ice layer
      st%enth_t_neu(kt,j,i)  = st%enth_t_neu(kt,jj,ii)
      st%omega_t_neu(kt,j,i) = st%omega_t_neu(kt,jj,ii)
      st%age_t_neu(kt,j,i)   = st%age_t_neu(kt,jj,ii)
   end do

   do kr=0, grd%KRMAX
      st%temp_r_neu(kr,j,i)  = st%temp_r_neu(kr,jj,ii)
   end do

   st%n_cts_neu(j,i)  = st%n_cts_neu(jj,ii)
   st%kc_cts_neu(j,i) = st%kc_cts_neu(jj,ii)
   st%zm_neu(j,i)     = st%zb(j,i)
   st%H_c_neu(j,i)    = st%H_c(j,i)
   st%H_t_neu(j,i)    = st%H_t(j,i)

else   ! st%maske(j,i) == 1,2 (ice-free land or sea point)

   st%n_cts_neu(j,i)  = -1
   st%kc_cts_neu(j,i) = 0
   st%zm_neu(j,i)     = st%zb(j,i)
   st%H_c_neu(j,i)    = st%H_c(j,i)
   st%H_t_neu(j,i)    = st%H_t(j,i)

   call calc_temp_enth_r(st,grd,par,tmr,i, j)

end if

!  ------ Lower right corner

i=grd%IMAX
j=0

if ( (st%maske(j,i) == 0).or.(st%maske(j,i) == 3) ) then
                              ! glaciated land or floating ice
   ii=i-1
   jj=j+1

   do kc=0, grd%KCMAX
      st%enth_c_neu(kc,j,i)  = st%enth_c_neu(kc,jj,ii)
      st%temp_c_neu(kc,j,i)  = st%temp_c_neu(kc,jj,ii)
      st%omega_c_neu(kc,j,i) = st%omega_c_neu(kc,jj,ii)
      st%age_c_neu(kc,j,i)   = st%age_c_neu(kc,jj,ii)
   end do

   do kt=0, grd%KTMAX
            ! redundant, lower (kt) ice layer
      st%enth_t_neu(kt,j,i)  = st%enth_t_neu(kt,jj,ii)
      st%omega_t_neu(kt,j,i) = st%omega_t_neu(kt,jj,ii)
      st%age_t_neu(kt,j,i)   = st%age_t_neu(kt,jj,ii)
   end do

   do kr=0, grd%KRMAX
      st%temp_r_neu(kr,j,i)  = st%temp_r_neu(kr,jj,ii)
   end do

   st%n_cts_neu(j,i)  = st%n_cts_neu(jj,ii)
   st%kc_cts_neu(j,i) = st%kc_cts_neu(jj,ii)
   st%zm_neu(j,i)     = st%zb(j,i)
   st%H_c_neu(j,i)    = st%H_c(j,i)
   st%H_t_neu(j,i)    = st%H_t(j,i)

else   ! st%maske(j,i) == 1,2 (ice-free land or sea point)

   st%n_cts_neu(j,i)  = -1
   st%kc_cts_neu(j,i) = 0
   st%zm_neu(j,i)     = st%zb(j,i)
   st%H_c_neu(j,i)    = st%H_c(j,i)
   st%H_t_neu(j,i)    = st%H_t(j,i)

   call calc_temp_enth_r(st,grd,par,tmr,i, j)

end if

!  ------ Upper left corner

i=0
j=grd%JMAX

if ( (st%maske(j,i) == 0).or.(st%maske(j,i) == 3) ) then
                              ! glaciated land or floating ice
   ii=i+1
   jj=j-1

   do kc=0, grd%KCMAX
      st%enth_c_neu(kc,j,i)  = st%enth_c_neu(kc,jj,ii)
      st%temp_c_neu(kc,j,i)  = st%temp_c_neu(kc,jj,ii)
      st%omega_c_neu(kc,j,i) = st%omega_c_neu(kc,jj,ii)
      st%age_c_neu(kc,j,i)   = st%age_c_neu(kc,jj,ii)
   end do

   do kt=0, grd%KTMAX
            ! redundant, lower (kt) ice layer
      st%enth_t_neu(kt,j,i)  = st%enth_t_neu(kt,jj,ii)
      st%omega_t_neu(kt,j,i) = st%omega_t_neu(kt,jj,ii)
      st%age_t_neu(kt,j,i)   = st%age_t_neu(kt,jj,ii)
   end do

   do kr=0, grd%KRMAX
      st%temp_r_neu(kr,j,i)  = st%temp_r_neu(kr,jj,ii)
   end do

   st%n_cts_neu(j,i)  = st%n_cts_neu(jj,ii)
   st%kc_cts_neu(j,i) = st%kc_cts_neu(jj,ii)
   st%zm_neu(j,i)     = st%zb(j,i)
   st%H_c_neu(j,i)    = st%H_c(j,i)
   st%H_t_neu(j,i)    = st%H_t(j,i)

else   ! st%maske(j,i) == 1,2 (ice-free land or sea point)

   st%n_cts_neu(j,i)  = -1
   st%kc_cts_neu(j,i) = 0
   st%zm_neu(j,i)     = st%zb(j,i)
   st%H_c_neu(j,i)    = st%H_c(j,i)
   st%H_t_neu(j,i)    = st%H_t(j,i)

   call calc_temp_enth_r(st,grd,par,tmr,i, j)

end if

!  ------ Upper right corner

i=grd%IMAX
j=grd%JMAX

if ( (st%maske(j,i) == 0).or.(st%maske(j,i) == 3) ) then
                              ! glaciated land or floating ice
   ii=i-1
   jj=j-1

   do kc=0, grd%KCMAX
      st%enth_c_neu(kc,j,i)  = st%enth_c_neu(kc,jj,ii)
      st%temp_c_neu(kc,j,i)  = st%temp_c_neu(kc,jj,ii)
      st%omega_c_neu(kc,j,i) = st%omega_c_neu(kc,jj,ii)
      st%age_c_neu(kc,j,i)   = st%age_c_neu(kc,jj,ii)
   end do

   do kt=0, grd%KTMAX
            ! redundant, lower (kt) ice layer
      st%enth_t_neu(kt,j,i)  = st%enth_t_neu(kt,jj,ii)
      st%omega_t_neu(kt,j,i) = st%omega_t_neu(kt,jj,ii)
      st%age_t_neu(kt,j,i)   = st%age_t_neu(kt,jj,ii)
   end do

   do kr=0, grd%KRMAX
      st%temp_r_neu(kr,j,i)  = st%temp_r_neu(kr,jj,ii)
   end do

   st%n_cts_neu(j,i)  = st%n_cts_neu(jj,ii)
   st%kc_cts_neu(j,i) = st%kc_cts_neu(jj,ii)
   st%zm_neu(j,i)     = st%zb(j,i)
   st%H_c_neu(j,i)    = st%H_c(j,i)
   st%H_t_neu(j,i)    = st%H_t(j,i)

else   ! st%maske(j,i) == 1,2 (ice-free land or sea point)

   st%n_cts_neu(j,i)  = -1
   st%kc_cts_neu(j,i) = 0
   st%zm_neu(j,i)     = st%zb(j,i)
   st%H_c_neu(j,i)    = st%H_c(j,i)
   st%H_t_neu(j,i)    = st%H_t(j,i)

   call calc_temp_enth_r(st,grd,par,tmr,i, j)

end if

!  ------ Lower and upper margins

do i=1, grd%IMAX-1

!    ---- Lower margin

   j=0

   if ( (st%maske(j,i) == 0).or.(st%maske(j,i) == 3) ) then
                                 ! glaciated land or floating ice
      ii=i
      jj=j+1

      do kc=0, grd%KCMAX
         st%enth_c_neu(kc,j,i)  = st%enth_c_neu(kc,jj,ii)
         st%temp_c_neu(kc,j,i)  = st%temp_c_neu(kc,jj,ii)
         st%omega_c_neu(kc,j,i) = st%omega_c_neu(kc,jj,ii)
         st%age_c_neu(kc,j,i)   = st%age_c_neu(kc,jj,ii)
      end do

      do kt=0, grd%KTMAX
               ! redundant, lower (kt) ice layer
         st%enth_t_neu(kt,j,i)  = st%enth_t_neu(kt,jj,ii)
         st%omega_t_neu(kt,j,i) = st%omega_t_neu(kt,jj,ii)
         st%age_t_neu(kt,j,i)   = st%age_t_neu(kt,jj,ii)
      end do

      do kr=0, grd%KRMAX
         st%temp_r_neu(kr,j,i)  = st%temp_r_neu(kr,jj,ii)
      end do

      st%n_cts_neu(j,i)  = st%n_cts_neu(jj,ii)
      st%kc_cts_neu(j,i) = st%kc_cts_neu(jj,ii)
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = st%H_t(j,i)

   else   ! st%maske(j,i) == 1,2 (ice-free land or sea point)

      st%n_cts_neu(j,i)  = -1
      st%kc_cts_neu(j,i) = 0
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = st%H_t(j,i)

      call calc_temp_enth_r(st,grd,par,tmr,i, j)

   end if

!    ---- Upper margin

   j=grd%JMAX

   if ( (st%maske(j,i) == 0).or.(st%maske(j,i) == 3) ) then
                                 ! glaciated land or floating ice
      ii=i
      jj=j-1

      do kc=0, grd%KCMAX
         st%enth_c_neu(kc,j,i)  = st%enth_c_neu(kc,jj,ii)
         st%temp_c_neu(kc,j,i)  = st%temp_c_neu(kc,jj,ii)
         st%omega_c_neu(kc,j,i) = st%omega_c_neu(kc,jj,ii)
         st%age_c_neu(kc,j,i)   = st%age_c_neu(kc,jj,ii)
      end do

      do kt=0, grd%KTMAX
               ! redundant, lower (kt) ice layer
         st%enth_t_neu(kt,j,i)  = st%enth_t_neu(kt,jj,ii)
         st%omega_t_neu(kt,j,i) = st%omega_t_neu(kt,jj,ii)
         st%age_t_neu(kt,j,i)   = st%age_t_neu(kt,jj,ii)
      end do

      do kr=0, grd%KRMAX
         st%temp_r_neu(kr,j,i)  = st%temp_r_neu(kr,jj,ii)
      end do

      st%n_cts_neu(j,i)  = st%n_cts_neu(jj,ii)
      st%kc_cts_neu(j,i) = st%kc_cts_neu(jj,ii)
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = st%H_t(j,i)

   else   ! st%maske(j,i) == 1,2 (ice-free land or sea point)

      st%n_cts_neu(j,i)  = -1
      st%kc_cts_neu(j,i) = 0
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = st%H_t(j,i)

      call calc_temp_enth_r(st,grd,par,tmr,i, j)

   end if

end do

!  ------ Left and right margins

do j=1, grd%JMAX-1

!    ---- Left margin

   i=0

   if ( (st%maske(j,i) == 0).or.(st%maske(j,i) == 3) ) then
                                 ! glaciated land or floating ice
      ii=i+1
      jj=j

      do kc=0, grd%KCMAX
         st%enth_c_neu(kc,j,i)  = st%enth_c_neu(kc,jj,ii)
         st%temp_c_neu(kc,j,i)  = st%temp_c_neu(kc,jj,ii)
         st%omega_c_neu(kc,j,i) = st%omega_c_neu(kc,jj,ii)
         st%age_c_neu(kc,j,i)   = st%age_c_neu(kc,jj,ii)
      end do

      do kt=0, grd%KTMAX
               ! redundant, lower (kt) ice layer
         st%enth_t_neu(kt,j,i)  = st%enth_t_neu(kt,jj,ii)
         st%omega_t_neu(kt,j,i) = st%omega_t_neu(kt,jj,ii)
         st%age_t_neu(kt,j,i)   = st%age_t_neu(kt,jj,ii)
      end do

      do kr=0, grd%KRMAX
         st%temp_r_neu(kr,j,i)  = st%temp_r_neu(kr,jj,ii)
      end do

      st%n_cts_neu(j,i)  = st%n_cts_neu(jj,ii)
      st%kc_cts_neu(j,i) = st%kc_cts_neu(jj,ii)
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = st%H_t(j,i)

   else   ! st%maske(j,i) == 1,2 (ice-free land or sea point)

      st%n_cts_neu(j,i)  = -1
      st%kc_cts_neu(j,i) = 0
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = st%H_t(j,i)

      call calc_temp_enth_r(st,grd,par,tmr,i, j)

   end if

!    ---- Right margin

   i=grd%IMAX

   if ( (st%maske(j,i) == 0).or.(st%maske(j,i) == 3) ) then
                                 ! glaciated land or floating ice
      ii=i-1
      jj=j

      do kc=0, grd%KCMAX
         st%enth_c_neu(kc,j,i)  = st%enth_c_neu(kc,jj,ii)
         st%temp_c_neu(kc,j,i)  = st%temp_c_neu(kc,jj,ii)
         st%omega_c_neu(kc,j,i) = st%omega_c_neu(kc,jj,ii)
         st%age_c_neu(kc,j,i)   = st%age_c_neu(kc,jj,ii)
      end do

      do kt=0, grd%KTMAX
               ! redundant, lower (kt) ice layer
         st%enth_t_neu(kt,j,i)  = st%enth_t_neu(kt,jj,ii)
         st%omega_t_neu(kt,j,i) = st%omega_t_neu(kt,jj,ii)
         st%age_t_neu(kt,j,i)   = st%age_t_neu(kt,jj,ii)
      end do

      do kr=0, grd%KRMAX
         st%temp_r_neu(kr,j,i)  = st%temp_r_neu(kr,jj,ii)
      end do

      st%n_cts_neu(j,i)  = st%n_cts_neu(jj,ii)
      st%kc_cts_neu(j,i) = st%kc_cts_neu(jj,ii)
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = st%H_t(j,i)

   else   ! st%maske(j,i) == 1,2 (ice-free land or sea point)

      st%n_cts_neu(j,i)  = -1
      st%kc_cts_neu(j,i) = 0
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = st%H_t(j,i)

      call calc_temp_enth_r(st,grd,par,tmr,i, j)

   end if

end do

end subroutine calc_temp_enth

!-------------------------------------------------------------------------------
!> Computation of temperature and age for a cold ice column with the
!! enthalpy method.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_1(st,grd,tmr,par, i, j)

use sico_maths_m,      only : tri_sle
use enth_temp_omega_m, only : temp_fct_enth

implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class),  intent(in)    :: grd
  type(sico_timer_class), intent(in)    :: tmr
  type(sico_par_class), intent(in)    :: par
integer, intent(in) :: i, j

integer :: kc, kt, kr
real(wp) :: ct1(0:100), ct2(0:100), ct3(0:100), ct4(0:100), &
            ce5(0:100), ce6(0:100), ce7(0:100), &
            ctr1, ccbe1, ccb2, ccb3, ccb4, clb1
real(wp) :: ct1_sg(0:100), ct2_sg(0:100), ct3_sg(0:100), &
            ct4_sg(0:100), adv_vert_sg(0:100), abs_adv_vert_sg(0:100)
real(wp) :: ci1(0:100), ci2(0:100)
real(wp) :: dtt_dxi, dtt_deta
real(wp) :: lgs_a0(0:200), &
            lgs_a1(0:200), &
            lgs_a2(0:200), &
            lgs_x(0:200), &
            lgs_b(0:200)

!-------- Check for boundary points --------

if ((i == 0).or.(i == grd%IMAX).or.(j == 0).or.(j == grd%JMAX)) &
   stop ' >>> calc_temp_enth_1: Boundary points not allowed.'

!-------- Abbreviations --------

call calc_temp_enth_1_a(st,grd,tmr,par,i, j, &
                        ct1, ct2, ct3, ct4, ce5, ce6, ce7, &
                        ctr1, ccbe1, ccb2, ccb3, ccb4, clb1, &
                        ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                        adv_vert_sg, abs_adv_vert_sg, &
                        ci1, ci2, dtt_dxi, dtt_deta)

!-------- Computation of the bedrock temperature
!         (upper boundary condition: old temperature at the ice base) --------

!  ------ Set-up of the the equations

call calc_temp_enth_1_b(st,grd,par,tmr,ctr1, clb1, i, j, &
                        lgs_a0, lgs_a1, lgs_a2, lgs_b)

!  ------ Solution of the system of linear equations

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KRMAX)

!  ------ Assignment of the result (predictor values)

do kr=0, grd%KRMAX
   st%temp_r_neu(kr,j,i) = lgs_x(kr)
end do

!-------- Computation of the ice enthalpy --------

!  ------ Set-up of the the equations

call calc_temp_enth_1_c(st,grd,par,tmr,ct1, ct2, ct3, ct4, ce5, ce6, ce7, &
                        ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                        ccbe1, ccb2, ccb3, ccb4, &
                        adv_vert_sg, abs_adv_vert_sg, &
                        dtt_dxi, dtt_deta, &
                        i, j, &
                        lgs_a0, lgs_a1, lgs_a2, lgs_b)

!  ------ Solution of the system of linear equations

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KCMAX)

!  ------ Assignment of the result

do kc=0, grd%KCMAX
   st%enth_c_neu(kc,j,i)  = lgs_x(kc)
   st%temp_c_neu(kc,j,i)  = temp_fct_enth(st%enth_c_neu(kc,j,i), st%temp_c_m(kc,j,i))
   st%omega_c_neu(kc,j,i) = 0.0_wp   ! solution is supposed to be for cold ice
end do

!-------- Water drainage from the non-existing temperate ice --------

st%Q_tld(j,i) = 0.0_wp

!-------- Set enthalpy and water content in the redundant,
!         lower (kt) ice layer to the value at the ice base --------

do kt=0, grd%KTMAX
   st%enth_t_neu(kt,j,i)  = st%enth_c_neu(0,j,i)
   st%omega_t_neu(kt,j,i) = st%omega_c_neu(0,j,i)
end do

!-------- Computation of the age of ice --------

!  ------ Set-up of the the equations

call calc_temp_enth_1_d(st,grd,par,tmr,ct1, ct2, ct3, ct4, ci1, ci2, &
                        ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                        adv_vert_sg, abs_adv_vert_sg, &
                        dtt_dxi, dtt_deta, &
                        i, j, &
                        lgs_a0, lgs_a1, lgs_a2, lgs_b)

!  ------ Solution of the system of linear equations

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KCMAX)

!  ------ Assignment of the result,
!         restriction to interval [0, AGE_MAX yr] --------

do kc=0, grd%KCMAX

   st%age_c_neu(kc,j,i) = lgs_x(kc)

   if (st%age_c_neu(kc,j,i) < (par%age_min*sec_year)) &
                           st%age_c_neu(kc,j,i) = 0.0_wp
   if (st%age_c_neu(kc,j,i) > (par%age_max*sec_year)) &
                           st%age_c_neu(kc,j,i) = par%age_max*sec_year

end do

!-------- Age of the ice in the redundant, lower (kt) ice layer --------

do kt=0, grd%KTMAX
   st%age_t_neu(kt,j,i) = st%age_c_neu(0,j,i)
end do

end subroutine calc_temp_enth_1

!-------------------------------------------------------------------------------
!> Computation of temperature and age for a cold ice column with the
!! enthalpy method: Abbreviations.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_1_a(st,grd,tmr,par,i, j, &
                              ct1, ct2, ct3, ct4, ce5, ce6, ce7, &
                              ctr1, ccbe1, ccb2, ccb3, ccb4, clb1, &
                              ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                              adv_vert_sg, abs_adv_vert_sg, &
                              ci1, ci2, dtt_dxi, dtt_deta)

use ice_material_properties_m, only : ratefac_c_t, kappa_val, c_val, &
                                      creep, viscosity

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class),  intent(in)    :: grd
type(sico_timer_class), intent(in)    :: tmr
type(sico_par_class), intent(in)    :: par
integer, intent(in) :: i, j

real(wp),    intent(out) :: ct1(0:100), ct2(0:100), ct3(0:100), &
                            ct4(0:100), ce5(0:100), ce6(0:100), &
                            ce7(0:100), &
                            ctr1, ccbe1, ccb2, ccb3, ccb4, clb1
real(wp),    intent(out) :: ct1_sg(0:100), ct2_sg(0:100), &
                            ct3_sg(0:100), ct4_sg(0:100), &
                            adv_vert_sg(0:100), abs_adv_vert_sg(0:100)
real(wp),    intent(out) :: ci1(0:100), ci2(0:100)
real(wp),    intent(out) :: dtt_dxi, dtt_deta

integer :: kc
real(wp) :: temp_c_help(0:100)

!-------- Initialisation --------

ct1             = 0.0_wp
ct2             = 0.0_wp
ct3             = 0.0_wp
ct4             = 0.0_wp
ce5             = 0.0_wp
ce6             = 0.0_wp
ce7             = 0.0_wp
ctr1            = 0.0_wp
clb1            = 0.0_wp
ct1_sg          = 0.0_wp
ct2_sg          = 0.0_wp
ct3_sg          = 0.0_wp
ct4_sg          = 0.0_wp
adv_vert_sg     = 0.0_wp
abs_adv_vert_sg = 0.0_wp
ci1             = 0.0_wp
ci2             = 0.0_wp
dtt_dxi         = 0.0_wp
dtt_deta        = 0.0_wp

!-------- Actual computation --------

ctr1 = grd%atr1

ccbe1 = grd%acb1 &
        *kappa_val(st%temp_c(0,j,i)) &
        /(c_val(st%temp_c(0,j,i))*st%H_c(j,i))
ccb2  = grd%acb2

if (par%dynamics==2) then

  if (.not.st%flag_shelfy_stream(j,i)) then

    ccb3 = grd%acb3*0.5_wp*(st%vx_t(0,j,i)+st%vx_t(0,j,i-1)) &
      *st%H_c(j,i)*st%dzs_dxi_g(j,i)
    ccb4 = grd%acb4*0.5_wp*(st%vy_t(0,j,i)+st%vy_t(0,j-1,i)) &
      *st%H_c(j,i)*st%dzs_deta_g(j,i)

  else   ! flag_shelfy_stream(j,i) == .true.

    if (par%hyb_mode==1) then
      ccb3 = -st%beta_drag(j,i) &
        * (st%vx_b_g(j,i)**2 + st%vy_b_g(j,i)**2)
    else
      ccb3 = -st%c_drag(j,i) &
        * sqrt(st%vx_b_g(j,i)**2  &
        +st%vy_b_g(j,i)**2) &
        **(1.0_wp+st%p_weert_inv(j,i))
    endif
    ccb4 = 0.0_wp

  end if

else    ! par%dynamics.ne.2

  ccb3 = grd%acb3*0.5_wp*(st%vx_t(0,j,i)+st%vx_t(0,j,i-1)) &
    *st%H_c(j,i)*st%dzs_dxi_g(j,i)
  ccb4 = grd%acb4*0.5_wp*(st%vy_t(0,j,i)+st%vy_t(0,j-1,i)) &
    *st%H_c(j,i)*st%dzs_deta_g(j,i)

endif

clb1 = grd%alb1*st%q_geo(j,i)

if (par%adv_vert==1) then

do kc=1, grd%KCMAX-1
   ct1(kc) = grd%at1(kc)/st%H_c(j,i)*0.5_wp*(st%vz_c(kc,j,i)+st%vz_c(kc-1,j,i))
end do

kc=0
ct1_sg(kc) = 0.5_wp*(grd%at1(kc)+grd%at1(kc+1))/st%H_c(j,i)*st%vz_c(kc,j,i)
             ! only needed for kc=0 ...
kc=grd%KCMAX-1
ct1_sg(kc) = 0.5_wp*(grd%at1(kc)+grd%at1(kc+1))/st%H_c(j,i)*st%vz_c(kc,j,i)
             ! ... and kc=grd%KCMAX-1

else if (par%adv_vert==2 .or. par%adv_vert==3) then

do kc=0, grd%KCMAX-1
   ct1_sg(kc) = 0.5_wp*(grd%at1(kc)+grd%at1(kc+1))/st%H_c(j,i)*st%vz_c(kc,j,i)
end do

endif

do kc=0, grd%KCMAX

   ct2(kc) = ( grd%at2_1(kc)*st%dzm_dtau(j,i) &
           +grd%at2_2(kc)*st%dH_c_dtau(j,i) )/st%H_c(j,i)
   ct3(kc) = ( grd%at3_1(kc)*st%dzm_dxi_g(j,i) &
           +grd%at3_2(kc)*st%dH_c_dxi_g(j,i) )/st%H_c(j,i) &
          *0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1)) *grd%insq_g11_g(j,i)
   ct4(kc) = ( grd%at4_1(kc)*st%dzm_deta_g(j,i) &
            +grd%at4_2(kc)*st%dH_c_deta_g(j,i) )/st%H_c(j,i) &
          *0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i)) *grd%insq_g22_g(j,i)
   ce5(kc) = grd%at5(kc)/st%H_c(j,i)

   if (par%dynamics==2) then

     if (.not.st%flag_shelfy_stream(j,i)) then
       ce7(kc) = grd%at7 &
         *st%enh_c(kc,j,i) &
         *ratefac_c_t(st%temp_c(kc,j,i), st%omega_c(kc,j,i), st%temp_c_m(kc,j,i)) &
         *creep(par%fin_visc,par%flow_law,st%sigma_c(kc,j,i)) &
         *st%sigma_c(kc,j,i)*st%sigma_c(kc,j,i)
     else
       ce7(kc) = 2.0_wp*grd%at7 &
         *viscosity(par%fin_visc,par%flow_law,st%de_c(kc,j,i), &
         st%temp_c(kc,j,i), st%temp_c_m(kc,j,i), 0.0_wp, &
         st%enh_c(kc,j,i), 0) &
         *st%de_c(kc,j,i)**2
     end if

   else   ! par%dynamics.ne.2

     ce7(kc) = grd%at7 &
       *st%enh_c(kc,j,i) &
       *ratefac_c_t(st%temp_c(kc,j,i), st%omega_c(kc,j,i), st%temp_c_m(kc,j,i)) &
       *creep(par%fin_visc,par%flow_law,st%sigma_c(kc,j,i)) &
       *st%sigma_c(kc,j,i)*st%sigma_c(kc,j,i)

   endif

   ci1(kc) = grd%ai1(kc)/st%H_c(j,i)

end do

if (par%adv_vert==1) then

kc=0
ct2_sg(kc) = 0.5_wp*(ct2(kc)+ct2(kc+1))
ct3_sg(kc) = 0.5_wp*(ct3(kc)+ct3(kc+1))
ct4_sg(kc) = 0.5_wp*(ct4(kc)+ct4(kc+1))
adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))   ! only needed for kc=0 ...
kc=grd%KCMAX-1
ct2_sg(kc) = 0.5_wp*(ct2(kc)+ct2(kc+1))
ct3_sg(kc) = 0.5_wp*(ct3(kc)+ct3(kc+1))
ct4_sg(kc) = 0.5_wp*(ct4(kc)+ct4(kc+1))
adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))   ! ... and kc=grd%KCMAX-1

else if (par%adv_vert==2 .or. par%adv_vert==3) then

do kc=0, grd%KCMAX-1
   ct2_sg(kc) = 0.5_wp*(ct2(kc)+ct2(kc+1))
   ct3_sg(kc) = 0.5_wp*(ct3(kc)+ct3(kc+1))
   ct4_sg(kc) = 0.5_wp*(ct4(kc)+ct4(kc+1))
   adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
   abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))
end do

endif

do kc=0, grd%KCMAX-1
   temp_c_help(kc) = 0.5_wp*(st%temp_c(kc,j,i)+st%temp_c(kc+1,j,i))
   ce6(kc) = grd%at6(kc) &
             *kappa_val(st%temp_c_help(kc))/(c_val(st%temp_c_help(kc))*st%H_c(j,i))
   ci2(kc) = grd%ai2(kc)/st%H_c(j,i)
end do

if (par%adv_hor==3) then
dtt_dxi  = 2.0_wp*grd%dtt_2dxi
dtt_deta = 2.0_wp*grd%dtt_2deta
endif

end subroutine calc_temp_enth_1_a

!-------------------------------------------------------------------------------
!> Computation of temperature and age for a cold ice column with the
!! enthalpy method:
!! Set-up of the equations for the bedrock temperature.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_1_b(st,grd,par,tmr,ctr1, clb1, i, j, &
                              lgs_a0, lgs_a1, lgs_a2, lgs_b)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class),  intent(in)    :: grd
type(sico_par_class),   intent(in)    :: par
type(sico_timer_class), intent(in)    :: tmr
integer, intent(in) :: i, j
real(wp),     intent(in) :: ctr1, clb1

real(wp),    intent(out) :: lgs_a0(0:200), &
                            lgs_a1(0:200), &
                            lgs_a2(0:200), &
                            lgs_b(0:200)

integer :: kr

!-------- Initialisation --------

lgs_a0 = 0.0_wp
lgs_a1 = 0.0_wp
lgs_a2 = 0.0_wp
lgs_b  = 0.0_wp

!-------- Actual computation --------

kr=0
lgs_a1(kr) =  1.0_wp
lgs_a2(kr) = -1.0_wp
lgs_b(kr)  = clb1

if (par%q_litho==1) then
!   (coupled heat-conducting bedrock)

do kr=1, grd%KRMAX-1
   lgs_a0(kr) = -ctr1
   lgs_a1(kr) = 1.0_wp + 2.0_wp*ctr1
   lgs_a2(kr) = -ctr1
   lgs_b(kr)  = st%temp_r(kr,j,i)
end do

else if (par%q_litho==0) then
!   (no coupled heat-conducting bedrock)

do kr=1, grd%KRMAX-1
   lgs_a0(kr) =  1.0_wp
   lgs_a1(kr) =  0.0_wp
   lgs_a2(kr) = -1.0_wp
   lgs_b(kr)  =  2.0_wp*clb1
end do

endif

kr=grd%KRMAX
lgs_a0(kr) = 0.0_wp
lgs_a1(kr) = 1.0_wp
lgs_b(kr)  = st%temp_c(0,j,i)

end subroutine calc_temp_enth_1_b

!-------------------------------------------------------------------------------
!> Computation of temperature and age for a cold ice column with the
!! enthalpy method:
!! Set-up of the equations for the ice enthalpy.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_1_c(st,grd,par,tmr,ct1, ct2, ct3, ct4, ce5, ce6, ce7, &
                              ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                              ccbe1, ccb2, ccb3, ccb4, &
                              adv_vert_sg, abs_adv_vert_sg, &
                              dtt_dxi, dtt_deta, &
                              i, j, &
                              lgs_a0, lgs_a1, lgs_a2, lgs_b)

use enth_temp_omega_m, only : enth_fct_temp_omega

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class),  intent(in)    :: grd
type(sico_par_class),   intent(in)    :: par
type(sico_timer_class), intent(in)    :: tmr
integer, intent(in) :: i, j
real(wp),     intent(in) :: ct1(0:100), ct2(0:100), ct3(0:100), &
                            ct4(0:100), ce5(0:100), ce6(0:100), &
                            ce7(0:100)
real(wp),     intent(in) :: ct1_sg(0:100), ct2_sg(0:100), &
                            ct3_sg(0:100), ct4_sg(0:100), &
                            ccbe1, ccb2, ccb3, ccb4, &
                            adv_vert_sg(0:100), abs_adv_vert_sg(0:100)
real(wp),     intent(in) :: dtt_dxi, dtt_deta

real(wp),    intent(out) :: lgs_a0(0:200), &
                            lgs_a1(0:200), &
                            lgs_a2(0:200), &
                            lgs_b(0:200)

integer :: kc, kr
real(wp) :: vx_c_help, vy_c_help
real(wp) :: adv_vert_help

!-------- Initialisation --------

lgs_a0 = 0.0_wp
lgs_a1 = 0.0_wp
lgs_a2 = 0.0_wp
lgs_b  = 0.0_wp

!-------- Actual computation --------

kr=grd%KRMAX
kc=0
lgs_a1(kc) = -ccbe1
lgs_a2(kc) =  ccbe1
lgs_b(kc)  =  ccb2*(st%temp_r_neu(kr,j,i)-st%temp_r_neu(kr-1,j,i)) + ccb3 + ccb4

do kc=1, grd%KCMAX-1

if (par%adv_vert==1) then

   lgs_a0(kc) = -0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) = 1.0_wp+ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) = 0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ce5(kc)*ce6(kc)

else if (par%adv_vert==2) then

   lgs_a0(kc) &
         = -0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) &
         = 1.0_wp &
           +0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           +ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) &
         =  0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           -ce5(kc)*ce6(kc)

else if (par%adv_vert==3) then

   adv_vert_help = 0.5_wp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) &
         = -max(adv_vert_help, 0.0_wp) &
           -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) &
         = 1.0_wp &
           +max(adv_vert_help, 0.0_wp)-min(adv_vert_help, 0.0_wp) &
           +ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) &
         =  min(adv_vert_help, 0.0_wp) &
           -ce5(kc)*ce6(kc)

endif

if (par%adv_hor==2) then

   lgs_b(kc) = st%enth_c(kc,j,i) + ce7(kc) &
       -grd%dtt_2dxi* &
          ( (st%vx_c(kc,j,i)-abs(st%vx_c(kc,j,i))) &
            *(st%enth_c(kc,j,i+1)-st%enth_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_c(kc,j,i-1)+abs(st%vx_c(kc,j,i-1))) &
            *(st%enth_c(kc,j,i)-st%enth_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_c(kc,j,i)-abs(st%vy_c(kc,j,i))) &
            *(st%enth_c(kc,j+1,i)-st%enth_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_c(kc,j-1,i)+abs(st%vy_c(kc,j-1,i))) &
            *(st%enth_c(kc,j,i)-st%enth_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

   vx_c_help = 0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1))
   vy_c_help = 0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i))

   lgs_b(kc) = st%enth_c(kc,j,i) + ce7(kc) &
       -dtt_dxi* &
          ( min(vx_c_help, 0.0_wp) &
            *(st%enth_c(kc,j,i+1)-st%enth_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_c_help, 0.0_wp) &
            *(st%enth_c(kc,j,i)-st%enth_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_c_help, 0.0_wp) &
            *(st%enth_c(kc,j+1,i)-st%enth_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_c_help, 0.0_wp) &
            *(st%enth_c(kc,j,i)-st%enth_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

end do

kc=grd%KCMAX
lgs_a0(kc) = 0.0_wp
lgs_a1(kc) = 1.0_wp
lgs_b(kc)  = enth_fct_temp_omega(st%temp_s(j,i), 0.0_wp)
                                ! zero water content at the ice surface

end subroutine calc_temp_enth_1_c

!-------------------------------------------------------------------------------
!> Computation of temperature and age for a cold ice column with the
!! enthalpy method:
!! Set-up of the equations for the age of ice.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_1_d(st,grd,par,tmr,ct1, ct2, ct3, ct4, ci1, ci2, &
                              ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                              adv_vert_sg, abs_adv_vert_sg, &
                              dtt_dxi, dtt_deta, &
                              i, j, &
                              lgs_a0, lgs_a1, lgs_a2, lgs_b)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class),  intent(in)    :: grd
type(sico_par_class),   intent(in)    :: par
type(sico_timer_class), intent(in)    :: tmr
integer, intent(in) :: i, j
real(wp),     intent(in) :: ct1(0:100), ct2(0:100), ct3(0:100), &
                            ct4(0:100), ci1(0:100), ci2(0:100)
real(wp),     intent(in) :: ct1_sg(0:100), ct2_sg(0:100), &
                            ct3_sg(0:100), ct4_sg(0:100), &
                            adv_vert_sg(0:100), abs_adv_vert_sg(0:100)
real(wp),     intent(in) :: dtt_dxi, dtt_deta

real(wp),    intent(out) :: lgs_a0(0:200), &
                            lgs_a1(0:200), &
                            lgs_a2(0:200), &
                            lgs_b(0:200)

integer :: kc
real(wp) :: vx_c_help, vy_c_help
real(wp) :: adv_vert_help

!-------- Initialisation --------

lgs_a0 = 0.0_wp
lgs_a1 = 0.0_wp
lgs_a2 = 0.0_wp
lgs_b  = 0.0_wp

!-------- Actual computation --------

kc=0                                                 ! adv_vert_sg(0) <= 0
lgs_a1(kc) = 1.0_wp - min(adv_vert_sg(kc), 0.0_wp)   ! (directed downward)
lgs_a2(kc) = min(adv_vert_sg(kc), 0.0_wp)            ! assumed/enforced

if (par%adv_hor==2) then

lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -grd%dtt_2dxi* &
          ( (st%vx_c(kc,j,i)-abs(st%vx_c(kc,j,i))) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_c(kc,j,i-1)+abs(st%vx_c(kc,j,i-1))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_c(kc,j,i)-abs(st%vy_c(kc,j,i))) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_c(kc,j-1,i)+abs(st%vy_c(kc,j-1,i))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

vx_c_help = 0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1))
vy_c_help = 0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i))

lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -dtt_dxi* &
          ( min(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

do kc=1, grd%KCMAX-1

if (par%adv_vert==1) then

   lgs_a0(kc) = -0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ci1(kc)*ci2(kc-1)
   lgs_a1(kc) = 1.0_wp+ci1(kc)*(ci2(kc)+ci2(kc-1))
   lgs_a2(kc) = 0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ci1(kc)*ci2(kc)

else if (par%adv_vert==2) then

   lgs_a0(kc) = -0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1))
   lgs_a1(kc) =  1.0_wp &
                +0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
                -0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  )
   lgs_a2(kc) =  0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  )

else if (par%adv_vert==3) then

   adv_vert_help = 0.5_wp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) = -max(adv_vert_help, 0.0_wp)
   lgs_a1(kc) =  1.0_wp &
                +max(adv_vert_help, 0.0_wp)-min(adv_vert_help, 0.0_wp)
   lgs_a2(kc) =  min(adv_vert_help, 0.0_wp)

endif

if (par%adv_hor==2) then

   lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -grd%dtt_2dxi* &
          ( (st%vx_c(kc,j,i)-abs(st%vx_c(kc,j,i))) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_c(kc,j,i-1)+abs(st%vx_c(kc,j,i-1))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_c(kc,j,i)-abs(st%vy_c(kc,j,i))) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_c(kc,j-1,i)+abs(st%vy_c(kc,j-1,i))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

   vx_c_help = 0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1))
   vy_c_help = 0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i))

   lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -dtt_dxi* &
          ( min(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

end do

kc=grd%KCMAX
if (st%as_perp(j,i) >= 0.0_wp) then
   lgs_a0(kc) = 0.0_wp
   lgs_a1(kc) = 1.0_wp
   lgs_b(kc)  = 0.0_wp
else
   lgs_a0(kc) = -max(adv_vert_sg(kc-1), 0.0_wp)
   lgs_a1(kc) = 1.0_wp + max(adv_vert_sg(kc-1), 0.0_wp)
                       ! adv_vert_sg(grd%KCMAX-1) >= 0 (directed upward)
                       ! assumed/enforced
if (par%adv_hor==2) then

   lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -grd%dtt_2dxi* &
          ( (st%vx_c(kc,j,i)-abs(st%vx_c(kc,j,i))) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_c(kc,j,i-1)+abs(st%vx_c(kc,j,i-1))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_c(kc,j,i)-abs(st%vy_c(kc,j,i))) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_c(kc,j-1,i)+abs(st%vy_c(kc,j-1,i))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

   vx_c_help = 0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1))
   vy_c_help = 0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i))

   lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -dtt_dxi* &
          ( min(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

end if

end subroutine calc_temp_enth_1_d

!-------------------------------------------------------------------------------
!> Computation of temperature and age for an ice column with a temperate base
!! with the enthalpy method.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_2(st,grd,tmr,par,i, j)

use sico_maths_m,      only : tri_sle
use enth_temp_omega_m, only : enth_fct_temp_omega, &
                              temp_fct_enth, omega_fct_enth

implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class),  intent(in)    :: grd
  type(sico_timer_class), intent(in)    :: tmr
  type(sico_par_class), intent(in)    :: par
integer, intent(in) :: i, j

integer :: kc, kt, kr
real(wp) :: ct1(0:100), ct2(0:100), ct3(0:100), ct4(0:100), &
            ce5(0:100), ce6(0:100), ce7(0:100), ctr1, clb1
real(wp) :: ct1_sg(0:100), ct2_sg(0:100), ct3_sg(0:100), &
            ct4_sg(0:100), adv_vert_sg(0:100), abs_adv_vert_sg(0:100)
real(wp) :: ci1(0:100), ci2(0:100)
real(wp) :: cqtlde(0:100), cm3(0:100)
real(wp) :: dtt_dxi, dtt_deta
real(wp) :: temp_c_val(0:100), omega_c_val(0:100)
real(wp) :: lgs_a0(0:200), &
            lgs_a1(0:200), &
            lgs_a2(0:200), &
            lgs_x(0:200), &
            lgs_b(0:200)

real(wp), parameter :: eps_omega=1.0e-12_wp

!-------- Check for boundary points --------

if ((i == 0).or.(i == grd%IMAX).or.(j == 0).or.(j == grd%JMAX)) &
   stop ' >>> calc_temp_enth_2: Boundary points not allowed.'

!-------- Abbreviations --------

call calc_temp_enth_2_a1(st,grd,par,tmr,i, j, &
                         ct1, ct2, ct3, ct4, ce5, ctr1, clb1, &
                         ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                         adv_vert_sg, abs_adv_vert_sg, &
                         ci1, cqtlde, dtt_dxi, dtt_deta)

do kc=0, grd%KCMAX
   temp_c_val(kc)  = st%temp_c(kc,j,i)
   omega_c_val(kc) = st%omega_c(kc,j,i)
end do

call calc_temp_enth_2_a2(st,grd,tmr,par, temp_c_val, omega_c_val, &
                         i, j, ce6, ce7, ci2, cm3)

!-------- Computation of the bedrock temperature --------

!  ------ Set-up of the the equations

call calc_temp_enth_2_b(st,grd,par,tmr,ctr1, clb1, i, j, &
                        lgs_a0, lgs_a1, lgs_a2, lgs_b)

!  ------ Solution of the system of linear equations

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KRMAX)

!  ------ Assignment of the result

do kr=0, grd%KRMAX
   st%temp_r_neu(kr,j,i) = lgs_x(kr)
end do

!-------- Computation of the ice enthalpy (predictor step) --------

!  ------ Set-up of the equations

call calc_temp_enth_2_c(st,grd,par,tmr,ct1, ct2, ct3, ct4, ce5, ce6, ce7, &
                        ct1_sg, ct2_sg, ct3_sg, ct4_sg, cm3, &
                        adv_vert_sg, abs_adv_vert_sg, &
                        dtt_dxi, dtt_deta, &
                        i, j, 0, &
                        lgs_a0, lgs_a1, lgs_a2, lgs_b)

!  ------ Solution of the system of linear equations

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KCMAX)

!  ------ Assignment of the result

do kc=0, grd%KCMAX
   st%enth_c_neu(kc,j,i)  = lgs_x(kc)
   st%temp_c_neu(kc,j,i)  = temp_fct_enth(st%enth_c_neu(kc,j,i), st%temp_c_m(kc,j,i))
   st%omega_c_neu(kc,j,i) = omega_fct_enth(st%enth_c_neu(kc,j,i), st%temp_c_m(kc,j,i))
end do

!  ------ Determination of the CTS

st%kc_cts_neu(j,i) = 0

do kc=1, grd%KCMAX-1
   if (st%omega_c_neu(kc,j,i) > eps_omega) then
      st%kc_cts_neu(j,i) = kc
   else
      exit
   end if
end do

!-------- Computation of the ice enthalpy
!         (corrector step for the cold-ice domain only
!         in order to fulfull the transition condition at the CTS) --------

if (par%calcmod==3) then  !/* ENTM scheme */

  if (st%kc_cts_neu(j,i) > 0) then

    !  ------ Update of the abbreviations where needed

    do kc=0, grd%KCMAX
      temp_c_val(kc)  = st%temp_c_neu(kc,j,i)
      omega_c_val(kc) = st%omega_c_neu(kc,j,i)
    end do

    call calc_temp_enth_2_a2(st,grd,tmr,par, temp_c_val, omega_c_val, &
      i, j, ce6, ce7, ci2, cm3)

    !  ------ Set-up of the equations

    call calc_temp_enth_2_c(st,grd,par,tmr,ct1, ct2, ct3, ct4, ce5, ce6, ce7, &
      ct1_sg, ct2_sg, ct3_sg, ct4_sg, cm3, &
      adv_vert_sg, abs_adv_vert_sg, &
      dtt_dxi, dtt_deta, &
      i, j, st%kc_cts_neu(j,i), &
      lgs_a0, lgs_a1, lgs_a2, lgs_b)

    !  ------ Solution of the system of linear equations

    call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KCMAX)

    !  ------ Assignment of the result

    kc=st%kc_cts_neu(j,i)
    st%enth_c_neu(kc,j,i)  = lgs_x(kc)
    st%temp_c_neu(kc,j,i)  = temp_fct_enth(st%enth_c_neu(kc,j,i), st%temp_c_m(kc,j,i))
    st%omega_c_neu(kc,j,i) = omega_fct_enth(st%enth_c_neu(kc,j,i), st%temp_c_m(kc,j,i))

    do kc=st%kc_cts_neu(j,i)+1, grd%KCMAX
      st%enth_c_neu(kc,j,i)  = lgs_x(kc)
      st%temp_c_neu(kc,j,i)  = temp_fct_enth(st%enth_c_neu(kc,j,i), st%temp_c_m(kc,j,i))
      st%omega_c_neu(kc,j,i) = 0.0_wp   ! cold-ice domain
    end do

  end if

else if (par%calcmod==2) then  ! /* ENTC scheme */

  !!! continue   ! no corrector step

else
  stop ' >>> calc_temp_enth_2: CALCMOD must be either 2 or 3!'
endif

!-------- Water drainage from temperate ice (if existing) --------

st%Q_tld(j,i) = 0.0_wp

do kc=0, st%kc_cts_neu(j,i)

   if (st%omega_c_neu(kc,j,i) > OMEGA_MAX) then

      st%Q_tld(j,i) = st%Q_tld(j,i) + cqtlde(kc)*(st%omega_c_neu(kc,j,i)-OMEGA_MAX)

      st%omega_c_neu(kc,j,i) = OMEGA_MAX
      st%enth_c_neu(kc,j,i)  = enth_fct_temp_omega(st%temp_c_neu(kc,j,i), OMEGA_MAX)

   end if

end do

!-------- Set enthalpy and water content in the redundant,
!         lower (kt) ice layer to the value at the ice base --------

do kt=0, grd%KTMAX
   st%enth_t_neu(kt,j,i)  = st%enth_c_neu(0,j,i)
   st%omega_t_neu(kt,j,i) = st%omega_c_neu(0,j,i)
end do

!-------- Computation of the age of ice --------

!  ------ Set-up of the the equations

call calc_temp_enth_2_d(st,grd,par,tmr,ct1, ct2, ct3, ct4, ci1, ci2, &
                        ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                        adv_vert_sg, abs_adv_vert_sg, &
                        dtt_dxi, dtt_deta, &
                        i, j, &
                        lgs_a0, lgs_a1, lgs_a2, lgs_b)

!  ------ Solution of the system of linear equations

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KCMAX)

!  ------ Assignment of the result
!         restriction to interval [0, AGE_MAX yr]

do kc=0, grd%KCMAX

   st%age_c_neu(kc,j,i) = lgs_x(kc)

   if (st%age_c_neu(kc,j,i) < (par%age_min*sec_year)) &
                           st%age_c_neu(kc,j,i) = 0.0_wp
   if (st%age_c_neu(kc,j,i) > (par%age_max*sec_year)) &
                           st%age_c_neu(kc,j,i) = par%age_max*sec_year

end do

!-------- Age of the ice in the redundant, lower (kt) ice layer --------

do kt=0, grd%KTMAX
   st%age_t_neu(kt,j,i) = st%age_c_neu(0,j,i)
end do

end subroutine calc_temp_enth_2

!-------------------------------------------------------------------------------
!> Computation of temperature and age for an ice column with a temperate base
!! with the enthalpy method: Abbreviations I.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_2_a1(st,grd,par,tmr,i, j, &
                               ct1, ct2, ct3, ct4, ce5, ctr1, clb1, &
                               ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                               adv_vert_sg, abs_adv_vert_sg, &
                               ci1, cqtlde, dtt_dxi, dtt_deta)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class),  intent(in)    :: grd
type(sico_par_class),   intent(in)    :: par
type(sico_timer_class), intent(in)    :: tmr
integer, intent(in) :: i, j

real(wp),    intent(out) :: ct1(0:100), ct2(0:100), ct3(0:100), &
                            ct4(0:100), ce5(0:100), &
                            ctr1, clb1
real(wp),    intent(out) :: ct1_sg(0:100), ct2_sg(0:100), &
                            ct3_sg(0:100), ct4_sg(0:100), &
                            adv_vert_sg(0:100), abs_adv_vert_sg(0:100)
real(wp),    intent(out) :: ci1(0:100), cqtlde(0:100)
real(wp),    intent(out) :: dtt_dxi, dtt_deta

integer :: kc

!-------- Initialisation --------

ct1             = 0.0_wp
ct2             = 0.0_wp
ct3             = 0.0_wp
ct4             = 0.0_wp
ce5             = 0.0_wp
ctr1            = 0.0_wp
clb1            = 0.0_wp
ct1_sg          = 0.0_wp
ct2_sg          = 0.0_wp
ct3_sg          = 0.0_wp
ct4_sg          = 0.0_wp
adv_vert_sg     = 0.0_wp
abs_adv_vert_sg = 0.0_wp
ci1             = 0.0_wp
cqtlde          = 0.0_wp
dtt_dxi         = 0.0_wp
dtt_deta        = 0.0_wp

!-------- Actual computation --------

ctr1 = grd%atr1
clb1 = grd%alb1*st%q_geo(j,i)

if (par%adv_vert==1) then

do kc=1, grd%KCMAX-1
   ct1(kc) = grd%at1(kc)/st%H_c(j,i)*0.5_wp*(st%vz_c(kc,j,i)+st%vz_c(kc-1,j,i))
end do

kc=0
ct1_sg(kc) = 0.5_wp*(grd%at1(kc)+grd%at1(kc+1))/st%H_c(j,i)*st%vz_c(kc,j,i)
             ! only needed for kc=0 ...
kc=grd%KCMAX-1
ct1_sg(kc) = 0.5_wp*(grd%at1(kc)+grd%at1(kc+1))/st%H_c(j,i)*st%vz_c(kc,j,i)
             ! ... and kc=grd%KCMAX-1

else if (par%adv_vert==2 .or. par%adv_vert==3) then

do kc=0, grd%KCMAX-1
   ct1_sg(kc) = 0.5_wp*(grd%at1(kc)+grd%at1(kc+1))/st%H_c(j,i)*st%vz_c(kc,j,i)
end do

endif

do kc=0, grd%KCMAX

   ct2(kc) = ( grd%at2_1(kc)*st%dzm_dtau(j,i) &
           +grd%at2_2(kc)*st%dH_c_dtau(j,i) )/st%H_c(j,i)
   ct3(kc) = ( grd%at3_1(kc)*st%dzm_dxi_g(j,i) &
           +grd%at3_2(kc)*st%dH_c_dxi_g(j,i) )/st%H_c(j,i) &
          *0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1)) *grd%insq_g11_g(j,i)
   ct4(kc) = ( grd%at4_1(kc)*st%dzm_deta_g(j,i) &
            +grd%at4_2(kc)*st%dH_c_deta_g(j,i) )/st%H_c(j,i) &
          *0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i)) *grd%insq_g22_g(j,i)
   ce5(kc) = grd%at5(kc)/st%H_c(j,i)
   ci1(kc) = grd%ai1(kc)/st%H_c(j,i)
   cqtlde(kc) = grd%aqtlde(kc)*st%H_c(j,i)

end do

if (par%adv_vert==1) then

kc=0
ct2_sg(kc) = 0.5_wp*(ct2(kc)+ct2(kc+1))
ct3_sg(kc) = 0.5_wp*(ct3(kc)+ct3(kc+1))
ct4_sg(kc) = 0.5_wp*(ct4(kc)+ct4(kc+1))
adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))   ! only needed for kc=0 ...
kc=grd%KCMAX-1
ct2_sg(kc) = 0.5_wp*(ct2(kc)+ct2(kc+1))
ct3_sg(kc) = 0.5_wp*(ct3(kc)+ct3(kc+1))
ct4_sg(kc) = 0.5_wp*(ct4(kc)+ct4(kc+1))
adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))   ! ... and kc=grd%KCMAX-1

else if (par%adv_vert==2 .or. par%adv_vert==3) then

do kc=0, grd%KCMAX-1
   ct2_sg(kc) = 0.5_wp*(ct2(kc)+ct2(kc+1))
   ct3_sg(kc) = 0.5_wp*(ct3(kc)+ct3(kc+1))
   ct4_sg(kc) = 0.5_wp*(ct4(kc)+ct4(kc+1))
   adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
   abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))
end do

endif

if (par%adv_hor==3) then
dtt_dxi  = 2.0_wp*grd%dtt_2dxi
dtt_deta = 2.0_wp*grd%dtt_2deta
endif

end subroutine calc_temp_enth_2_a1

!-------------------------------------------------------------------------------
!> Computation of temperature and age for an ice column with a temperate base
!! with the enthalpy method: Abbreviations II.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_2_a2(st,grd,tmr,par, temp_c_val, omega_c_val, &
                               i, j, ce6, ce7, ci2, cm3)

use ice_material_properties_m, only : ratefac_c_t, kappa_val, c_val, &
                                      creep, viscosity

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class),  intent(in)    :: grd
type(sico_timer_class), intent(in)    :: tmr
type(sico_par_class), intent(in)    :: par
integer, intent(in) :: i, j
real(wp),     intent(in) :: temp_c_val(0:100), omega_c_val(0:100)

real(wp),    intent(out) :: ce6(0:100), ce7(0:100), ci2(0:100), &
                            cm3(0:100)

integer :: kc
real(wp) :: temp_c_help(0:100)

!-------- Initialisation --------

ce6 = 0.0_wp
ce7 = 0.0_wp
ci2 = 0.0_wp
cm3 = 0.0_wp

!-------- Actual computation --------

do kc=0, grd%KCMAX

  if (par%dynamics==2) then

    if (.not.st%flag_shelfy_stream(j,i)) then
      ce7(kc) = grd%at7 &
        *st%enh_c(kc,j,i) &
        *ratefac_c_t(temp_c_val(kc), omega_c_val(kc), st%temp_c_m(kc,j,i)) &
        *creep(par%fin_visc,par%flow_law,st%sigma_c(kc,j,i)) &
        *st%sigma_c(kc,j,i)*st%sigma_c(kc,j,i)
    else
      ce7(kc) = 2.0_wp*grd%at7 &
        *viscosity(par%fin_visc,par%flow_law,st%de_c(kc,j,i), &
        temp_c_val(kc), st%temp_c_m(kc,j,i), omega_c_val(kc), &
        st%enh_c(kc,j,i), 2) &
        *st%de_c(kc,j,i)**2
    end if

  else   ! par%dynamics.ne.2

    ce7(kc) = grd%at7 &
      *st%enh_c(kc,j,i) &
      *ratefac_c_t(temp_c_val(kc), omega_c_val(kc), st%temp_c_m(kc,j,i)) &
      *creep(par%fin_visc,par%flow_law,st%sigma_c(kc,j,i)) &
      *st%sigma_c(kc,j,i)*st%sigma_c(kc,j,i)

  endif

   cm3(kc) = grd%am3(kc)*st%H_c(j,i)*c_val(temp_c_val(kc))

end do

do kc=0, st%kc_cts_neu(j,i)-1   ! temperate layer
   ce6(kc) = grd%at6(kc) &
             *NUE/st%H_c(j,i)
   ci2(kc) = grd%ai2(kc)/st%H_c(j,i)
end do

do kc=st%kc_cts_neu(j,i), grd%KCMAX-1   ! cold layer
   temp_c_help(kc) = 0.5_wp*(temp_c_val(kc)+temp_c_val(kc+1))
   ce6(kc) = grd%at6(kc) &
             *kappa_val(temp_c_help(kc))/(c_val(temp_c_help(kc))*st%H_c(j,i))
   ci2(kc) = grd%ai2(kc)/st%H_c(j,i)
end do

end subroutine calc_temp_enth_2_a2

!-------------------------------------------------------------------------------
!> Computation of temperature and age for an ice column with a temperate base
!! with the enthalpy method:
!! Set-up of the equations for the bedrock temperature.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_2_b(st,grd,par,tmr,ctr1, clb1, i, j, &
                              lgs_a0, lgs_a1, lgs_a2, lgs_b)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class),  intent(in)    :: grd
type(sico_par_class),   intent(in)    :: par
type(sico_timer_class), intent(in)    :: tmr
integer, intent(in) :: i, j
real(wp),     intent(in) :: ctr1, clb1

real(wp),    intent(out) :: lgs_a0(0:200), &
                            lgs_a1(0:200), &
                            lgs_a2(0:200), &
                            lgs_b(0:200)

integer :: kr

!-------- Initialisation --------

lgs_a0 = 0.0_wp
lgs_a1 = 0.0_wp
lgs_a2 = 0.0_wp
lgs_b  = 0.0_wp

!-------- Actual computation --------

kr=0
lgs_a1(kr) = 1.0_wp
lgs_a2(kr) = -1.0_wp
lgs_b(kr)    = clb1

if (par%q_litho==1) then
!   (coupled heat-conducting bedrock)

do kr=1, grd%KRMAX-1
   lgs_a0(kr) = - ctr1
   lgs_a1(kr) = 1.0_wp + 2.0_wp*ctr1
   lgs_a2(kr) = - ctr1
   lgs_b(kr)    = st%temp_r(kr,j,i)
end do

else if (par%q_litho==0) then
!   (no coupled heat-conducting bedrock)

do kr=1, grd%KRMAX-1
   lgs_a0(kr) = 1.0_wp
   lgs_a1(kr) = 0.0_wp
   lgs_a2(kr) = -1.0_wp
   lgs_b(kr)  = 2.0_wp*clb1
end do

endif

kr=grd%KRMAX
lgs_a0(kr) = 0.0_wp
lgs_a1(kr) = 1.0_wp
lgs_b(kr)  = st%temp_t_m(0,j,i)

end subroutine calc_temp_enth_2_b

!-------------------------------------------------------------------------------
!> Computation of temperature and age for an ice column with a temperate base
!! with the enthalpy method:
!! Set-up of the equations for the ice enthalpy.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_2_c(st,grd,par,tmr,ct1, ct2, ct3, ct4, ce5, ce6, ce7, &
                              ct1_sg, ct2_sg, ct3_sg, ct4_sg, cm3, &
                              adv_vert_sg, abs_adv_vert_sg, &
                              dtt_dxi, dtt_deta, &
                              i, j, kcmin, &
                              lgs_a0, lgs_a1, lgs_a2, lgs_b)

use enth_temp_omega_m, only : enth_fct_temp_omega

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class),  intent(in)    :: grd
type(sico_par_class),   intent(in)    :: par
type(sico_timer_class), intent(in)    :: tmr
integer, intent(in) :: i, j
integer, intent(in) :: kcmin
real(wp),     intent(in) :: ct1(0:100), ct2(0:100), ct3(0:100), &
                            ct4(0:100), ce5(0:100), ce6(0:100), &
                            ce7(0:100)
real(wp),     intent(in) :: ct1_sg(0:100), ct2_sg(0:100), &
                            ct3_sg(0:100), ct4_sg(0:100), &
                            adv_vert_sg(0:100), abs_adv_vert_sg(0:100)
real(wp),     intent(in) :: cm3(0:100)
real(wp),     intent(in) :: dtt_dxi, dtt_deta

real(wp),    intent(out) :: lgs_a0(0:200), &
                            lgs_a1(0:200), &
                            lgs_a2(0:200), &
                            lgs_b(0:200)

integer :: kc
real(wp) :: vx_c_help, vy_c_help
real(wp) :: adv_vert_help

!-------- Initialisation --------

lgs_a0 = 0.0_wp
lgs_a1 = 0.0_wp
lgs_a2 = 0.0_wp
lgs_b  = 0.0_wp

!-------- Actual computation --------

if (kcmin == 0) then   ! predictor step

   kc=0

   if (st%kc_cts_neu(j,i) == 0) then   ! temperate base without temperate layer

      lgs_a1(kc) = 1.0_wp
      lgs_a2(kc) = 0.0_wp
      lgs_b(kc)  = enth_fct_temp_omega(st%temp_c_m(kc,j,i), 0.0_wp)

   else   ! st%kc_cts_neu(j,i) > 0, temperate base with temperate layer

      lgs_a1(kc) =  1.0_wp
      lgs_a2(kc) = -1.0_wp
      lgs_b(kc)  =  0.0_wp

   end if

else   ! kcmin > 0, corrector step

   kc=0

   lgs_a1(kc) = 1.0_wp   ! dummy
   lgs_a2(kc) = 0.0_wp   ! setting,
   lgs_b(kc)  = 0.0_wp   ! not needed

   do kc=1, kcmin-1

      lgs_a0(kc) = 0.0_wp   ! dummy
      lgs_a1(kc) = 1.0_wp   ! setting,
      lgs_a2(kc) = 0.0_wp   ! not
      lgs_b(kc)  = 0.0_wp   ! needed

   end do

   kc=kcmin

   lgs_a0(kc) = 0.0_wp
   lgs_a1(kc) = 1.0_wp
   lgs_a2(kc) = -1.0_wp
   lgs_b(kc)  = -cm3(kc)

end if

do kc=kcmin+1, grd%KCMAX-1

if (par%adv_vert==1) then

   lgs_a0(kc) = -0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) = 1.0_wp+ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) = 0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ce5(kc)*ce6(kc)

else if (par%adv_vert==2) then

   lgs_a0(kc) &
         = -0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) &
         = 1.0_wp &
           +0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           +ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) &
         =  0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           -ce5(kc)*ce6(kc)

else if (par%adv_vert==3) then

   adv_vert_help = 0.5_wp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) &
         = -max(adv_vert_help, 0.0_wp) &
           -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) &
         = 1.0_wp &
           +max(adv_vert_help, 0.0_wp)-min(adv_vert_help, 0.0_wp) &
           +ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) &
         =  min(adv_vert_help, 0.0_wp) &
           -ce5(kc)*ce6(kc)

endif

if (par%adv_hor==2) then

   lgs_b(kc) = st%enth_c(kc,j,i) + ce7(kc) &
       -grd%dtt_2dxi* &
          ( (st%vx_c(kc,j,i)-abs(st%vx_c(kc,j,i))) &
            *(st%enth_c(kc,j,i+1)-st%enth_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_c(kc,j,i-1)+abs(st%vx_c(kc,j,i-1))) &
            *(st%enth_c(kc,j,i)-st%enth_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_c(kc,j,i)-abs(st%vy_c(kc,j,i))) &
            *(st%enth_c(kc,j+1,i)-st%enth_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_c(kc,j-1,i)+abs(st%vy_c(kc,j-1,i))) &
            *(st%enth_c(kc,j,i)-st%enth_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

   vx_c_help = 0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1))
   vy_c_help = 0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i))

   lgs_b(kc) = st%enth_c(kc,j,i) + ce7(kc) &
       -dtt_dxi* &
          ( min(vx_c_help, 0.0_wp) &
            *(st%enth_c(kc,j,i+1)-st%enth_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_c_help, 0.0_wp) &
            *(st%enth_c(kc,j,i)-st%enth_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_c_help, 0.0_wp) &
            *(st%enth_c(kc,j+1,i)-st%enth_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_c_help, 0.0_wp) &
            *(st%enth_c(kc,j,i)-st%enth_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

end do

kc=grd%KCMAX
lgs_a0(kc) = 0.0_wp
lgs_a1(kc) = 1.0_wp
lgs_b(kc)  = enth_fct_temp_omega(st%temp_s(j,i), 0.0_wp)
                                ! zero water content at the ice surface

end subroutine calc_temp_enth_2_c

!-------------------------------------------------------------------------------
!> Computation of temperature and age for an ice column with a temperate base
!! with the enthalpy method:
!! Set-up of the equations for the age of ice.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_2_d(st,grd,par,tmr,ct1, ct2, ct3, ct4, ci1, ci2, &
                              ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                              adv_vert_sg, abs_adv_vert_sg, &
                              dtt_dxi, dtt_deta, &
                              i, j, &
                              lgs_a0, lgs_a1, lgs_a2, lgs_b)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class),  intent(in)    :: grd
type(sico_par_class),   intent(in)    :: par
type(sico_timer_class), intent(in)    :: tmr
integer, intent(in) :: i, j
real(wp),     intent(in) :: ct1(0:100), ct2(0:100), ct3(0:100), &
                            ct4(0:100), ci1(0:100), ci2(0:100)
real(wp),     intent(in) :: ct1_sg(0:100), ct2_sg(0:100), &
                            ct3_sg(0:100), ct4_sg(0:100), &
                            adv_vert_sg(0:100), abs_adv_vert_sg(0:100)
real(wp),     intent(in) :: dtt_dxi, dtt_deta

real(wp),    intent(out) :: lgs_a0(0:200), &
                            lgs_a1(0:200), &
                            lgs_a2(0:200), &
                            lgs_b(0:200)

integer :: kc
real(wp) :: vx_c_help, vy_c_help
real(wp) :: adv_vert_help

!-------- Initialisation --------

lgs_a0 = 0.0_wp
lgs_a1 = 0.0_wp
lgs_a2 = 0.0_wp
lgs_b  = 0.0_wp

!-------- Actual computation --------

kc=0                                                 ! adv_vert_sg(0) <= 0
lgs_a1(kc) = 1.0_wp - min(adv_vert_sg(kc), 0.0_wp)   ! (directed downward)
lgs_a2(kc) = min(adv_vert_sg(kc), 0.0_wp)            ! assumed/enforced

if (par%adv_hor==2) then

lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -grd%dtt_2dxi* &
          ( (st%vx_c(kc,j,i)-abs(st%vx_c(kc,j,i))) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_c(kc,j,i-1)+abs(st%vx_c(kc,j,i-1))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_c(kc,j,i)-abs(st%vy_c(kc,j,i))) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_c(kc,j-1,i)+abs(st%vy_c(kc,j-1,i))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

vx_c_help = 0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1))
vy_c_help = 0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i))

lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -dtt_dxi* &
          ( min(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

do kc=1, grd%KCMAX-1

if (par%adv_vert==1) then

   lgs_a0(kc) = -0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ci1(kc)*ci2(kc-1)
   lgs_a1(kc) = 1.0_wp+ci1(kc)*(ci2(kc)+ci2(kc-1))
   lgs_a2(kc) = 0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ci1(kc)*ci2(kc)

else if (par%adv_vert==2) then

   lgs_a0(kc) = -0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1))
   lgs_a1(kc) =  1.0_wp &
                +0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
                -0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  )
   lgs_a2(kc) =  0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  )

else if (par%adv_vert==3) then

   adv_vert_help = 0.5_wp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) = -max(adv_vert_help, 0.0_wp)
   lgs_a1(kc) =  1.0_wp &
                +max(adv_vert_help, 0.0_wp)-min(adv_vert_help, 0.0_wp)
   lgs_a2(kc) =  min(adv_vert_help, 0.0_wp)

endif

if (par%adv_hor==2) then

   lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -grd%dtt_2dxi* &
          ( (st%vx_c(kc,j,i)-abs(st%vx_c(kc,j,i))) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_c(kc,j,i-1)+abs(st%vx_c(kc,j,i-1))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_c(kc,j,i)-abs(st%vy_c(kc,j,i))) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_c(kc,j-1,i)+abs(st%vy_c(kc,j-1,i))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

   vx_c_help = 0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1))
   vy_c_help = 0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i))

   lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -dtt_dxi* &
          ( min(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

end do

kc=grd%KCMAX
if (st%as_perp(j,i) >= 0.0_wp) then
   lgs_a0(kc) = 0.0_wp
   lgs_a1(kc) = 1.0_wp
   lgs_b(kc)  = 0.0_wp
else
   lgs_a0(kc) = -max(adv_vert_sg(kc-1), 0.0_wp)
   lgs_a1(kc) = 1.0_wp + max(adv_vert_sg(kc-1), 0.0_wp)
                       ! adv_vert_sg(grd%KCMAX-1) >= 0 (directed upward)
                       ! assumed/enforced
if (par%adv_hor==2) then

   lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -grd%dtt_2dxi* &
          ( (st%vx_c(kc,j,i)-abs(st%vx_c(kc,j,i))) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_c(kc,j,i-1)+abs(st%vx_c(kc,j,i-1))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_c(kc,j,i)-abs(st%vy_c(kc,j,i))) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_c(kc,j-1,i)+abs(st%vy_c(kc,j-1,i))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

   vx_c_help = 0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1))
   vy_c_help = 0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i))

   lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -dtt_dxi* &
          ( min(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

end if

end subroutine calc_temp_enth_2_d

!-------------------------------------------------------------------------------
!> Computation of temperature, age, water content and enthalpy for an
!! ice-free column.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_r(st,grd,par,tmr,i, j)

use sico_maths_m,      only : tri_sle
use enth_temp_omega_m, only : enth_fct_temp_omega

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class),  intent(in)    :: grd
type(sico_par_class),   intent(in)    :: par
type(sico_timer_class), intent(in)    :: tmr
integer, intent(in) :: i, j

integer :: kc, kt, kr
real(wp) :: ctr1, clb1
real(wp) :: lgs_a0(0:200), &
            lgs_a1(0:200), &
            lgs_a2(0:200), &
            lgs_x(0:200), &
            lgs_b(0:200)
real(wp) :: enth_val

!-------- Abbreviations --------

ctr1 = grd%atr1
clb1 = grd%alb1*st%q_geo(j,i)

!-------- Set up the equations for the bedrock temperature --------

kr=0
lgs_a1(kr) = 1.0_wp
lgs_a2(kr) = -1.0_wp
lgs_b(kr)    = clb1

if (par%q_litho==1) then
!   (coupled heat-conducting bedrock)

do kr=1, grd%KRMAX-1
   lgs_a0(kr) = - ctr1
   lgs_a1(kr) = 1.0_wp + 2.0_wp*ctr1
   lgs_a2(kr) = - ctr1
   lgs_b(kr)    = st%temp_r(kr,j,i)
end do

else if (par%q_litho==0) then
!   (no coupled heat-conducting bedrock)

do kr=1, grd%KRMAX-1
   lgs_a0(kr) = 1.0_wp
   lgs_a1(kr) = 0.0_wp
   lgs_a2(kr) = -1.0_wp
   lgs_b(kr)  = 2.0_wp*clb1
end do

endif

kr=grd%KRMAX
lgs_a0(kr) = 0.0_wp
lgs_a1(kr) = 1.0_wp
lgs_b(kr)   = st%temp_g(j,i)

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KRMAX)

!-------- Assign the result --------

do kr=0, grd%KRMAX
   st%temp_r_neu(kr,j,i) = lgs_x(kr)
end do

!-------- Water content, age and enthalpy
!                        in the non-existing lower (kt) ice layer --------

enth_val = enth_fct_temp_omega(st%temp_g(j,i), 0.0_wp)

do kt=0, grd%KTMAX
   st%omega_t_neu(kt,j,i) = 0.0_wp
   st%age_t_neu(kt,j,i)   = 0.0_wp
   st%enth_t_neu(kt,j,i)  = enth_val
end do

!-------- Temperature, age, water content and enthalpy
!                      in the non-existing upper (kc) ice layer --------

do kc=0, grd%KCMAX
   st%temp_c_neu(kc,j,i)  = st%temp_g(j,i)
   st%age_c_neu(kc,j,i)   = 0.0_wp
   st%omega_c_neu(kc,j,i) = 0.0_wp
   st%enth_c_neu(kc,j,i)  = enth_val
end do

end subroutine calc_temp_enth_r

!-------------------------------------------------------------------------------
!> Computation of temperature and age for ice shelves (floating ice)
!! with the enthalpy method.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_ssa(st,grd,tmr,par, i, j)

use ice_material_properties_m, only : kappa_val, c_val, viscosity
use sico_maths_m,              only : tri_sle
use enth_temp_omega_m,         only : enth_fct_temp_omega, &
                                      temp_fct_enth, omega_fct_enth

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class),  intent(in)    :: grd
type(sico_timer_class), intent(in)    :: tmr
type(sico_par_class), intent(in)      :: par
integer, intent(in) :: i, j

integer :: kc, kt, kr
real(wp) :: ct1(0:100), ct2(0:100), ct3(0:100), ct4(0:100), &
            ce5(0:100), ce6(0:100), ce7(0:100), ctr1, clb1
real(wp) :: ct1_sg(0:100), ct2_sg(0:100), ct3_sg(0:100), &
            ct4_sg(0:100), adv_vert_sg(0:100), abs_adv_vert_sg(0:100)
real(wp) :: ci1(0:100), ci2(0:100)
real(wp) :: temp_c_help(0:100)
real(wp) :: vx_c_help, vy_c_help
real(wp) :: adv_vert_help
real(wp) :: dtt_dxi, dtt_deta
real(wp) :: lgs_a0(0:200), &
            lgs_a1(0:200), &
            lgs_a2(0:200), &
            lgs_x(0:200), &
            lgs_b(0:200)
real(wp), parameter :: zero=0.0_wp

!-------- Check for boundary points --------

if ((i == 0).or.(i == grd%IMAX).or.(j == 0).or.(j == grd%JMAX)) &
   stop ' >>> calc_temp_enth_ssa: Boundary points not allowed.'

!-------- Abbreviations --------

ctr1 = grd%atr1
clb1 = grd%alb1*st%q_geo(j,i)

if (par%adv_vert==1) then

do kc=1, grd%KCMAX-1
   ct1(kc) = grd%at1(kc)/st%H_c(j,i)*0.5_wp*(st%vz_c(kc,j,i)+st%vz_c(kc-1,j,i))
end do

kc=0
ct1_sg(kc) = 0.5_wp*(grd%at1(kc)+grd%at1(kc+1))/st%H_c(j,i)*st%vz_c(kc,j,i)
             ! only needed for kc=0 ...
kc=grd%KCMAX-1
ct1_sg(kc) = 0.5_wp*(grd%at1(kc)+grd%at1(kc+1))/st%H_c(j,i)*st%vz_c(kc,j,i)
             ! ... and kc=grd%KCMAX-1

else if (par%adv_vert==2 .or. par%adv_vert==3) then

do kc=0, grd%KCMAX-1
   ct1_sg(kc) = 0.5_wp*(grd%at1(kc)+grd%at1(kc+1))/st%H_c(j,i)*st%vz_c(kc,j,i)
end do

endif

do kc=0, grd%KCMAX

   ct2(kc) = ( grd%at2_1(kc)*st%dzm_dtau(j,i) &
           +grd%at2_2(kc)*st%dH_c_dtau(j,i) )/st%H_c(j,i)
   ct3(kc) = ( grd%at3_1(kc)*st%dzm_dxi_g(j,i) &
           +grd%at3_2(kc)*st%dH_c_dxi_g(j,i) )/st%H_c(j,i) &
          *0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1)) *grd%insq_g11_g(j,i)
   ct4(kc) = ( grd%at4_1(kc)*st%dzm_deta_g(j,i) &
            +grd%at4_2(kc)*st%dH_c_deta_g(j,i) )/st%H_c(j,i) &
          *0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i)) *grd%insq_g22_g(j,i)
   ce5(kc) = grd%at5(kc)/st%H_c(j,i)
   ce7(kc) = 2.0_wp*grd%at7 &
             *viscosity(par%fin_visc,par%flow_law,st%de_ssa(j,i), &
                        st%temp_c(kc,j,i), st%temp_c_m(kc,j,i), 0.0_wp, &
                        st%enh_c(kc,j,i), 0) &
             *st%de_ssa(j,i)**2
   ci1(kc) = grd%ai1(kc)/st%H_c(j,i)

end do

if (par%adv_vert==1) then

kc=0
ct2_sg(kc) = 0.5_wp*(ct2(kc)+ct2(kc+1))
ct3_sg(kc) = 0.5_wp*(ct3(kc)+ct3(kc+1))
ct4_sg(kc) = 0.5_wp*(ct4(kc)+ct4(kc+1))
adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))   ! only needed for kc=0 ...
kc=grd%KCMAX-1
ct2_sg(kc) = 0.5_wp*(ct2(kc)+ct2(kc+1))
ct3_sg(kc) = 0.5_wp*(ct3(kc)+ct3(kc+1))
ct4_sg(kc) = 0.5_wp*(ct4(kc)+ct4(kc+1))
adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))   ! ... and kc=grd%KCMAX-1

else if (par%adv_vert==2 .or. par%adv_vert==3) then

do kc=0, grd%KCMAX-1
   ct2_sg(kc) = 0.5_wp*(ct2(kc)+ct2(kc+1))
   ct3_sg(kc) = 0.5_wp*(ct3(kc)+ct3(kc+1))
   ct4_sg(kc) = 0.5_wp*(ct4(kc)+ct4(kc+1))
   adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
   abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))
end do

endif

do kc=0, grd%KCMAX-1
   temp_c_help(kc) = 0.5_wp*(st%temp_c(kc,j,i)+st%temp_c(kc+1,j,i))
   ce6(kc) = grd%at6(kc) &
             *kappa_val(temp_c_help(kc))/(c_val(temp_c_help(kc))*st%H_c(j,i))
   ci2(kc) = grd%ai2(kc)/st%H_c(j,i)
end do

if (par%adv_hor==3) then
dtt_dxi  = 2.0_wp*grd%dtt_2dxi
dtt_deta = 2.0_wp*grd%dtt_2deta
endif

!-------- Set up the equations for the bedrock temperature --------

kr=0
lgs_a1(kr) = 1.0_wp
lgs_a2(kr) = -1.0_wp
lgs_b(kr)    = clb1

if (par%q_litho==1) then
!   (coupled heat-conducting bedrock)

do kr=1, grd%KRMAX-1
   lgs_a0(kr) = - ctr1
   lgs_a1(kr) = 1.0_wp + 2.0_wp*ctr1
   lgs_a2(kr) = - ctr1
   lgs_b(kr)    = st%temp_r(kr,j,i)
end do

else if (par%q_litho==0) then
!   (no coupled heat-conducting bedrock)

do kr=1, grd%KRMAX-1
   lgs_a0(kr) = 1.0_wp
   lgs_a1(kr) = 0.0_wp
   lgs_a2(kr) = -1.0_wp
   lgs_b(kr)  = 2.0_wp*clb1
end do

endif

kr=grd%KRMAX
lgs_a0(kr) = 0.0_wp
lgs_a1(kr) = 1.0_wp
lgs_b(kr)  = st%temp_c_m(0,j,i)-DELTA_TM_SW

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KRMAX)

!-------- Assign the result --------

do kr=0, grd%KRMAX
   st%temp_r_neu(kr,j,i) = lgs_x(kr)
end do

!-------- Set up the equations for the ice temperature --------

kc=0
lgs_a1(kc) = 1.0_wp
lgs_a2(kc) = 0.0_wp
lgs_b(kc)  = enth_fct_temp_omega(st%temp_c_m(kc,j,i)-DELTA_TM_SW, 0.0_wp)
                                                  ! zero water content assumed

do kc=1, grd%KCMAX-1

if (par%adv_vert==1) then

   lgs_a0(kc) = -0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) = 1.0_wp+ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) = 0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ce5(kc)*ce6(kc)

else if (par%adv_vert==2) then

   lgs_a0(kc) &
         = -0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) &
         = 1.0_wp &
           +0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           +ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) &
         =  0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           -ce5(kc)*ce6(kc)

else if (par%adv_vert==3) then

   adv_vert_help = 0.5_wp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) &
         = -max(adv_vert_help, 0.0_wp) &
           -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) &
         = 1.0_wp &
           +max(adv_vert_help, 0.0_wp)-min(adv_vert_help, 0.0_wp) &
           +ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) &
         =  min(adv_vert_help, 0.0_wp) &
           -ce5(kc)*ce6(kc)

endif

if (par%adv_hor==2) then

   lgs_b(kc) = st%enth_c(kc,j,i) + ce7(kc) &
       -grd%dtt_2dxi* &
          ( (st%vx_c(kc,j,i)-abs(st%vx_c(kc,j,i))) &
            *(st%enth_c(kc,j,i+1)-st%enth_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_c(kc,j,i-1)+abs(st%vx_c(kc,j,i-1))) &
            *(st%enth_c(kc,j,i)-st%enth_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_c(kc,j,i)-abs(st%vy_c(kc,j,i))) &
            *(st%enth_c(kc,j+1,i)-st%enth_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_c(kc,j-1,i)+abs(st%vy_c(kc,j-1,i))) &
            *(st%enth_c(kc,j,i)-st%enth_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

   vx_c_help = 0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1))
   vy_c_help = 0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i))

   lgs_b(kc) = st%enth_c(kc,j,i) + ce7(kc) &
       -dtt_dxi* &
          ( min(vx_c_help, 0.0_wp) &
            *(st%enth_c(kc,j,i+1)-st%enth_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_c_help, 0.0_wp) &
            *(st%enth_c(kc,j,i)-st%enth_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_c_help, 0.0_wp) &
            *(st%enth_c(kc,j+1,i)-st%enth_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_c_help, 0.0_wp) &
            *(st%enth_c(kc,j,i)-st%enth_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

end do

kc=grd%KCMAX
lgs_a0(kc) = 0.0_wp
lgs_a1(kc) = 1.0_wp
lgs_b(kc)  = enth_fct_temp_omega(st%temp_s(j,i), 0.0_wp)
             ! zero water content assumed

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KCMAX)

!-------- Assign the result --------

do kc=0, grd%KCMAX
   st%enth_c_neu(kc,j,i)  = lgs_x(kc)
   st%temp_c_neu(kc,j,i)  = temp_fct_enth(st%enth_c_neu(kc,j,i), st%temp_c_m(kc,j,i))
   st%omega_c_neu(kc,j,i) = omega_fct_enth(st%enth_c_neu(kc,j,i), st%temp_c_m(kc,j,i))
end do

!-------- Set enthalpy and water content in the redundant,
!         lower (kt) ice layer to the value at the ice base --------

do kt=0, grd%KTMAX
   st%enth_t_neu(kt,j,i)  = st%enth_c_neu(0,j,i)
   st%omega_t_neu(kt,j,i) = st%omega_c_neu(0,j,i)
end do

!-------- Water drainage from the non-existing temperate ice --------

st%Q_tld(j,i) = 0.0_wp

!-------- Set up the equations for the age of cold ice --------

kc=0                                                 ! adv_vert_sg(0) <= 0
lgs_a1(kc) = 1.0_wp - min(adv_vert_sg(kc), 0.0_wp)   ! (directed downward)
lgs_a2(kc) = min(adv_vert_sg(kc), 0.0_wp)            ! assumed/enforced

if (par%adv_hor==2) then

lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -grd%dtt_2dxi* &
          ( (st%vx_c(kc,j,i)-abs(st%vx_c(kc,j,i))) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_c(kc,j,i-1)+abs(st%vx_c(kc,j,i-1))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_c(kc,j,i)-abs(st%vy_c(kc,j,i))) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_c(kc,j-1,i)+abs(st%vy_c(kc,j-1,i))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

vx_c_help = 0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1))
vy_c_help = 0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i))

lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -dtt_dxi* &
          ( min(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

do kc=1, grd%KCMAX-1

if (par%adv_vert==1) then

   lgs_a0(kc) = -0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ci1(kc)*ci2(kc-1)
   lgs_a1(kc) = 1.0_wp+ci1(kc)*(ci2(kc)+ci2(kc-1))
   lgs_a2(kc) = 0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ci1(kc)*ci2(kc)

else if (par%adv_vert==2) then

   lgs_a0(kc) = -0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1))
   lgs_a1(kc) =  1.0_wp &
                +0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
                -0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  )
   lgs_a2(kc) =  0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  )

else if (par%adv_vert==3) then

   adv_vert_help = 0.5_wp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) = -max(adv_vert_help, 0.0_wp)
   lgs_a1(kc) =  1.0_wp &
                +max(adv_vert_help, 0.0_wp)-min(adv_vert_help, 0.0_wp)
   lgs_a2(kc) =  min(adv_vert_help, 0.0_wp)

endif

if (par%adv_hor==2) then

   lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -grd%dtt_2dxi* &
          ( (st%vx_c(kc,j,i)-abs(st%vx_c(kc,j,i))) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_c(kc,j,i-1)+abs(st%vx_c(kc,j,i-1))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_c(kc,j,i)-abs(st%vy_c(kc,j,i))) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_c(kc,j-1,i)+abs(st%vy_c(kc,j-1,i))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

   vx_c_help = 0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1))
   vy_c_help = 0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i))

   lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -dtt_dxi* &
          ( min(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

end do

kc=grd%KCMAX
if (st%as_perp(j,i) >= zero) then
   lgs_a0(kc) = 0.0_wp
   lgs_a1(kc) = 1.0_wp
   lgs_b(kc)  = 0.0_wp
else
   lgs_a0(kc) = -max(adv_vert_sg(kc-1), 0.0_wp)
   lgs_a1(kc) = 1.0_wp + max(adv_vert_sg(kc-1), 0.0_wp)
                       ! adv_vert_sg(grd%KCMAX-1) >= 0 (directed upward)
                       ! assumed/enforced
if (par%adv_hor==2) then

   lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -grd%dtt_2dxi* &
          ( (st%vx_c(kc,j,i)-abs(st%vx_c(kc,j,i))) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_c(kc,j,i-1)+abs(st%vx_c(kc,j,i-1))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_c(kc,j,i)-abs(st%vy_c(kc,j,i))) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_c(kc,j-1,i)+abs(st%vy_c(kc,j-1,i))) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

   vx_c_help = 0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1))
   vy_c_help = 0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i))

   lgs_b(kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
       -dtt_dxi* &
          ( min(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i+1)-st%age_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j+1,i)-st%age_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_c_help, 0.0_wp) &
            *(st%age_c(kc,j,i)-st%age_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

end if

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KCMAX)

!-------- Assign the result,
!         restriction to interval [0, AGE_MAX yr] --------

do kc=0, grd%KCMAX

   st%age_c_neu(kc,j,i) = lgs_x(kc)

   if (st%age_c_neu(kc,j,i) < (par%age_min*sec_year)) &
                           st%age_c_neu(kc,j,i) = 0.0_wp
   if (st%age_c_neu(kc,j,i) > (par%age_max*sec_year)) &
                           st%age_c_neu(kc,j,i) = par%age_max*sec_year

end do

!-------- Age of the ice in the redundant, lower (kt) ice layer --------

do kt=0, grd%KTMAX
   st%age_t_neu(kt,j,i) = st%age_c_neu(0,j,i)
end do

end subroutine calc_temp_enth_ssa

!-------------------------------------------------------------------------------

end module calc_temp_enth_m
!
