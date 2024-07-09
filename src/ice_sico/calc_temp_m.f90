!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ t e m p _ m
!
!> @file
!!
!! Computation of temperature, water content and age.
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
!> Computation of temperature, water content and age.
!<------------------------------------------------------------------------------
module calc_temp_m

  !$ use omp_lib
  use timer, only : sec_year

  use sico_types_m
  use sico_state
  use sico_grid_mod
  use sico_timer
  use sico_params, only : sico_par_class, delta_tm_sw, omega_max
  use ice_material_properties_m, only : ratefac_c, ratefac_t, kappa_val, c_val, &
                                       creep, viscosity
  use sico_maths_m,              only : tri_sle



  implicit none

  private
  public :: calc_temp_poly, calc_temp_cold, calc_temp_const

contains

!-------------------------------------------------------------------------------
!> Computation of temperature, water content and age in polythermal mode.
!<------------------------------------------------------------------------------
subroutine calc_temp_poly(st,grd,tmr,par)

implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class),  intent(in)    :: grd
  type(sico_timer_class), intent(in)    :: tmr
  type(sico_par_class), intent(in)    :: par

integer :: i, j, kc, kt, kr, ii, jj
real(wp) :: time_lag_cts
real(wp) :: Vol_t, Vol_t_smooth, korrfakt_t

!-------- Computation loop for temperature, water content and age --------

!$omp parallel do &
!$omp private (i,j,kc) 
do i=1, grd%IMAX-1   ! skipping domain margins
do j=1, grd%JMAX-1   ! skipping domain margins
      !!$ print *,i,j,omp_get_thread_num(),'/',omp_get_num_threads()

   if (st%maske(j,i)==0) then   ! glaciated land

!  ------ Old vertical column cold

      if (st%n_cts(j,i).eq.-1) then

         st%n_cts_neu(j,i) = st%n_cts(j,i)
         st%zm_neu(j,i)  = st%zb(j,i)
         st%H_c_neu(j,i) = st%H_c(j,i)
         st%H_t_neu(j,i) = st%H_t(j,i)

         call calc_temp1(st,grd,tmr,par,i,j) 

!    ---- Check whether base has become temperate

         if (st%temp_c_neu(0,j,i).gt.st%temp_c_m(0,j,i)) then

            st%n_cts_neu(j,i) = 0

            call calc_temp2(st,grd,tmr,par,i,j)

         end if

!    ---- Check whether even temperate layer has formed

         if ( &
              ( st%n_cts_neu(j,i).eq.0 ).and. &
              ( (st%temp_c_neu(1,j,i)-st%temp_c_neu(0,j,i)) &
                .gt.(grd%am1*st%H_c_neu(j,i)) ) &
            ) then

            st%n_cts_neu(j,i) = 1
            st%zm_neu(j,i)  = st%zb(j,i)+0.001_wp
            st%H_c_neu(j,i) = st%H_c(j,i)-0.001_wp
            st%H_t_neu(j,i) = st%H_t(j,i)+0.001_wp
!                 ! CTS preliminarily positioned 1 mm above ice base --------

            call calc_temp3(st,grd,tmr,par,i,j)

            call shift_cts_upward(st,grd,tmr,par,i,j)

         end if

!  ------ Old vertical column with temperate base

      else if (st%n_cts(j,i).eq.0) then

         st%n_cts_neu(j,i) = st%n_cts(j,i)
         st%zm_neu(j,i)  = st%zb(j,i)
         st%H_c_neu(j,i) = st%H_c(j,i)
         st%H_t_neu(j,i) = st%H_t(j,i)

         call calc_temp2(st,grd,tmr,par,i,j)

!    ---- Check whether temperate base becomes cold

         if ( (st%temp_c_neu(1,j,i)-st%temp_c_neu(0,j,i)) &
               .lt. (grd%am1*st%H_c(j,i)) ) then

            st%n_cts_neu(j,i) = -1

            call calc_temp1(st,grd,tmr,par,i,j)

            if (st%temp_c_neu(0,j,i).ge.st%temp_c_m(0,j,i)) then

               st%n_cts_neu(j,i) = 0

               call calc_temp2(st,grd,tmr,par,i,j)

            end if

         end if

!    ---- Check whether temperate layer has formed

         if ( &
              ( st%n_cts_neu(j,i).eq.0 ).and. &
              ( (st%temp_c_neu(1,j,i)-st%temp_c_neu(0,j,i)) &
                .gt.(grd%am1*st%H_c_neu(j,i)) ) &
            ) then

            st%n_cts_neu(j,i) = 1
            st%zm_neu(j,i)  = st%zb(j,i)+0.001_wp
            st%H_c_neu(j,i) = st%H_c(j,i)-0.001_wp
            st%H_t_neu(j,i) = st%H_t(j,i)+0.001_wp
!                 ! CTS preliminarily positioned 1 mm above ice base --------

            call calc_temp3(st,grd,tmr,par,i,j)

            call shift_cts_upward(st,grd,tmr,par,i,j)

         end if

!  ------ Old vertical column with temperate base and temperate layer

      else   ! n_cts(j,i).eq.1

         st%n_cts_neu(j,i) = st%n_cts(j,i)
         st%zm_neu(j,i)  = st%zm(j,i)
         st%H_c_neu(j,i) = st%H_c(j,i)
         st%H_t_neu(j,i) = st%H_t(j,i)

         call calc_temp3(st,grd,tmr,par,i,j)

         if ( (st%temp_c_neu(0,j,i)-(-BETA*st%H_c_neu(j,i))).gt.0.0_wp ) &
         then
            call shift_cts_upward(st,grd,tmr,par,i,j)
         else
            call shift_cts_downward(st,grd,tmr,par,i,j)
         end if

      end if

   else if (st%maske(j,i)==3 .and. par%margin==3) then   ! floating ice

      st%n_cts_neu(j,i) = -1
      st%zm_neu(j,i)  = st%zb(j,i)
      st%H_c_neu(j,i) = st%H_c(j,i) + st%H_t(j,i)
      st%H_t_neu(j,i) = 0.0_wp

      call calc_temp_ssa(st,grd,tmr,par,i,j)

!  ------ Reset temperatures above melting to the melting point
!         (should not occur, but just in case)

      do kc=0, grd%KCMAX
         if (st%temp_c_neu(kc,j,i) > st%temp_c_m(kc,j,i)) &
                    st%temp_c_neu(kc,j,i) = st%temp_c_m(kc,j,i)
      end do

   else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

      st%n_cts_neu(j,i) = -1
      st%zm_neu(j,i)  = st%zb(j,i)
      st%H_c_neu(j,i) = st%H_c(j,i)
      st%H_t_neu(j,i) = st%H_t(j,i)

      call calc_temp_r(st,grd,par,tmr,i,j)

endif

end do
end do   ! End of computation loop
!$omp end parallel do

!-------- Extrapolate values on margins --------

!  ------ Lower left corner

i=0
j=0

if ( (st%maske(j,i).eq.0).or.(st%maske(j,i).eq.3) ) then
                              ! glaciated land or floating ice
   ii=i+1
   jj=j+1

   do kc=0,grd%KCMAX
      st%temp_c_neu(kc,j,i) = st%temp_c_neu(kc,jj,ii)   ! set cold-ice temperature
      st%age_c_neu(kc,j,i)  = st%age_c_neu(kc,jj,ii)    ! set cold-ice age
   end do

   do kt=0,grd%KTMAX
      st%omega_t_neu(kt,j,i) = st%omega_t_neu(kt,jj,ii) ! set temp.-ice water content
      st%age_t_neu(kt,j,i)   = st%age_t_neu(kt,jj,ii)   ! set temp.-ice age
   end do

   do kr=0,grd%KRMAX
      st%temp_r_neu(kr,j,i) = st%temp_r_neu(kr,jj,ii)   ! set bedrock temperature
   end do

   st%n_cts_neu(j,i) = min(st%n_cts_neu(jj,ii),0)   ! temperate layer excluded
   st%H_c_neu(j,i)   = st%H_c(j,i)
   st%H_t_neu(j,i)   = st%H_t(j,i)
   st%zm_neu(j,i)  = st%zb(j,i)

else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

   st%n_cts_neu(j,i) = -1
   st%zm_neu(j,i)    = st%zb(j,i)
   st%H_c_neu(j,i)   = st%H_c(j,i)
   st%H_t_neu(j,i)   = st%H_t(j,i)

   call calc_temp_r(st,grd,par,tmr,i,j)

end if

!  ------ Lower right corner

i=grd%IMAX
j=0

if ( (st%maske(j,i).eq.0).or.(st%maske(j,i).eq.3) ) then
                              ! glaciated land or floating ice
   ii=i-1
   jj=j+1

   do kc=0,grd%KCMAX
      st%temp_c_neu(kc,j,i) = st%temp_c_neu(kc,jj,ii)   ! set cold-ice temperature
      st%age_c_neu(kc,j,i)  = st%age_c_neu(kc,jj,ii)    ! set cold-ice age
   end do

   do kt=0,grd%KTMAX
      st%omega_t_neu(kt,j,i) = st%omega_t_neu(kt,jj,ii) ! set temp.-ice water content
      st%age_t_neu(kt,j,i)   = st%age_t_neu(kt,jj,ii)   ! set temp.-ice age
   end do

   do kr=0,grd%KRMAX
      st%temp_r_neu(kr,j,i) = st%temp_r_neu(kr,jj,ii)   ! set bedrock temperature
   end do

   st%n_cts_neu(j,i) = min(st%n_cts_neu(jj,ii),0)   ! temperate layer excluded
   st%H_c_neu(j,i)   = st%H_c(j,i)
   st%H_t_neu(j,i)   = st%H_t(j,i)
   st%zm_neu(j,i)  = st%zb(j,i)

else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

   st%n_cts_neu(j,i) = -1
   st%zm_neu(j,i)    = st%zb(j,i)
   st%H_c_neu(j,i)   = st%H_c(j,i)
   st%H_t_neu(j,i)   = st%H_t(j,i)

   call calc_temp_r(st,grd,par,tmr,i,j)

end if

!  ------ Upper left corner

i=0
j=grd%JMAX

if ( (st%maske(j,i).eq.0).or.(st%maske(j,i).eq.3) ) then
                              ! glaciated land or floating ice
   ii=i+1
   jj=j-1

   do kc=0,grd%KCMAX
      st%temp_c_neu(kc,j,i) = st%temp_c_neu(kc,jj,ii)   ! set cold-ice temperature
      st%age_c_neu(kc,j,i)  = st%age_c_neu(kc,jj,ii)    ! set cold-ice age
   end do

   do kt=0,grd%KTMAX
      st%omega_t_neu(kt,j,i) = st%omega_t_neu(kt,jj,ii) ! set temp.-ice water content
      st%age_t_neu(kt,j,i)   = st%age_t_neu(kt,jj,ii)   ! set temp.-ice age
   end do

   do kr=0,grd%KRMAX
      st%temp_r_neu(kr,j,i) = st%temp_r_neu(kr,jj,ii)   ! set bedrock temperature
   end do

   st%n_cts_neu(j,i) = min(st%n_cts_neu(jj,ii),0)   ! temperate layer excluded
   st%H_c_neu(j,i)   = st%H_c(j,i)
   st%H_t_neu(j,i)   = st%H_t(j,i)
   st%zm_neu(j,i)  = st%zb(j,i)

else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

   st%n_cts_neu(j,i) = -1
   st%zm_neu(j,i)    = st%zb(j,i)
   st%H_c_neu(j,i)   = st%H_c(j,i)
   st%H_t_neu(j,i)   = st%H_t(j,i)

   call calc_temp_r(st,grd,par,tmr,i,j)

end if

!  ------ Upper right corner

i=grd%IMAX
j=grd%JMAX

if ( (st%maske(j,i).eq.0).or.(st%maske(j,i).eq.3) ) then
                              ! glaciated land or floating ice
   ii=i-1
   jj=j-1

   do kc=0,grd%KCMAX
      st%temp_c_neu(kc,j,i) = st%temp_c_neu(kc,jj,ii)   ! set cold-ice temperature
      st%age_c_neu(kc,j,i)  = st%age_c_neu(kc,jj,ii)    ! set cold-ice age
   end do

   do kt=0,grd%KTMAX
      st%omega_t_neu(kt,j,i) = st%omega_t_neu(kt,jj,ii) ! set temp.-ice water content
      st%age_t_neu(kt,j,i)   = st%age_t_neu(kt,jj,ii)   ! set temp.-ice age
   end do

   do kr=0,grd%KRMAX
      st%temp_r_neu(kr,j,i) = st%temp_r_neu(kr,jj,ii)   ! set bedrock temperature
   end do

   st%n_cts_neu(j,i) = min(st%n_cts_neu(jj,ii),0)   ! temperate layer excluded
   st%H_c_neu(j,i)   = st%H_c(j,i)
   st%H_t_neu(j,i)   = st%H_t(j,i)
   st%zm_neu(j,i)  = st%zb(j,i)

else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

   st%n_cts_neu(j,i) = -1
   st%zm_neu(j,i)    = st%zb(j,i)
   st%H_c_neu(j,i)   = st%H_c(j,i)
   st%H_t_neu(j,i)   = st%H_t(j,i)

   call calc_temp_r(st,grd,par,tmr,i,j)

end if

!  ------ Lower and upper margins

do i=1, grd%IMAX-1

!    ---- Lower margin

   j=0

   if ( (st%maske(j,i).eq.0).or.(st%maske(j,i).eq.3) ) then
                                 ! glaciated land or floating ice
      ii=i
      jj=j+1

      do kc=0,grd%KCMAX
         st%temp_c_neu(kc,j,i) = st%temp_c_neu(kc,jj,ii)   ! set cold-ice temperature
         st%age_c_neu(kc,j,i)  = st%age_c_neu(kc,jj,ii)    ! set cold-ice age
      end do

      do kt=0,grd%KTMAX
         st%omega_t_neu(kt,j,i) = st%omega_t_neu(kt,jj,ii) ! set temp.-ice water content
         st%age_t_neu(kt,j,i)   = st%age_t_neu(kt,jj,ii)   ! set temp.-ice age
      end do

      do kr=0,grd%KRMAX
         st%temp_r_neu(kr,j,i) = st%temp_r_neu(kr,jj,ii)   ! set bedrock temperature
      end do

      st%n_cts_neu(j,i) = min(st%n_cts_neu(jj,ii),0)   ! temperate layer excluded
      st%H_c_neu(j,i)   = st%H_c(j,i)
      st%H_t_neu(j,i)   = st%H_t(j,i)
      st%zm_neu(j,i)  = st%zb(j,i)

   else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

      st%n_cts_neu(j,i) = -1
      st%zm_neu(j,i)    = st%zb(j,i)
      st%H_c_neu(j,i)   = st%H_c(j,i)
      st%H_t_neu(j,i)   = st%H_t(j,i)

      call calc_temp_r(st,grd,par,tmr,i,j)

   end if

!    ---- Upper margin

   j=grd%JMAX

   if ( (st%maske(j,i).eq.0).or.(st%maske(j,i).eq.3) ) then
                                 ! glaciated land or floating ice
      ii=i
      jj=j-1

      do kc=0,grd%KCMAX
         st%temp_c_neu(kc,j,i) = st%temp_c_neu(kc,jj,ii)   ! set cold-ice temperature
         st%age_c_neu(kc,j,i)  = st%age_c_neu(kc,jj,ii)    ! set cold-ice age
      end do

      do kt=0,grd%KTMAX
         st%omega_t_neu(kt,j,i) = st%omega_t_neu(kt,jj,ii) ! set temp.-ice water content
         st%age_t_neu(kt,j,i)   = st%age_t_neu(kt,jj,ii)   ! set temp.-ice age
      end do

      do kr=0,grd%KRMAX
         st%temp_r_neu(kr,j,i) = st%temp_r_neu(kr,jj,ii)   ! set bedrock temperature
      end do

      st%n_cts_neu(j,i) = min(st%n_cts_neu(jj,ii),0)   ! temperate layer excluded
      st%H_c_neu(j,i)   = st%H_c(j,i)
      st%H_t_neu(j,i)   = st%H_t(j,i)
      st%zm_neu(j,i)  = st%zb(j,i)

   else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

      st%n_cts_neu(j,i) = -1
      st%zm_neu(j,i)    = st%zb(j,i)
      st%H_c_neu(j,i)   = st%H_c(j,i)
      st%H_t_neu(j,i)   = st%H_t(j,i)

      call calc_temp_r(st,grd,par,tmr,i,j)

   end if

end do

!  ------ Left and right margins

do j=1, grd%JMAX-1

!    ---- Left margin

   i=0

   if ( (st%maske(j,i).eq.0).or.(st%maske(j,i).eq.3) ) then
                                 ! glaciated land or floating ice
      ii=i+1
      jj=j

      do kc=0,grd%KCMAX
         st%temp_c_neu(kc,j,i) = st%temp_c_neu(kc,jj,ii)   ! set cold-ice temperature
         st%age_c_neu(kc,j,i)  = st%age_c_neu(kc,jj,ii)    ! set cold-ice age
      end do

      do kt=0,grd%KTMAX
         st%omega_t_neu(kt,j,i) = st%omega_t_neu(kt,jj,ii) ! set temp.-ice water content
         st%age_t_neu(kt,j,i)   = st%age_t_neu(kt,jj,ii)   ! set temp.-ice age
      end do

      do kr=0,grd%KRMAX
         st%temp_r_neu(kr,j,i) = st%temp_r_neu(kr,jj,ii)   ! set bedrock temperature
      end do

      st%n_cts_neu(j,i) = min(st%n_cts_neu(jj,ii),0)   ! temperate layer excluded
      st%H_c_neu(j,i)   = st%H_c(j,i)
      st%H_t_neu(j,i)   = st%H_t(j,i)
      st%zm_neu(j,i)  = st%zb(j,i)

   else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

      st%n_cts_neu(j,i) = -1
      st%zm_neu(j,i)    = st%zb(j,i)
      st%H_c_neu(j,i)   = st%H_c(j,i)
      st%H_t_neu(j,i)   = st%H_t(j,i)

      call calc_temp_r(st,grd,par,tmr,i,j)

   end if

!    ---- Right margin

   i=grd%IMAX

   if ( (st%maske(j,i).eq.0).or.(st%maske(j,i).eq.3) ) then
                                 ! glaciated land or floating ice
      ii=i-1
      jj=j

      do kc=0,grd%KCMAX
         st%temp_c_neu(kc,j,i) = st%temp_c_neu(kc,jj,ii)   ! set cold-ice temperature
         st%age_c_neu(kc,j,i)  = st%age_c_neu(kc,jj,ii)    ! set cold-ice age
      end do

      do kt=0,grd%KTMAX
         st%omega_t_neu(kt,j,i) = st%omega_t_neu(kt,jj,ii) ! set temp.-ice water content
         st%age_t_neu(kt,j,i)   = st%age_t_neu(kt,jj,ii)   ! set temp.-ice age
      end do

      do kr=0,grd%KRMAX
         st%temp_r_neu(kr,j,i) = st%temp_r_neu(kr,jj,ii)   ! set bedrock temperature
      end do

      st%n_cts_neu(j,i) = min(st%n_cts_neu(jj,ii),0)   ! temperate layer excluded
      st%H_c_neu(j,i)   = st%H_c(j,i)
      st%H_t_neu(j,i)   = st%H_t(j,i)
      st%zm_neu(j,i)  = st%zb(j,i)

   else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

      st%n_cts_neu(j,i) = -1
      st%zm_neu(j,i)    = st%zb(j,i)
      st%H_c_neu(j,i)   = st%H_c(j,i)
      st%H_t_neu(j,i)   = st%H_t(j,i)

      call calc_temp_r(st,grd,par,tmr,i,j)

   end if

end do

!-------- Dummy values for st%omega_c_neu and kc_cts_neu --------

st%omega_c_neu = 0.0_wp   ! not needed for
st%kc_cts_neu  = 0        ! the polythermal mode

!-------- Smoothing of H_t_neu with numerical diffusion --------

!  ------ Volume of temperate ice without smoothing

Vol_t = 0.0_wp
do i=0, grd%IMAX   ! extended to domain margins (22.1.02 -> V1.1)
do j=0, grd%JMAX   ! extended to domain margins (22.1.02 -> V1.1)
   if (st%n_cts_neu(j,i).eq.1) then
      Vol_t = Vol_t + st%H_t_neu(j,i)*grd%area(j,i)
   end if
end do
end do

!  ------ Smoothing

do i=1, grd%IMAX-1
do j=1, grd%JMAX-1
   if (st%n_cts_neu(j,i).ne.-1) then

      st%dH_t_smooth(j,i) = par%numdiff_h_t* ( -4.0_wp*st%H_t_neu(j,i) &
                             +st%H_t_neu(j,i+1)+st%H_t_neu(j,i-1) &
                             +st%H_t_neu(j+1,i)+st%H_t_neu(j-1,i) )
      if (st%dH_t_smooth(j,i).gt.0.001_wp) st%n_cts_neu(j,i) = 1

   end if
end do
end do

do i=1, grd%IMAX-1
do j=1, grd%JMAX-1
   if (st%n_cts_neu(j,i).eq.1) then
      st%H_t_neu(j,i) = st%H_t_neu(j,i) + st%dH_t_smooth(j,i)
   end if
end do
end do

!  ------ Volume of temperate ice with smoothing

Vol_t_smooth = 0.0_wp
do i=0, grd%IMAX   ! extended to domain margins (22.1.02 -> V1.1)
do j=0, grd%JMAX   ! extended to domain margins (22.1.02 -> V1.1)
   if (st%n_cts_neu(j,i).eq.1) then
      Vol_t_smooth = Vol_t_smooth + st%H_t_neu(j,i)*grd%area(j,i)
   end if
end do
end do

!  ------ Correction so that volume is not changed by the smoothing

if (Vol_t_smooth.gt.0.0_wp) then

   korrfakt_t = Vol_t/Vol_t_smooth
   do i=0, grd%IMAX   ! extended to domain margins (22.1.02 -> V1.1)
   do j=0, grd%JMAX   ! extended to domain margins (22.1.02 -> V1.1)
      if (st%n_cts_neu(j,i).eq.1) then
         st%H_t_neu(j,i) = st%H_t_neu(j,i)*korrfakt_t
!               st%zm_neu(j,i)  = st%zb(j,i) + st%H_t_neu(j,i)
!               st%H_c_neu(j,i) = st%zs(j,i) - st%zm_neu(j,i)
      end if
   end do
   end do

end if

!-------- Numerical time lag for evolution of H_t_neu --------

time_lag_cts = par%tau_cts*sec_year   ! yr --> s

do i=0, grd%IMAX   ! extended to domain margins (22.1.02 -> V1.1)
do j=0, grd%JMAX   ! extended to domain margins (22.1.02 -> V1.1)

   if (st%n_cts_neu(j,i).eq.1) then

      st%H_t_neu(j,i) = ( time_lag_cts*st%H_t(j,i) &
                       + tmr%dtime_temp*st%H_t_neu(j,i) ) &
                     /(time_lag_cts+tmr%dtime_temp)

      st%zm_neu(j,i)  = st%zb(j,i) + st%H_t_neu(j,i)
      st%H_c_neu(j,i) = st%zs(j,i) - st%zm_neu(j,i)

   end if

end do
end do

end subroutine calc_temp_poly

!-------------------------------------------------------------------------------
!> Computation of temperature and age in cold-ice mode.
!<------------------------------------------------------------------------------
subroutine calc_temp_cold(st,grd,tmr,par)

implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class),  intent(in)    :: grd
  type(sico_timer_class), intent(in)    :: tmr
  type(sico_par_class), intent(in)    :: par

integer :: i, j, kc, kr, ii, jj

!-------- Computation loop for temperature and age --------

!$omp parallel do &
!$omp private (i,j,kc) 
do i=1, grd%IMAX-1   ! skipping domain margins
do j=1, grd%JMAX-1   ! skipping domain margins

   if (st%maske(j,i)==0) then   ! glaciated land

      st%n_cts_neu(j,i) = -1
      st%zm_neu(j,i)  = st%zb(j,i)
      st%H_c_neu(j,i) = st%H_c(j,i)
      st%H_t_neu(j,i) = st%H_t(j,i)

      call calc_temp1(st,grd,tmr,par,i,j)

!  ------ Reset temperatures above melting to the melting point,
!         look for the CTS

      st%kc_cts_neu(j,i) = 0

      if (st%temp_c_neu(0,j,i).gt.st%temp_c_m(0,j,i)) then
         st%n_cts_neu(j,i)        = 0
         st%kc_cts_neu(j,i)       = 0
         st%temp_c_neu(0,j,i)     = st%temp_c_m(0,j,i)
         st%temp_r_neu(grd%KRMAX,j,i) = st%temp_c_m(0,j,i)
      end if

      do kc=1, grd%KCMAX
         if (st%temp_c_neu(kc,j,i).gt.st%temp_c_m(kc,j,i)) then
            st%kc_cts_neu(j,i)    = kc
            st%temp_c_neu(kc,j,i) = st%temp_c_m(kc,j,i)
         end if
      end do

   else if (st%maske(j,i)==3 .and. par%margin==3) then   ! floating ice

      st%n_cts_neu(j,i)  = -1
      st%kc_cts_neu(j,i) = 0
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = 0.0_wp

      call calc_temp_ssa(st,grd,tmr,par,i,j)

!  ------ Reset temperatures above melting to the melting point
!         (should not occur, but just in case)

      do kc=0, grd%KCMAX
         if (st%temp_c_neu(kc,j,i) > st%temp_c_m(kc,j,i)) &
                    st%temp_c_neu(kc,j,i) = st%temp_c_m(kc,j,i)
      end do

   else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

      st%n_cts_neu(j,i)  = -1
      st%kc_cts_neu(j,i) = 0
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = st%H_t(j,i)

      call calc_temp_r(st,grd,par,tmr,i,j)

   end if

end do
end do 
!$omp end parallel do

!-------- Extrapolate values on margins --------

!  ------ Lower left corner

i=0
j=0

if ( (st%maske(j,i).eq.0).or.(st%maske(j,i).eq.3) ) then
                              ! glaciated land or floating ice
   ii=i+1
   jj=j+1

   do kc=0,grd%KCMAX
      st%temp_c_neu(kc,j,i) = st%temp_c_neu(kc,jj,ii)   ! set cold-ice temperature
      st%age_c_neu(kc,j,i)  = st%age_c_neu(kc,jj,ii)    ! set cold-ice age
   end do

   do kr=0,grd%KRMAX
      st%temp_r_neu(kr,j,i) = st%temp_r_neu(kr,jj,ii)   ! set bedrock temperature
   end do

   st%n_cts_neu(j,i)  = st%n_cts_neu(jj,ii)
   st%kc_cts_neu(j,i) = st%kc_cts_neu(jj,ii)
   st%zm_neu(j,i)     = st%zb(j,i)
   st%H_c_neu(j,i)    = st%H_c(j,i)
   st%H_t_neu(j,i)    = st%H_t(j,i)

else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

   st%n_cts_neu(j,i)  = -1
   st%kc_cts_neu(j,i) = 0
   st%zm_neu(j,i)     = st%zb(j,i)
   st%H_c_neu(j,i)    = st%H_c(j,i)
   st%H_t_neu(j,i)    = st%H_t(j,i)

   call calc_temp_r(st,grd,par,tmr,i,j)

end if

!  ------ Lower right corner

i=grd%IMAX
j=0

if ( (st%maske(j,i).eq.0).or.(st%maske(j,i).eq.3) ) then
                              ! glaciated land or floating ice
   ii=i-1
   jj=j+1

   do kc=0,grd%KCMAX
      st%temp_c_neu(kc,j,i) = st%temp_c_neu(kc,jj,ii)   ! set cold-ice temperature
      st%age_c_neu(kc,j,i)  = st%age_c_neu(kc,jj,ii)    ! set cold-ice age
   end do

   do kr=0,grd%KRMAX
      st%temp_r_neu(kr,j,i) = st%temp_r_neu(kr,jj,ii)   ! set bedrock temperature
   end do

   st%n_cts_neu(j,i)  = st%n_cts_neu(jj,ii)
   st%kc_cts_neu(j,i) = st%kc_cts_neu(jj,ii)
   st%zm_neu(j,i)     = st%zb(j,i)
   st%H_c_neu(j,i)    = st%H_c(j,i)
   st%H_t_neu(j,i)    = st%H_t(j,i)

else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

   st%n_cts_neu(j,i)  = -1
   st%kc_cts_neu(j,i) = 0
   st%zm_neu(j,i)     = st%zb(j,i)
   st%H_c_neu(j,i)    = st%H_c(j,i)
   st%H_t_neu(j,i)    = st%H_t(j,i)

   call calc_temp_r(st,grd,par,tmr,i,j)

end if

!  ------ Upper left corner

i=0
j=grd%JMAX

if ( (st%maske(j,i).eq.0).or.(st%maske(j,i).eq.3) ) then
                              ! glaciated land or floating ice
   ii=i+1
   jj=j-1

   do kc=0,grd%KCMAX
      st%temp_c_neu(kc,j,i) = st%temp_c_neu(kc,jj,ii)   ! set cold-ice temperature
      st%age_c_neu(kc,j,i)  = st%age_c_neu(kc,jj,ii)    ! set cold-ice age
   end do

   do kr=0,grd%KRMAX
      st%temp_r_neu(kr,j,i) = st%temp_r_neu(kr,jj,ii)   ! set bedrock temperature
   end do

   st%n_cts_neu(j,i)  = st%n_cts_neu(jj,ii)
   st%kc_cts_neu(j,i) = st%kc_cts_neu(jj,ii)
   st%zm_neu(j,i)     = st%zb(j,i)
   st%H_c_neu(j,i)    = st%H_c(j,i)
   st%H_t_neu(j,i)    = st%H_t(j,i)

else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

   st%n_cts_neu(j,i)  = -1
   st%kc_cts_neu(j,i) = 0
   st%zm_neu(j,i)     = st%zb(j,i)
   st%H_c_neu(j,i)    = st%H_c(j,i)
   st%H_t_neu(j,i)    = st%H_t(j,i)

   call calc_temp_r(st,grd,par,tmr,i,j)

end if

!  ------ Upper right corner

i=grd%IMAX
j=grd%JMAX

if ( (st%maske(j,i).eq.0).or.(st%maske(j,i).eq.3) ) then
                              ! glaciated land or floating ice
   ii=i-1
   jj=j-1

   do kc=0,grd%KCMAX
      st%temp_c_neu(kc,j,i) = st%temp_c_neu(kc,jj,ii)   ! set cold-ice temperature
      st%age_c_neu(kc,j,i)  = st%age_c_neu(kc,jj,ii)    ! set cold-ice age
   end do

   do kr=0,grd%KRMAX
      st%temp_r_neu(kr,j,i) = st%temp_r_neu(kr,jj,ii)   ! set bedrock temperature
   end do

   st%n_cts_neu(j,i)  = st%n_cts_neu(jj,ii)
   st%kc_cts_neu(j,i) = st%kc_cts_neu(jj,ii)
   st%zm_neu(j,i)     = st%zb(j,i)
   st%H_c_neu(j,i)    = st%H_c(j,i)
   st%H_t_neu(j,i)    = st%H_t(j,i)

else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

   st%n_cts_neu(j,i)  = -1
   st%kc_cts_neu(j,i) = 0
   st%zm_neu(j,i)     = st%zb(j,i)
   st%H_c_neu(j,i)    = st%H_c(j,i)
   st%H_t_neu(j,i)    = st%H_t(j,i)

   call calc_temp_r(st,grd,par,tmr,i,j)

end if

!  ------ Lower and upper margins

do i=1, grd%IMAX-1

!    ---- Lower margin

   j=0

   if ( (st%maske(j,i).eq.0).or.(st%maske(j,i).eq.3) ) then
                                 ! glaciated land or floating ice
      ii=i
      jj=j+1

      do kc=0,grd%KCMAX
         st%temp_c_neu(kc,j,i) = st%temp_c_neu(kc,jj,ii)   ! set cold-ice temperature
         st%age_c_neu(kc,j,i)  = st%age_c_neu(kc,jj,ii)    ! set cold-ice age
      end do

      do kr=0,grd%KRMAX
         st%temp_r_neu(kr,j,i) = st%temp_r_neu(kr,jj,ii)   ! set bedrock temperature
      end do

      st%n_cts_neu(j,i)  = st%n_cts_neu(jj,ii)
      st%kc_cts_neu(j,i) = st%kc_cts_neu(jj,ii)
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = st%H_t(j,i)

   else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

      st%n_cts_neu(j,i)  = -1
      st%kc_cts_neu(j,i) = 0
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = st%H_t(j,i)

      call calc_temp_r(st,grd,par,tmr,i,j)

   end if

!    ---- Upper margin

   j=grd%JMAX

   if ( (st%maske(j,i).eq.0).or.(st%maske(j,i).eq.3) ) then
                                 ! glaciated land or floating ice
      ii=i
      jj=j-1

      do kc=0,grd%KCMAX
         st%temp_c_neu(kc,j,i) = st%temp_c_neu(kc,jj,ii)   ! set cold-ice temperature
         st%age_c_neu(kc,j,i)  = st%age_c_neu(kc,jj,ii)    ! set cold-ice age
      end do

      do kr=0,grd%KRMAX
         st%temp_r_neu(kr,j,i) = st%temp_r_neu(kr,jj,ii)   ! set bedrock temperature
      end do

      st%n_cts_neu(j,i)  = st%n_cts_neu(jj,ii)
      st%kc_cts_neu(j,i) = st%kc_cts_neu(jj,ii)
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = st%H_t(j,i)

   else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

      st%n_cts_neu(j,i)  = -1
      st%kc_cts_neu(j,i) = 0
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = st%H_t(j,i)

      call calc_temp_r(st,grd,par,tmr,i,j)

   end if

end do

!  ------ Left and right margins

do j=1, grd%JMAX-1

!    ---- Left margin

   i=0

   if ( (st%maske(j,i).eq.0).or.(st%maske(j,i).eq.3) ) then
                                 ! glaciated land or floating ice
      ii=i+1
      jj=j

      do kc=0,grd%KCMAX
         st%temp_c_neu(kc,j,i) = st%temp_c_neu(kc,jj,ii)   ! set cold-ice temperature
         st%age_c_neu(kc,j,i)  = st%age_c_neu(kc,jj,ii)    ! set cold-ice age
      end do

      do kr=0,grd%KRMAX
         st%temp_r_neu(kr,j,i) = st%temp_r_neu(kr,jj,ii)   ! set bedrock temperature
      end do

      st%n_cts_neu(j,i)  = st%n_cts_neu(jj,ii)
      st%kc_cts_neu(j,i) = st%kc_cts_neu(jj,ii)
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = st%H_t(j,i)

   else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

      st%n_cts_neu(j,i)  = -1
      st%kc_cts_neu(j,i) = 0
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = st%H_t(j,i)

      call calc_temp_r(st,grd,par,tmr,i,j)

   end if

!    ---- Right margin

   i=grd%IMAX

   if ( (st%maske(j,i).eq.0).or.(st%maske(j,i).eq.3) ) then
                                 ! glaciated land or floating ice
      ii=i-1
      jj=j

      do kc=0,grd%KCMAX
         st%temp_c_neu(kc,j,i) = st%temp_c_neu(kc,jj,ii)   ! set cold-ice temperature
         st%age_c_neu(kc,j,i)  = st%age_c_neu(kc,jj,ii)    ! set cold-ice age
      end do

      do kr=0,grd%KRMAX
         st%temp_r_neu(kr,j,i) = st%temp_r_neu(kr,jj,ii)   ! set bedrock temperature
      end do

      st%n_cts_neu(j,i)  = st%n_cts_neu(jj,ii)
      st%kc_cts_neu(j,i) = st%kc_cts_neu(jj,ii)
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = st%H_t(j,i)

   else   ! maske(j,i).eq.1,2 (ice-free land or sea point)

      st%n_cts_neu(j,i)  = -1
      st%kc_cts_neu(j,i) = 0
      st%zm_neu(j,i)     = st%zb(j,i)
      st%H_c_neu(j,i)    = st%H_c(j,i)
      st%H_t_neu(j,i)    = st%H_t(j,i)

      call calc_temp_r(st,grd,par,tmr,i,j)

   end if

end do

!-------- Dummy values for st%omega_c_neu --------

st%omega_c_neu = 0.0_wp   ! not computed in the cold-ice mode

end subroutine calc_temp_cold

!-------------------------------------------------------------------------------
!> Isothermal mode: Setting of the temperature and age to constant values.
!<------------------------------------------------------------------------------
subroutine calc_temp_const(st,grd,tmr,par)

implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class),  intent(in)    :: grd
  type(sico_timer_class), intent(in)    :: tmr
  type(sico_par_class), intent(in)    :: par

   st%temp_c_neu  = par%temp_const
   st%temp_r_neu  = par%temp_const

st%temp_c_neu = min(st%temp_c_neu, st%temp_c_m-eps)
             ! keep temperatures below the pressure melting point

st%omega_t_neu = 0.0_wp
st%omega_c_neu = 0.0_wp

st%Q_tld       = 0.0_wp

   st%age_c_neu   = par%age_const *sec_year   ! a --> s
   st%age_t_neu   = par%age_const *sec_year   ! a --> s

st%n_cts_neu   = -1
st%kc_cts_neu  = 0
st%zm_neu      = st%zb
st%H_c_neu     = st%H_c
st%H_t_neu     = 0.0_wp

end subroutine calc_temp_const

!-------------------------------------------------------------------------------
!> Computation of temperature and age for a cold ice column.
!<------------------------------------------------------------------------------
subroutine calc_temp1(st,grd,tmr,par,i,j)

implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class),  intent(in)    :: grd
  type(sico_timer_class), intent(in)    :: tmr
  type(sico_par_class), intent(in)    :: par

integer, intent(in) :: i, j

integer :: kc, kt, kr
real(wp) :: ct1(0:100), ct2(0:100), ct3(0:100), ct4(0:100), &
            ct5(0:100), ct6(0:100), ct7(0:100), ctr1, &
            ccb1, ccb2, ccb3, ccb4, clb1
real(wp) :: ct1_sg(0:100), ct2_sg(0:100), ct3_sg(0:100), &
            ct4_sg(0:100), adv_vert_sg(0:100), abs_adv_vert_sg(0:100)
real(wp) :: ci1(0:100), ci2(0:100)
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
   stop ' >>> calc_temp1: Boundary points not allowed.'

!-------- Abbreviations --------

ctr1 = grd%atr1

ccb1 = grd%acb1 &
   *kappa_val(st%temp_c(0,j,i)) &
   /st%H_c(j,i)
ccb2 = grd%acb2

if (par%dynamics==2) then

  if (.not.st%flag_shelfy_stream(j,i)) then

    ccb3 = grd%acb3*0.5_wp*(st%vx_t(0,j,i)+st%vx_t(0,j,i-1)) &
      *st%H_c(j,i)*st%dzs_dxi_g(j,i)
    ccb4 = grd%acb4*0.5_wp*(st%vy_t(0,j,i)+st%vy_t(0,j-1,i)) &
      *st%H_c(j,i)*st%dzs_deta_g(j,i)

  else   ! st%flag_shelfy_stream(j,i) == .true.

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
   ct5(kc) = grd%at5(kc) &
             /c_val(st%temp_c(kc,j,i)) &
             /st%H_c(j,i)

   if (par%dynamics==2) then
     if (.not.st%flag_shelfy_stream(j,i)) then
       ct7(kc) = grd%at7 &
         /c_val(st%temp_c(kc,j,i)) &
         *st%enh_c(kc,j,i) &
         *ratefac_c(st%temp_c(kc,j,i), st%temp_c_m(kc,j,i)) &
         *creep(par%fin_visc,par%flow_law,st%sigma_c(kc,j,i)) &
         *st%sigma_c(kc,j,i)*st%sigma_c(kc,j,i)
     else
       ct7(kc) = 2.0_wp*grd%at7 &
         /c_val(st%temp_c(kc,j,i)) &
         *viscosity(par%fin_visc,par%flow_law,st%de_c(kc,j,i), &
         st%temp_c(kc,j,i), st%temp_c_m(kc,j,i), 0.0_wp, &
         st%enh_c(kc,j,i), 0) &
         *st%de_c(kc,j,i)**2
     end if
   else ! par%dynamics.ne.2
     ct7(kc) = grd%at7 &
       /c_val(st%temp_c(kc,j,i)) &
       *st%enh_c(kc,j,i) &
       *ratefac_c(st%temp_c(kc,j,i), st%temp_c_m(kc,j,i)) &
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
   st%temp_c_help(kc) = 0.5_wp*(st%temp_c(kc,j,i)+st%temp_c(kc+1,j,i))
   ct6(kc) = grd%at6(kc) &
    *kappa_val(st%temp_c_help(kc)) &
    /st%H_c(j,i)
   ci2(kc) = grd%ai2(kc)/st%H_c(j,i)
end do

if (par%adv_hor==3) then
dtt_dxi  = 2.0_wp*grd%dtt_2dxi
dtt_deta = 2.0_wp*grd%dtt_2deta
endif

!-------- Set up the temperature equations (ice and bedrock
!         simultaneously) --------

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
kc=0
lgs_a0(kr) = ccb2
lgs_a1(kr) = -(ccb1+ccb2)
lgs_a2(kr) = ccb1
lgs_b(kr)  = ccb3+ccb4

do kc=1, grd%KCMAX-1

if (par%adv_vert==1) then

   lgs_a0(grd%KRMAX+kc) = -0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ct5(kc)*ct6(kc-1)
   lgs_a1(grd%KRMAX+kc) = 1.0_wp+ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(grd%KRMAX+kc) = 0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ct5(kc)*ct6(kc)

else if (par%adv_vert==2) then

   lgs_a0(grd%KRMAX+kc) &
         = -0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -ct5(kc)*ct6(kc-1)
   lgs_a1(grd%KRMAX+kc) &
         = 1.0_wp &
           +0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           +ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(grd%KRMAX+kc) &
         =  0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           -ct5(kc)*ct6(kc)

else if (par%adv_vert==3) then

   adv_vert_help = 0.5_wp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(grd%KRMAX+kc) &
         = -max(adv_vert_help, 0.0_wp) &
           -ct5(kc)*ct6(kc-1)
   lgs_a1(grd%KRMAX+kc) &
         = 1.0_wp &
           +max(adv_vert_help, 0.0_wp)-min(adv_vert_help, 0.0_wp) &
           +ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(grd%KRMAX+kc) &
         =  min(adv_vert_help, 0.0_wp) &
           -ct5(kc)*ct6(kc)

endif

if (par%adv_hor==2) then

   lgs_b(grd%KRMAX+kc) = st%temp_c(kc,j,i) + ct7(kc) &
       -grd%dtt_2dxi* &
          ( (st%vx_c(kc,j,i)-abs(st%vx_c(kc,j,i))) &
            *(st%temp_c(kc,j,i+1)-st%temp_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_c(kc,j,i-1)+abs(st%vx_c(kc,j,i-1))) &
            *(st%temp_c(kc,j,i)-st%temp_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_c(kc,j,i)-abs(st%vy_c(kc,j,i))) &
            *(st%temp_c(kc,j+1,i)-st%temp_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_c(kc,j-1,i)+abs(st%vy_c(kc,j-1,i))) &
            *(st%temp_c(kc,j,i)-st%temp_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

   vx_c_help = 0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1))
   vy_c_help = 0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i))

   lgs_b(grd%KRMAX+kc) = st%temp_c(kc,j,i) + ct7(kc) &
       -dtt_dxi* &
          ( min(vx_c_help, 0.0_wp) &
            *(st%temp_c(kc,j,i+1)-st%temp_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_c_help, 0.0_wp) &
            *(st%temp_c(kc,j,i)-st%temp_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_c_help, 0.0_wp) &
            *(st%temp_c(kc,j+1,i)-st%temp_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_c_help, 0.0_wp) &
            *(st%temp_c(kc,j,i)-st%temp_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

end do

kc=grd%KCMAX
lgs_a0(grd%KRMAX+kc) = 0.0_wp
lgs_a1(grd%KRMAX+kc) = 1.0_wp
lgs_b(grd%KRMAX+kc)  = st%temp_s(j,i)

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KCMAX+grd%KRMAX)

!-------- Assign the result --------

do kr=0, grd%KRMAX
   st%temp_r_neu(kr,j,i) = lgs_x(kr)
end do

do kc=0, grd%KCMAX
   st%temp_c_neu(kc,j,i) = lgs_x(grd%KRMAX+kc)
end do

!-------- Set water content in the non-existing temperate layer
!         to zero --------

do kt=0, grd%KTMAX
   st%omega_t_neu(kt,j,i) = 0.0_wp
end do

!-------- Water drainage from the non-existing temperate layer --------

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

!-------- Age of the ice in the non-existing temperate layer --------

do kt=0, grd%KTMAX
   st%age_t_neu(kt,j,i) = st%age_c_neu(0,j,i)
end do

end subroutine calc_temp1

!-------------------------------------------------------------------------------
!> Computation of temperature and age for an ice column with a temperate base
!! overlain by cold ice.
!<------------------------------------------------------------------------------
subroutine calc_temp2(st,grd,tmr,par,i,j)

implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class),  intent(in)    :: grd
  type(sico_timer_class), intent(in)    :: tmr
  type(sico_par_class), intent(in)    :: par

integer, intent(in) :: i, j

integer :: kc, kt, kr
real(wp) :: ct1(0:100), ct2(0:100), ct3(0:100), ct4(0:100), &
            ct5(0:100), ct6(0:100), ct7(0:100), ctr1, clb1
real(wp) :: ct1_sg(0:100), ct2_sg(0:100), ct3_sg(0:100), &
            ct4_sg(0:100), adv_vert_sg(0:100), abs_adv_vert_sg(0:100)
real(wp) :: ci1(0:100), ci2(0:100)
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
   stop ' >>> calc_temp2: Boundary points not allowed.'

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
   ct5(kc) = grd%at5(kc) &
             /c_val(st%temp_c(kc,j,i)) &
             /st%H_c(j,i)

   if (par%dynamics==2) then

     if (.not.st%flag_shelfy_stream(j,i)) then
       ct7(kc) = grd%at7 &
         /c_val(st%temp_c(kc,j,i)) &
         *st%enh_c(kc,j,i) &
         *ratefac_c(st%temp_c(kc,j,i), st%temp_c_m(kc,j,i)) &
         *creep(par%fin_visc,par%flow_law,st%sigma_c(kc,j,i)) &
         *st%sigma_c(kc,j,i)*st%sigma_c(kc,j,i)
     else
       ct7(kc) = 2.0_wp*grd%at7 &
         /c_val(st%temp_c(kc,j,i)) &
         *viscosity(par%fin_visc,par%flow_law,st%de_c(kc,j,i), &
         st%temp_c(kc,j,i), st%temp_c_m(kc,j,i), 0.0_wp, &
         st%enh_c(kc,j,i), 0) &
         *st%de_c(kc,j,i)**2
     end if

   else   ! par%dynamics.ne.2

     ct7(kc) = grd%at7 &
       /c_val(st%temp_c(kc,j,i)) &
       *st%enh_c(kc,j,i) &
       *ratefac_c(st%temp_c(kc,j,i), st%temp_c_m(kc,j,i)) &
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
   st%temp_c_help(kc) = 0.5_wp*(st%temp_c(kc,j,i)+st%temp_c(kc+1,j,i))
   ct6(kc) = grd%at6(kc) &
    *kappa_val(st%temp_c_help(kc)) &
    /st%H_c(j,i)
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
lgs_b(kr)   = st%temp_t_m(0,j,i)

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
lgs_b(kc)  = st%temp_c_m(0,j,i)

do kc=1, grd%KCMAX-1

if (par%adv_vert==1) then

   lgs_a0(kc) = -0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) = 1.0_wp+ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) = 0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ct5(kc)*ct6(kc)

else if (par%adv_vert==2) then

   lgs_a0(kc) &
         = -0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) &
         = 1.0_wp &
           +0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           +ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) &
         =  0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           -ct5(kc)*ct6(kc)

else if (par%adv_vert==3) then

   adv_vert_help = 0.5_wp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) &
         = -max(adv_vert_help, 0.0_wp) &
           -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) &
         = 1.0_wp &
           +max(adv_vert_help, 0.0_wp)-min(adv_vert_help, 0.0_wp) &
           +ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) &
         =  min(adv_vert_help, 0.0_wp) &
           -ct5(kc)*ct6(kc)

endif

if (par%adv_hor==2) then

   lgs_b(kc) = st%temp_c(kc,j,i) + ct7(kc) &
       -grd%dtt_2dxi* &
          ( (st%vx_c(kc,j,i)-abs(st%vx_c(kc,j,i))) &
            *(st%temp_c(kc,j,i+1)-st%temp_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_c(kc,j,i-1)+abs(st%vx_c(kc,j,i-1))) &
            *(st%temp_c(kc,j,i)-st%temp_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_c(kc,j,i)-abs(st%vy_c(kc,j,i))) &
            *(st%temp_c(kc,j+1,i)-st%temp_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_c(kc,j-1,i)+abs(st%vy_c(kc,j-1,i))) &
            *(st%temp_c(kc,j,i)-st%temp_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

   vx_c_help = 0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1))
   vy_c_help = 0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i))

   lgs_b(kc) = st%temp_c(kc,j,i) + ct7(kc) &
       -dtt_dxi* &
          ( min(vx_c_help, 0.0_wp) &
            *(st%temp_c(kc,j,i+1)-st%temp_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_c_help, 0.0_wp) &
            *(st%temp_c(kc,j,i)-st%temp_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_c_help, 0.0_wp) &
            *(st%temp_c(kc,j+1,i)-st%temp_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_c_help, 0.0_wp) &
            *(st%temp_c(kc,j,i)-st%temp_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

end do

kc=grd%KCMAX
lgs_a0(kc) = 0.0_wp
lgs_a1(kc) = 1.0_wp
lgs_b(kc)  = st%temp_s(j,i)

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KCMAX)

!-------- Assign the result --------

do kc=0, grd%KCMAX
   st%temp_c_neu(kc,j,i) = lgs_x(kc)
end do

!-------- Set water content in the non-existing temperate layer
!         to zero --------

do kt=0, grd%KTMAX
   st%omega_t_neu(kt,j,i) = 0.0_wp
end do

!-------- Water drainage from the non-existing temperate layer --------

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

!-------- Age of the ice in the non-existing temperate layer --------

do kt=0, grd%KTMAX
   st%age_t_neu(kt,j,i) = st%age_c_neu(0,j,i)
end do

end subroutine calc_temp2

!-------------------------------------------------------------------------------
!> Computation of temperature, water content and age for an ice column with a
!! temperate base overlain by a temperate-ice layer.
!<------------------------------------------------------------------------------
subroutine calc_temp3(st,grd,tmr,par,i,j)

implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class),  intent(in)    :: grd
  type(sico_timer_class), intent(in)    :: tmr
  type(sico_par_class), intent(in)    :: par

integer, intent(in) :: i, j

integer :: kc, kt, kr
real(wp) :: ct1(0:100), ct2(0:100), ct3(0:100), ct4(0:100), &
            ct5(0:100), ct6(0:100), ct7(0:100), ctr1, cm1, cm2, &
            clb1
real(wp) :: ct1_sg(0:100), ct2_sg(0:100), ct3_sg(0:100), &
            ct4_sg(0:100), adv_vert_sg(0:100), abs_adv_vert_sg(0:100)
real(wp) :: ci1(0:100), ci2(0:100), ci3
real(wp) :: cw1(0:100), cw2(0:100), cw3(0:100), cw4(0:100), &
            cw5, cw7(0:100), cw8, cw9(0:100)
real(wp) :: cw1_sg(0:100), cw2_sg(0:100), cw3_sg(0:100), &
            cw4_sg(0:100), adv_vert_w_sg(0:100), abs_adv_vert_w_sg(0:100)
real(wp) :: sigma_c_help(0:100), sigma_t_help(0:100)
real(wp) :: vx_c_help, vy_c_help, vx_t_help, vy_t_help
real(wp) :: adv_vert_help, adv_vert_w_help
real(wp) :: dtt_dxi, dtt_deta
real(wp) :: lgs_a0(0:200), &
            lgs_a1(0:200), &
            lgs_a2(0:200), &
            lgs_x(0:200), &
            lgs_b(0:200)
real(wp), parameter :: zero=0.0_wp

!-------- Check for boundary points --------

if ((i == 0).or.(i == grd%IMAX).or.(j == 0).or.(j == grd%JMAX)) &
   stop ' >>> calc_temp3: Boundary points not allowed.'

!-------- Abbreviations --------

ctr1 = grd%atr1
cm1  = grd%am1*st%H_c_neu(j,i)
clb1 = grd%alb1*st%q_geo(j,i)

if (par%adv_vert==1) then

do kc=1, grd%KCMAX-1
   ct1(kc) = grd%at1(kc)/st%H_c_neu(j,i)*0.5_wp*(st%vz_c(kc,j,i)+st%vz_c(kc-1,j,i))
end do

kc=0
ct1_sg(kc) = 0.5_wp*(grd%at1(kc)+grd%at1(kc+1))/st%H_c_neu(j,i)*st%vz_c(kc,j,i)
             ! only needed for kc=0 ...
kc=grd%KCMAX-1
ct1_sg(kc) = 0.5_wp*(grd%at1(kc)+grd%at1(kc+1))/st%H_c_neu(j,i)*st%vz_c(kc,j,i)
             ! ... and kc=grd%KCMAX-1

else if (par%adv_vert==2 .or. par%adv_vert==3) then

do kc=0, grd%KCMAX-1
   ct1_sg(kc) = 0.5_wp*(grd%at1(kc)+grd%at1(kc+1))/st%H_c_neu(j,i)*st%vz_c(kc,j,i)
end do

endif

do kc=0, grd%KCMAX

   ct2(kc) = ( grd%at2_1(kc)*st%dzm_dtau(j,i) &
           +grd%at2_2(kc)*st%dH_c_dtau(j,i) )/st%H_c_neu(j,i)
   ct3(kc) = ( grd%at3_1(kc)*st%dzm_dxi_g(j,i) &
           +grd%at3_2(kc)*st%dH_c_dxi_g(j,i) )/st%H_c_neu(j,i) &
          *0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1)) *grd%insq_g11_g(j,i)
   ct4(kc) = ( grd%at4_1(kc)*st%dzm_deta_g(j,i) &
            +grd%at4_2(kc)*st%dH_c_deta_g(j,i) )/st%H_c_neu(j,i) &
          *0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i)) *grd%insq_g22_g(j,i)
   ct5(kc) = grd%at5(kc) &
             /c_val(st%temp_c(kc,j,i)) &
             /st%H_c_neu(j,i)

   sigma_c_help(kc) &
           = RHO*G*st%H_c_neu(j,i)*(1.0_wp-grd%eaz_c_quotient(kc)) &
             *sqrt(st%dzs_dxi_g(j,i)**2+st%dzs_deta_g(j,i)**2)

   if (par%dynamics==2) then

     if (.not.st%flag_shelfy_stream(j,i)) then
       ct7(kc) = grd%at7 &
         /c_val(st%temp_c(kc,j,i)) &
         *st%enh_c(kc,j,i) &
         *ratefac_c(st%temp_c(kc,j,i), st%temp_c_m(kc,j,i)) &
         *creep(par%fin_visc,par%flow_law,sigma_c_help(kc)) &
         *sigma_c_help(kc)*sigma_c_help(kc)
     else
       ct7(kc) = 2.0_wp*grd%at7 &
         /c_val(st%temp_c(kc,j,i)) &
         *viscosity(par%fin_visc,par%flow_law,st%de_c(kc,j,i), &
         st%temp_c(kc,j,i), st%temp_c_m(kc,j,i), 0.0_wp, &
         st%enh_c(kc,j,i), 0) &
         *st%de_c(kc,j,i)**2
     end if

   else   ! par%dynamics.ne.2

     ct7(kc) = grd%at7 &
       /c_val(st%temp_c(kc,j,i)) &
       *st%enh_c(kc,j,i) &
       *ratefac_c(st%temp_c(kc,j,i), st%temp_c_m(kc,j,i)) &
       *creep(par%fin_visc,par%flow_law,sigma_c_help(kc)) &
       *sigma_c_help(kc)*sigma_c_help(kc)

   endif

   ci1(kc) = grd%ai1(kc)/st%H_c_neu(j,i)

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
   st%temp_c_help(kc) = 0.5_wp*(st%temp_c(kc,j,i)+st%temp_c(kc+1,j,i))
   ct6(kc) = grd%at6(kc) &
    *kappa_val(st%temp_c_help(kc)) &
    /st%H_c_neu(j,i)
   ci2(kc) = grd%ai2(kc)/st%H_c_neu(j,i)
end do

cw5 = grd%aw5/(st%H_t_neu(j,i)**2)
cw8 = grd%aw8
ci3 = grd%ai3/(st%H_t_neu(j,i)**2)

if (par%adv_vert==1) then

do kt=1, grd%KTMAX-1
   cw1(kt) = grd%aw1/st%H_t_neu(j,i)*0.5_wp*(st%vz_t(kt,j,i)+st%vz_t(kt-1,j,i))
end do

kt=grd%KTMAX
cw1(kt) = grd%aw1/st%H_t_neu(j,i)*0.5_wp*(st%vz_t(kt-1,j,i)+st%vz_c(0,j,i))

kt=0
cw1_sg(kt) = grd%aw1/st%H_t_neu(j,i)*st%vz_t(kt,j,i)
             ! only needed for kt=0 ...
kt=grd%KTMAX-1
cw1_sg(kt) = grd%aw1/st%H_t_neu(j,i)*st%vz_t(kt,j,i)
             ! ... and kt=grd%KTMAX-1

else if (par%adv_vert==2 .or. par%adv_vert==3) then

do kt=0, grd%KTMAX-1
   cw1_sg(kt) = grd%aw1/st%H_t_neu(j,i)*st%vz_t(kt,j,i)
end do

endif

do kt=1, grd%KTMAX-1
   cw9(kt) = grd%aw9 &
             *c_val(st%temp_t_m(kt,j,i)) &
    *( st%dzs_dtau(j,i) &
      +0.5_wp*(st%vx_t(kt,j,i)+st%vx_t(kt,j,i-1))*st%dzs_dxi_g(j,i) &
      +0.5_wp*(st%vy_t(kt,j,i)+st%vy_t(kt,j-1,i))*st%dzs_deta_g(j,i) &
      -0.5_wp*(st%vz_t(kt,j,i)+st%vz_t(kt-1,j,i)) )
end do

do kt=0, grd%KTMAX

   cw2(kt) = grd%aw2*(st%dzb_dtau(j,i)+grd%zeta_t(kt)*st%dH_t_dtau(j,i)) &
             /st%H_t_neu(j,i)
   cw3(kt) = grd%aw3*(st%dzb_dxi_g(j,i)+grd%zeta_t(kt)*st%dH_t_dxi_g(j,i)) &
             /st%H_t_neu(j,i) &
             *0.5_wp*(st%vx_t(kt,j,i)+st%vx_t(kt,j,i-1)) *grd%insq_g11_g(j,i)
   cw4(kt) = grd%aw4*(st%dzb_deta_g(j,i)+grd%zeta_t(kt)*st%dH_t_deta_g(j,i)) &
             /st%H_t_neu(j,i) &
             *0.5_wp*(st%vy_t(kt,j,i)+st%vy_t(kt,j-1,i)) *grd%insq_g22_g(j,i)
   sigma_t_help(kt) &
           = sigma_c_help(0) &
             + RHO*G*st%H_t_neu(j,i)*(1.0_wp-grd%zeta_t(kt)) &
               *sqrt(st%dzs_dxi_g(j,i)**2+st%dzs_deta_g(j,i)**2)

   if (par%dynamics==2) then

     if (.not.st%flag_shelfy_stream(j,i)) then
       cw7(kt) = grd%aw7 &
         *st%enh_t(kt,j,i) &
         *ratefac_t(st%omega_t(kt,j,i)) &
         *creep(par%fin_visc,par%flow_law,sigma_t_help(kt)) &
         *sigma_t_help(kt)*sigma_t_help(kt)
     else
       cw7(kt) = 2.0_wp*grd%aw7 &
         *viscosity(par%fin_visc,par%flow_law,st%de_t(kt,j,i), &
         st%temp_t_m(kt,j,i), st%temp_t_m(kt,j,i), &
         st%omega_t(kt,j,i), &
         st%enh_t(kt,j,i), 1) &
         *st%de_t(kt,j,i)**2
     end if

   else   ! par%dynamics.ne.2

     cw7(kt) = grd%aw7 &
       *st%enh_t(kt,j,i) &
       *ratefac_t(st%omega_t(kt,j,i)) &
       *creep(par%fin_visc,par%flow_law,sigma_t_help(kt)) &
       *sigma_t_help(kt)*sigma_t_help(kt)

   endif

end do

if (par%adv_vert==1) then

kt=0
cw2_sg(kt) = 0.5_wp*(cw2(kt)+cw2(kt+1))
cw3_sg(kt) = 0.5_wp*(cw3(kt)+cw3(kt+1))
cw4_sg(kt) = 0.5_wp*(cw4(kt)+cw4(kt+1))
adv_vert_w_sg(kt) = cw1_sg(kt)-cw2_sg(kt)-cw3_sg(kt)-cw4_sg(kt)
abs_adv_vert_w_sg(kt) = abs(adv_vert_w_sg(kt))   ! only needed for kt=0 ...
kt=grd%KTMAX-1
cw2_sg(kt) = 0.5_wp*(cw2(kt)+cw2(kt+1))
cw3_sg(kt) = 0.5_wp*(cw3(kt)+cw3(kt+1))
cw4_sg(kt) = 0.5_wp*(cw4(kt)+cw4(kt+1))
adv_vert_w_sg(kt) = cw1_sg(kt)-cw2_sg(kt)-cw3_sg(kt)-cw4_sg(kt)
abs_adv_vert_w_sg(kt) = abs(adv_vert_w_sg(kt))   ! ... and kt=grd%KTMAX-1

else if (par%adv_vert==2 .or. par%adv_vert==3) then

do kt=0, grd%KTMAX-1
   cw2_sg(kt) = 0.5_wp*(cw2(kt)+cw2(kt+1))
   cw3_sg(kt) = 0.5_wp*(cw3(kt)+cw3(kt+1))
   cw4_sg(kt) = 0.5_wp*(cw4(kt)+cw4(kt+1))
   adv_vert_w_sg(kt) = cw1_sg(kt)-cw2_sg(kt)-cw3_sg(kt)-cw4_sg(kt)
   abs_adv_vert_w_sg(kt) = abs(adv_vert_w_sg(kt))
end do

endif

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
lgs_b(kr)   = st%temp_t_m(0,j,i)

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KRMAX)

!-------- Assign the result --------

do kr=0, grd%KRMAX
   st%temp_r_neu(kr,j,i) = lgs_x(kr)
end do

!-------- Set up the equations for the water content in
!         temperate ice --------

kt=0
lgs_a1(kt) = 1.0_wp
lgs_a2(kt) = -1.0_wp
lgs_b(kt)  = 0.0_wp

do kt=1, grd%KTMAX-1

if (par%adv_vert==1) then

   lgs_a0(kt) = -0.5_wp*(cw1(kt)-cw2(kt)-cw3(kt)-cw4(kt)) - cw5
   lgs_a1(kt) = 1.0_wp + 2.0_wp*cw5
   lgs_a2(kt) = 0.5_wp*(cw1(kt)-cw2(kt)-cw3(kt)-cw4(kt)) - cw5

else if (par%adv_vert==2) then

   lgs_a0(kt) &
         = -0.5_wp*(adv_vert_w_sg(kt-1)+abs_adv_vert_w_sg(kt-1)) &
           -cw5
   lgs_a1(kt) &
         = 1.0_wp &
           +0.5_wp*(adv_vert_w_sg(kt-1)+abs_adv_vert_w_sg(kt-1)) &
           -0.5_wp*(adv_vert_w_sg(kt)  -abs_adv_vert_w_sg(kt)  ) &
           +2.0_wp*cw5
   lgs_a2(kt) &
         =  0.5_wp*(adv_vert_w_sg(kt)  -abs_adv_vert_w_sg(kt)  ) &
           -cw5

else if (par%adv_vert==3) then

   adv_vert_w_help = 0.5_wp*(adv_vert_w_sg(kt)+adv_vert_w_sg(kt-1))

   lgs_a0(kt) &
         = -max(adv_vert_w_help, 0.0_wp) &
           -cw5
   lgs_a1(kt) &
         = 1.0_wp &
           +max(adv_vert_w_help, 0.0_wp)-min(adv_vert_w_help, 0.0_wp) &
           +2.0_wp*cw5
   lgs_a2(kt) &
         =  min(adv_vert_w_help, 0.0_wp) &
           -cw5

endif

if (par%adv_hor==2) then

   lgs_b(kt) = st%omega_t(kt,j,i) + cw7(kt) + cw8 + cw9(kt) &
       -grd%dtt_2dxi* &
          ( (st%vx_t(kt,j,i)-abs(st%vx_t(kt,j,i))) &
            *(st%omega_t(kt,j,i+1)-st%omega_t(kt,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_t(kt,j,i-1)+abs(st%vx_t(kt,j,i-1))) &
            *(st%omega_t(kt,j,i)-st%omega_t(kt,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_t(kt,j,i)-abs(st%vy_t(kt,j,i))) &
            *(st%omega_t(kt,j+1,i)-st%omega_t(kt,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_t(kt,j-1,i)+abs(st%vy_t(kt,j-1,i))) &
            *(st%omega_t(kt,j,i)-st%omega_t(kt,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

   vx_t_help = 0.5_wp*(st%vx_t(kt,j,i)+st%vx_t(kt,j,i-1))
   vy_t_help = 0.5_wp*(st%vy_t(kt,j,i)+st%vy_t(kt,j-1,i))

   lgs_b(kt) = st%omega_t(kt,j,i) + cw7(kt) + cw8 + cw9(kt) &
       -dtt_dxi* &
          ( min(vx_t_help, 0.0_wp) &
            *(st%omega_t(kt,j,i+1)-st%omega_t(kt,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_t_help, 0.0_wp) &
            *(st%omega_t(kt,j,i)-st%omega_t(kt,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_t_help, 0.0_wp) &
            *(st%omega_t(kt,j+1,i)-st%omega_t(kt,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_t_help, 0.0_wp) &
            *(st%omega_t(kt,j,i)-st%omega_t(kt,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

end do

kt=grd%KTMAX

if (par%cts_melting_freezing==1) then

  if (st%am_perp(j,i) >= zero) then   ! melting condition
    lgs_a0(kt) = 0.0_wp
    lgs_a1(kt) = 1.0_wp
    lgs_b(kt)  = 0.0_wp
  else   ! st%am_perp(j,i) < 0.0, freezing condition
    lgs_a0(kt) = -1.0_wp
    lgs_a1(kt) = 1.0_wp
    lgs_b(kt)  = 0.0_wp
  end if

else if (par%cts_melting_freezing==2) then

  lgs_a0(kt) = 0.0_wp
  lgs_a1(kt) = 1.0_wp   ! melting condition assumed
  lgs_b(kt)  = 0.0_wp

else
  stop ' >>> calc_temp3: CTS_MELTING_FREEZING must be either 1 or 2!'
endif

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KTMAX)

!-------- Assign the result, compute the water drainage --------

st%Q_tld(j,i) = 0.0_wp

do kt=0, grd%KTMAX

   if (lgs_x(kt) < zero) then
      st%omega_t_neu(kt,j,i) = 0.0_wp   ! (as a precaution)
   else if (lgs_x(kt) < omega_MAX) then
      st%omega_t_neu(kt,j,i) = lgs_x(kt)
   else
      st%omega_t_neu(kt,j,i) = omega_MAX
      st%Q_tld(j,i) = st%Q_tld(j,i) &
                     +grd%aqtld*st%H_t_neu(j,i)*(lgs_x(kt)-omega_MAX)
   end if

end do

!-------- Set up the equations for the ice temperature --------

!  ------ Abbreviation for the jump of the temperature gradient with
!         the new omega

if (par%cts_melting_freezing==1) then

  if (st%am_perp(j,i) >= zero) then   ! melting condition
    cm2 = 0.0_wp
  else   ! st%am_perp(j,i) < 0.0, freezing condition
    cm2  = grd%am2*st%H_c_neu(j,i)*st%omega_t_neu(grd%KTMAX,j,i)*st%am_perp(j,i) &
      /kappa_val(st%temp_c(0,j,i))
  end if

else if (par%cts_melting_freezing==2) then

  cm2 = 0.0_wp   ! melting condition assumed

else
  stop ' >>> calc_temp3: CTS_MELTING_FREEZING must be either 1 or 2!'
endif

kc=0
lgs_a1(kc) = 1.0_wp
lgs_a2(kc) = -1.0_wp
lgs_b(kc)  = -cm1-cm2

do kc=1, grd%KCMAX-1

if (par%adv_vert==1) then

   lgs_a0(kc) = -0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) = 1.0_wp+ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) = 0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ct5(kc)*ct6(kc)

else if (par%adv_vert==2) then

   lgs_a0(kc) &
         = -0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) &
         = 1.0_wp &
           +0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           +ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) &
         =  0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           -ct5(kc)*ct6(kc)

else if (par%adv_vert==3) then

   adv_vert_help = 0.5_wp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) &
         = -max(adv_vert_help, 0.0_wp) &
           -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) &
         = 1.0_wp &
           +max(adv_vert_help, 0.0_wp)-min(adv_vert_help, 0.0_wp) &
           +ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) &
         =  min(adv_vert_help, 0.0_wp) &
           -ct5(kc)*ct6(kc)

endif

if (par%adv_hor==2) then

   lgs_b(kc) = st%temp_c(kc,j,i) + ct7(kc) &
       -grd%dtt_2dxi* &
          ( (st%vx_c(kc,j,i)-abs(st%vx_c(kc,j,i))) &
            *(st%temp_c(kc,j,i+1)-st%temp_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_c(kc,j,i-1)+abs(st%vx_c(kc,j,i-1))) &
            *(st%temp_c(kc,j,i)-st%temp_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_c(kc,j,i)-abs(st%vy_c(kc,j,i))) &
            *(st%temp_c(kc,j+1,i)-st%temp_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_c(kc,j-1,i)+abs(st%vy_c(kc,j-1,i))) &
            *(st%temp_c(kc,j,i)-st%temp_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

   vx_c_help = 0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1))
   vy_c_help = 0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i))

   lgs_b(kc) = st%temp_c(kc,j,i) + ct7(kc) &
       -dtt_dxi* &
          ( min(vx_c_help, 0.0_wp) &
            *(st%temp_c(kc,j,i+1)-st%temp_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_c_help, 0.0_wp) &
            *(st%temp_c(kc,j,i)-st%temp_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_c_help, 0.0_wp) &
            *(st%temp_c(kc,j+1,i)-st%temp_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_c_help, 0.0_wp) &
            *(st%temp_c(kc,j,i)-st%temp_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

end do

kc=grd%KCMAX
lgs_a0(kc) = 0.0_wp
lgs_a1(kc) = 1.0_wp
lgs_b(kc)  = st%temp_s(j,i)

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KCMAX)

!-------- Assign the result --------

do kc=0, grd%KCMAX
   st%temp_c_neu(kc,j,i) = lgs_x(kc)
end do

!-------- Set up the equations for the age (cold and temperate ice
!         simultaneously) --------

kt=0                                                  ! adv_vert_w_sg(0) <= 0
lgs_a1(kt) = 1.0_wp - min(adv_vert_w_sg(kt), 0.0_wp)  ! (directed downward)
lgs_a2(kt) = min(adv_vert_w_sg(kt), 0.0_wp)           ! assumed/enforced

if (par%adv_hor==2) then

lgs_b(kt) = st%age_t(kt,j,i) + tmr%dtime_temp &
       -grd%dtt_2dxi* &
          ( (st%vx_t(kt,j,i)-abs(st%vx_t(kt,j,i))) &
            *(st%age_t(kt,j,i+1)-st%age_t(kt,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_t(kt,j,i-1)+abs(st%vx_t(kt,j,i-1))) &
            *(st%age_t(kt,j,i)-st%age_t(kt,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_t(kt,j,i)-abs(st%vy_t(kt,j,i))) &
            *(st%age_t(kt,j+1,i)-st%age_t(kt,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_t(kt,j-1,i)+abs(st%vy_t(kt,j-1,i))) &
            *(st%age_t(kt,j,i)-st%age_t(kt,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

vx_t_help = 0.5_wp*(st%vx_t(kt,j,i)+st%vx_t(kt,j,i-1))
vy_t_help = 0.5_wp*(st%vy_t(kt,j,i)+st%vy_t(kt,j-1,i))

lgs_b(kt) = st%age_t(kt,j,i) + tmr%dtime_temp &
       -dtt_dxi* &
          ( min(vx_t_help, 0.0_wp) &
            *(st%age_t(kt,j,i+1)-st%age_t(kt,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_t_help, 0.0_wp) &
            *(st%age_t(kt,j,i)-st%age_t(kt,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_t_help, 0.0_wp) &
            *(st%age_t(kt,j+1,i)-st%age_t(kt,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_t_help, 0.0_wp) &
            *(st%age_t(kt,j,i)-st%age_t(kt,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

do kt=1, grd%KTMAX-1

if (par%adv_vert==1) then

   lgs_a0(kt) = -0.5_wp*(cw1(kt)-cw2(kt)-cw3(kt)-cw4(kt)) - ci3
   lgs_a1(kt) = 1.0_wp + 2.0_wp*ci3
   lgs_a2(kt) = 0.5_wp*(cw1(kt)-cw2(kt)-cw3(kt)-cw4(kt)) - ci3

else if (par%adv_vert==2) then

   lgs_a0(kt) = -0.5_wp*(adv_vert_w_sg(kt-1)+abs_adv_vert_w_sg(kt-1))
   lgs_a1(kt) = 1.0_wp &
               +0.5_wp*(adv_vert_w_sg(kt-1)+abs_adv_vert_w_sg(kt-1)) &
               -0.5_wp*(adv_vert_w_sg(kt)  -abs_adv_vert_w_sg(kt)  )
   lgs_a2(kt) =  0.5_wp*(adv_vert_w_sg(kt)  -abs_adv_vert_w_sg(kt)  )

else if (par%adv_vert==3) then

   adv_vert_w_help = 0.5_wp*(adv_vert_w_sg(kt)+adv_vert_w_sg(kt-1))

   lgs_a0(kt) = -max(adv_vert_w_help, 0.0_wp)
   lgs_a1(kt) = 1.0_wp &
               +max(adv_vert_w_help, 0.0_wp)-min(adv_vert_w_help, 0.0_wp)
   lgs_a2(kt) =  min(adv_vert_w_help, 0.0_wp)

endif

if (par%adv_hor==2) then

   lgs_b(kt) = st%age_t(kt,j,i) + tmr%dtime_temp &
       -grd%dtt_2dxi* &
          ( (st%vx_t(kt,j,i)-abs(st%vx_t(kt,j,i))) &
            *(st%age_t(kt,j,i+1)-st%age_t(kt,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_t(kt,j,i-1)+abs(st%vx_t(kt,j,i-1))) &
            *(st%age_t(kt,j,i)-st%age_t(kt,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_t(kt,j,i)-abs(st%vy_t(kt,j,i))) &
            *(st%age_t(kt,j+1,i)-st%age_t(kt,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_t(kt,j-1,i)+abs(st%vy_t(kt,j-1,i))) &
            *(st%age_t(kt,j,i)-st%age_t(kt,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

   vx_t_help = 0.5_wp*(st%vx_t(kt,j,i)+st%vx_t(kt,j,i-1))
   vy_t_help = 0.5_wp*(st%vy_t(kt,j,i)+st%vy_t(kt,j-1,i))

   lgs_b(kt) = st%age_t(kt,j,i) + tmr%dtime_temp &
       -dtt_dxi* &
          ( min(vx_t_help, 0.0_wp) &
            *(st%age_t(kt,j,i+1)-st%age_t(kt,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_t_help, 0.0_wp) &
            *(st%age_t(kt,j,i)-st%age_t(kt,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_t_help, 0.0_wp) &
            *(st%age_t(kt,j+1,i)-st%age_t(kt,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_t_help, 0.0_wp) &
            *(st%age_t(kt,j,i)-st%age_t(kt,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

end do

if (par%adv_vert==1) then

kt=grd%KTMAX
kc=0

lgs_a0(kt) = -0.5_wp*(cw1(kt)-cw2(kt)-cw3(kt)-cw4(kt)) - ci3
lgs_a1(kt) = 1.0_wp + 2.0_wp*ci3
lgs_a2(kt) = 0.5_wp*(cw1(kt)-cw2(kt)-cw3(kt)-cw4(kt)) - ci3

if (par%adv_hor==2) then

lgs_b(kt) = st%age_t(kt,j,i) + tmr%dtime_temp &
       -grd%dtt_2dxi* &
          ( (st%vx_t(kt,j,i)-abs(st%vx_t(kt,j,i))) &
            *(st%age_t(kt,j,i+1)-st%age_t(kt,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_t(kt,j,i-1)+abs(st%vx_t(kt,j,i-1))) &
            *(st%age_t(kt,j,i)-st%age_t(kt,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_t(kt,j,i)-abs(st%vy_t(kt,j,i))) &
            *(st%age_t(kt,j+1,i)-st%age_t(kt,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_t(kt,j-1,i)+abs(st%vy_t(kt,j-1,i))) &
            *(st%age_t(kt,j,i)-st%age_t(kt,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

vx_t_help = 0.5_wp*(st%vx_t(kt,j,i)+st%vx_t(kt,j,i-1))
vy_t_help = 0.5_wp*(st%vy_t(kt,j,i)+st%vy_t(kt,j-1,i))

lgs_b(kt) = st%age_t(kt,j,i) + tmr%dtime_temp &
       -dtt_dxi* &
          ( min(vx_t_help, 0.0_wp) &
            *(st%age_t(kt,j,i+1)-st%age_t(kt,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_t_help, 0.0_wp) &
            *(st%age_t(kt,j,i)-st%age_t(kt,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_t_help, 0.0_wp) &
            *(st%age_t(kt,j+1,i)-st%age_t(kt,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_t_help, 0.0_wp) &
            *(st%age_t(kt,j,i)-st%age_t(kt,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

else if (par%adv_vert==2 .or. par%adv_vert==3) then

kt=grd%KTMAX
kc=0

if (adv_vert_sg(kc) <= zero) then

   lgs_a0(grd%KTMAX+kc) = 0.0_wp
   lgs_a1(grd%KTMAX+kc) = 1.0_wp - adv_vert_sg(kc)
   lgs_a2(grd%KTMAX+kc) = adv_vert_sg(kc)

if (par%adv_hor==2) then

   lgs_b(grd%KTMAX+kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
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

   lgs_b(grd%KTMAX+kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
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

else if (adv_vert_w_sg(kt-1) >= zero) then

   lgs_a0(kt) = -adv_vert_w_sg(kt-1)
   lgs_a1(kt) = 1.0_wp + adv_vert_w_sg(kt-1)
   lgs_a2(kt) = 0.0_wp

if (par%adv_hor==2) then

   lgs_b(kt) = st%age_t(kt,j,i) + tmr%dtime_temp &
       -grd%dtt_2dxi* &
          ( (st%vx_t(kt,j,i)-abs(st%vx_t(kt,j,i))) &
            *(st%age_t(kt,j,i+1)-st%age_t(kt,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_t(kt,j,i-1)+abs(st%vx_t(kt,j,i-1))) &
            *(st%age_t(kt,j,i)-st%age_t(kt,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_t(kt,j,i)-abs(st%vy_t(kt,j,i))) &
            *(st%age_t(kt,j+1,i)-st%age_t(kt,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_t(kt,j-1,i)+abs(st%vy_t(kt,j-1,i))) &
            *(st%age_t(kt,j,i)-st%age_t(kt,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

   vx_t_help = 0.5_wp*(st%vx_t(kt,j,i)+st%vx_t(kt,j,i-1))
   vy_t_help = 0.5_wp*(st%vy_t(kt,j,i)+st%vy_t(kt,j-1,i))

   lgs_b(kt) = st%age_t(kt,j,i) + tmr%dtime_temp &
       -dtt_dxi* &
          ( min(vx_t_help, 0.0_wp) &
            *(st%age_t(kt,j,i+1)-st%age_t(kt,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_t_help, 0.0_wp) &
            *(st%age_t(kt,j,i)-st%age_t(kt,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_t_help, 0.0_wp) &
            *(st%age_t(kt,j+1,i)-st%age_t(kt,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_t_help, 0.0_wp) &
            *(st%age_t(kt,j,i)-st%age_t(kt,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

else

   lgs_a0(kt) = -0.5_wp
   lgs_a1(kt) = 1.0_wp
   lgs_a2(kt) = -0.5_wp
   lgs_b(kt)  = 0.0_wp
   ! Makeshift: Average of age_c(kc=1) and age_t(kt=grd%KTMAX-1)

end if

endif

do kc=1, grd%KCMAX-1

if (par%adv_vert==1) then

   lgs_a0(grd%KTMAX+kc) = -0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ci1(kc)*ci2(kc-1)
   lgs_a1(grd%KTMAX+kc) = 1.0_wp+ci1(kc)*(ci2(kc)+ci2(kc-1))
   lgs_a2(grd%KTMAX+kc) = 0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ci1(kc)*ci2(kc)

else if (par%adv_vert==2) then

   lgs_a0(grd%KTMAX+kc) = -0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1))
   lgs_a1(grd%KTMAX+kc) =  1.0_wp &
                      +0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
                      -0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  )
   lgs_a2(grd%KTMAX+kc) =  0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  )

else if (par%adv_vert==3) then

   adv_vert_help = 0.5_wp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(grd%KTMAX+kc) = -max(adv_vert_help, 0.0_wp)
   lgs_a1(grd%KTMAX+kc) =  1.0_wp &
                      +max(adv_vert_help, 0.0_wp)-min(adv_vert_help, 0.0_wp)
   lgs_a2(grd%KTMAX+kc) =  min(adv_vert_help, 0.0_wp)

endif

if (par%adv_hor==2) then

   lgs_b(grd%KTMAX+kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
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

   lgs_b(grd%KTMAX+kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
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
   lgs_a0(grd%KTMAX+kc) = 0.0_wp
   lgs_a1(grd%KTMAX+kc) = 1.0_wp
   lgs_b(grd%KTMAX+kc)  = 0.0_wp
else
   lgs_a0(grd%KTMAX+kc) = -max(adv_vert_sg(kc-1), 0.0_wp)
   lgs_a1(grd%KTMAX+kc) = 1.0_wp + max(adv_vert_sg(kc-1), 0.0_wp)
                             ! adv_vert_sg(grd%KCMAX-1) >= 0 (directed upward)
                             ! assumed/enforced
if (par%adv_hor==2) then

   lgs_b(grd%KTMAX+kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
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

   lgs_b(grd%KTMAX+kc) = st%age_c(kc,j,i) + tmr%dtime_temp &
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

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KCMAX+grd%KTMAX)

!-------- Assign the result,
!         restriction to interval [0, AGE_MAX yr] --------

do kt=0, grd%KTMAX

   st%age_t_neu(kt,j,i) = lgs_x(kt)

   if (st%age_t_neu(kt,j,i) < (par%age_min*sec_year)) &
                           st%age_t_neu(kt,j,i) = 0.0_wp
   if (st%age_t_neu(kt,j,i) > (par%age_max*sec_year)) &
                           st%age_t_neu(kt,j,i) = par%age_max*sec_year

end do

do kc=0, grd%KCMAX

   st%age_c_neu(kc,j,i) = lgs_x(grd%KTMAX+kc)

   if (st%age_c_neu(kc,j,i) < (par%age_min*sec_year)) &
                           st%age_c_neu(kc,j,i) = 0.0_wp
   if (st%age_c_neu(kc,j,i) > (par%age_max*sec_year)) &
                           st%age_c_neu(kc,j,i) = par%age_max*sec_year

end do

end subroutine calc_temp3

!-------------------------------------------------------------------------------
!> Computation of temperature and age for an ice-free column.
!<------------------------------------------------------------------------------
subroutine calc_temp_r(st,grd,par,tmr,i,j)

implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class),  intent(in)    :: grd
  type(sico_par_class), intent(in)    :: par
  type(sico_timer_class), intent(in)    :: tmr

integer, intent(in) :: i, j

integer :: kc, kt, kr
real(wp) :: ctr1, clb1
real(wp) :: lgs_a0(0:200), &
            lgs_a1(0:200), &
            lgs_a2(0:200), &
            lgs_x(0:200), &
            lgs_b(0:200)

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

!-------- Water content and age
!                       in the non-existing lower (kt) ice layer --------

do kt=0, grd%KTMAX
   st%omega_t_neu(kt,j,i) = 0.0_wp
   st%age_t_neu(kt,j,i)   = 0.0_wp
end do

!-------- Temperature and age
!                     in the non-existing upper (kc) ice layer --------

do kc=0, grd%KCMAX
   st%temp_c_neu(kc,j,i)  = st%temp_g(j,i)
   st%age_c_neu(kc,j,i)   = 0.0_wp
end do

end subroutine calc_temp_r

!-------------------------------------------------------------------------------
!> Upward shifting of the CTS.
!<------------------------------------------------------------------------------
subroutine shift_cts_upward(st,grd,tmr,par,i,j)

implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class),  intent(in)    :: grd
  type(sico_timer_class), intent(in)    :: tmr
  type(sico_par_class), intent(in)    :: par

integer, intent(in) :: i, j

real(wp) :: zm_shift
real(wp) :: difftemp_a, difftemp_b, interpol

zm_shift = 1.0_wp   ! CTS shift in intervals of 1 m

!-------- Temperature discrepancy from the computation of the main
!         program --------

difftemp_a = st%temp_c_neu(0,j,i)-(-BETA*st%H_c_neu(j,i))
if (difftemp_a <= 0.0_wp) return

!-------- Shift CTS upward until it is too high --------

   do while (difftemp_a > 0.0_dp)

   st%zm_neu(j,i)  = st%zm_neu(j,i)  + zm_shift
   if (st%zm_neu(j,i) >= st%zs(j,i)) then
      st%zm_neu(j,i)  = st%zm_neu(j,i) - zm_shift
      return
   end if
   st%H_c_neu(j,i) = st%H_c_neu(j,i) - zm_shift
   st%H_t_neu(j,i) = st%H_t_neu(j,i) + zm_shift

   st%dH_t_dtau(j,i) = (st%zm_neu(j,i)-st%zm(j,i))*tmr%dtime_temp_inv
   st%dzm_dtau(j,i)  = st%dzb_dtau(j,i)+st%dH_t_dtau(j,i)
   st%dH_c_dtau(j,i) = st%dzs_dtau(j,i)-st%dzm_dtau(j,i)

   st%am_perp(j,i) = st%am_perp_st(j,i) + st%dzm_dtau(j,i)

   call calc_temp3(st,grd,tmr,par,i,j)

   difftemp_b = difftemp_a
   difftemp_a = st%temp_c_neu(0,j,i)-(-BETA*st%H_c_neu(j,i))

 enddo

!-------- Interpolate the CTS position from the last (_a) and the
!         last but one (_b) value, weighed with the temperature
!         discrepancies at the CTS --------

interpol = difftemp_a/(difftemp_b-difftemp_a)*zm_shift

st%zm_neu(j,i)  = st%zm_neu(j,i)  + interpol
st%H_c_neu(j,i) = st%H_c_neu(j,i) - interpol
st%H_t_neu(j,i) = st%H_t_neu(j,i) + interpol

st%dH_t_dtau(j,i) = (st%zm_neu(j,i)-st%zm(j,i))*tmr%dtime_temp_inv
st%dzm_dtau(j,i)  = st%dzb_dtau(j,i)+st%dH_t_dtau(j,i)
st%dH_c_dtau(j,i) = st%dzs_dtau(j,i)-st%dzm_dtau(j,i)

st%am_perp(j,i) = st%am_perp_st(j,i) + st%dzm_dtau(j,i)

call calc_temp3(st,grd,tmr,par,i,j)

end subroutine shift_cts_upward

!-------------------------------------------------------------------------------
!> Downward shifting of the CTS.
!<------------------------------------------------------------------------------
subroutine shift_cts_downward(st,grd,tmr,par,i,j)

implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class),  intent(in)    :: grd
  type(sico_timer_class), intent(in)    :: tmr
  type(sico_par_class), intent(in)    :: par

integer, intent(in) :: i, j

real(wp) :: zm_shift
real(wp) :: difftemp_a, difftemp_b, interpol

zm_shift = 1.0_wp   ! CTS shift in intervals of 1 m

!-------- Temperature discrepancy from the computation of the main
!         program --------

difftemp_a = st%temp_c_neu(0,j,i)-(-BETA*st%H_c_neu(j,i))
if (difftemp_a >= 0.0_wp) return

!-------- Shift CTS downward until it is too low --------

do while (difftemp_a < 0.0_dp)

   st%zm_neu(j,i)  = st%zm_neu(j,i) - zm_shift

!  ------ Special case: CTS too close to the base

   if (st%zm_neu(j,i) <= st%zb(j,i)) then

      zm_shift = (st%zm_neu(j,i)+zm_shift)-(st%zb(j,i)+0.001_wp)
      st%zm_neu(j,i)  = st%zb(j,i)+0.001_wp
      st%H_c_neu(j,i) = st%H_c_neu(j,i)+st%H_t_neu(j,i)-0.001_wp
      st%H_t_neu(j,i) = 0.001_wp
!                   ! CTS positioned 1 mm above ice base --------

      st%dH_t_dtau(j,i) = (st%zm_neu(j,i)-st%zm(j,i))*tmr%dtime_temp_inv
      st%dzm_dtau(j,i)  = st%dzb_dtau(j,i)+st%dH_t_dtau(j,i)
      st%dH_c_dtau(j,i) = st%dzs_dtau(j,i)-st%dzm_dtau(j,i)

      st%am_perp(j,i) = st%am_perp_st(j,i) + st%dzm_dtau(j,i)

      call calc_temp3(st,grd,tmr,par,i,j)

      difftemp_b = difftemp_a
      difftemp_a = st%temp_c_neu(0,j,i)-(-BETA*st%H_c_neu(j,i))

      if (difftemp_a >= 0.0_wp) then ! CTS remains above the base

!    ---- Interpolate the CTS position from the last (_a) and the
!         last but one (_b) value, weighed with the temperature
!         discrepancies at the CTS --------

         interpol = difftemp_a/(difftemp_a-difftemp_b)*zm_shift

         st%zm_neu(j,i)  = st%zm_neu(j,i)  + interpol
         st%H_c_neu(j,i) = st%H_c_neu(j,i) - interpol
         st%H_t_neu(j,i) = st%H_t_neu(j,i) + interpol

         st%dH_t_dtau(j,i) = (st%zm_neu(j,i)-st%zm(j,i))*tmr%dtime_temp_inv
         st%dzm_dtau(j,i)  = st%dzb_dtau(j,i)+st%dH_t_dtau(j,i)
         st%dH_c_dtau(j,i) = st%dzs_dtau(j,i)-st%dzm_dtau(j,i)

         st%am_perp(j,i) = st%am_perp_st(j,i) + st%dzm_dtau(j,i)

         call calc_temp3(st,grd,tmr,par,i,j)

      else   ! CTS disappears

         st%n_cts_neu(j,i) = 0
         st%zm_neu(j,i) = st%zb(j,i)
         st%H_c_neu(j,i) = st%H_c_neu(j,i)+st%H_t_neu(j,i)
         st%H_t_neu(j,i) = 0.0_wp

         st%dH_t_dtau(j,i) = (st%zm_neu(j,i)-st%zm(j,i))*tmr%dtime_temp_inv
         st%dzm_dtau(j,i)  = st%dzb_dtau(j,i)+st%dH_t_dtau(j,i)
         st%dH_c_dtau(j,i) = st%dzs_dtau(j,i)-st%dzm_dtau(j,i)

         st%am_perp(j,i) = st%am_perp_st(j,i) + st%dzm_dtau(j,i)

         call calc_temp2(st,grd,tmr,par,i,j)

      end if

      return
   end if

!  ------ End of treatment of special case

   st%H_c_neu(j,i) = st%H_c_neu(j,i) + zm_shift
   st%H_t_neu(j,i) = st%H_t_neu(j,i) - zm_shift

   st%dH_t_dtau(j,i) = (st%zm_neu(j,i)-st%zm(j,i))*tmr%dtime_temp_inv
   st%dzm_dtau(j,i)  = st%dzb_dtau(j,i)+st%dH_t_dtau(j,i)
   st%dH_c_dtau(j,i) = st%dzs_dtau(j,i)-st%dzm_dtau(j,i)

   st%am_perp(j,i) = st%am_perp_st(j,i) + st%dzm_dtau(j,i)

   call calc_temp3(st,grd,tmr,par,i,j)

   difftemp_b = difftemp_a
   difftemp_a = st%temp_c_neu(0,j,i)-(-BETA*st%H_c_neu(j,i))

 enddo

!-------- Interpolate the CTS position from the last (_a) and the
!         last but one (_b) value, weighed with the temperature
!         discrepancies at the CTS --------

interpol = difftemp_a/(difftemp_a-difftemp_b)*zm_shift

st%zm_neu(j,i)  = st%zm_neu(j,i)  + interpol
st%H_c_neu(j,i) = st%H_c_neu(j,i) - interpol
st%H_t_neu(j,i) = st%H_t_neu(j,i) + interpol

st%dH_t_dtau(j,i) = (st%zm_neu(j,i)-st%zm(j,i))*tmr%dtime_temp_inv
st%dzm_dtau(j,i)  = st%dzb_dtau(j,i)+st%dH_t_dtau(j,i)
st%dH_c_dtau(j,i) = st%dzs_dtau(j,i)-st%dzm_dtau(j,i)

st%am_perp(j,i) = st%am_perp_st(j,i) + st%dzm_dtau(j,i)

call calc_temp3(st,grd,tmr,par,i,j)

end subroutine shift_cts_downward

!-------------------------------------------------------------------------------
!> Computation of temperature and age for ice shelves (floating ice).
!<------------------------------------------------------------------------------
subroutine calc_temp_ssa(st,grd,tmr,par,i,j)

implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class),  intent(in)    :: grd
  type(sico_timer_class), intent(in)    :: tmr
  type(sico_par_class), intent(in)    :: par

integer, intent(in) :: i, j

integer :: kc, kt, kr
real(wp) :: ct1(0:100), ct2(0:100), ct3(0:100), ct4(0:100), &
            ct5(0:100), ct6(0:100), ct7(0:100), ctr1, clb1
real(wp) :: ct1_sg(0:100), ct2_sg(0:100), ct3_sg(0:100), &
            ct4_sg(0:100), adv_vert_sg(0:100), abs_adv_vert_sg(0:100)
real(wp) :: ci1(0:100), ci2(0:100)
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
   stop ' >>> calc_temp_ssa: Boundary points not allowed.'

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
   ct5(kc) = grd%at5(kc) &
             /c_val(st%temp_c(kc,j,i)) &
             /st%H_c(j,i)
   ct7(kc) = 2.0_wp*grd%at7 &
             /c_val(st%temp_c(kc,j,i)) &
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
   st%temp_c_help(kc) = 0.5_wp*(st%temp_c(kc,j,i)+st%temp_c(kc+1,j,i))
   ct6(kc) = grd%at6(kc) &
    *kappa_val(st%temp_c_help(kc)) &
    /st%H_c(j,i)
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
lgs_b(kc)  = st%temp_c_m(0,j,i)-DELTA_TM_SW

do kc=1, grd%KCMAX-1

if (par%adv_vert==1) then

   lgs_a0(kc) = -0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) = 1.0_wp+ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) = 0.5_wp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ct5(kc)*ct6(kc)

else if (par%adv_vert==2) then

   lgs_a0(kc) &
         = -0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) &
         = 1.0_wp &
           +0.5_wp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           +ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) &
         =  0.5_wp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           -ct5(kc)*ct6(kc)

else if (par%adv_vert==3) then

   adv_vert_help = 0.5_wp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) &
         = -max(adv_vert_help, 0.0_wp) &
           -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) &
         = 1.0_wp &
           +max(adv_vert_help, 0.0_wp)-min(adv_vert_help, 0.0_wp) &
           +ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) &
         =  min(adv_vert_help, 0.0_wp) &
           -ct5(kc)*ct6(kc)

endif

if (par%adv_hor==2) then

   lgs_b(kc) = st%temp_c(kc,j,i) + ct7(kc) &
       -grd%dtt_2dxi* &
          ( (st%vx_c(kc,j,i)-abs(st%vx_c(kc,j,i))) &
            *(st%temp_c(kc,j,i+1)-st%temp_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +(st%vx_c(kc,j,i-1)+abs(st%vx_c(kc,j,i-1))) &
            *(st%temp_c(kc,j,i)-st%temp_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -grd%dtt_2deta* &
          ( (st%vy_c(kc,j,i)-abs(st%vy_c(kc,j,i))) &
            *(st%temp_c(kc,j+1,i)-st%temp_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +(st%vy_c(kc,j-1,i)+abs(st%vy_c(kc,j-1,i))) &
            *(st%temp_c(kc,j,i)-st%temp_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

else if (par%adv_hor==3) then

   vx_c_help = 0.5_wp*(st%vx_c(kc,j,i)+st%vx_c(kc,j,i-1))
   vy_c_help = 0.5_wp*(st%vy_c(kc,j,i)+st%vy_c(kc,j-1,i))

   lgs_b(kc) = st%temp_c(kc,j,i) + ct7(kc) &
       -dtt_dxi* &
          ( min(vx_c_help, 0.0_wp) &
            *(st%temp_c(kc,j,i+1)-st%temp_c(kc,j,i)) &
            *grd%insq_g11_sgx(j,i) &
           +max(vx_c_help, 0.0_wp) &
            *(st%temp_c(kc,j,i)-st%temp_c(kc,j,i-1)) &
            *grd%insq_g11_sgx(j,i-1) ) &
       -dtt_deta* &
          ( min(vy_c_help, 0.0_wp) &
            *(st%temp_c(kc,j+1,i)-st%temp_c(kc,j,i)) &
            *grd%insq_g22_sgy(j,i) &
           +max(vy_c_help, 0.0_wp) &
            *(st%temp_c(kc,j,i)-st%temp_c(kc,j-1,i)) &
            *grd%insq_g22_sgy(j-1,i) )

endif

end do

kc=grd%KCMAX
lgs_a0(kc) = 0.0_wp
lgs_a1(kc) = 1.0_wp
lgs_b(kc)  = st%temp_s(j,i)

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, grd%KCMAX)

!-------- Assign the result --------

do kc=0, grd%KCMAX
   st%temp_c_neu(kc,j,i) = lgs_x(kc)
end do

!-------- Set water content in the non-existing temperate layer
!         to zero --------

do kt=0, grd%KTMAX
   st%omega_t_neu(kt,j,i) = 0.0_wp
end do

!-------- Water drainage from the non-existing temperate layer --------

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

!-------- Age of the ice in the non-existing temperate layer --------

do kt=0, grd%KTMAX
   st%age_t_neu(kt,j,i) = st%age_c_neu(0,j,i)
end do

end subroutine calc_temp_ssa

!-------------------------------------------------------------------------------

end module calc_temp_m
!
