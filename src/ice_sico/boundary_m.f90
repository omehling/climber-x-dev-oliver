!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  b o u n d a r y _ m
!
!> @file
!!
!! Computation of the surface temperature (must be less than 0 deg C!)
!! and of the accumulation-ablation function.
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
!> Computation of the surface temperature (must be less than 0 deg C!)
!! and of the accumulation-ablation function.
!<------------------------------------------------------------------------------
module boundary_m

  use sico_types_m
  use sico_grid_mod, only : sico_grid_class
  use sico_state, only : sico_state_class
  use sico_params, only : sico_par_class, rho_sw, RHO
  use sico_timer, only : sico_timer_class
  use discharge_workers_m, only: discharge

  implicit none

  private
  public :: boundary

contains

!-------------------------------------------------------------------------------
!> Main routine of boundary_m:
!! Computation of the surface temperature (must be less than 0 deg C!)
!! and of the accumulation-ablation function.
!<------------------------------------------------------------------------------
subroutine boundary(st,grd,tmr,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_timer_class), intent(inout) :: tmr
type(sico_par_class), intent(in) :: par

integer :: i, j, n
real(wp) :: H_ice, H_sea
real(wp) :: calv_uw_coeff, r1_calv_uw, r2_calv_uw


!-------- Initialization of intent(out) variables --------

st%z_mar      = 0.0_wp

!  ------ Minimum bedrock elevation for extent of marine ice

if (par%margin==2) then

  if ( par%marine_ice_calving==2 .or. par%marine_ice_calving==3 ) then
    st%z_mar = par%z_mar
  else if ( par%marine_ice_calving==4 .or. par%marine_ice_calving==5 ) then
    stop 'marine_ice_calving==4 and marine_ice_calving==5 not implemented yet' 
    !st%z_mar = par%fact_z_mar*st%z_sl
  else if ( par%marine_ice_calving==6 .or. par%marine_ice_calving==7 ) then
    stop 'marine_ice_calving==6 and marine_ice_calving==7 not implemented yet' 
    !if (st%z_sl >= -80.0_wp) then
    !  st%z_mar = 2.5_wp*st%z_sl
    !else
    !  st%z_mar = 10.25_wp*(st%z_sl+80.0_wp)-200.0_wp
    !end if
    !st%z_mar = par%fact_z_mar*st%z_mar
  endif

endif

!  ------ Update of the mask according to the sea level

!    ---- Check all sea and floating-ice points and their direct
!         neighbours

do i=0, grd%IMAX
do j=0, grd%JMAX
   st%check_point(j,i) = .false.
end do
end do

do i=1, grd%IMAX-1
do j=1, grd%JMAX-1
   if (st%maske(j,i).ge.2) then
      st%check_point(j  ,i  ) = .true.
      st%check_point(j  ,i+1) = .true.
      st%check_point(j  ,i-1) = .true.
      st%check_point(j+1,i  ) = .true.
      st%check_point(j-1,i  ) = .true.
   end if
end do
end do

do i=1, grd%IMAX-1
do j=1, grd%JMAX-1
   if (st%check_point(j,i)) then
      !-------- Previously ice-free land point or sea point --------
      if ( (st%maske(j,i) == 1).or.(st%maske(j,i) == 2) ) then
        if (st%zl(j,i) > st%z_sl(j,i)) then
          st%maske_neu(j,i) = 1   ! now ice-free land
        else
          st%maske_neu(j,i) = 2   ! now sea point
        end if
      !-------- Previously grounded-ice or floating-ice point --------
      else   ! (maske(j,i) == 0, 3)
        if (st%zl(j,i) > st%z_sl(j,i)) then
          st%maske_neu(j,i) = 0   ! now grounded ice
        else
          H_ice = st%zs(j,i)-st%zb(j,i)   ! ice thickness
          H_sea = st%z_sl(j,i)-st%zl(j,i)   ! sea depth
          if ( H_ice < (RHO_SW/RHO*H_sea) ) then
            if (par%margin==1 .or. (par%margin==2 .and. par%marine_ice_formation==1)) then
              st%maske_neu(j,i) = 2     ! ice becomes floating, therefore now sea point (ice cut off)
            else if (par%margin==2 .and. par%marine_ice_formation==2) then
              st%maske_neu(j,i) = 0     ! now "underwater ice"
            else if (par%margin==3) then
              st%maske_neu(j,i) = 3     ! now floating ice
            endif
          else
            st%maske_neu(j,i) = 0     ! now grounded ice
          end if
        end if
      end if
   end if
end do
end do

!    ---- Assign new values of the mask

do i=1, grd%IMAX-1
do j=1, grd%JMAX-1
   if (st%check_point(j,i)) then
      st%maske(j,i) = st%maske_neu(j,i)
   end if
end do
end do

!-------- Calving rate of grounded ice --------

st%calving = 0.0_wp

if ((par%margin==2) .and. (par%marine_ice_formation==2) .and. (par%marine_ice_calving==9)) then

!-------- Ice thickness and sea depth --------

  !H     = max(H_c + H_t, 0.0_wp)   ! ice thickness
  !H_sea = max(z_sl - zl, 0.0_wp)   ! sea depth

!-------- Calving of "underwater ice" --------

  where ( (st%maske == 0).and.(max(st%H_c + st%H_t, 0.0_wp) < RHO_SW/RHO*max(st%z_sl - st%zl, 0.0_wp)) )
     st%calving = st%calving + par%calv_uw_coeff * max(st%H_c + st%H_t, 0.0_wp)**par%r1_calv_uw * max(st%z_sl - st%zl, 0.0_wp)**par%r2_calv_uw
  elsewhere
     st%calving = st%calving + 0.0_wp
  end where

endif

!-------- Ice discharge parameterization --------

if (par%i_disc.gt.0) then
  call discharge(st, grd, tmr, par)
  if (par%i_disc.eq.2) then
    ! apply discharge parameterisation only for Greenland
    where(st%id_mask.ne.2) 
      st%dis_perp = 0._wp
    endwhere
  endif
  st%calving = st%calving + st%dis_perp
endif

end subroutine boundary

!-------------------------------------------------------------------------------

end module boundary_m
!
