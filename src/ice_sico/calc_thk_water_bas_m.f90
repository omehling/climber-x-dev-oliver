!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ t h k _ w a t e r _ b a s _ m
!
!> @file
!!
!! Computation of the thickness of the water column under the ice base.
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
!> Computation of the thickness of the water column under the ice base.
!<------------------------------------------------------------------------------
module calc_thk_water_bas_m

  use timer, only : sec_year
  use sico_types_m
  use sico_state
  use sico_timer
  use sico_params, only : sico_par_class, RHO, RHO_w
  !use hydro_m

  implicit none

  private
  public :: calc_thk_water_bas

contains

!-------------------------------------------------------------------------------
!> Main subroutine of calc_thk_water_bas_m:
!! Computation of the thickness of the water column under the ice base.
!<------------------------------------------------------------------------------
  subroutine calc_thk_water_bas(st,par,tmr)

  implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_par_class), intent(in)      :: par
  type(sico_timer_class), intent(in)    :: tmr

  logical, save :: firstcall = .true.

!  real(wp), save                     :: rho_rho_w_ratio
!  integer , dimension(0:IMAX,0:JMAX) :: hydro_icemask
!  real(wp), dimension(0:IMAX,0:JMAX) :: hydro_topg, hydro_thk, &
!                                        hydro_temppabase, hydro_supply, &
!                                        hydro_sflux,  hydro_vflux, &
!                                        hydro_vfluxX, hydro_vfluxY, &
!                                        hydro_bwat
!  type(hydro_t), save                :: hydro
                                        !!! Does this need a save attribute?

!-------- Water column --------

if (par%basal_hydrology==1) then 

  stop 'basal_hydrology==1 not supported yet!'

!  if (firstcall) then
!
!     rho_rho_w_ratio = RHO/RHO_W
!
!     call hydro_init(hydro, xi, eta)
!     call hydro_gen_conf(hydro, &
!          & method='quinn', &
!          & avoid_frz=.false., &
!          & filter_len=0.0_wp, &
!          & rho_seawater=RHO_SW, &
!          & rho_freshwater=RHO_W, &
!          & rho_ice=RHO)
!
!  end if
!
!  hydro_topg       = transpose(st%zl)-transpose(st%z_sl)
!  hydro_temppabase = transpose(temph_b)
!
!  where (transpose(st%maske)==0)   ! grounded ice
!     hydro_icemask = 1
!     hydro_thk     = transpose(H_c+H_t)
!     hydro_supply  = rho_rho_w_ratio*transpose(Q_b_tot)
!  elsewhere
!     hydro_icemask = 0
!     hydro_thk     = 0.0_wp
!     hydro_supply  = 0.0_wp
!  end where
!
!  call hydro_set_topg(hydro, hydro_topg)
!  call hydro_set_thk(hydro, hydro_thk)
!  call hydro_set_temppabase(hydro, hydro_temppabase)
!  call hydro_set_supply(hydro, hydro_supply)
!  call hydro_set_mask(hydro, hydro_icemask)
!
!  call hydro_update(hydro)
!
!  call hydro_get_sflux(hydro, hydro_sflux)
!  call hydro_get_vflux(hydro, hydro_vflux, hydro_vfluxX, hydro_vfluxY)
!  call hydro_get_bwat(hydro, hydro_bwat)
!
!  st%q_w   = transpose(hydro_vflux)
!  st%q_w_x = transpose(hydro_vfluxX)
!  st%q_w_y = transpose(hydro_vfluxY)
!  st%H_w   = transpose(hydro_bwat)

else if (par%basal_hydrology==2) then 

    ! update basal water thickness
    where (st%maske==0) ! grounded ice
      st%H_w =  st%H_w + st%Q_b_tot*tmr%dtime - 0.001_wp/sec_year*tmr%dtime
      st%H_w = max(st%H_w,0._wp)
      st%H_w = min(st%H_w,par%H_w_max)
    endwhere

    st%q_w   = 0.0_wp
    st%q_w_x = 0.0_wp
    st%q_w_y = 0.0_wp

else
  where (st%maske==0) ! grounded ice
    st%q_w   = 0.0_wp
    st%q_w_x = 0.0_wp
    st%q_w_y = 0.0_wp
    st%H_w   = 0.0_wp   
  endwhere
endif

  where (st%maske==2)   ! ocean
     st%q_w   = 0.0_wp
     st%q_w_x = 0.0_wp
     st%q_w_y = 0.0_wp
     st%H_w   = st%z_sl-st%zl
     st%H_w = max(st%H_w,0._wp)
     st%H_w = min(st%H_w,par%H_w_max)
  elsewhere (st%maske==3)   ! floating ice
     st%q_w   = 0.0_wp
     st%q_w_x = 0.0_wp
     st%q_w_y = 0.0_wp
     st%H_w   = st%zb-st%zl
     st%H_w = max(st%H_w,0._wp)
     st%H_w = min(st%H_w,par%H_w_max)
  elsewhere (st%maske==1)   ! ice-free land
     st%q_w   = 0.0_wp
     st%q_w_x = 0.0_wp
     st%q_w_y = 0.0_wp
     st%H_w   = 0.0_wp
  end where

  if (firstcall) firstcall = .false.

  end subroutine calc_thk_water_bas

!-------------------------------------------------------------------------------

end module calc_thk_water_bas_m
!
