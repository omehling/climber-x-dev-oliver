!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  t o p o g r a d _ m
!
!> @file
!!
!! Calculation of topography gradients on the staggered grid and on the grid
!! points (including length rescaling with the corresponding components of the
!! metric tensor).
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
!> Calculation of topography gradients on the staggered grid and on the grid
!! points (including length rescaling with the corresponding components of the
!! metric tensor).
!<------------------------------------------------------------------------------
module topograd_m

  use sico_types_m
  use sico_state
  use sico_grid_mod
  use sico_timer

  implicit none

  public

contains

!-------------------------------------------------------------------------------
!> Calculation of topography gradients on the staggered grid and on the grid
!! points (the latter by second-order discretization).
!<------------------------------------------------------------------------------
  subroutine topograd_1(st,grd,n_switch)

  implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd

  integer, intent(in) :: n_switch

  integer                       :: i, j

!-------- Distinguish between old and new topography data --------

  if (n_switch == 1) then
     st%zs_aux = st%zs
     st%zm_aux = st%zm
     st%zb_aux = st%zb
  else if (n_switch == 2) then
     st%zs_aux = st%zs_neu
     st%zm_aux = st%zm_neu
     st%zb_aux = st%zb_neu
  else
     stop ' >>> topograd_1: Wrong value for n_switch!'
  end if

!-------- Topography gradients on the staggered grid --------

!  ------ x-derivatives

  do i=0, grd%IMAX-1
  do j=0, grd%JMAX
     st%dzs_dxi(j,i)  = (st%zs_aux(j,i+1)-st%zs_aux(j,i))*grd%dxi_inv &
                     *grd%insq_g11_sgx(j,i)
     st%dzm_dxi(j,i)  = (st%zm_aux(j,i+1)-st%zm_aux(j,i))*grd%dxi_inv &
                     *grd%insq_g11_sgx(j,i)
     st%dzb_dxi(j,i)  = (st%zb_aux(j,i+1)-st%zb_aux(j,i))*grd%dxi_inv &
                     *grd%insq_g11_sgx(j,i)
     st%dH_c_dxi(j,i) = st%dzs_dxi(j,i)-st%dzm_dxi(j,i)
     st%dH_t_dxi(j,i) = st%dzm_dxi(j,i)-st%dzb_dxi(j,i)
  end do
  end do

!  ------ y-derivatives

  do i=0, grd%IMAX
  do j=0, grd%JMAX-1
     st%dzs_deta(j,i)  = (st%zs_aux(j+1,i)-st%zs_aux(j,i))*grd%deta_inv &
                      *grd%insq_g22_sgy(j,i)
     st%dzm_deta(j,i)  = (st%zm_aux(j+1,i)-st%zm_aux(j,i))*grd%deta_inv &
                      *grd%insq_g22_sgy(j,i)
     st%dzb_deta(j,i)  = (st%zb_aux(j+1,i)-st%zb_aux(j,i))*grd%deta_inv &
                      *grd%insq_g22_sgy(j,i)
     st%dH_c_deta(j,i) = st%dzs_deta(j,i)-st%dzm_deta(j,i)
     st%dH_t_deta(j,i) = st%dzm_deta(j,i)-st%dzb_deta(j,i)
  end do
  end do

!-------- Topography gradients on the grid points --------

!  ------ x-derivatives

  do i=1, grd%IMAX-1
  do j=0, grd%JMAX
     st%dzs_dxi_g(j,i)  = (st%zs_aux(j,i+1)-st%zs_aux(j,i-1))*0.5_wp*grd%dxi_inv &
                       *grd%insq_g11_g(j,i)
     st%dzm_dxi_g(j,i)  = (st%zm_aux(j,i+1)-st%zm_aux(j,i-1))*0.5_wp*grd%dxi_inv &
                       *grd%insq_g11_g(j,i)
     st%dzb_dxi_g(j,i)  = (st%zb_aux(j,i+1)-st%zb_aux(j,i-1))*0.5_wp*grd%dxi_inv &
                       *grd%insq_g11_g(j,i)
     st%dH_c_dxi_g(j,i) = st%dzs_dxi_g(j,i)-st%dzm_dxi_g(j,i)
     st%dH_t_dxi_g(j,i) = st%dzm_dxi_g(j,i)-st%dzb_dxi_g(j,i)
  end do
  end do

  do j=0, grd%JMAX
     st%dzs_dxi_g(j,0)     = (st%zs_aux(j,1)-st%zs_aux(j,0))*grd%dxi_inv &
                          *grd%insq_g11_g(j,0)
     st%dzm_dxi_g(j,0)     = (st%zm_aux(j,1)-st%zm_aux(j,0))*grd%dxi_inv &
                          *grd%insq_g11_g(j,0)
     st%dzb_dxi_g(j,0)     = (st%zb_aux(j,1)-st%zb_aux(j,0))*grd%dxi_inv &
                          *grd%insq_g11_g(j,0)
     st%dH_c_dxi_g(j,0)    = st%dzs_dxi_g(j,0)-st%dzm_dxi_g(j,0)
     st%dH_t_dxi_g(j,0)    = st%dzm_dxi_g(j,0)-st%dzb_dxi_g(j,0)
     st%dzs_dxi_g(j,grd%IMAX)  = (st%zs_aux(j,grd%IMAX)-st%zs_aux(j,grd%IMAX-1)) &
                          *grd%dxi_inv &
                          *grd%insq_g11_g(j,grd%IMAX)
     st%dzm_dxi_g(j,grd%IMAX)  = (st%zm_aux(j,grd%IMAX)-st%zm_aux(j,grd%IMAX-1)) &
                          *grd%dxi_inv &
                          *grd%insq_g11_g(j,grd%IMAX)
     st%dzb_dxi_g(j,grd%IMAX)  = (st%zb_aux(j,grd%IMAX)-st%zb_aux(j,grd%IMAX-1)) &
                          *grd%dxi_inv &
                          *grd%insq_g11_g(j,grd%IMAX)
     st%dH_c_dxi_g(j,grd%IMAX) = st%dzs_dxi_g(j,grd%IMAX)-st%dzm_dxi_g(j,grd%IMAX)
     st%dH_t_dxi_g(j,grd%IMAX) = st%dzm_dxi_g(j,grd%IMAX)-st%dzb_dxi_g(j,grd%IMAX)
  end do

!  ------ y-derivatives

  do i=0, grd%IMAX
  do j=1, grd%JMAX-1
     st%dzs_deta_g(j,i)  = (st%zs_aux(j+1,i)-st%zs_aux(j-1,i)) &
                        *0.5_wp*grd%deta_inv &
                        *grd%insq_g22_g(j,i)
     st%dzm_deta_g(j,i)  = (st%zm_aux(j+1,i)-st%zm_aux(j-1,i)) &
                        *0.5_wp*grd%deta_inv &
                        *grd%insq_g22_g(j,i)
     st%dzb_deta_g(j,i)  = (st%zb_aux(j+1,i)-st%zb_aux(j-1,i)) &
                        *0.5_wp*grd%deta_inv &
                        *grd%insq_g22_g(j,i)
     st%dH_c_deta_g(j,i) = st%dzs_deta_g(j,i)-st%dzm_deta_g(j,i)
     st%dH_t_deta_g(j,i) = st%dzm_deta_g(j,i)-st%dzb_deta_g(j,i)
  end do
  end do

  do i=0, grd%IMAX
     st%dzs_deta_g(0,i)     = (st%zs_aux(1,i)-st%zs_aux(0,i))*grd%deta_inv &
                           *grd%insq_g22_g(0,i)
     st%dzm_deta_g(0,i)     = (st%zm_aux(1,i)-st%zm_aux(0,i))*grd%deta_inv &
                           *grd%insq_g22_g(0,i)
     st%dzb_deta_g(0,i)     = (st%zb_aux(1,i)-st%zb_aux(0,i))*grd%deta_inv &
                           *grd%insq_g22_g(0,i)
     st%dH_c_deta_g(0,i)    = st%dzs_deta_g(0,i)-st%dzm_deta_g(0,i)
     st%dH_t_deta_g(0,i)    = st%dzm_deta_g(0,i)-st%dzb_deta_g(0,i)
     st%dzs_deta_g(grd%JMAX,i)  = (st%zs_aux(grd%JMAX,i)-st%zs_aux(grd%JMAX-1,i)) &
                           *grd%deta_inv &
                           *grd%insq_g22_g(grd%JMAX,i)
     st%dzm_deta_g(grd%JMAX,i)  = (st%zm_aux(grd%JMAX,i)-st%zm_aux(grd%JMAX-1,i)) &
                           *grd%deta_inv &
                           *grd%insq_g22_g(grd%JMAX,i)
     st%dzb_deta_g(grd%JMAX,i)  = (st%zb_aux(grd%JMAX,i)-st%zb_aux(grd%JMAX-1,i)) &
                           *grd%deta_inv &
                           *grd%insq_g22_g(grd%JMAX,i)
     st%dH_c_deta_g(grd%JMAX,i) = st%dzs_deta_g(grd%JMAX,i)-st%dzm_deta_g(grd%JMAX,i)
     st%dH_t_deta_g(grd%JMAX,i) = st%dzm_deta_g(grd%JMAX,i)-st%dzb_deta_g(grd%JMAX,i)
  end do

  end subroutine topograd_1

!-------------------------------------------------------------------------------
!> Calculation of topography gradients on the staggered grid and on the grid
!! points (the latter by fourth-order discretization).
!<------------------------------------------------------------------------------
  subroutine topograd_2(st,grd,n_switch)

  implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd

  integer, intent(in) :: n_switch

  integer                       :: i, j

!-------- Distinguish between old and new topography data --------

  if (n_switch == 1) then
     st%zs_aux = st%zs
     st%zm_aux = st%zm
     st%zb_aux = st%zb
  else if (n_switch == 2) then
     st%zs_aux = st%zs_neu
     st%zm_aux = st%zm_neu
     st%zb_aux = st%zb_neu
  else
     stop ' >>> topograd_2: Wrong value for n_switch!'
  end if

!-------- Topography gradients on the staggered grid --------

!  ------ x-derivatives

  do i=0, grd%IMAX-1
  do j=0, grd%JMAX
     st%dzs_dxi(j,i)  = (st%zs_aux(j,i+1)-st%zs_aux(j,i))*grd%dxi_inv &
                     *grd%insq_g11_sgx(j,i)
     st%dzm_dxi(j,i)  = (st%zm_aux(j,i+1)-st%zm_aux(j,i))*grd%dxi_inv &
                     *grd%insq_g11_sgx(j,i)
     st%dzb_dxi(j,i)  = (st%zb_aux(j,i+1)-st%zb_aux(j,i))*grd%dxi_inv &
                     *grd%insq_g11_sgx(j,i)
     st%dH_c_dxi(j,i) = st%dzs_dxi(j,i)-st%dzm_dxi(j,i)
     st%dH_t_dxi(j,i) = st%dzm_dxi(j,i)-st%dzb_dxi(j,i)
  end do
  end do

!  ------ y-derivatives

  do i=0, grd%IMAX
  do j=0, grd%JMAX-1
     st%dzs_deta(j,i)  = (st%zs_aux(j+1,i)-st%zs_aux(j,i))*grd%deta_inv &
                      *grd%insq_g22_sgy(j,i)
     st%dzm_deta(j,i)  = (st%zm_aux(j+1,i)-st%zm_aux(j,i))*grd%deta_inv &
                      *grd%insq_g22_sgy(j,i)
     st%dzb_deta(j,i)  = (st%zb_aux(j+1,i)-st%zb_aux(j,i))*grd%deta_inv &
                      *grd%insq_g22_sgy(j,i)
     st%dH_c_deta(j,i) = st%dzs_deta(j,i)-st%dzm_deta(j,i)
     st%dH_t_deta(j,i) = st%dzm_deta(j,i)-st%dzb_deta(j,i)
  end do
  end do

!-------- Topography gradients on the grid points --------

!  ------ x-derivatives

  do i=2, grd%IMAX-2
  do j=0, grd%JMAX
     st%dzs_dxi_g(j,i) &
        = (        -st%zs_aux(j,i+2) + 8.0_wp*st%zs_aux(j,i+1) &
            -8.0_wp*st%zs_aux(j,i-1) +        st%zs_aux(j,i-2) ) &
          *grd%dxi12_inv &
          *grd%insq_g11_g(j,i)
     st%dzm_dxi_g(j,i) &
        = (      -st%zm_aux(j,i+2)   + 8.0_wp*st%zm_aux(j,i+1) &
            -8.0_wp*st%zm_aux(j,i-1) +        st%zm_aux(j,i-2) ) &
          *grd%dxi12_inv &
          *grd%insq_g11_g(j,i)
     st%dzb_dxi_g(j,i) &
        = (     -st%zb_aux(j,i+2)    + 8.0_wp*st%zb_aux(j,i+1) &
            -8.0_wp*st%zb_aux(j,i-1) +        st%zb_aux(j,i-2) ) &
          *grd%dxi12_inv &
          *grd%insq_g11_g(j,i)
     st%dH_c_dxi_g(j,i) = st%dzs_dxi_g(j,i)-st%dzm_dxi_g(j,i)
     st%dH_t_dxi_g(j,i) = st%dzm_dxi_g(j,i)-st%dzb_dxi_g(j,i)
  end do
  end do

  do j=0, grd%JMAX

     st%dzs_dxi_g(j,0)       = (st%zs_aux(j,1)-st%zs_aux(j,0))*grd%dxi_inv &
                            *grd%insq_g11_g(j,0)
     st%dzm_dxi_g(j,0)       = (st%zm_aux(j,1)-st%zm_aux(j,0))*grd%dxi_inv &
                            *grd%insq_g11_g(j,0)
     st%dzb_dxi_g(j,0)       = (st%zb_aux(j,1)-st%zb_aux(j,0))*grd%dxi_inv &
                            *grd%insq_g11_g(j,0)
     st%dH_c_dxi_g(j,0)      = st%dzs_dxi_g(j,0)-st%dzm_dxi_g(j,0)
     st%dH_t_dxi_g(j,0)      = st%dzm_dxi_g(j,0)-st%dzb_dxi_g(j,0)

     st%dzs_dxi_g(j,1)       = (st%zs_aux(j,2)-st%zs_aux(j,0)) &
                            *0.5_wp*grd%dxi_inv &
                            *grd%insq_g11_g(j,1)
     st%dzm_dxi_g(j,1)       = (st%zm_aux(j,2)-st%zm_aux(j,0)) &
                            *0.5_wp*grd%dxi_inv &
                            *grd%insq_g11_g(j,1)
     st%dzb_dxi_g(j,1)       = (st%zb_aux(j,2)-st%zb_aux(j,0)) &
                            *0.5_wp*grd%dxi_inv &
                            *grd%insq_g11_g(j,1)
     st%dH_c_dxi_g(j,1)      = st%dzs_dxi_g(j,1)-st%dzm_dxi_g(j,1)
     st%dH_t_dxi_g(j,1)      = st%dzm_dxi_g(j,1)-st%dzb_dxi_g(j,1)

     st%dzs_dxi_g(j,grd%IMAX-1)  = (st%zs_aux(j,grd%IMAX)-st%zs_aux(j,grd%IMAX-2)) &
                            *0.5_wp*grd%dxi_inv &
                            *grd%insq_g11_g(j,grd%IMAX-1)
     st%dzm_dxi_g(j,grd%IMAX-1)  = (st%zm_aux(j,grd%IMAX)-st%zm_aux(j,grd%IMAX-2)) &
                            *0.5_wp*grd%dxi_inv &
                            *grd%insq_g11_g(j,grd%IMAX-1)
     st%dzb_dxi_g(j,grd%IMAX-1)  = (st%zb_aux(j,grd%IMAX)-st%zb_aux(j,grd%IMAX-2)) &
                            *0.5_wp*grd%dxi_inv &
                            *grd%insq_g11_g(j,grd%IMAX-1)
     st%dH_c_dxi_g(j,grd%IMAX-1) = st%dzs_dxi_g(j,grd%IMAX-1) &
                            -st%dzm_dxi_g(j,grd%IMAX-1)
     st%dH_t_dxi_g(j,grd%IMAX-1) = st%dzm_dxi_g(j,grd%IMAX-1) &
                            -st%dzb_dxi_g(j,grd%IMAX-1)

     st%dzs_dxi_g(j,grd%IMAX)    = (st%zs_aux(j,grd%IMAX)-st%zs_aux(j,grd%IMAX-1)) &
                            *grd%dxi_inv &
                            *grd%insq_g11_g(j,grd%IMAX)
     st%dzm_dxi_g(j,grd%IMAX)    = (st%zm_aux(j,grd%IMAX)-st%zm_aux(j,grd%IMAX-1)) &
                            *grd%dxi_inv &
                            *grd%insq_g11_g(j,grd%IMAX)
     st%dzb_dxi_g(j,grd%IMAX)    = (st%zb_aux(j,grd%IMAX)-st%zb_aux(j,grd%IMAX-1)) &
                            *grd%dxi_inv &
                            *grd%insq_g11_g(j,grd%IMAX)
     st%dH_c_dxi_g(j,grd%IMAX)   = st%dzs_dxi_g(j,grd%IMAX)-st%dzm_dxi_g(j,grd%IMAX)
     st%dH_t_dxi_g(j,grd%IMAX)   = st%dzm_dxi_g(j,grd%IMAX)-st%dzb_dxi_g(j,grd%IMAX)

  end do

!  ------ y-derivatives

  do i=0, grd%IMAX
  do j=2, grd%JMAX-2
     st%dzs_deta_g(j,i) &
        = (     -st%zs_aux(j+2,i)    + 8.0_wp*st%zs_aux(j+1,i) &
            -8.0_wp*st%zs_aux(j-1,i) +        st%zs_aux(j-2,i) ) &
          *grd%deta12_inv &
          *grd%insq_g22_g(j,i)
     st%dzm_deta_g(j,i) &
        = (     -st%zm_aux(j+2,i)    + 8.0_wp*st%zm_aux(j+1,i) &
            -8.0_wp*st%zm_aux(j-1,i) +        st%zm_aux(j-2,i) ) &
          *grd%deta12_inv &
          *grd%insq_g22_g(j,i)
     st%dzb_deta_g(j,i) &
        = (     -st%zb_aux(j+2,i)    + 8.0_wp*st%zb_aux(j+1,i) &
            -8.0_wp*st%zb_aux(j-1,i) +        st%zb_aux(j-2,i) ) &
          *grd%deta12_inv &
          *grd%insq_g22_g(j,i)
     st%dH_c_deta_g(j,i) = st%dzs_deta_g(j,i)-st%dzm_deta_g(j,i)
     st%dH_t_deta_g(j,i) = st%dzm_deta_g(j,i)-st%dzb_deta_g(j,i)
  end do
  end do

  do i=0, grd%IMAX

     st%dzs_deta_g(0,i)       = (st%zs_aux(1,i)-st%zs_aux(0,i))*grd%deta_inv &
                             *grd%insq_g22_g(0,i)
     st%dzm_deta_g(0,i)       = (st%zm_aux(1,i)-st%zm_aux(0,i))*grd%deta_inv &
                             *grd%insq_g22_g(0,i)
     st%dzb_deta_g(0,i)       = (st%zb_aux(1,i)-st%zb_aux(0,i))*grd%deta_inv &
                             *grd%insq_g22_g(0,i)
     st%dH_c_deta_g(0,i)      = st%dzs_deta_g(0,i)-st%dzm_deta_g(0,i)
     st%dH_t_deta_g(0,i)      = st%dzm_deta_g(0,i)-st%dzb_deta_g(0,i)

     st%dzs_deta_g(1,i)       = (st%zs_aux(2,i)-st%zs_aux(0,i)) &
                             *0.5_wp*grd%deta_inv &
                             *grd%insq_g22_g(1,i)
     st%dzm_deta_g(1,i)       = (st%zm_aux(2,i)-st%zm_aux(0,i)) &
                             *0.5_wp*grd%deta_inv &
                             *grd%insq_g22_g(1,i)
     st%dzb_deta_g(1,i)       = (st%zb_aux(2,i)-st%zb_aux(0,i)) &
                             *0.5_wp*grd%deta_inv &
                             *grd%insq_g22_g(1,i)
     st%dH_c_deta_g(1,i)      = st%dzs_deta_g(1,i)-st%dzm_deta_g(1,i)
     st%dH_t_deta_g(1,i)      = st%dzm_deta_g(1,i)-st%dzb_deta_g(1,i)

     st%dzs_deta_g(grd%JMAX-1,i)  = (st%zs_aux(grd%JMAX,i)-st%zs_aux(grd%JMAX-2,i)) &
                             *0.5_wp*grd%deta_inv &
                             *grd%insq_g22_g(grd%JMAX-1,i)
     st%dzm_deta_g(grd%JMAX-1,i)  = (st%zm_aux(grd%JMAX,i)-st%zm_aux(grd%JMAX-2,i)) &
                             *0.5_wp*grd%deta_inv &
                             *grd%insq_g22_g(grd%JMAX-1,i)
     st%dzb_deta_g(grd%JMAX-1,i)  = (st%zb_aux(grd%JMAX,i)-st%zb_aux(grd%JMAX-2,i)) &
                             *0.5_wp*grd%deta_inv &
                             *grd%insq_g22_g(grd%JMAX-1,i)
     st%dH_c_deta_g(grd%JMAX-1,i) = st%dzs_deta_g(grd%JMAX-1,i) &
                             -st%dzm_deta_g(grd%JMAX-1,i)
     st%dH_t_deta_g(grd%JMAX-1,i) = st%dzm_deta_g(grd%JMAX-1,i) &
                             -st%dzb_deta_g(grd%JMAX-1,i)

     st%dzs_deta_g(grd%JMAX,i)    = (st%zs_aux(grd%JMAX,i)-st%zs_aux(grd%JMAX-1,i)) &
                             *grd%deta_inv &
                             *grd%insq_g22_g(grd%JMAX,i)
     st%dzm_deta_g(grd%JMAX,i)    = (st%zm_aux(grd%JMAX,i)-st%zm_aux(grd%JMAX-1,i)) &
                             *grd%deta_inv &
                             *grd%insq_g22_g(grd%JMAX,i)
     st%dzb_deta_g(grd%JMAX,i)    = (st%zb_aux(grd%JMAX,i)-st%zb_aux(grd%JMAX-1,i)) &
                             *grd%deta_inv &
                             *grd%insq_g22_g(grd%JMAX,i)
     st%dH_c_deta_g(grd%JMAX,i)   = st%dzs_deta_g(grd%JMAX,i)-st%dzm_deta_g(grd%JMAX,i)
     st%dH_t_deta_g(grd%JMAX,i)   = st%dzm_deta_g(grd%JMAX,i)-st%dzb_deta_g(grd%JMAX,i)

  end do

  end subroutine topograd_2

!-------------------------------------------------------------------------------

end module topograd_m

