!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ v z _ m
!
!> @file
!!
!! Computation of the vertical velocity vz.
!!
!! @section Copyright
!!
!! Copyright 2009-2017 Ralf Greve, Tatsuru Sato
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
!> Computation of the vertical velocity vz.
!<------------------------------------------------------------------------------
module calc_vz_m

  use sico_types_m
  use sico_state
  use sico_grid_mod

  implicit none

  private
  public :: calc_vz_grounded, calc_vz_floating, calc_vz_static

contains

!-------------------------------------------------------------------------------
!> Computation of the vertical velocity vz for grounded ice.
!<------------------------------------------------------------------------------
subroutine calc_vz_grounded(st,grd)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd

integer :: i, j, kc, kt
real(wp)     :: cvz00, cvz0(0:100), cvz1(0:100), cvz2(0:100), &
                cvz3(0:100), cvz4(0:100), cvz5(0:100)


do i=1, grd%IMAX-1
do j=1, grd%JMAX-1

   if (st%maske(j,i)==0) then   ! grounded ice

!-------- Abbreviations --------

      cvz00 = 0.5_wp*(st%vx_t(0,j,i)+st%vx_t(0,j,i-1))*st%dzb_dxi_g(j,i) &
              +0.5_wp*(st%vy_t(0,j,i)+st%vy_t(0,j-1,i))*st%dzb_deta_g(j,i) &
              +st%dzb_dtau(j,i)-st%Q_b_tot(j,i)

      kt=0
      cvz0(kt) = st%H_t(j,i)*grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i) &
              *( (st%vx_t(kt,j,i)*grd%sq_g22_sgx(j,i) &
                 -st%vx_t(kt,j,i-1)*grd%sq_g22_sgx(j,i-1))*grd%dxi_inv &
                +(st%vy_t(kt,j,i)*grd%sq_g11_sgy(j,i) &
                 -st%vy_t(kt,j-1,i)*grd%sq_g11_sgy(j-1,i))*grd%deta_inv ) &
              *grd%dzeta_t
      cvz1(kt) = (st%dzb_dxi_g(j,i)+grd%zeta_t(kt)*st%dH_t_dxi_g(j,i)) &
              *(st%vx_t(kt+1,j,i)+st%vx_t(kt+1,j,i-1) &
               -st%vx_t(kt,j,i)  -st%vx_t(kt,j,i-1))*0.5_wp
      cvz2(kt) = (st%dzb_deta_g(j,i)+grd%zeta_t(kt)*st%dH_t_deta_g(j,i)) &
              *(st%vy_t(kt+1,j,i)+st%vy_t(kt+1,j-1,i) &
               -st%vy_t(kt,j,i)  -st%vy_t(kt,j-1,i))*0.5_wp

      do kt=1, grd%KTMAX-1
         cvz0(kt) = st%H_t(j,i)*grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i) &
              *( (st%vx_t(kt,j,i)*grd%sq_g22_sgx(j,i) &
                 -st%vx_t(kt,j,i-1)*grd%sq_g22_sgx(j,i-1))*grd%dxi_inv &
                +(st%vy_t(kt,j,i)*grd%sq_g11_sgy(j,i) &
                 -st%vy_t(kt,j-1,i)*grd%sq_g11_sgy(j-1,i))*grd%deta_inv ) &
              *grd%dzeta_t
         cvz1(kt) = (st%dzb_dxi_g(j,i)+grd%zeta_t(kt)*st%dH_t_dxi_g(j,i)) &
              *(st%vx_t(kt+1,j,i)+st%vx_t(kt+1,j,i-1) &
               -st%vx_t(kt-1,j,i)-st%vx_t(kt-1,j,i-1))*0.25_wp
         cvz2(kt) = (st%dzb_deta_g(j,i)+grd%zeta_t(kt)*st%dH_t_deta_g(j,i)) &
              *(st%vy_t(kt+1,j,i)+st%vy_t(kt+1,j-1,i) &
               -st%vy_t(kt-1,j,i)-st%vy_t(kt-1,j-1,i))*0.25_wp
      end do

      kt=grd%KTMAX
      cvz0(kt) = st%H_t(j,i)*grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i) &
              *( (st%vx_t(kt,j,i)*grd%sq_g22_sgx(j,i) &
                 -st%vx_t(kt,j,i-1)*grd%sq_g22_sgx(j,i-1))*grd%dxi_inv &
                +(st%vy_t(kt,j,i)*grd%sq_g11_sgy(j,i) &
                 -st%vy_t(kt,j-1,i)*grd%sq_g11_sgy(j-1,i))*grd%deta_inv ) &
              *grd%dzeta_t
      cvz1(kt) = (st%dzb_dxi_g(j,i)+grd%zeta_t(kt)*st%dH_t_dxi_g(j,i)) &
              *(st%vx_t(kt,j,i)  +st%vx_t(kt,j,i-1) &
               -st%vx_t(kt-1,j,i)-st%vx_t(kt-1,j,i-1))*0.5_wp
      cvz2(kt) = (st%dzb_deta_g(j,i)+grd%zeta_t(kt)*st%dH_t_deta_g(j,i)) &
              *(st%vy_t(kt,j,i)  +st%vy_t(kt,j-1,i) &
               -st%vy_t(kt-1,j,i)-st%vy_t(kt-1,j-1,i))*0.5_wp

      kc=0
      cvz3(kc) = grd%avz3(kc)*st%H_c(j,i) &
                 *grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i) &
              *( (st%vx_c(kc,j,i)*grd%sq_g22_sgx(j,i) &
                 -st%vx_c(kc,j,i-1)*grd%sq_g22_sgx(j,i-1))*grd%dxi_inv &
                +(st%vy_c(kc,j,i)*grd%sq_g11_sgy(j,i) &
                 -st%vy_c(kc,j-1,i)*grd%sq_g11_sgy(j-1,i))*grd%deta_inv )
      cvz4(kc) = (st%dzm_dxi_g(j,i) &
              +grd%eaz_c_quotient(kc)*st%dH_c_dxi_g(j,i)) &
              *(st%vx_c(kc+1,j,i)+st%vx_c(kc+1,j,i-1) &
               -st%vx_c(kc,j,i)  -st%vx_c(kc,j,i-1))*0.5_wp
      cvz5(kc) = (st%dzm_deta_g(j,i) &
              +grd%eaz_c_quotient(kc)*st%dH_c_deta_g(j,i)) &
              *(st%vy_c(kc+1,j,i)+st%vy_c(kc+1,j-1,i) &
               -st%vy_c(kc,j,i)  -st%vy_c(kc,j-1,i))*0.5_wp

      do kc=1, grd%KCMAX-1
         cvz3(kc) = grd%avz3(kc)*st%H_c(j,i) &
                    *grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i) &
              *( (st%vx_c(kc,j,i)*grd%sq_g22_sgx(j,i) &
                 -st%vx_c(kc,j,i-1)*grd%sq_g22_sgx(j,i-1))*grd%dxi_inv &
                +(st%vy_c(kc,j,i)*grd%sq_g11_sgy(j,i) &
                 -st%vy_c(kc,j-1,i)*grd%sq_g11_sgy(j-1,i))*grd%deta_inv )
         cvz4(kc) = (st%dzm_dxi_g(j,i) &
              +grd%eaz_c_quotient(kc)*st%dH_c_dxi_g(j,i)) &
              *(st%vx_c(kc+1,j,i)+st%vx_c(kc+1,j,i-1) &
               -st%vx_c(kc-1,j,i)-st%vx_c(kc-1,j,i-1))*0.25_wp
         cvz5(kc) = (st%dzm_deta_g(j,i) &
              +grd%eaz_c_quotient(kc)*st%dH_c_deta_g(j,i)) &
              *(st%vy_c(kc+1,j,i)+st%vy_c(kc+1,j-1,i) &
               -st%vy_c(kc-1,j,i)-st%vy_c(kc-1,j-1,i))*0.25_wp
      end do

      kc=grd%KCMAX
      cvz3(kc) = grd%avz3(kc)*st%H_c(j,i) &
                 *grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i) &
              *( (st%vx_c(kc,j,i)*grd%sq_g22_sgx(j,i) &
                 -st%vx_c(kc,j,i-1)*grd%sq_g22_sgx(j,i-1))*grd%dxi_inv &
                +(st%vy_c(kc,j,i)*grd%sq_g11_sgy(j,i) &
                 -st%vy_c(kc,j-1,i)*grd%sq_g11_sgy(j-1,i))*grd%deta_inv )
      cvz4(kc) = (st%dzm_dxi_g(j,i) &
              +grd%eaz_c_quotient(kc)*st%dH_c_dxi_g(j,i)) &
              *(st%vx_c(kc,j,i)  +st%vx_c(kc,j,i-1) &
               -st%vx_c(kc-1,j,i)-st%vx_c(kc-1,j,i-1))*0.5_wp
      cvz5(kc) = (st%dzm_deta_g(j,i) &
              +grd%eaz_c_quotient(kc)*st%dH_c_deta_g(j,i)) &
              *(st%vy_c(kc,j,i)  +st%vy_c(kc,j-1,i) &
               -st%vy_c(kc-1,j,i)-st%vy_c(kc-1,j-1,i))*0.5_wp

!-------- Computation of vz_b --------

      st%vz_b(j,i) = cvz00

!-------- Computation of vz --------

      if ((st%n_cts(j,i) == -1).or.(st%n_cts(j,i) == 0)) then
                        ! cold ice base, temperate ice base

         do kt=0, grd%KTMAX-1
            st%vz_t(kt,j,i) = st%vz_b(j,i)
         end do

         st%vz_m(j,i) = st%vz_b(j,i)

         st%vz_c(0,j,i) = st%vz_m(j,i) &
                       +0.5_wp*(-cvz3(0)+cvz4(0)+cvz5(0))

         do kc=1, grd%KCMAX-1
            st%vz_c(kc,j,i) = st%vz_c(kc-1,j,i) &
                           -cvz3(kc)+cvz4(kc)+cvz5(kc)
         end do

      else   ! st%n_cts(j,i) == 1, temperate ice layer

         st%vz_t(0,j,i) = st%vz_b(j,i) &
                       +0.5_wp*(-cvz0(0)+cvz1(0)+cvz2(0))

         do kt=1, grd%KTMAX-1
            st%vz_t(kt,j,i) = st%vz_t(kt-1,j,i) &
                           -cvz0(kt)+cvz1(kt)+cvz2(kt)
         end do

         st%vz_m(j,i) = st%vz_t(grd%KTMAX-1,j,i) &
                     +0.5_wp*(-cvz0(grd%KTMAX)+cvz1(grd%KTMAX)+cvz2(grd%KTMAX))

         st%vz_c(0,j,i) = st%vz_m(j,i) &
                       +0.5_wp*(-cvz3(0)+cvz4(0)+cvz5(0))

         do kc=1, grd%KCMAX-1
            st%vz_c(kc,j,i) = st%vz_c(kc-1,j,i) &
                           -cvz3(kc)+cvz4(kc)+cvz5(kc)
         end do

      end if

!-------- Computation of vz_s --------

      kc=grd%KCMAX

      st%vz_s(j,i) = st%vz_c(kc-1,j,i) &
                  +0.5_wp*(-cvz3(kc)+cvz4(kc)+cvz5(kc))

   else   ! maske(j,i) /= 0 (not grounded ice)

      st%vz_b(j,i) = 0.0_wp

      do kt=0, grd%KTMAX-1
         st%vz_t(kt,j,i) = 0.0_wp
      end do

      st%vz_m(j,i) = 0.0_wp

      do kc=0, grd%KCMAX-1
         st%vz_c(kc,j,i) = 0.0_wp
      end do

      st%vz_s(j,i) = 0.0_wp

   end if

end do
end do

end subroutine calc_vz_grounded

!-------------------------------------------------------------------------------
!> Computation of the vertical velocity vz for floating ice.
!<------------------------------------------------------------------------------
subroutine calc_vz_floating(st,grd)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd

integer :: i, j, kt, kc
real(wp) :: dvx_dxi, dvy_deta


!-------- Computation of vz --------

do i=1, grd%IMAX-1
do j=1, grd%JMAX-1

   if (st%maske(j,i)==3) then   ! floating ice

!  ------ Derivatives of the horizontal velocity

      dvx_dxi =  grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i) &
                 *(st%vx_m(j,i)*grd%sq_g22_sgx(j,i)-st%vx_m(j,i-1)*grd%sq_g22_sgx(j,i-1)) &
                 *grd%dxi_inv

      dvy_deta = grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i) &
                 *(st%vy_m(j,i)*grd%sq_g11_sgy(j,i)-st%vy_m(j-1,i)*grd%sq_g11_sgy(j-1,i)) &
                 *grd%deta_inv

!  ------ Basal velocity vz_b

      st%vz_b(j,i) = 0.5_wp*(st%vx_m(j,i)+st%vx_m(j,i-1))*st%dzb_dxi_g(j,i) &
                  +0.5_wp*(st%vy_m(j,i)+st%vy_m(j-1,i))*st%dzb_deta_g(j,i) &
                  +st%dzb_dtau(j,i)-st%Q_b_tot(j,i)
                  ! kinematic boundary condition at the ice base

!  ------ Velocity at sea level vz_sl

      st%vz_sl(j,i) = st%vz_b(j,i) - (st%z_sl(j,i)-st%zb(j,i))*(dvx_dxi+dvy_deta)

!  ------ Surface velocity st%vz_s

      st%vz_s(j,i) = st%vz_sl(j,i) - (st%zs(j,i)-st%z_sl(j,i))*(dvx_dxi+dvy_deta)

!  ------ Velocity vz_m at the interface between
!                              the upper (kc) and the lower (kt) domain

      if ((st%n_cts(j,i) == -1).or.(st%n_cts(j,i) == 0)) then
                        ! cold ice base, temperate ice base

         st%vz_m(j,i) = st%vz_b(j,i)

      else   ! n_cts(j,i) == 1, temperate ice layer

         st%vz_m(j,i) = st%vz_b(j,i) - st%H_t(j,i)*(dvx_dxi+dvy_deta)

      end if

!  ------ 3-D velocity vz_c and vz_t

      do kc=0, grd%KCMAX-1
         st%vz_c(kc,j,i) = st%vz_sl(j,i) &
                        -(st%zm(j,i)+grd%eaz_c_quotient_sgz(kc)*st%H_c(j,i)-st%z_sl(j,i)) &
                         *(dvx_dxi+dvy_deta)
      end do

      if ((st%n_cts(j,i) == -1).or.(st%n_cts(j,i) == 0)) then
                        ! cold ice base, temperate ice base

         do kt=0, grd%KTMAX-1
            st%vz_t(kt,j,i) = st%vz_b(j,i)
         end do

      else   ! n_cts(j,i) == 1, temperate ice layer

         do kt=0, grd%KTMAX-1
            st%vz_t(kt,j,i) = st%vz_sl(j,i) &
                           -(st%zb(j,i) &
                             +0.5_wp*(grd%zeta_t(kt)+grd%zeta_t(kt+1))*st%H_t(j,i) &
                             -st%z_sl(j,i)) &
                            *(dvx_dxi+dvy_deta)
         end do

      end if

   end if

end do
end do

end subroutine calc_vz_floating

!-------------------------------------------------------------------------------
!> Computation of the vertical velocity vz for static ice.
!<------------------------------------------------------------------------------
subroutine calc_vz_static(st,grd)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd

integer :: i, j, kc, kt

do i=0, grd%IMAX
do j=0, grd%JMAX

   if ((st%maske(j,i)==0).or.(st%maske(j,i)==3)) then   ! grounded or floating ice

      st%vz_b(j,i) = st%dzb_dtau(j,i)-st%Q_b_tot(j,i)
                  ! kinematic boundary condition at the ice base

      do kt=0, grd%KTMAX-1
         st%vz_t(kt,j,i) = st%vz_b(j,i)
      end do

      st%vz_m(j,i) = st%vz_b(j,i)

      do kc=0, grd%KCMAX-1
         st%vz_c(kc,j,i) = st%vz_b(j,i)
      end do

      st%vz_s(j,i) = st%vz_b(j,i)

   else   ! maske(j,i) == (1 or 2)

      st%vz_b(j,i) = 0.0_wp

      do kt=0, grd%KTMAX-1
         st%vz_t(kt,j,i) = 0.0_wp
      end do

      st%vz_m(j,i) = 0.0_wp

      do kc=0, grd%KCMAX-1
         st%vz_c(kc,j,i) = 0.0_wp
      end do

      st%vz_s(j,i) = 0.0_wp

   end if

end do
end do

end subroutine calc_vz_static

!-------------------------------------------------------------------------------

end module calc_vz_m
!
