!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ d x y z _ m
!
!> @file
!!
!! Computation of all components of the strain-rate tensor, the full
!! effective strain rate and the shear fraction.
!!
!! @section Copyright
!!
!! Copyright 2014-2017 Ralf Greve
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
!> Computation of all components of the strain-rate tensor, the full
!! effective strain rate and the shear fraction.
!<------------------------------------------------------------------------------
module calc_dxyz_m

  use sico_types_m
  use sico_state
  use sico_grid_mod
  use sico_params, only : sico_par_class, eps_dp

  implicit none

  private
  public :: calc_dxyz

contains

!-------------------------------------------------------------------------------
!> Main subroutine of calc_dxyz_m:
!! Computation of all components of the strain-rate tensor, the full
!! effective strain rate and the shear fraction.
!<------------------------------------------------------------------------------
  subroutine calc_dxyz(st,grd,par)

  implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class), intent(in) :: grd
  type(sico_par_class), intent(in) :: par

  integer                       :: i, j, kc, kt
  real(wp)                           :: H_c_inv, H_t_inv

  real(wp), dimension(:), allocatable :: dxx_c !> dxx_c(kc): Strain rate dxx in the upper (kc) ice domain
  real(wp), dimension(:), allocatable :: dyy_c !> dyy_c(kc): Strain rate dyy in the upper (kc) ice domain
  real(wp), dimension(:), allocatable :: dxy_c !> dxy_c(kc): Strain rate dxy in the upper (kc) ice domain
  real(wp), dimension(:), allocatable :: dxz_c !> dxz_c(kc): Strain rate dxz in the upper (kc) ice domain
  real(wp), dimension(:), allocatable :: dyz_c !> dyz_c(kc): Strain rate dyz in the upper (kc) ice domain

  real(wp), dimension(:), allocatable :: dxx_t !> dxx_t(kt): Strain rate dxx in the lower (kt) ice domain
  real(wp), dimension(:), allocatable :: dyy_t !> dyy_t(kt): Strain rate dyy in the lower (kt) ice domain
  real(wp), dimension(:), allocatable :: dxy_t !> dxy_t(kt): Strain rate dxy in the lower (kt) ice domain
  real(wp), dimension(:), allocatable :: dxz_t !> dxz_t(kt): Strain rate dxz in the lower (kt) ice domain
  real(wp), dimension(:), allocatable :: dyz_t !> dyz_t(kt): Strain rate dyz in the lower (kt) ice domain

  real(wp), dimension(:), allocatable         :: lxy_c, lyx_c, lxz_c, lzx_c, lyz_c, lzy_c
  real(wp), dimension(:), allocatable         :: lxy_t, lyx_t, lxz_t, lzx_t, lyz_t, lzy_t
  real(wp), dimension(:), allocatable         :: shear_c_squared
  real(wp), dimension(:), allocatable         :: shear_t_squared
  real(wp)                           :: abs_v_ssa_inv, nx, ny
  real(wp)                           :: shear_x_help, shear_y_help
  real(wp)                           :: lambda_shear_help


  allocate(dxx_c(0:grd%KCMAX))
  allocate(dyy_c(0:grd%KCMAX))
  allocate(dxy_c(0:grd%KCMAX))
  allocate(dxz_c(0:grd%KCMAX))
  allocate(dyz_c(0:grd%KCMAX))

  allocate(dxx_t(0:grd%KTMAX)) 
  allocate(dyy_t(0:grd%KTMAX)) 
  allocate(dxy_t(0:grd%KTMAX)) 
  allocate(dxz_t(0:grd%KTMAX)) 
  allocate(dyz_t(0:grd%KTMAX)) 

  allocate(lxy_c(0:grd%KCMAX))
  allocate(lyx_c(0:grd%KCMAX))
  allocate(lxz_c(0:grd%KCMAX))
  allocate(lzx_c(0:grd%KCMAX))
  allocate(lyz_c(0:grd%KCMAX))
  allocate(lzy_c(0:grd%KCMAX))
  
  allocate(lxy_t(0:grd%KTMAX))
  allocate(lyx_t(0:grd%KTMAX))
  allocate(lxz_t(0:grd%KTMAX))
  allocate(lzx_t(0:grd%KTMAX))
  allocate(lyz_t(0:grd%KTMAX))
  allocate(lzy_t(0:grd%KTMAX))

  allocate(shear_c_squared(0:grd%KCMAX))
  allocate(shear_t_squared(0:grd%KTMAX))

!-------- Initialisation --------

  dxx_c          = 0.0_wp
  dyy_c          = 0.0_wp
  dxy_c          = 0.0_wp
  dxz_c          = 0.0_wp
  dyz_c          = 0.0_wp
  st%de_c           = 0.0_wp
  st%lambda_shear_c = 0.0_wp

  dxx_t          = 0.0_wp
  dyy_t          = 0.0_wp
  dxy_t          = 0.0_wp
  dxz_t          = 0.0_wp
  dyz_t          = 0.0_wp
  st%de_t           = 0.0_wp
  st%lambda_shear_t = 0.0_wp

!-------- Computation --------

  do i=1, grd%IMAX-1
  do j=1, grd%JMAX-1

     if ((st%maske(j,i) == 0).or.(st%maske(j,i) == 3)) then
                                                 ! grounded or floating ice

        H_c_inv = 1.0_wp/(abs(st%H_c(j,i))+eps_dp)

        kc=0

           dxx_c(kc) = (st%vx_c(kc,j,i)-st%vx_c(kc,j,i-1))*grd%fact_x(j,i)
           dyy_c(kc) = (st%vy_c(kc,j,i)-st%vy_c(kc,j-1,i))*grd%fact_y(j,i)

           lxy_c(kc) = (  (st%vx_c(kc,j+1,i)+st%vx_c(kc,j+1,i-1)) &
                        - (st%vx_c(kc,j-1,i)+st%vx_c(kc,j-1,i-1)) ) &
                       *0.25_wp*grd%fact_y(j,i)

           lyx_c(kc) = (  (st%vy_c(kc,j,i+1)+st%vy_c(kc,j-1,i+1)) &
                        - (st%vy_c(kc,j,i-1)+st%vy_c(kc,j-1,i-1)) ) &
                       *0.25_wp*grd%fact_x(j,i)

           lzx_c(kc) = (st%vz_m(j,i+1)-st%vz_m(j,i-1))*0.5_wp*grd%fact_x(j,i)
           lzy_c(kc) = (st%vz_m(j+1,i)-st%vz_m(j-1,i))*0.5_wp*grd%fact_y(j,i)

           lxz_c(kc) = (  (st%vx_c(kc+1,j,i)+st%vx_c(kc+1,j,i-1)) &
                        - (st%vx_c(kc  ,j,i)+st%vx_c(kc  ,j,i-1)) ) &
                       *0.5_wp*grd%fact_z_c(kc)*H_c_inv

           lyz_c(kc) = (  (st%vy_c(kc+1,j,i)+st%vy_c(kc+1,j-1,i)) &
                        - (st%vy_c(kc  ,j,i)+st%vy_c(kc  ,j-1,i)) ) &
                       *0.5_wp*grd%fact_z_c(kc)*H_c_inv

        ! end kc=0

        do kc=1, grd%KCMAX-1

           dxx_c(kc) = (st%vx_c(kc,j,i)-st%vx_c(kc,j,i-1))*grd%fact_x(j,i)
           dyy_c(kc) = (st%vy_c(kc,j,i)-st%vy_c(kc,j-1,i))*grd%fact_y(j,i)

           lxy_c(kc) = (  (st%vx_c(kc,j+1,i)+st%vx_c(kc,j+1,i-1)) &
                        - (st%vx_c(kc,j-1,i)+st%vx_c(kc,j-1,i-1)) ) &
                       *0.25_wp*grd%fact_y(j,i)

           lyx_c(kc) = (  (st%vy_c(kc,j,i+1)+st%vy_c(kc,j-1,i+1)) &
                        - (st%vy_c(kc,j,i-1)+st%vy_c(kc,j-1,i-1)) ) &
                       *0.25_wp*grd%fact_x(j,i)

           lzx_c(kc) = (  (st%vz_c(kc,j,i+1)+st%vz_c(kc-1,j,i+1)) &
                        - (st%vz_c(kc,j,i-1)+st%vz_c(kc-1,j,i-1)) ) &
                       *0.25_wp*grd%fact_x(j,i)

           lzy_c(kc) = (  (st%vz_c(kc,j+1,i)+st%vz_c(kc-1,j+1,i)) &
                        - (st%vz_c(kc,j-1,i)+st%vz_c(kc-1,j-1,i)) ) &
                       *0.25_wp*grd%fact_y(j,i)

           lxz_c(kc) = (  (st%vx_c(kc+1,j,i)+st%vx_c(kc+1,j,i-1)) &
                        - (st%vx_c(kc-1,j,i)+st%vx_c(kc-1,j,i-1)) ) &
                       *0.25_wp*grd%fact_z_c(kc)*H_c_inv

           lyz_c(kc) = (  (st%vy_c(kc+1,j,i)+st%vy_c(kc+1,j-1,i)) &
                        - (st%vy_c(kc-1,j,i)+st%vy_c(kc-1,j-1,i)) ) &
                       *0.25_wp*grd%fact_z_c(kc)*H_c_inv

        end do

        kc=grd%KCMAX

           dxx_c(kc) = (st%vx_c(kc,j,i)-st%vx_c(kc,j,i-1))*grd%fact_x(j,i)
           dyy_c(kc) = (st%vy_c(kc,j,i)-st%vy_c(kc,j-1,i))*grd%fact_y(j,i)

           lxy_c(kc) = (  (st%vx_c(kc,j+1,i)+st%vx_c(kc,j+1,i-1)) &
                        - (st%vx_c(kc,j-1,i)+st%vx_c(kc,j-1,i-1)) ) &
                       *0.25_wp*grd%fact_y(j,i)

           lyx_c(kc) = (  (st%vy_c(kc,j,i+1)+st%vy_c(kc,j-1,i+1)) &
                        - (st%vy_c(kc,j,i-1)+st%vy_c(kc,j-1,i-1)) ) &
                       *0.25_wp*grd%fact_x(j,i)

           lzx_c(kc) = (st%vz_s(j,i+1)-st%vz_s(j,i-1))*0.5_wp*grd%fact_x(j,i)
           lzy_c(kc) = (st%vz_s(j+1,i)-st%vz_s(j-1,i))*0.5_wp*grd%fact_y(j,i)

           lxz_c(kc) = (  (st%vx_c(kc  ,j,i)+st%vx_c(kc  ,j,i-1)) &
                        - (st%vx_c(kc-1,j,i)+st%vx_c(kc-1,j,i-1)) ) &
                       *0.5_wp*grd%fact_z_c(kc)*H_c_inv

           lyz_c(kc) = (  (st%vy_c(kc  ,j,i)+st%vy_c(kc  ,j-1,i)) &
                        - (st%vy_c(kc-1,j,i)+st%vy_c(kc-1,j-1,i)) ) &
                       *0.5_wp*grd%fact_z_c(kc)*H_c_inv

        ! end kc=grd%KCMAX

        dxy_c(:) = 0.5_wp*(lxy_c+lyx_c)
        dxz_c(:) = 0.5_wp*(lxz_c+lzx_c)
        dyz_c(:) = 0.5_wp*(lyz_c+lzy_c)

        do kc=0, grd%KCMAX

           shear_c_squared(kc) =   dxz_c(kc)*dxz_c(kc) &
                                 + dyz_c(kc)*dyz_c(kc)

           st%de_c(kc,j,i) =  sqrt(dxx_c(kc)*dxx_c(kc) &
                                 + dyy_c(kc)*dyy_c(kc) &
                                 + dxx_c(kc)*dyy_c(kc) &
                                 + dxy_c(kc)*dxy_c(kc) &
                                 + shear_c_squared(kc) + eps_dp**2)

           ! fixme, can this be 0??
           !if (st%de_c(kc,j,i).ne.0._wp) print *,'WARNING: de_c==0'
           if (st%de_c(kc,j,i).ne.0._wp) then
             st%lambda_shear_c(kc,j,i) = sqrt(shear_c_squared(kc))/(st%de_c(kc,j,i)+eps_dp)
           endif

        end do

        if (par%calcmod==1) then

          if ((st%n_cts(j,i) == -1).or.(st%n_cts(j,i) == 0)) then
                            ! cold ice base, temperate ice base

             dxx_t(:)          = dxx_c(0)
             dyy_t(:)          = dyy_c(0)
             dxy_t(:)          = dxy_c(0)
             dxz_t(:)          = dxz_c(0)
             dyz_t(:)          = dyz_c(0)
             st%de_t(:,j,i)           = st%de_c(0,j,i)
             st%lambda_shear_t(:,j,i) = st%lambda_shear_c(0,j,i)

          else   ! st%n_cts(j,i) == 1, temperate ice layer

             H_t_inv = 1.0_wp/(abs(st%H_t(j,i))+eps_dp)

             kt=0

                dxx_t(kt) = (st%vx_t(kt,j,i)-st%vx_t(kt,j,i-1))*grd%fact_x(j,i)
                dyy_t(kt) = (st%vy_t(kt,j,i)-st%vy_t(kt,j-1,i))*grd%fact_y(j,i)

                lxy_t(kt) = (  (st%vx_t(kt,j+1,i)+st%vx_t(kt,j+1,i-1)) &
                             - (st%vx_t(kt,j-1,i)+st%vx_t(kt,j-1,i-1)) ) &
                            *0.25_wp*grd%fact_y(j,i)

                lyx_t(kt) = (  (st%vy_t(kt,j,i+1)+st%vy_t(kt,j-1,i+1)) &
                             - (st%vy_t(kt,j,i-1)+st%vy_t(kt,j-1,i-1)) ) &
                            *0.25_wp*grd%fact_x(j,i)

                lzx_t(kt) = (st%vz_b(j,i+1)-st%vz_b(j,i-1))*0.5_wp*grd%fact_x(j,i)
                lzy_t(kt) = (st%vz_b(j+1,i)-st%vz_b(j-1,i))*0.5_wp*grd%fact_y(j,i)

                lxz_t(kt) = (  (st%vx_t(kt+1,j,i)+st%vx_t(kt+1,j,i-1)) &
                             - (st%vx_t(kt  ,j,i)+st%vx_t(kt  ,j,i-1)) ) &
                            *0.5_wp*grd%fact_z_t*H_t_inv

                lyz_t(kt) = (  (st%vy_t(kt+1,j,i)+st%vy_t(kt+1,j-1,i)) &
                             - (st%vy_t(kt  ,j,i)+st%vy_t(kt  ,j-1,i)) ) &
                            *0.5_wp*grd%fact_z_t*H_t_inv

             ! end kt=0

             do kt=1, grd%KTMAX-1

                dxx_t(kt) = (st%vx_t(kt,j,i)-st%vx_t(kt,j,i-1))*grd%fact_x(j,i)
                dyy_t(kt) = (st%vy_t(kt,j,i)-st%vy_t(kt,j-1,i))*grd%fact_y(j,i)

                lxy_t(kt) = (  (st%vx_t(kt,j+1,i)+st%vx_t(kt,j+1,i-1)) &
                             - (st%vx_t(kt,j-1,i)+st%vx_t(kt,j-1,i-1)) ) &
                            *0.25_wp*grd%fact_y(j,i)

                lyx_t(kt) = (  (st%vy_t(kt,j,i+1)+st%vy_t(kt,j-1,i+1)) &
                             - (st%vy_t(kt,j,i-1)+st%vy_t(kt,j-1,i-1)) ) &
                            *0.25_wp*grd%fact_x(j,i)

                lzx_t(kt) = (  (st%vz_t(kt,j,i+1)+st%vz_t(kt-1,j,i+1)) &
                             - (st%vz_t(kt,j,i-1)+st%vz_t(kt-1,j,i-1)) ) &
                            *0.25_wp*grd%fact_x(j,i)

                lzy_t(kt) = (  (st%vz_t(kt,j+1,i)+st%vz_t(kt-1,j+1,i)) &
                             - (st%vz_t(kt,j-1,i)+st%vz_t(kt-1,j-1,i)) ) &
                            *0.25_wp*grd%fact_y(j,i)

                lxz_t(kt) = (  (st%vx_t(kt+1,j,i)+st%vx_t(kt+1,j,i-1)) &
                             - (st%vx_t(kt-1,j,i)+st%vx_t(kt-1,j,i-1)) ) &
                            *0.25_wp*grd%fact_z_t*H_t_inv

                lyz_t(kt) = (  (st%vy_t(kt+1,j,i)+st%vy_t(kt+1,j-1,i)) &
                             - (st%vy_t(kt-1,j,i)+st%vy_t(kt-1,j-1,i)) ) &
                            *0.25_wp*grd%fact_z_t*H_t_inv

             end do

             kt=grd%KTMAX

                dxx_t(kt) = (st%vx_t(kt,j,i)-st%vx_t(kt,j,i-1))*grd%fact_x(j,i)
                dyy_t(kt) = (st%vy_t(kt,j,i)-st%vy_t(kt,j-1,i))*grd%fact_y(j,i)

                lxy_t(kt) = (  (st%vx_t(kt,j+1,i)+st%vx_t(kt,j+1,i-1)) &
                             - (st%vx_t(kt,j-1,i)+st%vx_t(kt,j-1,i-1)) ) &
                            *0.25_wp*grd%fact_y(j,i)

                lyx_t(kt) = (  (st%vy_t(kt,j,i+1)+st%vy_t(kt,j-1,i+1)) &
                             - (st%vy_t(kt,j,i-1)+st%vy_t(kt,j-1,i-1)) ) &
                            *0.25_wp*grd%fact_x(j,i)

                lzx_t(kt) = (st%vz_m(j,i+1)-st%vz_m(j,i-1))*0.5_wp*grd%fact_x(j,i)
                lzy_t(kt) = (st%vz_m(j+1,i)-st%vz_m(j-1,i))*0.5_wp*grd%fact_y(j,i)

                lxz_t(kt) = (  (st%vx_t(kt  ,j,i)+st%vx_t(kt  ,j,i-1)) &
                             - (st%vx_t(kt-1,j,i)+st%vx_t(kt-1,j,i-1)) ) &
                            *0.5_wp*grd%fact_z_t*H_t_inv

                lyz_t(kt) = (  (st%vy_t(kt  ,j,i)+st%vy_t(kt  ,j-1,i)) &
                             - (st%vy_t(kt-1,j,i)+st%vy_t(kt-1,j-1,i)) ) &
                            *0.5_wp*grd%fact_z_t*H_t_inv

             ! end kt=grd%KTMAX

             dxy_t(:) = 0.5_wp*(lxy_t+lyx_t)
             dxz_t(:) = 0.5_wp*(lxz_t+lzx_t)
             dyz_t(:) = 0.5_wp*(lyz_t+lzy_t)

             do kt=0, grd%KTMAX

                shear_t_squared(kt) =   dxz_t(kt)*dxz_t(kt) &
                                      + dyz_t(kt)*dyz_t(kt)

                st%de_t(kt,j,i)  = sqrt(dxx_t(kt)*dxx_t(kt) &
                                      + dyy_t(kt)*dyy_t(kt) &
                                      + dxx_t(kt)*dyy_t(kt) &
                                      + dxy_t(kt)*dxy_t(kt) &
                                      + shear_t_squared(kt) + eps_dp**2)

                st%lambda_shear_t(kt,j,i) = sqrt(shear_t_squared(kt))/(st%de_t(kt,j,i)+eps_dp)

             end do

          end if

        else if (par%calcmod==0 .or. par%calcmod==2 .or. par%calcmod==3 .or. par%calcmod==-1) then

          dxx_t(:)          = dxx_c(0)
          dyy_t(:)          = dyy_c(0)
          dxy_t(:)          = dxy_c(0)
          dxz_t(:)          = dxz_c(0)
          dyz_t(:)          = dyz_c(0)
          st%de_t(:,j,i)           = st%de_c(0,j,i)
          st%lambda_shear_t(:,j,i) = st%lambda_shear_c(0,j,i)

        else
          stop ' >>> calc_dxyz: Parameter CALCMOD must be -1, 0, 1, 2 or 3!'
        endif

!  ------ Modification of the shear fraction for floating ice (ice shelves)

        if (st%maske(j,i) == 3) then   ! floating ice

           abs_v_ssa_inv = 1.0_wp / &
                           sqrt( 0.25_wp*(st%vx_m(j,i)+st%vx_m(j,i-1))**2 &
                                +0.25_wp*(st%vy_m(j,i)+st%vy_m(j-1,i))**2 &
                                +eps_dp**2)

           nx = -0.5_wp*(st%vy_m(j,i)+st%vy_m(j-1,i)) * abs_v_ssa_inv
           ny =  0.5_wp*(st%vx_m(j,i)+st%vx_m(j,i-1)) * abs_v_ssa_inv

           shear_x_help =          dxx_c(grd%KCMAX)*nx       &
                          +        dxy_c(grd%KCMAX)*ny       &
                          -        dxx_c(grd%KCMAX)*nx**3    &
                          - 2.0_wp*dxy_c(grd%KCMAX)*nx**2*ny &
                          -        dyy_c(grd%KCMAX)*nx*ny**2
                       ! strain rate for ice shelves independent of depth,
                       ! thus surface values used here

           shear_y_help =          dxy_c(grd%KCMAX)*nx       &
                          +        dyy_c(grd%KCMAX)*ny       &
                          -        dxx_c(grd%KCMAX)*nx**2*ny &
                          - 2.0_wp*dxy_c(grd%KCMAX)*nx*ny**2 &
                          -        dyy_c(grd%KCMAX)*ny**3
                       ! strain rate for ice shelves independent of depth,
                       ! thus surface values used here

           lambda_shear_help = sqrt(shear_x_help**2 + shear_y_help**2 + eps_dp**2) &
                               / (st%de_ssa(j,i)+eps_dp)

           st%lambda_shear_c(:,j,i) = lambda_shear_help
           st%lambda_shear_t(:,j,i) = lambda_shear_help

        end if

!  ------ Constrain the shear fraction to reasonable [0,1] interval

        st%lambda_shear_c(:,j,i) = min(max(st%lambda_shear_c(:,j,i), 0.0_wp), 1.0_wp)
        st%lambda_shear_t(:,j,i) = min(max(st%lambda_shear_t(:,j,i), 0.0_wp), 1.0_wp)

     else   ! st%maske(j,i) == 1 or 2; ice-free land or ocean

        dxx_c(:)          = 0.0_wp
        dyy_c(:)          = 0.0_wp
        dxy_c(:)          = 0.0_wp
        dxz_c(:)          = 0.0_wp
        dyz_c(:)          = 0.0_wp
        st%de_c(:,j,i)           = 0.0_wp
        st%lambda_shear_c(:,j,i) = 0.0_wp

        dxx_t(:)          = 0.0_wp
        dyy_t(:)          = 0.0_wp
        dxy_t(:)          = 0.0_wp
        dxz_t(:)          = 0.0_wp
        dyz_t(:)          = 0.0_wp
        st%de_t(:,j,i)           = 0.0_wp
        st%lambda_shear_t(:,j,i) = 0.0_wp

     end if

  end do
  end do

  deallocate(dxx_c) 
  deallocate(dyy_c)
  deallocate(dxy_c)
  deallocate(dxz_c)
  deallocate(dyz_c)

  deallocate(dxx_t)
  deallocate(dyy_t)
  deallocate(dxy_t)
  deallocate(dxz_t)
  deallocate(dyz_t)

  deallocate(lxy_c)
  deallocate(lyx_c)
  deallocate(lxz_c)
  deallocate(lzx_c)
  deallocate(lyz_c)
  deallocate(lzy_c)

  deallocate(lxy_t)
  deallocate(lyx_t)
  deallocate(lxz_t)
  deallocate(lzx_t)
  deallocate(lyz_t)
  deallocate(lzy_t)

  deallocate(shear_c_squared)
  deallocate(shear_t_squared)

  end subroutine calc_dxyz

!-------------------------------------------------------------------------------

end module calc_dxyz_m
!
