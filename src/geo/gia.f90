!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : g i a _ m o d
!
!  Purpose : Computation of the glacial isostatic adjustment of the lithosphere
!            surface
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
module gia_mod
  
  use precision, only : wp
  use constants, only : rho_i, rho_as, rho_sw
  use timer, only : sec_year, dt_geo
  use geo_params, only : gia_time_lag

  implicit none

  private
  public :: gia_init, gia_update

contains

  subroutine gia_init(z_bed, h_ice, &
      z_bed_rel)

    implicit none

    real(wp), intent(in) :: z_bed(:,:)
    real(wp), intent(in) :: h_ice(:,:)
    real(wp), intent(out) :: z_bed_rel(:,:)


    ! relaxed bedrock topography
    ! Isostatic bedrock adjustment with local lithosphere and relaxing asthenosphere (LLRA model)
    z_bed_rel = z_bed + rho_i/rho_as * h_ice

    return

  end subroutine gia_init


  subroutine gia_update(z_bed_rel, h_ice, mask, sea_level, d_sea_level, &
      z_bed)

    implicit none

    real(wp), intent(in) :: z_bed_rel(:,:)
    real(wp), intent(in) :: h_ice(:,:)
    integer,  intent(in) :: mask(:,:)
    real(wp), intent(in) :: sea_level
    real(wp), intent(in) :: d_sea_level

    real(wp), intent(inout) :: z_bed(:,:)

 
    ! add sea level change
    z_bed = z_bed - d_sea_level

    ! convert bedrock elevation to be relative to present-day sea level 
    z_bed = z_bed + sea_level

    where (mask.eq.0 .or. mask.eq.1)  ! grounded ice or ice-free land 

      z_bed = 1.0_wp/(gia_time_lag*sec_year+dt_geo)*( gia_time_lag*sec_year*z_bed &
        + dt_geo*(z_bed_rel-rho_i/rho_as*h_ice) )

    elsewhere   ! ocean or shelf ice

      z_bed = 1.0_wp/(gia_time_lag*sec_year+dt_geo)*( gia_time_lag*sec_year*z_bed &
        + dt_geo*(z_bed_rel-rho_sw/rho_as*sea_level) )
      ! Water load relative to the present sea-level stand (0 m)
      ! -> can be positive or negative

    endwhere

    ! convert bedrock elevation back to be relative to current sea level 
    z_bed = z_bed - sea_level


    return

end subroutine gia_update

!!-------- Smoothing of the topography (Gaussian filter) --------
!
!filter_width = dx    ! half span of filtered area, in km
!                     ! (chosen as one grid spacing)
!
!sigma_filter = filter_width/dx   ! half span of filtered area,
!                                 ! in grid points
!
!n_filter     = ceiling(2.0_wp*sigma_filter)
!n_filter     = max(n_filter, 5)
!
!do i=0, grd%IMAX
!do j=0, grd%JMAX
!
!   sum_weigh = 0.0_wp
!   st%zl0_smoothed(j,i) = 0._wp
!
!   do m=-n_filter, n_filter
!   do n=-n_filter, n_filter
!
!      i_f = i+m
!      j_f = j+n
!
!      if (i_f <    0) i_f =    0
!      if (i_f > grd%IMAX) i_f = grd%IMAX
!
!      if (j_f <    0) j_f =    0
!      if (j_f > grd%JMAX) j_f = grd%JMAX
!
!      dist      = sqrt(real(m,dp)**2+real(n,dp)**2)
!      weigh     = exp(-(dist/sigma_filter)**2)
!      sum_weigh = sum_weigh + weigh
!
!      st%zl0_smoothed(j,i) = st%zl0_smoothed(j,i) + weigh*st%zl0_raw(j_f,i_f)
!
!   end do
!   end do
!
!   st%zl0_smoothed(j,i) = st%zl0_smoothed(j,i)/sum_weigh
!
!end do
!end do
!

end module gia_mod
!
