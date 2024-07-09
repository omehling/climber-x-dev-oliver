!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c o r a l s _ t o p o _ m o d
!
!  Purpose : topographic information for corals
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
module corals_topo_mod

  use precision, only : wp
  use constants, only : pi, r_earth
  use geo_grid, only : ni, nj, ni_topo, nj_topo, n_topo_sur
  use geo_grid, only : i0_topo, i1_topo, j0_topo, j1_topo, lon_topo, lat_topo

  implicit none

  private
  public :: corals_topo

contains


  subroutine corals_topo(z_topo, &
      coral_f_area, coral_f_topo)

    implicit none

    real(wp), intent(in) :: z_topo(:,:)
    real(wp), intent(out) :: coral_f_area(:,:,-250:)
    real(wp), intent(out) :: coral_f_topo(:,:,-250:)

    integer :: i, j, k, ii, jj, iii, jjj, nocn
    real(wp) :: phi1, phi2, theta1, theta2, tmp, dist
    real(wp) :: alpha
    real(wp), dimension(:,:), allocatable :: topography_factor
    real(wp), dimension(:), allocatable :: topography_factor_cell
    logical, dimension(:), allocatable :: mask_lev
    real(wp), dimension(:), allocatable :: z_topo_cell

    real(wp), parameter :: degrees_to_radians = pi/180._wp


    allocate( topography_factor(ni_topo,nj_topo) )
    allocate( topography_factor_cell(n_topo_sur) )
    allocate( mask_lev(n_topo_sur) )
    allocate( z_topo_cell(n_topo_sur) )

    do i=1,ni_topo
      do j=1,nj_topo
        alpha = 0._wp
        do ii=i-1,i+1
          do jj=j-1,j+1
            iii = ii
            if (iii.eq.0) iii = ni_topo
            if (iii.eq.ni_topo+1) iii = 1
            jjj = jj
            jjj = max(1,jjj)
            jjj = min(nj_topo,jjj)
            if (i.eq.iii .and. j.eq.jjj) then
              ! skip 
            else
              ! Compute spherical distance from spherical coordinates.
              ! phi = 90 - latitude
              phi1 = (90.0_wp - lat_topo(jjj))*degrees_to_radians
              phi2 = (90.0_wp - lat_topo(j))*degrees_to_radians
              ! theta = longitude
              theta1 = lon_topo(iii)*degrees_to_radians
              theta2 = lon_topo(i)*degrees_to_radians
              tmp =  sin(phi1)*sin(phi2)*cos(theta1-theta2) + cos(phi1)*cos(phi2)
              dist = acos( tmp )*R_earth
              if (dist.gt.0._wp) then
                alpha = alpha + atan((z_topo(iii,jjj)-z_topo(i,j))/dist)
                !print *,geo%hires%z_topo(iii,jjj)-geo%hires%z_topo(i,j),dist
                !print *,atan((geo%hires%z_topo(iii,jjj)-geo%hires%z_topo(i,j))/dist)
              endif
            endif
          enddo
        enddo
        alpha = min(alpha,1.7_wp)
        alpha = max(alpha,0.013_wp)
        topography_factor(i,j) = log(alpha*100._wp)/5._wp
      enddo
    enddo

    do j=1,nj
      do i=1,ni
        z_topo_cell = pack(z_topo(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)),.true.)
        topography_factor_cell = pack(topography_factor(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)),.true.)
        ! recompute only if coral file does not exist (has not been generated during previous runs)
        do k=-250,50
          mask_lev = z_topo_cell.gt.k .and. z_topo_cell.le.k+1._wp
          nocn = count(mask_lev) 
          ! hypsometry for corals
          coral_f_area(i,j,k) = real(nocn,wp)/real(n_topo_sur,wp)   ! fraction of area between k and k+1
          ! topography factor for corals
          if (nocn.gt.0) then
            coral_f_topo(i,j,k) = sum(topography_factor_cell,mask_lev)/real(nocn,wp)
          else
            coral_f_topo(i,j,k) = 0._wp
          endif
        enddo
      enddo
    enddo

    deallocate(topography_factor)
    deallocate(topography_factor_cell)
    deallocate(mask_lev)

  end subroutine corals_topo

end module corals_topo_mod
