!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c o n n e c t _ o c n _ m o d
!
!  Purpose : ensure connectivity of the ocean
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
module connect_ocn_mod

  use precision, only : wp
  use constants, only : pi
  use geo_params, only : f_crit, f_crit_eq
  use geo_params, only : l_close_panama, l_close_bering, l_fix_cell_stab
  use geo_params, only : l_ocn_below_shelf
  use geo_grid, only : ni, nj, n_topo_sur, i_topo_sur, j_topo_sur, i0_topo, i1_topo, j0_topo, j1_topo
  use fill_ocean_mod, only : fill_ocean_lowres

  implicit none

  private
  public :: connect_ocn 

contains

  subroutine connect_ocn(hires_mask, lon, lat, lon_ocn_origin, lat_ocn_origin, &
      f_ocn)

  implicit none

  integer, intent(inout) :: hires_mask(:,:)
  real(wp), intent(in) :: lat(:)
  real(wp), intent(in) :: lon(:)
  real(wp), intent(in) :: lat_ocn_origin  
  real(wp), intent(in) :: lon_ocn_origin

  real(wp), intent(out) :: f_ocn(:,:)

  integer :: i, j, nocn
  real(wp) :: fcrit

  integer, dimension(:), allocatable :: mask_cell
  integer, dimension(:,:), allocatable :: mask_ocn_connect


  allocate( mask_cell(n_topo_sur) )

  !$omp parallel do private(i,j,fcrit,nocn,mask_cell)
  do j=1,nj
    do i=1,ni

      mask_cell = pack(hires_mask(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)),.true.)

      ! ocean fraction
      if (l_ocn_below_shelf) then
        ! count number of ocean points and floating ice in grid cell
        nocn = count(mask_cell.eq.2 .or. mask_cell.eq.3) 
      else
        ! count number of ocean points in grid cell
        nocn = count(mask_cell.eq.2) 
      endif
      f_ocn(i,j) = real(nocn,wp)/real(n_topo_sur,wp)

      ! apply critical threshold for ocean/land
      if (j==1 .or. j.eq.nj) then
        fcrit = f_crit*2._wp  ! fixme?!
      else
        fcrit = f_crit + max(0._wp,(f_crit_eq-f_crit))*cos(lat(j)*pi/180._wp)**30
      endif
      if (f_ocn(i,j).lt.fcrit) then
        f_ocn(i,j) = 0._wp
      endif

      ! fix for Panama
      if (l_close_panama) then
        ! close Panama, ocean fraction=0
        f_ocn(16,23) = 0._wp
        f_ocn(17:18,22) = 0._wp
        f_ocn(19:20,21) = 0._wp
        f_ocn(21,20) = 0._wp
      endif
      ! option to close Bering 
      if (l_close_bering) then
        ! close Bering, ocean fraction=0
        f_ocn(1:4,32) = 0._wp
      endif
      ! fix to keep special cell causing instabilities in ocean model always wet 
      if (l_fix_cell_stab) then
        f_ocn(46,20) = max(0.10001_wp,f_ocn(46,20))
      endif

    enddo
  enddo
  !$omp end parallel do

  ! make sure that all ocean cells are connected on the low resolution grid
  ! where not, set ocean fraction to zero to avoid having 'ocean lakes'
  allocate( mask_ocn_connect(ni,nj))
  call fill_ocean_lowres(f_ocn,lon,lat,lon_ocn_origin,lat_ocn_origin, &
    mask_ocn_connect)
  where (mask_ocn_connect.eq.0)
    f_ocn = 0._wp
  endwhere

  !$omp parallel do private(i,j,nocn,mask_cell)
  do j=1,nj
    do i=1,ni

      mask_cell = pack(hires_mask(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)),.true.)

      ! count number of ocean points in grid cell
      if (l_ocn_below_shelf) then
        ! count number of ocean points and floating ice in grid cell
        nocn = count(mask_cell.eq.2 .or. mask_cell.eq.3) 
      else
        ! count number of ocean points in grid cell
        nocn = count(mask_cell.eq.2) 
      endif
      ! if number of ocean cells is larger than 0 but the ocean fraction on the
      ! low resolution grid is zero, then set high resolution ocean cells to land and floating ice cells to grounded ice
      if (nocn.gt.0 .and. mask_ocn_connect(i,j).eq.0) then
        where (hires_mask(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j))==2)  ! ocean
          hires_mask(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)) = 1   ! land
        endwhere
        where (hires_mask(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j))==3)  ! floating ice
          hires_mask(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)) = 0   ! grounded ice
        endwhere
      endif

    enddo
  enddo
  !$omp end parallel do

  deallocate( mask_cell )
  deallocate( mask_ocn_connect)


  end subroutine connect_ocn 

end module connect_ocn_mod
