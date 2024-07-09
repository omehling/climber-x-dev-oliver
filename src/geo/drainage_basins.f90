!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : d r a i n a g e _ b a s i n s _  m o d
!
!  Purpose : find drainage basins
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
module drainage_basins_mod

  use climber_grid, only : basin_mask, i_atlantic, i_pacific, i_indian, i_medi

  implicit none

  private
  public :: drainage_basins

contains

  subroutine drainage_basins(i_runoff_hires, j_runoff_hires, i_runoff, j_runoff, &
      drain_basins_ocn_hires, drain_basins_ocn, idivide_pac_atl, idivide_atl_indpac)

    implicit none

    integer, intent(in) :: i_runoff_hires(:,:)
    integer, intent(in) :: j_runoff_hires(:,:)
    integer, intent(in) :: i_runoff(:,:)
    integer, intent(in) :: j_runoff(:,:)
    integer, intent(out) :: drain_basins_ocn_hires(:,:)
    integer, intent(out) :: drain_basins_ocn(:,:)
    integer, intent(out) :: idivide_pac_atl(:)
    integer, intent(out) :: idivide_atl_indpac(:)

    integer :: i, j, ni, nj, ip1


    ! on low resolution grid

    ni = size(i_runoff,1)
    nj = size(j_runoff,2)

    do i=1,ni
      do j=1,nj
        if (i_runoff(i,j).eq.0) then
          drain_basins_ocn(i,j) = 0     ! ocean point
        else
          drain_basins_ocn(i,j) = basin_mask(i_runoff(i,j),j_runoff(i,j))
        endif
      enddo
    enddo

    ! derive index of longitude of continental divide between Pacific and Atlantic
    do j=1,nj
      ! Pacific -> Atlantic
      idivide_pac_atl(j) = 0
      do i=1,ni
        ip1 = i+1
        if (ip1.eq.ni+1) ip1 = 1
        if (drain_basins_ocn(i,j).eq.i_pacific .and. drain_basins_ocn(ip1,j).eq.i_atlantic) then
          idivide_pac_atl(j) = i
          exit
        endif
      enddo
      ! Atlantic -> Indo-Pacific
      idivide_atl_indpac(j) = 0
      do i=1,ni
        ip1 = i+1
        if (ip1.eq.ni+1) ip1 = 1
        if ((drain_basins_ocn(i,j).eq.i_atlantic .or. drain_basins_ocn(i,j).eq.i_medi) .and. (drain_basins_ocn(ip1,j).eq.i_pacific .or. drain_basins_ocn(ip1,j).eq.i_indian)) then
          idivide_atl_indpac(j) = i
          exit
        endif
      enddo
    enddo

    ! on high resolution grid

    ni = size(i_runoff_hires,1)
    nj = size(j_runoff_hires,2)

    do i=1,ni
      do j=1,nj
        if (i_runoff_hires(i,j).eq.0) then
          drain_basins_ocn_hires(i,j) = 0   ! ocean point
        else
          drain_basins_ocn_hires(i,j) = basin_mask(i_runoff_hires(i,j),j_runoff_hires(i,j))
        endif
      enddo
    enddo

    return
  end subroutine drainage_basins

end module drainage_basins_mod
