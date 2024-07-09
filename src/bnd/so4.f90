!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s o 4 _ m o d
!
!  Purpose : atmospheric sulfate aerosol load
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
module so4_mod

   use precision, only : wp
   use ncio

   implicit none

   real(wp), dimension(:), allocatable :: time
   real(wp), dimension(:,:,:), allocatable :: so4_data

contains

   subroutine so4_init(so4_file)

   implicit none

   character (len=*), intent(in), optional :: so4_file

   integer :: ntime, nlon, nlat


   ntime = nc_size(trim(so4_file),"time")
   nlon  = nc_size(trim(so4_file),"lon")
   nlat  = nc_size(trim(so4_file),"lat")
   allocate( time(ntime) )
   allocate( so4_data(nlon,nlat,ntime) )
   call nc_read(trim(so4_file),"time",time)    ! time in years BP
   call nc_read(trim(so4_file),"so4",so4_data) 


   end subroutine so4_init


    subroutine so4_update(iso4, time_now, so4)

    implicit none

    integer, intent(in) :: iso4
    real(wp), intent(in) :: time_now       ! current year BP
    real(wp), intent(inout) :: so4(:,:)  ! current atmospheric so4 load

    integer :: imin, i0, i1
    real(wp) :: w0, w1


    if (iso4.eq.1) then

      ! interpolate from so4_data to current year
      if (time_now.le.minval(time)) then
        so4(:,:) = so4_data(:,:,lbound(so4_data,3))
      elseif (time_now.ge.maxval(time)) then
        so4(:,:) = so4_data(:,:,ubound(so4_data,3))
      else
        imin = minloc(abs(time-time_now),1) 
        if (time(imin).lt.time_now) then
          i0 = imin
          i1 = imin+1
        else
          i0 = imin-1
          i1 = imin
        endif
        if (time(i1)-time(i0).eq.0._wp) then
          w0 = 1._wp
        else
          w0 = 1._wp - abs(time(i0)-time_now)/(time(i1)-time(i0))
        endif
        w1 = 1._wp - w0
        so4(:,:) = w0*so4_data(:,:,i0) + w1*so4_data(:,:,i1)
      endif

    endif


   return

  end subroutine so4_update

end module so4_mod
