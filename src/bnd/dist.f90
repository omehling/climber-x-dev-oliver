!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : d i s t _ m o d
!
!  Purpose : vegetation disturbance
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
module dist_mod

   use precision, only : wp
   use timer, only : sec_year
   use ncio

   implicit none

   real(wp), dimension(:), allocatable :: time
   real(wp), dimension(:,:,:,:), allocatable :: disturbance_data

contains

   subroutine dist_init(idist, dist_file)

   implicit none

   integer, intent(in) :: idist
   character (len=*), intent(in) :: dist_file

   integer :: ntime, nlon, nlat


   if (idist.eq.1) then

     nlon  = nc_size(trim(dist_file),"lon")
     nlat  = nc_size(trim(dist_file),"lat")
     allocate( disturbance_data(5,nlon,nlat,1) )
     call nc_read(trim(dist_file),"dist",disturbance_data) 

   else if (idist.eq.2) then

     ntime = nc_size(trim(dist_file),"time")
     nlon  = nc_size(trim(dist_file),"lon")
     nlat  = nc_size(trim(dist_file),"lat")
     allocate( time(ntime) )
     allocate( disturbance_data(5,nlon,nlat,ntime) )
     call nc_read(trim(dist_file),"time",time)    ! time in years BP
     call nc_read(trim(dist_file),"dist",disturbance_data) 

   endif

   end subroutine dist_init


    subroutine dist_update(idist, time_now, disturbance)

    implicit none

    integer, intent(in) :: idist
    real(wp), intent(in) :: time_now       ! current year BP
    real(wp), intent(inout) :: disturbance(:,:,:)  ! current crop fraction

    integer :: imin, i0, i1
    real(wp) :: w0, w1


    if (idist.eq.1) then

      disturbance(:,:,:)    = disturbance_data(:,:,:,1)

    else if (idist.eq.2) then

      ! interpolate to current year
      if (time_now.le.minval(time)) then
        disturbance(:,:,:)    = disturbance_data(:,:,:,lbound(disturbance_data,3))
      elseif (time_now.ge.maxval(time)) then
        disturbance(:,:,:)    = disturbance_data(:,:,:,ubound(disturbance_data,3))
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
        disturbance(:,:,:)    = w0*disturbance_data(:,:,:,i0)    + w1*disturbance_data(:,:,:,i1)
      endif

    else if (idist.eq.3) then
      ! manually set disturbance rate

      disturbance(:,:,:) = 1._wp/(1.e6_wp*sec_year)     ! 1/s
      if (time_now.gt.0._wp) then
        disturbance(1,:,:) = 1._wp/(50._wp*sec_year) * min(1._wp,time_now/200._wp)  ! tropical forest
      endif

    else if (idist.eq.4) then
      ! manually set disturbance rate

      disturbance(:,:,:) = 1._wp/(1.e6_wp*sec_year)     ! 1/s
      if (time_now.gt.0._wp) then
        disturbance(2,:,:) = 1._wp/(50._wp*sec_year) * min(1._wp,time_now/200._wp)  ! boreal forest
      endif

    else if (idist.eq.5) then
      ! manually set disturbance rate

      disturbance(:,:,:) = 1._wp/(1.e6_wp*sec_year)     ! 1/s
      if (time_now.gt.0._wp) then
        disturbance(1,:,:) = 1._wp/(50._wp*sec_year) * min(1._wp,time_now/200._wp)  ! tropical forest
        disturbance(2,:,:) = 1._wp/(50._wp*sec_year) * min(1._wp,time_now/200._wp)  ! boreal forest
      endif

    endif


   return

  end subroutine dist_update

end module dist_mod
