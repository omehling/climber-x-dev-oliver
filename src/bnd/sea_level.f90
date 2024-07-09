!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s e a _ l e v e l _ m o d
!
!  Purpose : sea level
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
module sea_level_mod

   use precision, only : wp
   use control, only : isea_level, sea_level_const
   use ncio

   implicit none

   real(wp), dimension(:), allocatable :: time, sea_level_data

contains

   subroutine sea_level_init(sea_level_file)

   implicit none

   character (len=*), intent(in) :: sea_level_file

   integer :: ntime


   ntime = nc_size(trim(sea_level_file),"time")
   allocate( time(ntime) )
   allocate( sea_level_data(ntime) )
   call nc_read(trim(sea_level_file),"time",time)    ! time in years BP
   call nc_read(trim(sea_level_file),"seal",sea_level_data) 
    

    end subroutine sea_level_init


    subroutine sea_level_update(time_now,sea_level)

    implicit none

    real(wp), intent(in) :: time_now          ! current year BP
    real(wp), intent(out) :: sea_level  ! current sea level

    integer :: imin, i0, i1
    real(wp) :: w0, w1


    if (isea_level.eq.1) then

      ! interpolate from sea_level_data to current year
      if (time_now.le.time(lbound(time,1))) then
        sea_level = sea_level_data(lbound(sea_level_data,1))
      elseif (time_now.ge.time(ubound(time,1))) then
        sea_level = sea_level_data(ubound(sea_level_data,1))
      else
        imin = minloc(abs(time-time_now),1) 
        if (time(imin).lt.time_now) then
          i0 = imin
          i1 = imin+1
        else
          i0 = imin-1
          i1 = imin
        endif
        w0 = 1._wp - abs(time(i0)-time_now)/(time(i1)-time(i0))
        w1 = 1._wp - w0
        sea_level = w0*sea_level_data(i0) + w1*sea_level_data(i1)
      endif

    else if (isea_level.eq.2) then

      sea_level = sea_level - 0.5_wp 
      sea_level = max(sea_level,-120._wp)

    endif

    return

    end subroutine sea_level_update

end module sea_level_mod

