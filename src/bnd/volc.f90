!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : v o l c _ m o d
!
!  Purpose : volcanic radiative forcing
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
module volc_mod

   use precision, only : wp
   use ncio

   implicit none

   real(wp), dimension(:), allocatable :: time, volc_data

contains

   subroutine volc_init(volc_file)

   implicit none

   character (len=*), intent(in), optional :: volc_file

   integer :: ntime


   ntime = nc_size(trim(volc_file),"time")
   allocate( time(ntime) )
   allocate( volc_data(ntime) )
   call nc_read(trim(volc_file),"time",time)    ! time in years BP
   call nc_read(trim(volc_file),"volc",volc_data) 


   end subroutine volc_init


    subroutine volc_update(ivolc, time_now, volc)

    implicit none

    integer, intent(in) :: ivolc
    real(wp), intent(in) :: time_now       ! current year BP
    real(wp), intent(inout) :: volc  ! current volcanic radiative forcing

    integer :: imin, i0, i1
    real(wp) :: w0, w1


    if (ivolc.eq.1) then

      ! interpolate from volc_data to current year
      if (time_now.le.minval(time)) then
        volc = volc_data(lbound(volc_data,1))
      elseif (time_now.ge.maxval(time)) then
        volc = volc_data(ubound(volc_data,1))
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
        volc = w0*volc_data(i0) + w1*volc_data(i1)
      endif

    endif


   return

  end subroutine volc_update

end module volc_mod
