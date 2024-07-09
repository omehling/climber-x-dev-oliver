!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : d 1 3 c _ a t m _ m o d
!
!  Purpose : atmospheric d13C
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
module d13c_atm_mod

   use precision, only : wp
   use ncio

   implicit none

   real(wp), dimension(:), allocatable :: time, d13c_atm_data

contains

   subroutine d13c_atm_init(d13c_atm_file)

   implicit none

   character (len=*), intent(in) :: d13c_atm_file

   integer :: ntime


   ntime = nc_size(d13c_atm_file,"time")
   allocate( d13c_atm_data(ntime) )
   allocate( time(ntime) )
   call nc_read(trim(d13c_atm_file),"time",time)    ! time in years BP
   call nc_read(trim(d13c_atm_file),"d13C",d13c_atm_data) 
    

    end subroutine d13c_atm_init


    subroutine d13c_atm_update(time_now,d13c_atm)

    implicit none

    real(wp), intent(in) :: time_now  ! current year BP
    real(wp), intent(out) :: d13c_atm

    integer :: imin, i0, i1
    real(wp) :: w0, w1


    ! interpolate from d13c_atm_data to current year
    if (time_now.le.time(lbound(time,1))) then
      d13c_atm = d13c_atm_data(lbound(d13c_atm_data,1))
    elseif (time_now.ge.time(ubound(time,1))) then
      d13c_atm = d13c_atm_data(ubound(d13c_atm_data,1))
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
      d13c_atm = w0*d13c_atm_data(i0) + w1*d13c_atm_data(i1)
    endif


    return

    end subroutine d13c_atm_update

end module d13c_atm_mod
