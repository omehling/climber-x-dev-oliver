!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : D 1 4 c _ a t m _ m o d
!
!  Purpose : atmospheric radiocarbon 
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
module D14c_atm_mod

   use precision, only : wp
   use ncio

   implicit none

   real(wp), dimension(:), allocatable :: time, D14c_atm_data

contains

   subroutine D14c_atm_init(D14c_atm_file)

   implicit none

   character (len=*), intent(in) :: D14c_atm_file

   integer :: ntime


   ntime = nc_size(D14c_atm_file,"time")
   allocate( time(ntime) )
   allocate( D14c_atm_data(ntime) )
   call nc_read(trim(D14c_atm_file),"time",time)    ! time in years BP
   call nc_read(trim(D14c_atm_file),"D14C",D14c_atm_data) 
    

    end subroutine D14c_atm_init


    subroutine D14c_atm_update(time_now,D14c_atm)

    implicit none

    real(wp), intent(in) :: time_now  ! current year BP
    real(wp), intent(out) :: D14c_atm

    integer :: imin, i0, i1
    real(wp) :: w0, w1


    ! interpolate from D14c_atm_data to current year
    if (time_now.le.time(lbound(time,1))) then
      D14c_atm = D14c_atm_data(lbound(D14c_atm_data,1))
    elseif (time_now.ge.time(ubound(time,1))) then
      D14c_atm = D14c_atm_data(ubound(D14c_atm_data,1))
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
      D14c_atm = w0*D14c_atm_data(i0) + w1*D14c_atm_data(i1)
    endif


    return

    end subroutine D14c_atm_update

end module D14c_atm_mod
