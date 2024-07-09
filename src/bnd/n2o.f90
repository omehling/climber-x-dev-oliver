!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : n 2 o _ m o d
!
!  Purpose : atmospheric N2O concentration
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
module n2o_mod

   use precision, only : wp
   use ncio

   implicit none

   real(wp), dimension(:), allocatable :: time, n2o_data

contains

   subroutine n2o_init(n2o_file)

   implicit none

   character (len=*), intent(in), optional :: n2o_file

   integer :: ntime


   ntime = nc_size(trim(n2o_file),"time")
   allocate( time(ntime) )
   allocate( n2o_data(ntime) )
   call nc_read(trim(n2o_file),"time",time)    ! time in years BP
   call nc_read(trim(n2o_file),"n2o",n2o_data) 


   end subroutine n2o_init


    subroutine n2o_update(in2o, time_now, n2o)

    implicit none

    integer, intent(in) :: in2o
    real(wp), intent(in) :: time_now       ! current year BP
    real(wp), intent(inout) :: n2o  ! current atmospheric n2o concentration

    integer :: imin, i0, i1
    real(wp) :: w0, w1


    if (in2o.eq.1 .or. in2o.eq.-1) then

      ! interpolate from n2o_data to current year
      if (time_now.le.minval(time)) then
        n2o = n2o_data(lbound(n2o_data,1))
      elseif (time_now.ge.maxval(time)) then
        n2o = n2o_data(ubound(n2o_data,1))
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
        n2o = w0*n2o_data(i0) + w1*n2o_data(i1)
      endif

    else if (in2o.eq.2) then

      ! 1 % increase per year
      n2o = n2o*1.01

    endif


   return

  end subroutine n2o_update

end module n2o_mod
