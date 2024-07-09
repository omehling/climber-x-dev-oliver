!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s o l _ m o d
!
!  Purpose : solar 'constant'
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
module sol_mod

   use precision, only : wp
   use ncio

   implicit none

   real(wp), dimension(:), allocatable :: time, sol_data

contains

   subroutine sol_init(sol_file)

   implicit none

   character (len=*), intent(in), optional :: sol_file

   integer :: ntime


   ntime = nc_size(trim(sol_file),"time")
   allocate( time(ntime) )
   allocate( sol_data(ntime) )
   call nc_read(trim(sol_file),"time",time)    ! time in years BP
   call nc_read(trim(sol_file),"sol",sol_data) 


   end subroutine sol_init


    subroutine sol_update(isol, time_now, sol)

    implicit none

    integer, intent(in) :: isol
    real(wp), intent(in) :: time_now       ! current year BP
    real(wp), intent(inout) :: sol  ! current atmospheric sol concentration

    integer :: imin, i0, i1
    real(wp) :: w0, w1


    if (isol.eq.1) then

      ! interpolate from sol_data to current year
      if (time_now.le.minval(time)) then
        sol = sol_data(lbound(sol_data,1))
      elseif (time_now.ge.maxval(time)) then
        sol = sol_data(ubound(sol_data,1))
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
        sol = w0*sol_data(i0) + w1*sol_data(i1)
      endif

    endif


   return

  end subroutine sol_update

end module sol_mod
