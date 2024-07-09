!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c o 2 _ r a d _ m o d
!
!  Purpose : atmospheric CO2 concentration for radiation
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
module co2_rad_mod

   use precision, only : wp
   use ncio

   implicit none

   real(wp), dimension(:), allocatable :: time, co2_rad_data

contains

   subroutine co2_rad_init(co2_rad_file)

   implicit none

   character (len=*), intent(in), optional :: co2_rad_file

   integer :: ntime


   ntime = nc_size(trim(co2_rad_file),"time")
   allocate( time(ntime) )
   allocate( co2_rad_data(ntime) )
   call nc_read(trim(co2_rad_file),"time",time)    ! time in years BP
   call nc_read(trim(co2_rad_file),"co2",co2_rad_data) 


   end subroutine co2_rad_init


    subroutine co2_rad_update(ico2_rad, time_now, co2_rad)

    implicit none

    integer, intent(in) :: ico2_rad
    real(wp), intent(in) :: time_now       ! current year BP
    real(wp), intent(inout) :: co2_rad  ! current atmospheric CO2 concentration

    integer :: imin, i0, i1
    real(wp) :: w0, w1


    if (ico2_rad.eq.2) then

      ! interpolate from co2_rad_data to current year
      if (time_now.le.minval(time)) then
        co2_rad = co2_rad_data(lbound(co2_rad_data,1))
      elseif (time_now.ge.maxval(time)) then
        co2_rad = co2_rad_data(ubound(co2_rad_data,1))
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
        co2_rad = w0*co2_rad_data(i0) + w1*co2_rad_data(i1)
      endif

    else if (ico2_rad.eq.3) then

      ! 1 % increase per year
      co2_rad = co2_rad*1.01

    else if (ico2_rad.eq.4) then

      if (time_now.lt.2000._wp) then
        co2_rad = co2_rad - 0.05_wp
      else
        co2_rad = co2_rad + 0.05_wp
      endif

    endif


   return

  end subroutine co2_rad_update

end module co2_rad_mod
