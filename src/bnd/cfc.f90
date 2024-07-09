!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c f c _ m o d
!
!  Purpose : atmospheric CFCs concentrations
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
module cfc_mod

   use precision, only : wp
   use ncio

   implicit none

   real(wp), dimension(:), allocatable :: time, cfc11_data, cfc12_data

   private
   public :: cfc_init, cfc_update

contains

   subroutine cfc_init(cfc_file)

   implicit none

   character (len=*), intent(in) :: cfc_file

   integer :: ntime


   ntime = nc_size(trim(cfc_file),"time")
   allocate( time(ntime) )
   allocate( cfc11_data(ntime) )
   allocate( cfc12_data(ntime) )
   call nc_read(trim(cfc_file),"time",time)    ! time in years BP
   call nc_read(trim(cfc_file),"cfc11",cfc11_data) 
   call nc_read(trim(cfc_file),"cfc12",cfc12_data) 
    

    end subroutine cfc_init


    subroutine cfc_update(time_now,cfc11,cfc12)

    implicit none

    real(wp), intent(in) :: time_now       ! current year BP
    real(wp), intent(out) :: cfc11  ! current atmospheric cfc11 concentration
    real(wp), intent(out) :: cfc12  ! current atmospheric cfc12 concentration

    integer :: imin, i0, i1
    real(wp) :: w0, w1


    ! interpolate from cfc_data to current year
    if (time_now.le.minval(time)) then
      cfc11 = cfc11_data(lbound(cfc11_data,1))
      cfc12 = cfc12_data(lbound(cfc12_data,1))
    elseif (time_now.ge.maxval(time)) then
      cfc11 = cfc11_data(ubound(cfc11_data,1))
      cfc12 = cfc12_data(ubound(cfc12_data,1))
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
      cfc11 = w0*cfc11_data(i0) + w1*cfc11_data(i1)
      cfc12 = w0*cfc12_data(i0) + w1*cfc12_data(i1)
    endif
   return

    end subroutine cfc_update

end module cfc_mod
