!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : r e s t o r e _ s a l i n i t y _ m o d
!
!  Purpose : Correction of salinity to globally averaged reference salinity
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Andrey Ganopolski and Matteo Willeit
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
module restore_salinity_mod

  use precision, only : wp
  use ocn_grid
  use ocn_params, only : saln0

  implicit none

  private
  public :: restore_salinity

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  r e s t o r e _ s a l i n i t y
  !   Purpose    :  restore global salinity to avoid drifts
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine restore_salinity(sal)

    implicit none

    real(wp), dimension(:,:,:) :: sal

    integer :: i, j, k
    real(wp) :: sglob, scorr


    sglob = 0._wp

    do k=1,maxk
       do j=1,maxj
          do i=1,maxi
             sglob = sglob + sal(i,j,k)*ocn_vol(i,j,k)
          enddo
       enddo
    enddo

    sglob = sglob/ocn_vol_tot

    scorr = saln0 - sglob

!      if (KSEM.ne.2) GLOB_SAL=34.7
!      GLOB_SAL=34.7+35.*(0.9*TOI_VOL/OVOL)  ! account for ice sheet volume (ocean volume)
!      SCORR=GLOB_SAL-SGLOB

    do k=1,maxk
       do j=1,maxj
          do i=1,maxi
             sal(i,j,k) = sal(i,j,k) + scorr
          enddo
       enddo
    enddo

    where (sal.lt.0._wp) sal = 0._wp


   return

  end subroutine restore_salinity

end module restore_salinity_mod
