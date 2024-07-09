!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : d u s t _ m o d
!
!  Purpose : dust cycle
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
module dust_mod

  use atm_params, only : wp
  use atm_grid, only : im, jm
  use atm_params, only : hatm, c_dust_dry, c_dust_wet, c_dust_mec

  implicit none

  private
  public :: dust

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  d u s t 
  !   Purpose    :  compute dust load, deposition and optical thickness
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine dust(dam, prc, ra2a, hdust, &
      dust_load, dust_dep_dry, dust_dep_wet, dust_dep, dust_ot)

    implicit none

    real(wp), intent(in ) :: dam(:,:)
    real(wp), intent(in ) :: prc(:,:)
    real(wp), intent(in ) :: ra2a(:,:)
    real(wp), intent(in ) :: hdust(:,:)

    real(wp), intent(out) :: dust_load(:,:)
    real(wp), intent(out) :: dust_dep_dry(:,:)
    real(wp), intent(out) :: dust_dep_wet(:,:)
    real(wp), intent(out) :: dust_dep(:,:)
    real(wp), intent(out) :: dust_ot(:,:)

    integer :: i, j
    real(wp) :: heff


    do j=1,jm
      do i=1,im

        !------------------------------------------------
        ! dust load
        heff = hdust(i,j)*hatm/(hdust(i,j)+hatm) ! m
        dust_load(i,j) = dam(i,j)*heff*ra2a(i,j) ! kg/m2

        !------------------------------------------------
        ! dry and wet dust deposition
        dust_dep_dry(i,j) = c_dust_dry / hdust(i,j) * dust_load(i,j)  ! kg/m2/s
        dust_dep_wet(i,j) = c_dust_wet * prc(i,j) * dust_load(i,j)    ! kg/m2/s
        dust_dep(i,j) = dust_dep_wet(i,j) + dust_dep_dry(i,j)         ! kg/m2/s

        !------------------------------------------------
        ! dust optical thickness
        dust_ot(i,j) = dust_load(i,j) * c_dust_mec

      enddo
    enddo


    return

  end subroutine dust

end module dust_mod

