!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : l a k e _ r h o _ m o d
!
!  Purpose : lake water density
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
module lake_rho_mod

  use precision, only : wp

  implicit none

  private
  public :: lake_rho

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l a k e _ r h o 
  !   Purpose    :  freshwater density as function of temperature
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function lake_rho(t) result(r)

    real(wp), intent(in) :: t
    real(wp) :: r

    r = 1000._wp*(1._wp-1.9549e-5_wp*abs(t-277._wp)**1.68)

    return

  end function lake_rho

end module lake_rho_mod
