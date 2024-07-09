!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : i c e _ p a r _ m o d
!
!  Purpose : ice sheet properies
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit and Andrey Ganopolski
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
module ice_par_mod

  use precision, only : wp
  use constants, only : lambda_i, cap_i
  use constants, only : rho_w, rho_i, T0
  use lnd_grid, only : z_int, z, nl
  use lnd_params, only : snow_par

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i c e _ p a r _ t h e r m a l
  !   Purpose    :  update ice thermal properties
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ice_par_thermal(h_snow, &
                             cap_ice,lambda_ice,lambda_int_ice)

    implicit none

    real(wp), intent(in) :: h_snow
    real(wp), dimension(0:), intent(out) :: cap_ice, lambda_ice, lambda_int_ice

    integer :: k
    real(wp), dimension(0:nl) :: z_loc


     z_loc(0) = -0.5_wp * h_snow
     z_loc(1:nl) = z(1:nl)

     ! snow heat capacity, J/m3/K
     cap_ice(0) = cap_i*snow_par%rho_snow

     ! snow thermal conductivity, W/m/K
     lambda_ice(0) = snow_par%lambda_snow

     do k=1,nl
      lambda_ice(k) = lambda_i
      cap_ice(k) = rho_i * cap_i
     enddo

     ! ice thermal conductivity at the soil levels (interfaces), W/m/K
     do k=0,nl-1
      lambda_int_ice(k) = lambda_ice(k) * lambda_ice(k+1) * (z_loc(k+1) - z_loc(k))  &
                    / ( lambda_ice(k) * (z_loc(k+1) - z_int(k)) + lambda_ice(k+1) * (z_int(k) - z_loc(k)) )
     enddo
     lambda_int_ice(nl) = 0._wp


    return

  end subroutine ice_par_thermal


end module ice_par_mod
