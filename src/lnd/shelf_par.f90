!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s h e l f _ p a r _ m o d
!
!  Purpose : shelf parameters
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
module shelf_par_mod

  use precision, only : wp
  use constants, only : lambda_i, lambda_w, cap_i, cap_w
  use constants, only : rho_w, rho_i
  use lnd_grid, only : z_int, z, nl
  use lnd_params, only : soil_par 

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s h e l f _ p a r _ t h e r m a l
  !   Purpose    :  update shelf soil thermal properties
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine shelf_par_thermal(theta_w_shelf,theta_i_shelf,theta_sat,lambda_s, &
                              cap_shelf,lambda_int_shelf)

    implicit none

    real(wp), dimension(:), intent(in) :: theta_w_shelf, theta_i_shelf, theta_sat, lambda_s
    real(wp), dimension(:), intent(out) :: cap_shelf
    real(wp), dimension(0:), intent(out) :: lambda_int_shelf

    integer :: k
    real(wp), dimension(nl) :: lambda_shelf
    real(wp), dimension(0:nl) :: z_loc


     z_loc(0) = 0._wp
     z_loc(1:nl) = z(1:nl)


     do k=1,nl

      ! thermal properties
      lambda_shelf(k) = lambda_s(k)**(1._wp - theta_sat(k)) * lambda_w**theta_w_shelf(k) * lambda_i**theta_i_shelf(k)
 
      ! heat capacity of shelf, J/m3/K
      cap_shelf(k) = (1._wp - theta_sat(k)) * soil_par%cap_s      &
                  + theta_w_shelf(k) * rho_w * cap_w  &
                  + theta_i_shelf(k) * rho_i * cap_i

     enddo

     ! shelf thermal conductivity at the soil levels (interfaces), W/m/K
     lambda_int_shelf(0) = lambda_shelf(1)
     do k=1,nl-1
      lambda_int_shelf(k) = lambda_shelf(k) * lambda_shelf(k+1) * (z_loc(k+1) - z_loc(k))  &
                    / ( lambda_shelf(k) * (z_loc(k+1) - z_int(k)) + lambda_shelf(k+1) * (z_int(k) - z_loc(k)) )
     enddo
     lambda_int_shelf(nl) = 0._wp


    return

  end subroutine shelf_par_thermal

end module shelf_par_mod
