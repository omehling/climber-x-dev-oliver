!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : w i n d _ m o d
!
!  Purpose : wind stress term for barotropic streamfunction
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Neil R. Edwards and Matteo Willeit
!
! This file is part of CLIMBER-X.
!
! This file was ported from the original c-GOLDSTEIN model,
! see Edwards and Marsh (2005)
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
module wind_mod

  use precision, only : wp
  use ocn_grid, only : maxi, maxj, maxk, k1, c, rcv, rdsv, rdphi, R_earth, rh

  implicit none

  private
  public :: wind

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  w i n d
  !   Purpose    :  wind stress term for barotropic streamfunction
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine wind(tau, &
                  ubar_wind) 

    implicit none

    real(wp), dimension(:,:,:), intent(in) :: tau
    real(wp), dimension(:,:), intent(out) :: ubar_wind

    integer i, j, ip1


    !$omp parallel do private(i,j,ip1)
    do j=1,maxj
      do i=1,maxi
          ip1=mod(i,maxi) + 1
          if (max(k1(i,j),k1(i+1,j),k1(i,j+1),k1(i+1,j+1)).le.maxk) then
            ubar_wind(i,j) = (tau(2,ip1,j)*rh(2,i+1,j) - tau(2,i,j)*rh(2,i,j)) &  ! N/m4 or kg/s2/m3
                   * rdphi*rcv(j)/R_earth &
                   - (tau(1,i,j+1)*c(j+1)*rh(1,i,j+1) - tau(1,i,j)*c(j)*rh(1,i,j)) &
                   * rdsv(j)/R_earth
          else
            ubar_wind(i,j) = 0._wp
          endif
       enddo
    enddo
    !$omp end parallel do

   return

  end subroutine wind

end module wind_mod
