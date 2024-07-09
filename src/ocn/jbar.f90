!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : j b a r _ m o d
!
!  Purpose : calculate jbar forcing for streamfunction
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
module jbar_mod

  use precision, only : wp
  use constants, only : g
  use ocn_grid, only : maxi, maxj, maxk, k1, dz, dza, ku, getj, rdphi, rdsv, R_earth, rh
  !$use omp_lib

  implicit none

  private
  public :: jbar

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  j b a r
  !   Purpose    :  calculate jbar forcing for streamfunction
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine jbar(rho, &
                  bp, ubar_jbar)
    
    implicit none

    real(wp), dimension(:,:,:), intent(in) :: rho

    real(wp), dimension(:,:,:), intent(out) :: bp
    real(wp), dimension(:,:), intent(out) :: ubar_jbar

    integer :: i, j, k, ip1
    real(wp) :: tv1, tv2, tv3, tv4


    ! calculate easy part of double p integral
    !$omp parallel do private(i,j,k)
    do j=1,maxj
       do i=1,maxi
          bp(i,j,:) = 0._wp
          if (k1(i,j).le.maxk) then
             do k=k1(i,j)+1,maxk
                bp(i,j,k) = bp(i,j,k-1) - g*(rho(i,j,k) + rho(i,j,k-1))*dza(k-1)*0.5_wp  ! kg/m/s2
             enddo
          endif
       enddo
       ! periodic b.c. 
       if (k1(1,j).lt.maxk) then
          do k=k1(1,j),maxk
             bp(maxi+1,j,k) = bp(1,j,k)
          enddo
       endif
    enddo
    !$omp end parallel do

    ubar_jbar(:,:) = 0._wp
    !$omp parallel do private(i,j,k,ip1,tv1,tv2,tv3,tv4)
    do j=1,maxj-1
      do i=1,maxi
        ip1=mod(i,maxi) + 1
        if (getj(i,j)) then
          tv1 = 0._wp
          do k=ku(2,ip1,j),maxk
            tv1 = tv1 + bp(ip1,j+1,k)*dz(k)
          enddo
          tv2 = 0._wp
          do k=ku(2,ip1,j),maxk
            tv2 = tv2 + bp(ip1,j,k)*dz(k)
          enddo
          tv3 = 0._wp
          do k=ku(2,i,j),maxk
            tv3 = tv3 + bp(i,j+1,k)*dz(k)
          enddo
          tv4 = 0._wp
          do k=ku(2,i,j),maxk
            tv4 = tv4 + bp(i,j,k)*dz(k)
          enddo
          ubar_jbar(i,j) = &
            + ((tv3 - tv4)*rh(2,i,j) &  ! kg/s2/m3
            - (tv1 - tv2)*rh(2,ip1,j)) &
            * rdphi*rdsv(j)/R_earth**2
          tv1 = 0._wp
          do k=ku(1,i,j+1),maxk
            tv1 = tv1 + bp(ip1,j+1,k)*dz(k)
          enddo
          tv2 = 0._wp
          do k=ku(1,i,j),maxk
            tv2 = tv2 + bp(ip1,j,k)*dz(k)
          enddo
          tv3 = 0._wp
          do k=ku(1,i,j+1),maxk
            tv3 = tv3 + bp(i,j+1,k)*dz(k)
          enddo
          tv4 = 0._wp
          do k=ku(1,i,j),maxk
            tv4 = tv4 + bp(i,j,k)*dz(k)
          enddo
          ubar_jbar(i,j) = ubar_jbar(i,j) &
            + ((tv1 - tv3)*rh(1,i,j+1) &
            - (tv2 - tv4)*rh(1,i,j)) &
            * rdphi*rdsv(j)/R_earth**2
        endif
      enddo
    enddo
    !$omp end parallel do

   return
         
  end subroutine jbar

end module jbar_mod
