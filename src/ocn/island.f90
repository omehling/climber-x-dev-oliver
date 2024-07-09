!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : i s l a n d _ m o d
!
!  Purpose : calculate path integral around islands
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
module island_mod

  use precision, only : wp
  use ocn_grid, only : maxi, maxj, maxk, c, dphi, rcv, dsv, dz, rh, R_earth, ku
  use ocn_grid, only : npi, lpisl, ipisl, jpisl
  use ocn_params, only : fcor, fcorv, drag, rho0
  !$use om_lib

  implicit none

  private
  public :: island

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i s l a n d
  !   Purpose    :  calculate path integral around island
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine island(ubloc,tau,bp, &
                   isl,indj, &
                   erisl1)

    implicit none

    real(wp), dimension(:,0:,0:), intent(in) :: ubloc
    real(wp), dimension(:,:,:), intent(in) :: tau
    real(wp), intent(in) :: bp(:,:,:)
    integer, intent(in) :: isl, indj
  
    real(wp), intent(out) :: erisl1

    integer :: i, k, lpi, ipi, jpi
    real(wp) :: cor, tv1, tv2


    erisl1 = 0._wp
    !$omp parallel do private (i,k,lpi,ipi,jpi,cor,tv1,tv2) reduction(+:erisl1)
    do i=1,npi(isl)
       lpi = lpisl(i,isl)
       ipi = ipisl(i,isl)
       jpi = jpisl(i,isl)
       if (abs(lpi).eq.1) then
          !cor = - fcor(jpi)*0.25_wp*(ubloc(2,ipi,jpi) + ubloc(2,ipi+1,jpi) &  ! m/s2
          !    + ubloc(2,ipi,jpi-1) + ubloc(2,ipi+1,jpi-1))
          cor = - 0.5_wp * (fcorv(jpi)*0.5_wp*(ubloc(2,ipi,jpi) + ubloc(2,ipi+1,jpi))  &  ! m/s2
              + fcorv(jpi-1)*0.5_wp*(ubloc(2,ipi,jpi-1) + ubloc(2,ipi+1,jpi-1)))
       else
          !cor = fcorv(jpi)*0.25_wp*(ubloc(1,ipi-1,jpi) + ubloc(1,ipi,jpi) &
          !    + ubloc(1,ipi-1,jpi+1) + ubloc(1,ipi,jpi+1))
          cor = 0.5_wp * (fcor(jpi)*0.5_wp*(ubloc(1,ipi-1,jpi) + ubloc(1,ipi,jpi)) &
              + fcor(jpi+1)*0.5_wp*(ubloc(1,ipi-1,jpi+1) + ubloc(1,ipi,jpi+1)))
       endif

       if (jpi.lt.maxj) then ! to avoid acessing rcv(maxj) and dsv(maxj), out of bound!
          erisl1 = erisl1 + sign(1,lpi)*(drag(abs(lpi),ipi,jpi) &  ! m/s2
              * ubloc(abs(lpi),ipi,jpi) + cor &
              - indj*tau(abs(lpi),ipi,jpi)/rho0*rh(abs(lpi),ipi,jpi)) &
              * (c(jpi)*dphi*(2._wp - abs(lpi)) &
              + rcv(jpi)*dsv(jpi)*(abs(lpi) - 1._wp))  
       else
          erisl1 = erisl1 + sign(1,lpi)*(drag(abs(lpi),ipi,jpi) &  ! m/s2
              * ubloc(abs(lpi),ipi,jpi) + cor &
              - indj*tau(abs(lpi),ipi,jpi)/rho0*rh(abs(lpi),ipi,jpi)) &
              * (c(jpi)*dphi*(2._wp - abs(lpi)))
       endif

       ! calc tricky bits and add to source term for path integral round 
       ! islands, all sums have at least one element
     
       if (indj.eq.1) then
          if (abs(lpi).eq.1) then
             tv1 = 0._wp
             do k=ku(1,ipi,jpi),maxk
                tv1 = tv1 + bp(ipi+1,jpi,k)*dz(k)  ! kg/s2
             enddo
             tv2 = 0._wp
             do k=ku(1,ipi,jpi),maxk
                tv2 = tv2 + bp(ipi,jpi,k)*dz(k)
             enddo
             erisl1 = erisl1 + (tv1-tv2) &  ! m/s2
                    * sign(1,lpi)*rh(1,ipi,jpi)/(rho0*R_earth)
          else
             tv1 = 0._wp
             do k=ku(2,ipi,jpi),maxk
                tv1 = tv1 + bp(ipi,jpi+1,k)*dz(k)
             enddo
             tv2 = 0._wp
             do k=ku(2,ipi,jpi),maxk
                tv2 = tv2 + bp(ipi,jpi,k)*dz(k)
             enddo
             erisl1 = erisl1 + (tv1 - tv2) &
                    * sign(1,lpi)*rh(2,ipi,jpi)/(rho0*R_earth)
          endif
       endif
    enddo
    !$omp end parallel do

   return

  end subroutine island

end module island_mod
