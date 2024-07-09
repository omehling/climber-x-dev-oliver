!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : i n v e r t _ m o d
!
!  Purpose : invert matrix for barotropic streamfunction
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
module invert_mod

  use precision, only : wp
  use ocn_grid, only : maxi, maxj, maxk, k1, c, cv, ds, dsv, dphi, rh, rcv, rdphi, R_earth
  use ocn_params, only : fcor, fcorv, drag, drhcor_max

  implicit none

  private
  public :: invert

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i n v e r t
  !   Purpose    :  invert matrix for barotropic streamfunction
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine invert(gap,ratm)

    implicit none

    real(wp), dimension(:,:), intent(out) :: gap
    real(wp), dimension(:,:), intent(out) :: ratm


    integer i, j, k, l, n, m, im
    real(wp) :: tv, tv1, rh1, rh2, drh, rat


    n = maxi
    m = maxj + 1

    ! Set equation at Psi points, assuming periodic b.c. in i.
    ! Cannot solve at both i=0 and i=maxi as periodicity => would
    ! have singular matrix. At dry points equation is trivial.

    gap(:,:) = 0._wp

    !$omp parallel do private(i,j,k,l,tv,tv1,rh1,rh2,drh)
    do i=1,maxi
       do j=0,maxj
          k=i + j*n
          if (max(k1(i,j),k1(i+1,j),k1(i,j+1),k1(i+1,j+1)).le.maxk) then
             if (j.eq.0) stop 'j==0'
             if (j.eq.maxj) stop 'j==maxj'
             ! Coriolis terms, limit topography gradient
             rh1 = rh(1,i,j+1)
             rh2 = rh(1,i,j)
             drh = abs((rh1-rh2)/dphi)
             if (drh.gt.drhcor_max) then
               drh = drhcor_max
               rh2 = rh1 - sign(drh,(rh1-rh2)/dphi)*dphi
             endif
             tv  = (fcor(j+1)*rh1 - fcor(j)*rh2) /(2._wp*dphi*dsv(j)*R_earth**2)  ! 1/s/m3
             rh1 = rh(2,i+1,j)
             rh2 = rh(2,i,j)
             drh = abs((rh1-rh2)/dsv(j))
             if (drh.gt.drhcor_max) then
               drh = drhcor_max
               rh2 = rh1 - sign(drh,(rh1-rh2)/dsv(j))*dsv(j)
             endif
             tv1 = (fcorv(j)*rh1  - fcorv(j)*rh2)/(2._wp*dphi*dsv(j)*R_earth**2)

             gap(k,2) = drag(1,i,j)*c(j)**2*rh(1,i,j)/(ds(j)*dsv(j)*R_earth**2) + tv1  ! 1/s/m3

             l = n+1
             ! for periodic boundary in i
             if (i.eq.1) l = 2*n+1

             gap(k,l) = drag(2,i,j)*rh(2,i,j)/(cv(j)*dphi*R_earth)**2 - tv

             gap(k,n+2) = - (drag(2,i,j)*rh(2,i,j) + drag(2,i+1,j)*rh(2,i+1,j)) &
                        /(cv(j)*dphi*R_earth)**2 &
                        - (drag(1,i,j)*c(j)**2*rh(1,i,j)/ds(j) + drag(1,i,j+1)*c(j+1)**2*rh(1,i,j+1)/ds(j+1)) &
                        /(dsv(j)*R_earth**2)

             l = n+3
             ! for periodic boundary in i
             if (i.eq.maxi) l=3

             gap(k,l) = drag(2,i+1,j)*rh(2,i+1,j)/(cv(j)*dphi*R_earth)**2 + tv
             gap(k,2*n+2) = drag(1,i,j+1)*c(j+1)**2*rh(1,i,j+1)/(ds(j+1)*dsv(j)*R_earth**2) - tv1

          else

             gap(k,n+2) = 1._wp 

          endif
       enddo
    enddo
    !$omp end parallel do

    ! now invert the thing

    do i=1,n*m-1
       im = min(i+n+1,n*m)
       do j=i+1,im
          rat = gap(j,n+2-j+i)/gap(i,n+2)
          ratm(j,j-i) = rat
          if (rat.ne.0) then
             do k=n+2-j+i,2*n+3-j+i
                gap(j,k)=gap(j,k) - rat*gap(i,k+j-i)
             enddo
          endif
       enddo
    enddo

   return

  end subroutine invert

end module invert_mod
