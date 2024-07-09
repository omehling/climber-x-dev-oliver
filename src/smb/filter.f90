!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : f i l t e r 
!
!  Purpose : Gaussian filter 
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
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
module filter

  use precision, only : wp

  implicit none

  private 
  public :: filter_gaussian

contains

  subroutine filter_gaussian(field_in, field_out)

    implicit none

    real(wp), intent(in) :: field_in(:,:)
    real(wp), intent(out) :: field_out(:,:)

    integer :: i, j, i, jj, i_f, j_f, ni, nj, n_filter
    real(wp) :: sigma_filter_pts
    real(wp) :: dist, weigh, sum_weight, field_tmp
    real(wp), allocatable :: weight(:,:)


    ni = size(field_in,1)
    nj = size(field_in,2)

    sigma_filter_pts = sigma_filter/dx   ! half span of filtered area, in grid points
    n_filter     = ceiling(2.0_wp*sigma_filter_pts)

    ! compute weights for n_filter x n_filter square
    allocate(weight(-n_filter:n_filter,-n_filter:n_filter)
    do ii=-n_filter, n_filter
      do jj=-n_filter, n_filter
        dist = sqrt(real(ii,wp)**2+real(jj,wp)**2)
        weight(ii,jj) = exp(-(dist/sigma_filter_pts)**2)
      enddo
    enddo

    do i=1,ni
      do j=1,nj
        sum_weight = 0._wp
        field_tmp  = 0._wp
        do ii=-n_filter, n_filter
          do jj=-n_filter, n_filter
            i_f = i+ii
            j_f = j+jj
            if (i_f < 1)  i_f = 1
            if (i_f > ni) i_f = ni
            if (j_f < 1)  j_f = 1
            if (j_f > nj) j_f = nj
            sum_weight = sum_weight + weight(ii,jj)
            field_tmp = field_tmp + weight(ii,jj)*field_in(i_f,j_f)
          enddo
        enddo
        field_out(i,j) = field_tmp/sum_weight
      enddo
    enddo

    return

  end subroutine filter_gaussian

end module filter
