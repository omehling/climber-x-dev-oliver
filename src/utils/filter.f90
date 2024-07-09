!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : f i l t e r
!
!  Purpose : filter routines
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
module filter

  use precision, only : wp

  implicit none

  public :: filter1d, smooth2

contains

  !***********************************************************************
  !      subroutine f i l t e r 1 d
  !***********************************************************************
  !     filtering nfil times, assuming periodic boundary conditions 
  !***********************************************************************
  subroutine filter1d(X,nfil)

    implicit none

    real(wp), intent(inout) :: X(:)
    integer, intent(in) :: nfil

    integer :: n, i, im, ipl, imi
    real(wp) :: xm
    real(wp), allocatable :: Y(:)


    im = size(X,1)

    if (nfil==-1) then

      ! fill with average value
      xm = 0._wp
      do i=1,im
        xm = xm + X(i)/real(im,wp)
      enddo
      X(:) = xm

    else

      allocate(Y(im))

      do n=1,nfil

        Y(:) = X(:)

        do i=1,im
          ipl=i+1
          if (ipl.gt.im) ipl=1
          imi=i-1
          if (imi.lt.1) imi=im 
          X(i)=(Y(imi)+Y(i)+Y(ipl))/3._wp
        enddo

      enddo

      deallocate(Y)

    endif


    return
  end subroutine filter1d


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s m o o t h 2
  !   Purpose    :  smoothing of 2D fields
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smooth2(X, niter)

    implicit none

    real(wp), intent(inout) :: X(:,:)
    integer, intent(in) :: niter

    integer :: iter, i, j, imi, ipl, im, jm
    real(wp), allocatable :: Y(:,:)


    im = size(X,1)
    jm = size(X,2)
    allocate(Y(im,jm))

    ! smoothing is repeated niter times
    do iter=1,niter

      ! save field in temporary variable
      do j=1,jm
        do i=1,im
          Y(i,j) = X(i,j)
        enddo
      enddo        

      ! smooth with 4-cell neighbors
      do i=1,im
        imi = i-1
        if (imi.eq.0) imi = im
        ipl = i+1
        if (ipl.gt.im) ipl = 1

        j = 1
        X(i,j) = 1._wp/3._wp*(Y(i,j)+Y(imi,j)+Y(ipl,j))
        do j=2,jm-1       
          X(i,j) = 0.2_wp*(Y(i,j)+Y(imi,j)+Y(ipl,j)+Y(i,j+1)+Y(i,j-1))
        enddo
        j = jm
        X(i,j) = 1._wp/3._wp*(Y(i,j)+Y(imi,j)+Y(ipl,j))
      enddo

    enddo       

    deallocate(Y)

    return

  end subroutine smooth2

end module filter
