!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s m o o t h _ a t m _ m o d
!
!  Purpose : smoothing functions for atmosphere
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
module smooth_atm_mod

  use atm_params, only : wp

  implicit none

  private
  public :: smooth2, smooth2_m, smooth2eq, zona, zofil

contains

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


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s m o o t h 2 _ m
  !   Purpose    :  smoothing of 2D fields
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smooth2_m(X,niter)

    implicit none

    real(wp), intent(inout) :: X(:,:)
    integer, intent(in) :: niter

    integer :: iter, i, j, imi, ipl, im, jm
    real(wp), allocatable :: Y(:,:)


    im = size(X,1)
    jm = size(X,2)
    allocate(Y(im,jm))

    do iter=1,niter

      do j=1,jm
        do i=1,im
          Y(i,j) = X(i,j)
        enddo
      enddo        

      !$omp parallel do private(i,j,imi,ipl) 
      do i=1,im
        imi = i-1
        if (imi.eq.0) imi = im
        ipl = i+1
        if (ipl.gt.im) ipl = 1

        j = 1
        X(i,j) = 1._wp/3._wp*(Y(i,j)+Y(imi,j)+Y(ipl,j))
        do j=2,jm-1       
          X(i,j) = 0.4_wp*Y(i,j) + 0.15_wp*(Y(imi,j)+Y(ipl,j)+Y(i,j+1)+Y(i,j-1))
        enddo
        j = jm
        X(i,j) = 1._wp/3._wp*(Y(i,j)+Y(imi,j)+Y(ipl,j))
      enddo
      !$omp end parallel do

    enddo       

    deallocate(Y)

    return

  end subroutine smooth2_m


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s m o o t h 2 e q
  !   Purpose    :  smoothing of 2D fields around the equator
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smooth2eq(X, nj, niter)

    implicit none

    real(wp), intent(inout) :: X(:,:)
    integer, intent(in) :: nj
    integer, intent(in) :: niter

    integer :: iter, i, j, imi, ipl, im, jm
    real(wp), allocatable :: Y(:,:)


    im = size(X,1)
    jm = size(X,2)
    allocate(Y(im,jm))

    do iter=1,niter

      do j=1,jm
        do i=1,im
          Y(i,j) = X(i,j)
        enddo
      enddo        

      do i=1,im
        imi = i-1
        if (imi.eq.0) imi = im
        ipl = i+1
        if (ipl.gt.im) ipl = 1

        do j=jm/2-nj,jm/2+1+nj       
          X(i,j) = 0.2_wp*(Y(i,j)+Y(imi,j)+Y(ipl,j)+Y(i,j+1)+Y(i,j-1))
        enddo
      enddo

    enddo       

    deallocate(Y)

    return

  end subroutine smooth2eq


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  z o n a
  !   Purpose    :  compute average of field over latitudinal belt
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function zona(X, cost, j1, j2)

    implicit none

    real(wp), intent(in) :: X(:)
    real(wp), intent(in) :: cost(:)
    integer, intent(in) :: j1, j2

    integer :: j
    real(wp) :: zona, weight


    zona = 0._wp
    weight = 0._wp

    do j=j1,j2
      zona = zona+cost(j)*X(j)  
      weight = weight+cost(j)
    enddo

    zona = zona/weight

    return

  end function zona


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  z o f i l
  !   Purpose    :  zonal filtering nite times for lat belt with index j 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine zofil(x,jf,nite)

    implicit none

    real(wp), intent(inout) :: X(:,:)
    integer, intent(in) :: jf, nite

    integer :: n, i, ipl, imi, im, jm
    real(wp), allocatable :: Y(:,:)


    im = size(X,1)
    jm = size(X,2)
    allocate(Y(im,jm))

    do n=1,nite

      do i=1,im
        Y(i,jf)=X(i,jf)
      enddo

      do i=1,im
        ipl=i+1
        if (ipl.gt.im) ipl=1
        imi=i-1
        if (imi.lt.1) imi=im 
        X(i,jf)=(Y(imi,jf)+Y(i,jf)+Y(ipl,jf))/3._wp
      enddo

    enddo

    deallocate(Y)

    return

  end subroutine zofil

end module smooth_atm_mod
