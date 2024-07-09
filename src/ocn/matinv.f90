!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : m a t i n v _ m o d
!
!  Purpose : matrix inversion
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
module matinv_mod
  ! includes subroutines to solve a set of n linear equations by direct inversion
  ! solves amat*x = rhs, putting result into rhs
  ! split in to two parts to invert and multiply separately
  ! to save a lot of cpu if matrix is constant in time

  use precision, only : wp

  implicit none

  private
  public :: matinv, matmult

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  m a t i n v
  !   Purpose    :  invert
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine matinv(nisl, &
                    amat)
    
    implicit none

    integer, intent(in) :: nisl
    real(wp), intent(inout) :: amat(:,:)

    integer i, j, k

    ! elimination

    do i=1,nisl-1
       do j=i+1,nisl
          do k=i+1,nisl
             amat(j,k) = amat(i,i)*amat(j,k) - amat(j,i)*amat(i,k)
             !if (j.eq.nisl .and. k.eq.nisl) then
             !  print *,i,amat(i,i),amat(j,k),amat(j,i),amat(i,k)
             !  print *,amat(j,k),amat(i,i)*amat(j,k),amat(j,i)*amat(i,k)
             !endif
          enddo
       enddo
    enddo

   return

  end subroutine matinv


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  m a t m u l t
  !   Purpose    :  multiply
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine matmult(nisl,amat, &
                     rhs)

    implicit none

    integer, intent(in) :: nisl
    real(wp), intent(in) :: amat(:,:)
    real(wp), intent(inout) :: rhs(:)

    integer i, j


    do i=1,nisl-1
       do j=i+1,nisl
          rhs(j) = amat(i,i)*rhs(j) - amat(j,i)*rhs(i)
       enddo
    enddo

    ! back substitution
    !print *,'nisl',nisl
    !print *,'rhs(nisl)',rhs(nisl)
    !print *,'amat(nisl,nisl)',amat(nisl,nisl)

    rhs(nisl) = rhs(nisl)/amat(nisl,nisl)
    do i=nisl-1,1,-1
       do j=i+1,nisl
          rhs(i) = rhs(i) - amat(i,j)*rhs(j)
       enddo
       rhs(i) = rhs(i)/amat(i,i)
    enddo

!c     print*,(rhs(i),i=1,nisl)
!c     print*
!c     print*,cmat(1,1)*rhs(1) + cmat(1,2)*rhs(2) + cmat(1,3)*rhs(3) &
!c           + cmat(1,4)*rhs(4) - orhs(1), &
!c           cmat(2,1)*rhs(1) + cmat(2,2)*rhs(2) + cmat(2,3)*rhs(3) &
!c           + cmat(2,4)*rhs(4) - orhs(2), &
!c           cmat(3,1)*rhs(1) + cmat(3,2)*rhs(2) + cmat(3,3)*rhs(3) &
!c           + cmat(3,4)*rhs(4) - orhs(3)
!c     print*,((amat(i,j),i=1,nisl),j=1,nisl)

   return

  end subroutine matmult

end module matinv_mod
