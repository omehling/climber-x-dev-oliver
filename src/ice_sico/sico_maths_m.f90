!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ m a t h s _ m
!
!> @file
!!
!! Several mathematical tools used by SICOPOLIS.
!!
!! @section Copyright
!!
!! Copyright 2009-2017 Ralf Greve
!!
!! @section License
!!
!! This file is part of SICOPOLIS.
!!
!! SICOPOLIS is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! SICOPOLIS is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with SICOPOLIS.  If not, see <http://www.gnu.org/licenses/>.
!<
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------------------------------------------------------------
!> Several mathematical tools used by SICOPOLIS.
!<------------------------------------------------------------------------------
module sico_maths_m

  use sico_types_m

  implicit none

  integer, parameter :: iter_max = 1000

  public

contains

!-------------------------------------------------------------------------------
!> SOR solver for a system of linear equations lgs_a*lgs_x=lgs_b
!! [matrix storage: compressed sparse row CSR,
!! represented by arrays lgs_a_value(values), lgs_a_index (indices)
!! and lgs_a_ptr (pointers)].
!<------------------------------------------------------------------------------
  subroutine sor_sprs(lgs_a_value, lgs_a_index, lgs_a_diag_index, lgs_a_ptr, &
                      lgs_b_value, &
                      nnz, nmax, omega, eps_sor, lgs_x_value, ierr)

  implicit none

  integer,                     intent(in) :: nnz, nmax
  real(wp),                         intent(in) :: omega, eps_sor
  integer, dimension(nmax+1),  intent(in) :: lgs_a_ptr
  integer, dimension(nnz),     intent(in) :: lgs_a_index
  integer, dimension(nmax),    intent(in) :: lgs_a_diag_index
  real(wp),     dimension(nnz),     intent(in) :: lgs_a_value
  real(wp),     dimension(nmax),    intent(in) :: lgs_b_value

  integer,                    intent(out) :: ierr
  real(wp),     dimension(nmax), intent(inout) :: lgs_x_value

  integer :: iter
  integer :: nr, k
  real(wp), allocatable, dimension(:) :: lgs_x_value_prev
  real(wp)     :: b_nr
  logical      :: flag_convergence


  allocate(lgs_x_value_prev(nmax))

  iter_loop : do iter=1, iter_max

     lgs_x_value_prev = lgs_x_value

     do nr=1, nmax

        b_nr = 0.0_wp 

        do k=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
           b_nr = b_nr + lgs_a_value(k)*lgs_x_value(lgs_a_index(k))
        end do

        lgs_x_value(nr) = lgs_x_value(nr) &
                          -omega*(b_nr-lgs_b_value(nr)) &
                                /lgs_a_value(lgs_a_diag_index(nr))

     end do

     flag_convergence = .true.
     do nr=1, nmax
        if (abs(lgs_x_value(nr)-lgs_x_value_prev(nr)) > eps_sor) then
           flag_convergence = .false.
           exit
        end if
     end do

     if (flag_convergence) then
        !write(6,'(10x,a,i0)') 'sor_sprs: iter = ', iter
        ierr = 0   ! convergence criterion fulfilled
        deallocate(lgs_x_value_prev)
        return
     end if

  end do iter_loop

  !write(6,'(10x,a,i0)') 'sor_sprs: iter = ', iter
  ierr = -1   ! convergence criterion not fulfilled
  deallocate(lgs_x_value_prev)

  end subroutine sor_sprs

!-------------------------------------------------------------------------------
!> Solution of a system of linear equations Ax=b with tridiagonal matrix A.
!! @param[in]  a0       a0(j) is element A_(j,j-1) of Matrix A
!! @param[in]  a1       a1(j) is element A_(j,j)   of Matrix A
!! @param[in]  a2       a2(j) is element A_(j,j+1) of Matrix A
!! @param[in]  b        inhomogeneity vector
!! @param[in]  nrows    size of matrix A (indices run from 0 (!!!) to nrows)
!! @param[out] x        Solution vector.
!<------------------------------------------------------------------------------
  subroutine tri_sle(a0, a1, a2, x, b, nrows)

  implicit none

  integer,             intent(in)    :: nrows
  real(wp), dimension(0:*), intent(in)    :: a0, a2
  real(wp), dimension(0:*), intent(inout) :: a1, b

  real(wp), dimension(0:*), intent(out)   :: x

  real(wp), allocatable, dimension(:) :: help_x
  integer :: n

!--------  Generate an upper triangular matrix
!                      ('obere Dreiecksmatrix') --------

  do n=1, nrows
     a1(n)   = a1(n) - a0(n)/a1(n-1)*a2(n-1)
  end do

  do n=1, nrows
     b(n)    = b(n) - a0(n)/a1(n-1)*b(n-1)
     ! a0(n)  = 0.0_wp , not needed in the following, therefore
     !                   not set
  end do

!-------- Iterative solution of the new system --------

  ! x(nrows) = b(nrows)/a1(nrows)

  ! do n=nrows-1, 0, -1
  !    x(n) = (b(n)-a2(n)*x(n+1))/a1(n)
  ! end do

  allocate(help_x(0:nrows))

  help_x(0) = b(nrows)/a1(nrows)

  do n=1, nrows
     help_x(n) = b(nrows-n)/a1(nrows-n) &
                -a2(nrows-n)/a1(nrows-n)*help_x(n-1)
  end do

  do n=0, nrows
     x(n) = help_x(nrows-n)
  end do

  !       (The trick with the help_x was introduced in order to avoid
  !        the negative step in the original, blanked-out loop.)

  deallocate(help_x)

  !  WARNING: Subroutine does not check for elements of the main
  !           diagonal becoming zero. In this case it crashes even
  !           though the system may be solvable. Otherwise ok.

  end subroutine tri_sle

!-------------------------------------------------------------------------------
!> Bilinear interpolation.
!<------------------------------------------------------------------------------
  function bilinint(x1, x2, y1, y2, z11, z12, z21, z22, x, y)

  implicit none

  real(wp), intent(in) :: x1, x2, y1, y2, z11, z12, z21, z22, x, y

  real(wp) :: t, u
  real(wp) :: bilinint

  real(wp), parameter :: I = 1.0_wp

  t = (x-x1)/(x2-x1)
  u = (y-y1)/(y2-y1)

  bilinint = (I-t)*(I-u)*z11 + (I-t)*u*z12 + t*(I-u)*z21 + t*u*z22

  end function bilinint

!-------------------------------------------------------------------------------

end module sico_maths_m
!
