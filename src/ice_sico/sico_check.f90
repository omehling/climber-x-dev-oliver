module sico_check

  use precision, only : wp
  use sico_state, only : sico_state_class
  use sico_grid_mod, only : sico_grid_class

  implicit none

  private
  public :: check_vel

contains


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c h e c k _ v e l
  !   Purpose    :  check velocity range 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine check_vel(st,grd, error)

    implicit none

    type(sico_state_class), intent(inout) :: st
    type(sico_grid_class), intent(in) :: grd

    logical, intent(inout) :: error

    integer :: i, j


    do i=0, grd%IMAX
      do j=0, grd%JMAX-1

        ! check for NaNs
        if (st%vx_m_sia(j,i).ne.st%vx_m_sia(j,i)) then
          print *
          print *,'ice velocity == NaN'
          print *,'(i,j)',i,j
          print *,'vx_m_sia',st%vx_m_sia(j,i)
          error = .true.
        endif
        if (st%vx_m_ssa(j,i).ne.st%vx_m_ssa(j,i)) then
          print *
          print *,'ice velocity == NaN'
          print *,'(i,j)',i,j
          print *,'vx_m_ssa',st%vx_m_ssa(j,i)
          error = .true.
        endif

      enddo
    enddo


    return

  end subroutine check_vel


end module sico_check
