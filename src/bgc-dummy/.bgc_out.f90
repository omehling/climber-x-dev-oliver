module bgc_out

  use precision, only : wp
  use bgc_def, only : bgc_class

  implicit none

  private
  public :: bgc_diag_init, bgc_diag

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b g c _ d i a g _ i n i t
  ! Purpose  :  Initialize netcdf output for ocean biogeochemistry
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bgc_diag_init(ni, nj, nk, z)

    implicit none

    integer, intent(in) :: ni, nj, nk
    real(wp), intent(in) :: z(:)

   return

  end subroutine bgc_diag_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  b g c _ d i a g
  !   Purpose    :  ocean biogeochemistry diagnostics
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bgc_diag(bgc)

    implicit none

    type(bgc_class), intent(in) :: bgc

   return

  end subroutine bgc_diag

end module bgc_out
