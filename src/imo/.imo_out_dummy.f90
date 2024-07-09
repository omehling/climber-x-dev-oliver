module imo_out

  use imo_def, only : imo_class

  implicit none

  private
  public :: imo_diag, imo_diag_init


contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  i m o _ d i a g _ i n i t
  ! Purpose  :  Initialize netcdf output for imo
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine imo_diag_init(imo)

    implicit none

    type(imo_class) :: imo

   return

  end subroutine imo_diag_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i m o _ d i a g
  !   Purpose    :  sea ice diagnostics
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine imo_diag(imo)

    implicit none

    type(imo_class) :: imo

   return

  end subroutine imo_diag

end module imo_out
