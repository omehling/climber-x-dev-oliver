module smb_out

  use smb_def, only : smb_class

  implicit none

  private
  public :: smb_diag, smb_diag_init


contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s m b _ d i a g _ i n i t
  ! Purpose  :  Initialize netcdf output for smb
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_diag_init(smb)

    implicit none

    type(smb_class) :: smb

   return

  end subroutine smb_diag_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s m b _ d i a g
  !   Purpose    :  sea ice diagnostics
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_diag(smb)

    implicit none

    type(smb_class) :: smb

   return

  end subroutine smb_diag

end module smb_out
