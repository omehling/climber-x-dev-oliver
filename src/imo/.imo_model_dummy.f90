module imo_model

  use precision, only : wp, dp
  use coord, only : grid_class
  use imo_def, only : imo_class

  implicit none

  private
  public :: imo_init, imo_update, imo_end, imo_write_restart

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i m o _ u p d a t e
  !   Purpose    :  update ice melt to ocean
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine imo_update(imo)

  !$  use omp_lib

  implicit none

  type(imo_class) :: imo

  return

  end subroutine imo_update
      

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i m o _ i n i t
  !   Purpose    :  initialize ice melt to ocean 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine imo_init(imo,grid,cmn_grid)

    implicit none

    type(imo_class) :: imo
    type(grid_class), intent(in) :: grid 
    type(grid_class), intent(in) :: cmn_grid 

  return

  end subroutine imo_init
      
     
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  i m o _ e n d 
  ! Purpose  :  end imo
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine imo_end(imo)

    implicit none

    type(imo_class) :: imo

    return

  end subroutine imo_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  i m o _ w r i t e _ r e s t a r t
  ! Purpose  :  Write restart file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine imo_write_restart(fnm,imo)

    implicit none

    character (len=*) :: fnm
    type(imo_class) :: imo

  end subroutine imo_write_restart

end module imo_model
