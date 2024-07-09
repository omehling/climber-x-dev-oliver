!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  f l a g _ u p d a t e _ g f _ g l _ c f _ m
!
!> @file
!!
!! Update of the flags for the land-terminating grounded front,
!! marine-terminating grounded front, grounding line and calving front.
!!
!! @section Copyright
!!
!! Copyright 2018-2019 Ralf Greve
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
!> Update of the flags for the land-terminating grounded front,
!! marine-terminating grounded front, grounding line and calving front.
!<------------------------------------------------------------------------------
module flag_update_gf_gl_cf_m

  use sico_types_m
  use sico_state
  use sico_grid_mod
  use sico_params, only : sico_par_class

  implicit none

  private
  public :: flag_update_gf_gl_cf

contains

!-------------------------------------------------------------------------------
!> Main subroutine of flag_update_gf_gl_cf_m:
!! Update of the flags for the land-terminating grounded front,
!! marine-terminating grounded front, grounding line and calving front.
!<------------------------------------------------------------------------------
  subroutine flag_update_gf_gl_cf(st,grd,par)

  implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class), intent(in) :: grd
  type(sico_par_class), intent(in) :: par

  integer :: i, j

  st%flag_grounded_front_a_1 = .false.
  st%flag_grounded_front_a_2 = .false.

  st%flag_grounded_front_b_1 = .false.
  st%flag_grounded_front_b_2 = .false.

  st%flag_grounding_line_1 = .false.
  st%flag_grounding_line_2 = .false.

  st%flag_calving_front_1 = .false.
  st%flag_calving_front_2 = .false.

  do i=1, grd%IMAX-1
    do j=1, grd%JMAX-1

      if ( (st%maske(j,i)==0) &   ! grounded ice point
        .and. &
        (    (st%maske(j,i+1)==1)   &   ! with
        .or.(st%maske(j,i-1)==1)   &   ! one
        .or.(st%maske(j+1,i)==1)   &   ! neighbouring
        .or.(st%maske(j-1,i)==1) ) &   ! ice-free land point
        ) &
        st%flag_grounded_front_a_1(j,i) = .true.

      if ( (st%maske(j,i)==1) &   ! ice-free land point
        .and. &
        (    (st%maske(j,i+1)==0)   &   ! with
        .or.(st%maske(j,i-1)==0)   &   ! one
        .or.(st%maske(j+1,i)==0)   &   ! neighbouring
        .or.(st%maske(j-1,i)==0) ) &   ! grounded ice point
        ) &
        st%flag_grounded_front_a_2(j,i) = .true.

      if ( (st%maske(j,i)==0) &   ! grounded ice point
        .and. &
        (    (st%maske(j,i+1)==2)   &   ! with
        .or.(st%maske(j,i-1)==2)   &   ! one
        .or.(st%maske(j+1,i)==2)   &   ! neighbouring
        .or.(st%maske(j-1,i)==2) ) &   ! ocean point
        ) &
        st%flag_grounded_front_b_1(j,i) = .true.

      if ( (st%maske(j,i)==2) &   ! ocean point
        .and. &
        (    (st%maske(j,i+1)==0)   &   ! with
        .or.(st%maske(j,i-1)==0)   &   ! one
        .or.(st%maske(j+1,i)==0)   &   ! neighbouring
        .or.(st%maske(j-1,i)==0) ) &   ! grounded ice point
        ) &
        st%flag_grounded_front_b_2(j,i) = .true.

      if (par%margin==3) then

        if ( (st%maske(j,i)==0) &   ! grounded ice point
          .and. &
          (    (st%maske(j,i+1)==3)   &   ! with
          .or.(st%maske(j,i-1)==3)   &   ! one
          .or.(st%maske(j+1,i)==3)   &   ! neighbouring
          .or.(st%maske(j-1,i)==3) ) &   ! floating ice point
          ) &
          st%flag_grounding_line_1(j,i) = .true.

        if ( (st%maske(j,i)==3) &   ! floating ice point
          .and. &
          (    (st%maske(j,i+1)==0)   &   ! with
          .or.(st%maske(j,i-1)==0)   &   ! one
          .or.(st%maske(j+1,i)==0)   &   ! neighbouring
          .or.(st%maske(j-1,i)==0) ) &   ! grounded ice point
          ) &
          st%flag_grounding_line_2(j,i) = .true.

        if ( (st%maske(j,i)==3) &   ! floating ice point
          .and. &
          (    (st%maske(j,i+1)==2)   &   ! with
          .or.(st%maske(j,i-1)==2)   &   ! one
          .or.(st%maske(j+1,i)==2)   &   ! neighbouring
          .or.(st%maske(j-1,i)==2) ) &   ! ocean point
          ) &
          st%flag_calving_front_1(j,i) = .true.

        if ( (st%maske(j,i)==2) &   ! ocean point
          .and. &
          (    (st%maske(j,i+1)==3)   &   ! with
          .or.(st%maske(j,i-1)==3)   &   ! one
          .or.(st%maske(j+1,i)==3)   &   ! neighbouring
          .or.(st%maske(j-1,i)==3) ) &   ! floating ice point
          ) &
          st%flag_calving_front_2(j,i) = .true.

      endif

    end do
  end do

  do i=1, grd%IMAX-1

    j=0

    if ( (st%maske(j,i)==1) &   ! ice-free land point
      .and. (st%maske(j+1,i)==0) &   ! with one neighbouring
      ) &                               ! grounded ice point
      st%flag_grounded_front_a_2(j,i) = .true.

    if ( (st%maske(j,i)==2) &   ! ocean point
      .and. (st%maske(j+1,i)==0) &   ! with one neighbouring
      ) &                               ! grounded ice point
      st%flag_grounded_front_b_2(j,i) = .true.

    if (par%margin==3) then

      if ( (st%maske(j,i)==2) &   ! ocean point
        .and. (st%maske(j+1,i)==3) &   ! with one neighbouring
        ) &                               ! floating ice point
        st%flag_calving_front_2(j,i) = .true.

    endif 

    j=grd%JMAX

    if ( (st%maske(j,i)==1) &   ! ice-free land point
      .and. (st%maske(j-1,i)==0) &   ! with one neighbouring
      ) &                               ! grounded ice point
      st%flag_grounded_front_a_2(j,i) = .true.

    if ( (st%maske(j,i)==2) &   ! ocean point
      .and. (st%maske(j-1,i)==0) &   ! with one neighbouring
      ) &                               ! grounded ice point
      st%flag_grounded_front_b_2(j,i) = .true.

    if (par%margin==3) then

      if ( (st%maske(j,i)==2) &   ! ocean point
        .and. (st%maske(j-1,i)==3) &   ! with one neighbouring
        ) &                               ! floating ice point
        st%flag_calving_front_2(j,i) = .true.

    endif 

  end do

  do j=1, grd%JMAX-1

    i=0

    if ( (st%maske(j,i)==1) &   ! ice-free land point
      .and. (st%maske(j,i+1)==0) &   ! with one neighbouring
      ) &                               ! grounded ice point
      st%flag_grounded_front_a_2(j,i) = .true.

    if ( (st%maske(j,i)==2) &   ! ocean point
      .and. (st%maske(j,i+1)==0) &   ! with one neighbouring
      ) &                               ! grounded ice point
      st%flag_grounded_front_b_2(j,i) = .true.

    if (par%margin==3) then

      if ( (st%maske(j,i)==2) &   ! ocean point
        .and. (st%maske(j,i+1)==3) &   ! with one neighbouring
        ) &                               ! floating ice point
        st%flag_calving_front_2(j,i) = .true.

    endif 

    i=grd%IMAX

    if ( (st%maske(j,i)==1) &   ! ice-free land point
      .and. (st%maske(j,i-1)==0) &   ! with one neighbouring
      ) &                               ! grounded ice point
      st%flag_grounded_front_a_2(j,i) = .true.

    if ( (st%maske(j,i)==2) &   ! ocean point
      .and. (st%maske(j,i-1)==0) &   ! with one neighbouring
      ) &                               ! grounded ice point
      st%flag_grounded_front_b_2(j,i) = .true.

    if (par%margin==3) then

      if ( (st%maske(j,i)==2) &   ! ocean point
        .and. (st%maske(j,i-1)==3) &   ! with one neighbouring
        ) &                               ! floating ice point
        st%flag_calving_front_2(j,i) = .true.

    endif   

  end do

  end subroutine flag_update_gf_gl_cf

!-------------------------------------------------------------------------------

end module flag_update_gf_gl_cf_m
!
