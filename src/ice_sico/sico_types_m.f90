!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ t y p e s _ m
!
!> @file
!!
!! Declarations of kind types for SICOPOLIS.
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
!> Declarations of kind types for SICOPOLIS.
!<------------------------------------------------------------------------------
module sico_types_m

implicit none
save

integer, parameter :: sp  = kind(1.0)              !< Single-precision reals
integer, parameter :: dp  = kind(1.0d0)            !< Double-precision reals

integer, parameter :: wp = dp

end module sico_types_m
!
