!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : p r e c i s i o n
!
!  Purpose : set precision 
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
module precision

  implicit none

  ! Floating point section
  
  integer,  parameter :: dp  = kind(1.d0)
  integer,  parameter :: sp  = kind(1.0)
  !INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)
  !INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)

  INTEGER, PARAMETER :: wp = dp   ! working precision

end module precision
