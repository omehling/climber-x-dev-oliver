!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  e n t h _ t e m p _ o m e g a _ m
!
!> @file
!!
!! Conversion from temperature (temp) and water content (omega) to enthalpy
!! (enth) and vice versa.
!!
!! @section Copyright
!!
!! Copyright 2013-2017 Ralf Greve, Heinz Blatter
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
!> Conversion from temperature (temp) and water content (omega) to enthalpy
!! (enth) and vice versa.
!<------------------------------------------------------------------------------
module enth_temp_omega_m

use sico_types_m

implicit none
save

!> c_int_table:     Temperature integral of the specific heat of ice.
!>                  Index is temperature in deg C.
   real(wp), dimension(-256:255), private       :: c_int_table

!> c_int_inv_table: Inverse of the temperature integral of the specific heat
!>                  of ice. Index is enthalpy in J/kg (zero for 0 deg C).
   real(wp), dimension(-524288:524287), private :: c_int_inv_table

!> n_temp_min: Lower index limit of properly defined values in c_int_table
!>             (n_temp_min >= -256).
   integer, private :: n_temp_min

!> n_temp_max: Upper index limit of properly defined values in c_int_table
!>             (n_temp_max <= 255).
   integer, private :: n_temp_max

!> n_enth_min: Lower index limit of properly defined values in c_int_inv_table
!>             (n_enth_min >= -524288).
   integer, private :: n_enth_min

!> n_enth_max: Upper index limit of properly defined values in c_int_inv_table
!>             (n_enth_max <= 524287).
   integer, private :: n_enth_max

!> L: Latent heat of ice
   real(wp), private :: L

!> L_inv: Inverse of the latent heat of ice
   real(wp), private :: L_inv

public  :: calc_c_int_table, calc_c_int_inv_table
public  :: enth_fct_temp_omega, temp_fct_enth, omega_fct_enth
private :: c_int_val, c_int_inv_val

contains

!-------------------------------------------------------------------------------
!> Computation of the temperature integral of the specific heat of ice as a
!! table (c_int_table). Further, definition of the latent heat of ice.
!<------------------------------------------------------------------------------
subroutine calc_c_int_table(c_table, n_tmp_min, n_tmp_max, L_val)

implicit none

integer,                             intent(in) :: n_tmp_min, n_tmp_max
real(wp), dimension(n_tmp_min:n_tmp_max), intent(in) :: c_table
real(wp), optional,                       intent(in) :: L_val

integer :: n
real(wp)     :: c_int_zero

!-------- Initialisation --------

c_int_table = 0.0_wp

n_temp_min = n_tmp_min
n_temp_max = n_tmp_max

if ((n_temp_min <= -256).or.(n_temp_max >= 255)) &
   stop ' >>> calc_c_int_table: Temperature indices out of allowed range!'

!-------- Numerical integration with the trapezoidal rule (spacing
!         of data in c_table and c_int_table assumed to be 1 deg C) --------

do n=n_temp_min+1, n_temp_max
   c_int_table(n) = c_int_table(n-1) + 0.5_wp*(c_table(n-1)+c_table(n))
                    ! that's the real stuff
end do

do n=n_temp_max+1, 255
   c_int_table(n) = c_int_table(n_temp_max)   ! dummy values
end do

!-------- Shift of the zero level to 0 deg C --------

c_int_zero = c_int_table(0)

do n=-256, 255
   c_int_table(n) = c_int_table(n) - c_int_zero
end do

!-------- Further stuff: Latent heat of ice --------

if ( present(L_val) ) then
   L = L_val
else
   L = 3.35e+05_wp   ! in J/kg
end if

L_inv = 1.0_wp/L

end subroutine calc_c_int_table

!-------------------------------------------------------------------------------
!> Computation of the inverse of the temperature integral of the specific heat
!! of ice as a table (c_int_inv_table).
!<------------------------------------------------------------------------------
subroutine calc_c_int_inv_table()

implicit none

integer :: n
integer :: n_temp_1, n_temp_2
real(wp)     :: enth_min, enth_max
real(wp)     :: enth_val, enth_1, enth_2

!-------- Initialisation --------

c_int_inv_table = 0.0_wp

enth_min = c_int_val(real(n_temp_min,wp))
enth_max = c_int_val(real(n_temp_max,wp))

n_enth_min = ceiling(enth_min)
n_enth_max = floor(enth_max)

if ((n_enth_min <= -524288).or.(n_enth_max >= 524287)) &
   stop ' >>> calc_c_int_inv_table: Enthalpy indices out of allowed range!'

!-------- Linear interpolation between adjacent enthalpy values --------

n_temp_1 = n_temp_min
n_temp_2 = n_temp_min+1

do n=n_enth_min, n_enth_max

   enth_val = real(n,wp)

   do

      if ((n_temp_1 > n_temp_max).or.(n_temp_2 > n_temp_max)) &
         stop ' >>> calc_c_int_inv_table: Temperature indices out of allowed range!'

      enth_1 = c_int_val(real(n_temp_1,wp))
      enth_2 = c_int_val(real(n_temp_2,wp))

      if ( (enth_1 <= enth_val).and.(enth_2 >= enth_val) ) exit

      n_temp_1 = n_temp_1+1
      n_temp_2 = n_temp_2+1

   end do

   c_int_inv_table(n) = real(n_temp_1,wp) &
                        + (real(n_temp_2,wp)-real(n_temp_1,wp)) &
                          * (enth_val-enth_1)/(enth_2-enth_1)
                                              ! Linear interpolation
end do

do n=-524288, n_enth_min-1
   c_int_inv_table(n) = c_int_inv_table(n_enth_min)   ! dummy values
end do

do n=n_enth_max+1, 524287
   c_int_inv_table(n) = c_int_inv_table(n_enth_max)   ! dummy values
end do

end subroutine calc_c_int_inv_table

!-------------------------------------------------------------------------------
!> Temperature integral of the specific heat of ice
!! (enthalpy as function of temperature).
!<------------------------------------------------------------------------------
function c_int_val(temp_val)

implicit none

real(wp)             :: c_int_val

real(wp), intent(in) :: temp_val

integer :: n_temp_1, n_temp_2

n_temp_1 = floor(temp_val)
n_temp_1 = max(min(n_temp_1, n_temp_max-1), n_temp_min)
n_temp_2 = n_temp_1 + 1

c_int_val = c_int_table(n_temp_1) &
            + (c_int_table(n_temp_2)-c_int_table(n_temp_1)) &
              * (temp_val-real(n_temp_1,wp))   ! Linear interpolation

end function c_int_val

!-------------------------------------------------------------------------------
!> Inverse function of c_int_val (temperature as function of enthalpy).
!<------------------------------------------------------------------------------
function c_int_inv_val(enth_val)

implicit none

real(wp)             :: c_int_inv_val

real(wp), intent(in) :: enth_val

integer :: n_enth_1, n_enth_2

n_enth_1 = floor(enth_val)
n_enth_1 = max(min(n_enth_1, n_enth_max-1), n_enth_min)
n_enth_2 = n_enth_1 + 1

c_int_inv_val = c_int_inv_table(n_enth_1) &
                + (c_int_inv_table(n_enth_2)-c_int_inv_table(n_enth_1)) &
                  * (enth_val-real(n_enth_1,wp))   ! Linear interpolation

end function c_int_inv_val

!-------------------------------------------------------------------------------
!> Enthalpy as a function of temperature and water content.
!<------------------------------------------------------------------------------
function enth_fct_temp_omega(temp_val, omega_val)

implicit none

real(wp)              :: enth_fct_temp_omega

real(wp), intent(in)  :: temp_val, omega_val

enth_fct_temp_omega = c_int_val(temp_val) + L*omega_val

end function enth_fct_temp_omega

!-------------------------------------------------------------------------------
!> Temperature as a function of enthalpy.
!<------------------------------------------------------------------------------
function temp_fct_enth(enth_val, temp_m_val)

implicit none

real(wp)             :: temp_fct_enth

real(wp), intent(in) :: enth_val
real(wp), intent(in) :: temp_m_val

real(wp) :: enth_i

enth_i = c_int_val(temp_m_val)   ! Enthalpy of pure ice at the melting point

if (enth_val < enth_i) then   ! cold ice
   temp_fct_enth = c_int_inv_val(enth_val)
else   ! temperate ice
   temp_fct_enth = temp_m_val
end if

end function temp_fct_enth

!-------------------------------------------------------------------------------
!> Water content as a function of enthalpy.
!<------------------------------------------------------------------------------
function omega_fct_enth(enth_val, temp_m_val)

implicit none

real(wp)             :: omega_fct_enth

real(wp), intent(in) :: enth_val
real(wp), intent(in) :: temp_m_val

real(wp) :: enth_i

enth_i = c_int_val(temp_m_val)   ! Enthalpy of pure ice at the melting point

omega_fct_enth = max((enth_val-enth_i)*L_inv, 0.0_wp)

end function omega_fct_enth

!-------------------------------------------------------------------------------

end module enth_temp_omega_m
!
