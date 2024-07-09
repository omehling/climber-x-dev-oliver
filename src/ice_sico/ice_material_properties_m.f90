!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  i c e _ m a t e r i a l _ p r o p e r t i e s _ m
!
!> @file
!!
!! Material properties of ice:
!! Rate factor, heat conductivity, specific heat (heat capacity),
!! creep function, viscosity.
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
!> Material properties of ice:
!! Rate factor, heat conductivity, specific heat (heat capacity),
!! creep function, viscosity.
!<------------------------------------------------------------------------------
module ice_material_properties_m

use sico_types_m, only : wp
use sico_params, only : RF, KAPPA, C, RHO, RHO_C, KAPPA_C, C_C, R_T

implicit none
save

!> n_temp_min: Lower index limit of properly defined values in RF, KAPPA and C
!>             (n_temp_min >= -256).
   integer, private :: n_temp_min = -190

!> n_temp_max: Upper index limit of properly defined values in RF, KAPPA and C
!>             (n_temp_max <= 255).
   integer, private :: n_temp_max = 10

   real(wp), parameter :: GR_SIZE = 1.0e-03_wp
!                         Average grain size (in m; only for FLOW_LAW==2)

   real(wp), parameter :: SIGMA_RES = 1.0e+04_wp
!                         Residual stress (finite-viscosity contribution)
!                         in the creep response function
!                         (in Pa; only for FLOW_LAW==1, 2, 3 and FIN_VISC==2)

private
public :: ratefac_c, ratefac_t, ratefac_c_t, kappa_val, c_val, &
          viscosity, creep

contains


!-------------------------------------------------------------------------------
!> Rate factor for cold ice:
!! Linear interpolation of tabulated values in RF(.).
!<------------------------------------------------------------------------------
function ratefac_c(temp_val, temp_m_val)

use sico_types_m

implicit none
real(wp)             :: ratefac_c
real(wp), intent(in) :: temp_val, temp_m_val

integer :: n_temp_1, n_temp_2
real(wp)     :: temp_h_val

temp_h_val = temp_val-temp_m_val

n_temp_1 = floor(temp_h_val)
n_temp_1 = max(min(n_temp_1, n_temp_max-1), n_temp_min)
n_temp_2 = n_temp_1 + 1

ratefac_c = RF(n_temp_1) &
            + (RF(n_temp_2)-RF(n_temp_1)) &
              * (temp_h_val-real(n_temp_1,wp))   ! Linear interpolation

end function ratefac_c

!-------------------------------------------------------------------------------
!> Rate factor for temperate ice.
!<------------------------------------------------------------------------------
function ratefac_t(omega_val)

use sico_types_m

implicit none
real(wp)             :: ratefac_t
real(wp), intent(in) :: omega_val

ratefac_t = RF(0)*(1.0_wp+R_T*(omega_val))

end function ratefac_t

!-------------------------------------------------------------------------------
!> Rate factor for cold and temperate ice:
!! Combination of ratefac_c and ratefac_t (only for the enthalpy method).
!<------------------------------------------------------------------------------
function ratefac_c_t(temp_val, omega_val, temp_m_val)

use sico_types_m

implicit none
real(wp)             :: ratefac_c_t
real(wp), intent(in) :: temp_val, temp_m_val, omega_val

integer :: n_temp_1, n_temp_2
real(wp)     :: temp_h_val

temp_h_val = temp_val-temp_m_val

n_temp_1 = floor(temp_h_val)
n_temp_1 = max(min(n_temp_1, n_temp_max-1), n_temp_min)
n_temp_2 = n_temp_1 + 1

ratefac_c_t = ( RF(n_temp_1) &
                + (RF(n_temp_2)-RF(n_temp_1)) &
                  * (temp_h_val-real(n_temp_1,wp)) ) &
              * (1.0_wp+R_T*(omega_val))

end function ratefac_c_t

!-------------------------------------------------------------------------------
!> Heat conductivity of ice:
!! Linear interpolation of tabulated values in KAPPA(.).
!<------------------------------------------------------------------------------
function kappa_val(temp_val)

use sico_types_m

implicit none
real(wp)             :: kappa_val
real(wp), intent(in) :: temp_val

integer :: n_temp_1, n_temp_2
real(wp)     :: kappa_ice

!-------- Heat conductivity of pure ice --------

n_temp_1 = floor(temp_val)
n_temp_1 = max(min(n_temp_1, n_temp_max-1), n_temp_min)
n_temp_2 = n_temp_1 + 1

#if defined(FRAC_DUST)
kappa_ice = KAPPA(n_temp_1) &
            + (KAPPA(n_temp_2)-KAPPA(n_temp_1)) &
              * (temp_val-real(n_temp_1,wp))
#else
kappa_val = KAPPA(n_temp_1) &
            + (KAPPA(n_temp_2)-KAPPA(n_temp_1)) &
              * (temp_val-real(n_temp_1,wp))
#endif

!-------- If dust is present (polar caps of Mars):
!         Heat conductivity of ice-dust mixture --------

#if defined(FRAC_DUST)
kappa_val = (1.0_wp-FRAC_DUST)*kappa_ice + FRAC_DUST*KAPPA_C
#endif

end function kappa_val

!-------------------------------------------------------------------------------
!> Specific heat of ice:
!! Linear interpolation of tabulated values in C(.).
!<------------------------------------------------------------------------------
function c_val(temp_val)

use sico_types_m

implicit none
real(wp)             :: c_val
real(wp), intent(in) :: temp_val

integer :: n_temp_1, n_temp_2
real(wp)     :: c_ice

!-------- Specific heat of pure ice --------

n_temp_1 = floor(temp_val)
n_temp_1 = max(min(n_temp_1, n_temp_max-1), n_temp_min)
n_temp_2 = n_temp_1 + 1

#if defined(FRAC_DUST)
c_ice = C(n_temp_1) &
        + (C(n_temp_2)-C(n_temp_1)) &
          * (temp_val-real(n_temp_1,wp))
#else
c_val = C(n_temp_1) &
        + (C(n_temp_2)-C(n_temp_1)) &
          * (temp_val-real(n_temp_1,wp))
#endif

!-------- If dust is present (polar caps of Mars):
!         Specific heat of ice-dust mixture --------

#if defined(FRAC_DUST)
c_val = rho_inv * ( (1.0_wp-FRAC_DUST)*RHO*c_ice + FRAC_DUST*RHO_C*C_C )
#endif

end function c_val

!-------------------------------------------------------------------------------
!> Creep response function for ice.
!<------------------------------------------------------------------------------
function creep(fin_visc,flow_law,sigma_val)

use sico_types_m

implicit none

real(wp)             :: creep
integer, intent(in) :: fin_visc
integer, intent(in) :: flow_law
real(wp), intent(in) :: sigma_val

real(wp), parameter :: sm_coeff_1 = 8.5112e-15_wp, &   ! s^-1 Pa^-1
                       sm_coeff_2 = 8.1643e-25_wp, &   ! s^-1 Pa^-3
                       sm_coeff_3 = 9.2594e-12_wp      ! Pa^-2


if (fin_visc==1) then

if (flow_law==1) then

creep = sigma_val*sigma_val
!       Glen's flow law (n=3)

else if (flow_law==2) then

creep = sigma_val**0.8_wp * GR_SIZE**(-1.4_wp)
!       Goldsby-Kohlstedt flow law (n=1.8, d=1.4)

else if (flow_law==3) then

creep = sigma_val*sigma_val*sigma_val
!       Durham's flow law (n=4)

endif

else if (fin_visc==2) then

if (flow_law==1) then

creep = sigma_val*sigma_val + SIGMA_RES*SIGMA_RES
!       Glen's flow (n=3) with additional finite viscosity

else if (flow_law==2) then

creep = (sigma_val**0.8_wp + SIGMA_RES**0.8_wp) * GR_SIZE**(-1.4_wp)
!       Goldsby-Kohlstedt flow law (n=1.8, d=1.4)
!       with additional finite viscosity

else if (flow_law==3) then

creep = sigma_val*sigma_val*sigma_val + SIGMA_RES*SIGMA_RES*SIGMA_RES
!       Durham's flow law (n=4) with additional finite viscosity

endif

endif

if (flow_law==4) then

creep = sm_coeff_1 &
        + sm_coeff_2 * (sigma_val*sigma_val) &
          * (1.0_wp + sm_coeff_3 * (sigma_val*sigma_val))
!       Smith-Morland (polynomial) flow law, normalised to a dimensionless
!       rate factor with A(-10C)=1.

endif

end function creep

!-------------------------------------------------------------------------------
!> Ice viscosity as a function of the effective strain rate and the temperature
!! (in cold ice) or the water content (in temperate ice) or both (for the
!! enthalpy method).
!<------------------------------------------------------------------------------
function viscosity(fin_visc, flow_law, de_val, temp_val, temp_m_val, omega_val, enh_val, &
                   i_flag_cold_temp)

use sico_types_m

implicit none

real(wp)                 :: viscosity
integer, intent(in) :: fin_visc
integer, intent(in) :: flow_law
real(wp)    , intent(in) :: de_val
real(wp)    , intent(in) :: temp_val, temp_m_val
real(wp)    , intent(in) :: omega_val
real(wp)    , intent(in) :: enh_val
integer, intent(in) :: i_flag_cold_temp

real(wp) :: ratefac_val
real(wp) :: de_val_m

real(wp) :: n_power_law, inv_n_power_law, n_grain_size

real(wp), parameter :: de_min = 1.0e-30_wp   ! minimum value for the
                                             ! effective strain rate

real(wp), parameter :: sm_coeff_1 = 8.5112e-15_wp, &   ! s^-1 Pa^-1
                       sm_coeff_2 = 8.1643e-25_wp, &   ! s^-1 Pa^-3
                       sm_coeff_3 = 9.2594e-12_wp      ! Pa^-2

!-------- Rate factor and effective strain rate --------

if (i_flag_cold_temp == 0) then   ! cold ice
   ratefac_val = ratefac_c(temp_val, temp_m_val)
else if (i_flag_cold_temp == 1) then   ! temperate ice
   ratefac_val = ratefac_t(omega_val)
else   ! enthalpy method
   ratefac_val = ratefac_c_t(temp_val, omega_val, temp_m_val)
end if

de_val_m = max(de_val, de_min)

if (fin_visc==1) then

if (flow_law==1) then

!-------- Glen's flow law (n=3) --------

inv_n_power_law = 0.333333333333333_wp   ! 1/3

viscosity = 0.5_wp * de_val_m**(inv_n_power_law-1.0_wp) &
                   * (enh_val*ratefac_val)**(-inv_n_power_law)

else if (flow_law==2) then

!-------- Goldsby-Kohlstedt flow law (n=1.8, d=1.4) --------

inv_n_power_law = 0.555555555555555_wp   ! 1/1.8
n_grain_size    = 1.4_wp

viscosity = 0.5_wp * de_val_m**(inv_n_power_law-1.0_wp) &
                   * GR_SIZE**(n_grain_size*inv_n_power_law) &
                   * (enh_val*ratefac_val)**(-inv_n_power_law)

else if (flow_law==3) then

!-------- Durham's flow law (n=4) --------

inv_n_power_law = 0.25_wp   ! 1/4

viscosity = 0.5_wp * de_val_m**(inv_n_power_law-1.0_wp) &
                   * (enh_val*ratefac_val)**(-inv_n_power_law)

endif

else if (fin_visc==2) then

if (flow_law==1) then

!-------- Glen's flow (n=3) with additional finite viscosity --------

n_power_law = 3.0_wp

viscosity = visc_iter(de_val_m, ratefac_val, enh_val, n_power_law, SIGMA_RES)

else if (flow_law==2) then

!-------- Goldsby-Kohlstedt flow law (n=1.8, d=1.4)
!                           with additional finite viscosity --------

n_power_law  = 1.8_wp
n_grain_size = 1.4_wp

write(6,'(a)') ' >>> viscosity: Computation of the viscosity as a function'
write(6,'(a)') '           of the effective strain rate not yet implemented'
write(6,'(a)') '           for grain-size-dependent finite viscosity flow laws!'
stop 

else if (flow_law==3) then

!-------- Durham's flow law (n=4) with additional finite viscosity --------

n_power_law = 4.0_wp

viscosity = visc_iter(de_val_m, ratefac_val, enh_val, n_power_law, SIGMA_RES)

endif

endif

if (flow_law==4) then

!-------- Smith-Morland (polynomial) flow law --------

viscosity = visc_iter_sm(de_val_m, ratefac_val, enh_val, &
                         sm_coeff_1, sm_coeff_2, sm_coeff_3)

endif

end function viscosity

!-------------------------------------------------------------------------------
!> Iterative computation of the viscosity by solving equation (4.28)
!! by Greve and Blatter (Springer, 2009).
!<------------------------------------------------------------------------------
function visc_iter(de_val_m, ratefac_val, enh_val, n_power_law, sigma_res)

use sico_types_m

implicit none

real(wp)             :: visc_iter
real(wp), intent(in) :: de_val_m
real(wp), intent(in) :: ratefac_val, enh_val
real(wp), intent(in) :: n_power_law, sigma_res


integer :: n
integer :: max_iters
real(wp)     :: visc_val, res
logical      :: flag_rescheck1, flag_rescheck2

real(wp), parameter :: eps = 1.0e-05_wp   ! convergence parameter

!-------- Determination of the order of magnitude --------

visc_val = 1.0e+10_wp   ! initial guess (very low value)

flag_rescheck1 = .false.
n              = 0
max_iters      = 30

do while ((.not.flag_rescheck1).and.(n <= max_iters))

   n = n+1

   res = fct_visc(de_val_m, ratefac_val, enh_val, visc_val, &
                  n_power_law, sigma_res)

   if (res < 0.0_wp) then
      visc_val = 10.0_wp*visc_val
   else
      flag_rescheck1 = .true.
   end if

end do

!-------- Newton's method --------

if (flag_rescheck1) then
   ! only if order of magnitude could be detected successfully

   flag_rescheck2 = .false.
   n              = 0
   max_iters      = 1000

   do while ((.not.flag_rescheck2).and.(n <= max_iters))

      n = n+1

      visc_val = visc_val &
                 - res &
                   /fct_visc_deriv(de_val_m, ratefac_val, enh_val, visc_val, &
                                   n_power_law, sigma_res)

      res = fct_visc(de_val_m, ratefac_val, enh_val, visc_val, &
                     n_power_law, sigma_res)

      if (abs(res) < eps) then 
         flag_rescheck2 = .true. 
      end if

   end do

end if

visc_iter = visc_val

end function visc_iter

!-------------------------------------------------------------------------------
!> Viscosity polynomial
!! [equation (4.28) by Greve and Blatter (Springer, 2009)].
!<------------------------------------------------------------------------------
function fct_visc(de_val_m, ratefac_val, enh_val, visc_val, &
                  n_power_law, sigma_res)

use sico_types_m

implicit none

real(wp)             :: fct_visc
real(wp), intent(in) :: de_val_m
real(wp), intent(in) :: ratefac_val, enh_val
real(wp), intent(in) :: visc_val
real(wp), intent(in) :: n_power_law, sigma_res

fct_visc = 2.0_wp**n_power_law &
             *enh_val*ratefac_val &
             *de_val_m**(n_power_law-1.0_wp) &
             *visc_val**n_power_law &
          + 2.0_wp*enh_val*ratefac_val &
             *sigma_res**(n_power_law-1.0_wp) &
             *visc_val &
          - 1.0_wp

end function fct_visc

!-------------------------------------------------------------------------------
!> Derivative of the viscosity polynomial
!! [equation (4.28) by Greve and Blatter (Springer, 2009)].
!<------------------------------------------------------------------------------
function fct_visc_deriv(de_val_m, ratefac_val, enh_val, visc_val, &
                        n_power_law, sigma_res)

use sico_types_m

implicit none

real(wp)             :: fct_visc_deriv
real(wp), intent(in) :: de_val_m
real(wp), intent(in) :: ratefac_val, enh_val
real(wp), intent(in) :: visc_val
real(wp), intent(in) :: n_power_law, sigma_res

fct_visc_deriv = 2.0_wp**n_power_law*n_power_law &
                   *enh_val*ratefac_val &
                   *de_val_m**(n_power_law-1.0_wp) &
                   *visc_val**(n_power_law-1.0_wp) &
                 + 2.0_wp*enh_val*ratefac_val &
                   *sigma_res**(n_power_law-1.0_wp)
         
end function fct_visc_deriv


!-------------------------------------------------------------------------------
!> Iterative computation of the viscosity by solving equation (4.33)
!! [analogous to (4.28)] by Greve and Blatter (Springer, 2009).
!<------------------------------------------------------------------------------
function visc_iter_sm(de_val_m, ratefac_val, enh_val, &
                      sm_coeff_1, sm_coeff_2, sm_coeff_3)

use sico_types_m

implicit none

real(wp)             :: visc_iter_sm
real(wp), intent(in) :: de_val_m
real(wp), intent(in) :: ratefac_val, enh_val
real(wp), intent(in) :: sm_coeff_1, sm_coeff_2, sm_coeff_3

integer :: n
integer :: max_iters
real(dp)     :: visc_val, res
logical      :: flag_rescheck1, flag_rescheck2

real(wp), parameter :: eps = 1.0e-05_wp   ! convergence parameter

!-------- Determination of the order of magnitude --------

visc_val = 1.0e+10_wp   ! initial guess (very low value)

flag_rescheck1 = .false.
n              = 0
max_iters      = 30

do while ((.not.flag_rescheck1).and.(n <= max_iters))

   n = n+1

   res = fct_visc_sm(de_val_m, ratefac_val, enh_val, visc_val, &
                     sm_coeff_1, sm_coeff_2, sm_coeff_3)

   if (res < 0.0_wp) then
      visc_val = 10.0_wp*visc_val
   else
      flag_rescheck1 = .true.
   end if

end do

!-------- Newton's method --------

if (flag_rescheck1) then
   ! only if order of magnitude could be detected successfully

   flag_rescheck2 = .false.
   n              = 0
   max_iters      = 1000

   do while ((.not.flag_rescheck2).and.(n <= max_iters))

      n = n+1

      visc_val = visc_val &
                 - res &
                   /fct_visc_sm_deriv(de_val_m, ratefac_val, &
                                      enh_val, visc_val, &
                                      sm_coeff_1, sm_coeff_2, sm_coeff_3)

      res = fct_visc_sm(de_val_m, ratefac_val, enh_val, visc_val, &
                        sm_coeff_1, sm_coeff_2, sm_coeff_3)

      if (abs(res) < eps) then
         flag_rescheck2 = .true. 
      end if

   end do

end if

visc_iter_sm = visc_val
         
end function visc_iter_sm

!-------------------------------------------------------------------------------
!> Viscosity polynomial
!! [equation (4.33) by Greve and Blatter (Springer, 2009)].
!<------------------------------------------------------------------------------
function fct_visc_sm(de_val_m, ratefac_val, enh_val, visc_val, &
                     sm_coeff_1, sm_coeff_2, sm_coeff_3)

use sico_types_m

implicit none

real(wp)             :: fct_visc_sm
real(wp), intent(in) :: de_val_m
real(wp), intent(in) :: ratefac_val, enh_val
real(wp), intent(in) :: visc_val
real(wp), intent(in) :: sm_coeff_1, sm_coeff_2, sm_coeff_3

real(wp) :: de_visc_factor

de_visc_factor = de_val_m*de_val_m*visc_val*visc_val
                   
fct_visc_sm = 2.0_wp*enh_val*ratefac_val*visc_val &
              * ( sm_coeff_1 &
                  + 4.0_wp*sm_coeff_2*de_visc_factor &
                    * ( 1.0_wp + 4.0_wp*sm_coeff_3*de_visc_factor ) ) &
              - 1.0_wp

end function fct_visc_sm

!-------------------------------------------------------------------------------
!> Derivative of the viscosity polynomial
!! [equation (4.33) by Greve and Blatter (Springer, 2009)].
!<------------------------------------------------------------------------------
function fct_visc_sm_deriv(de_val_m, ratefac_val, enh_val, visc_val, &
                           sm_coeff_1, sm_coeff_2, sm_coeff_3)

use sico_types_m

implicit none

real(wp)             :: fct_visc_sm_deriv
real(wp), intent(in) :: de_val_m
real(wp), intent(in) :: ratefac_val, enh_val
real(wp), intent(in) :: visc_val
real(wp), intent(in) :: sm_coeff_1, sm_coeff_2, sm_coeff_3

real(wp) :: de_visc_factor

real(wp), parameter :: twenty_over_three = 6.666666666666667_wp
                          
de_visc_factor = de_val_m*de_val_m*visc_val*visc_val
                          
fct_visc_sm_deriv = 2.0_wp*sm_coeff_1*enh_val*ratefac_val &
                   + 24.0_wp*sm_coeff_2*enh_val*ratefac_val*de_visc_factor &
                     * ( 1.0_wp + twenty_over_three*sm_coeff_3*de_visc_factor )
         
end function fct_visc_sm_deriv


!-------------------------------------------------------------------------------

end module ice_material_properties_m
!
