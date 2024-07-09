!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s t e r e o _ p r o j _ m
!
!> @file
!!
!! Computation of the forward or inverse stereographic projection,
!! alternatively for a spherical or an ellipsoidal planet.
!!
!! @section Copyright
!!
!! Copyright 2009-2017 Ralf Greve, Reinhard Calov, Alex Robinson
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
!> Computation of the forward or inverse stereographic projection,
!! alternatively for a spherical or an ellipsoidal planet.
!<------------------------------------------------------------------------------
module stereo_proj_m

  use sico_types_m,     only : wp
  use sico_params, only : pi, eps

  implicit none

  private
  public :: stereo_forw_ellipsoid, stereo_inv_ellipsoid, &
            stereo_forw_sphere, stereo_inv_sphere

contains

!-------------------------------------------------------------------------------
!> Forward stereographic projection for an ellipsoidal planet.
!<------------------------------------------------------------------------------
  subroutine stereo_forw_ellipsoid(lambda_val, phi_val, A, B, &
                                   lambda0, phi0, x_val, y_val)

  implicit none

  real(wp), intent(in)  :: lambda_val, phi_val, A, B, lambda0, phi0
  real(wp), intent(out) :: x_val, y_val

  integer :: l
  integer :: sign_phi0
  real(wp) :: phi_aux, phi0_aux
  real(wp) :: e, mc, t, tc, kp, rho, phi_p
  real(wp) :: sinphi0, sinlambda0, cosphi0, coslambda0

  if (phi0 > eps) then   ! for northern hemisphere
     sign_phi0 =  1
  else if (phi0 < (-eps)) then   ! for southern hemisphere
     sign_phi0 = -1
  else
     stop ' >>> stereo_forw_ellipsoid: phi0 must be different from zero!'
  end if

  phi_aux  = phi_val * sign_phi0
  phi0_aux = phi0    * sign_phi0

  e=sqrt((A**2-B**2)/(A**2))

  sinphi0    = sin(phi0_aux)
  sinlambda0 = sin(lambda0)
  cosphi0    = cos(phi0_aux)
  coslambda0 = cos(lambda0)
  
  mc=cosphi0/sqrt(1.0_wp-e*e*sinphi0*sinphi0)
  t=sqrt(((1.0_wp-sin(phi_aux))/(1.0_wp+sin(phi_aux)))* &
        ((1.0_wp+e*sin(phi_aux))/(1.0_wp-e*sin(phi_aux)))**e)
  tc=sqrt(((1.0_wp-sinphi0)/(1.0_wp+sinphi0))* &
         ((1.0_wp+e*sinphi0)/(1.0_wp-e*sinphi0))**e)
  rho=A*mc*t/tc

  x_val =              rho*sin(lambda_val-lambda0)
  y_val = -sign_phi0 * rho*cos(lambda_val-lambda0)

  end subroutine stereo_forw_ellipsoid

!-------------------------------------------------------------------------------
!> Inverse stereographic projection for an ellipsoidal planet.
!<------------------------------------------------------------------------------
  subroutine stereo_inv_ellipsoid(x_val, y_val, A, B, &
                                  lambda0, phi0, lambda_val, phi_val)

  implicit none

  real(wp), intent(in)  :: x_val, y_val, A, B, lambda0, phi0
  real(wp), intent(out) :: lambda_val, phi_val

  integer :: l
  integer :: sign_phi0
  real(wp) :: phi_aux, phi0_aux
  real(wp) :: e, mc, t, tc, kp, rho, phi_p, residual
  real(wp) :: sinphi0, sinlambda0, cosphi0, coslambda0

  real(wp), parameter :: eps_residual = 1.0e-09_wp

  if (phi0 > eps) then   ! for northern hemisphere
     sign_phi0 =  1
  else if (phi0 < (-eps)) then   ! for southern hemisphere
     sign_phi0 = -1
  else
     stop ' >>> stereo_inv_ellipsoid: phi0 must be different from zero!'
  end if

  phi0_aux = phi0 * sign_phi0

  e=sqrt((A**2-B**2)/(A**2))

  sinphi0    = sin(phi0_aux)
  sinlambda0 = sin(lambda0)
  cosphi0    = cos(phi0_aux)
  coslambda0 = cos(lambda0)
  
  tc=sqrt(((1.0_wp-sinphi0)/(1.0_wp+sinphi0))* &
         ((1.0_wp+e*sinphi0)/(1.0_wp-e*sinphi0))**e)
  mc=cosphi0/sqrt(1.0_wp-e*e*sinphi0*sinphi0)
  rho=sqrt(x_val*x_val+y_val*y_val)
  t=rho*tc/(A*mc)

  if ((x_val /= 0.0_wp).or.(y_val /= 0.0_wp)) then
     lambda_val = lambda0 + sign_phi0*atan2(y_val,x_val) + 0.5_wp*pi
  else
     lambda_val = lambda0 + 0.5_wp*pi
  end if  

  !  fix-point iteration

  phi_p=0.5_wp*pi-2.0_wp*atan(t)
  l=0
  residual=3600.0_wp
  do while(residual >= eps_residual)
     l=l+1
     phi_aux=0.5_wp*pi-2.0_wp*atan(t*((1.0_wp-e*sin(phi_p))/ &
             (1.0_wp+e*sin(phi_p)))**(0.5_wp*e))
     residual=abs(phi_aux-phi_p)
     phi_p=phi_aux
  end do

  phi_val = phi_aux * sign_phi0

  if (lambda_val < 0.0_wp) then
     lambda_val = lambda_val + 2.0_wp*pi
  else if (lambda_val >= (2.0_wp*pi)) then
     lambda_val = lambda_val - 2.0_wp*pi
  end if

  end subroutine stereo_inv_ellipsoid

!-------------------------------------------------------------------------------
!> Forward stereographic projection for a spherical planet.
!<------------------------------------------------------------------------------
  subroutine stereo_forw_sphere(lambda_val, phi_val, R, lambda0, phi0, &
                                x_val, y_val)

  implicit none

  real(wp), intent(in)  :: lambda_val, phi_val, R, lambda0, phi0
  real(wp), intent(out) :: x_val, y_val

  real(wp) :: K

  if (phi0 > eps) then   ! for northern hemisphere

     K = (cos(0.25_wp*pi-0.5_wp*phi0))**2

     x_val =  2.0_wp*R*K*tan(0.25_wp*pi-0.5_wp*phi_val) &
                        *sin(lambda_val-lambda0)
     y_val = -2.0_wp*R*K*tan(0.25_wp*pi-0.5_wp*phi_val) &
                        *cos(lambda_val-lambda0)

  else if (phi0 < (-eps)) then   ! for southern hemisphere

     K = (cos(0.25_wp*pi+0.5_wp*phi0))**2

     x_val =  2.0_wp*R*K*tan(0.25_wp*pi+0.5_wp*phi_val) &
                        *sin(lambda_val-lambda0)
     y_val =  2.0_wp*R*K*tan(0.25_wp*pi+0.5_wp*phi_val) &
                        *cos(lambda_val-lambda0)

  else

     stop ' >>> stereo_forw_sphere: phi0 must be different from zero!'

  end if

  end subroutine stereo_forw_sphere

!-------------------------------------------------------------------------------
!> Inverse stereographic projection for a spherical planet.
!<------------------------------------------------------------------------------
  subroutine stereo_inv_sphere(x_val, y_val, R, lambda0, phi0, &
                               lambda_val, phi_val)

  implicit none

  real(wp), intent(in)  :: x_val, y_val, R, lambda0, phi0
  real(wp), intent(out) :: lambda_val, phi_val

  real(wp) :: K

  if (phi0 > eps) then   ! for northern hemisphere

     K = (cos(0.25_wp*pi-0.5_wp*phi0))**2

     phi_val = 0.5_wp*pi &
               -2.0_wp*atan(sqrt(x_val**2+y_val**2)/(2.0_wp*R*K))
     if ((x_val /= 0.0_wp).or.(y_val /= 0.0_wp)) then
        lambda_val = lambda0 + atan2(y_val,x_val) + 0.5_wp*pi
     else
        lambda_val = lambda0 + 0.5_wp*pi
     end if

     if (lambda_val < 0.0_wp) then
        lambda_val = lambda_val + 2.0_wp*pi
     else if (lambda_val >= (2.0_wp*pi)) then
        lambda_val = lambda_val - 2.0_wp*pi
     end if

  else if (phi0 < (-eps)) then   ! for southern hemisphere

     K = (cos(0.25_wp*pi+0.5_wp*phi0))**2

     phi_val = -0.5_wp*pi &
                +2.0_wp*atan(sqrt(x_val**2+y_val**2)/(2.0_wp*R*K))

     if ((x_val /= 0.0_wp).or.(y_val /= 0.0_wp)) then
        lambda_val = lambda0 - atan2(y_val,x_val) + 0.5_wp*pi
     else
        lambda_val = lambda0 + 0.5_wp*pi
     end if

     if (lambda_val < 0.0_wp) then
        lambda_val = lambda_val + 2.0_wp*pi
     else if (lambda_val >= (2.0_wp*pi)) then
        lambda_val = lambda_val - 2.0_wp*pi
     end if

  else

     stop ' >>> stereo_inv_sphere: phi0 must be different from zero!'

  end if

  end subroutine stereo_inv_sphere

!-------------------------------------------------------------------------------

end module stereo_proj_m
!
