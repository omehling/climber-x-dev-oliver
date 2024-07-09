!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s m b _ s i m p l e _ m
!
!  Purpose : Simple surface mass balance scheme
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit and Reinhard Calov
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
module smb_simple_m

use precision, only : wp
use coord, only : grid_class


contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s m b _ s i m p l e
  !   Purpose    :  simple surface mass balance scheme
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_simple(grid, co2, Smax65N, mask_maxice, z_sur, dz_dx_sur, dz_dy_sur, &
    ann_smb)

  implicit none

  type(grid_class), intent(in) ::  grid                   ! grid information
  real(wp), intent(in) :: co2                             ! atmospheric CO2 concentration (ppm)
  real(wp), intent(in) :: Smax65N                         ! maximum summer insolation at 65N (W/m2)
  integer, intent(in),  dimension(:,:) :: mask_maxice     ! mask of maximum allowed ice sheet extent (1)
  real(wp), intent(in), dimension(:,:) :: z_sur           ! surface elevation (m)
  real(wp), intent(in), dimension(:,:) :: dz_dx_sur       ! zonal gradient of surface elevation (m/m)
  real(wp), intent(in), dimension(:,:) :: dz_dy_sur       ! meridional gradient of surface elevation (m/m)

  real(wp), intent(out), dimension(:,:) :: ann_smb        ! annual surface mass balance (kg/m2/a)

  integer :: i, j, m1, m2, imi, ipl, jmi, jpl 

  real(wp) rho_i, mean_Smax65N, mod_gradz_aux
  real(wp) mod_gradz(grid%G%nx,grid%G%ny), ice_mask_n(grid%G%nx,grid%G%ny)
  real(wp) a1, a2, a3, a4, a5, b1, b2, b3, b4
  real(wp) zsn, z_sur_aux, dz, smb_aux, mod_gradz_aux1, ice_mask_n_aux

  ! Calculate the Ice-mask considering neighbour points: mask_ice_n

  m1 = 0   ! counter
  do i = 1,grid%G%nx
    do j = 1,grid%G%ny

       if (mask_maxice(i,j).gt.0) then
       m1=m1+1

       imi=i-1
       if (imi.eq.0) imi=grid%G%nx       
       ipl=i+1
       if (ipl.gt.grid%G%nx) ipl=1          
       jmi=j-1
       if (jmi.eq.0) jmi=1       
       jpl=j+1
       if (jpl.gt.grid%G%ny) jpl=grid%G%ny

      ! mod_gradz_aux=(z_sur(ipl,j)-z_sur(imi,j))**2+ & 
      ! (z_sur(i,jpl)-z_sur(i,jmi))**2
      ! mod_gradz(i,j)=1.e-5*sqrt(mod_gradz_aux)
       mod_gradz_aux = (dz_dx_sur(i,j))**2 + (dz_dy_sur(i,j))**2
       mod_gradz(i,j)=sqrt(mod_gradz_aux) 

       ice_mask_n(i,j)= &
          (mask_maxice(imi,jmi)+mask_maxice(i,jmi)+mask_maxice(ipl,jmi)+ &      
           mask_maxice(imi,j)+mask_maxice(i,j)+mask_maxice(ipl,j)+ &
           mask_maxice(imi,jpl)+mask_maxice(i,jpl)+ &
           mask_maxice(ipl,jpl))/9.
        
       endif
       enddo
      enddo 
      


  ! Coefficients I need
  rho_i = 900  ! Ice Density kg/m^3
  mean_Smax65N = 494.5   ! mean Orbital Forcing in the last 1000 yr
  
  ! Coefficients from Optimisation
!  a1 = 9000
!  a1 = 8900
!  a1 = 8500
!  a1 = 8000
!  a1 = 7000
  a1 = 6500
!  a1 = 6000
  a2 = 200
  a3 = 19
  a4 = 20
  a5 = 800000
  b1 = 10.08
  b2 = 1.2e-5
  b3 = 227.6
  b4 = 3.3e-7
  
  ! For the moment I calculate the :atitude like this
  !do j = 1,grid%G%ny
  !    lat(j) = j-1
  !enddo


  m2 = 0   ! counter
  do i = 1,grid%G%nx
    do j = 1,grid%G%ny

       if (mask_maxice(i,j).gt.0) then
            
            m2=m2+1

            zsn = a1 - a2*grid%lat(i,j) + a3*co2 + a4*(Smax65N-mean_Smax65N) ! snowline height
            z_sur_aux = z_sur(i,j)
            dz = z_sur_aux - zsn
            mod_gradz_aux1=mod_gradz(i,j)
            ice_mask_n_aux=ice_mask_n(i,j)
            if (dz.gt.0) then
             smb_aux = b1+b4*mod_gradz_aux1*exp(-(dz/a5)) +b3*(ice_mask_n_aux-1) ! cm/year
            else
             smb_aux = -b2*dz**2 +b3*(ice_mask_n_aux-1) ! cm/year
            endif

            smb_aux = smb_aux * 0.01 * rho_i      ! kg/m2/year
            ann_smb(i,j)=smb_aux
       else
            ann_smb(i,j)=0
       endif


    enddo
  enddo


  return

end subroutine smb_simple

end module smb_simple_m
