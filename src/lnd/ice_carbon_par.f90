!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : i c e _ c a r b o n _ p a r _ m o d
!
!  Purpose : parameters for carbon in soil below ice sheets
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit and Andrey Ganopolski
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
module ice_carbon_par_mod

  use precision, only : wp
  use timer, only : mon
  use lnd_grid, only : z_c, z_int_c, nl, nlc
  use lnd_params, only : soilc_par

  implicit none

  private
  public :: ice_carbon_par

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  b u r i e d _ c a r b o n _ p a r
  !   Purpose    :  update ice carbon decomposition rate and diffusivity
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ice_carbon_par(k_litter_ice,k_fast_ice,k_slow_ice,diff_icec,adv_icec)

    implicit none

    real(wp), dimension(:), intent(inout) :: k_litter_ice, k_fast_ice, k_slow_ice, diff_icec, adv_icec

    integer :: k
    real(wp), dimension(nl) :: klitter_ice, kfast_ice, kslow_ice
    real(wp), dimension(nlc) :: diff, difficec


     ! decomposition rate for carbon ice below ice sheets, 1/s; there is some decomposition in reality (Wadham 2008)
     ! glacial erosion might be important in removing carbon from below ice sheets!!
     ! should there be a transfer between pools? e.g. Zeng 2003
     klitter_ice = soilc_par%k_ice
     kfast_ice   = soilc_par%k_ice
     kslow_ice   = soilc_par%k_ice

     diff = soilc_par%diff_ice  ! m2/s
     if( soilc_par%diff_ice .gt. 0._wp ) then
      ! carbon diffusivity at the soil levels (interfaces), m2/s
      do k=1,nl-1
       difficec(k) = diff(k) * diff(k+1) * (z_c(k+1) - z_c(k))  &
                   / ( diff(k) * (z_c(k+1) - z_int_c(k)) + diff(k+1) * (z_int_c(k) - z_c(k)) )
      enddo
      difficec(nl:nlc) = 0._wp
     else
      difficec = 0._wp
     endif

     ! cumulate
     if (mon.eq.1) then
       k_litter_ice = 0._wp 
       k_fast_ice   = 0._wp 
       k_slow_ice   = 0._wp 
       diff_icec    = 0._wp 
     endif
     k_litter_ice(1:nl) = k_litter_ice(1:nl) + klitter_ice
     k_fast_ice(1:nl)   = k_fast_ice(1:nl) + kfast_ice
     k_slow_ice(1:nl)   = k_slow_ice(1:nl) + kslow_ice
     ! no decomposition in iceial layer
     k_litter_ice(nlc) = 0._wp
     k_fast_ice(nlc)   = 0._wp
     k_slow_ice(nlc)   = 0._wp

     diff_icec    = diff_icec + difficec
     adv_icec     = 0._wp

    return

  end subroutine ice_carbon_par


end module ice_carbon_par_mod
