!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s h e l f _ c a r b o n _ p a r _ m o d
!
!  Purpose : shelf carbon parameters
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
module shelf_carbon_par_mod

  use precision, only : wp
  use timer, only : nstep_mon_lnd, mon
  use constants, only : T0, k_boltz
  use lnd_grid, only : z_int, z, dz, nl
  use lnd_grid, only : z_int_c, z_c, dz_c, nlc
  use lnd_params, only : soilc_par, ch4_par

  implicit none

  private
  public :: shelf_carbon_par

contains


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s h e l f _ c a r b o n _ p a r
  !   Purpose    :  update shelf carbon decomposition rate and diffusivity
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine shelf_carbon_par(t_shelf_cum,theta_w_shelf_cum,theta_i_shelf_cum, &
                              k_litter_shelf,k_fast_shelf,k_slow_shelf,diff_shelfc,adv_shelfc,ch4_frac_shelf)

    implicit none

    real(wp), dimension(:), intent(inout) :: t_shelf_cum, theta_w_shelf_cum, theta_i_shelf_cum 
    real(wp), dimension(:), intent(inout) :: k_litter_shelf, k_fast_shelf, k_slow_shelf, diff_shelfc, adv_shelfc
    real(wp), dimension(:), intent(inout) :: ch4_frac_shelf

    integer :: k
    real(wp), dimension(nl) :: t_shelf_mean, theta_w_shelf_mean, theta_i_shelf_mean, ftemp, fmoist, fdepth, ftemp_ch4
    real(wp), dimension(nl) :: klitter_shelf, kfast_shelf, kslow_shelf
    real(wp), dimension(nlc) :: diff, diffshelfc


    t_shelf_mean       = t_shelf_cum / nstep_mon_lnd
    theta_w_shelf_mean = theta_w_shelf_cum / nstep_mon_lnd
    theta_i_shelf_mean = theta_i_shelf_cum / nstep_mon_lnd

    ! temperature dependence of decomposition rate
    do k=1,nl
      if (soilc_par%iresp_temp.eq.1) then ! Lloyd & Taylor, 1994
        if( t_shelf_mean(k) .gt. 240._wp ) then
          ftemp(k) = exp(308.56_wp * (1._wp/56.02_wp - 1._wp/(46.02_wp+t_shelf_mean(k)-T0)) )
        else
          ftemp(k) = 0._wp
        endif
      elseif (soilc_par%iresp_temp.eq.2) then ! Arrhenius
        if( t_shelf_mean(k) .gt. 260._wp ) then
          ftemp(k) = exp(soilc_par%Ea*(1._wp/(k_boltz*283.15_wp)-1._wp/(k_boltz*t_shelf_mean(k))))
        else
          ftemp(k) = 0._wp
        endif
      elseif (soilc_par%iresp_temp.eq.3) then ! Q10
        ftemp(k) = soilc_par%q10_c**((t_shelf_mean(k)-283.15_wp)/10._wp)
      endif
    enddo

    ! temperature dependence of methane emission fraction of soil respiration
    do k=1,nl
      ftemp_ch4(k) = exp(ch4_par%Ea_ch4*(1._wp/(k_boltz*303.15_wp)-1._wp/(k_boltz*t_shelf_mean(k))))
    enddo
    ch4_frac_shelf = ch4_par%ch4_frac_shelf * ftemp_ch4

    ! soil moisture dependence of decomposition rate
    do k=1,nl
      ! Porporato 2003
      if( theta_w_shelf_mean(k) .le. 0.3_wp ) then  ! 0.3 = theta at field capacity
        fmoist(k) = max(0._wp,theta_w_shelf_mean(k)/0.3_wp) ! linear increase below field capacity
      else
        fmoist(k) = 0.3_wp/theta_w_shelf_mean(k) ! hyperbolic decrease above field capacity
      endif
    enddo

    ! soil depth dependence of decomposition rate, Koven 2013
    do k=1,nl
      fdepth(k) = exp(-z(k)/soilc_par%z_tau)
    enddo

    ! decomposition rate, 1/s
    klitter_shelf = soilc_par%k10_litter * ftemp*fmoist*fdepth 
    kfast_shelf   = soilc_par%k10_fast   * ftemp*fmoist*fdepth 
    kslow_shelf   = soilc_par%k10_slow   * ftemp*fmoist*fdepth 

    diff = soilc_par%diff_shelf  ! m2/s
    ! carbon diffusivity at the soil levels (interfaces), m2/s
    if( soilc_par%diff_shelf .gt. 0._wp ) then
      do k=1,nl-1
        diffshelfc(k) = diff(k) * diff(k+1) * (z_c(k+1) - z_c(k))  &
          / ( diff(k) * (z_c(k+1) - z_int_c(k)) + diff(k+1) * (z_int_c(k) - z_c(k)) )
      enddo
      diffshelfc(nl:nlc) = 0._wp
    else
      diffshelfc = 0._wp
    endif

    ! cumulate
    if (mon.eq.1) then
      k_litter_shelf = 0._wp 
      k_fast_shelf   = 0._wp 
      k_slow_shelf   = 0._wp 
      diff_shelfc    = 0._wp 
    endif

    k_litter_shelf(1:nl) = k_litter_shelf(1:nl) + klitter_shelf
    k_fast_shelf(1:nl)   = k_fast_shelf(1:nl) + kfast_shelf
    k_slow_shelf(1:nl)   = k_slow_shelf(1:nl) + kslow_shelf
    ! no decomposition in burial layer
    k_litter_shelf(nlc) = 0._wp
    k_fast_shelf(nlc)   = 0._wp
    k_slow_shelf(nlc)   = 0._wp

    diff_shelfc    = diff_shelfc + diffshelfc ! same value as for soil carbon in vegetated part
    adv_shelfc     = 0._wp

    ! reset cumulated variables
    t_shelf_cum  = 0._wp
    theta_w_shelf_cum = 0._wp
    theta_i_shelf_cum = 0._wp


    return

  end subroutine shelf_carbon_par


end module shelf_carbon_par_mod
