!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : e n d _ c e l l
!
!  Purpose : end cells
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
module end_cell

  use precision, only : wp
  use constants, only : rho_w, rho_i, T0
  use lnd_grid, only : dz, nl, nlc, nsurf, npft, i_ice, i_lake, i_bare, flag_veg
  use lnd_params, only : veg_par, pft_par, snow_par, surf_par 


  implicit none

  private
  public :: end_cell_veg, end_cell_ice, end_cell_ice_grd, end_cell_shelf, end_cell_lake

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  e n d _ c e l l _ v e g
  !   Purpose    :  End prognostic variables of vegetated grid cell 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine end_cell_veg(c13_c12_atm,c14_c_atm, &
                          pft_frac,alt,gdd5,gdd,phen,phen_acc,lai_bal,lai,sai, &
                          npp_ann, npp13_ann, npp14_ann, &
                          leaf_c,root_c, &
                          stem_c,veg_h,veg_c,veg_c13,veg_c14, &
                          w_can,s_can,t_skin,t_skin_veg,t_soil, &
                          w_snow,h_snow,mask_snow, &
                          theta,theta_w,theta_i,theta_sat,w_w,w_i, &
                          litter_c,litter_c13,litter_c14,fast_c,fast_c13,fast_c14,slow_c,slow_c13,slow_c14, &
                          cato_c,cato_c13,cato_c14,litter_c_peat,litter_c13_peat,litter_c14_peat, &
                          acro_c,acro_c13,acro_c14,f_peat,f_peat_pot,w_table_min,w_table_peat, &
                          soil_resp,soil_resp13,soil_resp14,soil_c_tot,soil_c13_tot,soil_c14_tot, &
                          ch4_emis_wetland,ch4_emis_peat,c13h4_emis_wetland,c13h4_emis_peat, &
                          z0m) 

    real(wp), intent(in) :: c13_c12_atm, c14_c_atm    
    real(wp), dimension(:), intent(inout) :: pft_frac
    real(wp) :: alt, gdd5, t_skin_veg
    real(wp), dimension(:) :: gdd, phen, phen_acc, lai_bal, lai, sai
    real(wp), dimension(:), intent(inout) :: npp_ann, npp13_ann, npp14_ann
    real(wp), dimension(:), intent(inout) :: leaf_c
    real(wp), dimension(:), intent(inout) :: root_c
    real(wp), dimension(:), intent(inout) :: stem_c 
    real(wp), dimension(:), intent(inout) :: veg_h, veg_c, veg_c13, veg_c14
    real(wp), dimension(:) :: w_can, s_can, t_skin
    real(wp) :: w_snow, h_snow
    integer :: mask_snow
    real(wp), dimension(0:) :: t_soil
    real(wp), dimension(:) :: theta, theta_w, theta_i, theta_sat, w_w, w_i
    real(wp), dimension(:) :: litter_c,litter_c13,litter_c14,fast_c,fast_c13,fast_c14,slow_c,slow_c13,slow_c14
    real(wp), dimension(:) :: cato_c, cato_c13, cato_c14
    real(wp) :: litter_c_peat, litter_c13_peat, litter_c14_peat, acro_c, acro_c13, acro_c14
    real(wp) :: f_peat, f_peat_pot, w_table_min, w_table_peat
    real(wp), dimension(:) :: soil_resp, soil_resp13, soil_resp14
    real(wp), dimension(:) :: soil_c_tot, soil_c13_tot, soil_c14_tot
    real(wp) :: ch4_emis_wetland, ch4_emis_peat, c13h4_emis_wetland, c13h4_emis_peat
    real(wp), dimension(:) :: z0m

    integer :: n

     ! end prognostic variables

     alt   = -1._wp

     ! vegetation properties
     gdd5       = 1000._wp ! check
     gdd        = 0._wp
     phen       = 0._wp
     phen_acc   = 0._wp
     lai_bal    = pft_par%lai_min
     lai        = lai_bal
     sai        = lai_bal * veg_par%sai_scale
     npp_ann    = 0._wp
     npp13_ann  = 0._wp
     npp14_ann  = 0._wp
     leaf_c     = lai_bal/pft_par%sla
     root_c     = leaf_c
     stem_c     = pft_par%awl*lai_bal**pft_par%bwl
     veg_h      = pft_par%awh * lai_bal
     veg_c      = leaf_c + root_c + stem_c
     veg_c13    = veg_c * c13_c12_atm 
     veg_c14    = veg_c * c14_c_atm
     pft_frac   = veg_par%seed_fraction

     w_can   = 0._wp
     s_can   = 0._wp

     do n=1,nsurf
      if( flag_veg(n) .eq. 1 ) t_skin(n) = T0 + 1._wp
     enddo
     t_skin_veg = T0 + 1._wp

     w_snow         = 0._wp
     h_snow         = 0._wp
     mask_snow      = 0

     t_soil(0)      = T0
     t_soil(1:nl)   = T0 + 1._wp 

     ! assume soil is saturated after ice/ocean shelf/lake are retreating
     ! if skin temperature below 0°C all frozen, otherwise all liquid
     theta_i    = 0.5*theta_sat
     theta_w    = 0._wp
     w_w        = theta_w * dz(1:nl) * rho_w
     w_i        = theta_i * dz(1:nl) * rho_i
     theta      = w_w/(rho_w*dz(1:nl)) + w_i/(rho_i*dz(1:nl))

     ! end soil carbon
     litter_c   = 0._wp
     litter_c13 = 0._wp
     litter_c14 = 0._wp
     fast_c     = 0._wp
     fast_c13   = 0._wp
     fast_c14   = 0._wp
     slow_c     = 0._wp
     slow_c13   = 0._wp
     slow_c14   = 0._wp

     ! end peatlands
     cato_c     = 0._wp
     cato_c13   = 0._wp
     cato_c14   = 0._wp
     litter_c_peat  = 0._wp
     acro_c         = 0._wp
     litter_c13_peat= 0._wp
     acro_c13       = 0._wp
     litter_c14_peat= 0._wp
     acro_c14       = 0._wp
     f_peat         = 0._wp
     f_peat_pot     = 0._wp
     w_table_min    = 0._wp
     w_table_peat   = 0._wp

     soil_resp     = 0._wp
     soil_resp13   = 0._wp
     soil_resp14   = 0._wp
     soil_c_tot           = 0._wp
     soil_c13_tot         = 0._wp
     soil_c14_tot         = 0._wp
     ch4_emis_wetland     = 0._wp
     ch4_emis_peat        = 0._wp
     c13h4_emis_wetland   = 0._wp
     c13h4_emis_peat      = 0._wp

     do n=1,npft
      z0m(n)  = 0.1_wp * veg_h(n)
     enddo
     z0m(i_bare) = surf_par%z0m_bare

     return

  end subroutine end_cell_veg 


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  e n d _ c e l l _ i c e
  !   Purpose    :  End prognostic variables for ice
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine end_cell_ice(t_skin,t_ice,w_snow,h_snow,mask_snow)

    real(wp), dimension(:) :: t_skin
    real(wp) :: w_snow, h_snow
    integer :: mask_snow
    real(wp), dimension(0:) :: t_ice


     ! end prognostic variables

     t_skin(i_ice) = T0

     w_snow         = 0._wp
     h_snow         = 0._wp
     mask_snow      = 0

     t_ice          = T0


     return

  end subroutine end_cell_ice


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  e n d _ c e l l _ i c e _ g r d
  !   Purpose    :  End prognostic variables for ice
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine end_cell_ice_grd(litter_c_ice,litter_c13_ice,litter_c14_ice, &
                          fast_c_ice,fast_c13_ice,fast_c14_ice,slow_c_ice,slow_c13_ice,slow_c14_ice, &
                          soil_resp,soil_resp13,soil_resp14,soil_c_tot,soil_c13_tot,soil_c14_tot)

    real(wp), dimension(:) :: litter_c_ice,litter_c13_ice,litter_c14_ice
    real(wp), dimension(:) :: fast_c_ice,fast_c13_ice,fast_c14_ice,slow_c_ice,slow_c13_ice,slow_c14_ice
    real(wp) :: soil_resp, soil_resp13, soil_resp14
    real(wp) :: soil_c_tot, soil_c13_tot, soil_c14_tot


     litter_c_ice   = 0._wp
     litter_c13_ice = 0._wp
     litter_c14_ice = 0._wp
     fast_c_ice     = 0._wp
     fast_c13_ice   = 0._wp
     fast_c14_ice   = 0._wp
     slow_c_ice     = 0._wp
     slow_c13_ice   = 0._wp
     slow_c14_ice   = 0._wp

     soil_resp     = 0._wp 
     soil_resp13   = 0._wp 
     soil_resp14   = 0._wp
     soil_c_tot           = 0._wp 
     soil_c13_tot         = 0._wp 
     soil_c14_tot         = 0._wp


     return

  end subroutine end_cell_ice_grd


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  e n d _ c e l l _ s h e l f
  !   Purpose    :  Initialize prognostic variables for shelf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine end_cell_shelf(litter_c_shelf,litter_c13_shelf,litter_c14_shelf, &
                            fast_c_shelf,fast_c13_shelf,fast_c14_shelf,slow_c_shelf,slow_c13_shelf,slow_c14_shelf, &
                            soil_resp,soil_resp13,soil_resp14,soil_c_tot,soil_c13_tot,soil_c14_tot, &
                            ch4_emis_shelf,c13h4_emis_shelf)

    real(wp), dimension(:) :: litter_c_shelf,litter_c13_shelf,litter_c14_shelf
    real(wp), dimension(:) :: fast_c_shelf,fast_c13_shelf,fast_c14_shelf,slow_c_shelf,slow_c13_shelf,slow_c14_shelf
    real(wp) :: soil_resp, soil_resp13, soil_resp14
    real(wp) :: soil_c_tot, soil_c13_tot, soil_c14_tot
    real(wp) :: ch4_emis_shelf, c13h4_emis_shelf


     ! end shelf carbon
     litter_c_shelf   = 0._wp
     litter_c13_shelf = 0._wp
     litter_c14_shelf = 0._wp
     fast_c_shelf     = 0._wp
     fast_c13_shelf   = 0._wp
     fast_c14_shelf   = 0._wp
     slow_c_shelf     = 0._wp
     slow_c13_shelf   = 0._wp
     slow_c14_shelf   = 0._wp

     soil_resp     = 0._wp 
     soil_resp13   = 0._wp 
     soil_resp14   = 0._wp
     soil_c_tot           = 0._wp 
     soil_c13_tot         = 0._wp 
     soil_c14_tot         = 0._wp
     ch4_emis_shelf       = 0._wp
     c13h4_emis_shelf     = 0._wp


     return

  end subroutine end_cell_shelf


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  e n d _ c e l l _ l a k e
  !   Purpose    :  end prognostic variables for lake
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine end_cell_lake(t_skin,t_lake,w_snow,h_snow,mask_snow, &
                            litter_c_lake,litter_c13_lake,litter_c14_lake, &
                            fast_c_lake,fast_c13_lake,fast_c14_lake,slow_c_lake,slow_c13_lake,slow_c14_lake, &
                            soil_resp,soil_resp13,soil_resp14,soil_c_tot,soil_c13_tot,soil_c14_tot, &
                            ch4_emis_lake,c13h4_emis_lake)

    real(wp), dimension(:) :: t_skin
    real(wp) :: w_snow, h_snow
    integer :: mask_snow
    real(wp), dimension(0:) :: t_lake
    real(wp), dimension(:) :: litter_c_lake,litter_c13_lake,litter_c14_lake
    real(wp), dimension(:) :: fast_c_lake,fast_c13_lake,fast_c14_lake,slow_c_lake,slow_c13_lake,slow_c14_lake
    real(wp) :: soil_resp, soil_resp13, soil_resp14
    real(wp) :: soil_c_tot, soil_c13_tot, soil_c14_tot
    real(wp) :: ch4_emis_lake, c13h4_emis_lake


     ! end prognostic variables

     t_skin(i_lake) = T0

     w_snow         = 0._wp
     h_snow         = 0._wp
     mask_snow      = 0

     t_lake          = T0

     ! end lake carbon
     litter_c_lake   = 0._wp
     litter_c13_lake = 0._wp
     litter_c14_lake = 0._wp
     fast_c_lake     = 0._wp
     fast_c13_lake   = 0._wp
     fast_c14_lake   = 0._wp
     slow_c_lake     = 0._wp
     slow_c13_lake   = 0._wp
     slow_c14_lake   = 0._wp

     soil_resp     = 0._wp 
     soil_resp13   = 0._wp 
     soil_resp14   = 0._wp
     soil_c_tot           = 0._wp 
     soil_c13_tot         = 0._wp 
     soil_c14_tot         = 0._wp
     ch4_emis_lake       = 0._wp
     c13h4_emis_lake     = 0._wp


     return

  end subroutine end_cell_lake
 
end module end_cell
