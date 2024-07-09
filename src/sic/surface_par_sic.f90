!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s u r f a c e _ p a r _ s i c
!
!  Purpose : sea ice surface parameters
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
module surface_par_sic

  use precision, only : wp
  use constants, only : T0, karman, g, z_sfl, frac_vu
  use sic_params, only : i_ice_albedo, alb_sic_vis, alb_sic_nir, c_alb_sic, snow_par
  use sic_params, only : i_cd_ocn, l_neutral_ocn, l_neutral_sic, z0m_ocn, z0m_sic, Cde0, Cdh0
  use sic_params, only : f_Ri_stab, f_Ri_unstab

  real(wp), parameter :: zm_to_zh_sic = exp(-2._wp)  
  real(wp), parameter :: zm_to_zh_ocn = exp(-2._wp)  

  private
  public :: snow_albedo, ocn_sic_albedo, cdrag

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s n o w _ a l b e d o
  !   Purpose    :  albedo of snow 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine snow_albedo(f_sic, w_snow, w_snow_max, &
                        t_skin_sic, snow, dust_dep, coszm, &
                        alb_snow_vis_dir, alb_snow_vis_dif, alb_snow_nir_dir, alb_snow_nir_dif, &
                        snow_grain, dust_con)

    implicit none

    real(wp), intent(in   ) :: f_sic
    real(wp), intent(in   ) :: w_snow, w_snow_max
    real(wp), intent(in   ) :: t_skin_sic
    real(wp), intent(in   ) :: snow, dust_dep, coszm
    real(wp), intent(inout) :: alb_snow_vis_dir, alb_snow_vis_dif, alb_snow_nir_dir, alb_snow_nir_dif
    real(wp), intent(inout) :: snow_grain
    real(wp), intent(inout) :: dust_con


    if (f_sic.gt.0._wp) then

      ! snow grain size
      if (snow_par%l_snow_aging) then
        call snow_grain_size(t_skin_sic, snow, &
                             snow_grain)
      else
        snow_grain = snow_par%snow_grain_fresh
      endif

      ! dust effect on snow albedo
      if (snow_par%l_snow_dust) then
        ! dust concentration in top snow layer
        call dust_in_snow(dust_dep, snow, w_snow, w_snow_max, &
                          dust_con)
      else
        dust_con = 0._wp
      endif

      ! compute snow albedo
      if (snow_par%i_snow_albedo.eq.1) then
        ! climber-2 snow albedo parameterisation, following Warren & Wiscombe 1980
        call snow_albedo_ww(snow_grain, dust_con, coszm, &
                            alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif)
      else if (snow_par%i_snow_albedo.eq.2) then
        ! snow albedo parameterisation following Dang et al 2015
        call snow_albedo_dang(snow_grain, dust_con, coszm, &
                              alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif)
      endif

    else

      alb_snow_vis_dif = 0._wp 
      alb_snow_nir_dif = 0._wp 
      alb_snow_vis_dir = 0._wp 
      alb_snow_nir_dir = 0._wp 

    endif


    return

    end subroutine snow_albedo


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s n o w _ a l b e d o _ w w
  !   Purpose    :  compute snow albedo following Warren & Wiscombe 1980
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine snow_albedo_ww(snow_grain, dust_con, coszm, &
                            alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif)

    implicit none

    real(wp), intent(in ) :: snow_grain, dust_con, coszm
    real(wp), intent(out) :: alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif

    integer :: k
    real(wp) :: f_cosz, f_age
    real(wp) :: d1, rint, c_dust_new_vis, c_dust_age_vis, c_dust_new_nir, c_dust_age_nir
    real(wp) :: c_age_vis, c_age_nir
    real(wp), parameter :: d_alb_age_vis = 0.05_wp
    real(wp), parameter :: d_alb_age_nir = 0.25_wp

    real(wp), dimension(4) :: tab0 = (/1.001_wp,10._wp,100._wp,1000._wp/)
    real(wp), dimension(4) :: tab1 = (/0.00_wp, 0.02_wp,0.10_wp,0.30_wp/)
    real(wp), dimension(4) :: tab2 = (/0.01_wp, 0.05_wp,0.15_wp,0.30_wp/)


    ! Clear sky snow albedo, zenit angle dependence
    ! zenith angle factor slightly modified from BATS
    ! 2 in the denominator instead of 4 and applied also to angles < 60
    ! => better agreement with Gardner & Sharp 2010
    f_cosz = 0.5_wp*(3._wp/(1._wp+2._wp*coszm)-1._wp)
    f_cosz = max(0._wp,f_cosz)

    ! compute effect of dust concentration on snow albedo after Warren and Wiscombe 1980, Figure 5
    d1 = min(dust_con*1.d6,999._wp)
    if (d1.gt.1.0001_wp) then
      if (d1.gt.1000._wp) d1=1000._wp
      if (d1.ge.1.0001_wp .and. d1.lt.10._wp) k=1
      if (d1.ge.10._wp .and. d1.lt.100._wp)   k=2
      if (d1.ge.100._wp .and. d1.le.1000._wp) k=3
      rint = (log(d1)-log(tab0(k)))/(log(tab0(k+1))-log(tab0(k)))
      c_dust_new_vis = (1._wp-rint)*tab1(k)+rint*(tab1(k+1))
      c_dust_age_vis = (1._wp-rint)*tab2(k)+rint*(tab2(k+1))
    else
      c_dust_new_vis = 0._wp
      c_dust_age_vis = 0._wp
    endif
    ! NIR reduction of snow albedo by dust largely overestimated by this scheme!
    c_dust_new_nir = 0.5_wp*c_dust_new_vis
    c_dust_age_nir = 0.5_wp*c_dust_age_vis
    c_age_vis = d_alb_age_vis + c_dust_age_vis
    c_age_nir = d_alb_age_nir + c_dust_age_nir

    !f_age = (snow_grain-snow_par%snow_grain_fresh)/(snow_par%snow_grain_old-snow_par%snow_grain_fresh)
    f_age = log10(1._wp+(snow_grain-snow_par%snow_grain_fresh)/200._wp)/log10(1._wp+(snow_par%snow_grain_old-snow_par%snow_grain_fresh)/200._wp)
    alb_snow_vis_dif = snow_par%alb_snow_vis_dif_new - f_age*c_age_vis - c_dust_new_vis
    alb_snow_nir_dif = snow_par%alb_snow_nir_dif_new - f_age*c_age_nir - c_dust_new_nir
    alb_snow_vis_dir = alb_snow_vis_dif + 0.4_wp*f_cosz*(1._wp-alb_snow_vis_dif)
    alb_snow_nir_dir = alb_snow_nir_dif + 0.4_wp*f_cosz*(1._wp-alb_snow_nir_dif)


    return

  end subroutine snow_albedo_ww


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s n o w _ a l b e d o _ d a n g
  !   Purpose    :  compute snow albedo following Dang et al 2015
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine snow_albedo_dang(snow_grain,dust_con,coszm, &
                              alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif)

    implicit none

    real(wp), intent(in) :: snow_grain, dust_con, coszm
    real(wp), intent(out) :: alb_snow_vis_dir, alb_snow_vis_dif, alb_snow_nir_dir, alb_snow_nir_dif

    real(wp) :: r, rn, c
    real(wp) :: x, f, H, p
    real(wp) :: alpha_age_vis_dif, alpha_age_nir_dif, alpha_age_vis_dir, alpha_age_nir_dir
    real(wp) :: dalpha_vis_dif, dalpha_nir_dif, dalpha_vis_dir, dalpha_nir_dir
    real(wp), parameter :: r0 = 100._wp ! um
    real(wp), parameter :: c0 = 1.e-6_wp ! kg/kg
    real(wp), parameter :: dust_min = 1.e-8_wp  ! kg/kg


    ! snow grain radius
    r = snow_grain
    rn = log10(r/r0)
    if (dust_con.gt.dust_min) then
      x = log10(dust_con*1.e6_wp)
    endif

    !----------------------------
    ! diffuse visible snow albedo

    ! diffuse vis albedo including snow aging effect, eq. 2
    alpha_age_vis_dif = 0.9856_wp - 0.0202_wp*rn - 0.0125_wp*rn**2 
    ! effect of black carbon (equivalent)
    if (dust_con.gt.dust_min) then
      ! black carbon equivalent, diffuse vis radiation, eq. 9
      f = 152._wp + 15.92_wp*x - 0.39_wp*x**2
      c = dust_con/f 
      H = c/c0*(r/r0)**0.73_wp
      p = log10(H)
      dalpha_vis_dif = 10._wp**(-0.050_wp*p**2+0.514_wp*p-0.890_wp)
    else
      dalpha_vis_dif = 0._wp
    endif
    ! visible diffuse snow albedo
    alb_snow_vis_dif = alpha_age_vis_dif - dalpha_vis_dif

    !----------------------------
    ! diffuse near-infrared snow albedo

    ! diffuse nir albedo including snow aging effect
    alpha_age_nir_dif = 0.7493_wp - 0.1820_wp*rn - 0.0388_wp*rn**2
    ! effect of dust in the infrared is small
    dalpha_nir_dif = 0._wp
    ! near-infrared diffuse snow albedo
    alb_snow_nir_dif = alpha_age_nir_dif - dalpha_nir_dif

    !----------------------------
    ! clear-sky visible snow albedo

    ! solar zenith angle effect
    r = snow_grain*(1._wp+0.781_wp*(coszm-0.65_wp)**2) ! effective grain size corrected for zenith angle
    rn = log10(r/r0)
    ! direct albedo including snow aging effect
    alpha_age_vis_dir = 0.9849_wp - 0.0215_wp*rn - 0.0132_wp*rn**2
    ! effect of black carbon (equivalent)
    if (dust_con.gt.dust_min) then
      ! black carbon equivalent, direct radiation, equation (9)
      f = 155._wp + 17.15_wp*x + 0.27_wp*x**2
      c = dust_con/f 
      H = c/c0*(r/r0)**0.73_wp
      p = log10(H)
      dalpha_vis_dir = 10._wp**(-0.049_wp*p**2+0.525_wp*p-0.893_wp)
    else
      dalpha_vis_dir = 0._wp
    endif
    ! allband direct snow albedo
    alb_snow_vis_dir = alpha_age_vis_dir - dalpha_vis_dir

    !----------------------------
    ! clear-sky near-infrared snow albedo

    ! solar zenith angle effect
    r = snow_grain*(1._wp+0.791_wp*(coszm-0.65_wp)**2) ! effective grain size corrected for zenith angle
    rn = log10(r/r0)
    ! direct albedo including snow aging effect
    alpha_age_nir_dir = 0.6596_wp - 0.1927_wp*rn - 0.0229_wp*rn**2
    ! effect of dust in the infrared is small
    dalpha_nir_dir = 0._wp
    ! allband direct snow albedo
    alb_snow_nir_dir = alpha_age_nir_dir - dalpha_nir_dir


    return

  end subroutine snow_albedo_dang


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s n o w _ g r a i n _ s i z e
  !   Purpose    :  compute snow grain size
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine snow_grain_size(t_skin,snow, &
                             snow_grain)

    implicit none

    real(wp), intent(in) :: t_skin, snow
    real(wp), intent(inout) :: snow_grain

    real(wp) :: f_age, f_tage1, f_tage2, f_tage, f_p


    ! CLIMBER-2 parameterisation of snow grain size (age), tuned to MARv3.6 (using
    ! CROCUS snow model) for Greenland

    ! dry snow temperature dependence for snow grain size
    f_tage1 = exp( snow_par%f_age_t * min(0._wp, (t_skin-T0)) )
    ! melting snow temperature dependence for snow grain size
    f_tage2 = exp( min(0._wp, t_skin-T0) )
    ! snow grain size temperature factor
    f_tage = f_tage1 + f_tage2
    ! averaged 'snow age'       
    f_p = f_tage * (snow_par%snow_0/max(1.e-20_wp,snow))**snow_par%snow_1
    f_age = 1._wp - log(1._wp+f_p)/f_p
    ! snow grain size
    snow_grain = snow_par%snow_grain_fresh + (snow_par%snow_grain_old-snow_par%snow_grain_fresh)*f_age 


    return

  end subroutine snow_grain_size


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  d u s t _ i n _ s n o w
  !   Purpose    :  compute dust concentration in top snow layer
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine dust_in_snow(dust_dep,snow,w_snow,w_snow_max, &
                          dust_con)

    implicit none

    real(wp), intent(in) :: dust_dep
    real(wp), intent(in) :: snow
    real(wp), intent(in) :: w_snow, w_snow_max
    real(wp), intent(inout) :: dust_con

    real(wp) :: dust_con_melt_fac

    real(wp), parameter :: dust_con_max = 1000._wp*1.e-6_wp ! kg/kg


    ! compute dust concentration in snowfall
    dust_con = dust_dep/max(1.e-7_wp,snow) ! kg(dust)/m2/s * kg(snow)/m2/s = kg(dust)/kg(snow)
    ! increase dust concentration when snow melts, assuming dust is not removed with melted water
    ! scavenging efficiency of dust with meltwater is 10-30% (Doherty 2013)
    if (w_snow.gt.1._wp) then 
      dust_con_melt_fac = 1._wp + (w_snow_max-w_snow)/10._wp  ! melt of 10 kg/m2 swe causes a doubling of the dust concentration
      dust_con_melt_fac = max(dust_con_melt_fac,1._wp)
      dust_con_melt_fac = min(dust_con_melt_fac,5._wp)  ! limit to a five-fold increase
    else
      dust_con_melt_fac = 1._wp
    endif
    dust_con = dust_con_melt_fac*dust_con

    ! avoid underflow and negative values
    if (dust_con.lt.1.e-15_wp) dust_con = 0._wp
    ! limit dust concentration in snow to dust_con_max
    dust_con = min(dust_con,dust_con_max)


    return

  end subroutine dust_in_snow


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  o c n _ s i c _ a l b e d o
  !   Purpose    :  compute ocean and sea ice albedo 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_sic_albedo(f_ocn, f_sic, h_snow, t_skin_sic, coszm, &
                           alb_snow_vis_dir, alb_snow_vis_dif, alb_snow_nir_dir, alb_snow_nir_dif, &
                           alb_ocn_vis_dir, alb_ocn_vis_dif, alb_ocn_nir_dir, alb_ocn_nir_dif, &
                           alb_sic_vis_dir, alb_sic_vis_dif, alb_sic_nir_dir, alb_sic_nir_dif, &
                           albedo_sic, albedo_ocn)

    implicit none

    real(wp), intent(in) :: f_ocn
    real(wp), intent(in) :: f_sic, h_snow, t_skin_sic
    real(wp), intent(in) :: coszm
    real(wp), intent(in) :: alb_snow_vis_dir, alb_snow_vis_dif, alb_snow_nir_dir, alb_snow_nir_dif
    real(wp), intent(out) :: alb_ocn_vis_dir, alb_ocn_vis_dif, alb_ocn_nir_dir, alb_ocn_nir_dif
    real(wp), intent(out) :: alb_sic_vis_dir, alb_sic_vis_dif, alb_sic_nir_dir, alb_sic_nir_dif
    real(wp), intent(out) :: albedo_sic, albedo_ocn

    real(wp) :: f_snow, fitemp
    real(wp), parameter :: h_snow_crit = 0.02_wp ! m, critical snow height for snow fraction computation
    real(wp) :: alb_ice_vis_dir, alb_ice_vis_dif, alb_ice_nir_dir, alb_ice_nir_dif


    !------------------------------------------------------------------
    ! ocean albedo, with zenith angle depenence from CLIMBER-2

    if (f_ocn.gt.0._wp) then

      if (coszm .gt. 1.e-4_wp) then
        alb_ocn_vis_dir = min(0.03_wp/coszm,0.2_wp)
      else
        alb_ocn_vis_dir = 0.03_wp
      endif
      alb_ocn_nir_dir = alb_ocn_vis_dir
      alb_ocn_vis_dif = 0.06_wp
      alb_ocn_nir_dif = 0.06_wp

    else  ! no ocean

      alb_ocn_vis_dir = 0._wp
      alb_ocn_nir_dir = 0._wp
      alb_ocn_vis_dif = 0._wp
      alb_ocn_nir_dif = 0._wp

    endif

    !------------------------------------------------------------------
    ! sea ice albedo, considering snow on top

    if (f_sic.gt.0._wp) then

      if (i_ice_albedo.eq.1) then
        ! sea ice albedo scheme accounting for snow 

        ! bare sea ice albedo, with temperature dependence from CLIMBER-2
        fitemp = 0.2_wp*(t_skin_sic-T0+5._wp)
        fitemp = max(0._wp,fitemp)
        fitemp = min(1._wp,fitemp)

        alb_ice_vis_dir = alb_sic_vis-c_alb_sic*fitemp
        alb_ice_nir_dir = alb_sic_nir-c_alb_sic*fitemp
        alb_ice_vis_dif = alb_ice_vis_dir  
        alb_ice_nir_dif = alb_ice_nir_dir  

      else if (i_ice_albedo.eq.2) then
        ! http://www.cesm.ucar.edu/models/atm-cam/docs/description/node35.html

        ! bare sea ice albedo
        fitemp = c_alb_sic*(t_skin_sic-T0+1._wp)
        fitemp = max(0._wp,fitemp)
        fitemp = min(1._wp,fitemp)

        alb_ice_vis_dir = alb_sic_vis-fitemp
        alb_ice_nir_dir = alb_sic_nir-fitemp
        alb_ice_vis_dif = alb_ice_vis_dir  
        alb_ice_nir_dif = alb_ice_nir_dir  

      endif

      ! weighting of snow and bare sea ice albedo
      f_snow = h_snow/(h_snow+h_snow_crit)   ! as in CCSM3
      alb_sic_vis_dir = f_snow * alb_snow_vis_dir + (1._wp-f_snow) * alb_ice_vis_dir  
      alb_sic_vis_dif = f_snow * alb_snow_vis_dif + (1._wp-f_snow) * alb_ice_vis_dif
      alb_sic_nir_dir = f_snow * alb_snow_nir_dir + (1._wp-f_snow) * alb_ice_nir_dir
      alb_sic_nir_dif = f_snow * alb_snow_nir_dif + (1._wp-f_snow) * alb_ice_nir_dif

    else    ! no sea ice

      alb_sic_vis_dir = 0._wp
      alb_sic_nir_dir = 0._wp
      alb_sic_vis_dif = 0._wp
      alb_sic_nir_dif = 0._wp

    endif

    !------------------------------------------------------------------
    ! composite albedo, assuming half cloud cover!
    ! used only in offline simulations, otherwise for diagnostics only
    if (f_ocn.gt.0._wp) then
      albedo_ocn = frac_vu * 0.5_wp*(alb_ocn_vis_dir + alb_ocn_vis_dif) &
        + (1._wp-frac_vu) * 0.5_wp*(alb_ocn_nir_dir + alb_ocn_nir_dif)
      albedo_sic = frac_vu * 0.5_wp*(alb_sic_vis_dir + alb_sic_vis_dif) &
        + (1._wp-frac_vu) * 0.5_wp*(alb_sic_nir_dir + alb_sic_nir_dif)
    endif


   return

  end subroutine ocn_sic_albedo


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c d r a g
  !   Purpose    :  compute drag coefficients
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cdrag(f_ocn, t_skin_ocn, t_air_ocn, t_skin_sic, t_air_sic, wind, &
                   rough_m_ocn, rough_h_ocn, Cde_ocn, Cdh_ocn, rough_m_sic, rough_h_sic, Cde_sic, Cdh_sic)

    implicit none

    real(wp), intent(in) :: f_ocn
    real(wp), intent(in) :: t_skin_ocn
    real(wp), intent(in) :: t_air_ocn
    real(wp), intent(in) :: t_skin_sic
    real(wp), intent(in) :: t_air_sic
    real(wp), intent(in) :: wind
    real(wp), intent(out) :: rough_m_ocn, rough_h_ocn, rough_m_sic, rough_h_sic
    real(wp), intent(out) :: Cde_ocn, Cdh_ocn, Cde_sic, Cdh_sic

    real(wp) :: Cd
    real(wp) :: richardson
    real(wp) :: log_m 
    real(wp) :: log_h 
    real(wp) :: Cdrag_h_neutral 


    if (f_ocn.gt.0._wp) then

      ! ocean drag coefficients, see Taylor 2000 for further details
      if (i_cd_ocn.eq.0) then

        rough_m_ocn = z0m_ocn
        rough_h_ocn = rough_m_ocn * zm_to_zh_ocn

        ! exchange coefficient for moisture
        Cde_ocn = Cde0

        ! exchange coefficient for sensible heat
        Cdh_ocn = Cdh0

      else if (i_cd_ocn.eq.1) then

        ! Large & Yaeger 2004, neutral 10 m coefficients

        rough_m_ocn = z0m_ocn
        rough_h_ocn = rough_m_ocn * zm_to_zh_ocn

        ! drag coefficient
        Cd = (2.7_wp/wind+0.142_wp+wind/13.09)*1.e-3_wp

        ! coefficient for moisture
        Cde_ocn = 34.6_wp*sqrt(Cd)*1.e-3_wp

        ! coefficient for heat, with stability dependence
        if (t_skin_ocn.gt.t_air_ocn) then
          ! unstable
          Cdh_ocn = 32.7_wp*sqrt(Cd)*1.e-3_wp
        else
          Cdh_ocn = 18.0_wp*sqrt(Cd)*1.e-3_wp
        endif

      else if (i_cd_ocn.eq.2) then

        rough_m_ocn = z0m_ocn
        rough_h_ocn = rough_m_ocn * zm_to_zh_ocn

        if (l_neutral_ocn) then
          ! neutral stratification
          Cde_ocn = Cde0
          Cdh_ocn = Cdh0
        else
          ! bulk richardson number
          richardson = g * z_sfl * (1._wp - t_skin_ocn/t_air_ocn) / wind**2
          ! account for atmospheric stability through a Richardson number dependence following BATS
          if (richardson.lt.0._wp) then ! "unstable" stratification
            Cde_ocn = Cde0 * (1._wp - f_Ri_unstab * richardson)
            Cdh_ocn = Cdh0 * (1._wp - f_Ri_unstab * richardson)
          else ! "stable" stratification
            Cde_ocn = Cde0 / (1._wp + f_Ri_stab * richardson)
            Cdh_ocn = Cdh0 / (1._wp + f_Ri_stab * richardson)
          endif
        endif

      endif

      ! sea ice drag coefficients
      rough_m_sic = z0m_sic
      rough_h_sic = rough_m_sic * zm_to_zh_sic
      log_m = karman/log(z_sfl/rough_m_sic)
      log_h = karman/log(z_sfl/rough_h_sic)
      Cdrag_h_neutral = log_m * log_h    ! neutral drag coefficient for sea ice
      if (l_neutral_sic) then
        ! neutral stratification
        Cde_sic = Cdrag_h_neutral
      else
        ! bulk richardson number
        richardson = g * z_sfl * (1._wp - t_skin_sic/t_air_sic) / wind**2
        ! account for atmospheric stability through a Richardson number dependence following BATS
        if (richardson.lt.0._wp) then ! "unstable" stratification
          Cde_sic = Cdrag_h_neutral * (1._wp - f_Ri_unstab * richardson)
        else ! "stable" stratification
          Cde_sic = Cdrag_h_neutral / (1._wp + f_Ri_stab * richardson)
        endif
      endif

      Cdh_sic = Cde_sic 
    endif


    return

  end subroutine cdrag

end module surface_par_sic
