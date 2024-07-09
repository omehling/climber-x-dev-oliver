module sico_state

  use sico_types_m, only : wp

  implicit none

  type sico_state_class
   real(wp) :: z_mar
   real(wp) :: mean_accum

    !-------- Field quantities --------

   integer, dimension(:,:), allocatable :: maske !! maske(j,i): Ice-land-ocean mask. 0: grounded ice, 1: ice-free land, 2: ocean, 3: floating ice
   integer, dimension(:,:), allocatable :: maske_old !! maske_old(j,i): Old value of maske (at the previous time step)
   integer, dimension(:,:), allocatable :: maske_neu !! maske_neu(j,i): New value of maske computed during an integration step
   integer, dimension(:,:), allocatable :: id_mask !! id_mask(j,i): ice id mask 
   real(wp), dimension(:,:), allocatable :: H_sed !! H_sed(j,i): sediment thickness
   real(wp), dimension(:,:), allocatable :: f_sed !! f_sed(j,i): sediment fraction
   integer, dimension(:,:), allocatable :: mask_ablation_type !! mask_ablation_type(j,i):
   integer, dimension(:,:), allocatable :: n_cts !! n_cts(j,i): Mask for thermal conditions. -1: cold ice base, 0: temperate ice base with cold ice above, 1: temperate ice base with temperate ice layer above (only for POLY)
   integer, dimension(:,:), allocatable :: n_cts_neu !! (.)_neu: New value of quantity (.) computed during an integration step
   integer, dimension(:,:), allocatable :: kc_cts !! kc_cts(j,i): Position kc of the CTS (for COLD, ENTC, ENTM)
   integer, dimension(:,:), allocatable :: kc_cts_neu !! (.)_neu: New value of quantity (.) computed during an integration step

   integer, dimension(:,:), allocatable :: mask_maxextent !! mask_maxextent(j,i): maximum ice extent mask
   integer, dimension(:,:), allocatable :: maske_target    !! maske_target(j,i): target mask
   real(wp), dimension(:,:), allocatable :: H_target  !! H_target(j,i): target ice thickness

   logical,  dimension(:,:), allocatable :: flag_ocean_point
!! flag_grounding_line_1(j,i): Grounding line flag.
!!                              .true.: grounding line point (grounded ice point with at least one floating ice neighbour),
!!                             .false.: otherwise
   logical, dimension(:,:), allocatable :: flag_grounding_line_1
!! flag_grounding_line_2(j,i): Grounding line flag.
!!                              .true.: grounding line point (floating ice point with at least one grounded ice neighbour),
!!                             .false.: otherwise
   logical, dimension(:,:), allocatable :: flag_grounding_line_2
!! flag_calving_front_1(j,i): Calving front flag.
!!                              .true.: calving front point (floating ice point with at least one ocean neighbour),
!!                             .false.: otherwise
   logical, dimension(:,:), allocatable :: flag_calving_front_1
!! flag_calving_front_2(j,i): Calving front flag.
!!                              .true.: calving front point
!!                                      (ocean point with at least one floating ice neighbour),
!!                             .false.: otherwise
   logical, dimension(:,:), allocatable :: flag_calving_front_2
!> flag_grounded_front_a_1(j,i): Land-terminating grounded front flag.
!>                                .true.: grounded front point
!>                                        (grounded ice point with at least
!>                                        one ice-free land neighbour),
!>                               .false.: otherwise
   logical, dimension(:,:), allocatable :: flag_grounded_front_a_1
!> flag_grounded_front_a_2(j,i): Land-terminating grounded front flag.
!>                                .true.: grounded front point
!>                                        (ice-free land point with at least
!>                                        one grounded ice neighbour),
!>                               .false.: otherwise
   logical, dimension(:,:), allocatable :: flag_grounded_front_a_2
!> flag_grounded_front_b_1(j,i): Marine-terminating grounded front flag.
!>                                .true.: grounded front point
!>                                        (grounded ice point with at least
!>                                        one ocean neighbour),
!>                               .false.: otherwise
   logical, dimension(:,:), allocatable :: flag_grounded_front_b_1
!> flag_grounded_front_b_2(j,i): Marine-terminating grounded front flag.
!>                                .true.: grounded front point
!>                                        (ocean point with at least
!>                                        one grounded ice neighbour),
!>                               .false.: otherwise
   logical, dimension(:,:), allocatable :: flag_grounded_front_b_2
!! flag_shelfy_stream_x(j,i): Shelfy stream flag in x-direction, at (i+1/2,j).
!!                             .true.: shelfy stream point
!!                            .false.: otherwise
   logical, dimension(:,:), allocatable :: flag_shelfy_stream_x
!! flag_shelfy_stream_y(j,i): Shelfy stream flag in y-direction, at (i,j+1/2).
!!                             .true.: shelfy stream point
!!                            .false.: otherwise
   logical, dimension(:,:), allocatable :: flag_shelfy_stream_y
!! flag_shelfy_stream(j,i):   Shelfy stream flag on the main grid.
!!                             .true.: grounded ice, and at least one neighbour on the staggered grid is a shelfy stream point
!!                            .false.: otherwise
   logical, dimension(:,:), allocatable :: flag_shelfy_stream
   logical, dimension(:,:), allocatable :: flag_calc_vxy_ssa_x 
   logical, dimension(:,:), allocatable :: flag_calc_vxy_ssa_y

   real(wp), dimension(:,:), allocatable :: z_sl    !! z_sl(j,i): sea level
   real(wp), dimension(:,:), allocatable :: zs !! zs(j,i): Coordinate z of the surface topography
   real(wp), dimension(:,:), allocatable :: zs_std !! zs_std(j,i): standard deviation of the surface topography
   real(wp), dimension(:,:), allocatable :: zl_std !! zl_std(j,i): standard deviation of the bedrock topography
   real(wp), dimension(:,:), allocatable :: zm !! zm(j,i): Coordinate z of the bottom of the upper (kc) ice domain = top of the lower (kt) ice domain (position of the CTS for POLY, equal to zb for ISOT, COLD, ENTC, ENTM)
   real(wp), dimension(:,:), allocatable :: zb !! zb(j,i): Coordinate z of the ice base
   real(wp), dimension(:,:), allocatable :: zl !! zl(j,i): Coordinate z of the lithosphere surface
   real(wp), dimension(:,:), allocatable :: zl_fil !! zl_fil(j,i): Coordinate z of the smoothed lithosphere surface
   real(wp), dimension(:,:), allocatable :: zl0 !! zl0(j,i): zl for isostatically relaxed ice-free conditions
   real(wp), dimension(:,:), allocatable :: H_c !! H_c(j,i): Thickness of ice in the upper (kc) domain (thickness of the cold-ice layer for POLY, entire ice thickness for ISOT, COLD, ENTC, ENTM)
   real(wp), dimension(:,:), allocatable :: H_t !! H_t(j,i): Thickness of ice in the lower (kt) domain (thickness of the temperate layer for POLY, redundant and thus set to zero for ISOT, COLD, ENTC, ENTM)
   real(wp), dimension(:,:), allocatable :: H
   real(wp), dimension(:,:), allocatable :: H_neu
   real(wp), dimension(:,:), allocatable :: H_neu_flow
   real(wp), dimension(:,:), allocatable :: H_neu_tmp
   real(wp), dimension(:,:), allocatable :: H_smb
   real(wp), dimension(:,:), allocatable :: H_adv
   real(wp), dimension(:,:), allocatable :: H_eff
   real(wp), dimension(:,:), allocatable :: H_calv
   real(wp), dimension(:,:), allocatable :: H_sea_neu
   real(wp), dimension(:,:), allocatable :: H_balance
   real(wp), dimension(:,:), allocatable :: dzs_dxi !! dzs_dxi(j,i): Derivative of zs by xi (at i+1/2,j), one-sided gradients at the grounding line
   real(wp), dimension(:,:), allocatable :: dzs_dx_aux !! dzs_dx_aux(j,i): Derivative of zs by xi (at i+1/2,j)
   real(wp), dimension(:,:), allocatable :: dzm_dxi !! dzm_dxi(j,i): Derivative of zm by xi (at i+1/2,j)
   real(wp), dimension(:,:), allocatable :: dzb_dxi !! dzb_dxi(j,i): Derivative of zb by xi (at i+1/2,j)
   real(wp), dimension(:,:), allocatable :: dH_c_dxi !! dH_c_dxi(j,i): Derivative of H_c by xi (at i+1/2,j)
   real(wp), dimension(:,:), allocatable :: dH_t_dxi !! dH_t_dxi(j,i): Derivative of H_t by xi (at i+1/2,j)
   real(wp), dimension(:,:), allocatable :: dzs_deta !! dzs_deta(j,i): Derivative of zs by eta (at i,j+1/2)
   real(wp), dimension(:,:), allocatable :: dzs_dy_aux !! dzs_dy_aux(j,i): Derivative of zs by eta (at i+1/2,j), one-sided gradients at the grounding line
   real(wp), dimension(:,:), allocatable :: dzm_deta !! dzm_deta(j,i): Derivative of zm by eta (at i,j+1/2)
   real(wp), dimension(:,:), allocatable :: dzb_deta !! dzb_deta(j,i): Derivative of zb by eta (at i,j+1/2)
   real(wp), dimension(:,:), allocatable :: dH_c_deta !! dH_c_deta(j,i): Derivative of H_c by eta (at i,j+1/2)
   real(wp), dimension(:,:), allocatable :: dH_t_deta !! dH_t_deta(j,i): Derivative of H_t by eta (at i,j+1/2)
   real(wp), dimension(:,:), allocatable :: dzs_dxi_g !! (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(wp), dimension(:,:), allocatable :: dzm_dxi_g !! (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(wp), dimension(:,:), allocatable :: dzb_dxi_g !! (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(wp), dimension(:,:), allocatable :: dH_c_dxi_g !! (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(wp), dimension(:,:), allocatable :: dH_t_dxi_g !! (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(wp), dimension(:,:), allocatable :: dzs_deta_g !! (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(wp), dimension(:,:), allocatable :: dzm_deta_g !! (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(wp), dimension(:,:), allocatable :: dzb_deta_g !! (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(wp), dimension(:,:), allocatable :: dH_c_deta_g !! (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(wp), dimension(:,:), allocatable :: dH_t_deta_g !! (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(wp), dimension(:,:), allocatable :: dzs_dtau !! dzs_dtau(j,i): Derivative of zs by tau (time)
   real(wp), dimension(:,:), allocatable :: dzm_dtau !! dzm_dtau(j,i): Derivative of zm by tau (time)
   real(wp), dimension(:,:), allocatable :: dzb_dtau !! dzb_dtau(j,i): Derivative of zb by tau (time)
   real(wp), dimension(:,:), allocatable :: dzl_dtau !! dzl_dtau(j,i): Derivative of zl by tau (time)
   real(wp), dimension(:,:), allocatable :: dH_c_dtau !! dH_c_dtau(j,i): Derivative of H_c by tau (time)
   real(wp), dimension(:,:), allocatable :: dH_t_dtau !! dH_t_dtau(j,i): Derivative of H_t by tau (time)
   real(wp), dimension(:,:), allocatable :: dH_t_smooth !! dH_t_dtau(j,i): Derivative of H_t by tau (time)
   integer, dimension(:,:), allocatable :: p_weert !! p_weert(j,i): Weertman exponent for the basal shear stress
   integer, dimension(:,:), allocatable :: q_weert !! q_weert(j,i): Weertman exponent for the basal pressure
   real(wp), dimension(:,:), allocatable :: p_weert_inv !! p_weert_inv(j,i): Inverse of p_weert
   real(wp), dimension(:,:), allocatable :: c_slide !! c_slide(j,i): Basal sliding coefficient
   real(wp), dimension(:,:), allocatable :: c_drag !! c_drag(j,i): Auxiliary quantity for the computation of the basal drag
   real(wp), dimension(:,:), allocatable :: c_fric !! c_fric(j,i): Basal friction coefficient
   real(wp), dimension(:,:), allocatable :: vb_t !! vb_t(j,i): threshold basal velocity for regularized Coulomb friction law
   real(wp), dimension(:,:), allocatable :: delta !! delta(j,i): Fraction of overburden pressure at saturation
   real(wp), dimension(:,:), allocatable :: p_b_w !! p_b_w(j,i): Basal water pressure
   real(wp), dimension(:,:), allocatable :: vx_b !! vx_b(j,i): Velocity in x-direction at the ice base, at (i+1/2,j)
   real(wp), dimension(:,:), allocatable :: vy_b !! vy_b(j,i): Velocity in y-direction at the ice base, at (i,j+1/2)
   real(wp), dimension(:,:), allocatable :: vx_m !! vx_m(j,i): Mean (depth-averaged) velocity in x-direction, at (i+1/2,j)
   real(wp), dimension(:,:), allocatable :: vy_m !! vy_m(j,i): Mean (depth-averaged) velocity in y-direction, at (i,j+1/2)
   real(wp), dimension(:,:), allocatable :: vx_m_1 
   real(wp), dimension(:,:), allocatable :: vy_m_1 
   real(wp), dimension(:,:), allocatable :: vx_m_2 
   real(wp), dimension(:,:), allocatable :: vy_m_2 
   real(wp), dimension(:,:), allocatable :: vx_m_sia, vy_m_sia
   real(wp), dimension(:,:), allocatable :: vx_m_ssa, vy_m_ssa
   real(wp), dimension(:,:), allocatable :: vx_m_prev, vy_m_prev
   real(wp), dimension(:,:), allocatable :: d_help_b !! d_help_b(j,i)
   real(wp), dimension(:,:), allocatable :: ratio_sl !! ratio_sl(j,i): Ratio of basal to surface velocity (slip ratio)
   real(wp), dimension(:,:), allocatable :: ratio_sl_x !! ratio_sl_x(j,i): Ratio of basal to surface velocity (slip ratio) in x-direction, at (i+1/2,j)
   real(wp), dimension(:,:), allocatable :: ratio_sl_y !! ratio_sl_y(j,i): Ratio of basal to surface velocity (slip ratio) in y-direction, at (i,j+1/2)
   real(wp), dimension(:,:), allocatable :: vx_b_g !! (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(wp), dimension(:,:), allocatable :: vy_b_g !! (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(wp), dimension(:,:), allocatable :: vz_b !! vz_b(j,i): Velocity in z-direction at the ice base, at (i,j)
   real(wp), dimension(:,:), allocatable :: vz_m !! vz_m(j,i): Velocity in z-direction at the position z=zm (interface between the upper (kc) and the lower (kt) domain), at (i,j)
   real(wp), dimension(:,:), allocatable :: vx_s_g !! vx_s_g(j,i): Velocity in x-direction at the ice surface, at (i,j)
   real(wp), dimension(:,:), allocatable :: vy_s_g !! vy_s_g(j,i): Velocity in x-direction at the ice surface, at (i,j)
   real(wp), dimension(:,:), allocatable :: vz_s !! vz_s(j,i): Velocity in z-direction at the ice surface, at (i,j)
   real(wp), dimension(:,:), allocatable :: vz_sl !! vz_sl(j,i): Velocity in z-direction at sea level, at (i,j)
   real(wp), dimension(:,:), allocatable :: flui_ave_sia !! flui_ave_sia(j,i): Depth-averaged fluidity of the SIA
   real(wp), dimension(:,:), allocatable :: h_diff !! h_diff(j,i): Diffusivity of the SIA evolution equation of the ice surface
   real(wp), dimension(:,:), allocatable :: qx !! qx(j,i): Volume flux in x-direction (depth-integrated vx, at (i+1/2,j))
   real(wp), dimension(:,:), allocatable :: qy !! qy(j,i): Volume flux in y-direction (depth-integrated vy, at (i,j+1/2))
   real(wp), dimension(:,:), allocatable :: q_gl_g !! q_gl_g(j,i): Volume flux across the grounding line, at (i,j)
   real(wp), dimension(:,:), allocatable :: q_geo !! q_geo(j,i): Geothermal heat flux
   real(wp), dimension(:,:), allocatable :: temp_b !! temp_b(j,i): Basal temperature
   real(wp), dimension(:,:), allocatable :: temph_b !! temph_b(j,i): Basal temperature relative to the pressure melting point
   real(wp), dimension(:,:), allocatable :: Q_bm !! Q_bm(j,i): Basal melting rate
   real(wp), dimension(:,:), allocatable :: Q_bm_float !! Q_bm_float(j,i): Basal melting rate below floating ice
   real(wp), dimension(:,:), allocatable :: Q_tld !! Q_tld(j,i): Water drainage rate from the temperate layer
   real(wp), dimension(:,:), allocatable :: Q_b_tot !! Q_b_tot(j,i): Sum of Q_bm and Q_tld
   real(wp), dimension(:,:), allocatable :: Q_b_apl !! Q_b_apl(j,i): 
   real(wp), dimension(:,:), allocatable :: H_w !! H_w(j,i): Thickness of the water column under the ice base
   real(wp), dimension(:,:), allocatable :: q_w !! q_w(j,i): 
   real(wp), dimension(:,:), allocatable :: q_w_x !! q_w_x(j,i): 
   real(wp), dimension(:,:), allocatable :: q_w_y !! q_w_y(j,i):
   real(wp), dimension(:,:), allocatable :: runoff !! runoff(j,i): Runoff rate at the ice surface
   real(wp), dimension(:,:), allocatable :: runoff_apl !! runoff(j,i): Runoff rate at the ice surface
   real(wp), dimension(:,:), allocatable :: as_perp !! as_perp(j,i): Accumulation-ablation function at the ice surface (SMB)
   real(wp), dimension(:,:), allocatable :: as_perp_apl !! as_perp_apl(j,i): Applied accumulation-ablation function (SMB)
   real(wp), dimension(:,:), allocatable :: mb_source
   real(wp), dimension(:,:), allocatable :: mb_source_apl !! mb_source_apl(j,i): Applied mass balance source (SMB, BMB, calving)
   real(wp), dimension(:,:), allocatable :: smb_corr
   real(wp), dimension(:,:), allocatable :: calving !! calving(j,i): Calving rate of grounded ice
   real(wp), dimension(:,:), allocatable :: calving_apl !! calving(j,i): Calving rate of grounded ice
   integer,  dimension(:,:), allocatable :: mask_mar
   real(wp), dimension(:,:), allocatable :: cst_dist
   integer, dimension(:,:), allocatable :: i_cst
   integer, dimension(:,:), allocatable :: j_cst
   real(wp), dimension(:,:), allocatable :: cos_grad_tc
   real(wp), dimension(:,:), allocatable :: dis_perp
   real(wp), dimension(:,:), allocatable :: t_ocn !! t_ocn(j,i): ocean temperature
   real(wp), dimension(:,:), allocatable :: s_ocn !! s_ocn(j,i): ocean salinity
   real(wp), dimension(:,:), allocatable :: accum !! accum(j,i): 
   real(wp), dimension(:,:), allocatable :: accum_apl !! accum(j,i): C
   real(wp), dimension(:,:), allocatable :: temp_s !! temp_s(j,i): Ice surface temperature
   real(wp), dimension(:,:), allocatable :: temp_g !! temp_g(j,i): Ground temperature
   real(wp), dimension(:,:), allocatable :: am_perp !! am_perp(j,i): Ice volume flux across the z=zm interface
   real(wp), dimension(:,:), allocatable :: am_perp_st !! am_perp_st(j,i): Steady-state part of am_perp (without contribution of dzm_dtau)

   real(wp), dimension(:,:), allocatable :: zs_neu !! (.)_neu: New value of quantity (.) computed during an integration step
   real(wp), dimension(:,:), allocatable :: zm_neu !! (.)_neu: New value of quantity (.) computed during an integration step
   real(wp), dimension(:,:), allocatable :: zb_neu !! (.)_neu: New value of quantity (.) computed during an integration step
   real(wp), dimension(:,:), allocatable :: zl_neu !! (.)_neu: New value of quantity (.) computed during an integration step
   real(wp), dimension(:,:), allocatable :: H_c_neu !! (.)_neu: New value of quantity (.) computed during an integration step
   real(wp), dimension(:,:), allocatable :: H_t_neu !! (.)_neu: New value of quantity (.) computed during an integration step
   real(wp), dimension(:,:), allocatable :: zs_aux 
   real(wp), dimension(:,:), allocatable :: zm_aux 
   real(wp), dimension(:,:), allocatable :: zb_aux 

   real(wp), dimension(:,:,:), allocatable :: vx_c !! vx_c(kc,j,i): Velocity in x-direction in the upper (kc) ice domain (at (i+1/2,j,kc))
   real(wp), dimension(:,:,:), allocatable :: vy_c !! vy_c(kc,j,i): Velocity in y-direction in the upper (kc) ice domain (at (i,j+1/2,kc))
   real(wp), dimension(:,:,:), allocatable :: vz_c !! vz_c(kc,j,i): Velocity in z-direction in the upper (kc) ice domain (at (i,j,kc+1/2))
   real(wp), dimension(:,:,:), allocatable :: temp_c !! temp_c(kc,j,i): Temperature in the upper (kc) ice domain
   real(wp), dimension(:,:,:), allocatable :: temp_c_neu !! (.)_neu: New value of quantity (.) computed during an integration step
   real(wp), dimension(:,:,:), allocatable :: temp_c_m !! temp_c_m(kc,j,i): Melting temperature in the upper (kc) ice domain
   real(wp), dimension(:), allocatable :: temp_c_help
   real(wp), dimension(:,:,:), allocatable :: age_c !! age_c(kc,j,i): Age in the upper (kc) ice domain
   real(wp), dimension(:,:,:), allocatable :: age_c_neu !! (.)_neu: New value of quantity (.) computed during an integration step
   real(wp), dimension(:,:,:), allocatable :: sigma_c !! sigma_c(kc,j,i): Effective stress in the upper (kc) ice domain
   real(wp), dimension(:,:,:), allocatable :: txz_c !! txz_c(kc,j,i): shear stress in the upper (kc) ice domain
   real(wp), dimension(:,:,:), allocatable :: tyz_c !! tyz_c(kc,j,i): shear stress in the upper (kc) ice domain
   real(wp), dimension(:,:,:), allocatable :: enh_c !! enh_c(kc,j,i): Flow enhancement factor in the upper (kc) ice domain

   real(wp), dimension(:,:), allocatable :: de_ssa !! de_ssa(j,i): Effective strain rate of the SSA, at (i,j)
   real(wp), dimension(:,:), allocatable :: vis_int_g !! vis_int_g(j,i): Depth-integrated viscosity of the SSA, at (i,j)
   real(wp), dimension(:,:), allocatable :: vis_ave_g !! vis_ave_g(j,i): Depth-averaged viscosity of the SSA, at (i,j)
   real(wp), dimension(:,:), allocatable :: vis_ave_g_smooth !! vis_ave_g_smooth(j,i): Depth-averaged smoothed viscosity of the SSA, at (i,j)
   real(wp), dimension(:,:), allocatable :: vis_int_sgxy !! vis_int_sgxy(j,i)
   real(wp), dimension(:,:), allocatable :: beta_drag
   real(wp), dimension(:,:), allocatable :: vx_g !! vx_g(j,i): Velocity in x-direction of the SSA, at (i,j)
   real(wp), dimension(:,:), allocatable :: vy_g !! vy_g(j,i): Velocity in y-direction of the SSA, at (i,j)

   real(wp), dimension(:,:,:), allocatable :: vx_t !! vx_t(kt,j,i): Velocity in x-direction in the lower (kt) ice domain (at (i+1/2,j,kt))
   real(wp), dimension(:,:,:), allocatable :: vy_t !! vy_t(kt,j,i): Velocity in y-direction in the lower (kt) ice domain (at (i,j+1/2,kt))
   real(wp), dimension(:,:,:), allocatable :: vz_t !! vz_t(kt,j,i): Velocity in z-direction in the lower (kt) ice domain (at (i,j,kt+1/2))
   real(wp), dimension(:,:,:), allocatable :: omega_t !! omega_t(kt,j,i): Water content in the lower (kt) ice domain
   real(wp), dimension(:,:,:), allocatable :: omega_t_neu !! (.)_neu: New value of quantity (.) computed during an integration step
   real(wp), dimension(:,:,:), allocatable :: temp_t_m !! temp_t_m(kt,j,i): Melting temperature in the lower (kt) ice domain
   real(wp), dimension(:,:,:), allocatable :: age_t !! age_t(kt,j,i): Age in the lower (kt) ice domain
   real(wp), dimension(:,:,:), allocatable :: age_t_neu !! (.)_neu: New value of quantity (.) computed during an integration step
   real(wp), dimension(:,:,:), allocatable :: sigma_t !! sigma_t(kt,j,i): Effective stress in the lower (kt) ice domain
   real(wp), dimension(:,:,:), allocatable :: txz_t !! txz_t(kc,j,i): shear stress in the lower (kt) ice domain
   real(wp), dimension(:,:,:), allocatable :: tyz_t !! tyz_t(kc,j,i): shear stress in the lower (kt) ice domain
   real(wp), dimension(:,:,:), allocatable :: enh_t !! enh_t(kt,j,i): Flow enhancement factor in the lower (kt) ice domain

   real(wp), dimension(:,:,:), allocatable :: temp_r !! temp_r(kr,j,i): Temperature in the bedrock
   real(wp), dimension(:,:,:), allocatable :: temp_r_neu !! (.)_neu: New value of quantity (.) computed during an integration step

   real(wp), dimension(:,:,:), allocatable :: enth_c !! enth_c(kc,j,i): Enthalpy in the upper (kc) ice domain
   real(wp), dimension(:,:,:), allocatable :: enth_c_neu !! (.)_neu: New value of quantity (.) computed during an integration step
   real(wp), dimension(:,:,:), allocatable :: omega_c !! omega_c(kc,j,i): Water content in the upper (kc) ice domain
   real(wp), dimension(:,:,:), allocatable :: omega_c_neu !! (.)_neu: New value of quantity (.) computed during an integration step

   real(wp), dimension(:,:,:), allocatable :: enth_t !! enth_t(kt,j,i): Enthalpy in the lower (kt) ice domain
   real(wp), dimension(:,:,:), allocatable :: enth_t_neu !! (.)_neu: New value of quantity (.) computed during an integration step

   real(wp), dimension(:,:,:), allocatable :: de_c !! de_c(kc,j,i): Full effective strain rate in the upper (kc) ice domain
   real(wp), dimension(:,:,:), allocatable :: lambda_shear_c !! lambda_shear_c(kc,j,i): Shear fraction in the upper (kc) ice domain

   real(wp), dimension(:,:,:), allocatable :: de_t !! de_t(kt,j,i): Full effective strain rate in the lower (kt) ice domain
   real(wp), dimension(:,:,:), allocatable :: lambda_shear_t !! lambda_shear_t(kt,j,i): Shear fraction in the lower (kt) ice domain

   logical,  dimension(:,:), allocatable   :: check_point
   real(wp), dimension(:,:), allocatable   :: weigh_ssta_sia_x, weigh_ssta_sia_y
   real(wp), dimension(:,:), allocatable   :: p_b, p_b_red, p_b_red_lim, tau_b

   real(wp), dimension(:,:), allocatable   :: gamma_slide_inv
   logical, dimension(:,:), allocatable   :: sub_melt_flag

   real(wp), dimension(:,:), allocatable :: upH_x_1, upH_x_2, upH_y_1, upH_y_2
   real(wp), dimension(:,:), allocatable :: czs2, czs3
 end type sico_state_class

  private
  public :: sico_state_class

  end module sico_state
