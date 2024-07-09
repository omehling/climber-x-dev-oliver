module sico_params

  use sico_types_m

  use nml
  use constants, only : g, rho_w, rho_sw, rho => rho_i, L => Lf
  use control, only : out_dir
  use timer, only : sec_year

  implicit none

  type sico_par_class
    logical :: l_restart
    integer :: i_tibet
    integer :: i_temp_init
    real(wp) :: temp_init
    integer :: grid
    real(wp) :: deform
    integer :: l_smooth_zl0
    real(wp) :: sigma_filter_zl0
    real(wp) :: s_filter_zl0
    real(wp) :: zl0_cont
    integer :: l_smooth_zl
    real(wp) :: sigma_filter_zl
    real(wp) :: s_filter_zl
    real(wp) :: zl_cont
    real(wp) :: dtime
    real(wp) :: dtime_temp
    integer :: topograd
    integer :: gl_surf_grad
    integer :: calcmod
    integer :: cts_melting_freezing
    real(wp) :: temp_const
    real(wp) :: age_const
    integer :: basal_hydrology
    real(wp) :: H_w_max
    integer :: dynamics
    integer :: hyb_mode
    integer :: ice_control
    character(len=256) :: lis_opts
    integer :: margin
    integer :: marine_ice_formation
    integer :: marine_ice_calving
    real(wp) :: z_mar
    real(wp) :: fact_z_mar
    real(wp) :: calv_uw_coeff
    real(wp) :: r1_calv_uw
    real(wp) :: r2_calv_uw
    integer :: ice_shelf_calving
    real(wp) :: h_calv_deep
    real(wp) :: h_calv_shallow
    real(wp) :: zl_deep
    real(wp) :: zl_shallow
    real(wp) :: tau_calv
    integer :: i_disc
    integer :: ii_disc
    real(wp) :: c_dis_bm
    real(wp) :: c_dis_clv
    real(wp) :: alpha_tc
    real(wp) :: zl_std_crit_disc
    real(wp) :: m_H
    real(wp) :: m_D
    real(wp) :: m_R
    real(wp) :: m_T
    real(wp) :: tau_dmb
    real(wp) :: r_mar_eff
    integer :: i_dist_coast
    real(wp) :: dist_coast_min
    real(wp) :: dtime_mar_coa
    integer :: fin_visc
    integer :: flow_law
    integer :: enhmod 
    real(wp) :: enh_fact 
    real(wp) :: enh_intg 
    real(wp) :: age_trans_0 
    real(wp) :: date_trans1_0
    real(wp) :: date_trans2_0
    real(wp) :: date_trans3_0
    real(wp) :: enh_compr 
    real(wp) :: enh_shear
    real(wp) :: enh_shelf
    integer :: i_advance
    integer :: thk_evol
    integer :: calcthk
    real(wp) :: ovi_weight 
    real(wp) :: omega_sor
    real(wp) :: mean_accum
    real(wp) :: time_target_topo_init 
    real(wp) :: time_target_topo_final 
    real(wp) :: target_topo_tau 
    character (len=256) :: target_topo_file
    character (len=256) :: mask_maxextent_file
    integer :: mb_account
    logical :: l_alpine
    logical :: l_exclude_grl
    real(wp) :: h_isol_max

    real(wp) :: z_abyss
    integer :: l_abyss
    integer :: marine_ice_basal_melting
    real(wp) :: qbm_marine
    real(wp) :: qbm_min
    real(wp) :: qbm_max

    integer  :: slide_law
    real(wp) :: c_slide
    real(wp) :: gamma_slide
    real(wp) :: p_weert
    real(wp) :: q_weert
    real(wp) :: c_slide_sedi
    integer :: p_weert_sedi
    integer :: q_weert_sedi    
    real(wp) :: q_coulomb
    real(wp) :: c_fric_rock
    real(wp) :: c_fric_sed
    real(wp) :: vb_t_rock
    real(wp) :: vb_t_sed
    real(wp) :: delta_rock
    real(wp) :: delta_sed
    integer :: i_slide_rock_sed
    real(wp) :: h_sed_thresh
    real(wp) :: h_sed_min
    real(wp) :: h_sed_max
    real(wp) :: enh_slide_marine_sed
    integer :: nyears_ramp_up_slide
    integer :: i_p_b_red
    real(wp) :: red_pres_limit_fact
    integer :: hydro_slide_sat_fct
    real(wp) :: c_hw_slide
    real(wp) :: hw0_slide
    real(wp) :: vh_max
    real(wp) :: visc_min
    real(wp) :: visc_max
    real(wp) :: tol_iter_ssa 
    integer  :: n_iter_ssa
    integer  :: iter_init_ssa
    real(wp) :: visc_init_ssa
    integer  :: n_visc_smooth
    real(wp) :: visc_smooth_diff
    real(wp) :: relax_fact_ssa
    real(wp) :: ratio_sl_thresh
    integer  :: ssta_sia_weigh_fct
    real(wp) :: hyb_ref_speed

    real(wp) :: hd_min
    real(wp) :: hd_max

    real(wp) :: numdiff_h_t
    real(wp) :: tau_cts

    real(wp) :: agediff
    real(wp) :: age_min
    real(wp) :: age_max

    integer :: q_litho
    integer :: adv_hor
    integer :: adv_vert

    logical :: l_write_timer
    
    logical :: l_output_3d
  end type sico_par_class

   real(wp),parameter :: NUE = 1.0d-06 !       Water diffusivity = 1.0e-06 kg/(m*s) 
   real(wp),parameter :: BETA = 8.7d-04 !       Clausius-Clapeyron gradient = 8.7e-04 K/m 
   real(wp),parameter :: DELTA_TM_SW = 1.85d+00 !       Melting point depression of sea water due to its average salinity = 1.85 deg C
   real(wp),parameter :: OMEGA_MAX = 1.0d-02 !       Threshold value for the water content of temperate ice = 1%
   real(wp),parameter :: RHO_C_R = 2.0d+06 !       Density times specific heat of the lithosphere = 2.0e+06 J/(m3*K) 
   real(wp),parameter :: KAPPA_R = 3.0d+00 !       Heat conductivity of the lithosphere = 3 W/(m*K) 
   real(wp),parameter :: R_T = 1.8125d+02 !       Coefficient of the water-content dependence in the rate factor for temperate ice = 181.25
   real(wp),parameter :: rhow_rho_ratio = rho_w/RHO
!! RHO_C: Density of crustal material (dust)
   real(wp)                     :: RHO_C
                                                ! only for the Martian ice caps
!! KAPPA_C: Heat conductivity of crustal material (dust)
   real(wp)                     :: KAPPA_C
                                                ! only for the Martian ice caps
!! C_C: Specific heat of crustal material (dust)
   real(wp)                     :: C_C
                                                ! only for the Martian ice caps


!-------- Mathematical constants -------- 

!! pi: Constant pi
   real(wp), parameter :: pi = 3.141592653589793_wp
!! pi_inv: Inverse of pi
   real(wp), parameter :: pi_inv = 1.0_wp/pi
!! pi_180: pi divided by 180 (-> deg to rad)
   real(wp), parameter :: pi_180 = pi/180.0_wp
!! pi_180_inv: 180 divided by pi (-> rad to deg)
   real(wp), parameter :: pi_180_inv = 180.0_wp/pi

!! euler: Euler number
   real(wp), parameter :: euler = 2.718281828459045_wp

!! eps: Small number
   real(wp), parameter :: eps = 1.0e-05_wp
!! epsi: Very small number
   real(wp), parameter :: epsi = 1.0e-12_wp

!! eps_sp: Small number to single-precision accuracy
   real(sp), parameter :: eps_sp = epsilon(1.0_sp)
!! eps_sp_dp: Small number to single-precision accuracy in double precision
   real(dp), parameter :: eps_sp_dp = real(eps_sp,dp)
!! eps_dp: Small number to double-precision accuracy
   real(dp), parameter :: eps_dp = epsilon(1.0_dp)
!! eps_dp: Small number to porking-precision accuracy
   real(dp), parameter :: eps_wp = epsilon(1.0_wp)
!! ice-thickness epsilon, If computations are in double precision, single precision is a good choice
   real(wp), parameter :: eps_H = eps_sp_dp

!! RF(n): Tabulated values for the rate factor of cold ice
   real(wp), dimension(-190:10) :: RF
!! KAPPA(n): Tabulated values for the heat conductivity of ice
   real(wp), dimension(-190:10) :: KAPPA
!! C(n): Tabulated values for the specific heat of ice
   real(wp), dimension(-190:10) :: C

!! kei(n): Tabulated values of the kei function (Kelvin function of zero order)
   real(wp), dimension(-10000:10000) :: kei
!! n_data_kei: Number of tabulated values of the kei function
   integer:: n_data_kei
!! kei_r_max: Maximum value of the argument r of the tabulated kei function
   real(wp) :: kei_r_max
!! kei_r_incr: Increment of the argument r of the tabulated kei function
   real(wp) :: kei_r_incr


contains

  subroutine sico_params_init(par)

    implicit none

    type(sico_par_class), intent(inout) :: par

    ! read model settings and parameters from namelist
    call sico_par_load(trim(out_dir)//"/ice_sico_par.nml",par)


!-------- Table for rate factor A(T)
!         (by Greve, Weis and Hutter, 1998, Paleoclimates 2 (2-3), 133-161)
!         [A] = 1/(s*Pa3), [T] = deg C --------
!
  RF(  10)    = 3.200e-24_wp     ! Dummy
  RF(   9)    = 3.200e-24_wp     ! Dummy
  RF(   8)    = 3.200e-24_wp     ! Dummy
  RF(   7)    = 3.200e-24_wp     ! Dummy
  RF(   6)    = 3.200e-24_wp     ! Dummy
  RF(   5)    = 3.200e-24_wp     ! Dummy
  RF(   4)    = 3.200e-24_wp     ! Dummy
  RF(   3)    = 3.200e-24_wp     ! Dummy
  RF(   2)    = 3.200e-24_wp     ! Dummy
  RF(   1)    = 3.200e-24_wp     ! Dummy
  RF(   0)    = 3.200e-24_wp
  RF(  -1)    = 2.800e-24_wp
  RF(  -2)    = 2.400e-24_wp
  RF(  -3)    = 2.133e-24_wp
  RF(  -4)    = 1.867e-24_wp
  RF(  -5)    = 1.600e-24_wp
  RF(  -6)    = 1.378e-24_wp
  RF(  -7)    = 1.156e-24_wp
  RF(  -8)    = 9.340e-25_wp
  RF(  -9)    = 7.120e-25_wp
  RF( -10)    = 4.900e-25_wp
  RF( -11)    = 4.413e-25_wp
  RF( -12)    = 3.971e-25_wp
  RF( -13)    = 3.571e-25_wp
  RF( -14)    = 3.209e-25_wp
  RF( -15)    = 2.880e-25_wp
  RF( -16)    = 2.584e-25_wp
  RF( -17)    = 2.316e-25_wp
  RF( -18)    = 2.074e-25_wp
  RF( -19)    = 1.855e-25_wp
  RF( -20)    = 1.658e-25_wp
  RF( -21)    = 1.481e-25_wp
  RF( -22)    = 1.321e-25_wp
  RF( -23)    = 1.178e-25_wp
  RF( -24)    = 1.049e-25_wp
  RF( -25)    = 9.337e-26_wp
  RF( -26)    = 8.300e-26_wp
  RF( -27)    = 7.372e-26_wp
  RF( -28)    = 6.541e-26_wp
  RF( -29)    = 5.798e-26_wp
  RF( -30)    = 5.134e-26_wp
  RF( -31)    = 4.542e-26_wp
  RF( -32)    = 4.014e-26_wp
  RF( -33)    = 3.544e-26_wp
  RF( -34)    = 3.125e-26_wp
  RF( -35)    = 2.753e-26_wp
  RF( -36)    = 2.423e-26_wp
  RF( -37)    = 2.130e-26_wp
  RF( -38)    = 1.870e-26_wp
  RF( -39)    = 1.641e-26_wp
  RF( -40)    = 1.438e-26_wp
  RF( -41)    = 1.258e-26_wp
  RF( -42)    = 1.100e-26_wp
  RF( -43)    = 9.603e-27_wp
  RF( -44)    = 8.375e-27_wp
  RF( -45)    = 7.295e-27_wp
  RF( -46)    = 6.347e-27_wp
  RF( -47)    = 5.515e-27_wp
  RF( -48)    = 4.786e-27_wp
  RF( -49)    = 4.148e-27_wp
  RF( -50)    = 3.591e-27_wp
  RF( -51)    = 3.104e-27_wp
  RF( -52)    = 2.680e-27_wp
  RF( -53)    = 2.311e-27_wp
  RF( -54)    = 1.990e-27_wp
  RF( -55)    = 1.711e-27_wp
  RF( -56)    = 1.469e-27_wp
  RF( -57)    = 1.260e-27_wp
  RF( -58)    = 1.079e-27_wp
  RF( -59)    = 9.224e-28_wp
  RF( -60)    = 7.875e-28_wp
  RF( -61)    = 6.714e-28_wp
  RF( -62)    = 5.715e-28_wp
  RF( -63)    = 4.857e-28_wp
  RF( -64)    = 4.121e-28_wp
  RF( -65)    = 3.492e-28_wp
  RF( -66)    = 2.954e-28_wp
  RF( -67)    = 2.494e-28_wp
  RF( -68)    = 2.103e-28_wp
  RF( -69)    = 1.770e-28_wp
  RF( -70)    = 1.488e-28_wp
  RF( -71)    = 1.248e-28_wp
  RF( -72)    = 1.045e-28_wp
  RF( -73)    = 8.734e-29_wp
  RF( -74)    = 7.288e-29_wp
  RF( -75)    = 6.070e-29_wp
  RF( -76)    = 5.046e-29_wp
  RF( -77)    = 4.187e-29_wp
  RF( -78)    = 3.467e-29_wp
  RF( -79)    = 2.866e-29_wp
  RF( -80)    = 2.364e-29_wp
  RF( -81)    = 1.947e-29_wp
  RF( -82)    = 1.599e-29_wp
  RF( -83)    = 1.311e-29_wp
  RF( -84)    = 1.073e-29_wp
  RF( -85)    = 8.760e-30_wp
  RF( -86)    = 7.136e-30_wp
  RF( -87)    = 5.801e-30_wp
  RF( -88)    = 4.705e-30_wp
  RF( -89)    = 3.808e-30_wp
  RF( -90)    = 3.074e-30_wp
  RF( -91)    = 2.476e-30_wp
  RF( -92)    = 1.990e-30_wp
  RF( -93)    = 1.595e-30_wp
  RF( -94)    = 1.275e-30_wp
  RF( -95)    = 1.017e-30_wp
  RF( -96)    = 8.093e-31_wp
  RF( -97)    = 6.422e-31_wp
  RF( -98)    = 5.083e-31_wp
  RF( -99)    = 4.012e-31_wp
  RF(-100)    = 3.158e-31_wp
  RF(-101)    = 2.479e-31_wp
  RF(-102)    = 1.940e-31_wp
  RF(-103)    = 1.514e-31_wp
  RF(-104)    = 1.179e-31_wp
  RF(-105)    = 9.145e-32_wp
  RF(-106)    = 7.074e-32_wp
  RF(-107)    = 5.455e-32_wp
  RF(-108)    = 4.194e-32_wp
  RF(-109)    = 3.213e-32_wp
  RF(-110)    = 2.454e-32_wp
  RF(-111)    = 1.868e-32_wp
  RF(-112)    = 1.417e-32_wp
  RF(-113)    = 1.072e-32_wp
  RF(-114)    = 8.074e-33_wp
  RF(-115)    = 6.062e-33_wp
  RF(-116)    = 4.534e-33_wp
  RF(-117)    = 3.379e-33_wp
  RF(-118)    = 2.508e-33_wp
  RF(-119)    = 1.855e-33_wp
  RF(-120)    = 1.366e-33_wp
  RF(-121)    = 1.002e-33_wp
  RF(-122)    = 7.325e-34_wp
  RF(-123)    = 5.330e-34_wp
  RF(-124)    = 3.861e-34_wp
  RF(-125)    = 2.785e-34_wp
  RF(-126)    = 2.000e-34_wp
  RF(-127)    = 1.430e-34_wp
  RF(-128)    = 1.018e-34_wp
  RF(-129)    = 7.209e-35_wp
  RF(-130)    = 5.081e-35_wp
  RF(-131)    = 3.564e-35_wp
  RF(-132)    = 2.487e-35_wp
  RF(-133)    = 1.727e-35_wp
  RF(-134)    = 1.193e-35_wp
  RF(-135)    = 8.195e-36_wp
  RF(-136)    = 5.599e-36_wp
  RF(-137)    = 3.804e-36_wp
  RF(-138)    = 2.570e-36_wp
  RF(-139)    = 1.726e-36_wp
  RF(-140)    = 1.152e-36_wp
  RF(-141)    = 7.647e-37_wp
  RF(-142)    = 5.043e-37_wp
  RF(-143)    = 3.304e-37_wp
  RF(-144)    = 2.151e-37_wp
  RF(-145)    = 1.391e-37_wp
  RF(-146)    = 8.930e-38_wp
  RF(-147)    = 5.695e-38_wp
  RF(-148)    = 3.605e-38_wp
  RF(-149)    = 2.266e-38_wp
  RF(-150)    = 1.413e-38_wp
  RF(-151)    = 8.747e-39_wp
  RF(-152)    = 5.371e-39_wp
  RF(-153)    = 3.272e-39_wp
  RF(-154)    = 1.976e-39_wp
  RF(-155)    = 1.184e-39_wp
  RF(-156)    = 7.027e-40_wp
  RF(-157)    = 4.135e-40_wp
  RF(-158)    = 2.410e-40_wp
  RF(-159)    = 1.392e-40_wp
  RF(-160)    = 7.961e-41_wp
  RF(-161)    = 4.508e-41_wp
  RF(-162)    = 2.527e-41_wp
  RF(-163)    = 1.401e-41_wp
  RF(-164)    = 7.689e-42_wp
  RF(-165)    = 4.172e-42_wp
  RF(-166)    = 2.238e-42_wp
  RF(-167)    = 1.187e-42_wp
  RF(-168)    = 6.217e-43_wp
  RF(-169)    = 3.216e-43_wp
  RF(-170)    = 1.643e-43_wp
  RF(-171)    = 8.283e-44_wp
  RF(-172)    = 4.120e-44_wp
  RF(-173)    = 2.020e-44_wp
  RF(-174)    = 9.768e-45_wp
  RF(-175)    = 4.653e-45_wp
  RF(-176)    = 2.183e-45_wp
  RF(-177)    = 1.008e-45_wp
  RF(-178)    = 4.581e-46_wp
  RF(-179)    = 2.047e-46_wp
  RF(-180)    = 8.989e-47_wp
  RF(-181)    = 3.878e-47_wp
  RF(-182)    = 1.642e-47_wp
  RF(-183)    = 6.824e-48_wp
  RF(-184)    = 2.780e-48_wp
  RF(-185)    = 1.110e-48_wp
  RF(-186)    = 4.338e-49_wp
  RF(-187)    = 1.659e-49_wp
  RF(-188)    = 6.202e-50_wp
  RF(-189)    = 2.265e-50_wp
  RF(-190)    = 8.076e-51_wp

!-------- Table for heat conductivity kappa(T)
!         (by Greve, Weis and Hutter, 1998, Paleoclimates 2 (2-3), 133-161)
!         [kappa] = W/(m*K), [T] = deg C --------
  KAPPA(  10) = 2.072_wp     ! Dummy
  KAPPA(   9) = 2.072_wp     ! Dummy
  KAPPA(   8) = 2.072_wp     ! Dummy
  KAPPA(   7) = 2.072_wp     ! Dummy
  KAPPA(   6) = 2.072_wp     ! Dummy
  KAPPA(   5) = 2.072_wp     ! Dummy
  KAPPA(   4) = 2.072_wp     ! Dummy
  KAPPA(   3) = 2.072_wp     ! Dummy
  KAPPA(   2) = 2.072_wp     ! Dummy
  KAPPA(   1) = 2.072_wp     ! Dummy
  KAPPA(   0) = 2.072_wp
  KAPPA(  -1) = 2.083_wp
  KAPPA(  -2) = 2.095_wp
  KAPPA(  -3) = 2.107_wp
  KAPPA(  -4) = 2.119_wp
  KAPPA(  -5) = 2.131_wp
  KAPPA(  -6) = 2.144_wp
  KAPPA(  -7) = 2.156_wp
  KAPPA(  -8) = 2.168_wp
  KAPPA(  -9) = 2.181_wp
  KAPPA( -10) = 2.193_wp
  KAPPA( -11) = 2.206_wp
  KAPPA( -12) = 2.218_wp
  KAPPA( -13) = 2.231_wp
  KAPPA( -14) = 2.244_wp
  KAPPA( -15) = 2.256_wp
  KAPPA( -16) = 2.269_wp
  KAPPA( -17) = 2.282_wp
  KAPPA( -18) = 2.295_wp
  KAPPA( -19) = 2.308_wp
  KAPPA( -20) = 2.322_wp
  KAPPA( -21) = 2.335_wp
  KAPPA( -22) = 2.348_wp
  KAPPA( -23) = 2.362_wp
  KAPPA( -24) = 2.375_wp
  KAPPA( -25) = 2.389_wp
  KAPPA( -26) = 2.402_wp
  KAPPA( -27) = 2.416_wp
  KAPPA( -28) = 2.430_wp
  KAPPA( -29) = 2.444_wp
  KAPPA( -30) = 2.458_wp
  KAPPA( -31) = 2.472_wp
  KAPPA( -32) = 2.486_wp
  KAPPA( -33) = 2.500_wp
  KAPPA( -34) = 2.515_wp
  KAPPA( -35) = 2.529_wp
  KAPPA( -36) = 2.543_wp
  KAPPA( -37) = 2.558_wp
  KAPPA( -38) = 2.573_wp
  KAPPA( -39) = 2.587_wp
  KAPPA( -40) = 2.602_wp
  KAPPA( -41) = 2.617_wp
  KAPPA( -42) = 2.632_wp
  KAPPA( -43) = 2.647_wp
  KAPPA( -44) = 2.662_wp
  KAPPA( -45) = 2.677_wp
  KAPPA( -46) = 2.693_wp
  KAPPA( -47) = 2.708_wp
  KAPPA( -48) = 2.723_wp
  KAPPA( -49) = 2.739_wp
  KAPPA( -50) = 2.755_wp
  KAPPA( -51) = 2.770_wp
  KAPPA( -52) = 2.786_wp
  KAPPA( -53) = 2.802_wp
  KAPPA( -54) = 2.818_wp
  KAPPA( -55) = 2.834_wp
  KAPPA( -56) = 2.850_wp
  KAPPA( -57) = 2.867_wp
  KAPPA( -58) = 2.883_wp
  KAPPA( -59) = 2.900_wp
  KAPPA( -60) = 2.916_wp
  KAPPA( -61) = 2.933_wp
  KAPPA( -62) = 2.950_wp
  KAPPA( -63) = 2.966_wp
  KAPPA( -64) = 2.983_wp
  KAPPA( -65) = 3.001_wp
  KAPPA( -66) = 3.018_wp
  KAPPA( -67) = 3.035_wp
  KAPPA( -68) = 3.052_wp
  KAPPA( -69) = 3.070_wp
  KAPPA( -70) = 3.087_wp
  KAPPA( -71) = 3.105_wp
  KAPPA( -72) = 3.123_wp
  KAPPA( -73) = 3.140_wp
  KAPPA( -74) = 3.158_wp
  KAPPA( -75) = 3.177_wp
  KAPPA( -76) = 3.195_wp
  KAPPA( -77) = 3.213_wp
  KAPPA( -78) = 3.231_wp
  KAPPA( -79) = 3.250_wp
  KAPPA( -80) = 3.268_wp
  KAPPA( -81) = 3.287_wp
  KAPPA( -82) = 3.306_wp
  KAPPA( -83) = 3.325_wp
  KAPPA( -84) = 3.344_wp
  KAPPA( -85) = 3.363_wp
  KAPPA( -86) = 3.382_wp
  KAPPA( -87) = 3.401_wp
  KAPPA( -88) = 3.421_wp
  KAPPA( -89) = 3.440_wp
  KAPPA( -90) = 3.460_wp
  KAPPA( -91) = 3.480_wp
  KAPPA( -92) = 3.500_wp
  KAPPA( -93) = 3.520_wp
  KAPPA( -94) = 3.540_wp
  KAPPA( -95) = 3.560_wp
  KAPPA( -96) = 3.580_wp
  KAPPA( -97) = 3.601_wp
  KAPPA( -98) = 3.621_wp
  KAPPA( -99) = 3.642_wp
  KAPPA(-100) = 3.663_wp
  KAPPA(-101) = 3.684_wp
  KAPPA(-102) = 3.705_wp
  KAPPA(-103) = 3.726_wp
  KAPPA(-104) = 3.747_wp
  KAPPA(-105) = 3.769_wp
  KAPPA(-106) = 3.790_wp
  KAPPA(-107) = 3.812_wp
  KAPPA(-108) = 3.834_wp
  KAPPA(-109) = 3.856_wp
  KAPPA(-110) = 3.878_wp
  KAPPA(-111) = 3.900_wp
  KAPPA(-112) = 3.922_wp
  KAPPA(-113) = 3.945_wp
  KAPPA(-114) = 3.967_wp
  KAPPA(-115) = 3.990_wp
  KAPPA(-116) = 4.013_wp
  KAPPA(-117) = 4.036_wp
  KAPPA(-118) = 4.059_wp
  KAPPA(-119) = 4.082_wp
  KAPPA(-120) = 4.105_wp
  KAPPA(-121) = 4.129_wp
  KAPPA(-122) = 4.152_wp
  KAPPA(-123) = 4.176_wp
  KAPPA(-124) = 4.200_wp
  KAPPA(-125) = 4.224_wp
  KAPPA(-126) = 4.248_wp
  KAPPA(-127) = 4.272_wp
  KAPPA(-128) = 4.297_wp
  KAPPA(-129) = 4.321_wp
  KAPPA(-130) = 4.346_wp
  KAPPA(-131) = 4.371_wp
  KAPPA(-132) = 4.396_wp
  KAPPA(-133) = 4.421_wp
  KAPPA(-134) = 4.446_wp
  KAPPA(-135) = 4.472_wp
  KAPPA(-136) = 4.497_wp
  KAPPA(-137) = 4.523_wp
  KAPPA(-138) = 4.549_wp
  KAPPA(-139) = 4.575_wp
  KAPPA(-140) = 4.601_wp
  KAPPA(-141) = 4.627_wp
  KAPPA(-142) = 4.654_wp
  KAPPA(-143) = 4.680_wp
  KAPPA(-144) = 4.707_wp
  KAPPA(-145) = 4.734_wp
  KAPPA(-146) = 4.761_wp
  KAPPA(-147) = 4.788_wp
  KAPPA(-148) = 4.816_wp
  KAPPA(-149) = 4.843_wp
  KAPPA(-150) = 4.871_wp
  KAPPA(-151) = 4.899_wp
  KAPPA(-152) = 4.927_wp
  KAPPA(-153) = 4.955_wp
  KAPPA(-154) = 4.983_wp
  KAPPA(-155) = 5.012_wp
  KAPPA(-156) = 5.040_wp
  KAPPA(-157) = 5.069_wp
  KAPPA(-158) = 5.098_wp
  KAPPA(-159) = 5.127_wp
  KAPPA(-160) = 5.157_wp
  KAPPA(-161) = 5.186_wp
  KAPPA(-162) = 5.216_wp
  KAPPA(-163) = 5.246_wp
  KAPPA(-164) = 5.276_wp
  KAPPA(-165) = 5.306_wp
  KAPPA(-166) = 5.336_wp
  KAPPA(-167) = 5.367_wp
  KAPPA(-168) = 5.397_wp
  KAPPA(-169) = 5.428_wp
  KAPPA(-170) = 5.459_wp
  KAPPA(-171) = 5.490_wp
  KAPPA(-172) = 5.522_wp
  KAPPA(-173) = 5.553_wp
  KAPPA(-174) = 5.585_wp
  KAPPA(-175) = 5.617_wp
  KAPPA(-176) = 5.649_wp
  KAPPA(-177) = 5.681_wp
  KAPPA(-178) = 5.714_wp
  KAPPA(-179) = 5.746_wp
  KAPPA(-180) = 5.779_wp
  KAPPA(-181) = 5.812_wp
  KAPPA(-182) = 5.846_wp
  KAPPA(-183) = 5.879_wp
  KAPPA(-184) = 5.913_wp
  KAPPA(-185) = 5.946_wp
  KAPPA(-186) = 5.980_wp
  KAPPA(-187) = 6.015_wp
  KAPPA(-188) = 6.049_wp
  KAPPA(-189) = 6.084_wp
  KAPPA(-190) = 6.118_wp

!-------- Table for specific heat c(T)
!         (by Greve, Weis and Hutter, 1998, Paleoclimates 2 (2-3), 133-161)
!         [c] = J/(kg*K), [T] = deg C --------
  C(  10)     = 2.127e+03_wp     !  Dummy
  C(   9)     = 2.127e+03_wp     !  Dummy
  C(   8)     = 2.127e+03_wp     !  Dummy
  C(   7)     = 2.127e+03_wp     !  Dummy
  C(   6)     = 2.127e+03_wp     !  Dummy
  C(   5)     = 2.127e+03_wp     !  Dummy
  C(   4)     = 2.127e+03_wp     !  Dummy
  C(   3)     = 2.127e+03_wp     !  Dummy
  C(   2)     = 2.127e+03_wp     !  Dummy
  C(   1)     = 2.127e+03_wp     !  Dummy
  C(   0)     = 2.127e+03_wp
  C(  -1)     = 2.120e+03_wp
  C(  -2)     = 2.113e+03_wp
  C(  -3)     = 2.106e+03_wp
  C(  -4)     = 2.098e+03_wp
  C(  -5)     = 2.091e+03_wp
  C(  -6)     = 2.084e+03_wp
  C(  -7)     = 2.077e+03_wp
  C(  -8)     = 2.069e+03_wp
  C(  -9)     = 2.062e+03_wp
  C( -10)     = 2.055e+03_wp
  C( -11)     = 2.048e+03_wp
  C( -12)     = 2.040e+03_wp
  C( -13)     = 2.033e+03_wp
  C( -14)     = 2.026e+03_wp
  C( -15)     = 2.019e+03_wp
  C( -16)     = 2.011e+03_wp
  C( -17)     = 2.004e+03_wp
  C( -18)     = 1.997e+03_wp
  C( -19)     = 1.990e+03_wp
  C( -20)     = 1.982e+03_wp
  C( -21)     = 1.975e+03_wp
  C( -22)     = 1.968e+03_wp
  C( -23)     = 1.961e+03_wp
  C( -24)     = 1.953e+03_wp
  C( -25)     = 1.946e+03_wp
  C( -26)     = 1.939e+03_wp
  C( -27)     = 1.932e+03_wp
  C( -28)     = 1.924e+03_wp
  C( -29)     = 1.917e+03_wp
  C( -30)     = 1.910e+03_wp
  C( -31)     = 1.903e+03_wp
  C( -32)     = 1.895e+03_wp
  C( -33)     = 1.888e+03_wp
  C( -34)     = 1.881e+03_wp
  C( -35)     = 1.874e+03_wp
  C( -36)     = 1.866e+03_wp
  C( -37)     = 1.859e+03_wp
  C( -38)     = 1.852e+03_wp
  C( -39)     = 1.845e+03_wp
  C( -40)     = 1.837e+03_wp
  C( -41)     = 1.830e+03_wp
  C( -42)     = 1.823e+03_wp
  C( -43)     = 1.816e+03_wp
  C( -44)     = 1.808e+03_wp
  C( -45)     = 1.801e+03_wp
  C( -46)     = 1.794e+03_wp
  C( -47)     = 1.787e+03_wp
  C( -48)     = 1.779e+03_wp
  C( -49)     = 1.772e+03_wp
  C( -50)     = 1.765e+03_wp
  C( -51)     = 1.758e+03_wp
  C( -52)     = 1.750e+03_wp
  C( -53)     = 1.743e+03_wp
  C( -54)     = 1.736e+03_wp
  C( -55)     = 1.729e+03_wp
  C( -56)     = 1.721e+03_wp
  C( -57)     = 1.714e+03_wp
  C( -58)     = 1.707e+03_wp
  C( -59)     = 1.700e+03_wp
  C( -60)     = 1.692e+03_wp
  C( -61)     = 1.685e+03_wp
  C( -62)     = 1.678e+03_wp
  C( -63)     = 1.671e+03_wp
  C( -64)     = 1.663e+03_wp
  C( -65)     = 1.656e+03_wp
  C( -66)     = 1.649e+03_wp
  C( -67)     = 1.642e+03_wp
  C( -68)     = 1.634e+03_wp
  C( -69)     = 1.627e+03_wp
  C( -70)     = 1.620e+03_wp
  C( -71)     = 1.612e+03_wp
  C( -72)     = 1.605e+03_wp
  C( -73)     = 1.598e+03_wp
  C( -74)     = 1.591e+03_wp
  C( -75)     = 1.583e+03_wp
  C( -76)     = 1.576e+03_wp
  C( -77)     = 1.569e+03_wp
  C( -78)     = 1.562e+03_wp
  C( -79)     = 1.554e+03_wp
  C( -80)     = 1.547e+03_wp
  C( -81)     = 1.540e+03_wp
  C( -82)     = 1.533e+03_wp
  C( -83)     = 1.525e+03_wp
  C( -84)     = 1.518e+03_wp
  C( -85)     = 1.511e+03_wp
  C( -86)     = 1.504e+03_wp
  C( -87)     = 1.496e+03_wp
  C( -88)     = 1.489e+03_wp
  C( -89)     = 1.482e+03_wp
  C( -90)     = 1.475e+03_wp
  C( -91)     = 1.467e+03_wp
  C( -92)     = 1.460e+03_wp
  C( -93)     = 1.453e+03_wp
  C( -94)     = 1.446e+03_wp
  C( -95)     = 1.438e+03_wp
  C( -96)     = 1.431e+03_wp
  C( -97)     = 1.424e+03_wp
  C( -98)     = 1.417e+03_wp
  C( -99)     = 1.409e+03_wp
  C(-100)     = 1.402e+03_wp
  C(-101)     = 1.395e+03_wp
  C(-102)     = 1.388e+03_wp
  C(-103)     = 1.380e+03_wp
  C(-104)     = 1.373e+03_wp
  C(-105)     = 1.366e+03_wp
  C(-106)     = 1.359e+03_wp
  C(-107)     = 1.351e+03_wp
  C(-108)     = 1.344e+03_wp
  C(-109)     = 1.337e+03_wp
  C(-110)     = 1.330e+03_wp
  C(-111)     = 1.322e+03_wp
  C(-112)     = 1.315e+03_wp
  C(-113)     = 1.308e+03_wp
  C(-114)     = 1.301e+03_wp
  C(-115)     = 1.293e+03_wp
  C(-116)     = 1.286e+03_wp
  C(-117)     = 1.279e+03_wp
  C(-118)     = 1.272e+03_wp
  C(-119)     = 1.264e+03_wp
  C(-120)     = 1.257e+03_wp
  C(-121)     = 1.250e+03_wp
  C(-122)     = 1.243e+03_wp
  C(-123)     = 1.235e+03_wp
  C(-124)     = 1.228e+03_wp
  C(-125)     = 1.221e+03_wp
  C(-126)     = 1.214e+03_wp
  C(-127)     = 1.206e+03_wp
  C(-128)     = 1.199e+03_wp
  C(-129)     = 1.192e+03_wp
  C(-130)     = 1.185e+03_wp
  C(-131)     = 1.177e+03_wp
  C(-132)     = 1.170e+03_wp
  C(-133)     = 1.163e+03_wp
  C(-134)     = 1.156e+03_wp
  C(-135)     = 1.148e+03_wp
  C(-136)     = 1.141e+03_wp
  C(-137)     = 1.134e+03_wp
  C(-138)     = 1.127e+03_wp
  C(-139)     = 1.119e+03_wp
  C(-140)     = 1.112e+03_wp
  C(-141)     = 1.105e+03_wp
  C(-142)     = 1.098e+03_wp
  C(-143)     = 1.090e+03_wp
  C(-144)     = 1.083e+03_wp
  C(-145)     = 1.076e+03_wp
  C(-146)     = 1.069e+03_wp
  C(-147)     = 1.061e+03_wp
  C(-148)     = 1.054e+03_wp
  C(-149)     = 1.047e+03_wp
  C(-150)     = 1.040e+03_wp
  C(-151)     = 1.032e+03_wp
  C(-152)     = 1.025e+03_wp
  C(-153)     = 1.018e+03_wp
  C(-154)     = 1.010e+03_wp
  C(-155)     = 1.003e+03_wp
  C(-156)     = 9.960e+02_wp
  C(-157)     = 9.887e+02_wp
  C(-158)     = 9.815e+02_wp
  C(-159)     = 9.742e+02_wp
  C(-160)     = 9.670e+02_wp
  C(-161)     = 9.597e+02_wp
  C(-162)     = 9.525e+02_wp
  C(-163)     = 9.452e+02_wp
  C(-164)     = 9.380e+02_wp
  C(-165)     = 9.307e+02_wp
  C(-166)     = 9.235e+02_wp
  C(-167)     = 9.162e+02_wp
  C(-168)     = 9.090e+02_wp
  C(-169)     = 9.017e+02_wp
  C(-170)     = 8.944e+02_wp
  C(-171)     = 8.872e+02_wp
  C(-172)     = 8.799e+02_wp
  C(-173)     = 8.727e+02_wp
  C(-174)     = 8.654e+02_wp
  C(-175)     = 8.582e+02_wp
  C(-176)     = 8.509e+02_wp
  C(-177)     = 8.437e+02_wp
  C(-178)     = 8.364e+02_wp
  C(-179)     = 8.292e+02_wp
  C(-180)     = 8.219e+02_wp
  C(-181)     = 8.147e+02_wp
  C(-182)     = 8.074e+02_wp
  C(-183)     = 8.002e+02_wp
  C(-184)     = 7.929e+02_wp
  C(-185)     = 7.857e+02_wp
  C(-186)     = 7.784e+02_wp
  C(-187)     = 7.711e+02_wp
  C(-188)     = 7.639e+02_wp
  C(-189)     = 7.566e+02_wp
  C(-190)     = 7.494e+02_wp
 
    return

    end subroutine sico_params_init

    subroutine sico_par_load(filename,par)

    implicit none

    character (len=*),   intent(in)    :: filename
    type(sico_par_class), intent(inout) :: par


    ! Read parameters from file
    write(*,*) "sicopolis domain parameters ==========="
    call nml_read(filename,"ice_sico_par","i_tibet",par%i_tibet)
    call nml_read(filename,"ice_sico_par","i_temp_init",par%i_temp_init)
    call nml_read(filename,"ice_sico_par","temp_init",par%temp_init)
    call nml_read(filename,"ice_sico_par","grid",par%grid)
    call nml_read(filename,"ice_sico_par","deform",par%deform)
    call nml_read(filename,"ice_sico_par","l_smooth_zl0",par%l_smooth_zl0)
    call nml_read(filename,"ice_sico_par","sigma_filter_zl0",par%sigma_filter_zl0)
    call nml_read(filename,"ice_sico_par","s_filter_zl0",par%s_filter_zl0)
    call nml_read(filename,"ice_sico_par","zl0_cont",par%zl0_cont)
    call nml_read(filename,"ice_sico_par","l_smooth_zl",par%l_smooth_zl)
    call nml_read(filename,"ice_sico_par","sigma_filter_zl",par%sigma_filter_zl)
    call nml_read(filename,"ice_sico_par","s_filter_zl",par%s_filter_zl)
    call nml_read(filename,"ice_sico_par","zl_cont",par%zl_cont)
    call nml_read(filename,"ice_sico_par","dtime",par%dtime)
    call nml_read(filename,"ice_sico_par","dtime_temp",par%dtime_temp)

    if ((par%l_smooth_zl0==1.or.par%l_smooth_zl0==2).and.par%sigma_filter_zl0<=0) then
      print *,'ERROR, sigma_filter_zl0 should be positive for l_smooth_zl0= ', par%l_smooth_zl0
      print *,'sigma_filter_zl0= ',par%sigma_filter_zl0
      stop
    endif
    if ((par%l_smooth_zl==1.or.par%l_smooth_zl==2).and.par%sigma_filter_zl<=0) then
      print *,'ERROR, sigma_filter_zl should be positive for l_smooth_zl= ', par%l_smooth_zl
      print *,'sigma_filter_zl= ',par%sigma_filter_zl
      stop
    endif
    if (abs(1._wp/par%dtime-nint(1._wp/par%dtime)).gt.1.e-10_wp) then
      print *,'ERROR, Sicopolis time step dtime has to be an integer fraction of 1 year'
      print *,'dtime',par%dtime
      print *,1._wp/par%dtime,nint(1._wp/par%dtime),1._wp/par%dtime-nint(1._wp/par%dtime)
      stop
    endif
    if (abs(par%dtime_temp/par%dtime-nint(par%dtime_temp/par%dtime)).gt.1.e-10_wp) then
      print *,'ERROR, Sicopolis time step dtime_temp has to be an integer multiple of dtime'
      print *,'dtime',par%dtime
      print *,'dtime_temp',par%dtime_temp
      print *,par%dtime_temp/par%dtime,nint(par%dtime_temp/par%dtime),par%dtime_temp/par%dtime-nint(par%dtime_temp/par%dtime)
      stop
    endif
    call nml_read(filename,"ice_sico_par","topograd",par%topograd)
    call nml_read(filename,"ice_sico_par","gl_surf_grad",par%gl_surf_grad)
    call nml_read(filename,"ice_sico_par","basal_hydrology",par%basal_hydrology)
    call nml_read(filename,"ice_sico_par","H_w_max",par%H_w_max)
    call nml_read(filename,"ice_sico_par","calcmod",par%calcmod)
    call nml_read(filename,"ice_sico_par","cts_melting_freezing",par%cts_melting_freezing)
    call nml_read(filename,"ice_sico_par","temp_const",par%temp_const)
    call nml_read(filename,"ice_sico_par","age_const",par%age_const)
    call nml_read(filename,"ice_sico_par","dynamics",par%dynamics)
    call nml_read(filename,"ice_sico_par","hyb_mode",par%hyb_mode)
    call nml_read(filename,"ice_sico_par","ice_control",par%ice_control)
    call nml_read(filename,"ice_sico_par","lis_opts",par%lis_opts)
    call nml_read(filename,"ice_sico_par","margin",par%margin)
    call nml_read(filename,"ice_sico_par","marine_ice_formation",par%marine_ice_formation)
    call nml_read(filename,"ice_sico_par","marine_ice_calving",par%marine_ice_calving)
    call nml_read(filename,"ice_sico_par","z_mar",par%z_mar)
    call nml_read(filename,"ice_sico_par","fact_z_mar",par%fact_z_mar)
    call nml_read(filename,"ice_sico_par","calv_uw_coeff",par%calv_uw_coeff)
    call nml_read(filename,"ice_sico_par","r1_calv_uw",par%r1_calv_uw)
    call nml_read(filename,"ice_sico_par","r2_calv_uw",par%r2_calv_uw)
    call nml_read(filename,"ice_sico_par","ice_shelf_calving",par%ice_shelf_calving)
    call nml_read(filename,"ice_sico_par","h_calv_deep",par%h_calv_deep)
    call nml_read(filename,"ice_sico_par","h_calv_shallow",par%h_calv_shallow)
    call nml_read(filename,"ice_sico_par","zl_deep",par%zl_deep)
    call nml_read(filename,"ice_sico_par","zl_shallow",par%zl_shallow)
    call nml_read(filename,"ice_sico_par","tau_calv",par%tau_calv)
    par%tau_calv = par%tau_calv*sec_year
    call nml_read(filename,"ice_sico_par","i_disc",par%i_disc)
    call nml_read(filename,"ice_sico_par","ii_disc",par%ii_disc)
    call nml_read(filename,"ice_sico_par","c_dis_bm",par%c_dis_bm)
    call nml_read(filename,"ice_sico_par","c_dis_clv",par%c_dis_clv)
    call nml_read(filename,"ice_sico_par","alpha_tc",par%alpha_tc)
    call nml_read(filename,"ice_sico_par","zl_std_crit_disc",par%zl_std_crit_disc)
    call nml_read(filename,"ice_sico_par","m_H",par%m_H)
    call nml_read(filename,"ice_sico_par","m_D",par%m_D)
    call nml_read(filename,"ice_sico_par","m_R",par%m_R)
    call nml_read(filename,"ice_sico_par","m_T",par%m_T)
    call nml_read(filename,"ice_sico_par","tau_dmb",par%tau_dmb)
    par%tau_dmb = par%tau_dmb*sec_year      ! yr -> sec
    call nml_read(filename,"ice_sico_par","r_mar_eff",par%r_mar_eff)
    par%r_mar_eff = par%r_mar_eff *1000.0_wp   ! km -> m
    call nml_read(filename,"ice_sico_par","i_dist_coast",par%i_dist_coast)
    call nml_read(filename,"ice_sico_par","dist_coast_min",par%dist_coast_min)
    par%dist_coast_min = par%dist_coast_min*1000._wp    ! m -> km
    call nml_read(filename,"ice_sico_par","dtime_mar_coa",par%dtime_mar_coa)
    if (abs(par%dtime_mar_coa/par%dtime-nint(par%dtime_mar_coa/par%dtime)).gt.1.e-10_wp) then
      print *,'ERROR, Sicopolis time step dtime_mar_coa has to be an integer multiple of dtime'
      print *,'dtime',par%dtime
      print *,'dtime_mar_coa',par%dtime_mar_coa
      stop
    endif
    call nml_read(filename,"ice_sico_par","fin_visc",par%fin_visc)
    call nml_read(filename,"ice_sico_par","flow_law",par%flow_law)
    call nml_read(filename,"ice_sico_par","enhmod",par%enhmod       )
    call nml_read(filename,"ice_sico_par","enh_fact",par%enh_fact     )
    call nml_read(filename,"ice_sico_par","enh_intg",par%enh_intg     )
    call nml_read(filename,"ice_sico_par","age_trans_0",par%age_trans_0  )
    call nml_read(filename,"ice_sico_par","date_trans1_0",par%date_trans1_0)
    call nml_read(filename,"ice_sico_par","date_trans2_0",par%date_trans2_0)
    call nml_read(filename,"ice_sico_par","date_trans3_0",par%date_trans3_0)
    call nml_read(filename,"ice_sico_par","enh_compr",par%enh_compr    )
    call nml_read(filename,"ice_sico_par","enh_shear",par%enh_shear    )
    call nml_read(filename,"ice_sico_par","enh_shelf",par%enh_shelf    )
    call nml_read(filename,"ice_sico_par","i_advance",par%i_advance)
    call nml_read(filename,"ice_sico_par","thk_evol",par%thk_evol)
    call nml_read(filename,"ice_sico_par","calcthk",par%calcthk)
    call nml_read(filename,"ice_sico_par","ovi_weight",par%ovi_weight)
    call nml_read(filename,"ice_sico_par","omega_sor",par%omega_sor)
    call nml_read(filename,"ice_sico_par","mean_accum",par%mean_accum)
    call nml_read(filename,"ice_sico_par","time_target_topo_init",par%time_target_topo_init)
    par%time_target_topo_init = par%time_target_topo_init*sec_year
    call nml_read(filename,"ice_sico_par","time_target_topo_final",par%time_target_topo_final)
    par%time_target_topo_final = par%time_target_topo_final*sec_year
    call nml_read(filename,"ice_sico_par","target_topo_tau",par%target_topo_tau)
    par%target_topo_tau = par%target_topo_tau*sec_year
    call nml_read(filename,"ice_sico_par","target_topo_file",par%target_topo_file)
    call nml_read(filename,"ice_sico_par","mask_maxextent_file",par%mask_maxextent_file)
    call nml_read(filename,"ice_sico_par","mb_account",par%mb_account)
    call nml_read(filename,"ice_sico_par","l_alpine",par%l_alpine)
    call nml_read(filename,"ice_sico_par","l_exclude_grl",par%l_exclude_grl)
    call nml_read(filename,"ice_sico_par","h_isol_max",par%h_isol_max)

    call nml_read(filename,"ice_sico_par","z_abyss",par%z_abyss)
    call nml_read(filename,"ice_sico_par","l_abyss",par%l_abyss)
    call nml_read(filename,"ice_sico_par","marine_ice_basal_melting",par%marine_ice_basal_melting)
    call nml_read(filename,"ice_sico_par","qbm_marine",par%qbm_marine)
    call nml_read(filename,"ice_sico_par","qbm_min",par%qbm_min)
    call nml_read(filename,"ice_sico_par","qbm_max",par%qbm_max)

    call nml_read(filename,"ice_sico_par","slide_law",par%slide_law       )
    call nml_read(filename,"ice_sico_par","c_slide",par%c_slide         )
    call nml_read(filename,"ice_sico_par","gamma_slide",par%gamma_slide     )
    call nml_read(filename,"ice_sico_par","p_weert",par%p_weert         )
    call nml_read(filename,"ice_sico_par","q_weert",par%q_weert         )
    call nml_read(filename,"ice_sico_par","c_slide_sedi",par%c_slide_sedi    )
    call nml_read(filename,"ice_sico_par","p_weert_sedi",par%p_weert_sedi    )
    call nml_read(filename,"ice_sico_par","q_weert_sedi",par%q_weert_sedi    )
    call nml_read(filename,"ice_sico_par","q_coulomb",par%q_coulomb)
    call nml_read(filename,"ice_sico_par","vb_t_rock",par%vb_t_rock)
    call nml_read(filename,"ice_sico_par","vb_t_sed",par%vb_t_sed)
    par%vb_t_rock = par%vb_t_rock / sec_year  ! m/a -> m/s
    par%vb_t_sed  = par%vb_t_sed  / sec_year  ! m/a -> m/s
    call nml_read(filename,"ice_sico_par","c_fric_rock",par%c_fric_rock)
    call nml_read(filename,"ice_sico_par","c_fric_sed",par%c_fric_sed)
    call nml_read(filename,"ice_sico_par","delta_rock",par%delta_rock)
    call nml_read(filename,"ice_sico_par","delta_sed",par%delta_sed)
    call nml_read(filename,"ice_sico_par","i_slide_rock_sed",par%i_slide_rock_sed    )
    call nml_read(filename,"ice_sico_par","h_sed_thresh",par%h_sed_thresh    )
    call nml_read(filename,"ice_sico_par","h_sed_min",par%h_sed_min    )
    call nml_read(filename,"ice_sico_par","h_sed_max",par%h_sed_max    )
    call nml_read(filename,"ice_sico_par","enh_slide_marine_sed",par%enh_slide_marine_sed )
    call nml_read(filename,"ice_sico_par","nyears_ramp_up_slide",par%nyears_ramp_up_slide)
    call nml_read(filename,"ice_sico_par","i_p_b_red",par%i_p_b_red)
    call nml_read(filename,"ice_sico_par","red_pres_limit_fact",par%red_pres_limit_fact)
    call nml_read(filename,"ice_sico_par","hydro_slide_sat_fct",par%hydro_slide_sat_fct)
    call nml_read(filename,"ice_sico_par","c_hw_slide",par%c_hw_slide)
    call nml_read(filename,"ice_sico_par","hw0_slide",par%hw0_slide)
    call nml_read(filename,"ice_sico_par","vh_max",par%vh_max)
    call nml_read(filename,"ice_sico_par","visc_min",par%visc_min)
    call nml_read(filename,"ice_sico_par","visc_max",par%visc_max)
    call nml_read(filename,"ice_sico_par","tol_iter_ssa ",par%tol_iter_ssa )
    call nml_read(filename,"ice_sico_par","n_iter_ssa   ",par%n_iter_ssa   )
    call nml_read(filename,"ice_sico_par","iter_init_ssa",par%iter_init_ssa)
    call nml_read(filename,"ice_sico_par","visc_init_ssa",par%visc_init_ssa)
    call nml_read(filename,"ice_sico_par","n_visc_smooth",par%n_visc_smooth)
    call nml_read(filename,"ice_sico_par","visc_smooth_diff  ",par%visc_smooth_diff  )
    call nml_read(filename,"ice_sico_par","relax_fact_ssa    ",par%relax_fact_ssa    )
    call nml_read(filename,"ice_sico_par","ratio_sl_thresh   ",par%ratio_sl_thresh   )
    call nml_read(filename,"ice_sico_par","ssta_sia_weigh_fct",par%ssta_sia_weigh_fct)
    call nml_read(filename,"ice_sico_par","hyb_ref_speed",par%hyb_ref_speed)

    call nml_read(filename,"ice_sico_par","hd_min",par%hd_min)
    call nml_read(filename,"ice_sico_par","hd_max",par%hd_max)

    call nml_read(filename,"ice_sico_par","numdiff_h_t",par%numdiff_h_t)
    call nml_read(filename,"ice_sico_par","tau_cts",par%tau_cts)

    call nml_read(filename,"ice_sico_par","agediff",par%agediff)
    call nml_read(filename,"ice_sico_par","age_min",par%age_min)
    call nml_read(filename,"ice_sico_par","age_max",par%age_max)

    call nml_read(filename,"ice_sico_par","q_litho",par%q_litho)
    call nml_read(filename,"ice_sico_par","adv_hor",par%adv_hor)
    call nml_read(filename,"ice_sico_par","adv_vert",par%adv_vert)

    call nml_read(filename,"ice_sico_par","l_write_timer",par%l_write_timer)
    call nml_read(filename,"ice_sico_par","l_output_3d",par%l_output_3d)

  par%calv_uw_coeff = par%calv_uw_coeff / sec_year

! regularized Coulomb friction law available only for dynamics==2 and hyb_mode==1
if (par%slide_law==4 .and. .not.(par%dynamics==2 .and. par%hyb_mode==1)) then 
  write(6, fmt='(a)') ' >>> sico_init: regularized Coulomb friction law '
  write(6, fmt='(a)') '                available only for dynamics==2 and hyb_mode==1'
  stop
endif

! pseudo-plastic power law friction law available only for dynamics==2 and hyb_mode==1
if (par%slide_law==5 .and. .not.(par%dynamics==2 .and. par%hyb_mode==1)) then 
  write(6, fmt='(a)') ' >>> sico_init: pseudo-plastic power law friction law '
  write(6, fmt='(a)') '                available only for dynamics==2 and hyb_mode==1'
  stop
endif

if ((par%margin==3 .and. par%dynamics==1) .or. par%dynamics==2) then !/* requires SSA or SStA */
if (par%grid.ne.0) then
  write(6, fmt='(a)') ' >>> sico_init: Distortion correction not implemented'
  write(6, fmt='(a)') '                for the shallow shelf approximation (SSA)'
  write(6, fmt='(a)') '                or the shelfy stream approximation (SStA)'
  write(6, fmt='(a)') '                -> GRID==0 required!'
  stop
endif
endif

if ((par%margin==3 .or. par%dynamics==2) .and. (par%calcthk==1 .or. par%calcthk==2 .or. par%calcthk==3)) then
    write(6, fmt='(a)') ' >>> sico_main_loop:'
    write(6, fmt='(a)') '          Non-SIA dynamics combined with'
    write(6, fmt='(a)') '          SIA ice thickness solver!'
    stop
endif

!-------- Compatibility check of discretization schemes for the horizontal and
!         vertical advection terms in the temperature and age equations --------

if (par%adv_hor==1) then
stop ' >>> sico_init: Option ADV_HOR==1 (central differences) not defined!'
endif

if (par%calcthk/=4 .and. par%i_advance==1) then
    write(6, fmt='(a)') ' >>> sico_par_load:'
    write(6, fmt='(a)') '          i_advance=1 only'
    write(6, fmt='(a)') '          for calcthk=4 implemented now!'
    stop
endif

   return

end subroutine sico_par_load

end module sico_params
