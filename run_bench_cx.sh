
############################################################################
# script to run a complete set of CLIMBER-X simulations for model benchmark
############################################################################

# NOTE: running this set of benchmark simulations requires restart from an existing preindustrial spinup experiment!

# NOTE: By default, job is called with `--omp 16` which means a parallel run with 16 processors

# specify desidered name of output directory
outdir=bench_v098


echo $outdir

# preindustrial control
#./runcx -s -w 24 -o output/$outdir/preind -q priority -p ctl.nyears=5000 ctl.flag_bgc=T bgc.l_spinup_bgc=F bgc.l_spinup_sed=F bgc.l_sediments=F lnd.l_weathering=F lnd.l_river_export=F lnd.l_burial=F

./runcx -s -q priority -w 24 -o output/$outdir/preind_bgc_open_nut -p ctl.nyears=5000 ctl.flag_bgc=T bgc.l_spinup_bgc=T bgc.l_spinup_sed=T bgc.l_sediments=T bgc.l_conserve_phos=T bgc.l_conserve_sil=T lnd.l_weathering=T lnd.l_river_export=F lnd.l_burial=F

# historical 
#./runcx -s -q priority -w 6 -o output/$outdir/hist       -p ctl.nyears=1015 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.icfc=1 ctl.io3=2 ctl.iso4=1 ctl.iluc=2 ctl.isol=1 ctl.ivolc=1 ctl.flag_ch4=T ctl.ich4_rad=2 ctl.ich4_emis=1 bgc.l_sediments=F lnd.l_weathering=F lnd.l_river_export=F lnd.l_burial=F ocn.l_cfc=T
#./runcx -s -q standby  -w 6 -o output/$outdir/hist-nat   -p ctl.nyears=1015 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2=0 ctl.ich4=0 ctl.in2o=0 ctl.icfc=0 ctl.io3=1 ctl.iso4=0 ctl.iluc=1 ctl.isol=1 ctl.ivolc=1 bgc.l_sediments=F lnd.l_weathering=F lnd.l_river_export=F lnd.l_burial=F
#./runcx -s -q priority -w 6 -o output/$outdir/hist-ghg   -p ctl.nyears=1015 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.icfc=1 ctl.io3=2 ctl.iso4=0 ctl.iluc=1 ctl.isol=0 ctl.ivolc=0 bgc.l_sediments=F lnd.l_weathering=F lnd.l_river_export=F lnd.l_burial=F
#./runcx -s -q standby  -w 6 -o output/$outdir/hist-aer   -p ctl.nyears=1015 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2=0 ctl.ich4=0 ctl.in2o=0 ctl.icfc=0 ctl.io3=1 ctl.iso4=1 ctl.iluc=1 ctl.isol=0 ctl.ivolc=0 bgc.l_sediments=F lnd.l_weathering=F lnd.l_river_export=F lnd.l_burial=F
#./runcx -s -q priority -w 6 -o output/$outdir/hist-noluc -p ctl.nyears=1015 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.icfc=1 ctl.io3=2 ctl.iso4=1 ctl.iluc=1 ctl.isol=1 ctl.ivolc=1 bgc.l_sediments=F lnd.l_weathering=F lnd.l_river_export=F lnd.l_burial=F
#./runcx -s -q priority -w 6 -o output/$outdir/hist-esm   -p ctl.nyears=1015 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.flag_co2=T ctl.ico2_emis=1 ctl.iorbit=1 ctl.ich4=1 ctl.in2o=1 ctl.icfc=1 ctl.io3=2 ctl.iso4=1 ctl.iluc=2 ctl.isol=1 ctl.ivolc=1

# 1/2xCO2
./runcx -s -q standby -w 48 -o output/$outdir/05xCO2  -p ctl.nyears=10000 ctl.co2_const=140
# 2xCO2
./runcx -s -q standby -w 48 -o output/$outdir/2xCO2   -p ctl.nyears=10000 ctl.co2_const=560
# 4xCO2
./runcx -s -q standby -w 48 -o output/$outdir/4xCO2   -p ctl.nyears=10000 ctl.co2_const=1120

# feedback analysis
./runcx -s -q standby -w 48 -o output/$outdir/fb/CO2_1x_2x  ctl.l_feedbacks=T ctl.nyears=10000 ctl.co2_const=280
./runcx -s -q standby -w 48 -o output/$outdir/fb/CO2_05x_1x ctl.l_feedbacks=T ctl.nyears=10000 ctl.co2_const=140
./runcx -s -q standby -w 48 -o output/$outdir/fb/CO2_2x_4x  ctl.l_feedbacks=T ctl.nyears=10000 ctl.co2_const=560

# abrupt 0p5xCO2
./runcx -s -q priority -w 1 -o output/$outdir/abrupt0p5xCO2 ctl.nyears=150 ctl.co2_const=140 ctl.nyout_atm=150 ctl.nyout_ocn=150 ctl.nyout_sic=150 ctl.nyout_lnd=150

# abrupt 2xCO2
./runcx -s -q priority -w 1 -o output/$outdir/abrupt2xCO2 ctl.nyears=150 ctl.co2_const=560 ctl.nyout_atm=150 ctl.nyout_ocn=150 ctl.nyout_sic=150 ctl.nyout_lnd=150

# abrupt 4xCO2
./runcx -s -q priority -w 1 -o output/$outdir/abrupt4xCO2 ctl.nyears=150 ctl.co2_const=1120 ctl.nyout_atm=150 ctl.nyout_ocn=150 ctl.nyout_sic=150 ctl.nyout_lnd=150

# abrupt 4% increase in solar constant
./runcx -s -q priority -w 1 -o output/$outdir/abruptsolp4p ctl.nyears=150 ctl.sol_const=1416.5 ctl.nyout_atm=150 ctl.nyout_ocn=150 ctl.nyout_sic=150 ctl.nyout_lnd=150

# 1%/yr CO2 increase
./runcx -s -q priority -w 1 -o output/$outdir/1pctCO2/cpl ctl.nyears=140 ctl.ico2=2 ctl.ico2_rad=0 ctl.co2_const=280 ctl.flag_bgc=T ctl.nyout_atm=140 ctl.nyout_ocn=140 ctl.nyout_sic=140 ctl.nyout_lnd=140
./runcx -s -q priority -w 1 -o output/$outdir/1pctCO2/bgc ctl.nyears=140 ctl.ico2=2 ctl.ico2_rad=1 ctl.co2_const=280 ctl.flag_bgc=T ctl.nyout_atm=140 ctl.nyout_ocn=140 ctl.nyout_sic=140 ctl.nyout_lnd=140
./runcx -s -q priority -w 1 -o output/$outdir/1pctCO2/rad ctl.nyears=140 ctl.ico2=0 ctl.ico2_rad=3 ctl.co2_const=280 ctl.flag_bgc=T ctl.nyout_atm=140 ctl.nyout_ocn=140 ctl.nyout_sic=140 ctl.nyout_lnd=140
#./runcx -s -q standby -w 24 -o output/$outdir/1pctCO2_long ctl.nyears=5000 ctl.ico2=2 ctl.ico2_rad=0 ctl.co2_const=280

# aquaplanet
./runcx -s -q priority -w 1 -o output/$outdir/aqua ctl.nyears=10 ctl.l_aquaplanet=T ctl.flag_ocn=F ctl.flag_sic=F ctl.flag_lnd=F ctl.flag_bgc=F ctl.nyout_atm=10

# LGM
./runcx -s -q priority -w 24 -o output/$outdir/lgm ctl.nyears=5000 ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc ctl.fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ctl.co2_const=190 ctl.ch4_const=375 ctl.n2o_const=200

./runcx -s -q standby -w 48 -o output/$outdir/lgm_bgc ctl.nyears=10000 ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 ctl.flag_co2=T ctl.flag_bgc=T ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc ctl.fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ctl.ico2_rad=1 ctl.co2_rad_const=190 ctl.ch4_const=375 ctl.n2o_const=200 bgc.l_spinup_bgc=F bgc.l_spinup_sed=F bgc.l_sediments=F lnd.l_weathering=F lnd.l_river_export=F lnd.l_burial=F

# mid-Holocene
./runcx -s -q standby -w 24 -o output/$outdir/midHolo ctl.nyears=5000 ctl.iorbit=1 ctl.ecc_const=0.018682 ctl.obl_const=24.105 ctl.per_const=0.87 ctl.co2_const=264.4 ctl.ch4_const=597 ctl.n2o_const=262

# last interglacial
./runcx -s -q standby -w 24 -o output/$outdir/lig127k ctl.nyears=5000 ctl.iorbit=1 ctl.ecc_const=0.039378 ctl.obl_const=24.04 ctl.per_const=275.41 ctl.co2_const=275 ctl.ch4_const=685 ctl.n2o_const=255

./runcx -s -q standby -w 48 -o output/$outdir/lig127k_bgc ctl.flag_bgc=T ctl.nyears=10000 ctl.iorbit=1 ctl.ecc_const=0.039378 ctl.obl_const=24.04 ctl.per_const=275.41 ctl.co2_const=275 ctl.ch4_const=685 ctl.n2o_const=255

# inception with interactive CO2 
#./runcx -s -q standby -o output/$outdir/incep_co2 ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=275 ctl.flag_co2=T ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_imo=T

# deglaciation
./runcx -s -q medium -w 168 -o output/$outdir/deglac ctl.nyears=25000 ctl.year_ini=-25000 ctl.ifake_ice=1 ctl.ifake_geo=1 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.iorbit=2

# hysteresis
./runcx -s -q medium -w 168 -o output/$outdir/hyst2050/up   ctl.nyears=20000 ocn.l_hosing=T ocn.hosing=-0.5 ocn.hosing_trend=0.05  ocn.lat_min_hosing=20 ocn.lat_max_hosing=50
./runcx -s -q medium -w 168 -o output/$outdir/hyst2050/down ctl.nyears=20000 ocn.l_hosing=T ocn.hosing=0.5  ocn.hosing_trend=-0.05 ocn.lat_min_hosing=20 ocn.lat_max_hosing=50
./runcx -s -q medium -w 168 -o output/$outdir/hyst5070/up   ctl.nyears=20000 ocn.l_hosing=T ocn.hosing=-0.5 ocn.hosing_trend=0.05  ocn.lat_min_hosing=50 ocn.lat_max_hosing=70
./runcx -s -q medium -w 168 -o output/$outdir/hyst5070/down ctl.nyears=20000 ocn.l_hosing=T ocn.hosing=0.5  ocn.hosing_trend=-0.05 ocn.lat_min_hosing=50 ocn.lat_max_hosing=70

# preind with interactive ice sheets
./runcx -s -q standby -w 24 -o output/$outdir/preind_ice_sia    ctl.nyears=30000 ctl.n_accel=100 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_imo=T ctl.n_ice_domain=1 ice.dynamics=1 ice.margin=2 ice.calcthk=2 smb.l_t2m_obs=F 
./runcx -s -q standby -w 24 -o output/$outdir/preind_ice_hybrid ctl.nyears=30000 ctl.n_accel=100 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_imo=T ctl.n_ice_domain=1 ice.dynamics=2 ice.margin=3 ice.calcthk=6 smb.l_t2m_obs=F 

# preind with semi 
#./runcx -s -q standby -w 10 -o output/$outdir/preind_smb ctl.nyears=100 ctl.flag_smb=T ctl.n_ice_domain=3 ctl.nyout_smb=100 smb.i_smb=1 smb.grid1_name=GRL-10KM smb.grid2_name=EU-10KM smb.grid3_name=ANT-20KM smb.l_monthly_output=T

# PaleoDEM
./runcx -s -q standby -w 24 -o output/$outdir/65Ma  ctl.nyears=5000 ctl.fake_geo_const_file=input/topog_065Ma_climberX.nc ctl.fake_ice_const_file=input/topog_065Ma_climberX.nc ctl.atm_restart=F ctl.ocn_restart=F ctl.sic_restart=F ctl.lnd_restart=F ocn.l_isl_auto=T ocn.i_init=1 geo.l_close_panama=F geo.geo_ref_file=input/topog_065Ma_climberX.nc 
./runcx -s -q standby -w 24 -o output/$outdir/250Ma ctl.nyears=5000 ctl.fake_geo_const_file=input/topog_250Ma_climberX.nc ctl.fake_ice_const_file=input/topog_250Ma_climberX.nc ctl.atm_restart=F ctl.ocn_restart=F ctl.sic_restart=F ctl.lnd_restart=F ocn.l_isl_auto=T ocn.i_init=1 geo.l_close_panama=F geo.geo_ref_file=input/topog_250Ma_climberX.nc 

# CO2 pulse 
./runcx -s -q standby -w 24 -o output/$outdir/co2_pulse/1000PgC      ctl.flag_co2=T ctl.flag_bgc=T ctl.nyears=2000 ctl.ico2_emis=2 ctl.co2_pulse=1000 ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100

./runcx -s -q standby -w 24 -o output/$outdir/co2_pulse_hist/1000PgC ctl.flag_co2=T ctl.nyears=3000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2_emis=1 ctl.ich4=1 ctl.in2o=1 ctl.icfc=1 ctl.io3=2 ctl.iso4=1 ctl.iluc=2 ctl.isol=1 ctl.ivolc=1 ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100











# Previously with `-a`, but set to `-o` with runcx for now:
#./runcx -s -q priority -w 24 -o output/cc_pi_lgm/pi  ctl.nyears=5000 ctl.flag_bgc=T bgc.l_spinup_bgc=T bgc.l_spinup_sed=T lnd.l_weathering=F lnd.l_river_export=F

#./runcx -s -q priority -w 24 -o output/cc_pi_lgm/lgm ctl.nyears=5000 ctl.year_ini=-21000 ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc ctl.fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ctl.flag_co2=T ctl.flag_bgc=T ctl.ocn_restart=T ctl.lnd_restart=T ctl.sic_restart=F ctl.bgc_restart=T ctl.co2_degas_const=0 bgc.l_spinup_bgc=T lnd.l_weathering=F lnd.l_river_export=F

# closed carbon cycle spinup
#./runcx -s -q priority -w 24 -o output/spin_cc/close ctl.l_spinup_cc=T ctl.nyears_spinup_cc=10000 ctl.flag_bgc=T ctl.nyears=20000 bgc.l_spinup_bgc=F bgc.l_spinup_sed=F bgc.l_sediments=F lnd.l_weathering=F lnd.l_river_export=F lnd.l_burial=F

# open carbon cycle spinup but with conserved nutrients (phosphate and silicate)
#./runcx -s -q priority -w 24 -o output/spin_cc/open_nut ctl.l_spinup_cc=T ctl.nyears_spinup_cc=10000 ctl.flag_bgc=T ctl.nyears=20000 bgc.l_spinup_bgc=T bgc.l_spinup_sed=T bgc.l_sediments=T bgc.l_conserve_phos=T bgc.l_conserve_sil=T lnd.l_weathering=T lnd.l_river_export=F lnd.l_burial=F

# open carbon cycle spinup but with conserved nutrients (phosphate and silicate) and alkalinity
#./runcx -s -q standby -w 72 -o output/spin_cc_open_nut/alkT_k10000 ctl.l_spinup_cc=T ctl.nyears_spinup_cc=15000 ctl.flag_bgc=T ctl.nyears=20000 bgc.l_spinup_bgc=T bgc.l_spinup_sed=T bgc.l_sediments=T bgc.l_conserve_phos=T bgc.l_conserve_sil=T bgc.l_conserve_alk=T lnd.l_weathering=T lnd.l_river_export=F lnd.l_burial=F lnd.k_min=10000

# open carbon cycle spinup
#./runcx -s -q priority -w 24 -a output/spin_cc/open ctl.l_spinup_cc=T ctl.nyears_spinup_cc=10000 ctl.flag_bgc=T ctl.nyears=20000 bgc.l_spinup_bgc=T bgc.l_spinup_sed=T bgc.l_sediments=T bgc.l_conserve_phos=F bgc.l_conserve_sil=F lnd.l_weathering=T lnd.l_river_export=T lnd.l_burial=F

# LGM closed carbon cycle
#./runcx -s-q standby -o output/cc_lgm ctl.nyears=10000 ctl.year_ini=-21000 ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc ctl.fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ctl.flag_co2=T ctl.flag_bgc=T ctl.co2_degas_const=0 ictl.co2_rad=1 ctl.co2_rad_const=180 bgc.l_spinup_bgc=F bgc.l_spinup_sed=F bgc.l_sediments=F lnd.l_weathering=F lnd.l_river_export=F lnd.l_burial=F

