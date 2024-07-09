!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c a r b o n _ t r a n s _ m o d
!
!  Purpose : transfer of carbon between pools when changing fractions
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
module carbon_trans_mod
  
  use precision, only : wp
  use timer, only : day_year
  use control, only : check_carbon
  use lnd_grid, only : dz_c, rdz_c, nlc
  use lnd_params, only : dt_c, dt_day_c
  use lnd_params, only : peat_par 

  implicit none

  private
  public :: carbon_trans

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c a r b o n _ t r a n s
  !   Purpose    :  transfer carbon between different pools
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine carbon_trans(f_veg, f_veg_old, f_ice_grd, f_ice_grd_old, f_lake, f_lake_old, f_shelf, f_shelf_old, &
      dCpeat_dt, f_peat_pot, f_peat,&
      litter_c, fast_c, slow_c, litter_c13, fast_c13, slow_c13, litter_c14, fast_c14, slow_c14, &
      litter_c_shelf, fast_c_shelf, slow_c_shelf, litter_c13_shelf, fast_c13_shelf, slow_c13_shelf, litter_c14_shelf, fast_c14_shelf, slow_c14_shelf, &
      litter_c_lake, fast_c_lake, slow_c_lake, litter_c13_lake, fast_c13_lake, slow_c13_lake, litter_c14_lake, fast_c14_lake, slow_c14_lake, &
      litter_c_ice, fast_c_ice, slow_c_ice, litter_c13_ice, fast_c13_ice, slow_c13_ice, litter_c14_ice, fast_c14_ice, slow_c14_ice, &
      litter_c_peat, acro_c, cato_c, litter_c13_peat, acro_c13, cato_c13, litter_c14_peat, acro_c14, cato_c14)

    implicit none

    real(wp), intent(in) :: f_veg, f_veg_old, f_ice_grd, f_ice_grd_old, f_lake, f_lake_old, f_shelf, f_shelf_old 
    real(wp), intent(in) :: dCpeat_dt, f_peat_pot
    real(wp), intent(inout) :: f_peat
    real(wp), dimension(:), intent(inout) :: litter_c, fast_c, slow_c, litter_c13, fast_c13, slow_c13, litter_c14, fast_c14, slow_c14
    real(wp), dimension(:), intent(inout) :: litter_c_shelf, fast_c_shelf, slow_c_shelf, litter_c13_shelf, fast_c13_shelf, slow_c13_shelf, litter_c14_shelf, fast_c14_shelf, slow_c14_shelf
    real(wp), dimension(:), intent(inout) :: litter_c_lake, fast_c_lake, slow_c_lake, litter_c13_lake, fast_c13_lake, slow_c13_lake, litter_c14_lake, fast_c14_lake, slow_c14_lake
    real(wp), dimension(:), intent(inout) :: litter_c_ice, fast_c_ice, slow_c_ice, litter_c13_ice, fast_c13_ice, slow_c13_ice, litter_c14_ice, fast_c14_ice, slow_c14_ice
    real(wp), dimension(:), intent(inout) :: cato_c, cato_c13, cato_c14
    real(wp), intent(inout) :: litter_c_peat, acro_c, litter_c13_peat, acro_c13, litter_c14_peat, acro_c14

    real(wp), dimension(nlc) :: litter_c_old, fast_c_old, slow_c_old, litter_c13_old, fast_c13_old, slow_c13_old, litter_c14_old, fast_c14_old, slow_c14_old
    real(wp), dimension(nlc) :: litter_c_shelf_old, fast_c_shelf_old, slow_c_shelf_old, litter_c13_shelf_old, fast_c13_shelf_old, slow_c13_shelf_old, litter_c14_shelf_old, fast_c14_shelf_old, slow_c14_shelf_old
    real(wp), dimension(nlc) :: litter_c_lake_old, fast_c_lake_old, slow_c_lake_old, litter_c13_lake_old, fast_c13_lake_old, slow_c13_lake_old, litter_c14_lake_old, fast_c14_lake_old, slow_c14_lake_old
    real(wp), dimension(nlc) :: litter_c_ice_old, fast_c_ice_old, slow_c_ice_old, litter_c13_ice_old, fast_c13_ice_old, slow_c13_ice_old, litter_c14_ice_old, fast_c14_ice_old, slow_c14_ice_old
    real(wp), dimension(nlc) :: cato_c_old, cato_c13_old, cato_c14_old
    real(wp) :: litter_c_peat_old, acro_c_old, litter_c13_peat_old, acro_c13_old, litter_c14_peat_old, acro_c14_old

    logical :: flag_done
    real(wp) :: peat_c, f_peat_old, d_f_peat
    real(wp) :: d_f_veg, d_f_ice, d_f_lake, d_f_shelf, fac, dt_p
    real(wp) :: carbon_tot_old, carbon13_tot_old, carbon14_tot_old, carbon_tot, carbon13_tot, carbon14_tot


    ! save old carbon values
    litter_c_old = litter_c
    fast_c_old   = fast_c
    slow_c_old   = slow_c
    litter_c13_old = litter_c13
    fast_c13_old   = fast_c13
    slow_c13_old   = slow_c13
    litter_c14_old = litter_c14
    fast_c14_old   = fast_c14
    slow_c14_old   = slow_c14

    litter_c_shelf_old = litter_c_shelf
    fast_c_shelf_old   = fast_c_shelf
    slow_c_shelf_old   = slow_c_shelf
    litter_c13_shelf_old = litter_c13_shelf
    fast_c13_shelf_old   = fast_c13_shelf
    slow_c13_shelf_old   = slow_c13_shelf
    litter_c14_shelf_old = litter_c14_shelf
    fast_c14_shelf_old   = fast_c14_shelf
    slow_c14_shelf_old   = slow_c14_shelf

    litter_c_lake_old = litter_c_lake
    fast_c_lake_old   = fast_c_lake
    slow_c_lake_old   = slow_c_lake
    litter_c13_lake_old = litter_c13_lake
    fast_c13_lake_old   = fast_c13_lake
    slow_c13_lake_old   = slow_c13_lake
    litter_c14_lake_old = litter_c14_lake
    fast_c14_lake_old   = fast_c14_lake
    slow_c14_lake_old   = slow_c14_lake

    litter_c_ice_old = litter_c_ice
    fast_c_ice_old   = fast_c_ice
    slow_c_ice_old   = slow_c_ice
    litter_c13_ice_old = litter_c13_ice
    fast_c13_ice_old   = fast_c13_ice
    slow_c13_ice_old   = slow_c13_ice
    litter_c14_ice_old = litter_c14_ice
    fast_c14_ice_old   = fast_c14_ice
    slow_c14_ice_old   = slow_c14_ice

    litter_c_peat_old = litter_c_peat
    acro_c_old        = acro_c
    cato_c_old        = cato_c
    litter_c13_peat_old = litter_c13_peat
    acro_c13_old        = acro_c13
    cato_c13_old        = cato_c13
    litter_c14_peat_old = litter_c14_peat
    acro_c14_old        = acro_c14
    cato_c14_old        = cato_c14


    ! 1) transfer carbon between peat and mineral soil carbon pools if peat fraction or vegetated fraction is changing

    f_peat_old = f_peat
    if( peat_par%peat_area ) then

      dt_p = dt_c*day_year/real(dt_day_c,wp)

     ! changes in peatland area
     peat_c = litter_c_peat+acro_c+sum(cato_c*dz_c(1:nlc)) ! kg/m2
     if( dCpeat_dt .gt. peat_par%dCpeat_dt_min  .or. peat_c .gt. peat_par%Cpeat_min ) then 
       ! peatland expansion
       f_peat = min((1._wp+peat_par%peat_ch_rate*dt_p)*f_peat,f_peat_pot*f_veg)
     else ! peatland contraction
       f_peat = max((1._wp-peat_par%peat_ch_rate*10._wp/max(1e-20_wp,peat_c)*dt_p)*f_peat,peat_par%f_peat_min*f_veg) ! contraction rate depends on peat carbon content
       if (f_veg.eq.0._wp) f_peat = 0._wp
     endif
     ! when ice or water area expanding reduce peatland fraction proportionally
     if( f_veg .lt. f_veg_old .and. f_veg_old .gt. 0._wp ) then
      if( f_veg_old .eq. 0._wp ) print *,'ATT f_veg_old<0',f_veg_old,f_veg
      f_peat = f_peat * f_veg/f_veg_old 
      !f_peat = max(f_peat,peat_par%f_peat_min)
     endif

     ! shift carbon between peat and mineral soil fraction to adjust for peatland area changes
     d_f_peat = f_peat-f_peat_old
     if( d_f_peat .gt. 0._wp ) then ! peat is expanding, transfer carbon to peat
      litter_c_peat = (f_peat_old*litter_c_peat + litter_c(1)*dz_c(1)*d_f_peat)/f_peat  ! kgC/m2
      acro_c        = (f_peat_old*acro_c        + (fast_c(1)+slow_c(1))*dz_c(1)*d_f_peat)/f_peat ! kgC/m2
      cato_c(2:nlc)  = (f_peat_old*cato_c(2:nlc)  + (litter_c(2:nlc)+fast_c(2:nlc)+slow_c(2:nlc))*d_f_peat)/f_peat  ! kgC/m3
      ! C13
      litter_c13_peat = (f_peat_old*litter_c13_peat + litter_c13(1)*dz_c(1)*d_f_peat)/f_peat  ! kgC/m2
      acro_c13        = (f_peat_old*acro_c13        + (fast_c13(1)+slow_c13(1))*dz_c(1)*d_f_peat)/f_peat ! kgC/m2
      cato_c13(2:nlc)  = (f_peat_old*cato_c13(2:nlc)  + (litter_c13(2:nlc)+fast_c13(2:nlc)+slow_c13(2:nlc))*d_f_peat)/f_peat  ! kgC/m3
      ! C14
      litter_c14_peat = (f_peat_old*litter_c14_peat + litter_c14(1)*dz_c(1)*d_f_peat)/f_peat  ! kgC/m2
      acro_c14        = (f_peat_old*acro_c14        + (fast_c14(1)+slow_c14(1))*dz_c(1)*d_f_peat)/f_peat ! kgC/m2
      cato_c14(2:nlc)  = (f_peat_old*cato_c14(2:nlc)  + (litter_c14(2:nlc)+fast_c14(2:nlc)+slow_c14(2:nlc))*d_f_peat)/f_peat  ! kgC/m3
     elseif( d_f_peat .lt. 0._wp ) then ! peat is contracting, transfer peat to mineral carbon pools
      litter_c(1)   = ((f_veg_old-f_peat_old)*litter_c(1)   - litter_c_peat*rdz_c(1)*d_f_peat)/(f_veg_old-f_peat)  ! kgC/m3
      fast_c(1)     = ((f_veg_old-f_peat_old)*fast_c(1)     - acro_c*rdz_c(1)*d_f_peat)/(f_veg_old-f_peat)  ! kgC/m3
      slow_c(2:nlc)  = ((f_veg_old-f_peat_old)*slow_c(2:nlc)  - 0.5_wp*cato_c(2:nlc)*d_f_peat)/(f_veg_old-f_peat)  ! kgC/m3
      litter_c(2:nlc)= (f_veg_old-f_peat_old)*litter_c(2:nlc)/(f_veg_old-f_peat)
      fast_c(2:nlc)  = ((f_veg_old-f_peat_old)*fast_c(2:nlc)  - 0.5_wp*cato_c(2:nlc)*d_f_peat)/(f_veg_old-f_peat)
      slow_c(1)     = (f_veg_old-f_peat_old)*slow_c(1)/(f_veg_old-f_peat)
      ! C13
      litter_c13(1)   = ((f_veg_old-f_peat_old)*litter_c13(1)   - litter_c13_peat*rdz_c(1)*d_f_peat)/(f_veg_old-f_peat)  ! kgC/m3
      fast_c13(1)     = ((f_veg_old-f_peat_old)*fast_c13(1)     - acro_c13*rdz_c(1)*d_f_peat)/(f_veg_old-f_peat)  ! kgC/m3
      slow_c13(2:nlc)  = ((f_veg_old-f_peat_old)*slow_c13(2:nlc)  - 0.5_wp*cato_c13(2:nlc)*d_f_peat)/(f_veg_old-f_peat)  ! kgC/m3
      litter_c13(2:nlc)= (f_veg_old-f_peat_old)*litter_c13(2:nlc)/(f_veg_old-f_peat)
      fast_c13(2:nlc)  = ((f_veg_old-f_peat_old)*fast_c13(2:nlc)  - 0.5_wp*cato_c13(2:nlc)*d_f_peat)/(f_veg_old-f_peat)
      slow_c13(1)     = (f_veg_old-f_peat_old)*slow_c13(1)/(f_veg_old-f_peat)
      ! C14
      litter_c14(1)   = ((f_veg_old-f_peat_old)*litter_c14(1)   - litter_c14_peat*rdz_c(1)*d_f_peat)/(f_veg_old-f_peat)  ! kgC/m3
      fast_c14(1)     = ((f_veg_old-f_peat_old)*fast_c14(1)     - acro_c14*rdz_c(1)*d_f_peat)/(f_veg_old-f_peat)  ! kgC/m3
      slow_c14(2:nlc)  = ((f_veg_old-f_peat_old)*slow_c14(2:nlc)  - 0.5_wp*cato_c14(2:nlc)*d_f_peat)/(f_veg_old-f_peat)  ! kgC/m3
      litter_c14(2:nlc)= (f_veg_old-f_peat_old)*litter_c14(2:nlc)/(f_veg_old-f_peat)
      fast_c14(2:nlc)  = ((f_veg_old-f_peat_old)*fast_c14(2:nlc)  - 0.5_wp*cato_c14(2:nlc)*d_f_peat)/(f_veg_old-f_peat)
      slow_c14(1)     = (f_veg_old-f_peat_old)*slow_c14(1)/(f_veg_old-f_peat)
     endif

     carbon_tot_old = (f_veg_old-f_peat_old) * (sum(litter_c_old*dz_c(1:nlc)) + sum(fast_c_old*dz_c(1:nlc)) + sum(slow_c_old*dz_c(1:nlc))) &
                    + f_peat_old * (litter_c_peat_old+acro_c_old+sum(cato_c_old*dz_c(1:nlc))) 
     carbon_tot = (f_veg_old-f_peat) * (sum(litter_c*dz_c(1:nlc)) + sum(fast_c*dz_c(1:nlc)) + sum(slow_c*dz_c(1:nlc))) &
                + f_peat * (litter_c_peat+acro_c+sum(cato_c*dz_c(1:nlc)))

     carbon13_tot_old = (f_veg_old-f_peat_old) * (sum(litter_c13_old*dz_c(1:nlc)) + sum(fast_c13_old*dz_c(1:nlc)) + sum(slow_c13_old*dz_c(1:nlc))) &
                    + f_peat_old * (litter_c13_peat_old+acro_c13_old+sum(cato_c13_old*dz_c(1:nlc)))
     carbon13_tot = (f_veg_old-f_peat) * (sum(litter_c13*dz_c(1:nlc)) + sum(fast_c13*dz_c(1:nlc)) + sum(slow_c13*dz_c(1:nlc))) &
                + f_peat * (litter_c13_peat+acro_c13+sum(cato_c13*dz_c(1:nlc)))

     carbon14_tot_old = (f_veg_old-f_peat_old) * (sum(litter_c14_old*dz_c(1:nlc)) + sum(fast_c14_old*dz_c(1:nlc)) + sum(slow_c14_old*dz_c(1:nlc))) &
                    + f_peat_old * (litter_c14_peat_old+acro_c14_old+sum(cato_c14_old*dz_c(1:nlc)))
     carbon14_tot = (f_veg_old-f_peat) * (sum(litter_c14*dz_c(1:nlc)) + sum(fast_c14*dz_c(1:nlc)) + sum(slow_c14*dz_c(1:nlc))) &
                + f_peat * (litter_c14_peat+acro_c14+sum(cato_c14*dz_c(1:nlc)))
 
     if (check_carbon) then
     if( abs(carbon_tot-carbon_tot_old) .gt. 1.e-5_wp ) then
      print *,'carbon conservation after PEAT transfer ',carbon_tot-carbon_tot_old
      print *,f_peat_old,f_peat
     endif

     if( abs(carbon13_tot-carbon13_tot_old) .gt. 1.e-5_wp ) then
      print *,'carbon conservation after PEAT transfer ',carbon_tot-carbon_tot_old,carbon_tot,carbon_tot_old
      print *,'carbon 13 conservation after PEAT transfer ',carbon13_tot-carbon13_tot_old,carbon13_tot,carbon13_tot_old
      print *,f_peat_old,f_peat
      print *,f_veg_old,f_veg
     endif

     if( abs(carbon14_tot-carbon14_tot_old) .gt. 1.e-5_wp ) then
      print *,'carbon 14 conservation after PEAT transfer ',carbon14_tot-carbon14_tot_old
      print *,f_peat_old,f_peat
     endif
     endif


    endif

     ! 2) if shelf or ice or lake area is changing, 

     d_f_veg = f_veg - f_veg_old  ! change in vegetated fraction
     d_f_ice = f_ice_grd - f_ice_grd_old  ! change in ice fraction
     d_f_lake = f_lake - f_lake_old  ! change in lake fraction
     d_f_shelf = f_shelf - f_shelf_old ! change in shelf fraction

     flag_done = .false.

     ! vegetation replacing ice and/or shelf and/or lake 
     if( (.not. flag_done) .and. d_f_veg .gt. 0._wp .and. d_f_shelf .le. 0._wp .and. d_f_ice .le. 0._wp .and. d_f_lake .le. 0._wp ) then 
       flag_done = .true.
      ! take carbon from ice and/or shelf and/or lake and add it to the vegetation part
      litter_c   = ((f_veg_old-f_peat)*litter_c   - d_f_ice*litter_c_ice   - d_f_shelf*litter_c_shelf   - d_f_lake*litter_c_lake  ) / (f_veg-f_peat)
      fast_c     = ((f_veg_old-f_peat)*fast_c     - d_f_ice*fast_c_ice     - d_f_shelf*fast_c_shelf     - d_f_lake*fast_c_lake    ) / (f_veg-f_peat)
      slow_c     = ((f_veg_old-f_peat)*slow_c     - d_f_ice*slow_c_ice     - d_f_shelf*slow_c_shelf     - d_f_lake*slow_c_lake    ) / (f_veg-f_peat)
      litter_c13 = ((f_veg_old-f_peat)*litter_c13 - d_f_ice*litter_c13_ice - d_f_shelf*litter_c13_shelf - d_f_lake*litter_c13_lake) / (f_veg-f_peat)
      fast_c13   = ((f_veg_old-f_peat)*fast_c13   - d_f_ice*fast_c13_ice   - d_f_shelf*fast_c13_shelf   - d_f_lake*fast_c13_lake  ) / (f_veg-f_peat)
      slow_c13   = ((f_veg_old-f_peat)*slow_c13   - d_f_ice*slow_c13_ice   - d_f_shelf*slow_c13_shelf   - d_f_lake*slow_c13_lake  ) / (f_veg-f_peat)
      litter_c14 = ((f_veg_old-f_peat)*litter_c14 - d_f_ice*litter_c14_ice - d_f_shelf*litter_c14_shelf - d_f_lake*litter_c14_lake) / (f_veg-f_peat)
      fast_c14   = ((f_veg_old-f_peat)*fast_c14   - d_f_ice*fast_c14_ice   - d_f_shelf*fast_c14_shelf   - d_f_lake*fast_c14_lake  ) / (f_veg-f_peat)
      slow_c14   = ((f_veg_old-f_peat)*slow_c14   - d_f_ice*slow_c14_ice   - d_f_shelf*slow_c14_shelf   - d_f_lake*slow_c14_lake  ) / (f_veg-f_peat)
     endif

     ! ice replacing vegetation and/or shelf and/or lake
     if( (.not. flag_done) .and. d_f_ice .gt. 0._wp .and. d_f_veg .le. 0._wp .and. d_f_shelf .le. 0._wp .and. d_f_lake .le. 0._wp ) then
       flag_done = .true.
      ! take carbon from vegetated part and shelf and lake and add it to ice carbon
      litter_c_ice   = (f_ice_grd_old*litter_c_ice   - d_f_veg*litter_c   - d_f_shelf*litter_c_shelf   - d_f_lake*litter_c_lake  ) / f_ice_grd
      fast_c_ice     = (f_ice_grd_old*fast_c_ice     - d_f_veg*fast_c     - d_f_shelf*fast_c_shelf     - d_f_lake*fast_c_lake    ) / f_ice_grd
      slow_c_ice     = (f_ice_grd_old*slow_c_ice     - d_f_veg*slow_c     - d_f_shelf*slow_c_shelf     - d_f_lake*slow_c_lake    ) / f_ice_grd
      litter_c13_ice = (f_ice_grd_old*litter_c13_ice - d_f_veg*litter_c13 - d_f_shelf*litter_c13_shelf - d_f_lake*litter_c13_lake) / f_ice_grd
      fast_c13_ice   = (f_ice_grd_old*fast_c13_ice   - d_f_veg*fast_c13   - d_f_shelf*fast_c13_shelf   - d_f_lake*fast_c13_lake  ) / f_ice_grd
      slow_c13_ice   = (f_ice_grd_old*slow_c13_ice   - d_f_veg*slow_c13   - d_f_shelf*slow_c13_shelf   - d_f_lake*slow_c13_lake  ) / f_ice_grd
      litter_c14_ice = (f_ice_grd_old*litter_c14_ice - d_f_veg*litter_c14 - d_f_shelf*litter_c14_shelf - d_f_lake*litter_c14_lake) / f_ice_grd
      fast_c14_ice   = (f_ice_grd_old*fast_c14_ice   - d_f_veg*fast_c14   - d_f_shelf*fast_c14_shelf   - d_f_lake*fast_c14_lake  ) / f_ice_grd
      slow_c14_ice   = (f_ice_grd_old*slow_c14_ice   - d_f_veg*slow_c14   - d_f_shelf*slow_c14_shelf   - d_f_lake*slow_c14_lake  ) / f_ice_grd
     endif

     ! shelf replacing vegetation and/or ice and/or lake
     if( (.not. flag_done) .and. d_f_shelf .gt. 0._wp .and. d_f_veg .le. 0._wp .and. d_f_ice .le. 0._wp .and. d_f_lake .le. 0._wp ) then
       flag_done = .true.
      ! take carbon from vegetated part and shelf and lake and add it to ice carbon
      litter_c_shelf   = (f_shelf_old*litter_c_shelf   - d_f_veg*litter_c   - d_f_ice*litter_c_ice   - d_f_lake*litter_c_lake  ) / f_shelf 
      fast_c_shelf     = (f_shelf_old*fast_c_shelf     - d_f_veg*fast_c     - d_f_ice*fast_c_ice     - d_f_lake*fast_c_lake    ) / f_shelf
      slow_c_shelf     = (f_shelf_old*slow_c_shelf     - d_f_veg*slow_c     - d_f_ice*slow_c_ice     - d_f_lake*slow_c_lake    ) / f_shelf
      litter_c13_shelf = (f_shelf_old*litter_c13_shelf - d_f_veg*litter_c13 - d_f_ice*litter_c13_ice - d_f_lake*litter_c13_lake) / f_shelf
      fast_c13_shelf   = (f_shelf_old*fast_c13_shelf   - d_f_veg*fast_c13   - d_f_ice*fast_c13_ice   - d_f_lake*fast_c13_lake  ) / f_shelf
      slow_c13_shelf   = (f_shelf_old*slow_c13_shelf   - d_f_veg*slow_c13   - d_f_ice*slow_c13_ice   - d_f_lake*slow_c13_lake  ) / f_shelf
      litter_c14_shelf = (f_shelf_old*litter_c14_shelf - d_f_veg*litter_c14 - d_f_ice*litter_c14_ice - d_f_lake*litter_c14_lake) / f_shelf
      fast_c14_shelf   = (f_shelf_old*fast_c14_shelf   - d_f_veg*fast_c14   - d_f_ice*fast_c14_ice   - d_f_lake*fast_c14_lake  ) / f_shelf
      slow_c14_shelf   = (f_shelf_old*slow_c14_shelf   - d_f_veg*slow_c14   - d_f_ice*slow_c14_ice   - d_f_lake*slow_c14_lake  ) / f_shelf
     endif

     ! lake replacing vegetation and/or ice and/or shelf
     if( (.not. flag_done) .and. d_f_lake .gt. 0._wp .and. d_f_veg .le. 0._wp .and. d_f_ice .le. 0._wp .and. d_f_shelf .le. 0._wp ) then
       flag_done = .true.
      ! take carbon from vegetated part and shelf and lake and add it to ice carbon
      litter_c_lake   = (f_lake_old*litter_c_lake   - d_f_veg*litter_c   - d_f_ice*litter_c_ice   - d_f_shelf*litter_c_shelf  ) / f_lake 
      fast_c_lake     = (f_lake_old*fast_c_lake     - d_f_veg*fast_c     - d_f_ice*fast_c_ice     - d_f_shelf*fast_c_shelf    ) / f_lake
      slow_c_lake     = (f_lake_old*slow_c_lake     - d_f_veg*slow_c     - d_f_ice*slow_c_ice     - d_f_shelf*slow_c_shelf    ) / f_lake
      litter_c13_lake = (f_lake_old*litter_c13_lake - d_f_veg*litter_c13 - d_f_ice*litter_c13_ice - d_f_shelf*litter_c13_shelf) / f_lake
      fast_c13_lake   = (f_lake_old*fast_c13_lake   - d_f_veg*fast_c13   - d_f_ice*fast_c13_ice   - d_f_shelf*fast_c13_shelf  ) / f_lake
      slow_c13_lake   = (f_lake_old*slow_c13_lake   - d_f_veg*slow_c13   - d_f_ice*slow_c13_ice   - d_f_shelf*slow_c13_shelf  ) / f_lake
      litter_c14_lake = (f_lake_old*litter_c14_lake - d_f_veg*litter_c14 - d_f_ice*litter_c14_ice - d_f_shelf*litter_c14_shelf) / f_lake
      fast_c14_lake   = (f_lake_old*fast_c14_lake   - d_f_veg*fast_c14   - d_f_ice*fast_c14_ice   - d_f_shelf*fast_c14_shelf  ) / f_lake
      slow_c14_lake   = (f_lake_old*slow_c14_lake   - d_f_veg*slow_c14   - d_f_ice*slow_c14_ice   - d_f_shelf*slow_c14_shelf  ) / f_lake
     endif

     ! ice and/or shelf and/or lake expanding into vegetated area
     if( (.not. flag_done) .and. d_f_veg .lt. 0._wp .and. d_f_shelf .ge. 0._wp .and. d_f_ice .ge. 0._wp .and. d_f_lake .ge. 0._wp) then
       flag_done = .true.
      ! take carbon from the vegetation part and bury it under the ice sheet
      if( d_f_ice .ne. 0._wp ) then
       litter_c_ice   = (f_ice_grd_old*litter_c_ice   + d_f_ice*litter_c  ) / f_ice_grd
       fast_c_ice     = (f_ice_grd_old*fast_c_ice     + d_f_ice*fast_c    ) / f_ice_grd
       slow_c_ice     = (f_ice_grd_old*slow_c_ice     + d_f_ice*slow_c    ) / f_ice_grd
       litter_c13_ice = (f_ice_grd_old*litter_c13_ice + d_f_ice*litter_c13) / f_ice_grd
       fast_c13_ice   = (f_ice_grd_old*fast_c13_ice   + d_f_ice*fast_c13  ) / f_ice_grd
       slow_c13_ice   = (f_ice_grd_old*slow_c13_ice   + d_f_ice*slow_c13  ) / f_ice_grd
       litter_c14_ice = (f_ice_grd_old*litter_c14_ice + d_f_ice*litter_c14) / f_ice_grd
       fast_c14_ice   = (f_ice_grd_old*fast_c14_ice   + d_f_ice*fast_c14  ) / f_ice_grd
       slow_c14_ice   = (f_ice_grd_old*slow_c14_ice   + d_f_ice*slow_c14  ) / f_ice_grd
      endif
      ! take carbon from the vegetation part and put it into shelf carbon pool
      if( d_f_shelf .ne. 0._wp ) then
       litter_c_shelf   = (f_shelf_old*litter_c_shelf   + d_f_shelf*litter_c  ) / f_shelf
       fast_c_shelf     = (f_shelf_old*fast_c_shelf     + d_f_shelf*fast_c    ) / f_shelf
       slow_c_shelf     = (f_shelf_old*slow_c_shelf     + d_f_shelf*slow_c    ) / f_shelf
       litter_c13_shelf = (f_shelf_old*litter_c13_shelf + d_f_shelf*litter_c13) / f_shelf
       fast_c13_shelf   = (f_shelf_old*fast_c13_shelf   + d_f_shelf*fast_c13  ) / f_shelf
       slow_c13_shelf   = (f_shelf_old*slow_c13_shelf   + d_f_shelf*slow_c13  ) / f_shelf
       litter_c14_shelf = (f_shelf_old*litter_c14_shelf + d_f_shelf*litter_c14) / f_shelf
       fast_c14_shelf   = (f_shelf_old*fast_c14_shelf   + d_f_shelf*fast_c14  ) / f_shelf
       slow_c14_shelf   = (f_shelf_old*slow_c14_shelf   + d_f_shelf*slow_c14  ) / f_shelf
      endif
      ! take carbon from the vegetation part and bury it under lake
      if( d_f_lake .ne. 0._wp ) then
       litter_c_lake   = (f_lake_old*litter_c_lake   + d_f_lake*litter_c  ) / f_lake
       fast_c_lake     = (f_lake_old*fast_c_lake     + d_f_lake*fast_c    ) / f_lake
       slow_c_lake     = (f_lake_old*slow_c_lake     + d_f_lake*slow_c    ) / f_lake
       litter_c13_lake = (f_lake_old*litter_c13_lake + d_f_lake*litter_c13) / f_lake
       fast_c13_lake   = (f_lake_old*fast_c13_lake   + d_f_lake*fast_c13  ) / f_lake
       slow_c13_lake   = (f_lake_old*slow_c13_lake   + d_f_lake*slow_c13  ) / f_lake
       litter_c14_lake = (f_lake_old*litter_c14_lake + d_f_lake*litter_c14) / f_lake
       fast_c14_lake   = (f_lake_old*fast_c14_lake   + d_f_lake*fast_c14  ) / f_lake
       slow_c14_lake   = (f_lake_old*slow_c14_lake   + d_f_lake*slow_c14  ) / f_lake
      endif
     endif
 
     ! vegetation and/or shelf and/or lake expanding into ice area
     if( (.not. flag_done) .and. d_f_ice .lt. 0._wp .and. d_f_veg .ge. 0._wp .and. d_f_shelf .ge. 0._wp .and. d_f_lake .ge. 0._wp) then
       flag_done = .true.
      ! take ice carbon and add it to the vegetation part
      if( d_f_veg .ne. 0._wp ) then
       litter_c   = ((f_veg_old-f_peat)*litter_c   + d_f_veg*litter_c_ice  ) / (f_veg-f_peat) 
       fast_c     = ((f_veg_old-f_peat)*fast_c     + d_f_veg*fast_c_ice    ) / (f_veg-f_peat)
       slow_c     = ((f_veg_old-f_peat)*slow_c     + d_f_veg*slow_c_ice    ) / (f_veg-f_peat)
       litter_c13 = ((f_veg_old-f_peat)*litter_c13 + d_f_veg*litter_c13_ice) / (f_veg-f_peat)
       fast_c13   = ((f_veg_old-f_peat)*fast_c13   + d_f_veg*fast_c13_ice  ) / (f_veg-f_peat)
       slow_c13   = ((f_veg_old-f_peat)*slow_c13   + d_f_veg*slow_c13_ice  ) / (f_veg-f_peat)
       litter_c14 = ((f_veg_old-f_peat)*litter_c14 + d_f_veg*litter_c14_ice) / (f_veg-f_peat)
       fast_c14   = ((f_veg_old-f_peat)*fast_c14   + d_f_veg*fast_c14_ice  ) / (f_veg-f_peat)
       slow_c14   = ((f_veg_old-f_peat)*slow_c14   + d_f_veg*slow_c14_ice  ) / (f_veg-f_peat)
      endif
      ! take ice carbon and put it into shelf carbon pool
      if( d_f_shelf .ne. 0._wp ) then
       litter_c_shelf   = (f_shelf_old*litter_c_shelf   + d_f_shelf*litter_c_ice  ) / f_shelf
       fast_c_shelf     = (f_shelf_old*fast_c_shelf     + d_f_shelf*fast_c_ice    ) / f_shelf
       slow_c_shelf     = (f_shelf_old*slow_c_shelf     + d_f_shelf*slow_c_ice    ) / f_shelf
       litter_c13_shelf = (f_shelf_old*litter_c13_shelf + d_f_shelf*litter_c13_ice) / f_shelf
       fast_c13_shelf   = (f_shelf_old*fast_c13_shelf   + d_f_shelf*fast_c13_ice  ) / f_shelf
       slow_c13_shelf   = (f_shelf_old*slow_c13_shelf   + d_f_shelf*slow_c13_ice  ) / f_shelf
       litter_c14_shelf = (f_shelf_old*litter_c14_shelf + d_f_shelf*litter_c14_ice) / f_shelf
       fast_c14_shelf   = (f_shelf_old*fast_c14_shelf   + d_f_shelf*fast_c14_ice  ) / f_shelf
       slow_c14_shelf   = (f_shelf_old*slow_c14_shelf   + d_f_shelf*slow_c14_ice  ) / f_shelf
      endif
      ! take ice carbon and put it into lake carbon pool
      if( d_f_lake .ne. 0._wp ) then
       litter_c_lake   = (f_lake_old*litter_c_lake   + d_f_lake*litter_c_ice  ) / f_lake
       fast_c_lake     = (f_lake_old*fast_c_lake     + d_f_lake*fast_c_ice    ) / f_lake
       slow_c_lake     = (f_lake_old*slow_c_lake     + d_f_lake*slow_c_ice    ) / f_lake
       litter_c13_lake = (f_lake_old*litter_c13_lake + d_f_lake*litter_c13_ice) / f_lake
       fast_c13_lake   = (f_lake_old*fast_c13_lake   + d_f_lake*fast_c13_ice  ) / f_lake
       slow_c13_lake   = (f_lake_old*slow_c13_lake   + d_f_lake*slow_c13_ice  ) / f_lake
       litter_c14_lake = (f_lake_old*litter_c14_lake + d_f_lake*litter_c14_ice) / f_lake
       fast_c14_lake   = (f_lake_old*fast_c14_lake   + d_f_lake*fast_c14_ice  ) / f_lake
       slow_c14_lake   = (f_lake_old*slow_c14_lake   + d_f_lake*slow_c14_ice  ) / f_lake
      endif
     endif
 
     ! vegetation and/or ice and/or lake expanding into shelf area
     if( (.not. flag_done) .and. d_f_shelf .lt. 0._wp .and. d_f_veg .ge. 0._wp .and. d_f_ice .ge. 0._wp .and. d_f_lake .ge. 0._wp) then
       flag_done = .true.
      ! take shelf carbon and add it to the vegetation part
      if( d_f_veg .ne. 0._wp ) then
       litter_c   = ((f_veg_old-f_peat)*litter_c   + d_f_veg*litter_c_shelf  ) / (f_veg-f_peat) 
       fast_c     = ((f_veg_old-f_peat)*fast_c     + d_f_veg*fast_c_shelf    ) / (f_veg-f_peat)
       slow_c     = ((f_veg_old-f_peat)*slow_c     + d_f_veg*slow_c_shelf    ) / (f_veg-f_peat)
       litter_c13 = ((f_veg_old-f_peat)*litter_c13 + d_f_veg*litter_c13_shelf) / (f_veg-f_peat)
       fast_c13   = ((f_veg_old-f_peat)*fast_c13   + d_f_veg*fast_c13_shelf  ) / (f_veg-f_peat)
       slow_c13   = ((f_veg_old-f_peat)*slow_c13   + d_f_veg*slow_c13_shelf  ) / (f_veg-f_peat)
       litter_c14 = ((f_veg_old-f_peat)*litter_c14 + d_f_veg*litter_c14_shelf) / (f_veg-f_peat)
       fast_c14   = ((f_veg_old-f_peat)*fast_c14   + d_f_veg*fast_c14_shelf  ) / (f_veg-f_peat)
       slow_c14   = ((f_veg_old-f_peat)*slow_c14   + d_f_veg*slow_c14_shelf  ) / (f_veg-f_peat)
      endif
      ! take shelf carbon and put it into ice carbon pool
      if( d_f_ice .ne. 0._wp ) then
       litter_c_ice   = (f_ice_grd_old*litter_c_ice   + d_f_ice*litter_c_shelf  ) / f_ice_grd
       fast_c_ice     = (f_ice_grd_old*fast_c_ice     + d_f_ice*fast_c_shelf    ) / f_ice_grd
       slow_c_ice     = (f_ice_grd_old*slow_c_ice     + d_f_ice*slow_c_shelf    ) / f_ice_grd
       litter_c13_ice = (f_ice_grd_old*litter_c13_ice + d_f_ice*litter_c13_shelf) / f_ice_grd
       fast_c13_ice   = (f_ice_grd_old*fast_c13_ice   + d_f_ice*fast_c13_shelf  ) / f_ice_grd
       slow_c13_ice   = (f_ice_grd_old*slow_c13_ice   + d_f_ice*slow_c13_shelf  ) / f_ice_grd
       litter_c14_ice = (f_ice_grd_old*litter_c14_ice + d_f_ice*litter_c14_shelf) / f_ice_grd
       fast_c14_ice   = (f_ice_grd_old*fast_c14_ice   + d_f_ice*fast_c14_shelf  ) / f_ice_grd
       slow_c14_ice   = (f_ice_grd_old*slow_c14_ice   + d_f_ice*slow_c14_shelf  ) / f_ice_grd
      endif
      ! take shelf carbon and put it into lake carbon pool
      if( d_f_lake .ne. 0._wp ) then
       litter_c_lake   = (f_lake_old*litter_c_lake   + d_f_lake*litter_c_shelf  ) / f_lake
       fast_c_lake     = (f_lake_old*fast_c_lake     + d_f_lake*fast_c_shelf    ) / f_lake
       slow_c_lake     = (f_lake_old*slow_c_lake     + d_f_lake*slow_c_shelf    ) / f_lake
       litter_c13_lake = (f_lake_old*litter_c13_lake + d_f_lake*litter_c13_shelf) / f_lake
       fast_c13_lake   = (f_lake_old*fast_c13_lake   + d_f_lake*fast_c13_shelf  ) / f_lake
       slow_c13_lake   = (f_lake_old*slow_c13_lake   + d_f_lake*slow_c13_shelf  ) / f_lake
       litter_c14_lake = (f_lake_old*litter_c14_lake + d_f_lake*litter_c14_shelf) / f_lake
       fast_c14_lake   = (f_lake_old*fast_c14_lake   + d_f_lake*fast_c14_shelf  ) / f_lake
       slow_c14_lake   = (f_lake_old*slow_c14_lake   + d_f_lake*slow_c14_shelf  ) / f_lake
      endif
     endif
 
     ! vegetation and/or ice and/or lake expanding into lake area
     if( (.not. flag_done) .and. d_f_lake .lt. 0._wp .and. d_f_veg .ge. 0._wp .and. d_f_ice .ge. 0._wp .and. d_f_shelf .ge. 0._wp) then
       flag_done = .true.
      ! take lake carbon and add it to the vegetation part
      if( d_f_veg .ne. 0._wp ) then
       litter_c   = ((f_veg_old-f_peat)*litter_c   + d_f_veg*litter_c_lake  ) / (f_veg-f_peat) 
       fast_c     = ((f_veg_old-f_peat)*fast_c     + d_f_veg*fast_c_lake    ) / (f_veg-f_peat)
       slow_c     = ((f_veg_old-f_peat)*slow_c     + d_f_veg*slow_c_lake    ) / (f_veg-f_peat)
       litter_c13 = ((f_veg_old-f_peat)*litter_c13 + d_f_veg*litter_c13_lake) / (f_veg-f_peat)
       fast_c13   = ((f_veg_old-f_peat)*fast_c13   + d_f_veg*fast_c13_lake  ) / (f_veg-f_peat)
       slow_c13   = ((f_veg_old-f_peat)*slow_c13   + d_f_veg*slow_c13_lake  ) / (f_veg-f_peat)
       litter_c14 = ((f_veg_old-f_peat)*litter_c14 + d_f_veg*litter_c14_lake) / (f_veg-f_peat)
       fast_c14   = ((f_veg_old-f_peat)*fast_c14   + d_f_veg*fast_c14_lake  ) / (f_veg-f_peat)
       slow_c14   = ((f_veg_old-f_peat)*slow_c14   + d_f_veg*slow_c14_lake  ) / (f_veg-f_peat)
      endif
      ! take lake carbon and put it into ice carbon pool
      if( d_f_ice .ne. 0._wp ) then
       litter_c_ice   = (f_ice_grd_old*litter_c_ice   + d_f_ice*litter_c_lake  ) / f_ice_grd
       fast_c_ice     = (f_ice_grd_old*fast_c_ice     + d_f_ice*fast_c_lake    ) / f_ice_grd
       slow_c_ice     = (f_ice_grd_old*slow_c_ice     + d_f_ice*slow_c_lake    ) / f_ice_grd
       litter_c13_ice = (f_ice_grd_old*litter_c13_ice + d_f_ice*litter_c13_lake) / f_ice_grd
       fast_c13_ice   = (f_ice_grd_old*fast_c13_ice   + d_f_ice*fast_c13_lake  ) / f_ice_grd
       slow_c13_ice   = (f_ice_grd_old*slow_c13_ice   + d_f_ice*slow_c13_lake  ) / f_ice_grd
       litter_c14_ice = (f_ice_grd_old*litter_c14_ice + d_f_ice*litter_c14_lake) / f_ice_grd
       fast_c14_ice   = (f_ice_grd_old*fast_c14_ice   + d_f_ice*fast_c14_lake  ) / f_ice_grd
       slow_c14_ice   = (f_ice_grd_old*slow_c14_ice   + d_f_ice*slow_c14_lake  ) / f_ice_grd
      endif
      ! take lake carbon and put it into shelf carbon pool
      if( d_f_shelf .ne. 0._wp ) then
       litter_c_shelf   = (f_shelf_old*litter_c_shelf   + d_f_shelf*litter_c_lake  ) / f_shelf
       fast_c_shelf     = (f_shelf_old*fast_c_shelf     + d_f_shelf*fast_c_lake    ) / f_shelf
       slow_c_shelf     = (f_shelf_old*slow_c_shelf     + d_f_shelf*slow_c_lake    ) / f_shelf
       litter_c13_shelf = (f_shelf_old*litter_c13_shelf + d_f_shelf*litter_c13_lake) / f_shelf
       fast_c13_shelf   = (f_shelf_old*fast_c13_shelf   + d_f_shelf*fast_c13_lake  ) / f_shelf
       slow_c13_shelf   = (f_shelf_old*slow_c13_shelf   + d_f_shelf*slow_c13_lake  ) / f_shelf
       litter_c14_shelf = (f_shelf_old*litter_c14_shelf + d_f_shelf*litter_c14_lake) / f_shelf
       fast_c14_shelf   = (f_shelf_old*fast_c14_shelf   + d_f_shelf*fast_c14_lake  ) / f_shelf
       slow_c14_shelf   = (f_shelf_old*slow_c14_shelf   + d_f_shelf*slow_c14_lake  ) / f_shelf
      endif
     endif
 
     ! vegetation and shelf replacing ice and/or lake 
     if( (.not. flag_done) .and. d_f_veg.gt.0._wp .and. d_f_shelf.gt.0._wp .and. d_f_ice.le.0._wp .and. d_f_lake.le.0._wp .and. (-d_f_ice-d_f_lake).gt.0._wp ) then 
       flag_done = .true.
       ! take carbon from ice and/or lake and add it to the vegetation part
       fac = d_f_veg/(-d_f_ice-d_f_lake)
       litter_c   = ((f_veg_old-f_peat)*litter_c   - d_f_ice*fac*litter_c_ice   - d_f_lake*fac*litter_c_lake  ) / (f_veg-f_peat)
       fast_c     = ((f_veg_old-f_peat)*fast_c     - d_f_ice*fac*fast_c_ice     - d_f_lake*fac*fast_c_lake    ) / (f_veg-f_peat)
       slow_c     = ((f_veg_old-f_peat)*slow_c     - d_f_ice*fac*slow_c_ice     - d_f_lake*fac*slow_c_lake    ) / (f_veg-f_peat)
       litter_c13 = ((f_veg_old-f_peat)*litter_c13 - d_f_ice*fac*litter_c13_ice - d_f_lake*fac*litter_c13_lake) / (f_veg-f_peat)
       fast_c13   = ((f_veg_old-f_peat)*fast_c13   - d_f_ice*fac*fast_c13_ice   - d_f_lake*fac*fast_c13_lake  ) / (f_veg-f_peat)
       slow_c13   = ((f_veg_old-f_peat)*slow_c13   - d_f_ice*fac*slow_c13_ice   - d_f_lake*fac*slow_c13_lake  ) / (f_veg-f_peat)
       litter_c14 = ((f_veg_old-f_peat)*litter_c14 - d_f_ice*fac*litter_c14_ice - d_f_lake*fac*litter_c14_lake) / (f_veg-f_peat)
       fast_c14   = ((f_veg_old-f_peat)*fast_c14   - d_f_ice*fac*fast_c14_ice   - d_f_lake*fac*fast_c14_lake  ) / (f_veg-f_peat)
       slow_c14   = ((f_veg_old-f_peat)*slow_c14   - d_f_ice*fac*slow_c14_ice   - d_f_lake*fac*slow_c14_lake  ) / (f_veg-f_peat)
       ! take carbon from ice and/or lake and add it to the shelf part
       fac = d_f_shelf/(-d_f_ice-d_f_lake)
       litter_c_shelf   = (f_shelf_old*litter_c_shelf   - d_f_ice*fac*litter_c_ice   - d_f_lake*fac*litter_c_lake  ) / f_shelf 
       fast_c_shelf     = (f_shelf_old*fast_c_shelf     - d_f_ice*fac*fast_c_ice     - d_f_lake*fac*fast_c_lake    ) / f_shelf
       slow_c_shelf     = (f_shelf_old*slow_c_shelf     - d_f_ice*fac*slow_c_ice     - d_f_lake*fac*slow_c_lake    ) / f_shelf
       litter_c13_shelf = (f_shelf_old*litter_c13_shelf - d_f_ice*fac*litter_c13_ice - d_f_lake*fac*litter_c13_lake) / f_shelf
       fast_c13_shelf   = (f_shelf_old*fast_c13_shelf   - d_f_ice*fac*fast_c13_ice   - d_f_lake*fac*fast_c13_lake  ) / f_shelf
       slow_c13_shelf   = (f_shelf_old*slow_c13_shelf   - d_f_ice*fac*slow_c13_ice   - d_f_lake*fac*slow_c13_lake  ) / f_shelf
       litter_c14_shelf = (f_shelf_old*litter_c14_shelf - d_f_ice*fac*litter_c14_ice - d_f_lake*fac*litter_c14_lake) / f_shelf
       fast_c14_shelf   = (f_shelf_old*fast_c14_shelf   - d_f_ice*fac*fast_c14_ice   - d_f_lake*fac*fast_c14_lake  ) / f_shelf
       slow_c14_shelf   = (f_shelf_old*slow_c14_shelf   - d_f_ice*fac*slow_c14_ice   - d_f_lake*fac*slow_c14_lake  ) / f_shelf
     endif

     ! vegetation and ice replacing shelf and/or lake 
     if( (.not. flag_done) .and. d_f_veg.gt.0._wp .and. d_f_ice.gt.0._wp .and. d_f_shelf.le.0._wp .and. d_f_lake.le.0._wp .and. (-d_f_shelf-d_f_lake).gt.0._wp ) then 
       flag_done = .true.
       ! take carbon from shelf and/or lake and add it to the vegetation part
       fac = d_f_veg/(-d_f_shelf-d_f_lake)
       litter_c   = ((f_veg_old-f_peat)*litter_c   - d_f_shelf*fac*litter_c_shelf   - d_f_lake*fac*litter_c_lake  ) / (f_veg-f_peat)
       fast_c     = ((f_veg_old-f_peat)*fast_c     - d_f_shelf*fac*fast_c_shelf     - d_f_lake*fac*fast_c_lake    ) / (f_veg-f_peat)
       slow_c     = ((f_veg_old-f_peat)*slow_c     - d_f_shelf*fac*slow_c_shelf     - d_f_lake*fac*slow_c_lake    ) / (f_veg-f_peat)
       litter_c13 = ((f_veg_old-f_peat)*litter_c13 - d_f_shelf*fac*litter_c13_shelf - d_f_lake*fac*litter_c13_lake) / (f_veg-f_peat)
       fast_c13   = ((f_veg_old-f_peat)*fast_c13   - d_f_shelf*fac*fast_c13_shelf   - d_f_lake*fac*fast_c13_lake  ) / (f_veg-f_peat)
       slow_c13   = ((f_veg_old-f_peat)*slow_c13   - d_f_shelf*fac*slow_c13_shelf   - d_f_lake*fac*slow_c13_lake  ) / (f_veg-f_peat)
       litter_c14 = ((f_veg_old-f_peat)*litter_c14 - d_f_shelf*fac*litter_c14_shelf - d_f_lake*fac*litter_c14_lake) / (f_veg-f_peat)
       fast_c14   = ((f_veg_old-f_peat)*fast_c14   - d_f_shelf*fac*fast_c14_shelf   - d_f_lake*fac*fast_c14_lake  ) / (f_veg-f_peat)
       slow_c14   = ((f_veg_old-f_peat)*slow_c14   - d_f_shelf*fac*slow_c14_shelf   - d_f_lake*fac*slow_c14_lake  ) / (f_veg-f_peat)
       ! take carbon from shelf and/or lake and bury it under ice
       fac = d_f_ice/(-d_f_shelf-d_f_lake)
       litter_c_ice   = (f_ice_grd_old*litter_c_ice   - d_f_shelf*fac*litter_c_shelf   - d_f_lake*fac*litter_c_lake  ) / f_ice_grd 
       fast_c_ice     = (f_ice_grd_old*fast_c_ice     - d_f_shelf*fac*fast_c_shelf     - d_f_lake*fac*fast_c_lake    ) / f_ice_grd
       slow_c_ice     = (f_ice_grd_old*slow_c_ice     - d_f_shelf*fac*slow_c_shelf     - d_f_lake*fac*slow_c_lake    ) / f_ice_grd
       litter_c13_ice = (f_ice_grd_old*litter_c13_ice - d_f_shelf*fac*litter_c13_shelf - d_f_lake*fac*litter_c13_lake) / f_ice_grd
       fast_c13_ice   = (f_ice_grd_old*fast_c13_ice   - d_f_shelf*fac*fast_c13_shelf   - d_f_lake*fac*fast_c13_lake  ) / f_ice_grd
       slow_c13_ice   = (f_ice_grd_old*slow_c13_ice   - d_f_shelf*fac*slow_c13_shelf   - d_f_lake*fac*slow_c13_lake  ) / f_ice_grd
       litter_c14_ice = (f_ice_grd_old*litter_c14_ice - d_f_shelf*fac*litter_c14_shelf - d_f_lake*fac*litter_c14_lake) / f_ice_grd
       fast_c14_ice   = (f_ice_grd_old*fast_c14_ice   - d_f_shelf*fac*fast_c14_shelf   - d_f_lake*fac*fast_c14_lake  ) / f_ice_grd
       slow_c14_ice   = (f_ice_grd_old*slow_c14_ice   - d_f_shelf*fac*slow_c14_shelf   - d_f_lake*fac*slow_c14_lake  ) / f_ice_grd
     endif

     ! vegetation and lake replacing ice and/or shelf
     if( (.not. flag_done) .and. d_f_veg.gt.0._wp .and. d_f_lake.gt.0._wp .and. d_f_ice.le.0._wp .and. d_f_shelf.le.0._wp .and. (-d_f_ice-d_f_shelf).gt.0._wp ) then 
       flag_done = .true.
       ! take carbon from ice and/or shelf and add it to the vegetation part
       fac = d_f_veg/(-d_f_ice-d_f_shelf)
       litter_c   = ((f_veg_old-f_peat)*litter_c   - d_f_ice*fac*litter_c_ice   - d_f_shelf*fac*litter_c_shelf  ) / (f_veg-f_peat)
       fast_c     = ((f_veg_old-f_peat)*fast_c     - d_f_ice*fac*fast_c_ice     - d_f_shelf*fac*fast_c_shelf    ) / (f_veg-f_peat)
       slow_c     = ((f_veg_old-f_peat)*slow_c     - d_f_ice*fac*slow_c_ice     - d_f_shelf*fac*slow_c_shelf    ) / (f_veg-f_peat)
       litter_c13 = ((f_veg_old-f_peat)*litter_c13 - d_f_ice*fac*litter_c13_ice - d_f_shelf*fac*litter_c13_shelf) / (f_veg-f_peat)
       fast_c13   = ((f_veg_old-f_peat)*fast_c13   - d_f_ice*fac*fast_c13_ice   - d_f_shelf*fac*fast_c13_shelf  ) / (f_veg-f_peat)
       slow_c13   = ((f_veg_old-f_peat)*slow_c13   - d_f_ice*fac*slow_c13_ice   - d_f_shelf*fac*slow_c13_shelf  ) / (f_veg-f_peat)
       litter_c14 = ((f_veg_old-f_peat)*litter_c14 - d_f_ice*fac*litter_c14_ice - d_f_shelf*fac*litter_c14_shelf) / (f_veg-f_peat)
       fast_c14   = ((f_veg_old-f_peat)*fast_c14   - d_f_ice*fac*fast_c14_ice   - d_f_shelf*fac*fast_c14_shelf  ) / (f_veg-f_peat)
       slow_c14   = ((f_veg_old-f_peat)*slow_c14   - d_f_ice*fac*slow_c14_ice   - d_f_shelf*fac*slow_c14_shelf  ) / (f_veg-f_peat)
       ! take carbon from ice and/or shelf and add it to the lake part
       fac = d_f_lake/(-d_f_ice-d_f_shelf)
       litter_c_lake   = (f_lake_old*litter_c_lake   - d_f_ice*fac*litter_c_ice   - d_f_shelf*fac*litter_c_shelf  ) / f_lake 
       fast_c_lake     = (f_lake_old*fast_c_lake     - d_f_ice*fac*fast_c_ice     - d_f_shelf*fac*fast_c_shelf    ) / f_lake
       slow_c_lake     = (f_lake_old*slow_c_lake     - d_f_ice*fac*slow_c_ice     - d_f_shelf*fac*slow_c_shelf    ) / f_lake
       litter_c13_lake = (f_lake_old*litter_c13_lake - d_f_ice*fac*litter_c13_ice - d_f_shelf*fac*litter_c13_shelf) / f_lake
       fast_c13_lake   = (f_lake_old*fast_c13_lake   - d_f_ice*fac*fast_c13_ice   - d_f_shelf*fac*fast_c13_shelf  ) / f_lake
       slow_c13_lake   = (f_lake_old*slow_c13_lake   - d_f_ice*fac*slow_c13_ice   - d_f_shelf*fac*slow_c13_shelf  ) / f_lake
       litter_c14_lake = (f_lake_old*litter_c14_lake - d_f_ice*fac*litter_c14_ice - d_f_shelf*fac*litter_c14_shelf) / f_lake
       fast_c14_lake   = (f_lake_old*fast_c14_lake   - d_f_ice*fac*fast_c14_ice   - d_f_shelf*fac*fast_c14_shelf  ) / f_lake
       slow_c14_lake   = (f_lake_old*slow_c14_lake   - d_f_ice*fac*slow_c14_ice   - d_f_shelf*fac*slow_c14_shelf  ) / f_lake
     endif

     ! ice and shelf replacing vegetation and/or lake 
     if( (.not. flag_done) .and. d_f_ice.gt.0._wp .and. d_f_shelf.gt.0._wp .and. d_f_veg.le.0._wp .and. d_f_lake.le.0._wp .and. (-d_f_veg-d_f_lake).gt.0._wp ) then 
       flag_done = .true.
       ! take carbon from vegetation and/or lake and add it to the ice part
       fac = d_f_ice/(-d_f_veg-d_f_lake)
       litter_c_ice   = (f_ice_grd_old*litter_c_ice   - d_f_veg*fac*litter_c   - d_f_lake*fac*litter_c_lake  ) / f_ice_grd 
       fast_c_ice     = (f_ice_grd_old*fast_c_ice     - d_f_veg*fac*fast_c     - d_f_lake*fac*fast_c_lake    ) / f_ice_grd
       slow_c_ice     = (f_ice_grd_old*slow_c_ice     - d_f_veg*fac*slow_c     - d_f_lake*fac*slow_c_lake    ) / f_ice_grd
       litter_c13_ice = (f_ice_grd_old*litter_c13_ice - d_f_veg*fac*litter_c13 - d_f_lake*fac*litter_c13_lake) / f_ice_grd
       fast_c13_ice   = (f_ice_grd_old*fast_c13_ice   - d_f_veg*fac*fast_c13   - d_f_lake*fac*fast_c13_lake  ) / f_ice_grd
       slow_c13_ice   = (f_ice_grd_old*slow_c13_ice   - d_f_veg*fac*slow_c13   - d_f_lake*fac*slow_c13_lake  ) / f_ice_grd
       litter_c14_ice = (f_ice_grd_old*litter_c14_ice - d_f_veg*fac*litter_c14 - d_f_lake*fac*litter_c14_lake) / f_ice_grd
       fast_c14_ice   = (f_ice_grd_old*fast_c14_ice   - d_f_veg*fac*fast_c14   - d_f_lake*fac*fast_c14_lake  ) / f_ice_grd
       slow_c14_ice   = (f_ice_grd_old*slow_c14_ice   - d_f_veg*fac*slow_c14   - d_f_lake*fac*slow_c14_lake  ) / f_ice_grd
       ! take carbon from vegetation and/or lake and add it to the shelf part
       fac = d_f_shelf/(-d_f_veg-d_f_lake)
       litter_c_shelf   = (f_shelf_old*litter_c_shelf   - d_f_veg*fac*litter_c   - d_f_lake*fac*litter_c_lake  ) / f_shelf 
       fast_c_shelf     = (f_shelf_old*fast_c_shelf     - d_f_veg*fac*fast_c     - d_f_lake*fac*fast_c_lake    ) / f_shelf
       slow_c_shelf     = (f_shelf_old*slow_c_shelf     - d_f_veg*fac*slow_c     - d_f_lake*fac*slow_c_lake    ) / f_shelf
       litter_c13_shelf = (f_shelf_old*litter_c13_shelf - d_f_veg*fac*litter_c13 - d_f_lake*fac*litter_c13_lake) / f_shelf
       fast_c13_shelf   = (f_shelf_old*fast_c13_shelf   - d_f_veg*fac*fast_c13   - d_f_lake*fac*fast_c13_lake  ) / f_shelf
       slow_c13_shelf   = (f_shelf_old*slow_c13_shelf   - d_f_veg*fac*slow_c13   - d_f_lake*fac*slow_c13_lake  ) / f_shelf
       litter_c14_shelf = (f_shelf_old*litter_c14_shelf - d_f_veg*fac*litter_c14 - d_f_lake*fac*litter_c14_lake) / f_shelf
       fast_c14_shelf   = (f_shelf_old*fast_c14_shelf   - d_f_veg*fac*fast_c14   - d_f_lake*fac*fast_c14_lake  ) / f_shelf
       slow_c14_shelf   = (f_shelf_old*slow_c14_shelf   - d_f_veg*fac*slow_c14   - d_f_lake*fac*slow_c14_lake  ) / f_shelf
     endif

     ! lake and shelf replacing vegetation and/or ice 
     if( (.not. flag_done) .and. d_f_lake.gt.0._wp .and. d_f_shelf.gt.0._wp .and. d_f_veg.le.0._wp .and. d_f_ice.le.0._wp .and. (-d_f_veg-d_f_ice).gt.0._wp ) then 
       flag_done = .true.
       ! take carbon from vegetation and/or ice and add it to the lake part
       fac = d_f_lake/(-d_f_veg-d_f_ice)
       litter_c_lake   = (f_lake_old*litter_c_lake   - d_f_veg*fac*litter_c   - d_f_ice*fac*litter_c_ice   ) / f_lake 
       fast_c_lake     = (f_lake_old*fast_c_lake     - d_f_veg*fac*fast_c     - d_f_ice*fac*fast_c_ice     ) / f_lake
       slow_c_lake     = (f_lake_old*slow_c_lake     - d_f_veg*fac*slow_c     - d_f_ice*fac*slow_c_ice     ) / f_lake
       litter_c13_lake = (f_lake_old*litter_c13_lake - d_f_veg*fac*litter_c13 - d_f_ice*fac*litter_c13_ice ) / f_lake
       fast_c13_lake   = (f_lake_old*fast_c13_lake   - d_f_veg*fac*fast_c13   - d_f_ice*fac*fast_c13_ice   ) / f_lake
       slow_c13_lake   = (f_lake_old*slow_c13_lake   - d_f_veg*fac*slow_c13   - d_f_ice*fac*slow_c13_ice   ) / f_lake
       litter_c14_lake = (f_lake_old*litter_c14_lake - d_f_veg*fac*litter_c14 - d_f_ice*fac*litter_c14_ice ) / f_lake
       fast_c14_lake   = (f_lake_old*fast_c14_lake   - d_f_veg*fac*fast_c14   - d_f_ice*fac*fast_c14_ice   ) / f_lake
       slow_c14_lake   = (f_lake_old*slow_c14_lake   - d_f_veg*fac*slow_c14   - d_f_ice*fac*slow_c14_ice   ) / f_lake
       ! take carbon from vegetation and/or ice and add it to the shelf part
       fac = d_f_shelf/(-d_f_veg-d_f_ice)
       litter_c_shelf   = (f_shelf_old*litter_c_shelf   - d_f_veg*fac*litter_c   - d_f_ice*fac*litter_c_ice  ) / f_shelf 
       fast_c_shelf     = (f_shelf_old*fast_c_shelf     - d_f_veg*fac*fast_c     - d_f_ice*fac*fast_c_ice    ) / f_shelf
       slow_c_shelf     = (f_shelf_old*slow_c_shelf     - d_f_veg*fac*slow_c     - d_f_ice*fac*slow_c_ice    ) / f_shelf
       litter_c13_shelf = (f_shelf_old*litter_c13_shelf - d_f_veg*fac*litter_c13 - d_f_ice*fac*litter_c13_ice) / f_shelf
       fast_c13_shelf   = (f_shelf_old*fast_c13_shelf   - d_f_veg*fac*fast_c13   - d_f_ice*fac*fast_c13_ice  ) / f_shelf
       slow_c13_shelf   = (f_shelf_old*slow_c13_shelf   - d_f_veg*fac*slow_c13   - d_f_ice*fac*slow_c13_ice  ) / f_shelf
       litter_c14_shelf = (f_shelf_old*litter_c14_shelf - d_f_veg*fac*litter_c14 - d_f_ice*fac*litter_c14_ice) / f_shelf
       fast_c14_shelf   = (f_shelf_old*fast_c14_shelf   - d_f_veg*fac*fast_c14   - d_f_ice*fac*fast_c14_ice  ) / f_shelf
       slow_c14_shelf   = (f_shelf_old*slow_c14_shelf   - d_f_veg*fac*slow_c14   - d_f_ice*fac*slow_c14_ice  ) / f_shelf
     endif

     ! lake and ice replacing vegetation and/or shelf
     if( (.not. flag_done) .and. d_f_lake.gt.0._wp .and. d_f_ice.gt.0._wp .and. d_f_veg.le.0._wp .and. d_f_shelf.le.0._wp .and. (-d_f_veg-d_f_shelf).gt.0._wp ) then 
       flag_done = .true.
       ! take carbon from vegetation and/or shelf and add it to the lake part
       fac = d_f_lake/(-d_f_veg-d_f_shelf)
       litter_c_lake   = (f_lake_old*litter_c_lake   - d_f_veg*fac*litter_c   - d_f_shelf*fac*litter_c_shelf  ) / f_lake 
       fast_c_lake     = (f_lake_old*fast_c_lake     - d_f_veg*fac*fast_c     - d_f_shelf*fac*fast_c_shelf    ) / f_lake
       slow_c_lake     = (f_lake_old*slow_c_lake     - d_f_veg*fac*slow_c     - d_f_shelf*fac*slow_c_shelf    ) / f_lake
       litter_c13_lake = (f_lake_old*litter_c13_lake - d_f_veg*fac*litter_c13 - d_f_shelf*fac*litter_c13_shelf) / f_lake
       fast_c13_lake   = (f_lake_old*fast_c13_lake   - d_f_veg*fac*fast_c13   - d_f_shelf*fac*fast_c13_shelf  ) / f_lake
       slow_c13_lake   = (f_lake_old*slow_c13_lake   - d_f_veg*fac*slow_c13   - d_f_shelf*fac*slow_c13_shelf  ) / f_lake
       litter_c14_lake = (f_lake_old*litter_c14_lake - d_f_veg*fac*litter_c14 - d_f_shelf*fac*litter_c14_shelf) / f_lake
       fast_c14_lake   = (f_lake_old*fast_c14_lake   - d_f_veg*fac*fast_c14   - d_f_shelf*fac*fast_c14_shelf  ) / f_lake
       slow_c14_lake   = (f_lake_old*slow_c14_lake   - d_f_veg*fac*slow_c14   - d_f_shelf*fac*slow_c14_shelf  ) / f_lake
       ! take carbon from vegetation and/or ice and add it to the ice part
       fac = d_f_ice/(-d_f_veg-d_f_shelf)
       litter_c_ice   = (f_ice_grd_old*litter_c_ice   - d_f_veg*fac*litter_c   - d_f_shelf*fac*litter_c_shelf  ) / f_ice_grd 
       fast_c_ice     = (f_ice_grd_old*fast_c_ice     - d_f_veg*fac*fast_c     - d_f_shelf*fac*fast_c_shelf    ) / f_ice_grd
       slow_c_ice     = (f_ice_grd_old*slow_c_ice     - d_f_veg*fac*slow_c     - d_f_shelf*fac*slow_c_shelf    ) / f_ice_grd
       litter_c13_ice = (f_ice_grd_old*litter_c13_ice - d_f_veg*fac*litter_c13 - d_f_shelf*fac*litter_c13_shelf) / f_ice_grd
       fast_c13_ice   = (f_ice_grd_old*fast_c13_ice   - d_f_veg*fac*fast_c13   - d_f_shelf*fac*fast_c13_shelf  ) / f_ice_grd
       slow_c13_ice   = (f_ice_grd_old*slow_c13_ice   - d_f_veg*fac*slow_c13   - d_f_shelf*fac*slow_c13_shelf  ) / f_ice_grd
       litter_c14_ice = (f_ice_grd_old*litter_c14_ice - d_f_veg*fac*litter_c14 - d_f_shelf*fac*litter_c14_shelf) / f_ice_grd
       fast_c14_ice   = (f_ice_grd_old*fast_c14_ice   - d_f_veg*fac*fast_c14   - d_f_shelf*fac*fast_c14_shelf  ) / f_ice_grd
       slow_c14_ice   = (f_ice_grd_old*slow_c14_ice   - d_f_veg*fac*slow_c14   - d_f_shelf*fac*slow_c14_shelf  ) / f_ice_grd
     endif

     carbon_tot_old = (f_veg_old - f_peat_old) * (sum(litter_c_old*dz_c(1:nlc)) + sum(fast_c_old*dz_c(1:nlc)) + sum(slow_c_old*dz_c(1:nlc))) &
                    + f_peat_old * (litter_c_peat_old+acro_c_old+sum(cato_c_old*dz_c(1:nlc))) &
                    + f_shelf_old * (sum(litter_c_shelf_old*dz_c(1:nlc)) + sum(fast_c_shelf_old*dz_c(1:nlc)) + sum(slow_c_shelf_old*dz_c(1:nlc))) &
                    + f_lake_old * (sum(litter_c_lake_old*dz_c(1:nlc)) + sum(fast_c_lake_old*dz_c(1:nlc)) + sum(slow_c_lake_old*dz_c(1:nlc))) &
                    + f_ice_grd_old * (sum(litter_c_ice_old*dz_c(1:nlc)) + sum(fast_c_ice_old*dz_c(1:nlc)) + sum(slow_c_ice_old*dz_c(1:nlc))) 
     carbon13_tot_old = (f_veg_old - f_peat_old) * (sum(litter_c13_old*dz_c(1:nlc)) + sum(fast_c13_old*dz_c(1:nlc)) + sum(slow_c13_old*dz_c(1:nlc))) &
                    + f_peat_old * (litter_c13_peat_old+acro_c13_old+sum(cato_c13_old*dz_c(1:nlc))) &
                    + f_shelf_old * (sum(litter_c13_shelf_old*dz_c(1:nlc)) + sum(fast_c13_shelf_old*dz_c(1:nlc)) + sum(slow_c13_shelf_old*dz_c(1:nlc))) &
                    + f_lake_old * (sum(litter_c13_lake_old*dz_c(1:nlc)) + sum(fast_c13_lake_old*dz_c(1:nlc)) + sum(slow_c13_lake_old*dz_c(1:nlc))) &
                    + f_ice_grd_old * (sum(litter_c13_ice_old*dz_c(1:nlc)) + sum(fast_c13_ice_old*dz_c(1:nlc)) + sum(slow_c13_ice_old*dz_c(1:nlc)))  
     carbon14_tot_old = (f_veg_old - f_peat_old) * (sum(litter_c14_old*dz_c(1:nlc)) + sum(fast_c14_old*dz_c(1:nlc)) + sum(slow_c14_old*dz_c(1:nlc))) &
                    + f_peat_old * (litter_c14_peat_old+acro_c14_old+sum(cato_c14_old*dz_c(1:nlc))) &
                    + f_shelf_old * (sum(litter_c14_shelf_old*dz_c(1:nlc)) + sum(fast_c14_shelf_old*dz_c(1:nlc)) + sum(slow_c14_shelf_old*dz_c(1:nlc))) &
                    + f_lake_old * (sum(litter_c14_lake_old*dz_c(1:nlc)) + sum(fast_c14_lake_old*dz_c(1:nlc)) + sum(slow_c14_lake_old*dz_c(1:nlc))) &
                    + f_ice_grd_old * (sum(litter_c14_ice_old*dz_c(1:nlc)) + sum(fast_c14_ice_old*dz_c(1:nlc)) + sum(slow_c14_ice_old*dz_c(1:nlc)))  

     carbon_tot = (f_veg - f_peat) * (sum(litter_c*dz_c(1:nlc)) + sum(fast_c*dz_c(1:nlc)) + sum(slow_c*dz_c(1:nlc))) &
                + f_peat * (litter_c_peat+acro_c+sum(cato_c*dz_c(1:nlc))) &
                + f_shelf * (sum(litter_c_shelf*dz_c(1:nlc)) + sum(fast_c_shelf*dz_c(1:nlc)) + sum(slow_c_shelf*dz_c(1:nlc))) &
                + f_lake * (sum(litter_c_lake*dz_c(1:nlc)) + sum(fast_c_lake*dz_c(1:nlc)) + sum(slow_c_lake*dz_c(1:nlc))) &
                + f_ice_grd * (sum(litter_c_ice*dz_c(1:nlc)) + sum(fast_c_ice*dz_c(1:nlc)) + sum(slow_c_ice*dz_c(1:nlc)))         
     carbon13_tot = (f_veg - f_peat) * (sum(litter_c13*dz_c(1:nlc)) + sum(fast_c13*dz_c(1:nlc)) + sum(slow_c13*dz_c(1:nlc))) &
                  + f_peat * (litter_c13_peat+acro_c13+sum(cato_c13*dz_c(1:nlc))) &
                  + f_shelf * (sum(litter_c13_shelf*dz_c(1:nlc)) + sum(fast_c13_shelf*dz_c(1:nlc)) + sum(slow_c13_shelf*dz_c(1:nlc))) &
                  + f_lake * (sum(litter_c13_lake*dz_c(1:nlc)) + sum(fast_c13_lake*dz_c(1:nlc)) + sum(slow_c13_lake*dz_c(1:nlc))) &
                  + f_ice_grd * (sum(litter_c13_ice*dz_c(1:nlc)) + sum(fast_c13_ice*dz_c(1:nlc)) + sum(slow_c13_ice*dz_c(1:nlc)))
     carbon14_tot = (f_veg - f_peat) * (sum(litter_c14*dz_c(1:nlc)) + sum(fast_c14*dz_c(1:nlc)) + sum(slow_c14*dz_c(1:nlc))) &
                  + f_peat * (litter_c14_peat+acro_c14+sum(cato_c14*dz_c(1:nlc))) &
                  + f_shelf * (sum(litter_c14_shelf*dz_c(1:nlc)) + sum(fast_c14_shelf*dz_c(1:nlc)) + sum(slow_c14_shelf*dz_c(1:nlc))) &
                  + f_lake * (sum(litter_c14_lake*dz_c(1:nlc)) + sum(fast_c14_lake*dz_c(1:nlc)) + sum(slow_c14_lake*dz_c(1:nlc))) &
                  + f_ice_grd * (sum(litter_c14_ice*dz_c(1:nlc)) + sum(fast_c14_ice*dz_c(1:nlc)) + sum(slow_c14_ice*dz_c(1:nlc)))

     if( abs(carbon_tot-carbon_tot_old) .gt. 1.e-5_wp ) then
      print *,''
      print *,'carbon conservation after carbon transfer ',carbon_tot-carbon_tot_old
      print *,'carbon before ',carbon_tot_old
      print *,'carbon after ',carbon_tot
      print *,'dveg,dice,dlake,dshelf,dpeat',d_f_veg,d_f_ice,d_f_lake,d_f_shelf,d_f_peat
      print *,'veg,veg_old',f_veg,f_veg_old
      print *,'ice,ice_old',f_ice_grd,f_ice_grd_old
      print *,'shelf,shelf_ild',f_shelf,f_shelf_old
      print *,'lake,lake_old',f_lake,f_lake_old
      print *,'peat,peat_old',f_peat,f_peat_old
      !stop
     endif
     if( abs(carbon13_tot-carbon13_tot_old) .gt. 1.e-5_wp ) print *,'carbon 13 conservation after carbon transfer ',carbon13_tot-carbon13_tot_old
     if( abs(carbon14_tot-carbon14_tot_old) .gt. 1.e-5_wp ) print *,'carbon 14 conservation after carbon transfer ',carbon14_tot-carbon14_tot_old


     return

  end subroutine carbon_trans

end module carbon_trans_mod
