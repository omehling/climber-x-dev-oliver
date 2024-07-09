!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ b a s _ m e l t _ m
!
!> @file
!!
!! Computation of the basal melting rate.
!!
!! @section Copyright
!!
!! Copyright 2009-2017 Ralf Greve, Ben Galton-Fenzi, Tatsuru Sato
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
!> Computation of the basal melting rate.
!<------------------------------------------------------------------------------

module calc_bas_melt_m

  use timer, only : sec_year, sec_year_inv
  use sico_types_m
  use sico_state
  use sico_grid_mod
  use sico_timer
  use sico_params, only : sico_par_class, DELTA_TM_SW, rhow_rho_ratio, RHO_SW
  use ice_material_properties_m, only : kappa_val


  implicit none

  private
  public :: calc_qbm

contains

!-------------------------------------------------------------------------------
!> Computation of the basal melting rate Q_bm.
!! Summation of Q_bm and Q_tld (water drainage rate from the temperate layer).
!<------------------------------------------------------------------------------
subroutine calc_qbm(st,grd,tmr,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_timer_class), intent(in) :: tmr
type(sico_par_class), intent(in) :: par

integer :: i, j
integer :: n_ocean
logical :: flag_float
real(wp) :: frictional_heating
real(wp) :: Q_bm_grounded, Q_bm_marine
real(wp) :: qbm_min, qbm_max


!-------- Computation of Q_bm --------

st%Q_bm = 0.0_wp   ! initialisation

do i=1, grd%IMAX-1
do j=1, grd%JMAX-1

  if (par%margin==1 .or. par%margin==2) then

    ! flag for 'floating ice'
    flag_float = st%H_c(j,i)+st%H_t(j,i) < (st%z_sl(j,i)-st%zl(j,i))*RHO_SW/RHO

    if (st%maske(j,i)==0 .and. .not.flag_float) then   ! grounded ice

      if (st%n_cts(j,i)==-1) then

        frictional_heating = 0.0_wp
        st%Q_bm(j,i)       = 0.0_wp

      else if (st%n_cts(j,i)==0) then

        if (par%dynamics==2) then

          if (.not.st%flag_shelfy_stream(j,i)) then

            frictional_heating &
              = -grd%aqbm3a*st%H_c(j,i) &
              *0.5_wp*(st%vx_t(0,j,i)+st%vx_t(0,j,i-1)) &
              *st%dzs_dxi_g(j,i) &
              -grd%aqbm3a*st%H_c(j,i) &
              *0.5_wp*(st%vy_t(0,j,i)+st%vy_t(0,j-1,i)) &
              *st%dzs_deta_g(j,i)

          else   ! flag_shelfy_stream(j,i) == .true.

            if (par%hyb_mode==1) then
              frictional_heating &
                = grd%aqbm3b &
                * st%beta_drag(j,i) &
                * (st%vx_b_g(j,i)**2 + st%vy_b_g(j,i)**2)
            else
              frictional_heating &
                = grd%aqbm3b &
                * st%c_drag(j,i) &
                * sqrt(st%vx_b_g(j,i)**2  &
                +st%vy_b_g(j,i)**2) &
                **(1.0_wp+st%p_weert_inv(j,i))
            endif

          end if

        else    ! par%dynamics.ne.2

          frictional_heating &
            = -grd%aqbm3a*st%H_c(j,i) &
            *0.5_wp*(st%vx_t(0,j,i)+st%vx_t(0,j,i-1)) &
            *st%dzs_dxi_g(j,i) &
            -grd%aqbm3a*st%H_c(j,i) &
            *0.5_wp*(st%vy_t(0,j,i)+st%vy_t(0,j-1,i)) &
            *st%dzs_deta_g(j,i)

        endif

        st%Q_bm(j,i) =   grd%aqbm1*(st%temp_c(1,j,i)-st%temp_c(0,j,i))/st%H_c(j,i) &
          *kappa_val(st%temp_c(0,j,i)) &
          - grd%aqbm2*(st%temp_r(grd%KRMAX,j,i)-st%temp_r(grd%KRMAX-1,j,i)) &
          + frictional_heating

      else   ! n_cts(j,i)==1

        if (par%dynamics==2) then

          if (.not.st%flag_shelfy_stream(j,i)) then

            frictional_heating &
              = -grd%aqbm3a*(st%H_c(j,i)+st%H_t(j,i)) &
              *0.5_wp*(st%vx_t(0,j,i)+st%vx_t(0,j,i-1)) &
              *st%dzs_dxi_g(j,i) &
              -grd%aqbm3a*(st%H_c(j,i)+st%H_t(j,i)) &
              *0.5_wp*(st%vy_t(0,j,i)+st%vy_t(0,j-1,i)) &
              *st%dzs_deta_g(j,i)

          else   ! flag_shelfy_stream(j,i) == .true.

            if (par%hyb_mode==1) then
              frictional_heating &
                = grd%aqbm3b &
                * st%beta_drag(j,i) &
                * (st%vx_b_g(j,i)**2 + st%vy_b_g(j,i)**2)
            else
              frictional_heating &
                = grd%aqbm3b &
                * st%c_drag(j,i) &
                * sqrt(st%vx_b_g(j,i)**2  &
                +st%vy_b_g(j,i)**2) &
                **(1.0_wp+st%p_weert_inv(j,i))
            endif

          end if

        else ! par%dynamics.ne.2

          frictional_heating &
            = -grd%aqbm3a*(st%H_c(j,i)+st%H_t(j,i)) &
            *0.5_wp*(st%vx_t(0,j,i)+st%vx_t(0,j,i-1)) &
            *st%dzs_dxi_g(j,i) &
            -grd%aqbm3a*(st%H_c(j,i)+st%H_t(j,i)) &
            *0.5_wp*(st%vy_t(0,j,i)+st%vy_t(0,j-1,i)) &
            *st%dzs_deta_g(j,i)

        endif

        st%Q_bm(j,i) =   grd%aqbm4*kappa_val(st%temp_t_m(0,j,i)) &
          - grd%aqbm2*(st%temp_r(grd%KRMAX,j,i)-st%temp_r(grd%KRMAX-1,j,i)) &
          + frictional_heating

      end if

      if (par%marine_ice_basal_melting==1) then

        !!continue   ! do nothing

      else if (par%marine_ice_basal_melting==2) then

        if ( (st%zb(j,i) < st%z_sl(j,i)) &          ! marine ice
          .and. &
          (     (st%maske(j,i+1)==2) &   ! at least one
          .or.(st%maske(j,i-1)==2) &   ! nearest neighbour
          .or.(st%maske(j+1,i)==2) &   ! is
          .or.(st%maske(j-1,i)==2) &   ! ocean
          ) &
          ) then

          st%Q_bm(j,i) = par%qbm_marine *sec_year_inv*rhow_rho_ratio
          ! m/a water equiv. -> m/s ice equiv.

        end if

      else if (par%marine_ice_basal_melting==3) then

        if ( (st%zb(j,i) < st%z_sl(j,i)) &          ! marine ice
          .and. &
          (     (st%maske(j,i+1)==2) &   ! at least one
          .or.(st%maske(j,i-1)==2) &   ! nearest neighbour
          .or.(st%maske(j+1,i)==2) &   ! is
          .or.(st%maske(j-1,i)==2) &   ! ocean
          ) &
          ) then

          n_ocean = 0
          if (st%maske(j,i+1)==2) n_ocean = n_ocean+1
          if (st%maske(j,i-1)==2) n_ocean = n_ocean+1
          if (st%maske(j+1,i)==2) n_ocean = n_ocean+1
          if (st%maske(j-1,i)==2) n_ocean = n_ocean+1

          if ( n_ocean > 0 ) then

            Q_bm_grounded = st%Q_bm(j,i)
            Q_bm_marine   = par%qbm_marine *sec_year_inv*rhow_rho_ratio
            ! m/a water equiv. -> m/s ice equiv.

            st%Q_bm(j,i) = (1.0_wp-0.25_wp*real(n_ocean,dp)) * Q_bm_grounded &
              +0.25_wp*real(n_ocean,dp)  * Q_bm_marine
            ! weighed average of grounded ice melting (computed)
            ! and marine ice melting (prescribed)
          else
            write(6,'(a)') ' >>> calc_qbm: Marine ice margin point does not'
            write(6,'(a)') '               have an ocean neighbour!'
            stop
          end if

        end if

      endif

    else if ((st%maske(j,i)==0 .and. flag_float) .or. st%maske(j,i)==2) then   ! floating ice

      ! MW, apply basal melt to ice that would be floating according to floatation criterium
      st%Q_bm(j,i) = st%Q_bm_float(j,i) ! m/s ice equiv. 

    end if


  else if (par%margin==3) then

    if (st%maske(j,i)==0) then   ! grounded ice

      if (st%n_cts(j,i)==-1) then

        frictional_heating = 0.0_wp
        st%Q_bm(j,i)       = 0.0_wp

      else if (st%n_cts(j,i)==0) then

        if (par%dynamics==2) then

          if (.not.st%flag_shelfy_stream(j,i)) then

            frictional_heating &
              = -grd%aqbm3a*st%H_c(j,i) &
              *0.5_wp*(st%vx_t(0,j,i)+st%vx_t(0,j,i-1)) &
              *st%dzs_dxi_g(j,i) &
              -grd%aqbm3a*st%H_c(j,i) &
              *0.5_wp*(st%vy_t(0,j,i)+st%vy_t(0,j-1,i)) &
              *st%dzs_deta_g(j,i)

          else   ! flag_shelfy_stream(j,i) == .true.

            if (par%hyb_mode==1) then
              frictional_heating &
                = grd%aqbm3b &
                * st%beta_drag(j,i) &
                * (st%vx_b_g(j,i)**2 + st%vy_b_g(j,i)**2)
            else
              frictional_heating &
                = grd%aqbm3b &
                * st%c_drag(j,i) &
                * sqrt(st%vx_b_g(j,i)**2  &
                +st%vy_b_g(j,i)**2) &
                **(1.0_wp+st%p_weert_inv(j,i))
            endif

          end if

        else    ! par%dynamics.ne.2

          frictional_heating &
            = -grd%aqbm3a*st%H_c(j,i) &
            *0.5_wp*(st%vx_t(0,j,i)+st%vx_t(0,j,i-1)) &
            *st%dzs_dxi_g(j,i) &
            -grd%aqbm3a*st%H_c(j,i) &
            *0.5_wp*(st%vy_t(0,j,i)+st%vy_t(0,j-1,i)) &
            *st%dzs_deta_g(j,i)

        endif

        st%Q_bm(j,i) =   grd%aqbm1*(st%temp_c(1,j,i)-st%temp_c(0,j,i))/st%H_c(j,i) &
          *kappa_val(st%temp_c(0,j,i)) &
          - grd%aqbm2*(st%temp_r(grd%KRMAX,j,i)-st%temp_r(grd%KRMAX-1,j,i)) &
          + frictional_heating

      else   ! n_cts(j,i)==1

        if (par%dynamics==2) then

          if (.not.st%flag_shelfy_stream(j,i)) then

            frictional_heating &
              = -grd%aqbm3a*(st%H_c(j,i)+st%H_t(j,i)) &
              *0.5_wp*(st%vx_t(0,j,i)+st%vx_t(0,j,i-1)) &
              *st%dzs_dxi_g(j,i) &
              -grd%aqbm3a*(st%H_c(j,i)+st%H_t(j,i)) &
              *0.5_wp*(st%vy_t(0,j,i)+st%vy_t(0,j-1,i)) &
              *st%dzs_deta_g(j,i)

          else   ! flag_shelfy_stream(j,i) == .true.

            if (par%hyb_mode==1) then
              frictional_heating &
                = grd%aqbm3b &
                * st%beta_drag(j,i) &
                * (st%vx_b_g(j,i)**2 + st%vy_b_g(j,i)**2)
            else
              frictional_heating &
                = grd%aqbm3b &
                * st%c_drag(j,i) &
                * sqrt(st%vx_b_g(j,i)**2  &
                +st%vy_b_g(j,i)**2) &
                **(1.0_wp+st%p_weert_inv(j,i))
            endif

          end if

        else ! par%dynamics.ne.2

          frictional_heating &
            = -grd%aqbm3a*(st%H_c(j,i)+st%H_t(j,i)) &
            *0.5_wp*(st%vx_t(0,j,i)+st%vx_t(0,j,i-1)) &
            *st%dzs_dxi_g(j,i) &
            -grd%aqbm3a*(st%H_c(j,i)+st%H_t(j,i)) &
            *0.5_wp*(st%vy_t(0,j,i)+st%vy_t(0,j-1,i)) &
            *st%dzs_deta_g(j,i)

        endif

        st%Q_bm(j,i) =   grd%aqbm4*kappa_val(st%temp_t_m(0,j,i)) &
          - grd%aqbm2*(st%temp_r(grd%KRMAX,j,i)-st%temp_r(grd%KRMAX-1,j,i)) &
          + frictional_heating

      end if

      if (st%flag_grounding_line_1(j,i)) then
                                        ! grounding line (grounded-ice side)
        !!continue   ! do nothing

      else if ( (st%zb(j,i) < st%z_sl(j,i)) &          ! marine ice margin
        .and. &
        (     (st%maske(j,i+1)>=2) &   !  (at least one
        .or.(st%maske(j,i-1)>=2) &   !   nearest neighbour
        .or.(st%maske(j+1,i)>=2) &   !   is
        .or.(st%maske(j-1,i)>=2) &   !   ocean)
        ) &
        ) then

        if (par%marine_ice_basal_melting==1) then

          !!continue   ! do nothing

        else if (par%marine_ice_basal_melting==2) then

          st%Q_bm(j,i) = par%qbm_marine *sec_year_inv*rhow_rho_ratio
          ! m/a water equiv. -> m/s ice equiv.

        else if (par%marine_ice_basal_melting==3) then

          n_ocean = 0
          if (st%maske(j,i+1)>=2) n_ocean = n_ocean+1
          if (st%maske(j,i-1)>=2) n_ocean = n_ocean+1
          if (st%maske(j+1,i)>=2) n_ocean = n_ocean+1
          if (st%maske(j-1,i)>=2) n_ocean = n_ocean+1

          if ( n_ocean > 0 ) then

            Q_bm_grounded = st%Q_bm(j,i)
            Q_bm_marine   = par%qbm_marine *sec_year_inv*rhow_rho_ratio
            ! m/a water equiv. -> m/s ice equiv.

            st%Q_bm(j,i) = (1.0_wp-0.25_wp*real(n_ocean,dp)) * Q_bm_grounded &
              +0.25_wp*real(n_ocean,dp)  * Q_bm_marine
            ! weighed average of grounded ice melting (computed)
            ! and marine ice melting (prescribed)
          else
            write(6,'(a)') ' >>> calc_qbm: Marine ice margin point does not'
            write(6,'(a)') '          have a floating ice or ocean neighbour!'
            stop
          end if

        endif

      end if

    else if ( (st%maske(j,i)==2).or.(st%maske(j,i)==3) ) then   ! floating ice or ocean

          st%Q_bm(j,i) = st%Q_bm_float(j,i) ! m/s ice equiv. 

    end if

  end if

end do
end do

!-------- Limitation of Q_bm, Q_tld and Q_b_tot
!                       (only for grounded ice) --------

qbm_min = par%qbm_min *sec_year_inv*rhow_rho_ratio
! m/a water equiv. -> m/s ice equiv.

qbm_max = par%qbm_max *sec_year_inv*rhow_rho_ratio
! m/a water equiv. -> m/s ice equiv.

do i=0, grd%IMAX
  do j=0, grd%JMAX
    if (st%Q_tld(j,i)   < qbm_min) st%Q_tld(j,i)   = 0.0_wp
    if (st%Q_tld(j,i)   > qbm_max) st%Q_tld(j,i)   = qbm_max
  end do
end do

!-------- Sum of Q_bm and Q_tld --------

st%Q_b_tot = st%Q_bm + st%Q_tld


if (par%marine_ice_basal_melting==2 .or. par%marine_ice_basal_melting==3) then

  if (par%qbm_marine*sec_year_inv*rhow_rho_ratio > qbm_max) then
    write (6,'(a)') ' >>> calc_qbm: QBM_MARINE'
    write (6,'(a)') '               (basal melting rate at the ice front)'
    write (6,'(a)') '               is larger than the limiter qbm_max!'
    stop
  end if

endif

do i=0, grd%IMAX
  do j=0, grd%JMAX
    if (st%maske(j,i)==0) then   ! grounded ice
      if (st%Q_bm(j,i)    < qbm_min) st%Q_bm(j,i)    = 0.0_wp
      if (st%Q_bm(j,i)    > qbm_max) st%Q_bm(j,i)    = qbm_max
      if (st%Q_tld(j,i)   < qbm_min) st%Q_tld(j,i)   = 0.0_wp
      if (st%Q_tld(j,i)   > qbm_max) st%Q_tld(j,i)   = qbm_max
      if (st%Q_b_tot(j,i) < qbm_min) st%Q_b_tot(j,i) = 0.0_wp
      if (st%Q_b_tot(j,i) > qbm_max) st%Q_b_tot(j,i) = qbm_max
    end if
  end do
end do

end subroutine calc_qbm

end module calc_bas_melt_m
!
