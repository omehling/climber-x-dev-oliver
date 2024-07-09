!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : i c e _ m o d 
!
!  Purpose : Ice-sheet related quantities for SMB
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
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
module ice_mod

  use precision, only : wp
  use timer, only : n_year_smb, sec_year
  use smb_params, only : surf_par

  implicit none

  private
  public :: frac_ice, albedo_ice, margin_ice

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  f r a c _ i c e
  !   Purpose    :  compute ice fraction for background albedo
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine frac_ice(mask_ice, mask_margin, h_ice, z_sur_std, f_ice)

  implicit none

  integer, intent(in) :: mask_ice(:,:)
  integer, intent(in) :: mask_margin(:,:)
  real(wp), intent(in) :: h_ice(:,:)
  real(wp), intent(in) :: z_sur_std(:,:)

  real(wp), intent(out) :: f_ice(:,:)

  integer :: i, j, ii, jj, ni, nj, n
  real(wp) :: hice
  real(wp), dimension(:,:), allocatable :: h_ice_tmp

  real(wp), parameter :: eps = 1.e-10_wp


  ! subgrid ice fraction
  if (surf_par%i_f_ice.eq.0) then

    ! ice everywhere
    f_ice(:,:) = 1._wp

  else if (surf_par%i_f_ice.eq.1) then

    ! ice only over actual ice points
    where (mask_ice.eq.1) 
      f_ice = 1._wp
    elsewhere
      f_ice = 0._wp
    endwhere

  else if (surf_par%i_f_ice.eq.2) then

    ni = size(mask_ice,1)
    nj = size(mask_ice,2)

    ! compute ice thickness at the margin using neighboring points
    allocate(h_ice_tmp(ni,nj))
    h_ice_tmp = h_ice

    do j=2,nj-1
      do i=2,ni-1
        if (mask_margin(i,j).eq.1) then
          n = 0
          hice = 0._wp
          do jj=j-1,j+1
            do ii=i-1,i+1
              n = n+1
              hice = hice + h_ice_tmp(ii,jj)
            enddo
          enddo
          h_ice_tmp(i,j) = hice/real(n,wp)
        endif
      enddo
    enddo

    ! compute ice cover fraction from ice thickness and topographic roughness
    f_ice(:,:) = min(1._wp,h_ice_tmp(:,:)/(surf_par%h_ice_crit+eps)) &
               * tanh(h_ice_tmp(:,:)/(surf_par%c_fice*max(0._wp,z_sur_std(:,:)-surf_par%z_sur_std_crit)+eps))

    deallocate(h_ice_tmp)

  endif

  return

  end subroutine frac_ice


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  a l b e d o _ i c e
  !   Purpose    :  compute albedo of bare ice sheet
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine albedo_ice(mask_ice, mask_margin, f_ice, f_ice_old, ann_icemelt, ann_smb, dt_snowfree, &
                        alb_ice)

    implicit none

    integer, intent(in) :: mask_ice(:,:)
    integer, intent(in) :: mask_margin(:,:)
    real(wp), intent(in) :: f_ice(:,:)
    real(wp), intent(in) :: f_ice_old(:,:)
    real(wp), intent(in) :: ann_icemelt(:,:)
    real(wp), intent(in) :: ann_smb(:,:)
    real(wp), intent(in) :: dt_snowfree(:,:)

    real(wp), intent(inout) :: alb_ice(:,:)

    integer :: i, j, n, ii, jj, ni, nj
    real(wp) :: delta_t
    real(wp) :: tau_alb_firn_clean
    real(wp) :: aice

    real(wp), parameter :: eps = 1.e-10_wp


    ! ice albedo
    if (surf_par%i_alb_ice.eq.0) then

      ! constant ice albedo
      alb_ice(:,:) = surf_par%alb_ice_const

    else if (surf_par%i_alb_ice.eq.1) then

      ! variable ice albedo accounting for transition from firn to ice and for
      ! ice becoming dirty when exposed to the atmosphere 

      ni = size(f_ice,1)
      nj = size(f_ice,2)

      do j=1,nj
        do i=1,ni

          if (mask_margin(i,j).eq.0) then
            ! for non-marginal points

            if (f_ice_old(i,j).eq.0._wp .and. f_ice(i,j).eq.0._wp) then
              ! no ice, dummy value
              alb_ice(i,j) = surf_par%alb_firn
            else if (f_ice_old(i,j).eq.0._wp .and. f_ice(i,j).gt.0._wp) then
              ! newly formed ice
              alb_ice(i,j) = surf_par%alb_firn
            else if (f_ice(i,j).gt.0._wp) then
              delta_t = sec_year*real(n_year_smb,wp)  ! s, accounting for frequency of SEMI calls
              ! decrease albedo towards the value for clean ice when firn is melting
              if (alb_ice(i,j).gt.surf_par%alb_ice_clean .and. ann_icemelt(i,j).gt.eps) then
                tau_alb_firn_clean = max(1._wp,surf_par%w_firn/ann_icemelt(i,j))*sec_year*real(n_year_smb,wp)  ! s, time it takes to melt firn layer
                ! explicit
                !alb_ice(i,j) = alb_ice(i,j) - (alb_ice(i,j)-surf_par%alb_ice_clean)/tau_alb_firn_clean * delta_t
                ! implicit
                alb_ice(i,j) = (alb_ice(i,j)/delta_t+surf_par%alb_ice_clean/tau_alb_firn_clean) / (1._wp/delta_t+1._wp/tau_alb_firn_clean)
              endif
              ! reset to firn albedo if positive mass balance 
              if (ann_smb(i,j).gt.eps) then
                alb_ice(i,j) = surf_par%alb_firn 
              endif
              ! decrease albedo towards the value for dirty ice when bare ice is exposed, account for frequency of SEMI calls
              delta_t = max(1._wp,dt_snowfree(i,j)*real(n_year_smb,wp))  ! s, time snowfree, accounting for frequency of SEMI calls
              ! explicit
              !alb_ice(i,j) = alb_ice(i,j) - (alb_ice(i,j)-surf_par%alb_ice_dirty)/surf_par%tau_alb_ice_dirty * delta_t
              ! implicit
              alb_ice(i,j) = (alb_ice(i,j)/delta_t+surf_par%alb_ice_dirty/surf_par%tau_alb_ice_dirty) / (1._wp/delta_t+1._wp/surf_par%tau_alb_ice_dirty)
            endif

          endif

        enddo
      enddo

      if (surf_par%i_alb_ice_margin.eq.0) then
        ! constant ice margin albedo 

        where (mask_margin.eq.1) 
          alb_ice = surf_par%alb_ice_margin_const
        endwhere

      else if (surf_par%i_alb_ice_margin.eq.1) then
        ! for marginal points, compute ice albedo as average over ice cell neighbors

        do j=2,nj-1
          do i=2,ni-1
            if (mask_margin(i,j).eq.1) then
              alb_ice(i,j) = 0._wp
              n = 0
              aice = 0._wp
              do jj=j-1,j+1
                do ii=i-1,i+1
                  if (mask_ice(ii,jj).eq.1) then
                    n = n+1
                    aice = aice + alb_ice(ii,jj)
                  endif
                enddo
              enddo
              if (n.gt.0) then
                alb_ice(i,j) = aice/real(n,wp)
              endif
            endif
          enddo
        enddo

      endif

    endif

    return

  end subroutine albedo_ice


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  m a r g i n _ i c e
  !   Purpose    :  derive ice margin mask
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine margin_ice(mask_ice, mask_margin)

  implicit none

  integer, intent(in) :: mask_ice(:,:)
  integer, intent(out) :: mask_margin(:,:)

  integer :: i, j, ii, jj, ni, nj


  ni = size(mask_ice,1)
  nj = size(mask_ice,2)
  do j=1,nj
    do i=1,ni
      mask_margin(i,j) = 0
      ! look for ice neighbors
      if (mask_ice(i,j).eq.0 .and. j>1 .and. j<nj .and. i>1 .and. i<ni) then
        do jj=j-1,j+1
          do ii=i-1,i+1
            if (mask_ice(ii,jj).eq.1) then
              mask_margin(i,j) = 1
            endif
          enddo
        enddo
      endif
    enddo
  enddo

  return

  end subroutine margin_ice

end module ice_mod
