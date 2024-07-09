!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : i m o _ m o d e l
!
!  Purpose : main IMO model
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
module imo_model

  use nml
  use ncio
  use precision, only : wp, dp
  use constants, only : rho_i, rho_sw, Lf, cap_w
  use timer, only : nmon_year, nstep_year_imo, time_soy_imo, time_eoy_imo
  use control, only : out_dir, restart_in_dir, imo_restart, i_map
  use imo_grid, only : imo_grid_init
  use imo_params, only : imo_params_init, dt, i_imo, imo_const, k_1, k_2, l_fix_depth, fix_depth, depth_disc, l_depth_scale, zl_ref, l_bm_lake
  use imo_def, only : imo_class
  use coord, only : grid_allocate, grid_class
  use coord, only : map_init, map_class, map_field
  use coord, only : map_scrip_init, map_scrip_class, map_scrip_field

  implicit none


  private
  public :: imo_init, imo_update, imo_end, imo_write_restart

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i m o _ u p d a t e
  !   Purpose    :  update ice melt to ocean
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine imo_update(imo)

  !$  use omp_lib

  implicit none

  type(imo_class) :: imo

  integer :: i, j, im1, ip1, jm1, jp1, n, k, k1, k2
  integer :: nk_ocn, k_fix_ocn, k_disc_ocn
  integer :: nk_lake, k_fix_lake, k_disc_lake
  logical :: flag
  real(wp) :: t_sum, s_sum
  real(wp) :: t_npol, s_npol, t_spol, s_spol
  real(wp) :: w1, w2


  !---------------------------
  ! ocean

  nk_ocn = size(imo%z_ocn_in)

  ! extrapolate ocean temperature and salinity to cells where they are undefined
  do k = 1,nk_ocn
    flag = .true.
    imo%mask_ocn_in_tmp = imo%mask_ocn_in(:,:,k)
    do while (flag)
      flag = .false.
      do j = 1,imo%grid_in%G%ny
        do i = 1,imo%grid_in%G%nx
          if (imo%mask_ocn_in(i,j,k)==0) then
            n = 0
            t_sum = 0._wp
            s_sum = 0._wp
            ip1 = i+1; if (ip1.eq.imo%grid_in%G%nx+1) ip1 = 1
            im1 = i-1; if (im1.eq.0) im1 = imo%grid_in%G%nx
            jp1 = min(j+1,imo%grid_in%G%ny)
            jm1 = max(j-1,1)
            if (imo%mask_ocn_in(im1,j,k)==1) then
              n = n+1
              t_sum = t_sum + imo%t_ocn_in(im1,j,k)
              s_sum = s_sum + imo%s_ocn_in(im1,j,k)
            endif
            if (imo%mask_ocn_in(ip1,j,k)==1) then
              n = n+1
              t_sum = t_sum + imo%t_ocn_in(ip1,j,k)
              s_sum = s_sum + imo%s_ocn_in(ip1,j,k)
            endif
            if (imo%mask_ocn_in(i,jm1,k)==1) then
              n = n+1
              t_sum = t_sum + imo%t_ocn_in(i,jm1,k)
              s_sum = s_sum + imo%s_ocn_in(i,jm1,k)
            endif
            if (imo%mask_ocn_in(i,jp1,k)==1) then
              n = n+1
              t_sum = t_sum + imo%t_ocn_in(i,jp1,k)
              s_sum = s_sum + imo%s_ocn_in(i,jp1,k)
            endif
            if (n>0) then
              flag = .true.
              imo%t_ocn_in(i,j,k) = t_sum/real(n,wp)
              imo%s_ocn_in(i,j,k) = s_sum/real(n,wp)
              ! update ocean mask 
              imo%mask_ocn_in_tmp(i,j) = 1
            endif
          endif
        enddo
      enddo
      imo%mask_ocn_in(:,:,k) = imo%mask_ocn_in_tmp
    enddo
  enddo

  ! polar smoothing
  do k = 1,nk_ocn
    t_npol = 0._wp
    s_npol = 0._wp
    t_spol = 0._wp
    s_spol = 0._wp
    do i = 1,imo%grid_in%G%nx
      t_npol = t_npol + imo%t_ocn_in(i,imo%grid_in%G%ny,k) / imo%grid_in%G%nx
      s_npol = s_npol + imo%s_ocn_in(i,imo%grid_in%G%ny,k) / imo%grid_in%G%nx
      t_spol = t_spol + imo%t_ocn_in(i,1,k) / imo%grid_in%G%nx
      s_spol = s_spol + imo%s_ocn_in(i,1,k) / imo%grid_in%G%nx
    enddo
    imo%t_ocn_in(:,1,k) = t_spol
    imo%s_ocn_in(:,1,k) = s_spol
    imo%t_ocn_in(:,imo%grid_in%G%ny,k) = t_npol
    imo%s_ocn_in(:,imo%grid_in%G%ny,k) = s_npol
  enddo

  if (.not.allocated(imo%t_ocn)) allocate(imo%t_ocn(imo%grid%G%nx,imo%grid%G%ny,nk_ocn))
  if (.not.allocated(imo%s_ocn)) allocate(imo%s_ocn(imo%grid%G%nx,imo%grid%G%ny,nk_ocn))

  ! map shelf temperature and salinity to ice sheet grid
  if (i_map==1) then
    !$omp parallel do private(k)
    do k=1,nk_ocn
      call map_field(imo%map_cmn_to_ice,"t_ocn",imo%t_ocn_in(:,:,k), imo%t_ocn(:,:,k),method="bilinear")
      call map_field(imo%map_cmn_to_ice,"s_ocn",imo%s_ocn_in(:,:,k), imo%s_ocn(:,:,k),method="bilinear")
    enddo
    !$omp end parallel do 
  else if (i_map==2) then
    !$omp parallel do private(k)
    do k=1,nk_ocn
      call map_scrip_field(imo%maps_cmn_to_ice,"t_ocn",imo%t_ocn_in(:,:,k), imo%t_ocn(:,:,k),method="mean",missing_value=-9999._dp)
      call map_scrip_field(imo%maps_cmn_to_ice,"s_ocn",imo%s_ocn_in(:,:,k), imo%s_ocn(:,:,k),method="mean",missing_value=-9999._dp)
    enddo
    !$omp end parallel do 
  endif


  !---------------------------
  ! lake

  nk_lake = size(imo%z_lake_in)

  ! extrapolate lake temperature and salinity to cells where they are undefined
  flag = .true.
  imo%mask_lake_in_tmp = imo%mask_lake_in
  do while (flag)
    flag = .false.
    do j = 1,imo%grid_in%G%ny
      do i = 1,imo%grid_in%G%nx
        if (imo%mask_lake_in(i,j)==0) then
          ip1 = i+1; if (ip1.eq.imo%grid_in%G%nx+1) ip1 = 1
          im1 = i-1; if (im1.eq.0) im1 = imo%grid_in%G%nx
          jp1 = min(j+1,imo%grid_in%G%ny)
          jm1 = max(j-1,1)
          do k = 1,nk_lake
            n = 0
            t_sum = 0._wp
            s_sum = 0._wp
            if (imo%mask_lake_in(im1,j)==1) then
              n = n+1
              t_sum = t_sum + imo%t_lake_in(im1,j,k)
              s_sum = s_sum + imo%s_lake_in(im1,j,k)
            endif
            if (imo%mask_lake_in(ip1,j)==1) then
              n = n+1
              t_sum = t_sum + imo%t_lake_in(ip1,j,k)
              s_sum = s_sum + imo%s_lake_in(ip1,j,k)
            endif
            if (imo%mask_lake_in(i,jm1)==1) then
              n = n+1
              t_sum = t_sum + imo%t_lake_in(i,jm1,k)
              s_sum = s_sum + imo%s_lake_in(i,jm1,k)
            endif
            if (imo%mask_lake_in(i,jp1)==1) then
              n = n+1
              t_sum = t_sum + imo%t_lake_in(i,jp1,k)
              s_sum = s_sum + imo%s_lake_in(i,jp1,k)
            endif
            if (n>0) then
              imo%t_lake_in(i,j,k) = t_sum/real(n,wp)
              imo%s_lake_in(i,j,k) = s_sum/real(n,wp)
              ! update ocean mask 
              imo%mask_lake_in_tmp(i,j) = 1
            endif
          enddo
          flag = .true.
        endif
      enddo
    enddo
    imo%mask_lake_in = imo%mask_lake_in_tmp
  enddo

  if (.not.allocated(imo%t_lake)) allocate(imo%t_lake(imo%grid%G%nx,imo%grid%G%ny,nk_lake))
  if (.not.allocated(imo%s_lake)) allocate(imo%s_lake(imo%grid%G%nx,imo%grid%G%ny,nk_lake))

  ! map lake temperature and salinity to ice sheet grid
  if (i_map==1) then
    !$omp parallel do private(k)
    do k=1,nk_lake
      call map_field(imo%map_cmn_to_ice,"t_lake",imo%t_lake_in(:,:,k), imo%t_lake(:,:,k),method="bilinear")
      call map_field(imo%map_cmn_to_ice,"s_lake",imo%s_lake_in(:,:,k), imo%s_lake(:,:,k),method="bilinear")
    enddo
    !$omp end parallel do 
  else if (i_map==2) then
    !$omp parallel do private(k)
    do k=1,nk_lake
      call map_scrip_field(imo%maps_cmn_to_ice,"t_lake",imo%t_lake_in(:,:,k), imo%t_lake(:,:,k),method="mean",missing_value=-9999._dp)
      call map_scrip_field(imo%maps_cmn_to_ice,"s_lake",imo%s_lake_in(:,:,k), imo%s_lake(:,:,k),method="mean",missing_value=-9999._dp)
    enddo
    !$omp end parallel do 
  endif

  ! index of fix_depth 
  k_fix_ocn = minloc(abs(imo%z_ocn_in-fix_depth),1)
  k_disc_ocn = minloc(abs(imo%z_ocn_in-depth_disc),1)
  k_fix_lake = minloc(abs(imo%z_lake_in-fix_depth),1)
  k_disc_lake = minloc(abs(imo%z_lake_in-depth_disc),1)

  !$omp parallel do collapse(2) private(i,j,k,k1,k2,w1,w2)
  do i = 1,imo%grid%G%nx
    do j = 1,imo%grid%G%ny

      if (time_soy_imo) then
        imo%imo_ann(i,j) = 0._wp
        imo%t_disc(i,j) = 0._wp
        imo%s_disc(i,j) = 0._wp
      endif

      if (imo%mask_ocn_lake(i,j).eq.0) then
        ! use ocean values

        ! annual mean values for small-scale basal melt in ice sheet model
        imo%t_disc(i,j) = imo%t_disc(i,j) + imo%t_ocn(i,j,k_disc_ocn)/real(nstep_year_imo,wp)
        imo%s_disc(i,j) = imo%s_disc(i,j) + imo%s_ocn(i,j,k_disc_ocn)/real(nstep_year_imo,wp) 

        ! select which level to use
        if (l_fix_depth) then
          ! use temperature and salinity at fixed depth
          imo%t_imo(i,j) = imo%t_ocn(i,j,k_fix_ocn)
          imo%s_imo(i,j) = imo%s_ocn(i,j,k_fix_ocn) 
        else
          ! interpolate temperature and salinity to depth of ice shelf base
          k = 1
          do while (imo%z_ocn_in(k).lt.-imo%zb(i,j) .and. k.le.nk_ocn)
            k = k+1
          enddo
          k1 = max(1,k-1)
          k2 = min(nk_ocn,k)
          if (k1.eq.k2) then
            w1 = 1._wp
          else
            w1  = 1._wp - (-imo%zb(i,j)-imo%z_ocn_in(k1)) / (imo%z_ocn_in(k2)-imo%z_ocn_in(k1))
          endif
          w2  = 1._wp-w1
          !print *
          !print *,'zb',-imo%zb(i,j)
          !print *,'z_ocn_in',imo%z_ocn_in
          !print *,'k1,k2',k1,k2
          !print *,'w1,w2',w1,w2
          !print *,'t_ocn',imo%t_ocn(i,j,:)
          !print *,'s_ocn',imo%s_ocn(i,j,:)
          imo%t_imo(i,j) = w1*imo%t_ocn(i,j,k1) + w2*imo%t_ocn(i,j,k2)
          imo%s_imo(i,j) = w1*imo%s_ocn(i,j,k1) + w2*imo%s_ocn(i,j,k2)
          !print *,'t_imo',imo%t_imo(i,j)
          !print *,'s_imo',imo%s_imo(i,j)
        endif

      else if (imo%mask_ocn_lake(i,j).eq.1) then
        ! use lake values

        ! annual mean values for small-scale basal melt in ice sheet model
        imo%t_disc(i,j) = imo%t_disc(i,j) + imo%t_lake(i,j,k_disc_lake)/real(nstep_year_imo,wp)
        imo%s_disc(i,j) = imo%s_disc(i,j) + imo%s_lake(i,j,k_disc_lake)/real(nstep_year_imo,wp) 

        ! select which level to use
        if (l_fix_depth) then
          ! use temperature and salinity at fixed depth
          imo%t_imo(i,j) = imo%t_lake(i,j,k_fix_lake)
          imo%s_imo(i,j) = imo%s_lake(i,j,k_fix_lake) 
        else
          ! interpolate temperature and salinity to depth of ice shelf base
          k = 1
          do while (imo%z_lake_in(k).lt.-imo%zb(i,j) .and. k.le.nk_lake)
            k = k+1
          enddo
          k1 = max(1,k-1)
          k2 = min(nk_lake,k)
          if (k1.eq.k2) then
            w1 = 1._wp
          else
            w1  = 1._wp - (-imo%zb(i,j)-imo%z_lake_in(k1)) / (imo%z_lake_in(k2)-imo%z_lake_in(k1))
          endif
          w2  = 1._wp-w1
          !print *
          !print *,'zb',-imo%zb(i,j)
          !print *,'z_lake_in',imo%z_lake_in
          !print *,'k1,k2',k1,k2
          !print *,'w1,w2',w1,w2
          !print *,'t_lake',imo%t_lake(i,j,:)
          !print *,'s_lake',imo%s_lake(i,j,:)
          imo%t_imo(i,j) = w1*imo%t_lake(i,j,k1) + w2*imo%t_lake(i,j,k2)
          imo%s_imo(i,j) = w1*imo%s_lake(i,j,k1) + w2*imo%s_lake(i,j,k2)
          !print *,'t_imo',imo%t_imo(i,j)
          !print *,'s_imo',imo%s_imo(i,j)
        endif

      endif

      ! freezing temperature after Beckmann & Goose 2003 (eq. 2), zb<0
      imo%t_freeze(i,j) = 0.0939_wp - 0.057_wp*imo%s_imo(i,j) + 7.64e-4_wp*imo%zb(i,j) ! degC

      if (i_imo==0) then
        ! constant and uniform basal melt rate
        imo%imo(i,j) = imo_const  ! kg/m2/s
      else if (i_imo==1) then
        ! Beckmann & Goose 2003
        imo%imo(i,j) = k_1 * rho_sw*cap_w/(rho_i*Lf) * (imo%t_imo(i,j)-imo%t_freeze(i,j))  ! kg/m2/s
      else if (i_imo==2) then
        ! Pollard & DeConto 2012, eq 17
        imo%imo(i,j) = k_2 * rho_sw*cap_w/(rho_i*Lf) * abs(imo%t_imo(i,j)-imo%t_freeze(i,j))*(imo%t_imo(i,j)-imo%t_freeze(i,j))  ! kg/m2/s
      endif

      if (imo%mask_ocn_lake(i,j).eq.1 .and. .not.l_bm_lake) then
        ! suppress basal melt in lakes
        imo%imo(i,j) = 0._wp
      endif

      ! additional dependence on ocean depth, reduce melt in shallow water and enhance in deep water
      if (l_depth_scale) then
        if (imo%imo(i,j).gt.0._wp) then ! only for melt case, exclude refreezing
          imo%imo(i,j) = imo%imo(i,j) * (1._wp+max(0._wp,max(0._wp,-imo%zl_fil(i,j))-zl_ref)/zl_ref)
        endif
      endif

      ! exclude melt in Hudson Bay
!      where (imo%grid%lat.gt.50._wp .and. imo%grid%lat.lt.70._wp .and. imo%grid%lon.gt.-100._wp .and. imo%grid%lon.lt.-60._wp)
!        imo%imo = 0._wp
!      endwhere

      imo%imo_ann(i,j) = imo%imo_ann(i,j) + imo%imo(i,j)*dt   ! imo integrated over the year, kg/m2

    enddo
  enddo
  !$omp end parallel do


  return

  end subroutine imo_update
      

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i m o _ i n i t
  !   Purpose    :  initialize ice melt to ocean 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine imo_init(imo,grid,cmn_grid)

    implicit none

    type(imo_class) :: imo
    type(grid_class), intent(in) :: grid 
    type(grid_class), intent(in) :: cmn_grid 


    call imo_params_init
    call imo_grid_init

    ! store grid information
    imo%grid = grid

    ! Generate mapping from cmn to imo/ice
    if (i_map==1) then
      call map_init(imo%map_cmn_to_ice,cmn_grid,imo%grid,max_neighbors=4,lat_lim=10._dp,dist_max=1.e6_dp)
    else if (i_map==2) then
      call map_scrip_init(imo%maps_cmn_to_ice,cmn_grid,imo%grid,method="bil",fldr="maps",load=.TRUE.,clean=.FALSE.)
    endif

    imo%grid_in = cmn_grid

    ! Allocate arrays 
    call grid_allocate(imo%grid_in, imo%mask_ocn_in_tmp )
    call grid_allocate(imo%grid_in, imo%mask_lake_in_tmp )
    call grid_allocate(imo%grid, imo%mask_ocn_lake)
    call grid_allocate(imo%grid, imo%mask_ice_shelf )
    call grid_allocate(imo%grid, imo%zb )
    call grid_allocate(imo%grid, imo%t_imo )
    call grid_allocate(imo%grid, imo%s_imo )
    call grid_allocate(imo%grid, imo%t_disc )
    call grid_allocate(imo%grid, imo%s_disc )
    call grid_allocate(imo%grid, imo%t_freeze )
    call grid_allocate(imo%grid, imo%imo )
    call grid_allocate(imo%grid, imo%imo_ann )

    if (imo_restart) then

      ! read restart file 
      call imo_read_restart("restart/"//trim(restart_in_dir)//"/imo_"//trim(imo%grid%name)//"_restart.nc",imo)
      print *,'read restart file ',"restart/"//trim(restart_in_dir)//"/imo_"//trim(imo%grid%name)//"_restart.nc"

    else

      ! initialize for new run
      imo%imo     = 0._wp
      imo%imo_ann = 0._wp

    endif

    print*
    print*,'======================================================='
    print*,' Initialisation of imo ', trim(imo%grid%name),' complete'
    print*,'======================================================='
    print*

  return

  end subroutine imo_init
      
     
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  i m o _ e n d 
  ! Purpose  :  end imo
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine imo_end(imo)

    implicit none

    type(imo_class) :: imo


    ! Deallocate all state variables to free memory
    call imo_dealloc(imo)


    return

  end subroutine imo_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  i m o _ d e a l l o c 
  ! Purpose  :  deallocate imo variables
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine imo_dealloc(imo)

    implicit none 

    type(imo_class) :: imo


    ! Deallocate all state variables to free memory
    deallocate(imo%mask_ocn_in)
    deallocate(imo%mask_ocn_in_tmp)
    deallocate(imo%t_ocn_in )
    deallocate(imo%s_ocn_in )
    deallocate(imo%mask_lake_in)
    deallocate(imo%mask_lake_in_tmp)
    deallocate(imo%t_lake_in )
    deallocate(imo%s_lake_in )
    deallocate(imo%t_ocn )
    deallocate(imo%s_ocn )
    deallocate(imo%t_lake )
    deallocate(imo%s_lake )
    deallocate(imo%mask_ocn_lake)
    deallocate(imo%mask_ice_shelf )
    deallocate(imo%zb )
    deallocate(imo%t_imo )
    deallocate(imo%s_imo )
    deallocate(imo%t_disc )
    deallocate(imo%s_disc )
    deallocate(imo%t_freeze )
    deallocate(imo%imo )
    deallocate(imo%imo_ann )

    return

  end subroutine imo_dealloc 


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  i m o _ w r i t e _ r e s t a r t
  ! Purpose  :  Write restart file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine imo_write_restart(fnm,imo)

    use dim_name, only: dim_x, dim_y

    implicit none

    type(imo_class) :: imo

    character (len=*) :: fnm
    integer :: ncid
    integer :: nx, ny


    nx = imo%grid%G%nx
    ny = imo%grid%G%ny

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,dim_x,x=imo%grid%G%x0,dx=imo%grid%G%dx,nx=nx,axis="x",units="m",ncid=ncid)
    call nc_write_dim(fnm,dim_y,x=imo%grid%G%y0,dx=imo%grid%G%dy,nx=ny,axis="y",units="m",ncid=ncid)

    call nc_write(fnm,"imo_ann", imo%imo_ann, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="floating ice basal melt",grid_mapping="polar_stereographic",units="kg/m2/a",ncid=ncid)    

    call nc_close(ncid)


  end subroutine imo_write_restart


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  i m o _ r e a d _ r e s t a r t
  ! Purpose  :  read imo restart file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine imo_read_restart(fnm,imo)

    implicit none

    type(imo_class), intent(inout) :: imo

    character (len=*) :: fnm


    call nc_read(fnm,"imo_ann", imo%imo_ann)


    return

  end subroutine imo_read_restart

end module imo_model
