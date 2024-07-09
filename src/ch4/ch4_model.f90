!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c h 4 _ m o d e l
!
!  Purpose : atmospheric CH4 model
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
module ch4_model

    use precision, only : wp
    use control, only: ch4_restart, restart_in_dir, ch4_const, i_ch4_tau, ch4_tau_const, ch4_ref
    use control, only : ich4_emis, ch4_emis_file, ch4_emis_const
    use timer, only : sec_year
    use ch4_def, only : ch4_class

    use ncio

    implicit none

    ! conversion factor from kgCH4 to ppb CH4
    real(wp), parameter :: kgCH4_to_ppb = 1._wp/2.7476_wp*1e-9_wp  ! ppb/TgCH4 * TgCH4/kgCH4 
    ! http://how-it-looks.blogspot.com/2011/07/petagrams-of-carbon.html
    ! The mean mass of the atmosphere is 5.1480e18 kg. The molar mass of air is 28.966 g/mol. 
    ! So the atmosphere contains 5.1480e21 g / 28.966 g/mol = 1.7773e20 moles of air.
    ! Mole fractions of carbon dioxide are expressed in ppm and directly convertible from parts-per-million by volume.
    ! So 1 ppb CH4 = 1.7773e20 / 1e9 = 1.7773e11 moles CH4
    ! The molar mass of CH4 is 16.04 g/mol.
    ! So 1 ppb CH4 = 1.7773e11 mol * 16.04 g/mol = 2.851 Tg CH4.
    ! 2.7476 TgCH4 per ppb used in IPCC (Prather et al., 2012)

    real(wp), dimension(:), allocatable :: ch4_emis_time, ch4_emis_data, f_ch4_emis_agro_data

    private
    public :: ch4_init, ch4_update, ch4_end
    public :: ch4_read_restart, ch4_write_restart


contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c h 4 _ u p d a t e
  ! Purpose  :  update atmospheric ch4
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ch4_update(ch4,time_now)

    implicit none

    type(ch4_class) :: ch4
    real(wp), intent(in) :: time_now

    integer :: imin, i0, i1
    real(wp) :: w0, w1
    real(wp) :: ch4_tau


    ! update ch4 emission rate if needed    
    if (ich4_emis.eq.1) then
      ! interpolate from ch4_emis_data to current year
      if (time_now.le.ch4_emis_time(lbound(ch4_emis_time,1))) then
        ch4%dch4emis_dt = ch4_emis_data(lbound(ch4_emis_data,1)) ! kgCH4
        ch4%f_ch4emis_agro = f_ch4_emis_agro_data(lbound(f_ch4_emis_agro_data,1)) 
      elseif (time_now.ge.ch4_emis_time(ubound(ch4_emis_time,1))) then
        ch4%dch4emis_dt = ch4_emis_data(ubound(ch4_emis_data,1)) ! kgCH4
        ch4%f_ch4emis_agro = f_ch4_emis_agro_data(ubound(f_ch4_emis_agro_data,1)) 
      else
        imin = minloc(abs(ch4_emis_time-time_now),1) 
        if (ch4_emis_time(imin).lt.time_now) then
          i0 = imin
          i1 = imin+1
        else
          i0 = imin-1
          i1 = imin
        endif
        w0 = 1._wp - abs(ch4_emis_time(i0)-time_now)/(ch4_emis_time(i1)-ch4_emis_time(i0))
        w1 = 1._wp - w0
        ch4%dch4emis_dt = (w0*ch4_emis_data(i0) + w1*ch4_emis_data(i1)) !  kgCH4/yr
        ch4%f_ch4emis_agro = (w0*f_ch4_emis_agro_data(i0) + w1*f_ch4_emis_agro_data(i1)) 
      endif
    endif

    if (i_ch4_tau.eq.1) then
      ! constant lifetime
      ch4%tau = ch4_tau_const
    else if (i_ch4_tau.eq.2) then
      ! lifetime dependent on CH4 load, approximate first order relation
      ! derived from simulations with MPI-ESM for deglaciation and SSP scenarios provided by Thomas Kleinen
      ch4%tau = ch4_tau_const*(ch4%ch4/ch4_ref)**0.24
    else if (i_ch4_tau.eq.3) then
      ! lifetime dependent on CH4 load, approximate first order relation
      ! derived from simulations with MPI-ESM for deglaciation and SSP scenarios provided by Thomas Kleinen
      ch4%tau = ch4_tau_const*(ch4%ch4_slow/ch4_ref)**0.24
    endif
    ch4%dch4ox_dt = ch4%ch4m/ch4%tau   ! kgCH4/yr

    ! add fluxes from ocean, land, anthropogenic emissions and decay (oxidation)
    ch4%ch4m   = ch4%ch4m   + ch4%dch4ocn_dt   + ch4%dch4lnd_dt   + ch4%dch4emis_dt  - ch4%dch4ox_dt ! kgCH4

    ! CH4 concentration
    ch4%ch4 = ch4%ch4m * kgCH4_to_ppb   ! ppb

    !ch4%ch4_slow = 0.99_wp*ch4%ch4_slow + 0.01_wp*ch4%ch4
    ch4%ch4_slow = 0.995_wp*ch4%ch4_slow + 0.005_wp*ch4%ch4

   return

  end subroutine ch4_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c h 4 _ i n i t
  ! Purpose  :  initialize atmospheric ch4
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ch4_init(ch4)

    implicit none

    type(ch4_class) :: ch4

    integer :: ntime


    if (ch4_restart) then

      ! read restart file
      call ch4_read_restart("restart/"//trim(restart_in_dir)//"/ch4_restart.nc",ch4)

    else

      ! initialise
      ch4%ch4 = ch4_const  ! ppb

    endif

    ch4%ch4m = ch4%ch4 / kgCH4_to_ppb   ! kgCH4

    ch4%ch4_slow = ch4%ch4

    ! emissions
    if (ich4_emis.eq.0) then

      ! constant ch4 emissions
      ch4%dch4emis_dt = ch4_emis_const * 1.e9_wp ! kgCH4/yr
      ch4%f_ch4emis_agro = 0._wp

    else if (ich4_emis.eq.1) then

      ! read ch4 emissions file
      ntime = nc_size(trim(ch4_emis_file),"time")
      allocate( ch4_emis_time(ntime) )
      allocate( ch4_emis_data(ntime) )
      allocate( f_ch4_emis_agro_data(ntime) )
      call nc_read(trim(ch4_emis_file),"time",ch4_emis_time)    ! time in years BP
      call nc_read(trim(ch4_emis_file),"ch4_emis",ch4_emis_data) ! TgCH4/yr
      ch4_emis_data = ch4_emis_data * 1.e9_wp  ! kgCH4/yr
      call nc_read(trim(ch4_emis_file),"f_ch4_emis_agro",f_ch4_emis_agro_data)

    endif

  return

  end subroutine ch4_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c h 4 _ e n d 
  ! Purpose  :  end atmospheric ch4
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ch4_end(ch4)

    implicit none

    type(ch4_class) :: ch4



   return

  end subroutine ch4_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c h 4 _ w r i t e _ r e s t a r t
  ! Purpose  :  Write restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ch4_write_restart(fnm,ch4)

    implicit none

    character (len=*) :: fnm
    type(ch4_class) :: ch4


    call nc_create(fnm)
    call nc_write_dim(fnm,"x",x=[1])
    call nc_write(fnm,"ch4",   ch4%ch4,     dim1="x",  long_name="atmospheric CH4 concentration",units="ppb")

   return

  end subroutine ch4_write_restart


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c h 4 _ r e a d _ r e s t a r t
  ! Purpose  :  read restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ch4_read_restart(fnm,ch4)

    implicit none

    character (len=*) :: fnm
    type(ch4_class) :: ch4


    call nc_read(fnm,"ch4",     ch4%ch4)

   return

  end subroutine ch4_read_restart


end module ch4_model

