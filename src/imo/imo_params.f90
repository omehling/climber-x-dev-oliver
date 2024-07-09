!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : i m o _ p a r a m s
!
!  Purpose : IMO model parameters
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
module imo_params

  use precision, only : wp
  use nml
  use timer, only : dt_imo, sec_year
  use control, only : out_dir

  implicit none

  real(wp) :: dt, rdt

  integer :: i_imo
  real(wp) :: imo_const
  real(wp) :: k_1
  real(wp) :: k_2

  logical :: l_fix_depth
  real(wp) :: fix_depth
  real(wp) :: depth_disc

  logical :: l_bm_lake

  logical :: l_depth_scale
  real(wp) :: zl_ref

  logical :: l_monthly_output

contains

    subroutine imo_params_init

    implicit none


    ! time step
    dt = dt_imo
    rdt = 1._wp/dt

    ! read imo parameter file
    call imo_par_load(trim(out_dir)//"/imo_par.nml")


    return

    end subroutine imo_params_init


subroutine imo_par_load(filename)

    implicit none

    character (len=*) :: filename


    ! Read parameters from file
    write(*,*) "imo parameters ==========="
    call nml_read(filename,"imo_par","i_imo",i_imo)
    call nml_read(filename,"imo_par","imo_const",imo_const)
    imo_const = imo_const/sec_year  ! kg/m2/a -> kg/m2/s
    call nml_read(filename,"imo_par","k_1",k_1)
    call nml_read(filename,"imo_par","k_2",k_2)
    call nml_read(filename,"imo_par","l_fix_depth",l_fix_depth)
    call nml_read(filename,"imo_par","fix_depth",fix_depth)
    call nml_read(filename,"imo_par","depth_disc",depth_disc)
    call nml_read(filename,"imo_par","l_bm_lake",l_bm_lake)
    call nml_read(filename,"imo_par","l_depth_scale",l_depth_scale)
    call nml_read(filename,"imo_par","zl_ref",zl_ref)
    call nml_read(filename,"imo_par","l_monthly_output",l_monthly_output)

   return

end subroutine imo_par_load


end module imo_params
