!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c r i s a _ m o d
!
!  Purpose : computation of cross-isobar angle
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Andrey Ganopolski and Matteo Willeit
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
module crisa_mod

  use atm_params, only : wp
  use constants, only : pi, karman
  use atm_params, only : cd0_ocn, cd0_sic, i_acbar, acbar_max, acbar_scale, nsmooth_acbar
  use atm_grid, only : im, jm, nm, fcorta_sqrt, i_ocn, i_sic, i_ice, i_lake
  use smooth_atm_mod, only : smooth2
  !$ use omp_lib

  implicit none

  private
  public :: crisa

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c r i s a
  !   Purpose    :  computation of cross-isobar angle
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine crisa(frst, z0m, zoro, &
        cd, cda, cd0, cd0a, acbar, sin_cos_acbar, cos_acbar, sin_acbar, epsa)

    implicit none

    real(wp), intent(in )  :: frst(:,:,:)
    real(wp), intent(in )  :: z0m(:,:,:)
    real(wp), intent(in )  :: zoro(:,:)

    real(wp), intent(out) :: cd(:,:,:)
    real(wp), intent(out) :: cda(:,:)
    real(wp), intent(out) :: cd0(:,:,:)
    real(wp), intent(out) :: cd0a(:,:)
    real(wp), intent(out) :: acbar(:,:)
    real(wp), intent(out) :: sin_cos_acbar(:,:)
    real(wp), intent(out) :: cos_acbar(:,:,:)
    real(wp), intent(out) :: sin_acbar(:,:,:)
    real(wp), intent(out) :: epsa(:,:,:)

    integer :: i, j, n, iter
    real(wp) :: rhs, rhsn, alfa0, alfa1, alfa, acbarn

    real(wp), parameter :: z_ref = 100._wp      ! m, reference height


    !$omp parallel do collapse(2) private(i, j, n, iter, rhs, rhsn, alfa0, alfa1, alfa, acbarn)
    do j=1,jm
      do  i=1,im

        !------------------------------------------------
        ! drag coefficient 
        !------------------------------------------------

        do n=1,nm

          !------------------------------------------------
          ! neutral drag coefficient

          if (n.eq.i_ocn) then
            ! ocean

            cd(i,j,n) = cd0_ocn
            cd0(i,j,n) = cd(i,j,n)

          else if (n.eq.i_sic) then
            ! sea ice

            cd(i,j,n) = cd0_sic
            cd0(i,j,n) = cd(i,j,n)

          else if (n.eq.i_lake) then
            ! lake

            cd(i,j,n) = cd0_ocn
            cd0(i,j,n) = cd(i,j,n)

          else
            ! land or ice sheet

            if (frst(i,j,n).gt.0._wp) then
              cd(i,j,n) = (karman/log(z_ref/(z0m(i,j,n)+zoro(i,j))))**2
              cd0(i,j,n)   = (karman/log(z_ref/z0m(i,j,n)))**2
            else
              cd(i,j,n) = 0.01_wp
              cd0(i,j,n) = 0.01_wp
            endif

          endif

        enddo

        ! grid cell average
        cda(i,j) = sum(cd(i,j,:)*frst(i,j,:)) 

        ! drag coefficient without orographic component
        cd0a(i,j) = sum(cd0(i,j,:)*frst(i,j,:)) 


        !------------------------------------------------
        ! solve of the equation for cross-isobar angle   
        !------------------------------------------------

        if (i_acbar.eq.1) then
          rhs = cda(i,j)/fcorta_sqrt(j) 
        else if (i_acbar.eq.2) then
          rhs = (frst(i,j,i_ice)*cd0a(i,j)+(1._wp-frst(i,j,i_ice))*cda(i,j)) /fcorta_sqrt(j) 
        endif
        alfa0 = 0._wp
        alfa1 = pi/4._wp
        ! Iteration loop
        do iter=1,10
          alfa = 0.5_wp*(alfa0+alfa1)
          rhsn = sin(alfa)/sqrt(1._wp-sin(2._wp*alfa))
          if (rhsn.gt.rhs) then
            alfa1 = alfa
          else
            alfa0 = alfa
          endif
        enddo          
        alfa = alfa*acbar_scale
        alfa = max(alfa,0.05_wp)
        alfa = min(alfa,acbar_max)
        acbar(i,j) = alfa
        sin_cos_acbar(i,j) = sin(alfa)*cos(alfa)

        ! Solving of the equation for cross-isobar angle for each surface type
        do n=1,nm
          rhs = cd(i,j,n)/fcorta_sqrt(j) 
          alfa0 = 0._wp
          alfa1 = pi/4._wp
          ! Iteration loop
          do iter=1,10
            alfa = 0.5_wp*(alfa0+alfa1)
            rhsn = sin(alfa)/sqrt(1._wp-sin(2._wp*alfa))
            if (rhsn.gt.rhs) then
              alfa1 = alfa
            else
              alfa0 = alfa
            endif
          enddo          
          alfa = alfa*acbar_scale
          alfa = max(alfa,0.05_wp)
          alfa = min(alfa,acbar_max)
          acbarn = alfa
          cos_acbar(i,j,n) = cos(acbarn)
          sin_acbar(i,j,n) = sin(acbarn)
          epsa(i,j,n) = sqrt(1._wp-sin(2._wp*acbarn))

        enddo

      enddo
    enddo
    !$omp end parallel do

    ! smooth in space
    call smooth2(acbar,nsmooth_acbar)
    call smooth2(sin_cos_acbar,nsmooth_acbar)

    return

  end subroutine crisa

end module crisa_mod
