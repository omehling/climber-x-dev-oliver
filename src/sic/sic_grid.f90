!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s i c _ g r i d
!
!  Purpose : sea ice grid
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
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
module sic_grid

  use precision, only : wp
  use constants, only : pi, r_earth, omega
  use climber_grid, only: ni, nj

  implicit none

  integer :: maxi     !! number of gridcells in longitudinal direction []
  integer :: maxj     !! number of gridcells in latitudinal direction []

  real(wp) :: phi0   !! easternmost longitude [radians] 
  real(wp) :: dphi   !! longitudinal resolution [radians] 
  real(wp) :: rdphi  !! reverse of dphi [1/radians] 
  real(wp) :: dtheta !! latitudinal resolution [radians]
  real(wp) :: dy     !! latitudinal grid-cell width [m]
  real(wp) :: rdy    !! reverse of dy [1/m]
  real(wp) :: h1     !! thickness of top ocean layer [m]
  real(wp), dimension(:), allocatable :: s   !! sine of latitude at cell centers []
  real(wp), dimension(:), allocatable :: sv  !! sine of latitude at cell edges []
  real(wp), dimension(:), allocatable :: ds  !! Delta(s) []
  real(wp), dimension(:), allocatable :: dsv !! Delta(sv) []
  real(wp), dimension(:), allocatable :: rds !! reverse of ds
  real(wp), dimension(:), allocatable :: rdsv!! reverse of dsv
  real(wp), dimension(:), allocatable :: rds2!! 2.0/(dsv(j)+dsv(j-1))
  real(wp), dimension(:), allocatable :: c   !! cosine of latitude at cell centers []
  real(wp), dimension(:), allocatable :: cv  !! cosine of latitude at cell edges []
  real(wp), dimension(:), allocatable :: cv2 !! cv(j)*cv(j)*rdsv(j)
  real(wp), dimension(:), allocatable :: rc  !! reverse of c
  real(wp), dimension(:), allocatable :: rc2 !! rc*rc*rdphi
  real(wp), dimension(:), allocatable :: rcv !! reverse of cv
  real(wp), dimension(:), allocatable :: dx  !! longitudinal grid-cell width at cell centers [m] 
  real(wp), dimension(:), allocatable :: rdx !! reverse of dx
  real(wp), dimension(:), allocatable :: dxv !! longitudinal grid-cell width at cell edges [m]
  real(wp), dimension(:), allocatable :: rdxv!! reverse of dxv
  real(wp), dimension(:), allocatable :: fcor!! coriolis parameter []
  real(wp), dimension(:), allocatable :: fcorv!! coriolis parameter on v-grid []

  real(wp), dimension(:,:), allocatable :: area  !! horizonzal area of ocean cell fraction [m2]
  real(wp), dimension(:,:), allocatable :: area_old  !! horizonzal area of ocean cell fraction from previous year [m2]
  real(wp), dimension(:,:), allocatable :: area_full  !! horizonzal area of full grid cells [m2]

  real(wp), dimension(:,:), allocatable :: mask_u  !! mask on u-grid []
  real(wp), dimension(:,:), allocatable :: mask_v  !! mask on v-grid []
  real(wp), dimension(:,:), allocatable :: mask_q  !! mask on q-grid []

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s i c _ g r i d _ i n i t 
  !   Purpose    :  initialize sea ice grid
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_grid_init(f_ocn, dz1)

    implicit none

    real(wp), dimension(:,:), intent(in) :: f_ocn
    real(wp), intent(in) :: dz1

    integer :: i, j, ip1
    real(wp) :: th0, th1, s0, s1, theta, thv, dth, deg_to_rad


    maxi = ni   
    maxj = nj 
    h1 = dz1
    
    ! allocate variables
    allocate(s(0:maxj))
    allocate(sv(0:maxj))
    allocate(ds(maxj))
    allocate(dsv(1:maxj-1))
    allocate(rds(maxj))
    allocate(rds2(2:maxj-1))
    allocate(rdsv(1:maxj-1))
    allocate(c(0:maxj))
    allocate(cv(0:maxj))
    allocate(cv2(1:maxj-1))
    allocate(rc(0:maxj))
    allocate(rc2(0:maxj)) 
    allocate(rcv(1:maxj-1))
    allocate(dx(maxj))
    allocate(rdx(maxj))
    allocate(dxv(0:maxj))
    allocate(rdxv(0:maxj))
    allocate(area(maxi,maxj))
    allocate(area_old(maxi,maxj))
    allocate(area_full(maxi,maxj))
    allocate(mask_u(maxi,maxj))
    allocate(mask_v(maxi,maxj))
    allocate(mask_q(maxi,maxj))
    allocate(fcor(maxj))
    allocate(fcorv(0:maxj))

    ! parameters for setting up grid
    ! coords are latitude theta and longitude phi
    th0 = - pi/2._wp
    th1 = pi/2._wp
    s0 = sin(th0)
    s1 = sin(th1)
    deg_to_rad = pi/180._wp

    phi0 = -180.0*deg_to_rad
    dphi = 2._wp*pi/maxi
    rdphi = 1._wp/dphi

    ! set up horizontal grid: sin and cos factors at rho and v points (c grid)
    ! fix for global domain although only cv and cv2 are referred to at or beyond
    ! limits 24/6/2 if no flow out of N + S boundaries.

    sv(0) = s0
    cv(0) = cos(th0)
    ! set up const dlat grid
    dth = (th1 - th0)/maxj
    dtheta = dth
    do j=1,maxj
     thv = th0 + j*dth
     theta = thv - 0.5_wp*dth
     sv(j) = sin(thv)
     s(j) = sin(theta)
     cv(j) = cos(thv)
    enddo

    do j=1,maxj
     ds(j) = sv(j) - sv(j-1)
     rds(j) = 1._wp/ds(j)
     c(j) = sqrt(1._wp - s(j)*s(j))
     rc(j) = 1._wp/c(j)
     rc2(j) = rc(j)*rc(j)*rdphi
     if (j.lt.maxj) then
      dsv(j) = s(j+1) - s(j)
      rdsv(j) = 1._wp/dsv(j)
      cv2(j) = cv(j)*cv(j)*rdsv(j)
      rcv(j) = 1._wp/cv(j)
      if(j.gt.1) rds2(j) = 2._wp/(dsv(j)+dsv(j-1))
     endif
    enddo

    ! meridional grid resolution
    dy = dtheta*R_earth ! m
    rdy = 1._wp/dy
    ! zonal grid resolution at cell center (tracer points)
    do j=1,maxj
       dx(j) = dphi*R_earth*c(j)
       rdx(j) = 1._wp/dx(j)
    enddo
    ! zonal grid resolution on v-grid (meridional cell edges)
    do j=0,maxj
       dxv(j) = dphi*R_earth*cv(j)
       rdxv(j) = 1._wp/dxv(j)
    enddo

    ! area
    do j=1,maxj
      do i=1,maxi
        area_full(i,j) = dx(j)*dy
        if (f_ocn(i,j).gt.0._wp) then
          area(i,j) = dx(j)*dy*f_ocn(i,j)
        else
          area(i,j) = 0._wp
        endif
      enddo
    enddo

    ! u- and v-grid masks
    do j=1,maxj
      do i=1,maxi
        ! mask on u-grid
        ip1 = i+1
        if (ip1.eq.maxi+1) ip1=1
        if (f_ocn(i,j).gt.0._wp .and. f_ocn(ip1,j).gt.0._wp) then
          mask_u(i,j) = 1._wp
        else
          mask_u(i,j) = 0._wp
        endif
        if (j.eq.maxj) then
          mask_v(i,j) = 0._wp
          mask_q(i,j) = 0._wp
        else
          ! ocean mask on v-grid
          if (f_ocn(i,j).gt.0._wp .and. f_ocn(i,j+1).gt.0._wp) then
            mask_v(i,j) = 1._wp
          else
            mask_v(i,j) = 0._wp
          endif
          ! ocean mask on q-grid
          if (f_ocn(i,j).gt.0._wp .and. f_ocn(i,j+1).gt.0._wp .and. f_ocn(ip1,j).gt.0._wp .and. f_ocn(ip1,j+1).gt.0._wp) then
            mask_q(i,j) = 1._wp
          else
            mask_q(i,j) = 0._wp
          endif
        endif
      enddo
    enddo

    ! Coriolis parameter
    do j=1,maxj
      fcor(j) = sign(max(5.e-5_wp,abs(2._wp*omega*s(j))),s(j))
    enddo
    do j=0,maxj
      fcorv(j) = sign(max(5.e-5_wp,abs(2._wp*omega*sv(j))),sv(j))
    enddo


    return

  end subroutine sic_grid_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s i c _ g r i d _ u p d a t e
  !   Purpose    :  update sea ice grid
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_grid_update(f_ocn)

    implicit none

    real(wp), dimension(:,:), intent(in) :: f_ocn

    integer :: i, j, ip1


    ! save old area
    area_old = area

    ! update area
    do j=1,maxj
      do i=1,maxi
        area(i,j) = dx(j)*dy*f_ocn(i,j)
      enddo
    enddo

    ! update u- and v-grid masks
    do j=1,maxj
      do i=1,maxi
        ! mask on u-grid
        ip1 = i+1
        if (ip1.eq.maxi+1) ip1=1
        if (f_ocn(i,j).gt.0._wp .and. f_ocn(ip1,j).gt.0._wp) then
          mask_u(i,j) = 1._wp
        else
          mask_u(i,j) = 0._wp
        endif
        if (j.eq.maxj) then
          mask_v(i,j) = 0._wp
          mask_q(i,j) = 0._wp
        else
          ! update ocean mask on v-grid
          if (f_ocn(i,j).gt.0._wp .and. f_ocn(i,j+1).gt.0._wp) then
            mask_v(i,j) = 1._wp
          else
            mask_v(i,j) = 0._wp
          endif
          ! update ocean mask on q-grid
          if (f_ocn(i,j).gt.0._wp .and. f_ocn(i,j+1).gt.0._wp .and. f_ocn(ip1,j).gt.0._wp .and. f_ocn(ip1,j+1).gt.0._wp) then
            mask_q(i,j) = 1._wp
          else
            mask_q(i,j) = 0._wp
          endif
        endif
      enddo
    enddo

    return

  end subroutine sic_grid_update

end module sic_grid
