!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : b i a s _ c o r r _ m o d 
!
!  Purpose : climate bias correction
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
module bias_corr_mod

   use precision, only : wp
   use timer, only: doy, nday_year, monthly2daily
   use ncio

   implicit none

   real(wp), dimension(:,:), allocatable :: t2m_bias_ann
   real(wp), dimension(:,:,:), allocatable :: t2m_bias_mon
   real(wp), dimension(:,:,:), allocatable :: prc_bias_mon
   integer, dimension(nday_year) :: m0, m1
   real(wp), dimension(nday_year) :: wtm0, wtm1

   private
   public :: t2m_bias_corr, prc_bias_corr

contains

   subroutine t2m_bias_corr(i_t2m_bias_corr, bias_corr_file, t2m_bias)

   implicit none

   integer, intent(in) :: i_t2m_bias_corr
   character (len=*), intent(in) :: bias_corr_file

   real(wp), intent(out) :: t2m_bias(:,:,:)

   integer :: d, nmon, nlon, nlat


   if (i_t2m_bias_corr.eq.1) then

     ! read temperature bias file
     nlon = nc_size(trim(bias_corr_file),"lon")
     nlat = nc_size(trim(bias_corr_file),"lat")
     allocate( t2m_bias_ann(nlon,nlat) )
     call nc_read(trim(bias_corr_file),"t2m_bias",t2m_bias_ann) 

     do d=1,nday_year
       t2m_bias(:,:,d) = t2m_bias_ann 
     enddo

   else if (i_t2m_bias_corr.eq.2) then

     ! get weights for interpolation from monthly to daily
     call monthly2daily(m0,m1,wtm0,wtm1)

     ! read temperature bias file
     nmon = nc_size(trim(bias_corr_file),"mon")
     if (nmon.gt.13) stop 'bias_corr ERROR: i_bias_corr==2 requires input file with monthly values' 
     nlon = nc_size(trim(bias_corr_file),"lon")
     nlat = nc_size(trim(bias_corr_file),"lat")
     allocate( t2m_bias_mon(nlon,nlat,0:nmon+1) )
     call nc_read(trim(bias_corr_file),"t2m_bias",t2m_bias_mon(:,:,1:12),start=[1,1,1],count=[nlon,nlat,12]) 

     ! expand month for easier interpolation at borders
     t2m_bias_mon(:,:,0)  = t2m_bias_mon(:,:,12)
     t2m_bias_mon(:,:,13) = t2m_bias_mon(:,:,1)

     do d=1,nday_year
       t2m_bias(:,:,d) = wtm0(d)*t2m_bias_mon(:,:,m0(d)) + wtm1(d)*t2m_bias_mon(:,:,m1(d))
     enddo

   endif

   return

   end subroutine t2m_bias_corr


   subroutine prc_bias_corr(i_prc_bias_corr, bias_corr_file, prc_bias)

   implicit none

   integer, intent(in) :: i_prc_bias_corr
   character (len=*), intent(in) :: bias_corr_file

   real(wp), intent(out) :: prc_bias(:,:,:)

   integer :: d, nmon, nlon, nlat


   if (i_prc_bias_corr.eq.1) then

     ! read temperature bias file
     nlon = nc_size(trim(bias_corr_file),"lon")
     nlat = nc_size(trim(bias_corr_file),"lat")
     nmon = nc_size(trim(bias_corr_file),"mon")
     allocate( prc_bias_mon(nlon,nlat,nmon) )
     call nc_read(trim(bias_corr_file),"prc_bias",prc_bias_mon) 

     do d=1,nday_year
       prc_bias(:,:,d) = sum(prc_bias_mon(:,:,5:9),3)/5._wp 
     enddo

   else if (i_prc_bias_corr.eq.2) then

     ! get weights for interpolation from monthly to daily
     call monthly2daily(m0,m1,wtm0,wtm1)

     ! read temperature bias file
     nmon = nc_size(trim(bias_corr_file),"mon")
     if (nmon.gt.13) stop 'bias_corr ERROR: i_bias_corr==2 requires input file with monthly values' 
     nlon = nc_size(trim(bias_corr_file),"lon")
     nlat = nc_size(trim(bias_corr_file),"lat")
     allocate( prc_bias_mon(nlon,nlat,0:nmon+1) )
     call nc_read(trim(bias_corr_file),"prc_bias",prc_bias_mon(:,:,1:12),start=[1,1,1],count=[nlon,nlat,12]) 

     ! expand month for easier interpolation at borders
     prc_bias_mon(:,:,0)  = prc_bias_mon(:,:,12)
     prc_bias_mon(:,:,13) = prc_bias_mon(:,:,1)

     do d=1,nday_year
       prc_bias(:,:,d) = wtm0(d)*prc_bias_mon(:,:,m0(d)) + wtm1(d)*prc_bias_mon(:,:,m1(d))
     enddo

   endif

   return

   end subroutine prc_bias_corr

end module bias_corr_mod
