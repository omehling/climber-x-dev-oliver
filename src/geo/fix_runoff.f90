!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : f i x _ r u n o f f _ m o d
!
!  Purpose : make sure runoff end up in ocean cell
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
module fix_runoff_mod

  use precision, only : wp

  implicit none

  private
  public :: fix_runoff, fix_runoff_lake

contains


  subroutine fix_runoff(f_ocn, &
      i_runoff, j_runoff)

    implicit none

    real(wp), intent(in) :: f_ocn(:,:)
    integer, intent(inout) :: i_runoff(:,:)
    integer, intent(inout) :: j_runoff(:,:)

    integer :: i, j, ii, jj, n, k, ni, nj
    integer :: i0, j0, ip1, im1, jp1, jm1, ip2, im2, jp2, jm2  
    integer, allocatable, dimension(:) :: idx, jdx


    ni = size(f_ocn,1)
    nj = size(f_ocn,2)

    ! make sure runoff does not end up in dry cells

    allocate(idx(25))
    allocate(jdx(25))

    do i=1,ni
      do j=1,nj
        if (i_runoff(i,j).ne.0 .and. j_runoff(i,j).ne.0) then
          if (f_ocn(i_runoff(i,j),j_runoff(i,j)).eq.0._wp) then ! runoff to dry point, move to neighbor ocean point
            !print *
            !print '(5F7.1)', f_ocn(max(1,i_runoff(i,j)-2):min(ni,i_runoff(i,j)+2),min(nj,j_runoff(i,j)+2))
            !print '(5F7.1)', f_ocn(max(1,i_runoff(i,j)-2):min(ni,i_runoff(i,j)+2),min(nj,j_runoff(i,j)+1))
            !print '(5F7.1)', f_ocn(max(1,i_runoff(i,j)-2):min(ni,i_runoff(i,j)+2),j_runoff(i,j))
            !print '(5F7.1)', f_ocn(max(1,i_runoff(i,j)-2):min(ni,i_runoff(i,j)+2),max(1,j_runoff(i,j)-1))
            !print '(5F7.1)', f_ocn(max(1,i_runoff(i,j)-2):min(ni,i_runoff(i,j)+2),max(1,j_runoff(i,j)-2))
            !print *
            !print *,'dry point in',i_runoff(i,j),j_runoff(i,j),f_ocn(i_runoff(i,j),j_runoff(i,j))
            i0 = i_runoff(i,j)
            j0 = j_runoff(i,j)
            ip1 = i0+1
            if (ip1.eq.ni+1) ip1 = 1
            im1 = i0-1
            if (im1.eq.0) im1 = ni
            jp1 = j0+1
            if (jp1.eq.nj+1) jp1 = nj
            jm1 = j0-1
            if (jm1.eq.0) jm1 = 1
            ip2 = i0+2
            if (ip2.eq.ni+1) ip2 = 1
            if (ip2.eq.ni+2) ip2 = 2
            im2 = i0-2
            if (im2.eq.0) im2 = ni
            if (im2.eq.-1) im2 = ni-1
            jp2 = j0+2
            if (jp2.gt.nj) jp2 = nj
            jm2 = j0-2
            if (jm2.lt.1) jm2 = 1
            idx = (/i0, ip1, im1, i0,  i0,  im1, im1, ip1, ip1, ip2, im2, i0,  i0,  im2, im2, ip2, ip2, im1, ip1, im1, ip1, im2, im2, ip2, ip2/)
            jdx = (/j0, j0,  j0,  jp1, jm1, jm1, jp1, jm1, jp1, j0,  j0,  jp2, jm2, jm1, jp1, jm1, jp1, jm2, jm2, jp2, jp2, jm2, jp2, jm2, jp2/)
            n = 0
            do while (f_ocn(i_runoff(i,j),j_runoff(i,j)).eq.0._wp)
              n = n+1
              if (n.le.size(idx)) then
                i_runoff(i,j) = idx(n)
                j_runoff(i,j) = jdx(n)
              else
                ! still dry point, put runoff 'somewhere close' 
                print *,'WARNING, runoff to dry point',i,j,i0,j0
                loopj : do k=0,nj
                  do n=1,2
                    jj=j0+k*(-1.)**n
                    if (jj.lt.1 .or. jj.gt.nj) cycle
                    do ii=i0,ni
                      if (f_ocn(ii,jj).gt.0._wp) then
                        i_runoff(i,j) = ii
                        j_runoff(i,j) = jj
                        exit loopj
                      endif
                    enddo
                    do ii=i0,1,-1
                      if (f_ocn(ii,jj).gt.0._wp) then
                        i_runoff(i,j) = ii
                        j_runoff(i,j) = jj
                        exit loopj
                      endif
                    enddo
                  enddo
                enddo loopj
                print *,'runoff from dry point to',i_runoff(i,j),j_runoff(i,j)
              endif
            enddo
          endif
        endif
      enddo
    enddo

    deallocate(idx)
    deallocate(jdx)

  end subroutine fix_runoff


  subroutine fix_runoff_lake(f_ocn, &
      i_runoff, j_runoff)

    implicit none

    real(wp), intent(in) :: f_ocn(:,:)
    integer, intent(inout) :: i_runoff
    integer, intent(inout) :: j_runoff

    integer :: n, ni, nj, k
    integer :: i0, j0, ii, jj, ip1, im1, jp1, jm1, ip2, im2, jp2, jm2  
    integer, dimension(25) :: idx, jdx


    ni = size(f_ocn,1)
    nj = size(f_ocn,2)

    if (f_ocn(i_runoff,j_runoff).eq.0._wp) then ! runoff to dry point, move to neighbor ocean point
      i0 = i_runoff
      j0 = j_runoff
      ip1 = i0+1
      if (ip1.eq.ni+1) ip1 = 1
      im1 = i0-1
      if (im1.eq.0) im1 = ni
      jp1 = j0+1
      if (jp1.eq.nj+1) jp1 = nj
      jm1 = j0-1
      if (jm1.eq.0) jm1 = 1
      ip2 = i0+2
      if (ip2.eq.ni+1) ip2 = 1
      if (ip2.eq.ni+2) ip2 = 2
      im2 = i0-2
      if (im2.eq.0) im2 = ni
      if (im2.eq.-1) im2 = ni-1
      jp2 = j0+2
      if (jp2.gt.nj) jp2 = nj
      jm2 = j0-2
      if (jm2.lt.1) jm2 = 1
      idx = (/i0, ip1, im1, i0,  i0,  im1, im1, ip1, ip1, ip2, im2, i0,  i0,  im2, im2, ip2, ip2, im1, ip1, im1, ip1, im2, im2, ip2, ip2/)
      jdx = (/j0, j0,  j0,  jp1, jm1, jm1, jp1, jm1, jp1, j0,  j0,  jp2, jm2, jm1, jp1, jm1, jp1, jm2, jm2, jp2, jp2, jm2, jp2, jm2, jp2/)
      n = 0
      do while (f_ocn(i_runoff,j_runoff).eq.0._wp)
        n = n+1
        if (n.le.size(idx)) then
          i_runoff = idx(n)
          j_runoff = jdx(n)
        else
          ! still dry point, put runoff 'somewhere close' 
          print *,'WARNING, runoff to dry point (lake)',i0,j0
          loopj : do k=0,nj
            do n=1,2
              jj=j0+k*(-1.)**n
              if (jj.lt.1 .or. jj.gt.nj) cycle
              do ii=i0,ni
                if (f_ocn(ii,jj).gt.0._wp) then
                  i_runoff = ii
                  j_runoff = jj
                  exit loopj
                endif
              enddo
              do ii=i0,1,-1
                if (f_ocn(ii,jj).gt.0._wp) then
                  i_runoff = ii
                  j_runoff = jj
                  exit loopj
                endif
              enddo
            enddo
          enddo loopj
          print *,'runoff from dry point (lake) to',i_runoff,j_runoff
        endif
      enddo
    endif

  end subroutine fix_runoff_lake

end module fix_runoff_mod

