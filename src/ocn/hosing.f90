module hosing_mod

  use precision, only : wp, sp
  use constants, only : pi
  use timer, only: sec_year, time_soy_ocn
  use climber_grid, only : lat, lon, basin_mask, i_atlantic, i_pacific, i_southern
  use ocn_grid, only : maxi, maxj, dx, dy
  use ocn_params, only: dt, rho0, i_hosing, hosing_basin, hosing_ini, hosing_trend, hosing_sigma
  use ocn_params, only : year_hosing_ini, year_hosing_end, year_hosing_ramp, lat_min_hosing, lat_max_hosing, lon_min_hosing, lon_max_hosing
  use ocn_params, only : hosing_comp_basin, lat_min_hosing_comp, lat_max_hosing_comp, lon_min_hosing_comp, lon_max_hosing_comp

!  use hyster 
  
  implicit none

  integer :: i_hosing_1, i_hosing_2
  integer :: i_hosing_comp_1, i_hosing_comp_2
  integer , allocatable :: j_hosing(:,:)
  integer , allocatable :: j_hosing_comp(:,:)
  real(wp), allocatable :: rhosing(:,:)

!  type(hyster_class) :: hyst1 

  private
  public :: hosing_init, hosing_update

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  h o s i n g _ u p d a t e
  !   Purpose    :  update freshwater forcing
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine hosing_update(time,f_ocn,amoc,hosing,fw_hosing) 

    implicit none

    real(wp), intent(in) :: time                          ! [yr] Current model time 
    real(wp), dimension(:,:), intent(in) :: f_ocn
    real(wp), intent(in) :: amoc                          ! [Sv] current amoc strength
    real(wp), intent(out) :: hosing
    real(wp), dimension(:,:), intent(inout) :: fw_hosing

    integer :: i, j 
    real(wp) :: area_hosing
    real(wp) :: area_hosing_comp
    real(wp), dimension(2) :: xx
    real(wp) :: norm


    rhosing = 0._wp
    area_hosing = 0._wp

    ! compute total hosing area
    do i=i_hosing_1,i_hosing_2
      do j=j_hosing(i,1),j_hosing(i,2)
        if (hosing_basin.eq.1) then
          ! Atlantic
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_atlantic) area_hosing = area_hosing + dx(j)*dy*f_ocn(i,j)
        else if (hosing_basin.eq.2) then
          ! Pacific
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_pacific) area_hosing = area_hosing + dx(j)*dy*f_ocn(i,j)
        else if (hosing_basin.eq.3) then
          ! Southern ocean
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_southern) area_hosing = area_hosing + dx(j)*dy*f_ocn(i,j)
        endif
      enddo
    enddo

    ! derive hosing weight for each cell
    do i=i_hosing_1,i_hosing_2
      do j=j_hosing(i,1),j_hosing(i,2)
        if (hosing_basin.eq.1) then
          ! Atlantic
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_atlantic) rhosing(i,j) = 1.e6_wp/area_hosing  ! m3/s/Sv / m2 = m/s/Sv
        else if (hosing_basin.eq.2) then
          ! Pacific
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_pacific) rhosing(i,j) = 1.e6_wp/area_hosing
        else if (hosing_basin.eq.3) then
          ! Southern ocean
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_southern) rhosing(i,j) = 1.e6_wp/area_hosing  ! m3/s/Sv / m2 = m/s/Sv
        endif
      enddo
    enddo

    ! compute total hosing compensation area
    do i=1,maxi
      do j=1,maxj
        if (hosing_comp_basin.eq.0) then
          ! Global ocean
          if (f_ocn(i,j).gt.0._wp) area_hosing_comp = area_hosing_comp + dx(j)*dy*f_ocn(i,j)
        endif
      enddo
    enddo
    do i=i_hosing_comp_1,i_hosing_comp_2
      do j=j_hosing_comp(i,1),j_hosing_comp(i,2)
        if (hosing_comp_basin.eq.1) then
          ! Atlantic
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_atlantic) area_hosing_comp = area_hosing_comp + dx(j)*dy*f_ocn(i,j)
        else if (hosing_comp_basin.eq.2) then
          ! Pacific
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_pacific) area_hosing_comp = area_hosing_comp + dx(j)*dy*f_ocn(i,j)
        else if (hosing_comp_basin.eq.3) then
          ! Southern ocean
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_southern) area_hosing_comp = area_hosing_comp + dx(j)*dy*f_ocn(i,j)
        endif
      enddo
    enddo
    if (hosing_comp_basin.eq.0) then
      area_hosing_comp = area_hosing_comp - area_hosing
    endif

    ! derive hosing compensation weight for each cell
    do i=i_hosing_comp_1,i_hosing_comp_2
      do j=j_hosing_comp(i,1),j_hosing_comp(i,2)
        if (hosing_comp_basin.eq.1) then
          ! Atlantic
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_atlantic) rhosing(i,j) = -1.e6_wp/area_hosing_comp  ! m3/s/Sv / m2 = m/s/Sv
        else if (hosing_comp_basin.eq.2) then
          ! Pacific
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_pacific) rhosing(i,j) = -1.e6_wp/area_hosing_comp
        else if (hosing_comp_basin.eq.3) then
          ! Southern ocean
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_southern) rhosing(i,j) = -1.e6_wp/area_hosing_comp  ! m3/s/Sv / m2 = m/s/Sv
        endif
      enddo
    enddo
    if (hosing_comp_basin.eq.0) then
      where (rhosing.eq.0._wp)
        rhosing = -1e6_wp/area_hosing_comp ! m3/s/Sv / m2 = m/s/Sv
      endwhere
    endif


if (.true.) then 
  ! Old hosing method with constant hosing_trend value 

    if (hosing_sigma.gt.0._wp) then
      ! Get 2 numbers from uniform distribution between 0 and 1
      call random_number(xx)
      ! Convert to normal distribution using the Box-Mueller algorithm
      norm = sqrt(-2._wp*log(xx(1))) * cos(2._wp*pi*xx(2))
      ! keep within 3 sigma
      norm = min(3._wp,norm)
      norm = max(-3._wp,norm)
    else
      norm = 0._wp
    endif

    ! TODO add option to convert to red noise?
    ! https://atmos.washington.edu/~breth/classes/AM582/lect/lect8-notes.pdf
    if (i_hosing.eq.0) then

      if (time.ge.year_hosing_ini .and. time.le.year_hosing_end) then
        ! update total extra freshwater forcing
        hosing = min(1._wp,time/max(1,year_hosing_ramp))*hosing_ini + hosing_trend*sec_year*(time-year_hosing_ini)  + hosing_sigma*norm 
      else 
        hosing = 0._wp
      endif

    else if (i_hosing.eq.1) then

      if (time.lt.year_hosing_ini) then
        hosing = 0._wp
      else if (time.ge.year_hosing_ini .and. time.lt.year_hosing_end) then
        hosing = hosing_ini
      else if (time.ge.year_hosing_end .and. time.lt.(year_hosing_end+(year_hosing_end-year_hosing_ini))) then
        hosing = -hosing_ini
      else
        hosing = 0._wp
      endif

    else if (i_hosing.eq.2) then

      if (time.lt.year_hosing_ini) then
        hosing = 0._wp
      else if (time.ge.year_hosing_ini .and. time.lt.year_hosing_end) then
        hosing = hosing_ini
      else if (time.ge.year_hosing_end .and. time.lt.(year_hosing_end+(year_hosing_end-year_hosing_ini))) then
        hosing = -0.5_wp*hosing_ini
      else
        hosing = 0._wp
      endif

    else if (i_hosing.eq.3) then

      if (time.lt.year_hosing_ini) then
        hosing = 0._wp
      else if (time.ge.year_hosing_ini .and. time.lt.year_hosing_end) then
        hosing = hosing_ini
      else if (time.ge.year_hosing_end .and. time.lt.(year_hosing_end+(year_hosing_end-year_hosing_ini))) then
        hosing = -2._wp*hosing_ini
      else
        hosing = 0._wp
      endif

    else if (i_hosing.eq.4) then

      if (time.lt.year_hosing_ini) then
        hosing = 0._wp
      else if (time.ge.year_hosing_ini .and. time.lt.year_hosing_end) then
        hosing = hosing_ini
      else if (time.ge.year_hosing_end .and. time.lt.(year_hosing_end+(year_hosing_end-year_hosing_ini))) then
        hosing = -hosing_ini
      else
        hosing = 0._wp
      endif

    endif

else
!  ! New hosing method using hyster object 
!
!    ! Update hysteresis variable 
!    call hyster_calc_forcing(hyst1,time=real(time,sp),var=real(amoc,sp))
!
!    ! update total extra freshwater forcing [Sv]
!    hosing = hyst1%f_now 
!
end if 

    ! update 2D extra freshwater forcing
    fw_hosing = hosing * rhosing*rho0  ! Sv * m/s/Sv * kg/m3 = kg/m2/s


   return

  end subroutine hosing_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  h o s i n g _ i n i t
  !   Purpose    :  initialize freshwater forcing
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine hosing_init(f_ocn,time)

    implicit none

    real(wp), intent(in) :: f_ocn(:,:)
    real(wp), intent(in) :: time

    integer :: i, j
    logical :: flag
    real(wp) :: area_hosing
    real(wp) :: area_hosing_comp


    allocate(rhosing(maxi,maxj))
    allocate(j_hosing(maxi,2))
    allocate(j_hosing_comp(maxi,2))

    ! find index of western and eastern boundary of box where hosing should be applied
    do i=2,maxi
     ! western boundary
     if ((lon(i-1).le.lon_min_hosing) .and. (lon(i).ge.lon_min_hosing)) then
        i_hosing_1 = i-1
     endif
     ! eastern boundary
     if ((lon(i-1).le.lon_max_hosing) .and. (lon(i).ge.lon_max_hosing)) then
        i_hosing_2 = i
     endif
    enddo

    ! find index of western and eastern boundary of box where hosing compensation should be applied
    do i=2,maxi
     ! western boundary
     if ((lon(i-1).le.lon_min_hosing_comp) .and. (lon(i).ge.lon_min_hosing_comp)) then
        i_hosing_comp_1 = i-1
     endif
     ! eastern boundary
     if ((lon(i-1).le.lon_max_hosing_comp) .and. (lon(i).ge.lon_max_hosing_comp)) then
        i_hosing_comp_2 = i
     endif
    enddo

    ! find index of northern and southern boundary of box where hosing should be applied
    do j=2,maxj
     ! southern boundary
     if ((lat(j-1).le.lat_min_hosing) .and. (lat(j).ge.lat_min_hosing)) then
        j_hosing(:,1) = j
     endif
     ! northern boundary
     if ((lat(j-1).le.lat_max_hosing) .and. (lat(j).ge.lat_max_hosing)) then
        j_hosing(:,2) = j-1
     endif
    enddo

    ! find index of northern and southern boundary of box where hosing compensation should be applied
    do j=2,maxj
     ! southern boundary
     if ((lat(j-1).le.lat_min_hosing_comp) .and. (lat(j).ge.lat_min_hosing_comp)) then
        j_hosing_comp(:,1) = j
     endif
     ! northern boundary
     if ((lat(j-1).le.lat_max_hosing_comp) .and. (lat(j).ge.lat_max_hosing_comp)) then
        j_hosing_comp(:,2) = j-1
     endif
    enddo

    ! special case, first ocean cell around Antarctica, if lat_min_hosing==0 and lat_max_hosing==0
    if (lat_min_hosing==0 .and. lat_max_hosing==0) then
      do i=1,maxi
        flag = .false.
        do j=1,maxj
          if ( (.not. flag) .and. f_ocn(i,j).gt.0._wp) then  ! first ocean cell from south
            j_hosing(i,1:2) = j
            flag = .true.
          endif
        enddo
      enddo
    endif

    rhosing = 0._wp
    area_hosing = 0._wp

    ! compute total hosing area
    do i=i_hosing_1,i_hosing_2
      do j=j_hosing(i,1),j_hosing(i,2)
        if (hosing_basin.eq.1) then
          ! Atlantic
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_atlantic) area_hosing = area_hosing + dx(j)*dy*f_ocn(i,j)
        else if (hosing_basin.eq.2) then
          ! Pacific
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_pacific) area_hosing = area_hosing + dx(j)*dy*f_ocn(i,j)
        else if (hosing_basin.eq.3) then
          ! Southern ocean
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_southern) area_hosing = area_hosing + dx(j)*dy*f_ocn(i,j)
        endif
      enddo
    enddo

    ! derive hosing weight for each cell
    do i=i_hosing_1,i_hosing_2
      do j=j_hosing(i,1),j_hosing(i,2)
        if (hosing_basin.eq.1) then
          ! Atlantic
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_atlantic) rhosing(i,j) = 1.e6_wp/area_hosing  ! m3/s/Sv / m2 = m/s/Sv
        else if (hosing_basin.eq.2) then
          ! Pacific
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_pacific) rhosing(i,j) = 1.e6_wp/area_hosing
        else if (hosing_basin.eq.3) then
          ! Southern ocean
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_southern) rhosing(i,j) = 1.e6_wp/area_hosing  ! m3/s/Sv / m2 = m/s/Sv
        endif
      enddo
    enddo

    ! compute total hosing compensation area
    do i=1,maxi
      do j=1,maxj
        if (hosing_comp_basin.eq.0) then
          ! Global ocean
          if (f_ocn(i,j).gt.0._wp) area_hosing_comp = area_hosing_comp + dx(j)*dy*f_ocn(i,j)
        endif
      enddo
    enddo
    do i=i_hosing_comp_1,i_hosing_comp_2
      do j=j_hosing_comp(i,1),j_hosing_comp(i,2)
        if (hosing_comp_basin.eq.1) then
          ! Atlantic
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_atlantic) area_hosing_comp = area_hosing_comp + dx(j)*dy*f_ocn(i,j)
        else if (hosing_comp_basin.eq.2) then
          ! Pacific
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_pacific) area_hosing_comp = area_hosing_comp + dx(j)*dy*f_ocn(i,j)
        else if (hosing_comp_basin.eq.3) then
          ! Southern ocean
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_southern) area_hosing_comp = area_hosing_comp + dx(j)*dy*f_ocn(i,j)
        endif
      enddo
    enddo
    if (hosing_comp_basin.eq.0) then
      area_hosing_comp = area_hosing_comp - area_hosing
    endif

    ! derive hosing compensation weight for each cell
    do i=i_hosing_comp_1,i_hosing_comp_2
      do j=j_hosing_comp(i,1),j_hosing_comp(i,2)
        if (hosing_comp_basin.eq.1) then
          ! Atlantic
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_atlantic) rhosing(i,j) = -1.e6_wp/area_hosing_comp  ! m3/s/Sv / m2 = m/s/Sv
        else if (hosing_comp_basin.eq.2) then
          ! Pacific
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_pacific) rhosing(i,j) = -1.e6_wp/area_hosing_comp
        else if (hosing_comp_basin.eq.3) then
          ! Southern ocean
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_southern) rhosing(i,j) = -1.e6_wp/area_hosing_comp  ! m3/s/Sv / m2 = m/s/Sv
        endif
      enddo
    enddo
    if (hosing_comp_basin.eq.0) then
      where (rhosing.eq.0._wp)
        rhosing = -1e6_wp/area_hosing_comp ! m3/s/Sv / m2 = m/s/Sv
      endwhere
    endif


    ! convert hosing_trend from Sv/ky to Sv/s
    hosing_trend = hosing_trend/(1.e3_wp*sec_year)

    
!    ! Initialize hysteresis module for transient forcing experiments 
!    call hyster_init(hyst1,"hyster_ctrl.nml",real(time,sp),"hosing") 

    return

  end subroutine hosing_init

end module hosing_mod
