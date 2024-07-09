module noise_mod

  use precision, only : wp
  use constants, only : pi
  use climber_grid, only : lat, lon, basin_mask, i_atlantic, i_pacific, i_southern
  use ocn_grid, only : maxi, maxj
  use ocn_params, only: i_noise, noise_basin, noise_amp_fw, noise_amp_flx, noise_period, noise_autocorr, l_noise_only_coast
  use ocn_params, only : lat_min_noise, lat_max_noise, lon_min_noise, lon_max_noise

  implicit none

  integer :: i_noise_1, i_noise_2
  integer :: i_noise_comp_1, i_noise_comp_2
  integer, allocatable :: j_noise(:,:)
  integer, allocatable :: mask_noise(:,:)

  private
  public :: noise_init, noise_fw_update, noise_flx_update

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  n o i s e _ f w _ u p d a t e
  !   Purpose    :  update freshwater noise
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine noise_fw_update(time,f_ocn,mask_coast,noise_fw,fw_noise) 

    implicit none

    real(wp), intent(in) :: time                          ! [yr] Current model time 
    real(wp), dimension(:,:), intent(in) :: f_ocn
    integer, intent(in) :: mask_coast(:,:)
    real(wp), intent(inout) :: noise_fw
    real(wp), dimension(:,:), intent(inout) :: fw_noise

    integer :: i, j 
    real(wp), dimension(2) :: xx
    real(wp) :: norm
    real(wp) :: noise
    real(wp), save :: noise_old = 0._wp


    if (i_noise.eq.1) then
      ! Gaussian white noise

      ! Get 2 numbers from uniform distribution between 0 and 1
      call random_number(xx)
      ! Convert to normal distribution using the Box-Mueller algorithm
      norm = sqrt(-2._wp*log(xx(1))) * cos(2._wp*pi*xx(2))
      ! keep within 3 sigma
      norm = min(3._wp,norm)
      norm = max(-3._wp,norm)

      noise = noise_amp_fw*norm

    else if (i_noise.eq.2) then
      ! red noise

      ! Get 2 numbers from uniform distribution between 0 and 1
      call random_number(xx)
      ! Convert to normal distribution using the Box-Mueller algorithm
      norm = sqrt(-2._wp*log(xx(1))) * cos(2._wp*pi*xx(2))
      ! keep within 3 sigma
      norm = min(3._wp,norm)
      norm = max(-3._wp,norm)

      ! convert white noise to red noise following:
      ! https://atmos.washington.edu/~breth/classes/AM582/lect/lect8-notes.pdf
      norm = noise_autocorr*noise_old + sqrt(1._wp-noise_autocorr**2)*norm
      noise_old = norm

      noise = noise_amp_fw*norm

    else if (i_noise.eq.3) then
      ! periodic 'noise'

      noise = noise_amp_fw*cos(2._wp*pi*time/noise_period)

    endif

    noise_fw = noise

    ! update 2D freshwater noise
    do j=1,maxj
      do i=1,maxi
        if (mask_noise(i,j).eq.1) then
          if (.not.l_noise_only_coast .or. mask_coast(i,j).eq.1) then
            fw_noise(i,j) = noise  ! kg/m2/s
          else
            fw_noise(i,j) = 0._wp 
          endif
        else
          fw_noise(i,j) = 0._wp 
        endif
      enddo
    enddo


   return

  end subroutine noise_fw_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  n o i s e _ f l x _ u p d a t e
  !   Purpose    :  update heat flux noise
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine noise_flx_update(time,f_ocn,mask_coast,noise_flx,flx_noise) 

    implicit none

    real(wp), intent(in) :: time                          ! [yr] Current model time 
    real(wp), dimension(:,:), intent(in) :: f_ocn
    integer, intent(in) :: mask_coast(:,:)
    real(wp), intent(inout) :: noise_flx
    real(wp), dimension(:,:), intent(inout) :: flx_noise

    integer :: i, j 
    real(wp), dimension(2) :: xx
    real(wp) :: norm
    real(wp) :: noise
    real(wp), save :: noise_old = 0._wp


    if (i_noise.eq.1) then
      ! Gaussian white noise

      ! Get 2 numbers from uniform distribution between 0 and 1
      call random_number(xx)
      ! Convert to normal distribution using the Box-Mueller algorithm
      norm = sqrt(-2._wp*log(xx(1))) * cos(2._wp*pi*xx(2))
      ! keep within 3 sigma
      norm = min(3._wp,norm)
      norm = max(-3._wp,norm)

      noise = noise_amp_flx*norm

    else if (i_noise.eq.2) then
      ! red noise

      ! Get 2 numbers from uniform distribution between 0 and 1
      call random_number(xx)
      ! Convert to normal distribution using the Box-Mueller algorithm
      norm = sqrt(-2._wp*log(xx(1))) * cos(2._wp*pi*xx(2))
      ! keep within 3 sigma
      norm = min(3._wp,norm)
      norm = max(-3._wp,norm)

      ! convert white noise to red noise following:
      ! https://atmos.washington.edu/~breth/classes/AM582/lect/lect8-notes.pdf
      norm = noise_autocorr*noise_old + sqrt(1._wp-noise_autocorr**2)*norm
      noise_old = norm

      noise = noise_amp_flx*norm

    else if (i_noise.eq.3) then
      ! periodic 'noise'

      noise = noise_amp_flx*cos(2._wp*pi*time/noise_period)

    endif

    noise_flx = noise

    ! update 2D freshwater noise
    do j=1,maxj
      do i=1,maxi
        if (mask_noise(i,j).eq.1) then
          if (.not.l_noise_only_coast .or. mask_coast(i,j).eq.1) then
            flx_noise(i,j) = noise  ! W/m2
          else
            flx_noise(i,j) = 0._wp 
          endif
        else
          flx_noise(i,j) = 0._wp 
        endif
      enddo
    enddo


   return

  end subroutine noise_flx_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  n o i s e _ i n i t
  !   Purpose    :  initialize noise
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine noise_init(f_ocn)

    implicit none

    real(wp), intent(in) :: f_ocn(:,:)

    integer :: i, j


    allocate(mask_noise(maxi,maxj))
    allocate(j_noise(maxi,2))

    ! find index of western and eastern boundary of box where noise should be applied
    do i=2,maxi
     ! western boundary
     if ((lon(i-1).le.lon_min_noise) .and. (lon(i).ge.lon_min_noise)) then
        i_noise_1 = i
     endif
     ! eastern boundary
     if ((lon(i-1).le.lon_max_noise) .and. (lon(i).ge.lon_max_noise)) then
        i_noise_2 = i-1
     endif
    enddo

    ! find index of northern and southern boundary of box where noise should be applied
    do j=2,maxj
     ! southern boundary
     if ((lat(j-1).le.lat_min_noise) .and. (lat(j).ge.lat_min_noise)) then
        j_noise(:,1) = j
     endif
     ! northern boundary
     if ((lat(j-1).le.lat_max_noise) .and. (lat(j).ge.lat_max_noise)) then
        j_noise(:,2) = j-1
     endif
    enddo

    ! derive noise mask
    mask_noise(:,:) = 0
    do i=i_noise_1,i_noise_2
      do j=j_noise(i,1),j_noise(i,2)
        if (noise_basin.eq.1) then
          ! Atlantic
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_atlantic) mask_noise(i,j) = 1
        else if (noise_basin.eq.2) then
          ! Pacific
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_pacific) mask_noise(i,j) = 1
        else if (noise_basin.eq.3) then
          ! Southern ocean
          if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_southern) mask_noise(i,j) = 1
        endif
      enddo
    enddo


    return

  end subroutine noise_init

end module noise_mod
