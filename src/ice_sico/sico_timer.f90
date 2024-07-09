module sico_timer

  use sico_types_m
  use timer, only : sec_year

  implicit none

  type sico_timer_class
    real(wp) :: dtime
    real(wp) :: dtime_a
    real(wp) :: dtime_temp
    real(wp) :: dtime_mar_coa
    real(wp) :: time
    real(wp) :: time_init
    real(wp) :: dtime_inv
    real(wp) :: dtime_temp_inv
    integer :: itercount
    integer :: iter_temp
    integer :: n_discharge_call
    integer :: iter_mar_coa
  end type sico_timer_class

contains

  subroutine sico_timer_init(dtime, dtime_temp, dtime_mar_coa, timer)

    real(wp), intent(in) :: dtime
    real(wp), intent(in) :: dtime_temp
    real(wp), intent(in) :: dtime_mar_coa
    type(sico_timer_class) :: timer

    ! time step
    timer%dtime_a = dtime
    timer%dtime = dtime*sec_year
    timer%dtime_temp = dtime_temp*sec_year
    timer%dtime_mar_coa = dtime_mar_coa * sec_year ! a -> s

    timer%dtime_inv = 1.0_wp/timer%dtime
    timer%dtime_temp_inv = 1.0_wp/timer%dtime_temp

    timer%itercount = 0
    timer%iter_temp = nint(timer%dtime_temp/timer%dtime)

    timer%time = timer%time_init

    timer%n_discharge_call = -1
    timer%iter_mar_coa = nint(timer%dtime_mar_coa/timer%dtime)


    return

  end subroutine sico_timer_init

end module sico_timer

