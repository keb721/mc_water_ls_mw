! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: Git -*-
!=============================================================================!
!                            T I M E R                                        !
!=============================================================================!
! Contains routines to check time relative to a queue slot length             !
!=============================================================================!
module timer 

  use constants, only : dp
  implicit none

  private                ! Everything private

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
                                                      !... unless exposed here.
  public :: timer_init
  public :: timer_check_runtime
  public :: timer_elapsed_time

  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  public :: timer_qtime
  public :: timer_closetime
  real(kind=dp),save :: timer_qtime = 432000
  real(kind=dp),save :: timer_closetime = 3600

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  real(kind=dp),save :: last_time,current_time,elapsed_time
  integer,save       :: start_day
  logical,save       :: timer_initialised = .false.

contains
  
  subroutine timer_init()
    !-------------------------------------------------------------------------!
    ! Initialises the timer                                                   !
    !-------------------------------------------------------------------------!
    ! D.Quigley March 2010                                                    !
    !-------------------------------------------------------------------------!
    implicit none

    character(12) :: dat,tim,zon
    integer,dimension(8) :: info

    call date_and_time(dat,tim,zon,info)
    
    last_time = 3600_dp*real(info(5),kind=dp)+60.0_dp*real(info(6),kind=dp) &
               + real(info(7),kind=dp)+0.001_dp*real(info(8),kind=dp)

    start_day  = info(3)
    
!!$    write(*,*)
!!$    write(*,'("=====================")')
!!$    write(*,'("| Timer initialised |")')              
!!$    write(*,'("=====================")')
!!$    write(*,*)

    elapsed_time = 0.0_dp
    timer_initialised = .true.

    return

  end subroutine timer_init

  real(kind=dp) function timer_elapsed_time()
    !-------------------------------------------------------------------------!
    ! Returns time since timer was initialised                                !
    !-------------------------------------------------------------------------!
    ! D.Quigley March 2010                                                    !
    !-------------------------------------------------------------------------!
    implicit none

    character(12) :: dat,tim,zon
    integer,dimension(8) :: info

    if (.not.timer_initialised) then
       stop 'Called timer_elasped time before initialising timer!'
    end if

    call date_and_time(dat,tim,zon,info)
    
    current_time = 3600_dp*real(info(5),kind=dp)+60.0_dp*real(info(6),kind=dp) &
                   + real(info(7),kind=dp)+0.001_dp*real(info(8),kind=dp)
    
    ! if the day has changed....
    if (start_day/=info(3)) then
        last_time = last_time - 86400_dp
       start_day  = info(3)
    end if

    ! Compute the elapsed time
    elapsed_time = elapsed_time + current_time - last_time
    last_time    = current_time
    timer_elapsed_time = elapsed_time
    
    return

  end function timer_elapsed_time

  subroutine timer_check_runtime(safe)

    implicit none
    logical,intent(out) :: safe
    real(kind=dp)       :: tmptime
    logical,save        :: lwarn = .false.

    safe = .false.

    tmptime = timer_elapsed_time()

    if (timer_qtime - tmptime>timer_closetime) then
       safe = .true.
    elseif (.not.lwarn) then
       write(0,'("! ====================================================================== ! ")')
       write(0,'("! WARNING less than timer_closetime remaining before timer_qtime reached ! ")')
       write(0,'("! ====================================================================== ! ")')
       lwarn = .true.
    end if

    return

  end subroutine timer_check_runtime


end module timer
