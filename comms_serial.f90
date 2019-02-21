! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: Git -*-
!=============================================================================!
!                             C O M M S                                       !
!=============================================================================!
! Abstracted comms layer for implementing parallel tempering Monte-Carlo      !
! This is the dummy (stub) serial version.                                    !
!=============================================================================!
module comms

  use constants 
  use model
  implicit none

  integer :: mpi_status_size = 1
  integer :: mpi_success = 0

  ! buffers for sending and recieving all molecule data during parallel
  ! tempering moves.
  real(kind=dp),allocatable,dimension(:),save :: Rbuffer,Sbuffer

  integer :: myrank,size

contains


  subroutine Comms_Initialise()   
    !=========================================================================!
    ! Initialises MPI and gets the rank of this task, and the size of the     !
    ! global communicator.                                                    !
    !-------------------------------------------------------------------------!
    ! David Quigley June 2003                                                 !
    !=========================================================================!
    implicit none



    return

  end subroutine Comms_Initialise

  subroutine comms_allocate()
    !=========================================================================!
    ! Allocates arrays used in parallel tempering send/receive moves          !
    ! Must be called after water data structures are created.                 !
    ! Only water array is sent
    !-------------------------------------------------------------------------!
    ! David Quigley September 2006                                            !
    !=========================================================================!
    use userparams, only : nwater
    implicit none


    return

  end subroutine comms_allocate


  subroutine Comms_BcastReal(Rvalue,Length)

    !=========================================================================!
    ! Routine to broadcast a single real value from the root node (task 0)    !
    ! to all other nodes.                                                     !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley June 2003                                      !
    !=========================================================================!
    implicit none
    real(kind=dp),intent(INOUT) :: Rvalue
    integer,intent(IN)          :: Length


 
    return

  end subroutine Comms_BcastReal

  subroutine Comms_BcastChar(Cvalue,Length)

    !=========================================================================!
    ! Routine to broadcast a single character from the root node (task 0)     !
    ! to all other nodes.                                                     !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley June 2003                                      !
    !=========================================================================!
    implicit none
    character,intent(inout) :: Cvalue
    integer,intent(in)      :: length

    return

  end subroutine Comms_BcastChar

  subroutine Comms_BcastLog(Lvalue,Length)

    !=========================================================================!
    ! Routine to broadcast a single logical   from the root node (task 0)     !
    ! to all other nodes.                                                     !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley June 2003                                      !
    !=========================================================================!
    implicit none
    logical,intent(inout) :: Lvalue
    integer,intent(in)    :: length 

    return

  end subroutine Comms_BcastLog

  subroutine comms_p2preal(Rvalue,length,snode,rnode)
    !=========================================================================!
    ! Sends real data from snode to rnode. Must be called by both snode and   !
    ! rnode.                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley June 2003                                      !
    !=========================================================================!
    implicit none
    real(kind=dp),intent(INOUT) :: Rvalue
    integer,intent(IN)          :: Length,snode,rnode


    return

  end subroutine comms_p2preal

  subroutine comms_p2pint(Ivalue,length,snode,rnode)
    !=========================================================================!
    ! Sends integer data from snode to rnode. Must be called by both snode    !
    ! rnode.                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley June 2003                                      !
    !=========================================================================!
    implicit none
    integer,intent(INOUT) :: Ivalue
    integer,intent(IN)    :: Length,snode,rnode

    return

  end subroutine comms_p2pint

  subroutine Comms_BcastInt(Ivalue,Length)

    !=========================================================================!
    ! Routine to broadcast a single real value from the root node (task 0)    !
    ! to all other nodes.                                                     !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley June 2003                                      !
    !=========================================================================!
    implicit none
    integer,intent(INOUT) :: Ivalue
    integer,intent(IN)    :: Length
   

    return

  end subroutine Comms_BcastInt

  subroutine comms_allreduce_eta(weight,Length)
    !=========================================================================!
    ! Wrapper for mpi_allreduce (real - double precision)                     !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley June 2003                                      !
    !=========================================================================!
    implicit none
    integer,intent(in) :: Length
    real(kind=dp),intent(inout),dimension(1:Length) :: weight

    return

  end subroutine comms_allreduce_eta

  subroutine comms_allreduce_hist(histogram,Length)
    !=========================================================================!
    ! Wrapper for mpi_allreduce (real - double precision)                     !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley June 2003                                      !
    !=========================================================================!
    implicit none
    integer,intent(in)   :: Length
    real(kind=dp),intent(inout),dimension(Length)   :: histogram


    return

  end subroutine comms_allreduce_hist

  subroutine comms_set_histogram(hist_in,length)
    !=========================================================================!
    ! Zero's out the internal array hist_last_sync. To be used when zeroing   !
    ! the histogram across all nodes.                                         !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley June 2003                                      !
    !=========================================================================!
    implicit none
    integer,intent(in) :: length
    real(kind=dp),dimension(length),intent(in) :: hist_in


    return

  end subroutine comms_set_histogram

  subroutine Comms_Finalise()

    !=========================================================================!
    ! Initialises MPI and gets the rank of this task, and the size of the     !
    ! global communicator.                                                    !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    ! MPI must have been initialised                                          !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley June 2003                                      !
    !=========================================================================!
    implicit none

    return

  end subroutine Comms_Finalise

  subroutine comms_barrier()
    !=========================================================================!
    ! Initialises MPI and gets the rank of this task, and the size of the     !
    ! global communicator.                                                    !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    ! MPI must have been initialised                                          !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley May 2011                                       !
    !=========================================================================!
    implicit none



    return

  end subroutine comms_barrier



end module comms
