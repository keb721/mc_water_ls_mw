! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: Git -*-
!=============================================================================!
!                             C O M M S                                       !
!=============================================================================!
! Abstracted comms layer for implementing parallel tempering Monte-Carlo      !
! This is the MPI version.                                                    !
!=============================================================================!
module comms

  use constants, only : dp
  implicit none
  include 'mpif.h'  


  ! buffers for sending and recieving all molecule data during parallel
  ! tempering moves.
  real(kind=dp),allocatable,dimension(:),save :: Rbuffer,Sbuffer
  real(kind=dp),allocatable,dimension(:),save :: hist_last_sync,eta_last_sync
  real(kind=dp),allocatable,dimension(:),save :: uhist_last_sync
  

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

    integer :: ierr !,provided

    ! Initialise MPI
    call mpi_init(ierr)


!!$    ! Initialise MPI - checking thread support
!!$    call mpi_init_thread(MPI_THREAD_FUNNELED,provided,ierr)
!!$    
!!$    if (provided==MPI_THREAD_FUNNELED) then
!!$       write(0,'("MPI provided MPI_THREAD_FUNNELED : OK")')
!!$    else
!!$       write(0,'("-------- W A R N I N G ------------------")')
!!$       write(0,'(" MPI did not provide MPI_THREAD_FUNNELED ")')
!!$       if (provided==MPI_THREAD_SINGLE) then
!!$          write(0,'(" MPI provided MPI_THREAD_SINGLE")')
!!$       else if (provided==MPI_THREAD_SERIALIZED) then
!!$          write(0,'(" MPI provided MPI_THREAD_SERIALIZED")')
!!$       else if (provided==MPI_THREAD_MULTIPLE) then
!!$          write(0,'(" MPI provided MPI_THREAD_MULTIPLE")')
!!$       end if
!!$       write(0,'("-----------------------------------------")')
!!$    end if
!!$
!!$    if (ierr.ne.MPI_SUCCESS) stop 'Error initialising MPI - exiting!'

    ! Get the size of the communicator
    call mpi_comm_size(MPI_COMM_WORLD,size,ierr)
    if (ierr.ne.MPI_SUCCESS) stop 'Could not get size of communicator'

    ! Get my rank in the communicator
    call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)
    if (ierr.ne.MPI_SUCCESS) stop 'Could not get rank'

    return

  end subroutine Comms_Initialise

  subroutine comms_allocate()
    !=========================================================================!
    ! Allocates arrays used in parallel tempering send/receive moves          !
    ! Must be called after water data structures are created.                 !
    !-------------------------------------------------------------------------!
    ! David Quigley September 2006                                            !
    !=========================================================================!
    use userparams, only : nwater,mc_ensemble,nbins,samplerun
    use model,      only : cspm,ljspm
    implicit none
    integer :: dim,ierr

    ! 3 reals for each site + 4 for quaternions
    dim = nwater* (3*cspm + 3*ljspm + 4)
    if ( mc_ensemble == 'npt' ) dim = dim + 9
    allocate(Rbuffer(1:dim),stat=ierr)
    allocate(Sbuffer(1:dim),stat=ierr)
    if (ierr/=0) stop 'Error allocating send and recieve buffers'
    
    allocate(eta_last_sync(1:nbins),stat=ierr)
    allocate(hist_last_sync(1:nbins),stat=ierr)
    if (samplerun) allocate(uhist_last_sync(1:nbins),stat=ierr)

    if (ierr/=0) stop 'Error allocating bin sync arrays'

    eta_last_sync = 0.0_dp
    hist_last_sync = 0.0_dp
    if (samplerun) uhist_last_sync = 0.0_dp

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
    integer :: ierr

    call MPI_BCAST(Rvalue,Length,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if (ierr.ne.MPI_SUCCESS) stop 'Error in Comms_BcastReal - comms operation failed'


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
    integer :: ierr

    call MPI_BCAST(Cvalue,Length,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    if (ierr.ne.MPI_SUCCESS) stop 'Error in Comms_BcastChar - comms operation failed'

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
    integer :: ierr

    call MPI_BCAST(Lvalue,Length,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (ierr.ne.MPI_SUCCESS) stop 'Error in Comms_BcastLog - comms operation failed'

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
    integer :: ierr

    integer,dimension(1:MPI_STATUS_SIZE) :: status


    if (myrank==snode) then
       call MPI_Send(Rvalue,length,MPI_DOUBLE_PRECISION,rnode,25,MPI_COMM_WORLD,ierr)
       if (ierr/=0) stop 'Error during send'
    elseif(myrank==rnode) then
       call MPI_Recv(Rvalue,length,MPI_DOUBLE_PRECISION,snode,25,MPI_COMM_WORLD,status,ierr)
       if (ierr/=0) stop 'Error during receive'
    end if


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
    integer :: ierr

    integer,dimension(1:MPI_STATUS_SIZE) :: status


    if (myrank==snode) then
       call MPI_Send(Ivalue,length,MPI_DOUBLE_PRECISION,rnode,25,MPI_COMM_WORLD,ierr)
       if (ierr/=0) stop 'Error during send'
    elseif(myrank==rnode) then
       call MPI_Recv(Ivalue,length,MPI_DOUBLE_PRECISION,snode,25,MPI_COMM_WORLD,status,ierr)
       if (ierr/=0) stop 'Error during receive'
    end if


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
    integer :: ierr

    call MPI_BCAST(Ivalue,Length,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if (ierr.ne.MPI_SUCCESS) stop 'Error in Comms_BcastInt - comms operation failed'

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
    real(kind=dp),dimension(:),allocatable :: buff
    integer :: ierr

    weight = weight - eta_last_sync

    ! sanity check
!!$    do i = 1,2
!!$       do j = 1,length
!!$          if (weight(j,i)<0.0_dp) stop 'Negative growth of eta since last sync'
!!$       end do
!!$    end do

    allocate(buff(1:length))
    call MPI_ALLreduce(weight(1),buff(1),length,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    if (ierr.ne.MPI_SUCCESS) stop 'Error in comms_allreducereal - comms operation failed'

    weight = buff + eta_last_sync
    eta_last_sync = weight


    deallocate(buff)

    return

  end subroutine comms_allreduce_eta

  subroutine comms_get_max(myval,maxval)
    !=========================================================================!
    ! Find the maximum of a double over all MPI ranks.                        !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley February 2019                                  !
    !=========================================================================!
    implicit none
    real(kind=dp),intent(in)  :: myval
    real(kind=dp),intent(out) :: maxval
    
    integer :: ierr

    call MPI_Allreduce(myval,maxval,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (ierr/=0) stop 'Error finding maximum value over all nodes'

    return


  end subroutine comms_get_max

subroutine comms_join_uhist(uhist,length,overlap,joined)
    !=========================================================================!
    ! Join unbiased histogram together from multiple overlapping windows.     !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley February 2019                                  !
    !=========================================================================!

    implicit none
    integer,intent(in) :: length,overlap
    real(kind=dp),intent(in),dimension(1:length) :: uhist
    real(kind=dp),intent(out),dimension(1:length) :: joined

    real(kind=dp),dimension(1:length) :: recvbuff

    real(kind=dp) :: shift,myave,nextav

    integer,dimension(MPI_STATUS_SIZE) :: status
    
    integer :: bins_per_window,irank,ierr,k,my_start_bin,my_end_bin

    ! Clunky - replace with something useing MBAR?

    bins_per_window = length/size 
    
    if (myrank==0) joined = uhist

    do irank = 1,size-1

       if (myrank == 0) then

          ! Receive from irank at every iteration
          call MPI_Recv(recvbuff,length,MPI_DOUBLE_PRECISION,irank,99,MPI_COMM_WORLD,status,ierr)
          if (ierr/=0) stop 'Error receiving from partner in comms_join_uhist'


       elseif (myrank==irank) then
       
          ! Send to rank zero if iteration counter matches current rank
          call MPI_Send(uhist,length,MPI_DOUBLE_PRECISION,0,99,MPI_COMM_WORLD,ierr)
          if (ierr/=0) stop 'Error sending to partner in comms_join_uhist'

       end if

       if (myrank==0) then


          ! Range of bins rank 0 already has data for
          my_start_bin = 1
          my_end_bin = irank*bins_per_window 

          myave = 0.0_dp
          do k = my_end_bin-overlap,my_end_bin+overlap
             myave = myave + log(joined(k))
          end do
          myave = myave/real(2*overlap+1,kind=dp)

          ! Average from the overlapping bins in the next window
          nextav = 0.0_dp
          do k = my_end_bin-overlap,my_end_bin+overlap
             nextav = nextav + log(recvbuff(k))
          end do
          nextav = nextav/real(2*overlap+1,kind=dp)

          shift = myave-nextav
          if (isnan(shift)) shift = 0.0_dp

          do k = my_end_bin+1,length
             joined(k) = recvbuff(k)*exp(shift)
          end do
         
       end if

    end do ! end loop over ranks

    ! Broadcast joined to everyone
    call MPI_Bcast(joined,length,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if (ierr/=0) stop 'Error broadcasting joined unbiased histogram'

    return

  end subroutine comms_join_uhist

  subroutine comms_join_eta(weight,Length,overlap,joined)
    !=========================================================================!
    ! Join  weight functions together from multiple overlapping windows.      !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley February 2019                                  !
    !=========================================================================!
    implicit none
    integer,intent(in) :: length,overlap
    real(kind=dp),intent(in),dimension(1:length) :: weight
    real(kind=dp),intent(out),dimension(1:length) :: joined

    real(kind=dp),dimension(1:length) :: recvbuff

    real(kind=dp) :: shift,myave,nextav

    integer,dimension(MPI_STATUS_SIZE) :: status
    
    integer :: bins_per_window,irank,ierr,k,my_start_bin,my_end_bin

    ! Clunky - replace with something useing MBAR?

    bins_per_window = length/size 
    
    if (myrank==0) joined = weight

    do irank = 1,size-1

       if (myrank == 0) then

          ! Receive from irank at every iteration
          call MPI_Recv(recvbuff,length,MPI_DOUBLE_PRECISION,irank,99,MPI_COMM_WORLD,status,ierr)
          if (ierr/=0) stop 'Error receiving from partner in comms_join_eta'


       elseif (myrank==irank) then
       
          ! Send to rank zero if iteration counter matches current rank
          call MPI_Send(weight,length,MPI_DOUBLE_PRECISION,0,99,MPI_COMM_WORLD,ierr)
          if (ierr/=0) stop 'Error sending to partner in comms_join_eta'

       end if

       if (myrank==0) then


          ! Range of bins rank 0 already has data for
          my_start_bin = 1
          my_end_bin = irank*bins_per_window 

          myave = sum(joined(my_end_bin-overlap:my_end_bin+overlap))
          myave = myave/real(2*overlap+1,kind=dp)

          ! Average from the overlapping bins in the next window
          nextav = sum(recvbuff(my_end_bin-overlap:my_end_bin+overlap))
          nextav = nextav/real(2*overlap+1,kind=dp)

          shift = myave-nextav

          do k = my_end_bin+1,length
             joined(k) = recvbuff(k) + shift
          end do
         
       end if

    end do ! end loop over ranks

    shift = joined(length/2+1)
    do k = 1,length
       joined(k) = joined(k)-shift
    end do

    ! Broadcast joined to everyone
    call MPI_Bcast(joined,length,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if (ierr/=0) stop 'Error broadcasting joined weight function'


    return

  end subroutine comms_join_eta

  subroutine comms_allreduce_hist(histogram,Length)
    !=========================================================================!
    ! Wrapper for mpi_allreduce (real - double precision)                     !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley June 2003                                      !
    !=========================================================================!
    implicit none
    integer,intent(in)   :: Length
    real(kind=dp),intent(inout),dimension(Length)   :: histogram
    real(kind=dp),dimension(:),allocatable :: buff
    integer :: ierr

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    histogram = histogram - hist_last_sync

    ! sanity check
!!$    do j = 1,length
!!$       if (histogram(j)<0) stop 'Negative growth of histogram since last sync'
!!$    end do


    allocate(buff(1:length))
    call MPI_ALLreduce(histogram,buff,length,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    if (ierr.ne.MPI_SUCCESS) stop 'Error in comms_allreducereal - comms operation failed'

    histogram = buff + hist_last_sync
    hist_last_sync = histogram

    deallocate(buff)

    return

  end subroutine comms_allreduce_hist

  subroutine comms_allreduce_uhist(histogram,Length)
    !=========================================================================!
    ! Wrapper for mpi_allreduce (real - double precision)                     !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley June 2003                                      !
    ! Unbiased histogram version February 2019                                !
    !=========================================================================!
    implicit none
    integer,intent(in)   :: Length
    real(kind=dp),intent(inout),dimension(Length)   :: histogram
    real(kind=dp),dimension(:),allocatable :: buff
    integer :: ierr

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    histogram = histogram - uhist_last_sync

    ! sanity check
!!$    do j = 1,length
!!$       if (histogram(j)<0) stop 'Negative growth of histogram since last sync'
!!$    end do


    allocate(buff(1:length))
    call MPI_ALLreduce(histogram,buff,length,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    if (ierr.ne.MPI_SUCCESS) stop 'Error in comms_allreducereal - comms operation failed'

    histogram = buff + uhist_last_sync
    uhist_last_sync = histogram

    deallocate(buff)

    return

  end subroutine comms_allreduce_uhist


  subroutine comms_set_histogram(hist_in,length)
    !=========================================================================!
    ! Sets the internal array hist_last_sync. To be used when zeroing         ! 
    ! the histogram across all nodes.                                         !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley June 2003                                      !
    !=========================================================================!
    implicit none
    integer,intent(in) :: length
    real(kind=dp),dimension(length),intent(in) :: hist_in

    hist_last_sync = hist_in

    return

  end subroutine comms_set_histogram

  subroutine comms_set_uhistogram(hist_in,length)
    !=========================================================================!
    ! Sets the internal array uhist_last_sync. To be used when zeroing        !
    ! the histogram across all nodes.                                         !
    !-------------------------------------------------------------------------!
    ! Written by David Quigley June 2003                                      !
    ! Unbiased histogram version February 2019                                !
    !=========================================================================!
    implicit none
    integer,intent(in) :: length
    real(kind=dp),dimension(length),intent(in) :: hist_in

    uhist_last_sync = hist_in

    return

  end subroutine comms_set_uhistogram


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
    use userparams, only : samplerun
    implicit none
    integer :: ierr

    !Finalise comms
    call MPI_FINALIZE(ierr)
    if (ierr.ne.MPI_SUCCESS) stop 'Error finalising MPI!'

    deallocate(RBuffer,stat=ierr)
    deallocate(SBuffer,stat=ierr)
    if (ierr/=0) stop 'Error deallocating send and recieve buffer'

    deallocate(hist_last_sync,stat=ierr)
    deallocate(eta_last_sync,stat=ierr)
    if (samplerun) deallocate(uhist_last_sync,stat=ierr)
    if (ierr/=0) stop 'Error deallocating bin arrays'

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
    integer :: ierr

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    return

  end subroutine comms_barrier

end module comms
