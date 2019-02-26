! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: Git -*-
!=============================================================================!
!                             M C _ M O V E S                                 !
!=============================================================================!
! Contains to implement various Monte-Carlo trial moves and accept/reject     !
! based on the change in energy.                                              !
!=============================================================================!
module mc_moves


  use constants, only : dp

  implicit none                                       !Impose strong typing
  private   

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
                                                      !... unless exposed here.

  public :: mc_init
  public :: mc_deinit
  public :: mc_cycle
  public :: mc_checkpoint_write

  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  public :: ls_mu
  public :: ref_enthalpy

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!


  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  integer,save :: mc_cycle_num = 0

  ! Running counters of moves acceptance 
  integer,save :: mc_accepted_rsteps = 0
  integer,save :: mc_accepted_vsteps = 0
  integer,save :: mc_accepted_swtch  = 0

  ! Running counters of move attempts
  integer,save :: mc_attempted_rsteps = 0
  integer,save :: mc_attempted_vsteps = 0
  integer,save :: mc_attempted_swtch  = 0

  ! Per-molecule counters of move attempts                                                                                 
  integer,allocatable,dimension(:),save :: mc_translations

  ! allocated in mc_init - backups in case of mc move rejection
  real(kind=dp),dimension(2),save   :: backup_model_energy
  
  real(kind=dp),allocatable,dimension(:,:),save :: backup_local_real_energy

  ! lattice-switch order parameter
  real(kind=dp),save :: ls_mu

  real(kind=dp),allocatable,dimension(:),save   :: histogram
  real(kind=dp),allocatable,dimension(:),save   :: weight
  real(kind=dp),allocatable,dimension(:),save   :: unbiased_hist
  real(kind=dp),allocatable,dimension(:),save   :: mu_bin
  real(kind=dp),allocatable,dimension(:),save   :: binwidth
  real(kind=dp) :: av_binwidth
  real(kind=dp) :: log_unbiased_norm
  
  
  ! variables controlling the non-uniform bin-grid
  real(kind=dp),save :: a_pos,a_neg  ! Initial value of GPs
  real(kind=dp),save :: r_pos,r_neg  ! Common ratios of GPs
  real(kind=dp),save :: s_pos,s_neg  ! GP sums 

  ! modification factor as read from the input file
  real(kind=dp),save :: orig_wl_factor
  logical      ,save :: wl_invt_active = .false. 


  integer,save :: startcycle = 0      ! zero if not restart
  logical,save :: firstcycle = .true. ! is this the original wl_factor

  ! Energy checks
  real(kind=dp),dimension(2),save :: ref_enthalpy,average_energy=0.0_dp

  ! Checks for range of overlap parameter change in a single move
  real(kind=dp),save :: max_dmu=0.0_dp,min_dmu=huge(1.0_dp),dmu 
  
  ! Extra variables for Swetnam-style Wang-Landau
  real(kind=dp),save:: sumhist = 0.0

  ! Unit numberes for dumping of histograms etc
  integer,parameter :: ubh = 61    ! unbiased histogram
  integer,parameter :: his = 62    ! raw histogram
  integer,parameter :: wgt = 63    ! multicanonical weights
  integer,parameter :: wlf = 64    ! record of f values
  integer,parameter :: chk = 65    ! checkpoint file
  

  ! Lookup list for "other" lattice
  integer,dimension(2) :: partner_lattice

  ! Local (per MPI rank) domain variables
  real(kind=dp),save :: my_mu_max, my_mu_min
  integer,save :: my_start_bin, my_end_bin

  ! Has the walker on this MPI rank reached its window
  logical,save :: walker_in_window = .false.


contains

  subroutine mc_cycle()
    !------------------------------------------------------------------------------!
    ! Performs an MC cycle and updates parameters toward target acceptance ratio   !
    ! if during equilibration phase and is samplerun is not set true.              !
    !------------------------------------------------------------------------------!
    ! D.Quigley January 2007                                                       !
    !           Februray 2019 : Updated for domain decomposition option.           !
    !------------------------------------------------------------------------------!
    use comms,      only : comms_allreduce_eta,comms_allreduce_hist, &
                           comms_allreduce_uhist,myrank,size
    use random,     only : random_uniform_random
    use userparams, only : nwater,mc_ensemble,nbins,allow_switch,  &
         allow_vol,allow_trans,pressure,num_lattices, &
         mc_always_switch,mc_trans_prob,mc_vol_prob,mc_switch_prob, &
         deltaG_int,chkpt_dump_int,list_update_int, &
         flat_chk_int,latt_sync_int,monitor_int,mpi_sync_int, &
         parallel_strategy, eq_mc_cycles, samplerun
    use model,      only : volume
    use util,       only : util_determinant,util_recipmatrix
    use energy

    implicit none

    real(kind=dp),save :: sum_prob,transP,swP,volP
    real(kind=dp)      :: xi
    integer            :: imove,ils,irank
    logical            :: in_window
    logical,save       :: firstpass = .true.


    ! Increment MC cycle number
    mc_cycle_num = mc_cycle_num + 1

    !-------------------------------------------------------!
    ! Compute probability for performing each MC move type. !
    !-------------------------------------------------------!
    if (firstpass) then

       firstpass = .false.

       ! Compute move probabilities
       if (mc_always_switch)   mc_switch_prob = 0.0_dp
       if (.not.allow_switch)  mc_switch_prob = 0.0_dp
       if (mc_ensemble=='nvt') mc_vol_prob    = 0.0_dp
       if (.not.allow_vol)     mc_vol_prob    = 0.0_dp
       if (.not.allow_trans)   mc_trans_prob  = 0.0_dp

       sum_prob = mc_trans_prob + mc_vol_prob + mc_switch_prob

       transP = mc_trans_prob/sum_prob
       volP   = mc_vol_prob/sum_prob
       swP    = mc_switch_prob/sum_prob

       ! Store as cumulative probability
       volP   = volP  + transP
       swP    = swP   + volP

       if (swP < 0.999_dp ) stop 'Cumulative move type probability error'

    end if

    !-------------------------------------!
    ! Sanity check on equilibration       !
    !-------------------------------------!
    if (parallel_strategy == 'dd') then      

       if (mc_cycle_num < eq_mc_cycles) then
          
          walker_in_window = ( ( ls_mu > my_mu_min ) .and. ( ls_mu < my_mu_max ) )

       else if (mc_cycle_num == eq_mc_cycles) then

          ! Check that the walker on this MPI task has made its way into 
          ! its window. Otherwise we've got a problem and need to abort.
          if ( .not. walker_in_window) then
             write(0,'("Error : Not all walkers have reached their designated window ")')
             write(0,'("        after ",I10," MC cycles.")')eq_mc_cycles
             write(0,*)
             do irank = 0,size-1
                if (myrank==irank) then
                   write(0,'("Rank ",I2, " has overlap parameter mu = ",F12.6)')irank,ls_mu
                end if
             end do
             stop
          end if

       else

          ! This must be a restart
          walker_in_window = .true.

       end if

    end if


    !-------------------------------------!
    ! Perform nwater trial moves          !
    !-------------------------------------!

    ! Neighbour list update
    if ( mod(mc_cycle_num,list_update_int)==0 ) then
       do ils = 1,num_lattices
          call compute_neighbours(ils)
       end do
    end if
 
    do imove = 1,nwater

       xi = random_uniform_random()

       if (xi<transP) then
          call mc_water_translation()
          call mc_update_wl_bins()
          mc_attempted_rsteps = mc_attempted_rsteps + 1       
       elseif (xi<volP) then
          call mc_volume()
          call mc_update_wl_bins()
          mc_attempted_vsteps = mc_attempted_vsteps + 1
       elseif (xi<swP) then
          if (.not.((parallel_strategy=='dd').and.(mc_cycle_num<eq_mc_cycles))) then
             call mc_lattice_switch()
             mc_attempted_swtch = mc_attempted_swtch + 1
          end if
       end if
       
       if (mc_always_switch) then
          if (.not.((parallel_strategy=='dd').and.(mc_cycle_num<eq_mc_cycles))) then
             call mc_lattice_switch()
             mc_attempted_swtch = mc_attempted_swtch + 1
          end if
       end if

    end do

    ! accumulate the average energy in the current sample block
    average_energy(1:num_lattices) = average_energy(1:num_lattices) + model_energy(1:num_lattices)
    if ( mc_ensemble == 'npt' ) & 
         average_energy(1:num_lattices) = average_energy(1:num_lattices) + pressure*volume(1:num_lattices)

    ! If doing lattice switching with multiple walkers, sychronise weights and histograms
    if ((num_lattices == 2)) then

       ! MPI - synchronise histogram and weights at this interval
       if ( mod(mc_cycle_num,mpi_sync_int)==0 ) then
          !       write(mylog,'(" DEBUG - syncing weights and histograms at cycle ",I10)')mc_cycle_num

          if (parallel_strategy=='mw') then
             ! Construct global histogram as sum over individual walker histograms
             call comms_allreduce_eta(weight,nbins)
             call comms_allreduce_hist(histogram,nbins)
             if (samplerun) call comms_allreduce_uhist(unbiased_hist,nbins)

          elseif (parallel_strategy=='dd') then
             ! Build global histogram from individual windows ?

          endif


       end if
       
    end if ! end if lattice-switching
    
    ! Adjust move sizes, report acceptance rates and 
    ! output histogram and weight functions.
    if ( mod(mc_cycle_num,monitor_int)==0 ) then
!       write(mylog,'(" DEBUG - monitoring stats at cycle ",I10)')mc_cycle_num
       call mc_monitor_stats() 
    end if

    if (num_lattices == 2) then

       ! Check for flatness of histogram and take action
       ! accordingly.
       if ( mod(mc_cycle_num,flat_chk_int)==0 ) then
          !       write(mylog,'(" DEBUG - checking flatness  at cycle ",I10)')mc_cycle_num
          call mc_check_flatness()
       end if
       
       ! Enforce synchoronisation of the two Markov chains
       if (mod(mc_cycle_num,latt_sync_int)==0 ) then
          !       write(mylog,'(" DEBUG - checking synchronisation  at cycle ",I10)')mc_cycle_num
          call mc_check_chain_synchronisation()
       end if
       
       if (mod(mc_cycle_num,deltaG_int)==0) then
          !       write(mylog,'(" DEBUG - computing delta G at cycle ",I10)')mc_cycle_num
          !call mc_compute_deltaG_from_eta()
          if (samplerun) call mc_compute_deltaG_from_hist()
       end if

    end if

    !--------------------------------------------------------!
    ! Checkpointing                                          !
    !--------------------------------------------------------!
    if ( mod(mc_cycle_num,chkpt_dump_int)==0 ) then
!       write(mylog,'(" DEBUG - writing checkpoint at cycle ",I10)')mc_cycle_num
       call mc_checkpoint_write()
    end if


    return

  end subroutine mc_cycle


  subroutine mc_checkpoint_write()
    !------------------------------------------------------------------------------!
    ! writes a checkpoint file used for restarting the calculation from the        !
    ! current point                                                                !
    !------------------------------------------------------------------------------!
    ! D.Quigley January 2007                                                       !
    !------------------------------------------------------------------------------!
    use comms,      only : myrank,comms_barrier
    use userparams, only : nbins,mc_max_trans,wl_factor,nwater,mc_dv_max,samplerun
    use model,      only : hmatrix,ljr,ref_ljr,ls
    implicit none
    integer :: i,ierr
    integer,save :: fn = 1
    character(1) :: fns
    character(3) :: myrstr
    
    write(fns,'(I1)')fn
    write(myrstr,'(I3.3)')myrank
    open(unit=chk,file='checkpoint'//myrstr//'.dat.'//fns,status='replace', &
         form='unformatted',iostat=ierr)
    if ( ierr/=0 ) stop 'Error opening checkpoint.dat'

    if (fn==1) then
       fn = 2
    else
       fn = 1
    end if

    ! write number of units
    write(chk)nwater

    ! current mc cycle
    write(chk)mc_cycle_num

    ! write move generation data
    write(chk)mc_max_trans,mc_dv_max
    
    ! write wang landau data
    write(chk)wl_factor
    write(chk)(histogram(i),i=1,nbins)
    write(chk)(weight(i),i=1,nbins)
    write(chk)wl_invt_active

    if (samplerun) write(chk)(unbiased_hist(i),i=1,nbins)

    ! write cell data
    write(chk)hmatrix

    ! write reference positions
    write(chk)ref_ljr

    ! write positions
    write(chk)ljr

    ! write the currently active lattice
    write(chk)ls

    close(chk)

    ! Make sure everyone has caught up otherwise may end up with the guy who
    ! finishes first here writing another checkpoint before the slowest guy
    ! has written this one.
    call comms_barrier()

    return

  end subroutine mc_checkpoint_write


  subroutine mc_checkpoint_load()
    !------------------------------------------------------------------------------!
    ! Reads a checkpoint file used for restarting the calculation  from the        !
    !------------------------------------------------------------------------------!
    ! D.Quigley January 2007                                                       !
    !------------------------------------------------------------------------------!
    use comms,      only : myrank,comms_set_histogram,comms_bcastint, &
                           comms_allreduce_uhist,comms_set_uhistogram
    use io,         only : mylog
    use userparams, only : nbins,mc_max_trans,wl_factor,mc_dv_max,samplerun, &
                           parallel_strategy
    use model,      only : hmatrix,ljr,ref_ljr,ls

    implicit none

    integer       :: i
    real(kind=dp) :: tmp_wl_factor
    integer       :: fn,fnu
    integer,dimension(2) :: ferr,fcyc,tmpint
    character(1)  :: fs
    character(3)  :: rankstr

    write(rankstr,'(I3.3)')myrank

    do fn = 1,2
       
       write(fs,'(I1)')fn
       open(unit=chk+fn,file='checkpoint'//rankstr//'.dat.'//fs,status='old',form='unformatted',iostat=ferr(fn))

       ! read number of units
       read(chk+fn,iostat=ferr(fn))tmpint(fn)

       ! read the current MC cycle number
       read(chk+fn,iostat=ferr(fn))fcyc(fn)

    end do

    ! Decide which checkpoint file to use
    if ((ferr(1)/=0).and.(ferr(2)/=0)) then
       ! Both useless - give up
       write(mylog,'("Error on rank ",I5," no valid checkpoint file.")')myrank
       stop 'Error using checkpoint files - check logs'
    else if (ferr(1)/=0) then
       ! Use 2
       fnu = 2
    else if (ferr(2)/=0) then
       ! Use 1
       fnu = 1
    else
       ! Both are useful - use most recent
       if (fcyc(1)>fcyc(2)) then
          fnu = 1
       else
          fnu = 2
       end if
    end if

    startcycle = fcyc(fnu)
    mc_cycle_num = startcycle

    ! Synchronise
    call comms_bcastint(mc_cycle_num,1)

    ! read move generation data
    read(chk+fnu)mc_max_trans,mc_dv_max
    
    ! read Wang-Landau data
    tmp_wl_factor = wl_factor
    read(chk+fnu)wl_factor
    read(chk+fnu)(histogram(i),i=1,nbins)
    read(chk+fnu)(weight(i),i=1,nbins)
    read(chk+fnu)wl_invt_active

    if (samplerun) read(chk+fnu)(unbiased_hist(i),i=1,nbins)

    ! Make sure the comms module knows what's going on
    if (parallel_strategy=='mw') then
       call comms_set_histogram(histogram,nbins)
       if (samplerun) call comms_set_uhistogram(unbiased_hist,nbins)
    end if

    ! compute sum over all histogram bins
    sumhist = real(sum(histogram),kind=dp)

    ! check for updates of wl_factor
    if ( wl_factor < orig_wl_factor ) then
       !write(*,*)"wl_factor smaller than original : firstcycle = false"
       firstcycle = .false.
    end if

    ! read cell data
    read(chk+fnu)hmatrix

    ! read reference positions
    read(chk+fnu)ref_ljr

    ! read positions
    read(chk+fnu)ljr

    ! read current lattice if present
    read(chk+fnu,end=10)ls


    10 close(chk+1,iostat=ferr(1))
       close(chk+2,iostat=ferr(2))
     
    return

  end subroutine mc_checkpoint_load


  subroutine mc_init(outcycle)
    !------------------------------------------------------------------------------!
    ! Computes a trial translation of a single molecule and computes the resulting !
    ! change in energy. The move is accepted or rejected appropriately.            !
    !------------------------------------------------------------------------------!
    ! D.Quigley January 2007                                                       !
    !           Feburary 2019 : Updated for domain decomposition.                  !
    !------------------------------------------------------------------------------!
    use constants, only : Kb
    use comms,     only : myrank,comms_bcastreal,comms_allreduce_eta,size
    use io,        only : glog,mylog
    use model,     only : volume,hmatrix,recip_matrix,ls
    use energy,    only : model_energy,compute_model_energy,compute_ivects
    use userparams,only : temperature,pressure,nbins,wl_factor, &
                          nwater,mu_max,mu_min,num_lattices,leshift, &
                          parallel_strategy, window_overlap,max_mc_cycles, &
                          eq_mc_cycles,samplerun
    use util,      only : util_determinant,util_recipmatrix
    implicit none
    integer,intent(out) :: outcycle ! pass first cycle num back up to caller

    integer :: ierr,k,ils,Ns,ibin,irank,bins_per_window
    logical :: lexist1,lexist2,lexist
    character(3)  :: rankstring
    character(29) :: dumchar29
    real(kind=dp) :: tmp_wl_factor,tmpsum,r_new,mu_u,mu_l,dumreal,beta,dumreal2
    real(kind=dp) :: hits_per_bin,unbiased_norm,incr

    write(mylog,'("#                                                              #")')
    write(mylog,'("# Initialising Monte-Carlo calculation                         #")')
    write(mylog,'("# ------------------------------------                         #")')
    write(mylog,'("#                                                              #")')

    !----------------------------------------!
    ! Allocate arrays for resetting moves    !
    !----------------------------------------!
    allocate(backup_local_real_energy(1:nwater,1:2),stat=ierr)
    if (ierr/=0) stop 'Error allocating backup local real energies'

    !-------------------------------------------!
    ! Arrays for monitoring moves per molecule  !
    !-------------------------------------------!
    allocate(mc_translations(1:nwater),stat=ierr)
    if (ierr/=0) stop 'Error allocating move counter arrays'
    mc_translations = 0


    !--------------------------------------------------!
    ! Set up a non-uniform grid of bins using a GP     !
    !--------------------------------------------------!

    ! First allocate bin arrays, making sure that nbins
    ! is odd
    if (mod(nbins,2)==0) nbins = nbins + 1

    allocate(mu_bin(1:nbins)       ,stat=ierr)
    allocate(binwidth(1:nbins)     ,stat=ierr)
    allocate(weight(1:nbins)       ,stat=ierr)  
    allocate(histogram(1:nbins)    ,stat=ierr)
    allocate(unbiased_hist(1:nbins),stat=ierr)
    weight = 0.0_dp
    histogram = 0.0_dp

    ! Unbiased_hist 
    unbiased_hist = 0.0_dp

    ! Initialise geometric series of grid spacings
    s_pos = abs(mu_max) - 0.5_dp
    s_neg = abs(mu_min) - 0.5_dp

    a_pos = 1.0_dp ! First positive bin width
    a_neg = 1.0_dp ! First negative bin width

    ! nbins/2 bin spaces between mu=+0.5 and mu=mu_max
    ! nbins/2 bin spaces between mu=-0.5 and mu=mu_min
    Ns = nbins/2

!!$    write(0,'("Rank ",I5," has Ns = ",I5)')myrank,Ns

    ! Compute r_pos - common ratio of positive series
    r_pos = 1.1_dp
    k = 0
    do
       k = k + 1
       tmpsum = a_pos*(1.0_dp-r_pos**Ns)/(1.0_dp-r_pos)
       r_new = r_pos*(s_pos/tmpsum)**(1.0_dp/real(Ns,kind=dp))
       if (abs(r_new-r_pos)<=2.0_dp*epsilon(1.0_dp)) exit
       if (k>1000000) exit
!!$       write(0,'(3E15.6)')abs(r_new-r_pos),tiny(1.0_dp),epsilon(1.0_dp)
       r_pos = r_new  ! converged
    end do

    if (num_lattices==2) then
       write(mylog,'("#                                                              #")')
       write(mylog,'("# Using common ratio of ",F12.6," for positive GP           #")')r_pos
       write(mylog,'("# Error in common ratio: ",E12.6,"                          #")')abs(r_new-r_pos)
       write(mylog,'("#                                                              #")')
    end if

    ! Compute r_new - common ratio of negative series
    r_neg = 1.1_dp
    k = 0
    do
       k = k + 1
       tmpsum = a_neg*(1.0_dp-r_neg**Ns)/(1.0_dp-r_neg)
       r_new = r_neg*(s_neg/tmpsum)**(1.0_dp/real(Ns,kind=dp))
       if (abs(r_new-r_neg)<=2.0_dp*epsilon(1.0_dp)) exit
       if (k>1000000) exit
       r_neg = r_new  ! converged
    end do

    if (num_lattices==2) then
       write(mylog,'("#                                                              #")')
       write(mylog,'("# Using common ratio of ",F12.6," for negative GP           #")')r_neg
       write(mylog,'("# Error in common ratio: ",E12.6,"                          #")')abs(r_new-r_neg)
       write(mylog,'("#                                                              #")')
    end if

    ! Compute mid-bin values and binwidths for negative series
    mu_u = -0.5_dp
    k = 0
    do ibin = nbins/2,1,-1
       mu_l = mu_u - a_neg*r_neg**k
       mu_bin(ibin)   = 0.5_dp*(mu_u+mu_l)
       binwidth(ibin) = mu_u-mu_l
       !write(*,'(2F15.6)')mu_l,mu_u
       mu_u = mu_l
       k = k + 1
    end do

    ! Middle bin is special case
    mu_bin(nbins/2+1)   = 0.0_dp
    binwidth(nbins/2+1) = 1.0_dp

    ! Min-bin values
    mu_l = 0.5_dp      ! upper limit of middle bin
    k = 0
    do ibin = nbins/2+2,nbins
       mu_u = mu_l + a_pos*r_pos**k
       mu_bin(ibin)   = 0.5_dp*(mu_u+mu_l)
       binwidth(ibin) = mu_u - mu_l 
       !write(*,'(2F15.6)')mu_l,mu_u
       mu_l = mu_u
       k = k + 1
    end do

    ! Compute average binwidth
    k = 0
    av_binwidth = 0.0_dp
    do ibin = 1,nbins
       av_binwidth = av_binwidth + binwidth(ibin)
    end do
    av_binwidth = av_binwidth/real(nbins,kind=dp)

    ! Parallel strategy
    select case (parallel_strategy)
       case('dd')

          ! Domain decomposition with overlap at domain boundaries
          bins_per_window = nbins/size 

          if (myrank==0) then

             my_start_bin = 1
             my_end_bin   = bins_per_window + window_overlap

             my_mu_min = mu_min
             my_mu_max = mu_min+sum(binwidth(1:my_end_bin))

          end if

          if (size > 1) then

             do irank = 1,size-2
                if (irank==myrank) then
                   
                   my_start_bin = irank * bins_per_window - window_overlap
                   my_end_bin   = (irank+1) * bins_per_window + window_overlap
                   
                   my_mu_min = mu_min+sum(binwidth(1:my_start_bin-1))
                   my_mu_max = mu_min+sum(binwidth(1:my_end_bin))
                
                end if
             end do
             
             if (myrank==size-1) then
                
                my_start_bin = myrank * bins_per_window - window_overlap
                my_end_bin = nbins             
                
                my_mu_min = mu_min+sum(binwidth(1:my_start_bin-1))
                my_mu_max = mu_max
                
             end if

          end if

          ! Identify which lattice should be active
          if (my_mu_max < 0.0_dp) ls = 1
          if (my_mu_min > 0.0_dp) ls = 2

          write(mylog,'("#                                                              #")')
          write(mylog,'("# This rank will use bin ",I6," to ",I6,"                      #")')my_start_bin,my_end_bin
          write(mylog,'("# Lower limit of mu: ",F12.6,"                              #")')my_mu_min
          write(mylog,'("# Upper limit of mu: ",F12.6,"                              #")')my_mu_max
          write(mylog,'("#                                                              #")')

       case('mw')
          
          ! Every MPI rank works on the same domain
          my_start_bin = 1
          my_end_bin = nbins

          my_mu_min = mu_min
          my_mu_max = mu_max

       case default
          stop 'Unknown parallel strategy'
       end select



    !--------------------------------------------------!
    ! Read any existing set of multicanonical weights. !
    ! These will be subsequently replaced if this is a !
    ! restart and we read checkpoint.dat               !
    !--------------------------------------------------!
    ! Store wl_factor from input file to compare with 
    ! that in existing eta_weights.dat. Use that in 
    ! eta_weights.dat if smaller.
    orig_wl_factor = wl_factor ! as read from input file

    if (num_lattices==2) then

       ! Rank 0 - reads weights and broadcasts them to all other ranks
       if (myrank==0) then
          
          inquire(file='eta_weights.dat',exist=lexist)
          if (lexist) then
             
             write(glog,'("#                                                              #")')
             write(glog,'("# Found existing eta_weights.dat in current directory.         #")')
             write(glog,'("#                                                              #")')
             
             open(unit=wgt,file='eta_weights.dat',status='old',iostat=ierr)
             if (ierr/=0) stop 'Error opening eta_weights.dat'
             
             read(wgt,'(A29,E20.12)')dumchar29,tmp_wl_factor
             if (ierr/=0) then
                backspace(wgt)
                read(wgt,'(A29,F20.12)',iostat=ierr)dumchar29,tmp_wl_factor
             end if
             
             if ( tmp_wl_factor > 1e-10 ) then
                wl_factor = min(wl_factor,tmp_wl_factor)
                if (samplerun) wl_factor = 0.0_dp
             end if
             k = 1
             do 
                read(wgt,*,end=10)dumreal,dumreal2
                weight(k) = dumreal2
                k = k + 1
             end do
             
10           close(wgt)
             
          end if
          
       end if

       ! MPI - make sure everyone has the same weights
       call comms_bcastreal(wl_factor,1)
       call comms_allreduce_eta(weight,nbins)

       ! Estimate the sum over unbiased counts in all bins in a fashion
       ! which is resistant to overflow. For use later when incrementing
       ! unbiased histogram.
       hits_per_bin = real(max_mc_cycles,kind=dp)-real(eq_mc_cycles,kind=dp)
       hits_per_bin = hits_per_bin*real(size*nwater,kind=dp)/real(nbins,kind=dp)

       ! First bin
       incr = hits_per_bin*av_binwidth

       ! log of expected first bin hits when 
       !                 = log(incr*exp(eta_weight(mu_bin(1))
       log_unbiased_norm = log(incr) + weight(1)


       do k = 2,nbins

          incr = hits_per_bin*av_binwidth

          if ( log_unbiased_norm > weight(k) + log(incr) ) then
             ! Use log(a) --> log(a+b) = log(a) + log(1+b/a)
             log_unbiased_norm = log_unbiased_norm + &
                  log(1.0_dp+incr*exp(weight(k)-log_unbiased_norm))             
          else
             ! Use log(a) --> log(a+b) = log(b) + log(1+a/b)
             log_unbiased_norm = log(incr) + weight(k) + &
                  log(1.0_dp+exp(log_unbiased_norm-weight(k))/incr)
          end if

       end do

       if (parallel_strategy=='dd') then

          ! Keep only the portion of the weights that relate to my window.
          weight(1:my_start_bin-1) = 0.0_dp
          weight(my_end_bin+1:nbins) = 0.0_dp

       end if
       
       if ( wl_factor < orig_wl_factor ) then
          write(mylog,'("#                                                              #")')
          write(mylog,'("# Using smaller weight increment read from eta_weights.dat     #")')
          write(mylog,'("#                                                              #")')
          firstcycle = .false.
       end if
       
    end if

    ! if this is a restart then call mc_checkpoint_load and update 
    ! all properties of the energy module
    write(rankstring,'(I3.3)')myrank
    inquire(file='checkpoint'//rankstring//'.dat.1',exist=lexist1)
    inquire(file='checkpoint'//rankstring//'.dat.2',exist=lexist2)

    if ( lexist1.or.lexist2 ) then

       write(mylog,'("#                                                              #")')
       write(mylog,'("# Found a checkpoint file on rank :",I6,"                      #")')myrank

       call mc_checkpoint_load()     

       outcycle = startcycle
       write(mylog,'("#  --restarting from cycle ",I10,"                          #")')startcycle
       write(mylog,'("#                                                              #")')

       do ils = 1,num_lattices
          volume(ils) = abs(util_determinant(hmatrix(:,:,ils)))
          call util_recipmatrix(hmatrix(:,:,ils),recip_matrix(:,:,ils))
          call compute_ivects(ils)
       end do

       if (num_lattices==2) call mc_check_chain_synchronisation()

       do ils = 1,num_lattices
          call compute_model_energy(ils)
       end do

    end if

    ! Compute the initial lattice-switching overlap parameter
    beta  = 1.0_dp/(kB*temperature)
    if (num_lattices==2) then
       ls_mu = model_energy(1) + pressure*volume(1) - model_energy(2) - pressure*volume(2) 
       if (leshift) ls_mu = ls_mu - ref_enthalpy(1) + ref_enthalpy(2)
       ls_mu = ls_mu*beta - real(Nwater,kind=dp)*log(volume(1)/volume(2))
    end if

    ! Create lookup array for 'other lattice'
    if (num_lattices==2) then
       partner_lattice = (/2,1/)
    else
       partner_lattice = (/1,1/)
    end if

    ! all walkers are in window if not domain decomposed
    if (parallel_strategy/='dd') walker_in_window=.true.


    return

  end subroutine mc_init

  subroutine mc_deinit()

    implicit none
    integer :: ierr

    deallocate(weight,histogram,binwidth,mu_bin,unbiased_hist,stat=ierr)
    deallocate(backup_local_real_energy)

    if (ierr/=0) stop 'Error releasing memory in mc module.'

    return

  end subroutine mc_deinit

  real(kind=dp) function eta_weight(mu_in)
    !------------------------------------------------------------------------------!
    ! Computes the current value of the weight function eta as a function of the   !
    ! overlap parameter mu.                                                        !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------!
    use userparams, only : nbins,eta_interp,eq_mc_cycles
    implicit none
   
    real(kind=dp),intent(in) :: mu_in
    integer             :: k

    real(kind=dp) :: gradient


    ! Never apply any weight during equilibration?
    !if (mc_cycle_num < eq_mc_cycles) return

    ! TODO - don't want to penalise walkers for not having reached their window yet
    if (.not.walker_in_window) return

    if (mu_in < my_mu_min ) then 
       eta_weight = huge(1.0_dp)
       return
    endif 
    if (mu_in > my_mu_max ) then
       eta_weight = huge(1.0_dp)
       return
    end if
    
    ! k is the bin number in which the current mu lies
    k = mu_to_bin(mu_in)

    if (eta_interp) then

       if (k==my_start_bin) then ! first bin special case

          gradient = 2.0_dp*(weight(k+1)-weight(k))/ (binwidth(k)+binwidth(k+1))
          eta_weight = weight(k) + (mu_in-mu_bin(k))*gradient
          
       elseif (k==my_end_bin) then ! last bin special case
       
          gradient = 2.0_dp*(weight(k)-weight(k-1))/(binwidth(k)+binwidth(k-1))
          eta_weight = weight(k) + (mu_in-mu_bin(k))*gradient

       else

          if (mu_in>mu_bin(k)) then

             gradient = 2.0_dp*(weight(k+1)-weight(k))/(binwidth(k)+binwidth(k+1))
             eta_weight =weight(k) + (mu_in-mu_bin(k))*gradient

          else

             gradient = 2.0_dp*(weight(k)-weight(k-1))/(binwidth(k)+binwidth(k-1))
             eta_weight =weight(k-1) + (mu_in-mu_bin(k-1))*gradient

          end if

       end if


    else

       eta_weight = weight(k)    

    end if
  
    return

  end function eta_weight

  subroutine mc_water_translation()
    !------------------------------------------------------------------------------!
    ! Computes a trial translation of a single molecule and computes the resulting !
    ! change in energy. The move is accepted or rejected appropriately.            !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------!
    use constants,  only : kb,invpi
    use random
    use userparams, only : nwater,mc_max_trans,temperature,num_lattices, &
         mc_ensemble,pressure,leshift
    use model,      only : ljspm,ljr,ls,volume,hmatrix,recip_matrix
    use energy,     only : compute_local_real_energy,model_energy,compute_model_energy
    implicit none

    ! local variables
    real(kind=dp) :: beta,x,y,z,norm,r,sx,sy,sz,diffkT
    real(kind=dp) :: zeta,eta_old,eta_new
    real(kind=dp),dimension(2)   :: old_energy,new_energy,deltaE
    real(kind=dp),dimension(2)   :: new_local_energy
    real(kind=dp),dimension(3,2) :: transvec

    ! loop counters/array indices
    integer :: imol,ilj,ils,lsn,idim

#ifdef DEBUG
    real(kind=dp),dimension(2) :: energy_dbg
#endif
    
    lsn = partner_lattice(ls)

    ! beta on this rank
    beta = 1.0_dp/(Kb*temperature)

    ! select a molecule at random
    x = random_uniform_random()
    imol = min(int(x*real(nwater,kind=dp)) + 1,nwater)

    mc_translations(imol) = mc_translations(imol) + 1

    !$omp parallel do default(shared) private(ils,ic)
    do ils = 1,num_lattices

       ! store old energy of this molecule
       old_energy(ils)    = compute_local_real_energy(imol,ils)
       
       !properties of energy module in case of rejection
       backup_model_energy(ils) = model_energy(ils)

       ! update model energy
       model_energy(ils) = model_energy(ils) - old_energy(ils)

    end do

    ! generate a random displacement
    x = random_uniform_random()
    y = random_uniform_random()
    z = random_uniform_random()

    x = 2.0_dp*x - 1.0_dp
    y = 2.0_dp*y - 1.0_dp
    z = 2.0_dp*z - 1.0_dp

    norm = 1.0_dp/sqrt(x*x+y*y+z*z)

    x = x*norm
    y = y*norm
    z = z*norm

    r = random_uniform_random()*2.0_dp - 1.0_dp

    x = x * mc_max_trans*r
    y = y * mc_max_trans*r
    z = z * mc_max_trans*r

    ! Size of move in scaled coordinates
    sx = recip_matrix(1,1,ls)*x + &
         recip_matrix(2,1,ls)*y + &
         recip_matrix(3,1,ls)*z
    sy = recip_matrix(1,2,ls)*x + &
         recip_matrix(2,2,ls)*y + &
         recip_matrix(3,2,ls)*z  
    sz = recip_matrix(1,3,ls)*x + &
         recip_matrix(2,3,ls)*y + &
         recip_matrix(3,3,ls)*z 

    sx = sx*0.5_dp*invPi 
    sy = sy*0.5_dp*invPi 
    sz = sz*0.5_dp*invPi 

    transvec(1,ls) = x
    transvec(2,ls) = y
    transvec(3,ls) = z

    ! move to make in the non-active lattice
    if (num_lattices==2) then
       do idim=1,3
          transvec(idim,lsn) = hmatrix(idim,1,lsn)*sx + &
                               hmatrix(idim,2,lsn)*sy + &
                               hmatrix(idim,3,lsn)*sz
       end do
    end if

    ! debug check for when cell vectors are the same
!!$    write(0,'("Move in lattice 1 :",3F15.6)')transvec(:,1)
!!$    write(0,'("Move in lattice 2 :",3F15.6)')transvec(:,2)


    ! translate molecule imol
    !$omp parallel do default(shared) private(ils,ilj,ic,ics)
    do ils = 1,num_lattices

       do ilj = 1,ljspm
          ljr(:,ilj,imol,ils) = ljr(:,ilj,imol,ils)  + transvec(:,ils)
       end do

       ! add in the real part
       new_local_energy(ils)  = compute_local_real_energy(imol,ils)
       new_energy(ils)        = new_local_energy(ils)

       ! update model energy
       model_energy(ils)  = model_energy(ils) + new_energy(ils)

       ! compute energy difference
       deltaE(ils)  = new_energy(ils) - old_energy(ils)

    end do

#ifdef DEBUG
    energy_dbg = model_energy
    call compute_model_energy(1)
    call compute_model_energy(2)
    energy_dbg = energy_dbg - model_energy
    if ( any(energy_dbg>1d-10) ) then
       write(0,'("DEBUG - translate 1, energy diff = ",2F15.6)')energy_dbg*hart_to_ev
    end if
#endif

    if (num_lattices == 1) then
       ! Single box branch
       diffkT = beta*deltaE(1)
    else

       ! Lattice switching branch

       ! update the order parameter
       eta_old = eta_weight(ls_mu)
       ls_mu   = ls_mu + (deltaE(1) - deltaE(2))*beta
       eta_new = eta_weight(ls_mu)
       
       diffkT = deltaE(ls)*beta + eta_new - eta_old


#ifdef MINU

       ! In which reference lattice does this move result in the lower energy?
       if (leshift) then
          lsn = minloc(model_energy+pressure*volume-ref_enthalpy,1)          
       else          
          lsn = minloc(model_energy+pressure*volume,1)          
       end if

       ! Is this is different to the current lattice then incorporate a switch into the move
       if (lsn/=ls) then
          if (mc_ensemble == 'npt' ) then
             diffkT = beta*model_energy(lsn) - beta*backup_model_energy(ls) + beta*pressure*(volume(lsn)-volume(ls)) &
                  -real(Nwater,kind=dp)*log(volume(lsn)/volume(ls)) + eta_new - eta_old      
             if (leshift) diffkT = diffkT - beta*ref_enthalpy(lsn) + beta*ref_enthalpy(ls)    
          else
             diffkT  = beta*model_energy(lsn) - beta*backup_model_energy(ls) + eta_new - eta_old
             if (leshift) diffkT = diffkT - beta*ref_enthalpy(lsn) + beta*ref_enthalpy(ls)
          end if
       end if
       
#endif
       
    end if ! num_lattices

    ! accept or reject 
    zeta = random_uniform_random()
    if (zeta<min(1.0_dp,exp(-diffkT))) then
              
       !===============!
       ! Move accepted !
       !===============!

#ifdef DEBUG
       !write(0,'("DEBUG - Accepted a translation with prob :",2F15.6)')min(1.0_dp,exp(-diffkT)),zeta
#endif

       mc_accepted_rsteps = mc_accepted_rsteps + 1
       dmu = abs(deltaE(1) - deltaE(2))*beta
       if ( dmu < min_dmu ) min_dmu = dmu
       if ( dmu > max_dmu ) max_dmu = dmu

            
!!$       if ( dmu < 0.0001 ) then
!!$          write(*,'("Warning in mc_water_translation : delta mu = ",F15.6)')dmu
!!$          write(*,'("DeltaE for this move was        : delta E  = ",F15.6)')deltaE(ls)
!!$       end if

       ! Update lattice
#ifdef MINU
       ls = lsn
#endif

    else

       !===============!
       ! Move rejected !
       !===============!

#ifdef DEBUG
       !write(0,'("DEBUG - Rejected a translation with prob :",2F15.6)')min(1.0_dp,exp(-diffkT)),zeta
#endif

       do ils = 1,num_lattices

          ! translate molecule imol back to where it came from
          do ilj = 1,ljspm
             ljr(:,ilj,imol,ils) = ljr(:,ilj,imol,ils) - transvec(:,ils)
          end do

          ! reset energy module
          model_energy(ils) = backup_model_energy(ils)

       end do

       ! update the order parameter
       if (num_lattices==2) ls_mu = ls_mu - (deltaE(1) - deltaE(2))*beta

#ifdef DEBUG
       energy_dbg = model_energy
       do ils = 1,num_lattices
          call compute_model_energy(ils)
       end do
       energy_dbg = energy_dbg - model_energy
       !if ( any(energy_dbg>1d-10) ) then
          write(0,'("DEBUG - translate 2, energy diff = ",2F15.6)')energy_dbg*hart_to_ev
       !end if
       !stop
#endif

    end if

    return

  end subroutine mc_water_translation


  subroutine mc_volume()
    !------------------------------------------------------------------------------!
    ! Computes a trial change in the cell vectors resulting from either an         !
    ! isotropic or anisotropic expansion/contraction. The fractional co-ordinates  !
    ! of the molecule centre-of-masses are maintained during the trial move.       !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------!
    use random
    use constants,  only : invpi,kb
    use userparams, only : nwater,mc_dv_max,pressure,nwater,num_lattices, &
                           temperature,leshift
    use util,       only : util_determinant,util_recipmatrix
    use energy,     only : model_energy,compute_model_energy,compute_ivects
    use model,      only : ljspm,ljr,volume,hmatrix,recip_matrix,ls,ref_ljr
    implicit none


    ! local variables
    real(kind=dp) :: beta,x,diffkT
    real(kind=dp) :: compare,old_eta=0.0_dp,new_eta=0.0_dp,old_mu
    real(kind=dp),dimension(2) :: old_energy,new_energy,deltaE,old_volume

    real(kind=dp),dimension(3,3,2) :: old_hmatrix,old_recip_matrix
    real(kind=dp),dimension(3,3)   :: delta_hmatrix

    real(kind=dp),dimension(3)   :: old_pos,new_pos,transvec

    integer :: lsn

    ! loop counters/array indices
    integer :: idim,jdim,ils,ilj,imol


    beta = 1.0_dp/(kB*temperature)

    lsn = partner_lattice(ls)

    ! store properties of energy module in case of rejection
    backup_model_energy(1:num_lattices)  = model_energy(1:num_lattices)

    ! store old model energy
    old_energy(1:num_lattices) = model_energy(1:num_lattices)

    do ils = 1,num_lattices
       call util_recipmatrix(hmatrix(:,:,ils),recip_matrix(:,:,ils))
    end do

    old_hmatrix(:,:,1:num_lattices)      = hmatrix(:,:,1:num_lattices)
    old_recip_matrix(:,:,1:num_lattices) = recip_matrix(:,:,1:num_lattices)
    old_volume(1:num_lattices)           = volume(1:num_lattices)

    ! now perform a volume step
    x = random_uniform_random()
    idim = int(x*3.0_dp)+1
    x = random_uniform_random()
    jdim = int(x*3.0_dp)+1

    x = random_uniform_random()

    delta_hmatrix = 0.0_dp

    delta_hmatrix(idim,jdim) =  (2.0_dp*x-1.0_dp)*mc_dv_max
    delta_hmatrix(jdim,idim) =  delta_hmatrix(idim,jdim)

    hmatrix(:,:,1) = hmatrix(:,:,1) + delta_hmatrix
    if (num_lattices==2) hmatrix(:,:,2) = hmatrix(:,:,2) + delta_hmatrix

    !$omp parallel do private(ils,old_pos,new_pos,imol,ilj,ics,idim,transvec)
    do ils = 1,num_lattices

       ! update positions of molecules centers due to box scaling
       do imol = 1,nwater

          old_pos(:) = ljr(:,1,imol,ils)

          ! compute fractional co-ordinates
          new_pos(1) = recip_matrix(1,1,ils)*old_pos(1) + &
                       recip_matrix(2,1,ils)*old_pos(2) + &
                       recip_matrix(3,1,ils)*old_pos(3)
          new_pos(2) = recip_matrix(1,2,ils)*old_pos(1) + &
                       recip_matrix(2,2,ils)*old_pos(2) + &
                       recip_matrix(3,2,ils)*old_pos(3)  
          new_pos(3) = recip_matrix(1,3,ils)*old_pos(1) + &
                       recip_matrix(2,3,ils)*old_pos(2) + &
                       recip_matrix(3,3,ils)*old_pos(3)        

          new_pos = new_pos*0.5_dp*invPi 

          do idim=1,3
             transvec(idim) = hmatrix(idim,1,ils)*new_pos(1) + &
                              hmatrix(idim,2,ils)*new_pos(2) + &
                              hmatrix(idim,3,ils)*new_pos(3)
          end do

          transvec   = transvec - old_pos

          do ilj = 1,ljspm
             ljr(:,ilj,imol,ils) = ljr(:,ilj,imol,ils) + transvec(:)
          end do

       end do

       ! update positions of reference centers due to box scaling
       do imol = 1,nwater

          old_pos(:) = ref_ljr(:,1,imol,ils)

          ! compute fractional co-ordinates
          new_pos(1) = recip_matrix(1,1,ils)*old_pos(1) + &
                       recip_matrix(2,1,ils)*old_pos(2) + &
                       recip_matrix(3,1,ils)*old_pos(3)
          new_pos(2) = recip_matrix(1,2,ils)*old_pos(1) + &
                       recip_matrix(2,2,ils)*old_pos(2) + &
                       recip_matrix(3,2,ils)*old_pos(3)  
          new_pos(3) = recip_matrix(1,3,ils)*old_pos(1) + &
                       recip_matrix(2,3,ils)*old_pos(2) + &
                       recip_matrix(3,3,ils)*old_pos(3)        

          new_pos = new_pos*0.5_dp*invPi 

          do idim=1,3
             transvec(idim) = hmatrix(idim,1,ils)*new_pos(1) + &
                              hmatrix(idim,2,ils)*new_pos(2) + &
                              hmatrix(idim,3,ils)*new_pos(3)
          end do

          transvec   = transvec - old_pos

          do ilj = 1,ljspm
             ref_ljr(:,ilj,imol,ils) = ref_ljr(:,ilj,imol,ils) + transvec(:)
          end do

       end do


       ! update volume and anything dependent on
       ! reciprocal lattice.
       volume(ils) = abs(util_determinant(hmatrix(:,:,ils)))
       call util_recipmatrix(hmatrix(:,:,ils),recip_matrix(:,:,ils))
       call compute_ivects(ils) 
       call compute_model_energy(ils)

       new_energy(ils) = model_energy(ils)

    end do

    ! change in energy
    deltaE(1:num_lattices) = new_energy(1:num_lattices) - old_energy(1:num_lattices)

    if (num_lattices==2) then

       old_eta = eta_weight(ls_mu)
       old_mu  = ls_mu
       ls_mu   = ( model_energy(1)+pressure*volume(1) ) - ( model_energy(2)+pressure*volume(2) )  
       if (leshift) ls_mu = ls_mu - ref_enthalpy(1) + ref_enthalpy(2)
       ls_mu   = ls_mu*beta - real(Nwater,kind=dp)*log(volume(1)/volume(2))
       new_eta = eta_weight(ls_mu)

    end if

    ! spin new random number
    x = random_uniform_random()

    diffkT = beta*deltaE(ls) + new_eta - old_eta + beta*Pressure*(volume(ls)-old_volume(ls)) &
             - real(Nwater,kind=dp)*log(volume(ls)/old_volume(ls))

    if (num_lattices==2) then
           
#ifdef MINU
    
       ! In which reference lattice does this move result in the lower energy?
       if (leshift) then          
          lsn     = minloc(model_energy+pressure*volume-ref_enthalpy,1)          
       else          
          lsn     = minloc(model_energy+pressure*volume,1)          
       end if

       ! Is this is different to the current lattice then incorporate a switch into the move
       if (lsn/=ls) then
          diffkT = beta*model_energy(lsn) - beta*backup_model_energy(ls) + beta*pressure*(volume(lsn)-old_volume(ls)) &
               -real(Nwater,kind=dp)*log(volume(lsn)/old_volume(ls)) + new_eta - old_eta     
          if (leshift) diffkT = diffkT - beta*ref_enthalpy(lsn) + beta*ref_enthalpy(ls)  
       end if
       
#endif


    end if ! lattices
              
    ! calculate acceptance probability based on the current lattice
    compare = exp(-diffkT)
    compare = min(1.0_dp,compare)

    if (x < compare) then
       mc_accepted_vsteps = mc_accepted_vsteps + 1

       if (num_lattices == 2) then
          dmu = abs(old_mu - ls_mu)       
          if ( dmu < min_dmu ) min_dmu = dmu
          if ( dmu > max_dmu ) max_dmu = dmu

!!$       if ( dmu < 0.0001 ) then
!!$          write(*,'("Warning in mc_volume            : delta mu = ",F15.6)')dmu
!!$          write(*,'("DeltaE for this move was        : delta E  = ",F15.6)')deltaE(ls)
!!$       end if

       end if


#ifdef MINU
       ! Update active lattice
       ls = lsn
#endif

    else

       ! reject
       volume(1:num_lattices)      = old_volume(1:num_lattices)
       hmatrix(:,:,1:num_lattices) = old_hmatrix(:,:,1:num_lattices)

       !$omp parallel do private(ils,old_pos,new_pos,imol,ilj,ics,idim,transvec)
       do ils = 1,num_lattices

          ! update positions of molecules centers due to box scaling
          do imol = 1,nwater

             old_pos(:) = ljr(:,1,imol,ils)

             ! compute fractional co-ordinates
             new_pos(1) = recip_matrix(1,1,ils)*old_pos(1) + &
                          recip_matrix(2,1,ils)*old_pos(2) + &
                          recip_matrix(3,1,ils)*old_pos(3)
             new_pos(2) = recip_matrix(1,2,ils)*old_pos(1) + &
                          recip_matrix(2,2,ils)*old_pos(2) + &
                          recip_matrix(3,2,ils)*old_pos(3)  
             new_pos(3) = recip_matrix(1,3,ils)*old_pos(1) + &
                          recip_matrix(2,3,ils)*old_pos(2) + &
                          recip_matrix(3,3,ils)*old_pos(3)        

             new_pos = new_pos*0.5_dp*invPi 

             do idim=1,3
                transvec(idim) = hmatrix(idim,1,ils)*new_pos(1) + &
                                 hmatrix(idim,2,ils)*new_pos(2) + &
                                 hmatrix(idim,3,ils)*new_pos(3)
             end do

             transvec   = transvec - old_pos

             do ilj = 1,ljspm
                ljr(:,ilj,imol,ils) = ljr(:,ilj,imol,ils) + transvec(:)
             end do

          end do

          ! update positions of reference centers due to box scaling
          do imol = 1,nwater

             old_pos(:) = ref_ljr(:,1,imol,ils)

             ! compute fractional co-ordinates
             new_pos(1) = recip_matrix(1,1,ils)*old_pos(1) + &
                          recip_matrix(2,1,ils)*old_pos(2) + &
                          recip_matrix(3,1,ils)*old_pos(3)
             new_pos(2) = recip_matrix(1,2,ils)*old_pos(1) + &
                          recip_matrix(2,2,ils)*old_pos(2) + &
                          recip_matrix(3,2,ils)*old_pos(3)  
             new_pos(3) = recip_matrix(1,3,ils)*old_pos(1) + &
                          recip_matrix(2,3,ils)*old_pos(2) + &
                          recip_matrix(3,3,ils)*old_pos(3)        

             new_pos = new_pos*0.5_dp*invPi 

             do idim=1,3
                transvec(idim) = hmatrix(idim,1,ils)*new_pos(1) + &
                                 hmatrix(idim,2,ils)*new_pos(2) + &
                                 hmatrix(idim,3,ils)*new_pos(3)
             end do

             transvec   = transvec - old_pos

             do ilj = 1,ljspm
                ref_ljr(:,ilj,imol,ils) = ref_ljr(:,ilj,imol,ils) + transvec(:)
             end do

          end do

       end do


       recip_matrix(:,:,1:num_lattices) = old_recip_matrix(:,:,1:num_lattices)

!!$       kx            = backup_kx
!!$       ky            = backup_ky
!!$       kz            = backup_kz

       do ils=1,num_lattices
          call compute_ivects(ils)
       end do

!!$       rho_k         = backup_rho_k
!!$       expor         = backup_expor
!!$       recip_energy  = backup_recip_energy
       model_energy(1:num_lattices)  = backup_model_energy(1:num_lattices)

!!$       local_real_energy = backup_local_real_energy

       if (num_lattices==2) then
          ls_mu  = ( model_energy(1) + pressure*volume(1) ) - ( model_energy(2) + pressure*volume(2) )  
          if (leshift) ls_mu = ls_mu - ref_enthalpy(1) + ref_enthalpy(2)
          ls_mu  = ls_mu*beta - real(Nwater,kind=dp)*log(volume(1)/volume(2))
       end if

    end if

    return

  end subroutine mc_volume

  subroutine mc_lattice_switch()

    use constants, only  : kB
    use model,  only     : ls,volume
    use energy, only     : model_energy
    use random, only     : random_uniform_random
    use userparams, only : pressure,nwater,mc_ensemble,temperature,leshift,num_lattices

    implicit none
    real(kind=dp) :: beta,compare,x,old_eta,new_eta,diffkT
    integer :: lsn

    if (num_lattices==1) then
       stop 'Called mc_lattice switch with only 1 lattice!'
    end if
    
    !return
    beta  = 1.0_dp/(kB*temperature)

    lsn = partner_lattice(ls)


    old_eta = eta_weight(ls_mu)
    new_eta = eta_weight(ls_mu)


    if (mc_ensemble == 'npt' ) then

       diffkT = beta*model_energy(lsn) - beta*model_energy(ls) + beta*pressure*(volume(lsn)-volume(ls)) &
            -real(Nwater,kind=dp)*log(volume(lsn)/volume(ls)) + new_eta - old_eta

       if (leshift) diffkT = diffkT - beta*ref_enthalpy(lsn) + beta*ref_enthalpy(ls)

       compare = min(1.0_dp,exp(-diffkT))
    else
       diffkT  = beta*model_energy(lsn) - beta*model_energy(ls) + new_eta - old_eta
       if (leshift) diffkT = diffkT - beta*ref_enthalpy(lsn) + beta*ref_enthalpy(ls)
       compare = min(1.0_dp,exp(-diffkT)) 
    end if

    x = random_uniform_random()

    if (x < compare) then

      
       mc_accepted_swtch = mc_accepted_swtch + 1
    
       ls_mu  = ( model_energy(1) + pressure*volume(1) ) - ( model_energy(2) + pressure*volume(2) )  
       if (leshift) ls_mu = ls_mu - ref_enthalpy(1) + ref_enthalpy(2)
       ls_mu  = ls_mu*beta - real(Nwater,kind=dp)*log(volume(1)/volume(2))

       ls = lsn

    end if
    
    return


  end subroutine mc_lattice_switch


  subroutine mc_update_wl_bins()
    !--------------------------------------------------------------------!
    ! Updated the Wang-Landau histogram and the current estimate of the  !
    ! multicanonical weight function. Called after every MC move.        !
    !--------------------------------------------------------------------!
    ! D.Quigley July 2010                                                !
    !--------------------------------------------------------------------!
    use userparams,   only : nbins,mu_min,mu_max,wl_swetnam,wl_factor, &
                             wl_alpha,nwater,samplerun,eq_mc_cycles
    use random,       only : random_uniform_random
    use util,         only : util_determinant
    use model,        only : ls
    implicit none

    real(kind=dp) :: incr,minbin,binfrac
    integer       :: k,i,lsn

    ! Give the system a chance to equilibrate
    if ( mc_cycle_num < eq_mc_cycles ) return

    ! Update the histogram
    k = mu_to_bin(ls_mu)
    if ( (k<1).or.(k>nbins) ) return

    histogram(k)  = histogram(k)  + av_binwidth/binwidth(k)

    ! If this is a sample run then also accumulate an unbiased histogram
    ! but then stop and do nothing to the weights
    if (samplerun) then

       incr = av_binwidth/binwidth(k)
       unbiased_hist(k) = unbiased_hist(k) + incr*exp(eta_weight(ls_mu)-log_unbiased_norm)
       !write(0,*)incr*exp(eta_weight(ls_mu)),incr*exp(eta_weight(ls_mu)-log_unbiased_norm)
       return

    end if

    ! Use Adam Swetnam's method for computing the optimal wl_factor
    ! based on the current state of the histogram.
    if (wl_swetnam) then

       sumhist = sumhist + 1.0_dp ! total probability 

       ! Compute the optimal increment at the current flatness
       wl_factor = 0.0_dp
       do i = 1,nbins
          binfrac   = binwidth(i)/(mu_max-mu_min-1.0_dp)
          wl_factor = wl_factor + (histogram(i)*binwidth(i)/sumhist - binfrac)**2
       end do
       wl_factor = sqrt(wl_factor/real(nbins,kind=dp)) 

       ! wl_factor now contains the R.M.S. deviation from a
       ! flat probability distribution, defining a suitable
       ! energy increment
       wl_factor = log(wl_factor)
       wl_factor = wl_factor*wl_alpha*real(nbins,kind=dp)
       wl_factor = min(wl_factor,orig_wl_factor)

    else if (wl_invt_active) then

       wl_factor = min(wl_factor,real(nbins,kind=dp)/real((mc_cycle_num)*nwater,kind=dp))

    end if

    ! Define 'other' phase
    lsn = partner_lattice(ls)

    ! Set increment
    incr = wl_factor ! default

!!$ More trouble than its worth.
!!$    if ( weight(k)<tiny(1.0_dp) ) then
!!$       ! If this bin hasn't been visited use smallest non-zero weight
!!$       incr = maxval(weight(:),1)
!!$       do i = 1,nbins
!!$          if ( (weight(i)>tiny(1.0_dp)).and.(weight(i)<incr) ) then
!!$             incr = weight(i)
!!$          end if
!!$       end do
!!$       incr = max(incr,wl_factor)
!!$    end if

    ! Increment the current bin only
    weight(k)  = weight(k)  + av_binwidth*incr/binwidth(k)

    minbin = minval(weight(my_start_bin:my_end_bin),1)
    do i = my_start_bin,my_end_bin
       weight(i) = weight(i) - minbin
    end do

    return

  end subroutine mc_update_wl_bins

  subroutine mc_monitor_stats
    !------------------------------------------------------------------------------------!
    ! Adjusts (if appropriate) the translation step size e.t.c. toward a target          !
    ! acceptance ratio of 50%, and reports to mylog various statistics about the various !
    ! trial move types.                                                                  !
    !------------------------------------------------------------------------------------!
    ! D.Quigley July 2010                                                                !
    !------------------------------------------------------------------------------------!
    use comms,      only : myrank,comms_allreduce_hist,comms_allreduce_eta, &
                           comms_join_eta, comms_get_max,comms_allreduce_uhist, &
                           comms_join_uhist
    use io,         only : mylog
    use constants,  only : bohr_to_ang,kb,hart_to_ev
    use userparams, only : eq_adjust_mc,mc_max_trans,mc_dv_max, &
                           mc_target_ratio,mc_ensemble,nwater,temperature,nbins, &
                           mu_min,mu_max,wl_factor,allow_trans, parallel_strategy, &
                           allow_switch,eq_mc_cycles,monitor_int,num_lattices, &
                           window_overlap,samplerun
    use energy,     only : model_energy,compute_model_energy
    implicit none

    real(kind=dp)              :: atr,avr,alr,max_wl_factor
    real(kind=dp),dimension(2) :: orig_energy
    integer                    :: ierr,k,ils

    real(kind=dp),allocatable,dimension(:) :: joined

    character(19) :: wgtstring
    character(17) :: hisstring
    character(26) :: ubhstring

    ! Compute acceptance ratios for each move type
    atr = real(mc_accepted_rsteps,kind=dp)/real(mc_attempted_rsteps,kind=dp)
    avr = real(mc_accepted_vsteps,kind=dp)/real(mc_attempted_vsteps,kind=dp)
    alr = real(mc_accepted_swtch,kind=dp)/real(mc_attempted_swtch,kind=dp)

    ! If adjusting toward a target acceptance ratio then do so.
    ! (Shouldn't be doing this for production runs)
    if (eq_adjust_mc.and.(mc_cycle_num<eq_mc_cycles)) then
       mc_max_trans = max(mc_max_trans * atr/mc_target_ratio,0.1_dp)
       mc_dv_max = max(mc_dv_max * avr/mc_target_ratio,0.0001_dp)
    end if

    write(mylog,'("#                                                              #")')   
    write(mylog,'("# Monte-Carlo statistics at cycle ",I10,"                   #")')mc_cycle_num
    write(mylog,'("#--------------------------------------------------------------#")')
    write(mylog,'("#                                                              #")')   
    write(mylog,'("#                                                              #")')   
    if (allow_trans) then
       write(mylog,'("# Accepted ",F8.2," percent of ",I8," translation moves      #")')atr*100.0_dp,mc_attempted_rsteps
       write(mylog,'("#                                                              #")')  
       write(mylog,'("# Translations attempted per molecule :  average = ",I10,"  #")')int(real(sum(mc_translations),kind=dp) & 
            /real(nwater,kind=dp))
       write(mylog,'("#                                            min = ",I10,"  #")')minval(mc_translations,1)
       write(mylog,'("#                                            max = ",I10,"  #")')maxval(mc_translations,1)
       write(mylog,'("#                                                              #")')   
    end if
    if (mc_ensemble=='npt' ) then
       write(mylog,'("# Accepted ",F10.2," percent of ",I8," box moves            #")')avr*100.0_dp,mc_attempted_vsteps
    end if
    if (allow_switch) then
       write(mylog,'("# Accepted ",F10.2," percent of ",I8," lattice switches     #")')alr*100.0_dp,mc_attempted_swtch
    end if
    write(mylog,'("#                                                              #")')  

    if (eq_adjust_mc.and.(mc_cycle_num<eq_mc_cycles)) then
       write(mylog,'("# Adjusting move parameters toward a target acceptance ratio   #")')
       write(mylog,'("#                                                              #")')  
       write(mylog,'("# Maximum molecule translation :", F15.6," Ang            #")')mc_max_trans*bohr_to_ang
       write(mylog,'("# Maximum cell displacement    :", F15.6," Ang            #")')mc_dv_max*bohr_to_ang
    else
       write(mylog,'("# Adjustment of move parameters is disabled                    #")')
    end if
       write(mylog,'("#                                                              #")')  

    write(mylog,'("# Average excitation energy:                                   #")')
    write(mylog,'("# --------------------------                                   #")')
       write(mylog,'("#                                                              #")')  
    do ils = 1,num_lattices
       write(mylog,'("# lattice ",I6," = ",F15.6," kT per D.O.F.               #")') &
           ils, (average_energy(ils)/real(monitor_int,kind=dp)-ref_enthalpy(ils))/(Kb*temperature*3*nwater)
    end do

    write(mylog,'("#                                                              #")')  

    write(mylog,'("# Overlap bins spanned per MC move : min = ",F15.6,"     #")')min_dmu*real(nbins,kind=dp)/(mu_max - mu_min)
    write(mylog,'("# Overlap bins spanned per MC move : max = ",F15.6,"     #")')max_dmu*real(nbins,kind=dp)/(mu_max - mu_min)
    write(mylog,'("#                                                              #")')  


    write(mylog,'("# Checking accumulated energies                                #")')
    write(mylog,'("# -----------------------------                                #")')
    write(mylog,'("#  latt   Stored (eV)        Computed (eV)       drift (eV)    #")')
    write(mylog,'("#                                                              #")')  

    orig_energy(1:num_lattices) = model_energy(1:num_lattices)*hart_to_ev

    do ils = 1,num_lattices
       call compute_model_energy(ils)
       write(mylog,'("#  ",I1,3x,F12.6,9x,F12.6,6x,F12.6,4x," #")')ils,orig_energy(ils),model_energy(ils)*hart_to_ev, &
            orig_energy(ils)-model_energy(ils)*hart_to_ev
    end do

    write(mylog,'("#                                                              #")')  


    ! Reset accumulators
    mc_accepted_rsteps = 0
    mc_accepted_vsteps = 0
    mc_accepted_swtch  = 0

    mc_attempted_rsteps = 0
    mc_attempted_vsteps = 0
    mc_attempted_swtch  = 0

    mc_translations = 0
    average_energy = 0.0_dp

    max_dmu = 0.0_dp
    min_dmu = huge(1.0_dp)


    if (num_lattices == 2 ) then

       if ( parallel_strategy == 'mw' ) then

          ! MPI - Synchronise histogram and weights here
          call comms_allreduce_hist(histogram,nbins)
          call comms_allreduce_eta(weight,nbins)

          if (samplerun) call comms_allreduce_uhist(unbiased_hist,nbins)

          if (myrank==0) then

             ! Update multicanonical weight and histogram files
             if (.not.samplerun) then
                open(unit=wgt,file='eta_weights.dat',status='replace',iostat=ierr)
                if (ierr/=0) stop 'Error opening eta_weights.dat for output'
                write(wgt,'("#Current energy increment = ",E20.12)')wl_factor
             end if

             open(unit=his,file='histogram.dat',status='replace',iostat=ierr)
             write(his,'("#Current energy increment = ",E20.12)')wl_factor

             open(unit=ubh,file='unbiased_histogram.dat',status='replace',iostat=ierr)
             write(his,'("#Current energy increment = ",E20.12)')wl_factor

             
             do k=1,nbins
                if(.not.samplerun) write(wgt,*)mu_bin(k),weight(k)
                write(his,*)mu_bin(k),histogram(k)
                write(ubh,*)mu_bin(k),unbiased_hist(k)
             end do
             
             if (.not.samplerun) close(wgt)
             close(his)
             close(ubh)
             
          end if
          
       else if (parallel_strategy == 'dd') then
          
          ! Write local weights to seperate files
          write(wgtstring,'("eta_weights_",I3.3,".dat")')myrank
          open(unit=wgt,file=wgtstring,status='replace',iostat=ierr)
          if (ierr/=0) stop 'Error opening weights file for output'
          write(wgt,'("#Current energy increment = ",E20.12)')wl_factor          

          write(hisstring,'("histogram_",I3.3,".dat")')myrank
          open(unit=his,file=hisstring,status='replace',iostat=ierr)
          if (ierr/=0) stop 'Error opening histogram file for output'
          write(his,'("#Current energy increment = ",E20.12)')wl_factor

          write(ubhstring,'("unbiased_histogram_",I3.3,".dat")')myrank
          open(unit=ubh,file=ubhstring,status='replace',iostat=ierr)
          if (ierr/=0) stop 'Error opening unbiased histogram file for output'
          write(ubh,'("#Current energy increment = ",E20.12)')wl_factor


          do k=my_start_bin, my_end_bin
             write(wgt,*)mu_bin(k),weight(k)
             write(his,*)mu_bin(k),histogram(k)
             write(ubh,*)mu_bin(k),unbiased_hist(k)
          end do          

          close(wgt)
          close(his)
          close(ubh)

          ! Fit all the windows together
          allocate(joined(1:nbins),stat=ierr)
          if (ierr/=0) stop 'Error in mc_monitor_stats - could not allocate joined weights'
          call comms_join_eta(weight,nbins,window_overlap,joined)

          ! Should write the maximum wl_factor over all ranks here?
          call comms_get_max(wl_factor,max_wl_factor)

          if (myrank==0) then
             
             open(unit=wgt,file='eta_weights.dat',status='replace',iostat=ierr)
             if (ierr/=0) stop 'Error opening eta_weights.dat for output'
             write(wgt,'("#Current energy increment = ",E20.12)')max_wl_factor

             do k=1,nbins
                write(wgt,*)mu_bin(k),joined(k)
             end do
             
             close(wgt)

          end if

          ! Also join the unbiased histogram
          if (samplerun) then

             call comms_join_uhist(unbiased_hist,nbins,window_overlap,joined)
             
             if (myrank==0) then

                open(unit=ubh,file='unbiased_histogram.dat',status='replace',iostat=ierr)
                if (ierr/=0) stop 'Error opening unbiased_histogram.dat for output'
                write(ubh,'("#Current energy increment = ",E20.12)')0.0
                
                do k=1,nbins
                   write(ubh,*)mu_bin(k),joined(k)
                end do
                
                close(ubh)

             end if

          end if

          deallocate(joined,stat=ierr)
          if (ierr/=0) stop 'Error deallocating joined weights'

       else
          ! Uknown parallelisation
          stop 'Unknown parallel_strategy in mc_monitor_stats'
       end if

    end if


  end subroutine mc_monitor_stats

  subroutine mc_check_flatness()
    !--------------------------------------------------------!
    ! Check for a flat histogram and update the Wang-Landau  !
    ! modification factor accordingly                        !
    !--------------------------------------------------------!
    ! D.Quigley July 2010                                    !
    !--------------------------------------------------------!
    use comms,      only : myrank,comms_allreduce_hist, &
                           comms_set_histogram,comms_bcastlog
    use io,         only : glog,mylog
    use userparams, only : wl_schedule,nbins,wl_flattol, &
                           wl_minhist,wl_factor,wl_swetnam, &
                           nwater,wl_useinvt,samplerun,invt_dump_int, &
                           parallel_strategy
    implicit none

    character(20) :: wlstring
    real(kind=dp) :: av,wl_invt

    integer :: mini,count,k,ierr
    logical :: flat,lexist

    logical, save :: histogram_reset = .false.

    ! Never check for flatness if this is a sample run, i.e.
    ! a run with fixed weights and no histogram resets.
    if (samplerun.or.(sum(histogram)<tiny(1.0_dp))) return

    ! MPI - synchronise histogram
    if (parallel_strategy=='mw') then ! multiple walkers
       call comms_allreduce_hist(histogram,nbins)
    end if

    ! If the Wang-Landau f has never been modified
    ! reset the histogram once we've visited every
    ! bin at wl_minhist times. 
    mini = nint(minval(histogram))
    if ( (firstcycle).and.(.not.histogram_reset ) )then
       if (mini>wl_minhist) then
          histogram_reset = .true.
          histogram = 0.0_dp
          call comms_set_histogram(histogram,nbins)
          return
       end if
    end if
   
    ! Compute mean histogram value
    av = 0.0_dp
    count = 0
    do k = my_start_bin,my_end_bin 
       av = av + real(histogram(k),kind=dp)
       count = count + 1
    end do
    av = av / real(count,kind=dp)

    if (parallel_strategy=='mw') then

       if (myrank==0) then
          write(glog,'("# Checking flatness of histogram at cycle ",I10,"           #")')mc_cycle_num
          write(glog,'("# --------------------------------------------------           #")')   
          write(glog,'("#                                                              #")')   
          write(glog,'("# Most  populated histogram bin = ",F10.4," % of mean         #")')100.0_dp*maxval(histogram,1)/av
          write(glog,'("# Least populated histogram bin = ",F10.4," % of mean         #")')100.0_dp*minval(histogram,1)/av
          write(glog,'("#                                                              #")')   
       end if

    elseif (parallel_strategy=='dd') then

          write(mylog,'("# Checking flatness of histogram at cycle ",I10,"           #")')mc_cycle_num
          write(mylog,'("# --------------------------------------------------           #")')   
          write(mylog,'("#                                                              #")')   
          write(mylog,'("# Most  populated histogram bin = ",F10.4," % of mean         #")')100.0_dp*maxval(histogram,1)/av
          write(mylog,'("# Least populated histogram bin = ",F10.4," % of mean         #")')100.0_dp*minval(histogram,1)/av
          write(mylog,'("#                                                              #")')   


    end if



    ! Using standard Wang-Landau with f only updated
    ! every time the histogram is flat.
    if (.not.(wl_invt_active.or.wl_swetnam) ) then

       !==========================================!
       ! Check for flatness of current histogram  !
       !==========================================!

       if ( wl_schedule == 0 ) then
          
          ! Flatness criterion is that all histogram bins
          ! must be within 100*wl_flattol % of the mean.
          flat = .true.
          do k = my_start_bin,my_end_bin 
             if (abs(real(histogram(k),kind=dp) - av)/av > wl_flattol ) flat=.false.
          end do

       elseif (wl_schedule == 1) then
          
          ! Flatness criterion is that all histogram bins
          ! must have been visited wl_minhist times.
          mini = nint(minval(histogram(my_start_bin:my_end_bin)))
          flat = .true.
          if ( mini<wl_minhist ) flat=.false.
          
       elseif ( wl_schedule == 2 ) then
          
          ! Flatness criterion is that all histogram bins
          ! must be above 100*wl_flattol % of the mean.
          flat = .true.
          do k = my_start_bin,my_end_bin 
             if ( real(histogram(k),kind=dp) < (1.0_dp-wl_flattol)*av ) flat=.false.
          end do

       else
          stop 'Error - unknown wl_schedule value'
       end if
       
       ! Make sure everyone agrees that the histogram is flat
       if (parallel_strategy=='mw') call comms_bcastlog(flat,1)


       if (flat) then
          
          if (parallel_strategy=='mw') then

             ! Shift everything down such that minimum is at zero. 
             av = weight(nbins/2+1)
             do k = 1,nbins
                weight(k) = weight(k) - av
             end do
          
             if (myrank==0) then
             
                inquire(file='wlf.dat',exist=lexist)
                if (.not.lexist.or.firstcycle) then
                   ! Open a new file 
                   open(unit=wlf,file='wlf.dat',status='replace',iostat=ierr)
                else
                   ! Append to existing file
                   open(unit=wlf,file='wlf.dat',status='old',position='append',iostat=ierr)
                end if
                if (ierr/=0) stop 'Error opening wlf.dat'
                
                write(wlf,'(I10,E20.12)')mc_cycle_num,wl_factor
                write(wlf,'(I10,E20.12)')mc_cycle_num,0.5_dp*wl_factor
                close(wlf)
                
                ! Write histogram and weights tagged with the current f
                write(wlstring,'(F20.12)')wl_factor
                
                open(unit=wgt,file="eta_weights.dat_"//adjustl(wlstring),status='replace',iostat=ierr)
                if (ierr/=0) stop 'Error opening eta_weights.dat for output'
                write(wgt,'("#Current energy increment = ",E20.12)')wl_factor
                
                open(unit=his,file="histogram.dat_"//adjustl(wlstring),status='replace',iostat=ierr)
                if (ierr/=0) stop 'Error opening histogram.dat for output'
                write(his,'("#Current energy increment = ",E20.12)')wl_factor
                
                do k=1,nbins
                   write(wgt,*)mu_bin(k),weight(k)
                   write(his,*)mu_bin(k),histogram(k)
                end do
                close(wgt)
                close(his)
                
             end if ! end if myrank==0
             
             ! Reset histogram and update wl_factor
             histogram = 0.0_dp
             call comms_set_histogram(histogram,nbins)
             wl_factor = wl_factor*0.5_dp
             if (myrank==0) write(glog,'("#                                                              #")') 
             if (myrank==0) write(glog,'("# Flatness criterion satisfied - updating wl_factor            #")')
             if (myrank==0) write(glog,'("#                                                              #")') 
             
             firstcycle = .false.
 
          elseif(parallel_strategy=='dd') then

             ! Write current weights and histogram for the present rank only

             ! Reset histogram and update wl_factor
             histogram = 0.0_dp
             wl_factor = wl_factor*0.5_dp
             write(mylog,'("#                                                              #")') 
             write(mylog,'("# Flatness criterion satisfied - updating wl_factor            #")')
             write(mylog,'("#                                                              #")') 

             firstcycle = .false.
            
          else
             stop 'Unknown parallel strategy in mc_check_flatness'             
          end if


          
       end if ! end if flat
       
       ! switch to 1/t if appropriate
       wl_invt = real(nbins,kind=dp)/real(mc_cycle_num*nwater,kind=dp)
       if ( (wl_factor < wl_invt).and.(wl_factor>tiny(1.0_dp)) ) then
          if ( wl_useinvt ) then
             wl_invt_active = .true.
             if (myrank==0) write(glog,'("# Switching to 1/t method                     #")')
             wl_factor = wl_invt
          end if
       end if
         
    else
       
       ! If computing f at every step, having either switched to /1t
       ! or if using Adam Swetnam's formula.
       if ( (myrank==0).and.(mod(mc_cycle_num,invt_dump_int)==0) ) then

             ! Write the current wl_factor to wlf.dat
             inquire(file='wlf.dat',exist=lexist)
             if (lexist) then
                open(unit=wlf,file='wlf.dat',status='old',position='append')
             else
                open(unit=wlf,file='wlf.dat',status='new')
             end if

             write(wlf,'(I10,E20.12)')mc_cycle_num,wl_factor

             close(wlf)

             write(wlstring,'(I20.20)')mc_cycle_num
             open(unit=wgt,file="eta_weights.dat_"//adjustl(wlstring),status='replace',iostat=ierr)
             if (ierr/=0) stop 'Error opening eta_weights.dat for output'
             write(wgt,'("#Current energy increment = ",E20.12)')wl_factor

             open(unit=his,file="histogram.dat_"//adjustl(wlstring),status='replace',iostat=ierr)
             if (ierr/=0) stop 'Error opening histogram.dat for output'
             write(his,'("#Current energy increment = ",E20.12)')wl_factor

             do k=1,nbins
                write(wgt,*)mu_bin(k),weight(k)
                write(his,*)mu_bin(k),histogram(k)
             end do
             close(wgt)
             close(his)

          end if ! if (mod_cycle_num,invt_dump_int)

       end if !  if invT active of Swetnam

       return

  end subroutine mc_check_flatness

  integer function mu_to_bin(mu_in)
    !--------------------------------------!
    ! Takes a value of mu and returns the  !
    ! index of the bin in which it resides !
    !--------------------------------------!
    use userparams, only : nbins
    implicit none
    real(kind=dp),intent(in) :: mu_in
    real(kind=dp) :: arg
    integer       :: Ns
    
    ! Central bin in middle
    if (abs(mu_in)<=0.5_dp) then
       Ns = nbins/2+1
       mu_to_bin = Ns
       return
    end if

    if (mu_in>0.0_dp) then
       arg   = 1.0_dp-(mu_in-0.5_dp)*(1.0_dp-r_pos)/a_pos
       Ns    = nbins/2 + 2 + int(log(arg)/log(r_pos))
    else
       arg   = 1.0_dp-(abs(mu_in)-0.5_dp)*(1.0_dp-r_neg)/a_neg
       Ns    = nbins/2 - int(log(arg)/log(r_neg))
    end if

    mu_to_bin = Ns

  end function mu_to_bin

  subroutine mc_check_chain_synchronisation()
    !--------------------------------------------------------!
    ! Make sure that the two Markov chains are perfectly in  !
    ! sync, and correct that in lattice 2 to match that in   !
    ! lattice 1.                                             !
    !--------------------------------------------------------!
    ! D. Quigley February 2011                               !
    !--------------------------------------------------------!
    use io,         only : mylog
    use constants,  only : invPi,hart_to_ev,kB
    use userparams, only : nwater,pressure,num_lattices,temperature,leshift
    use util,       only : util_recipmatrix,util_determinant
    use model,      only : hmatrix,ref_hmatrix,ljr,ref_ljr,recip_matrix,volume
    use energy,     only : model_energy,compute_model_energy,compute_ivects
    implicit none

    real(kind=dp),dimension(2)     :: current_energy
    real(kind=dp),dimension(3,3,2) :: hmat_diff
    real(kind=dp),allocatable,dimension(:,:,:) :: spos_diff
    real(kind=dp),dimension(3,2)   :: svect,ref_svect
    real(kind=dp) :: beta

    integer :: ils,iwater,ierr

    beta = 1.0_dp/(kB*temperature)

    write(mylog,'("#                                                              #")')  
    write(mylog,'("# Enforcing chain synchronisation                              #")')
    write(mylog,'("#--------------------------------------------------------------#")')  
    write(mylog,'("#                                                              #")')  


    ! Compute the energy using stored coordinates
    call compute_model_energy(1)
    call compute_model_energy(2)
    current_energy = model_energy

    ! Compute the current overlap parameter using these energies
    ls_mu     = model_energy(1) + pressure*volume(1) - model_energy(2) - pressure*volume(2) 
    if (leshift) ls_mu = ls_mu - ref_enthalpy(1) + ref_enthalpy(2)
    ls_mu     = ls_mu*beta - real(Nwater,kind=dp)*log(volume(1)/volume(2))
    write(mylog,'("# Original  overlap parameter : ",F15.6,"                #")')ls_mu

    ! Compare lattice vectors to reference
    do ils = 1,num_lattices ! loop over lattices
       hmat_diff(:,:,ils) = hmatrix(:,:,ils) - ref_hmatrix(:,:,ils)
    end do

!!$    write(mylog,*)
!!$    write(mylog,'("Cell vector displacements")')
!!$    write(mylog,'("=========================")')
!!$    write(mylog,*)
!!$    write(mylog,'(22X,"Lattice 1",36X,"Lattice 2",32X,"Difference")')
!!$    write(mylog,'(22X,"---------",36X,"---------",32X,"----------")')
!!$    write(mylog,*)
!!$    write(mylog,'(6F15.6,3E15.6)')hmat_diff(:,1,1),hmat_diff(:,1,2),hmat_diff(:,1,2)-hmat_diff(:,1,1)
!!$    write(mylog,'(6F15.6,3E15.6)')hmat_diff(:,2,1),hmat_diff(:,2,2),hmat_diff(:,2,2)-hmat_diff(:,2,1)
!!$    write(mylog,'(6F15.6,3E15.6)')hmat_diff(:,3,1),hmat_diff(:,3,2),hmat_diff(:,3,2)-hmat_diff(:,3,1)

    ! Correct hmatrix in 2nd lattice to have same hmat_diff as first.
    hmatrix(:,:,2) = ref_hmatrix(:,:,2) + hmat_diff(:,:,1)

    ! Update matrix of reciprocal lattice vectors
    call util_recipmatrix(hmatrix(:,:,1),recip_matrix(:,:,1))
    call util_recipmatrix(hmatrix(:,:,2),recip_matrix(:,:,2))

    ! Allocate array to hold displacements from reference
    ! scaled molecular positions.
    allocate(spos_diff(1:3,1:nwater,1:2),stat=ierr)
    if (ierr/=0) stop 'Error allocating spos_diff'

    ! Loop over water molecules
    do iwater = 1,nwater

       ! Loop over lattices
       do ils = 1,num_lattices

          ! Current scaled coordinates
          svect(1,ils) = recip_matrix(1,1,ils)*ljr(1,1,iwater,ils) + &
                         recip_matrix(2,1,ils)*ljr(2,1,iwater,ils) + &
                         recip_matrix(3,1,ils)*ljr(3,1,iwater,ils)
          svect(2,ils) = recip_matrix(1,2,ils)*ljr(1,1,iwater,ils) + &
                         recip_matrix(2,2,ils)*ljr(2,1,iwater,ils) + &
                         recip_matrix(3,2,ils)*ljr(3,1,iwater,ils)  
          svect(3,ils) = recip_matrix(1,3,ils)*ljr(1,1,iwater,ils) + &
                         recip_matrix(2,3,ils)*ljr(2,1,iwater,ils) + &
                         recip_matrix(3,3,ils)*ljr(3,1,iwater,ils) 
          
          svect(1,ils) = svect(1,ils)*0.5_dp*invPi 
          svect(2,ils) = svect(2,ils)*0.5_dp*invPi 
          svect(3,ils) = svect(3,ils)*0.5_dp*invPi 

          ! Current scaled coordinates
          ref_svect(1,ils) = recip_matrix(1,1,ils)*ref_ljr(1,1,iwater,ils) + &
                             recip_matrix(2,1,ils)*ref_ljr(2,1,iwater,ils) + &
                             recip_matrix(3,1,ils)*ref_ljr(3,1,iwater,ils)
          ref_svect(2,ils) = recip_matrix(1,2,ils)*ref_ljr(1,1,iwater,ils) + &
                             recip_matrix(2,2,ils)*ref_ljr(2,1,iwater,ils) + &
                             recip_matrix(3,2,ils)*ref_ljr(3,1,iwater,ils)  
          ref_svect(3,ils) = recip_matrix(1,3,ils)*ref_ljr(1,1,iwater,ils) + &
                             recip_matrix(2,3,ils)*ref_ljr(2,1,iwater,ils) + &
                             recip_matrix(3,3,ils)*ref_ljr(3,1,iwater,ils) 
          
          ref_svect(1,ils) = ref_svect(1,ils)*0.5_dp*invPi 
          ref_svect(2,ils) = ref_svect(2,ils)*0.5_dp*invPi 
          ref_svect(3,ils) = ref_svect(3,ils)*0.5_dp*invPi 

          ! Displacement of scaled coordinates from reference scaled coordinates
          spos_diff(:,iwater,ils) = svect(:,ils) - ref_svect(:,ils)

       end do

       ! Force scaled displacements in lattice 2 to match those in lattice 1
       svect(:,2) = ref_svect(:,2) + spos_diff(:,iwater,1)
       ljr(:,1,iwater,2) = matmul(hmatrix(:,:,2),svect(:,2))

    end do

!!$    write(mylog,*)
!!$    write(mylog,'("Fractional position displacements")')
!!$    write(mylog,'("=================================")')
!!$    write(mylog,*)
!!$    write(mylog,'(22X,"Lattice 1",36X,"Lattice 2",32X,"Difference")')
!!$    write(mylog,'(22X,"---------",36X,"---------",32X,"----------")')
!!$    write(mylog,*)
!!$    do iwater = 1,nwater
!!$       write(*,'(6F15.6,3E15.6)')spos_diff(:,iwater,1),spos_diff(:,iwater,2), &
!!$                                 spos_diff(:,iwater,2) - spos_diff(:,iwater,1)
!!$    end do

    ! Release memory
    deallocate(spos_diff,stat=ierr)
    if (ierr/=0) stop 'Error deallocating spos_diff'

!!$    write(mylog,*)
!!$    write(mylog,'("RMS Rotational error")')
!!$    write(mylog,'("====================")')
!!$    write(mylog,*)
!!$    write(mylog,'("    Lattice 1        Lattice 2")')
!!$    write(mylog,'("    =========        =========")')

!!$    ! Compute the difference between the stored charge site coordinates and those
!!$    ! obtained by applying the stored quaternion to the reference charge sites.
!!$    do iwater = 1,nwater
!!$       do ils = 1,num_lattices 
!!$          rot_err(ils) = 0.0_dp
!!$          do ics = 1,cspm
!!$             tmpvect(:) = quatconj(quat(:,iwater),ref_csr(:,ics,iwater,ils)-ref_ljr(:,1,iwater,ils))
!!$             tmpvect(:) = tmpvect(:) - (csr(:,ics,iwater,ils)-ljr(:,1,iwater,ils))
!!$             rot_err(ils) = rot_err(ils) + dot_product(tmpvect,tmpvect)
!!$          end do
!!$          rot_err(ils) = rot_err(ils)/3.0_dp
!!$          rot_err(ils) = sqrt(rot_err(ils))
!!$       end do

!!$       write(mylog,'(3E15.6)')rot_err(:)

       ! Enforce synchoronisation by making sure this difference is zero in both lattices
!!$       do ils = 1,num_lattices
!!$          do ics = 1,cspm
!!$             csr(:,ics,iwater,ils) = quatconj(quat(:,iwater),ref_csr(:,ics,iwater,ils)-ref_ljr(:,1,iwater,ils))
!!$             csr(:,ics,iwater,ils) =  csr(:,ics,iwater,ils) + ljr(:,1,iwater,ils)
!!$          end do
!!$       end do

!!$    end do

    ! Compute reciprocal and real space lattice vectors, charge density
    ! and cached exponentials.
    do ils = 1,num_lattices
       volume(ils) = abs(util_determinant(hmatrix(:,:,ils)))
!       call compute_kvects(ils)
       call compute_ivects(ils)
!       call compute_rho_k(ils)
!       call compute_expor(ils)
    end do
    
    ! Compute new model energy
    call compute_model_energy(1)
    call compute_model_energy(2)


    ! Compute new overlap parameter
    ls_mu     = model_energy(1) + pressure*volume(1) - model_energy(2) - pressure*volume(2) 
    if (leshift) ls_mu = ls_mu - ref_enthalpy(1) + ref_enthalpy(2)
    ls_mu     = ls_mu*beta - real(Nwater,kind=dp)*log(volume(1)/volume(2))

    write(mylog,'("# Corrected overlap parameter : ",F15.6,"                #")')ls_mu

    write(mylog,'("#                                                              #")')  
    write(mylog,'("# Change in energy after synchronisation                       #")')
    write(mylog,'("# --------------------------------------                       #")')
    write(mylog,'("#                                                              #")')  
    write(mylog,'("# lattice ",I6," : ",E15.6,"                             #")')1,(model_energy(1)-current_energy(1))*hart_to_ev
    write(mylog,'("# lattice ",I6," : ",E15.6,"                             #")')2,(model_energy(2)-current_energy(2))*hart_to_ev
    write(mylog,'("#                                                              #")')  

    return

  end subroutine mc_check_chain_synchronisation

!!$  subroutine mc_compute_deltaG_from_eta()
!!$    !--------------------------------------------------------!
!!$    ! Assumes that the current set of multicanonical weights !
!!$    ! are perfect, i.e. exactly -G(mu), and computes the     !
!!$    ! free energy difference between the two lattices.       !
!!$    !--------------------------------------------------------!
!!$    ! D. Quigley March 2011                                  !
!!$    !--------------------------------------------------------!
!!$    use comms,      only : myrank,comms_allreduce_eta
!!$    use io,         only : glog
!!$    use constants,  only : hart_to_kJpm,Kb,hart_to_eV
!!$    use userparams, only : temperature,nbins,nwater,leshift, &
!!$                           parallel_strategy
!!$    implicit none
!!$    real(kind=dp),allocatable,dimension(:) :: tmpP
!!$    real(kind=dp) :: Pnorm,pA,pB,deltaG,beta
!!$    integer :: i,ierr
!!$
!!$    ! MPI - allreduce on eta
!!$    if ( parallel_strategy == 'mw' ) then
!!$       call comms_allreduce_eta(weight,nbins)
!!$    elseif ( parallel_strategy == 'dd' ) then
!!$      
!!$    else
!!$       ! Unknown parallel strategy
!!$    end if
!!$
!!$    allocate(tmpP(1:nbins),stat=ierr)
!!$    if (ierr/=0) stop 'Error allocating tmpP in mc_compute_deltaG_from_eta'
!!$
!!$    ! Compute an estimate for P(mu)
!!$    Pnorm = 0.0_dp
!!$    do i = 1,nbins
!!$       tmpP(i) = exp(weight(i))
!!$       Pnorm   = Pnorm + tmpP(i)*binwidth(i)
!!$    end do
!!$    do i = 1,nbins
!!$       tmpP(i) = tmpP(i)/Pnorm
!!$    end do
!!$
!!$    ! Integrate up to nbins/2
!!$    pA = 0.0_dp
!!$    do i = 1,nbins/2
!!$       pA = pA + tmpP(i)*0.5_dp*(binwidth(i)+binwidth(i+1))
!!$       pA = pA + 0.5_dp*binwidth(i)*(tmpP(i+1)-tmpP(i))
!!$    end do
!!$
!!$    ! Integrate the rest of the way
!!$    pB = 0.0_dp
!!$    do i = nbins/2+1,nbins
!!$       pB = pB + 0.5_dp*binwidth(i-1)*(tmpP(i-1)-tmpP(i))  
!!$       pB = pB + tmpP(i)*0.5_dp*(binwidth(i-1)+binwidth(i))
!!$    end do
!!$
!!$    ! Free energy difference in units of kT
!!$    ! negative if pA > pB, so energy A->B
!!$    ! A is lattice 1, B is lattice 2
!!$    deltaG = log(pA/pB)
!!$    beta = 1.0_dp/(Kb*temperature)
!!$    if (leshift) deltaG = deltaG + beta*ref_enthalpy(2) - beta*ref_enthalpy(1)
!!$
!!$    if (myrank==0) then
!!$       write(glog,'("#                                                              #")')  
!!$       write(glog,'("# Estimate of delta G direct from weights at cycle ",I10,"  #")')mc_cycle_num
!!$       write(glog,'("#------------------------------------------------------------  #")')
!!$       write(glog,'("#                                                              #")')  
!!$       write(glog,'("# G(lattice2) - G(lattice1) = ",F15.8, " kT/molecule      #")')deltaG/real(nwater,kind=dp)
!!$       write(glog,'("# G(lattice2) - G(lattice1) = ",F15.8, " J/mole           #")') &
!!$             Kb*temperature*hart_to_kJpm*1000.0_dp*deltaG/real(nwater,kind=dp)
!!$       write(glog,'("# G(lattice2) - G(lattice1) = ",F15.8, " meV/molecule     #")') &
!!$             Kb*temperature*hart_to_eV*1000.0_dp*deltaG/real(nwater,kind=dp) 
!!$       write(glog,'("#                                                              #")')  
!!$    end if
!!$
!!$    deallocate(tmpP)
!!$
!!$    return
!!$
!!$  end subroutine mc_compute_deltaG_from_eta

  subroutine mc_compute_deltaG_from_hist()
    !--------------------------------------------------------!
    ! Compute the current estimate of delta G by unweighting !
    ! the current histogram using the current weights, and   !
    ! computing the area under each segment. Most useful for !
    ! long runs in which the weights are static and the      !
    ! histogram never resets.                                !
    !--------------------------------------------------------!
    ! D. Quigley March 2011                                  !
    !--------------------------------------------------------!
    use comms,      only : myrank,comms_allreduce_hist,comms_allreduce_eta, &
                           comms_allreduce_uhist, comms_join_uhist
    use io,         only : glog 
    use constants,  only : hart_to_kJpm,Kb,hart_to_ev
    use userparams, only : temperature,nbins,nwater,leshift,samplerun, &
                           parallel_strategy,window_overlap
    implicit none
    real(kind=dp),allocatable,dimension(:) :: normP
    real(kind=dp) :: Pnorm,pA,pB,deltaG,beta
    integer :: i,ierr
    character(10) :: cyclestring
    character(33) :: filename

    real(kind=dp),allocatable,dimension(:) :: joined

!!$    ! MPI - all reduce on eta and histogram
!!$    call comms_allreduce_eta(weight,nbins)
!!$    call comms_allreduce_hist(histogram,nbins)

    allocate(joined(1:nbins),stat=ierr)
    if (ierr/=0) stop 'Error allocating joined histogram in mc_compute_deltaG_from_hist'
    
    if (samplerun) then
       if (parallel_strategy=='mw') then
          call comms_allreduce_uhist(unbiased_hist,nbins)
          joined = unbiased_hist
       elseif (parallel_strategy=='dd') then
          call comms_join_uhist(unbiased_hist,nbins,window_overlap,joined)          
       else
          stop 'Error in mc_compute_deltaG_from_hist - urecognised parallel strategy'
       end if
    end if

!!$    allocate(tmpP(1:nbins),stat=ierr)
!!$    if (ierr/=0) stop 'Error allocating tmpP in mc_compute_deltaG_from_hist'
    allocate(normP(1:nbins),stat=ierr)
    if (ierr/=0) stop 'Error allocating normP in mc_compute_deltaG_from_hist'
    
    ! Compute normalised P as the current unbiased histogram
    Pnorm = 0.0_dp
    do i = 1,nbins
       Pnorm = Pnorm + joined(i)*binwidth(i)
    end do
    do i = 1,nbins
       normP(i) = joined(i)/Pnorm
    end do    

!!$    ! Compute an estimate for P(mu)
!!$    Pnorm = 0.0_dp
!!$    do i = 1,nbins
!!$       tmpP(i) = normP(i)*exp(weight(i))
!!$       Pnorm   = Pnorm + tmpP(i)*binwidth(i)
!!$    end do
!!$    do i = 1,nbins
!!$       tmpP(i) = tmpP(i)/Pnorm
!!$    end do

    ! Sum up to nbins/2
    pA = 0.0_dp
    do i = 1,nbins/2
!!$       pA = pA + P(i)*0.5_dp*(binwidth(i)+binwidth(i+1))
!!$       pA = pA + 0.5_dp*binwidth(i)*(tmpP(i+1)-tmpP(i))
       pA = pA + normP(i)*binwidth(i)
    end do

    ! Sum the rest of the way
    pB = 0.0_dp
    do i = nbins/2+1,nbins
!!$       pB = pB + 0.5_dp*binwidth(i-1)*(tmpP(i-1)-tmpP(i))
!!$       pB = pB + tmpP(i)*0.5_dp*(binwidth(i-1)+binwidth(i))
       pB = pB + normP(i)*binwidth(i)
    end do

    ! Free energy difference in units of kT
    ! negative if pA > pB, so energy A->B
    ! A is lattice 1, B is lattice 2
    deltaG = log(pA/pB)
    beta = 1.0_dp/(Kb*temperature)
    if (leshift) deltaG = deltaG + beta*ref_enthalpy(2) - beta*ref_enthalpy(1)    

    if (myrank == 0) then

       write(cyclestring,'(I10.10)')mc_cycle_num
       filename = 'unbiased_histogram_'//cyclestring//'.dat'
       
       write(glog,'("#                                                              #")')  
       write(glog,'("# Estimate of delta G from histogram at cycle      ",I10,"  #")')mc_cycle_num
       write(glog,'("#------------------------------------------------------------  #")')
       write(glog,'("#                                                              #")')  
       write(glog,'("# G(lattice2) - G(lattice1) = ",F15.8, " kT/molecule      #")')deltaG/real(nwater,kind=dp)
       write(glog,'("# G(lattice2) - G(lattice1) = ",F15.8, " J/mole           #")') &
             Kb*temperature*hart_to_kJpm*1000.0_dp*deltaG/real(nwater,kind=dp)
       write(glog,'("# G(lattice2) - G(lattice1) = ",F15.8, " meV/molecule     #")') &
             Kb*temperature*hart_to_eV*1000.0_dp*deltaG/real(nwater,kind=dp) 
       write(glog,'("#                                                              #")')  
       
       call flush(glog)

       open(unit=ubh,file=trim(filename),status='replace',iostat=ierr)
       if (ierr/=0) stop 'Error opening file for unweighted histogram'
       
       do i=1,nbins
          write(ubh,*)mu_bin(i),normP(i)
       end do
       
       close(ubh)

    end if

    deallocate(normP,joined)

    return

  end subroutine mc_compute_deltaG_from_hist

  function quatmul_norm(a,b)

    implicit none
    real(kind=dp),dimension(4) :: quatmul_norm
    real(kind=dp),dimension(4),intent(in) :: a,b
    real(kind=dp),dimension(4) :: tmpquat
    real(kind=dp) :: tnorm
    
    tmpquat(1)   = a(1)*b(1) - a(2)*b(2) - a(3)*b(3) - a(4)*b(4) 
    tmpquat(2:4) = a(1)*b(2:4) + b(1)*a(2:4) + cross_product(a(2:4),b(2:4)) 

    tnorm   = sqrt(dot_product(tmpquat,tmpquat))
    quatmul_norm =  tmpquat/tnorm

    return

  end function quatmul_norm

  function quatmul(a,b)

    implicit none
    real(kind=dp),dimension(4) :: quatmul
    real(kind=dp),dimension(4),intent(in) :: a,b
    real(kind=dp),dimension(4) :: tmpquat
    
    tmpquat(1)   = a(1)*b(1) - a(2)*b(2) - a(3)*b(3) - a(4)*b(4) 
    tmpquat(2:4) = a(1)*b(2:4) + b(1)*a(2:4) + cross_product(a(2:4),b(2:4)) 

    quatmul = tmpquat

    return

  end function quatmul


  function quatinv(a)

    implicit none
    real(kind=dp),dimension(4) :: quatinv
    real(kind=dp),dimension(4),intent(in) :: a
    real(kind=dp),dimension(4) :: tmpquat
    
    tmpquat    = -a
    tmpquat(1) = -tmpquat(1)
    quatinv    = tmpquat

    ! N.B. assumes a normalised

    return

  end function quatinv

  function quatconj(a,v)

    implicit none
    real(kind=dp),dimension(3) :: quatconj
    real(kind=dp),dimension(4),intent(in) :: a
    real(kind=dp),dimension(3),intent(in) :: v
    real(kind=dp),dimension(4) :: q,b

    q(1)   = 0.0_dp
    q(2:4) = v 

    b = quatinv(a)
    q = quatmul(q,b)
    q = quatmul(a,q)
    
    quatconj = q(2:4)

    return

  end function quatconj

  function cross_product(a,b)

    implicit none
    real(kind=dp),dimension(3),intent(in) :: a,b
    real(kind=dp),dimension(3) :: cross_product


    cross_product(1) = a(2)*b(3) - a(3)*b(2)
    cross_product(2) = -(a(1)*b(3)-a(3)*b(1))
    cross_product(3) = a(1)*b(2)-b(1)*a(2)

  end function cross_product


end module mc_moves
