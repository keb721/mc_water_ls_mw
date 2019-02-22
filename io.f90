! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: Git -*-
!=============================================================================!
!                            I O                                              !
!=============================================================================!
! Handles the reading of the input file and the dumping of configurations     !
! to DLPOLY HISTORY files.                                                    !
!=============================================================================!

module io

  Use constants, Only : dp,int32                !Minimal useage where practical
  use omp_lib

  Implicit None                                 !Impose strong typing
  Private   

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
                                                      !... unless exposed here.
  public :: io_read_input
  public :: io_hist_header
  public :: io_hist_append
  public :: io_write_global_header
  public :: io_write_rank_header

  public :: io_write_psf
  public :: io_write_dcd_header
  public :: io_write_dcd_snapshot

  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  public :: glog,mylog,mytherm              ! global and per rank log files


  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!


  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  character(10),save :: hist_filename        ! filename for history file
  integer,save       :: hist  = 33           ! unit number for history file
  integer,save       :: inp   = 34           ! unit number for input file
  integer,save       :: glog  = 35           ! unit number for global log file
  integer,save       :: mylog = 36           ! unit number for rank log file
  integer,save       :: mytherm = 71         ! unit number for therm file
  integer,save       :: psf   = 81           ! unit number for psf file
  integer,save       :: dcd   = 82           ! unit number for dcd file

contains
 

  subroutine io_read_input()
    !------------------------------------------------------------------------------!
    ! Reads the input file specified on the command line and checks that all inputs!
    ! are valid. To be called after initialisation of communications.              !
    !------------------------------------------------------------------------------!
    ! D.Quigley October 2006                                                       !
    !------------------------------------------------------------------------------!
    use constants, only : ang_to_bohr
    use comms,     only : myrank,comms_bcastint,comms_bcastchar,comms_bcastreal, &
         comms_bcastlog,size
    use model,     only : ls
    use userparams
    use energy
    use timer
    implicit none

    ! command line data
    integer                       :: iarg,idata,ierr
    integer                       :: num_args
    character(30), dimension(0:10):: command_line
    character(30)                 :: file_name
    integer(kind=int32)           :: last_dot=30
!!$    integer,  external ::  iargc
!!$    external  getarg

    ! declare namelists
    namelist/potential/model_type,mw_sigma,mw_epsilon,mw_lambda,sw_bigA,sw_B,sw_gamma,sw_a,sw_p,sw_q,cos0

    namelist/thermal/temperature,pressure

    namelist/MonteCarlo/mc_ensemble,mc_max_trans,mc_dv_max,mc_target_ratio,   &
         wl_factor,wl_schedule,wl_flattol,wl_minhist,wl_useinvt,wl_swetnam, &
         wl_alpha,eta_interp,samplerun, &
         nbins,mu_max,mu_min,&
         allow_switch,allow_vol,allow_trans, &
         mc_trans_prob,mc_vol_prob,mc_switch_prob,mc_always_switch,leshift, &
         input_ref_enthalpy

    namelist/config/nwater,num_lattices,method,r_overlap,ls

    namelist/bookkeeping/list_update_int,traj_output_int,file_output_int,latt_sync_int,mpi_sync_int, &
         chkpt_dump_int,monitor_int,flat_chk_int,invt_dump_int,eq_adjust_mc,deltaG_int, &
         max_mc_cycles,eq_mc_cycles,timer_qtime,timer_closetime

    namelist/parallelisation/parallel_strategy,window_overlap

       
    
    if ( myrank == 0 ) then

       ! check that there is only one argument.
       num_args = iargc()

       if (num_args<1) then
          write(0,*)
          write(0,*) '           M C _ W A T E R _ L S _ M W     '
          write(0,*)
          write(0,*) '         Usage: mw_water_ls <input file>   '
          write(0,*)
          write(0,*) '        D. Quigley - University of Warwick '
          write(0,*)
          stop
       end if

       do iarg = 1, num_args
          call getarg(iarg,command_line(iarg))
       end do

       ! find the last dot in the filename.
       do idata = 1, len(seedname)
          if(command_line(1)(idata:idata)=='.') last_dot = idata
       end do

       ! set the seedname.
       seedname = command_line(1)(1:last_dot-1)  ! index() intrinsic replaced.

       !-----------------------------------------------------------!
       ! Open the input file.                                      !
       !-----------------------------------------------------------!

       ! set the file name
       file_name = trim(command_line(1))

       ! open it
       open (inp,file=file_name,status='old',iostat=ierr)
       if(ierr/=0) then
          stop 'Unable to open input file.'
       end if

       !-----------------------------------------------------------!
       ! Potential                                                 !
       !-----------------------------------------------------------!

       ! read the potential namelists
       read(inp,nml=potential,iostat=ierr)
       if (ierr/=0) stop 'Error reading potential namelist'

       !-------------------------!
       !   Thermal Parameters    !
       !-------------------------!    
       read(inp,nml=thermal,iostat=ierr)
       if (ierr/=0) stop 'Error reading thermal namelist'

       ! sanity cheks
       if ( temperature < 0.0_dp ) stop 'Error temperature must be positive'

       ! convert into internal representation
       pressure = pressure/aup_to_atm

       !-------------------------!
       !  MC Parameters          !
       !-------------------------!    

       mc_switch_prob  = 0.1_dp
       mc_vol_prob     = 1.0_dp/real(nwater,kind=dp)

       read(inp,nml=MonteCarlo,iostat=ierr)    
       if (ierr/=0) stop 'Error reading MonteCarlo namelist'

       ! sanity checks
       if ( (mc_ensemble/='nvt').and.(mc_ensemble)/='npt' ) then
          write(0,*)'Error - unrecognised ensemble.'
          write(0,*)' Please choose either npt or nvt. '
          stop
       end if

       ! convert to internal units
       mc_max_trans = mc_max_trans*ang_to_bohr
       mc_dv_max    = mc_dv_max*ang_to_bohr

       !-------------------------!
       !  Configuration          !
       !-------------------------! 
       read(inp,nml=config,iostat=ierr)
       if (ierr/=0) stop 'Error reading config namelist'

       ! sanity checks
       if (nwater<1 ) stop 'Error - invalid number of waters'
       if (r_overlap < 0.0_dp ) stop 'Error - invalid r_overlap'
       if ( (trim(method) /= 'xmol') ) then
          write(0,*)'Invalid initialisation option in config namelist'
          write(0,*)
          write(0,*)"           use method=\'xmol\' only             "
          stop
       end if


       ! convert to internal representation
       r_overlap = r_overlap*ang_to_bohr

       select case(num_lattices)
          case(1)
             ! Single box mode - override all lattice switch settings
             allow_switch = .false.
             mc_switch_prob = 0.0_dp
             mc_always_switch = .false.
             ls = 1 ! Start (and stay for all time) in lattice 1
          case(2)
             ! Do nothing
          case default
             write(0,*)'Error num_lattices must equal 1 or 2!'
             stop             
       end select

       !-------------------------!
       ! Bookkeeping             !
       !-------------------------!
       read(inp,nml=bookkeeping,iostat=ierr)
       if (ierr/=0) stop 'Error reading bookkeeping namelist'

       ! sanity checks
       if ( list_update_int < 1 ) stop 'Error - list_update_int must be > 0 '
       if ( traj_output_int < 1 ) stop 'Error - traj_output_int must be > 0 '
       if ( file_output_int < 1 ) stop 'Error - file_output_int must be > 0 '
       if ( max_mc_cycles   < 1 ) stop 'Error - max_mc_cycles must be   > 0 '
       if ( eq_mc_cycles    < 1 ) stop 'Error - eq_mc_cycles must be    > 0'

       !-------------------------!
       ! Parallelisation         !
       !-------------------------!
       read(inp,nml=parallelisation,iostat=ierr)
       if ((ierr/=0).and.(size>1)) then
          !write(0,'("Parallelisation namelist missing or contains invalid entries")')
          write(0,'("Using default multi-walker parallelisation.")')
       end if

       close(inp)

    end if

    ! Broadcast all input parameters
    ! grep '::' userparams.f90 | sed 's/,save//' | sed 's/(kind=dp)//' | awk '{print $1,$3}' > dv.dat

    ! awk '{if ($1=="integer") print "call comms_bcastint("$2",1)"}' dv.dat 
    call comms_bcastint(nwater,1)                         
    call comms_bcastint(num_lattices,1)
    call comms_bcastint(wl_schedule,1)
    call comms_bcastint(wl_minhist,1)
    call comms_bcastint(nbins,1)
    call comms_bcastint(list_update_int,1)
    call comms_bcastint(traj_output_int,1)
    call comms_bcastint(file_output_int,1)
    call comms_bcastint(latt_sync_int,1)
    call comms_bcastint(mpi_sync_int,1)
    call comms_bcastint(chkpt_dump_int,1)
    call comms_bcastint(monitor_int,1)
    call comms_bcastint(flat_chk_int,1)
    call comms_bcastint(invt_dump_int,1)
    call comms_bcastint(deltaG_int,1)
    call comms_bcastint(max_mc_cycles,1)
    call comms_bcastint(eq_mc_cycles,1)
    call comms_bcastint(num_lattices,1)

    !awk '{if ($1=="real") print "call comms_bcastreal("$2",1)"}' dv.dat
    call comms_bcastreal(r_overlap,1)
    call comms_bcastreal(pressure,1)
    call comms_bcastreal(temperature,1)
    call comms_bcastreal(mc_max_trans,1)
    call comms_bcastreal(mc_target_ratio,1)
    call comms_bcastreal(mc_dv_max,1)
    call comms_bcastreal(wl_factor,1)
    call comms_bcastreal(wl_flattol,1)
    call comms_bcastreal(wl_alpha,1)
    call comms_bcastreal(mu_min,1)
    call comms_bcastreal(mu_max,1)
    call comms_bcastreal(mc_trans_prob,1)
    call comms_bcastreal(mc_vol_prob,1)
    call comms_bcastreal(mc_switch_prob,1)
    call comms_bcastreal(input_ref_enthalpy(1),2)

    !awk '{if ($1=="logical") print "call comms_bcastlog("$2",1)"}' dv.dat
    call comms_bcastlog(wl_useinvt,1)
    call comms_bcastlog(wl_swetnam,1)
    call comms_bcastlog(eta_interp,1)
    call comms_bcastlog(samplerun,1)
    call comms_bcastlog(allow_switch,1)
    call comms_bcastlog(allow_vol,1)
    call comms_bcastlog(allow_trans,1)
    call comms_bcastlog(mc_always_switch,1)
    call comms_bcastlog(eq_adjust_mc,1)
    call comms_bcastlog(leshift,1)

    ! grep character dv.dat | sed 's/(/ /g' | sed 's/)/ /g' | awk '{print "call comms_bcastchar("$3","$2")"}'
    call comms_bcastchar(model_type,12)
    call comms_bcastchar(method,6)
    call comms_bcastchar(mc_ensemble,3)
    call comms_bcastchar(seedname,30)


    call comms_bcastreal(mw_sigma,1)
    call comms_bcastreal(mw_epsilon,1)
    call comms_bcastreal(mw_lambda,1)

    call comms_bcastreal(sw_bigA,1)
    call comms_bcastreal(sw_B,1)
    call comms_bcastreal(sw_gamma,1)
    call comms_bcastreal(sw_a,1)
    call comms_bcastint(sw_p,1)
    call comms_bcastint(sw_q,1)
    call comms_bcastreal(cos0,1)

    call comms_bcastchar(parallel_strategy,2)
    call comms_bcastint(window_overlap,1)
    
    return

  end subroutine io_read_input

!!$  subroutine io_hist_header()
!!$    !------------------------------------------------------------------------------!
!!$    ! Writes the header of a DLPOLY HISTORY file. This is appended to at intervals !
!!$    ! in the calculation.                                                          !
!!$    !------------------------------------------------------------------------------!
!!$    ! D.Quigley October 2006                                                       !
!!$    !------------------------------------------------------------------------------!
!!$    use userparams, only : model_type,nwater
!!$    use model,      only : mass,charge
!!$    implicit none
!!$    character(80) :: header
!!$    character(3)  :: rankstr
!!$    real(8)       :: natms
!!$    character(8),allocatable,dimension(:)  :: atname
!!$    real(8),allocatable,dimension(:)       :: weight
!!$    real(8),allocatable,dimension(:)       :: ocharge
!!$
!!$    integer :: ierr,imol,n,k
!!$    integer :: myrank = 0
!!$
!!$    write(hist_filename,'("HISTORY",I3.3)')myrank
!!$    open(unit=hist,file=hist_filename,status='replace',form='unformatted',iostat=ierr)
!!$    if (ierr/=0) stop 'Error opening HISTORY file'
!!$
!!$    write(rankstr,'(I3)')myrank
!!$    header="Metadynamics + Parallel tempering Monte-Carlo. HISTORY file for rank "//rankstr
!!$
!!$    select case (model_type)
!!$    case ("tip4p/2005")
!!$       natms = real(4*nwater,kind=8)
!!$       !print*,myrank,nwater
!!$    case ("tip4p")
!!$       natms = real(4*nwater,kind=8)
!!$    case ("tip4p/Ice")
!!$       natms = real(4*nwater,kind=8)
!!$    case default
!!$       write(0,*)'HISTORY file output for selected model not implemented'
!!$    end select
!!$
!!$    n = anint(natms)
!!$
!!$    allocate(atname(1:n),stat=ierr)
!!$    allocate(weight(1:n),stat=ierr)
!!$    allocate(ocharge(1:n),stat=ierr)
!!$    if (ierr/=0) stop 'Error allocating memory in io_hist_header'
!!$
!!$    k = 1
!!$    do imol = 1,nwater
!!$
!!$       select case (trim(model_type))
!!$       case ("tip4p/2005")
!!$
!!$          atname(k)   = 'O_tip4p'
!!$          atname(k+1) = 'H_tip4p'
!!$          atname(k+2) = 'H_tip4p'
!!$          atname(k+3) = 'M_tip4p'
!!$
!!$          ocharge(k)   = 0.0_dp
!!$          ocharge(k+1) = charge(2)
!!$          ocharge(k+2) = charge(3)
!!$          ocharge(k+3) = charge(1)
!!$
!!$
!!$          weight(k)    = 15.9998_dp
!!$          weight(k+1)  = 1.0080_dp
!!$          weight(k+2)  = 1.0080_dp
!!$          weight(k+3)  = 0.0_dp
!!$
!!$          k = k + 4
!!$
!!$       case("tip4p")
!!$
!!$          atname(k)   = 'O_tip4p'
!!$          atname(k+1) = 'H_tip4p'
!!$          atname(k+2) = 'H_tip4p'
!!$          atname(k+3) = 'M_tip4p'
!!$
!!$          ocharge(k)   = 0.0_dp
!!$          ocharge(k+1) = charge(2)
!!$          ocharge(k+2) = charge(3)
!!$          ocharge(k+3) = charge(1)
!!$
!!$
!!$          weight(k)    = 15.9998_dp
!!$          weight(k+1)  = 1.0080_dp
!!$          weight(k+2)  = 1.0080_dp
!!$          weight(k+3)  = 0.0_dp
!!$
!!$          k = k + 4
!!$
!!$       case("tip4p/Ice")
!!$
!!$          atname(k)   = 'O_tip4p'
!!$          atname(k+1) = 'H_tip4p'
!!$          atname(k+2) = 'H_tip4p'
!!$          atname(k+3) = 'M_tip4p'
!!$
!!$          ocharge(k)   = 0.0_dp
!!$          ocharge(k+1) = charge(2)
!!$          ocharge(k+2) = charge(3)
!!$          ocharge(k+3) = charge(1)
!!$
!!$
!!$          weight(k)    = 15.9998_dp
!!$          weight(k+1)  = 1.0080_dp
!!$          weight(k+2)  = 1.0080_dp
!!$          weight(k+3)  = 0.0_dp
!!$
!!$          k = k + 4
!!$
!!$       case default
!!$
!!$          write(0,*)'HISTORY file output for selected model not implemented'
!!$
!!$       end select
!!$
!!$    end do
!!$
!!$
!!$    ! write the HISTORY FILE HEADER
!!$    write(hist)header
!!$    write(hist)natms
!!$    write(hist)atname
!!$    write(hist)weight
!!$    write(hist)ocharge
!!$
!!$    deallocate(atname,weight,ocharge,stat=ierr)
!!$    if (ierr/=0) stop 'Error releasing memory in io_hist_header'
!!$
!!$    close(hist)
!!$
!!$
!!$  end subroutine io_hist_header
!!$
!!$  subroutine io_hist_append(icyc)
!!$    !------------------------------------------------------------------------------!
!!$    ! Appends the current configuration to a DLPOLY HISTORY file.                  !
!!$    !------------------------------------------------------------------------------!
!!$    ! D.Quigley October 2006                                                       !
!!$    !------------------------------------------------------------------------------!
!!$    use util,       only : util_images
!!$    use constants,  only : bohr_to_ang
!!$    use userparams, only : model_type,nwater
!!$    use model,    only   : ljr,csr,hmatrix,recip_matrix,ls
!!$    implicit none
!!$
!!$    integer,intent(in) :: icyc
!!$
!!$    real(8) :: nstep
!!$    real(8) :: natms
!!$    real(8) :: keytrj
!!$    real(8) :: imcon
!!$    real(8) :: tstep
!!$
!!$    real(8),allocatable,dimension(:) :: xxx,yyy,zzz
!!$
!!$    real(8),dimension(3) :: vec_o,vec_m,vec_h1,vec_h2
!!$
!!$    integer :: imol,ierr,n,k
!!$
!!$    integer :: myrank=0
!!$
!!$    write(hist_filename,'("HISTORY",I3.3)')myrank
!!$    open(unit=hist,file=hist_filename,status='old',form='unformatted',position='append',iostat=ierr)
!!$    if (ierr/=0) stop 'Error opening HISTORY file'
!!$
!!$
!!$    select case (trim(model_type))
!!$    case ("tip4p/2005")
!!$       natms = real(4*nwater,kind=8)
!!$    case ("tip4p")
!!$       natms = real(4*nwater,kind=8)
!!$    case ("tip4p/Ice")
!!$       natms = real(4*nwater,kind=8)
!!$    case default
!!$       write(0,*)'HISTORY file output for selected model not implemented'
!!$    end select
!!$
!!$
!!$    nstep  = real(icyc,kind=8)
!!$    keytrj = 0.0_dp
!!$    imcon  = 3.0_dp
!!$    tstep  = 1.0_dp
!!$
!!$    n = anint(natms)
!!$    allocate(xxx(1:n),stat=ierr)
!!$    allocate(yyy(1:n),stat=ierr)
!!$    allocate(zzz(1:n),stat=ierr)
!!$    if (ierr/=0) stop 'Error allocating coords in io_hist_append'
!!$
!!$    k = 1
!!$    do imol = 1,nwater
!!$
!!$       select case (trim(model_type))
!!$       case ("tip4p/2005")
!!$
!!$
!!$          vec_o(1) = ljr(1,1,imol,ls)
!!$          vec_o(2) = ljr(2,1,imol,ls)
!!$          vec_o(3) = ljr(3,1,imol,ls)
!!$
!!$          vec_m(1) = csr(1,1,imol,ls)
!!$          vec_m(2) = csr(2,1,imol,ls)
!!$          vec_m(3) = csr(3,1,imol,ls)
!!$
!!$          vec_m(:) = vec_m(:) - vec_o(:)
!!$
!!$          vec_h1(1) = csr(1,2,imol,ls)
!!$          vec_h1(2) = csr(2,2,imol,ls)
!!$          vec_h1(3) = csr(3,2,imol,ls)
!!$
!!$          vec_h1(:) = vec_h1(:) - vec_o(:)          
!!$
!!$          vec_h2(1) = csr(1,3,imol,ls)
!!$          vec_h2(2) = csr(2,3,imol,ls)
!!$          vec_h2(3) = csr(3,3,imol,ls)
!!$
!!$          vec_h2(:) = vec_h2(:) - vec_o(:) 
!!$
!!$
!!$          call util_images(vec_o,hmatrix(:,:,ls),recip_matrix(:,:,ls))
!!$
!!$          vec_m(:) = vec_o(:) + vec_m(:)
!!$          
!!$          vec_h1(:) = vec_h1(:) + vec_o(:)
!!$          vec_h2(:) = vec_h2(:) + vec_o(:)          
!!$
!!$          xxx(k) = vec_o(1)*bohr_to_ang
!!$          yyy(k) = vec_o(2)*bohr_to_ang
!!$          zzz(k) = vec_o(3)*bohr_to_ang
!!$
!!$          xxx(k+1) = vec_h1(1)*bohr_to_ang
!!$          yyy(k+1) = vec_h1(2)*bohr_to_ang
!!$          zzz(k+1) = vec_h1(3)*bohr_to_ang
!!$
!!$          xxx(k+2) = vec_h2(1)*bohr_to_ang
!!$          yyy(k+2) = vec_h2(2)*bohr_to_ang
!!$          zzz(k+2) = vec_h2(3)*bohr_to_ang
!!$
!!$          xxx(k+3) = vec_m(1)*bohr_to_ang
!!$          yyy(k+3) = vec_m(2)*bohr_to_ang
!!$          zzz(k+3) = vec_m(3)*bohr_to_ang
!!$
!!$          k = k + 4
!!$
!!$
!!$       case ("tip4p/Ice")
!!$
!!$
!!$          vec_o(1) = ljr(1,1,imol,ls)
!!$          vec_o(2) = ljr(2,1,imol,ls)
!!$          vec_o(3) = ljr(3,1,imol,ls)
!!$
!!$          vec_m(1) = csr(1,1,imol,ls)
!!$          vec_m(2) = csr(2,1,imol,ls)
!!$          vec_m(3) = csr(3,1,imol,ls)
!!$
!!$          vec_m(:) = vec_m(:) - vec_o(:)
!!$
!!$          vec_h1(1) = csr(1,2,imol,ls)
!!$          vec_h1(2) = csr(2,2,imol,ls)
!!$          vec_h1(3) = csr(3,2,imol,ls)
!!$
!!$          vec_h1(:) = vec_h1(:) - vec_o(:)          
!!$
!!$          vec_h2(1) = csr(1,3,imol,ls)
!!$          vec_h2(2) = csr(2,3,imol,ls)
!!$          vec_h2(3) = csr(3,3,imol,ls)
!!$
!!$          vec_h2(:) = vec_h2(:) - vec_o(:) 
!!$
!!$          call util_images(vec_o,hmatrix(:,:,ls),recip_matrix(:,:,ls))
!!$
!!$          vec_m(:) = vec_o(:) + vec_m(:)
!!$          
!!$          vec_h1(:) = vec_h1(:) + vec_o(:)
!!$          vec_h2(:) = vec_h2(:) + vec_o(:)          
!!$
!!$          xxx(k) = vec_o(1)*bohr_to_ang
!!$          yyy(k) = vec_o(2)*bohr_to_ang
!!$          zzz(k) = vec_o(3)*bohr_to_ang
!!$
!!$          xxx(k+1) = vec_h1(1)*bohr_to_ang
!!$          yyy(k+1) = vec_h1(2)*bohr_to_ang
!!$          zzz(k+1) = vec_h1(3)*bohr_to_ang
!!$
!!$          xxx(k+2) = vec_h2(1)*bohr_to_ang
!!$          yyy(k+2) = vec_h2(2)*bohr_to_ang
!!$          zzz(k+2) = vec_h2(3)*bohr_to_ang
!!$
!!$          xxx(k+3) = vec_m(1)*bohr_to_ang
!!$          yyy(k+3) = vec_m(2)*bohr_to_ang
!!$          zzz(k+3) = vec_m(3)*bohr_to_ang
!!$
!!$          k = k + 4
!!$
!!$       case("tip4p")
!!$
!!$          vec_o(1) = ljr(1,1,imol,ls)
!!$          vec_o(2) = ljr(2,1,imol,ls)
!!$          vec_o(3) = ljr(3,1,imol,ls)
!!$
!!$          vec_m(1) = csr(1,1,imol,ls)
!!$          vec_m(2) = csr(2,1,imol,ls)
!!$          vec_m(3) = csr(3,1,imol,ls)
!!$
!!$          vec_m(:) = vec_m(:) - vec_o(:)
!!$
!!$          vec_h1(1) = csr(1,2,imol,ls)
!!$          vec_h1(2) = csr(2,2,imol,ls)
!!$          vec_h1(3) = csr(3,2,imol,ls)
!!$
!!$          vec_h1(:) = vec_h1(:) - vec_o(:)          
!!$
!!$          vec_h2(1) = csr(1,3,imol,ls)
!!$          vec_h2(2) = csr(2,3,imol,ls)
!!$          vec_h2(3) = csr(3,3,imol,ls)
!!$
!!$          vec_h2(:) = vec_h2(:) - vec_o(:) 
!!$
!!$          call util_images(vec_o,hmatrix(:,:,ls),recip_matrix(:,:,ls))
!!$
!!$          vec_m(:) = vec_o(:) + vec_m(:)
!!$          
!!$          vec_h1(:) = vec_h1(:) + vec_o(:)
!!$          vec_h2(:) = vec_h2(:) + vec_o(:)          
!!$
!!$          xxx(k) = vec_o(1)*bohr_to_ang
!!$          yyy(k) = vec_o(2)*bohr_to_ang
!!$          zzz(k) = vec_o(3)*bohr_to_ang
!!$
!!$          xxx(k+1) = vec_h1(1)*bohr_to_ang
!!$          yyy(k+1) = vec_h1(2)*bohr_to_ang
!!$          zzz(k+1) = vec_h1(3)*bohr_to_ang
!!$
!!$          xxx(k+2) = vec_h2(1)*bohr_to_ang
!!$          yyy(k+2) = vec_h2(2)*bohr_to_ang
!!$          zzz(k+2) = vec_h2(3)*bohr_to_ang
!!$
!!$          xxx(k+3) = vec_m(1)*bohr_to_ang
!!$          yyy(k+3) = vec_m(2)*bohr_to_ang
!!$          zzz(k+3) = vec_m(3)*bohr_to_ang
!!$
!!$          k = k + 4
!!$
!!$       case default
!!$
!!$          write(0,*)'HISTORY file output for selected model not implemented'
!!$
!!$       end select
!!$
!!$    end do
!!$
!!$    write(hist)nstep,natms,keytrj,imcon,tstep
!!$    write(hist)hmatrix(:,:,ls)*bohr_to_ang
!!$    write(hist)xxx
!!$    write(hist)yyy
!!$    write(hist)zzz
!!$
!!$
!!$    deallocate(xxx,yyy,zzz,stat=ierr)
!!$    if (ierr/=0) stop 'Error releasing memory in io_hist_append'
!!$
!!$    close(hist)
!!$
!!$  end subroutine io_hist_append

  subroutine io_write_psf()
    !------------------------------------------------------!
    ! Writes a VMD compatible psf file for a collection of !
    ! non-bonded beads. The symbol for oxygen is used for  !
    ! each bead, this is an arbitary choice. The radius    !
    ! of the beads can be changed once loaded into vmd.    !
    !------------------------------------------------------!
    ! Adapted from earlier code, D. Quigley April 2015     !
    !------------------------------------------------------!
    use userparams, only : nwater, num_lattices
    implicit none    
    integer :: i,ierr,k

    ! open the psf file
    open(unit=psf,file='mW.psf',status='replace',iostat=ierr)
    if (ierr/=0) stop 'Error opening mW.psf for output'

    ! write the header
    write(psf,'(A3)')'PSF'
    write(psf,'("         1 !NTITLE")')
    write(psf,*)
    write(psf,'(I8,1x,"!NATOM")')nwater*num_lattices

    ! format for x-plor style psf
10  format(I8,1x,A4,1x,I4,1x,A4,1x,A4,1x,A5,1x,F10.6,1x,5x,F8.4,10x,"0")

    ! Write the atoms information. The last entry is the mass
    ! which will be ignored by VMD anyway.
    k = 1
    do i = 1,nwater*num_lattices
       write(psf,10)k,"BULK",i,"UNK ","O ","O ",0.0,1.0
       k = k + 1
    end do
    
    ! Write the number of bonds
    write(psf,*)
    write(psf,'(I8,1x,"!NBOND: bonds")')0

    ! write irrelevant stuff - chimera will check for it and
    ! complain if it isn't there.
    write(psf,'(I8,1x,"!NTHETA: angles")')0
    write(psf,'(I8,1x,"!NPHI: torsions")')0
    write(psf,'(I8,1x,"!NIMPHI: torsions")')0    
    write(psf,'(I8,1x,"!NDON: donors")')0 
    write(psf,'(I8,1x,"!NACC: acceptors")')0 
    
    close(psf)
    
    return

  end subroutine io_write_psf

  subroutine io_write_dcd_header()
    !------------------------------------------------------!
    ! Writes the header of a VMD compatible dcd file for a !
    ! linear chain of bonded beads.                        !
    !------------------------------------------------------!
    ! Adapted from earlier code, D. Quigley April 2015     !
    !------------------------------------------------------!
    use userparams, only : nwater,num_lattices
    implicit none

    ! arrays for header
    integer,dimension(20) :: icntrl
    character(4) :: hdr='CORD'
    character*80,dimension(32) :: dcdtitle

    integer :: i,ierr

    open(unit=dcd,file='mW.dcd',status='replace',iostat=ierr,form='unformatted')
    if (ierr/=0) stop 'Error opening chain.dcd file - quitting'

    ! write the dcd header - most of this will be ignored
    icntrl(1)     = 1000                 ! number of snapshots in history file
    icntrl(2)     = 0
    icntrl(3)     = 100                  ! gap in steps between snapshots (doesn't matter)
    icntrl(4)     = 100*1000             ! total numbe of steps (VMD ignores this)
    icntrl(5:7)   = 0
    icntrl(8)     = nwater*num_lattices  ! Ndeg 
    icntrl(9)     = 0                    ! no fixed atoms
    icntrl(10)    = 0
    icntrl(11)    = 1                    ! 1/0 for unit cell presence
    icntrl(12:19) = 0
    icntrl(20)    = 24                   ! Charmm version number (fixes dcd format)

    write(dcd)hdr,icntrl
    write(dcd)1,(dcdtitle(i),i=1,1)
    write(dcd)nwater*num_lattices

    close(dcd)

    return

  end subroutine io_write_dcd_header

  subroutine io_write_dcd_snapshot()
    !======================================================!
    ! Writes a snapshot of the current positions to a the  !
    ! dcd file. Expects a 2D array r(1:3,1:nchain) in      !
    ! double precision which holds the coordinates.        !
    !======================================================!
    use constants,  only : bohr_to_ang  
    use userparams, only : nwater,num_lattices
    use model,      only : ljr,hmatrix,ls
    implicit none
    
    real(kind=dp),allocatable,dimension(:,:) :: rcopy
    real(kind=dp),dimension(3)   :: unita,unitb,unitc
    real(kind=dp),parameter :: invPi = 1.0_dp/3.141592653589793238462643383279502884197_dp

    ! charmm style cell vector array
    real(kind=dp),dimension(6) :: xtlabc   

    integer :: ierr,i,os

    allocate(rcopy(1:3,nwater*num_lattices),stat=ierr)

    ! Active lattice
    rcopy(:,1:nwater) = ljr(:,1,1:nwater,ls)*bohr_to_ang

    ! Other lattice if present
    if (num_lattices==2) then

       if (ls==1) os = 2
       if (ls==2) os = 1
       rcopy(:,nwater+1:nwater*num_lattices) = ljr(:,1,1:nwater,os)*bohr_to_ang

    end if

       
    open(unit=dcd,file='mW.dcd',status='old',position='append',iostat=ierr,form='unformatted')
    if (ierr/=0) stop 'Error opening mW.dcd file - quitting'
    
    xtlabc(1) = sqrt(dot_product(hmatrix(:,1,ls),hmatrix(:,1,ls)))*bohr_to_ang
    xtlabc(3) = sqrt(dot_product(hmatrix(:,2,ls),hmatrix(:,2,ls)))*bohr_to_ang
    xtlabc(6) = sqrt(dot_product(hmatrix(:,3,ls),hmatrix(:,3,ls)))*bohr_to_ang

    unita(:)  = hmatrix(:,1,ls)/xtlabc(1)
    unitb(:)  = hmatrix(:,2,ls)/xtlabc(3)
    unitc(:)  = hmatrix(:,3,ls)/xtlabc(6)

    xtlabc(2) = acos(dot_product(unita,unitb))*180.0*invPi
    xtlabc(4) = acos(dot_product(unita,unitc))*180.0*invPi
    xtlabc(5) = acos(dot_product(unitb,unitc))*180.0*invPi

    ! Write the information to file, note conversion to single
    ! precision to save on file size.
    write(dcd)xtlabc
    write(dcd)(-real(rcopy(1,i),kind=4),i=1,nwater*num_lattices)
    write(dcd)(-real(rcopy(2,i),kind=4),i=1,nwater*num_lattices)
    write(dcd)( real(rcopy(3,i),kind=4),i=1,nwater*num_lattices)

    !call flush(dcd)
    close(dcd)

    deallocate(rcopy,stat=ierr)
    
  end subroutine io_write_dcd_snapshot


  subroutine io_hist_header()
    !------------------------------------------------------------------------------!
    ! Writes the header of a DLPOLY HISTORY file. This is appended to at intervals !
    ! in the calculation.                                                          !
    !------------------------------------------------------------------------------!
    ! D.Quigley October 2006                                                       !
    !------------------------------------------------------------------------------!
    use userparams, only : model_type,nwater,num_lattices
    implicit none
    character(80) :: header
    character(3)  :: rankstr
    real(8)       :: natms
    character(8),allocatable,dimension(:)  :: atname
    real(8),allocatable,dimension(:)       :: weight
    real(8),allocatable,dimension(:)       :: ocharge

    integer :: ierr,imol,n,k
    integer :: myrank = 0

    write(hist_filename,'("HISTORY",I3.3)')myrank
    open(unit=hist,file=hist_filename,status='replace',form='unformatted',iostat=ierr)
    if (ierr/=0) stop 'Error opening HISTORY file'

    write(rankstr,'(I3)')myrank
    header="HISTORY file for rank "//rankstr

    select case (model_type)
    case ("mW")
       natms = real(num_lattices*nwater,kind=8)
    case default
       write(0,*)'HISTORY file output for selected model not implemented'
    end select

    n = anint(natms,kind=dp)

    allocate(atname(1:n),stat=ierr)
    allocate(weight(1:n),stat=ierr)
    allocate(ocharge(1:n),stat=ierr)
    if (ierr/=0) stop 'Error allocating memory in io_hist_header'

    k = 1
    do imol = 1,num_lattices*nwater

       select case (trim(model_type))
       case ("mW")

          atname(k)   = 'O_mW'
          ocharge(k)  = 0.0_dp
          weight(k)   = 15.9998_dp+2.0_dp * 1.0080_dp

          k = k + 1

       case default

          write(0,*)'HISTORY file output for selected model not implemented'

       end select

    end do


    ! write the HISTORY FILE HEADER
    write(hist)header
    write(hist)natms
    write(hist)atname
    write(hist)weight
    write(hist)ocharge

    deallocate(atname,weight,ocharge,stat=ierr)
    if (ierr/=0) stop 'Error releasing memory in io_hist_header'

    close(hist)


  end subroutine io_hist_header

  subroutine io_hist_append(icyc)
    !------------------------------------------------------------------------------!
    ! Appends the current configuration to a DLPOLY HISTORY file.                  !
    !------------------------------------------------------------------------------!
    ! D.Quigley October 2006                                                       !
    !------------------------------------------------------------------------------!
    use util,       only : util_images
    use constants,  only : bohr_to_ang
    use userparams, only : model_type,nwater,num_lattices
    use model,    only   : ljr,hmatrix,ls
    implicit none

    integer,intent(in) :: icyc

    real(8) :: nstep
    real(8) :: natms
    real(8) :: keytrj
    real(8) :: imcon
    real(8) :: tstep

    real(8),allocatable,dimension(:) :: xxx,yyy,zzz

    real(8),dimension(3) :: vec_o

    integer :: imol,ierr,n,k,os,ll

    integer :: myrank=0

    write(hist_filename,'("HISTORY",I3.3)')myrank
    open(unit=hist,file=hist_filename,status='old',form='unformatted',position='append',iostat=ierr)
    if (ierr/=0) stop 'Error opening HISTORY file'


    select case (trim(model_type))
    case ("mW")
       natms = real(num_lattices*nwater,kind=8)
    case default
       write(0,*)'HISTORY file output for selected model not implemented'
    end select


    nstep  = real(icyc,kind=8)
    keytrj = 0.0_dp
    imcon  = 3.0_dp
    tstep  = 1.0_dp

    n = anint(natms,kind=dp)
    allocate(xxx(1:n),stat=ierr)
    allocate(yyy(1:n),stat=ierr)
    allocate(zzz(1:n),stat=ierr)
    if (ierr/=0) stop 'Error allocating coords in io_hist_append'

    ll = ls
    if (num_lattices==1) ll = 1
    k = 1
    do imol = 1,nwater

       select case (trim(model_type))
       case ("mW")

          vec_o(1) = ljr(1,1,imol,ll)
          vec_o(2) = ljr(2,1,imol,ll)
          vec_o(3) = ljr(3,1,imol,ll)

          xxx(k) = vec_o(1)*bohr_to_ang +2.0_dp*hmatrix(1,1,ll)*bohr_to_ang
          yyy(k) = vec_o(2)*bohr_to_ang +2.0_dp*hmatrix(2,1,ll)*bohr_to_ang
          zzz(k) = vec_o(3)*bohr_to_ang +2.0_dp*hmatrix(3,1,ll)*bohr_to_ang

          k = k + 1

       case default

          write(0,*)'HISTORY file output for selected model not implemented'

       end select
      
    end do

    if (num_lattices==2) then

       if (ll==1) os = 2
       if (ll==2) os = 1

       do imol = 1,nwater

          select case (trim(model_type))
          case ("mW")
     
             vec_o(1) = ljr(1,1,imol,os)
             vec_o(2) = ljr(2,1,imol,os)
             vec_o(3) = ljr(3,1,imol,os)    
 
             xxx(k) = vec_o(1)*bohr_to_ang +2.0_dp*hmatrix(1,1,os)*bohr_to_ang
             yyy(k) = vec_o(2)*bohr_to_ang +2.0_dp*hmatrix(2,1,os)*bohr_to_ang
             zzz(k) = vec_o(3)*bohr_to_ang +2.0_dp*hmatrix(3,1,os)*bohr_to_ang

             k = k + 1

          case default
             
             write(0,*)'HISTORY file output for selected model not implemented'
             
          end select
          
       end do
             
    end if ! if num_lattices==2

    
    write(hist)nstep,natms,keytrj,imcon,tstep
    write(hist)hmatrix(:,:,ls)*bohr_to_ang
    write(hist)xxx
    write(hist)yyy
    write(hist)zzz


    deallocate(xxx,yyy,zzz,stat=ierr)
    if (ierr/=0) stop 'Error releasing memory in io_hist_append'

    close(hist)

  end subroutine io_hist_append

  subroutine io_write_global_header(restart)
    !------------------------------------------------------------------------------!
    ! Writes initialisation info / input and basic information about this job to   !
    ! the global output file.                                                      !
    !------------------------------------------------------------------------------!
    ! D.Quigley April 2011                                                         !
    !------------------------------------------------------------------------------!
    use userparams
    use comms, only : myrank,size
!$    use omp_lib   
    implicit none
    logical,intent(in) :: restart
    integer :: ierr

    if (myrank/=0) return ! only rank 0 writes global log

    ! Open the main log file
    if (restart) then
       open(unit=glog,file='mc.log',status='old',position='append',iostat=ierr)
       if (ierr/=0) stop 'Error opening mc.log for append'
    else
       open(unit=glog,file='mc.log',status='replace',iostat=ierr)
       if (ierr/=0) stop 'Error creating new mc.log'
    end if
    
    
    write(glog,'("#==============================================================#")')
    write(glog,'("#      Lattice-switching MC code for mW water molecules        #")')
    write(glog,'("#                                                              #")')
    write(glog,'("#                         D. Quigley                           #")')
    write(glog,'("#                 University of Warwick, UK                    #")')
    write(glog,'("#==============================================================#")')
    write(glog,'("#                                                              #")')
    write(glog,'("# Number of MPI tasks               : ",I6,"                   #")')size
!$  write(glog,'("# Number of OpenMP threads per task : ",I6,"                   #")')omp_get_num_threads()
    write(glog,'("#                                                              #")')   
!$    if ( num_lattices < omp_get_num_threads() ) then
!$    write(glog,'("# Warning - running more threads than lattices!                #")')
!$    write(glog,'("#                                                              #")')   
!$    end if
    if (restart) then

    end if

    call flush(glog)

    return

  end subroutine io_write_global_header

  subroutine io_write_rank_header(restart)
    !------------------------------------------------------------------------------!
    ! Writes initialisation info / input and basic information about this job to   !
    ! the global output file.                                                      !
    !------------------------------------------------------------------------------!
    ! D.Quigley April 2011                                                         !
    !------------------------------------------------------------------------------!
    use comms, only : myrank,size
!$  use omp_lib
    implicit none
    logical,intent(in) :: restart

    character(3) :: myrstr
    integer :: ierr

    write(myrstr,'(I3.3)')myrank

    if (restart) then
       open(unit=mylog,file='node'//myrstr//'.log',status='old',position='append',iostat=ierr)
       if (ierr/=0) then
          write(0,'("Error on rank ",I3," could not open node",I3.3,".log for append")')myrank,myrank
       end if
    else
       open(unit=mylog,file='node'//myrstr//'.log',status='replace',iostat=ierr)
       if (ierr/=0) then
          write(0,'("Error on rank ",I3," could not create node",I3.3,".log")')myrank,myrank
       end if
    end if

    write(mylog,'("#==============================================================#")')
    write(mylog,'("# Log file for rank ",I6,"    of ",I6,"                        #")')myrank+1,size
    write(mylog,'("#==============================================================#")')
!$  write(mylog,'("# Number of threads : ",I6,"                                   #")')omp_get_max_threads()
    write(mylog,'("#                                                              #")')   

    call flush(mylog)

    return

  end subroutine io_write_rank_header

end module io
