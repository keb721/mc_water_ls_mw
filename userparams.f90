! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: Git -*-
!=============================================================================!
!                             U S E R P A R A M S                             !
!=============================================================================!
! Holds parameters specified by the user which define the calculation.        !
! Sensible defaults are specified in most cases.                              !
!=============================================================================!
module userparams

  use constants, only : dp,aup_to_atm,ang_to_bohr,hart_to_dlpol
  implicit none

  ! Model parameters
  character(12),save :: model_type  = "mW"               ! water model to use

  ! Configuration parameters
  integer,save        :: nwater = 768                    ! number of water molecules 
  integer,save        :: num_lattices = 2                ! should be two for lattice-switching
  character(6),save   :: method = 'xmol'                 ! initialisation method
  real(kind=dp),save  :: r_overlap = 1.7_dp*ang_to_bohr  ! min overlap if random

  ! Thermal parameters
  real(kind=dp),save :: pressure    = 1.0_dp/aup_to_atm  ! pressure in a.u.p.
  real(kind=dp),save :: temperature = 240.0_dp           ! temperature in Kelvin

  ! Monte-Carlo parameters
  character(3),save  :: mc_ensemble      = 'npt'         ! ensemble to sample
  real(kind=dp),save :: mc_max_trans     = 0.6_dp        ! maximum translation
  real(kind=dp),save :: mc_target_ratio  = 0.50_dp       ! target move acceptance ratio
  real(kind=dp),save :: mc_dv_max        = 0.1_dp        ! maximum cell vector change

  real(kind=dp),save :: wl_factor    = 0.05              ! Wang-Landau modification factor
  integer,save       :: wl_schedule  = 0                 ! 0 = flatness, 1 = min visits
  integer,save       :: wl_minhist   = 20                ! minimum number of hits per bin
  real(kind=dp),save :: wl_flattol   = 0.05              ! flatness tolerance 
  logical,save       :: wl_useinvt   = .false.           ! use 1/t method
  logical,save       :: wl_swetnam   = .false.           ! use Adam Swetnam's formula for f
  real(kind=dp),save :: wl_alpha     = 1.0               ! prefactor in above formula
  logical,save       :: eta_interp   = .true.            ! inerpolate between grid values of eta
  logical,save       :: samplerun    = .false.           ! sampling run with constant eta
  logical,save       :: leshift      = .false.           ! shift by reference energy config and correct in output

  integer,save       :: nbins = 200                      ! number of bins to use
  real(kind=dp),save :: mu_min = -8000.0_dp              ! minimum overlap parameter
  real(kind=dp),save :: mu_max = +8000.0_dp              ! maximum overlap parameter

  logical,save       :: allow_switch   = .true.          ! allow lattice switches
  logical,save       :: allow_vol      = .true.          ! allow volume moves
  logical,save       :: allow_trans    = .true.          ! allow translation moves

  real(kind=dp),save :: mc_trans_prob    = 0.5_dp        ! Relative mol translation P
  real(kind=dp),save :: mc_vol_prob      = 0.01_dp       ! Relative volume move P
  real(kind=dp),save :: mc_switch_prob   = 0.00_dp       ! Relative switch P
  logical,save       :: mc_always_switch = .true.        ! Attempt switch after every move

  ! Override reference energies
  real(kind=dp),dimension(2),save :: input_ref_enthalpy = (/0.0_dp,0.0_dp/)

  ! book-keeping parameters
  integer,save        :: list_update_int = 50            ! Verlet list update interval
  integer,save        :: traj_output_int = 5000000       ! HISTORY file dump interval
  integer,save        :: file_output_int = 5             ! sample interval
  integer,save        :: latt_sync_int   = 10000         ! internal at which to enforce synchronisation
  integer,save        :: mpi_sync_int    = 250           ! interval at which to reduce bins across all MPI tasks
  integer,save        :: chkpt_dump_int  = 1000          ! interval at which to write checkpoint file
  integer,save        :: monitor_int     = 1000          ! interval at which to report to log file
  integer,save        :: flat_chk_int    = 10000         ! interval at which to check histogram flatness
  integer,save        :: invt_dump_int   = 500000        ! interval at which to dump histograms in invt mode
  logical,save        :: eq_adjust_mc    = .false.       ! adjust mc parameters in equil
  integer,save        :: deltaG_int      = 100000        ! interval at which to compute deltaG
  integer,save        :: max_mc_cycles   = 1000          ! total number of mc cycles
  integer,save        :: eq_mc_cycles    = 25000         ! number of equilibration cycles

  ! parallelisation
  character(2),save   :: parallel_strategy = 'mw'        ! Multiple walkers of domain decomposition (dd)
  integer,save        :: window_overlap = 2              ! Number of bins to overlap at window boundaries

  ! seedname for calculation
  character(30),save :: seedname


end module userparams


