! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: Git -*-
!=============================================================================!
!                             M C _ W A T E R                                 !
!=============================================================================!
program mc_water

  use constants,   only : dp,kb,aud_to_kgm3,water_mass,hart_to_ev,bohr_to_ang
  use random,      only : random_set_random_seed,random_uniform_random,random_test_uniform
  use io,          only : io_read_input,io_write_dcd_header,io_write_dcd_snapshot,mylog,glog, &
                          io_write_rank_header,io_write_global_header,mytherm,io_write_psf
  use init,        only : read_xmol
  use comms,       only : myrank,size,comms_initialise,comms_allocate,comms_finalise, &
                          comms_bcastlog,comms_allocate
  use userparams,  only : seedname,max_mc_cycles,eq_mc_cycles,temperature, &
                          file_output_int,nwater,leshift,input_ref_enthalpy, &
                           traj_output_int,pressure,mc_ensemble,wl_factor,num_lattices,samplerun
  use mc_moves,    only : mc_cycle,mc_init,mc_deinit,ls_mu,ref_enthalpy,mc_checkpoint_write
  use energy,      only : model_energy,compute_model_energy,energy_init,energy_deinit
  use model,       only : create_model,destroy_model,volume,ls,hmatrix
  use util,        only : util_hmatrix_to_abc
  use timer

  implicit none

  real(kind=dp) :: density
  integer       :: icyc,ierr,ils
  logical       :: restart,equil,safe,restart1,restart2
  integer       :: start_cycle = 0
  character(3)  :: myrstr

  real(kind=dp) :: alength,blength,clength,alpha,beta,gamma

  logical :: RNG_test = .true.

  external :: cleanexit

  !---------------------------------!
  ! Set comms going                 !
  !---------------------------------! 
  call comms_initialise()
  write(myrstr,'(I3.3)')myrank

  !---------------------------------!
  ! Set timer going                 !
  !---------------------------------! 
  call timer_init()

  !---------------------------------!
  ! Seed random number generator    !
  !---------------------------------!
  ! N.B. extra two optional arguments ensure
  ! different sequence of RNs on each rank
  call random_set_random_seed(0,dumrank=myrank,dumsize=size)

  !---------------------------------!
  ! Read input files                !
  !---------------------------------!
  call io_read_input()

  !---------------------------------!
  ! Is this a restart, if so count  !
  ! the number of cycles written    !
  ! so far.                         !
  !---------------------------------!
  inquire(file='checkpoint'//myrstr//'.dat.1',exist=restart1)
  inquire(file='checkpoint'//myrstr//'.dat.2',exist=restart2)
  restart = restart1.or.restart2
  call comms_bcastlog(restart,1)

  !---------------------------------!
  ! Open log files                  !
  !---------------------------------!
  call io_write_global_header(restart)
  call io_write_rank_header(restart)

  !---------------------------------!
  ! Test random number generator    !
  !---------------------------------!
  if (RNG_test) then
     call random_test_uniform(mylog)
  end if


  !---------------------------------!
  ! Open sample files for output    !
  !---------------------------------!
  if ( restart ) then
     open(unit=mytherm,file=trim(seedname)//myrstr//'_therm.dat',status='old',position='append',iostat=ierr)
     if (ierr/=0) stop 'Error opening therm.dat'
  else
     open(unit=mytherm,file=trim(seedname)//myrstr//'_therm.dat',status='replace',iostat=ierr)
     if (ierr/=0) stop 'Error opening therm.dat'
  end if

  !---------------------------------!
  ! Create initial system           !
  !---------------------------------!
  call create_model()

  !---------------------------------!
  ! psf file and dcd header         !
  !---------------------------------!
  !call io_hist_header
  call io_write_psf()
  call io_write_dcd_header()
  
  !---------------------------------!
  ! Initial structure               !
  !---------------------------------!
  call read_xmol()

  !---------------------------------!
  ! Compute initial energy          !
  !---------------------------------!
  call energy_init()

  !---------------------------------!
  ! Allocate buffers.               !
  !---------------------------------!
  call comms_allocate()

  !---------------------------------!
  ! Report initial energy           !
  !---------------------------------!
  do ils = 1,num_lattices
 !    write(0,'("Computing energy for lattice ",I3)')ils
     call compute_model_energy(ils)      ! Recompute energies  
  end do
  write(mylog,'("#                                                              #")') 
  write(mylog,'("# Reference lattices                                           #")')
  write(mylog,'("# -------------------------------------------------------------#")')
  if (num_lattices == 2) then
     write(mylog,'("# Computed energies   = ",2F15.6,"  eV     #")')model_energy*hart_to_ev
     write(mylog,'("# Computed enthalpies = ",2F15.6,"  eV     #")')(model_energy+pressure*volume)*hart_to_ev 
  else
     write(mylog,'("# Computed energy     = ",F15.6,"  eV                    #")')model_energy(1)*hart_to_ev
     write(mylog,'("# Computed enthalpy   = ",F15.6,"  eV                    #")')(model_energy(1)+pressure*volume(1))*hart_to_ev 
  end if
  write(mylog,'("#                                                              #")')   
  call flush(mylog)

  !---------------------------------!
  ! Shift energies of two lattices  !
  ! to match?                       !
  !---------------------------------!
  ref_enthalpy(1:num_lattices) = model_energy(1:num_lattices) 
  if (mc_ensemble=='npt') ref_enthalpy(1:num_lattices) = ref_enthalpy(1:num_lattices) + pressure*volume(1:num_lattices)

  ! Override this with input_ref_enthalpy if non-zero
  if (any(abs(input_ref_enthalpy)>tiny(1.0_dp))) ref_enthalpy = input_ref_enthalpy

  !---------------------------------!
  ! Set up mc module                !
  !---------------------------------!
  call mc_init(start_cycle)
  do ils = 1,num_lattices
     call compute_model_energy(ils)      ! Recompute energies  
  end do

  !------------------------------------!
  ! Setup a signal handler for SIGTERM !   
  ! which is usually signal no. 15.    !
  !------------------------------------!
  call signal('SIGTERM',cleanexit)
  ! KEB: I think this might just make it more generally portable but it now runs/compiles locally for me
  ! See hackpad for comments

  !--------------------------------------!
  ! Computing initial overlap parameter  !
  !--------------------------------------!
  if (num_lattices==2)  then
     beta  = 1.0_dp/(kB*temperature)
     ls_mu = model_energy(1) + pressure*volume(1) - model_energy(2) - pressure*volume(2) 
     if (leshift) ls_mu = ls_mu - ref_enthalpy(1) + ref_enthalpy(2)
     ls_mu = ls_mu*beta - real(Nwater,kind=dp)*log(volume(1)/volume(2))
  end if

  !-------------------------------------!
  ! Loop over Monte-Carlo cycles        !
  !-------------------------------------!
  do icyc = start_cycle+1,start_cycle+max_mc_cycles

     !-------------------------------------!
     ! Are we equilibrated?                !
     !-------------------------------------!
     if (mod(icyc,eq_mc_cycles)==0) equil = .false.

     !-------------------------------------!
     ! Perform a single MC cycle           !
     !-------------------------------------!
     call mc_cycle()

     !-------------------------------------!
     ! Write configuration to xmol file.   !
     !-------------------------------------!
     if ( (mod(icyc,traj_output_int)==0) ) call io_write_dcd_snapshot()

     !-------------------------------------!
     ! Write energy and density to unit 26 !
     !-------------------------------------!
     if ( (mod(icyc,file_output_int)==0) ) then

        density = real(nwater,kind=dp)*water_mass/volume(ls)

        if (num_lattices==1) then

           call util_hmatrix_to_abc(hmatrix(:,:,1),alength,blength,clength,alpha,beta,gamma)

           write(mytherm,'(I8,E15.6,5x,F15.6,6F15.6)')icyc,model_energy(1)*hart_to_ev,volume(1)*bohr_to_ang**3, &
                alength*bohr_to_ang,blength*bohr_to_ang,clength*bohr_to_ang,alpha,beta,gamma

        else

           if ( (wl_factor < tiny(1.0_dp)).or.(samplerun) ) then
              ! This is a sample run, need volume for histogram reweighting
              write(mytherm,'(I8,E15.6,5x,3F15.6,1x,I1)')icyc,model_energy(ls)*hart_to_ev,ls_mu,volume*bohr_to_ang**3,ls
           else
              ! Not a sample run
              write(mytherm,'(I8,E15.6,5x,2F15.6,1x,I1)')icyc,model_energy(ls)*hart_to_ev,ls_mu,density*aud_to_kgm3,ls
           end if

        end if

     end if

     !------------------------------------!
     ! Check that enough time remains in  !
     ! the queue to continue.             !
     !------------------------------------!
     call timer_check_runtime(safe)
     call comms_bcastlog(safe,1)

     if ( (.not.safe).and.(myrank==0) ) then
        write(glog,*)
        write(glog,'("!============================================!")')
        write(glog,'("! Approaching end of queue time - stopping   !")')
        write(glog,'("!============================================!")')
        write(glog,*)
     end if
     if (.not.safe) exit

  end do

  !-------------------------!
  ! Write final checkpoint  !
  !-------------------------!
  call mc_checkpoint_write()

  !-------------------------!
  ! Relase memory           !
  !-------------------------!
  call energy_deinit()
  call mc_deinit()
  call destroy_model()
  call comms_finalise()

  ! Close units?
  close(mylog)
  close(mytherm)
  if (myrank==0) close(glog)

end program mc_water

subroutine cleanexit()
  !-------------------------------------------------------------------------!
  ! Routine to handle clean shutdown of code in the event of SIGTERM being  !
  ! sent to the process by PBS or Condor etc.                               !
  !-------------------------------------------------------------------------!
  ! D.Quigley August 2011                                                   !
  !-------------------------------------------------------------------------!
  use mc_moves, only : mc_checkpoint_write,mc_deinit
  use energy,   only : energy_deinit
  use model,    only : destroy_model
  use comms,    only : comms_finalise,myrank
  use io,       only : mylog,mytherm,glog,mytherm
  implicit none
  
  !-------------------------!
  ! Write final checkpoint  !
  !-------------------------!
  call mc_checkpoint_write()

  !-------------------------!
  ! Relase memory           !
  !-------------------------!
  call energy_deinit()
  call mc_deinit()
  call destroy_model()
  call comms_finalise()

  ! Close units?
  close(mylog)
  close(mytherm)
  if (myrank==0) close(glog)

end subroutine cleanexit
