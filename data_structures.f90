! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: Git -*-
!=============================================================================!
!                        D A T A _ S T R U C T U R E S                        !
!=============================================================================!
! Defines the data structures used to hold the molecular positions and        !
! orientations as well as the interaction parameters for the model in use.    !
!=============================================================================!

module model

  use constants, only : dp,int32                ! Minimal useage where practical
!  use hbond_module                              ! G. Tribello's H-bond module

  implicit none                                 !Impose strong typing
  private   

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
                                                      !... unless exposed here.
  
  public :: create_model      ! create  water data structures
  public :: destroy_model     ! destroy water data structures

  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  public :: ljspm,cspm
  public :: mass
  public :: ljr,ref_ljr
  public :: hmatrix,ref_hmatrix,recip_matrix,ref_recipmatrix,volume,ls


  ! masses of each site
  integer,save :: nsites
  real(kind=dp),allocatable,dimension(:),save :: mass

  ! array of lennard-jones sites and reference state for quaternions
  real(kind=dp),allocatable,dimension(:,:,:,:),save :: ljr,ref_ljr

  ! matrix of cell vectors
  real(kind=dp),allocatable,dimension(:,:,:),save :: hmatrix,ref_hmatrix

  ! matrix of reciprocal cell vectors
  real(kind=dp),allocatable,dimension(:,:,:),save :: recip_matrix,ref_recipmatrix

  !system volume
  real(kind=dp),allocatable,dimension(:),save :: volume

  ! current lattice-switching reference system
  integer,save :: ls = 1

  integer :: ljspm,cspm

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!


  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

contains

  subroutine create_model
    !------------------------------------------------------------------------------!
    ! Allocates the array of water moleclues after assiging interaction parameters !
    ! and reference geometries. Note that if the initial configuration  is to be   !
    ! read from an xmol file the reference molecules are overwritten with that     !
    ! data.                                                                        !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------!    
    use constants, only : ang_to_bohr
    use userparams,only : model_type,nwater,num_lattices
    implicit none

    integer :: imol    ! loop counter
    integer :: ierr ! error flag


    ! Allocate matricies and volume
    allocate(hmatrix(1:3,1:3,1:num_lattices),stat=ierr)
    if (ierr/=0) stop 'Error allocating memory for hmatrices'

    allocate(ref_hmatrix(1:3,1:3,1:num_lattices),stat=ierr)
    if (ierr/=0) stop 'Error allocating memory for ref_hmatrices'

    allocate(recip_matrix(1:3,1:3,1:num_lattices),stat=ierr)
    if (ierr/=0) stop 'Error allocating memory for recip_matrices'

    allocate(ref_recipmatrix(1:3,1:3,1:num_lattices),stat=ierr)
    if (ierr/=0) stop 'Error allocating memory for ref_recip_matrices'


    allocate(volume(1:num_lattices),stat=ierr)
    if (ierr/=0) stop 'Error allocating memory for volumes'

    select case(trim(model_type))

    case ('mW')

       ljspm = 1 ! one Lennard-Jones site on the oxygen
       cspm =  0 ! no charge sites

       nsites = 1

       allocate(mass(1:cspm+ljspm),stat=ierr)
       if (ierr/=0) stop 'error allocating masses'

       mass(1) = 15.9998_dp+2.0_dp*1.0080_dp


       ! allocate reference water arrays
       allocate(ref_ljr(1:3,1:1,1:nwater,1:num_lattices),stat=ierr)

       if (ierr/=0) stop 'Error allocating reference arrays'

       do imol = 1,nwater

          ! oxygen Lennard-Jones site at the origin
          ref_ljr(:,1,imol,:) = 0.0_dp

       end do

       ! units
       ref_ljr = ref_ljr*ang_to_bohr

    case default

       write(*,*)'Unrecognised water model'
       stop

    end select


    ! allocate working arrays
    allocate(ljr(1:3,1:ljspm,1:nwater,1:num_lattices),stat=ierr)
    if (ierr/=0) stop 'Error allocating working position arrays'

    return

  end subroutine create_model

  subroutine destroy_model()
    !------------------------------------------------------------------------------!
    ! Release all memory used by the data structures in this module.               !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------!
    implicit none

    integer :: ierr

    ! deallocate matrices and volumes
    deallocate(hmatrix,ref_hmatrix,stat=ierr)
    deallocate(recip_matrix,stat=ierr)
    deallocate(ref_recipmatrix,stat=ierr)
    deallocate(volume,stat=ierr)
    if (ierr/=0) stop 'Error releasing matrix/volume memory'

    ! deallocate interaction parameters
    deallocate(mass,stat=ierr)
    if (ierr/=0) stop 'Error releasing model parameters'

    ! deallocate positions and quaternions
    deallocate(ref_ljr,stat=ierr)
    deallocate(ljr,stat=ierr)
    if (ierr/=0) stop 'Error releasing postions and quaternions'

    return

  end subroutine destroy_model

end module model
