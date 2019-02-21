! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: Git -*-
!=============================================================================!
!                               E N E R G Y                                   !
!=============================================================================!
! Contains all subroutines associated with caculating the energy of the       !
! simulation cell. Where possible, terms which depend on the cell-vectors     !
! are only recalculated following a volume change. Routines are provided for  !
! computing the change in energy following a single molecule move.            !
!=============================================================================!
module energy

  Use constants, Only : dp,int32,ang_to_bohr    ! Minimal useage where practical

  Implicit None                                 ! Impose strong typing
  Private   

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
                                                      !... unless exposed here.

  public :: energy_init
  public :: energy_deinit
  public :: compute_ivects
  public :: compute_model_energy
  public :: compute_local_real_energy
  public :: compute_neighbours

  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  public :: model_energy
  public :: nivect,ivect

  ! current energy due to the model Hamiltonian
  real(kind=dp),allocatable,dimension(:),save :: model_energy

  ! lattice translation vectors
  integer,allocatable,dimension(:) :: nivect
  real(kind=dp),allocatable,dimension(:,:,:) :: ivect

  public :: mw_sigma,mw_epsilon,mw_lambda
  public :: sw_bigA,sw_B,sw_gamma,sw_a,sw_p,sw_q,cos0
  
  ! Constants defining the mW water model
  real(kind=dp) :: mw_sigma   = 2.3925_dp*ang_to_bohr     ! In Bohr
  real(kind=dp) :: mw_epsilon = 6.189_dp/627.509469_dp    ! In Hartree
  real(kind=dp) :: mw_lambda  = 23.15_dp                  ! Dimensionless

  ! General Stillinger Weber constants
  real(kind=dp) :: sw_bigA = 7.049556277_dp
  real(kind=dp) :: sw_B = 0.6022245584_dp
  real(kind=dp) :: sw_gamma = 1.2_dp
  real(kind=dp) :: sw_a = 1.8_dp
  integer :: sw_p=4,sw_q=0
  real(kind=dp) :: cos0 = -0.33331324756

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  ! Neighbour list
  integer,parameter :: maxneigh = 50
  integer,allocatable,dimension(:,:),save :: nn
  integer,allocatable,dimension(:,:,:),save ::jn,vn
  
contains

  subroutine energy_init
    !------------------------------------------------------------------------------!
    ! Initialises all arrays/caches used in calculcating water-water interaction   !
    ! energies.                                                                    !
    ! N.B. must be called after initial positions have been created. The arrays    !
    ! allocated here should be de-allocated in energy_deinit.                      !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------!
    use util,       only : util_determinant,util_images,util_recipmatrix
    use userparams, only : num_lattices,nwater
    use model,      only : hmatrix,volume
    implicit none

    integer :: ils,im,jm,km    ! loop counters
    integer :: ierr ! error flag

    allocate(model_energy(1:num_lattices),stat=ierr)
    if (ierr/=0) stop 'Error allocating model and recip energy arrays'

    allocate(nivect(1:num_lattices),stat=ierr)
    if (ierr/=0) stop 'Error allocating nivect'


    do ils = 1,num_lattices

       im = floor((sw_a*mw_sigma+1.0_dp)/sqrt(dot_product(hmatrix(:,1,ils),hmatrix(:,1,ils))))+1
       jm = floor((sw_a*mw_sigma+1.0_dp)/sqrt(dot_product(hmatrix(:,2,ils),hmatrix(:,2,ils))))+1
       km = floor((sw_a*mw_sigma+1.0_dp)/sqrt(dot_product(hmatrix(:,3,ils),hmatrix(:,3,ils))))+1

       nivect(ils) = (2*im+1)*(2*jm+1)*(2*km+1)
#ifdef DEBUG
       write(0,'("DEBUG - number of ivects for lattice ",I2," : ",I4)')ils,nivect(ils)
#endif
       volume(ils) = abs(util_determinant(hmatrix(:,:,ils)))

    end do

    allocate(ivect(1:3,1:maxval(nivect,1),1:num_lattices),stat=ierr)
    if (ierr/=0) stop 'Error allocating ivect'

    ! compute current image translation vectors
    do ils = 1,num_lattices
       call compute_ivects(ils)
    end do

    ! Neighbours
    allocate(nn(1:nwater,1:num_lattices),stat=ierr)
    if (ierr/=0) stop 'Error allocating nn array in molint.F90'
    allocate(vn(1:maxneigh,1:nwater,1:num_lattices),stat=ierr)
    if (ierr/=0) stop 'Error allocating vn array in molint.F90'
    allocate(jn(1:maxneigh,1:nwater,1:num_lattices),stat=ierr)
    if (ierr/=0) stop 'Error allocating jn array in molint.F90'


    do ils = 1,num_lattices
       call compute_neighbours(ils)
       call compute_model_energy(ils)
    end do

    return

  end subroutine energy_init

  subroutine energy_deinit()
    !------------------------------------------------------------------------------!
    ! Releases all memory allocated in energy_init                                 !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------!
    implicit none

    integer :: ierr

    ! deallocate image vectors
    deallocate(ivect,stat=ierr)
    if (ierr/=0) stop 'Error deallocating ivect'

    return

  end subroutine energy_deinit


  subroutine compute_ivects(ils)
    !------------------------------------------------------------------------------!
    ! Computes the translation vectors needed to include all images such that each !
    ! atom sees all the images (including those of itself) within the cut-off      !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2009                                                     !
    !------------------------------------------------------------------------------!
    use model,      only : hmatrix
    implicit none
    integer,intent(in) :: ils
 
    !loop counters
    integer :: icell,jcell,kcell,im,jm,km,k
    real(kind=dp),dimension(3) :: sx,sy,sz

    im = floor(sw_a*mw_sigma/sqrt(dot_product(hmatrix(:,1,ils),hmatrix(:,1,ils))))+1
    jm = floor(sw_a*mw_sigma/sqrt(dot_product(hmatrix(:,2,ils),hmatrix(:,2,ils))))+1
    km = floor(sw_a*mw_sigma/sqrt(dot_product(hmatrix(:,3,ils),hmatrix(:,3,ils))))+1
    
    nivect(ils) = (2*im+1)*(2*jm+1)*(2*km+1)

    ! we'd like the central cell to be entry 0
    ! we can flag it as non-self interacting
    ivect(:,1,ils) = 0.0_dp
    
    k = 2
    do icell = -im,im
       sx = real(icell,kind=dp)*hmatrix(:,1,ils)
       do jcell = -jm,jm
          sy = real(jcell,kind=dp)*hmatrix(:,2,ils) 
          do kcell = -km,km
             sz = real(kcell,kind=dp)*hmatrix(:,3,ils)
             
             if ( abs(icell)+abs(jcell)+abs(kcell) == 0 ) cycle
             ivect(:,k,ils)  = sx + sy + sz
             k = k + 1

          end do
       end do
    end do

    return

  end subroutine compute_ivects


  real(kind=dp) function compute_local_real_energy(imol,ils)
    !------------------------------------------------------------------------------!
    ! Calculates the real-space contribution to the energy due water number imol.  !
    ! To be used when computing the changes in real-space energy due to a trial    !
    ! single molecule move.                                                        !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------!
    use model,      only : ljr
    implicit none
    integer,intent(in) :: imol,ils

    ! local variables
    real(kind=dp)                :: Evdw,rcsq,ctheta,exp1,exp2,csq
    real(kind=dp)                :: r1_ij,r2_ij,r1_ik,r2_ik,tmpE,r1_jk,r2_jk
    real(kind=dp),dimension(3)   :: tmpvect,tmpvect2
    real(kind=dp),dimension(3)   :: ilj,jlj,klj

    integer :: jmol,ln,ln2 ! loop counters
    integer :: ji,ki,kmol

    Evdw  = 0.0_dp

    rcsq = mw_sigma*sw_a*mw_sigma*sw_a

    ilj(:) = ljr(:,1,imol,ils)

    do ln = 1,nn(imol,ils) ! loop over other molecule jmol
       
       jmol = jn(ln,imol,ils)   ! molecule
       ji   = vn(ln,imol,ils)   ! image
       
       jlj(:) = ljr(:,1,jmol,ils) + ivect(:,ji,ils)

       ! compute separation vector
       tmpvect = jlj(:) - ilj(:)
       r2_ij   = tmpvect(1)**2 + tmpvect(2)**2 + tmpvect(3)**2
             
       ! compute interactions if in range with the current image
       if ( r2_ij < rcsq ) then
          
          r1_ij = sqrt(r2_ij)
          
          ! Pair interaction
          exp2 = exp(mw_sigma/(r1_ij-mw_sigma*sw_a))
          tmpE = sw_bigA*mw_epsilon*(sw_B*(mw_sigma*mw_sigma/r2_ij)**2-1.0_dp)
          tmpE = tmpE*exp2

          exp2 = exp(sw_gamma*mw_sigma/(r1_ij-mw_sigma*sw_a))
          
          Evdw  = Evdw + tmpE
             
          ! Three body interactions like jmol--imol--kmol
          do ln2 = ln+1,nn(imol,ils) ! other neighbours of imol as 3rd body
             
             kmol = jn(ln2,imol,ils)
             ki   = vn(ln2,imol,ils)
             
             klj(:) = ljr(:,1,kmol,ils)  + ivect(:,ki,ils)    
             
             tmpvect2(:) = klj(:) - ilj(:)
             r2_ik      = tmpvect2(1)**2 + tmpvect2(2)**2 + tmpvect2(3)**2
             
             if ( r2_ik < rcsq ) then
                r1_ik = sqrt(r2_ik)
                
                ctheta = dot_product(tmpvect,tmpvect2)/(r1_ik*r1_ij)
                csq  = (ctheta - cos0)**2
                exp1 = exp(sw_gamma*mw_sigma/(r1_ik-mw_sigma*sw_a))
                Evdw = Evdw + mw_lambda*mw_epsilon*exp1*exp2*csq
                
             end if ! if kmol within range inside image ki
             
          end do ! 3rd body

          tmpvect = -tmpvect ! now vector from j to i

          ! Three body interactions like imol--jmol--kmol
          do ln2 = 1,nn(jmol,ils) ! other neighbours of jmol as 3rd body
             
             kmol = jn(ln2,jmol,ils)
             ki   = vn(ln2,jmol,ils)

             ! This is the information of neighbour ln2 of jmol in the central
             ! image, which needs translating by ivect(:,ji,ils) to get 
             ! neighbours of jmol in the relevant image.             
             klj(:) = ljr(:,1,kmol,ils)  + ivect(:,ki,ils) + ivect(:,ji,ils)    
             
             tmpvect2(:) = klj(:) - jlj(:)
             r2_jk      = tmpvect2(1)**2 + tmpvect2(2)**2 + tmpvect2(3)**2
             
             if ( r2_jk < rcsq ) then
                r1_jk= sqrt(r2_jk)
                
                ctheta = dot_product(tmpvect,tmpvect2)/(r1_jk*r1_ij)
                if (ctheta > 0.99_dp) then
                   !if (kmol/=imol) stop 'Indexing error'
                   cycle ! skip reverse vector back to imol
                end if
                csq  = (ctheta - cos0)**2
                exp1 = exp(sw_gamma*mw_sigma/(r1_jk-mw_sigma*sw_a))
                Evdw = Evdw + mw_lambda*mw_epsilon*exp1*exp2*csq
                
             end if ! if kmol within range inside image ki
             
          end do ! 3rd body

          
       end if ! jmol within range of jmol inside image ji
       
    end do  ! end loop over neighbours of imol
    
    


    ! set return value
    compute_local_real_energy = Evdw
    
    return 

  end function compute_local_real_energy


  subroutine compute_model_energy(ils)
    !------------------------------------------------------------------------------!
    ! Computes the energy of the entire simulation cell, using Lennard-Jones plus  !
    ! Ewald summation. To be used when sampling the energy of the entire system or !
    ! when computing the energy change from a volume/cell trial move.              !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------!
    use userparams, only : nwater
    use model,      only : ljr

    implicit none
    integer,intent(in) :: ils
    real(kind=dp),dimension(3)   :: tmpvect,tmpvect2
    real(kind=dp)                :: r1_ij,r2_ij,r1_ik,r2_ik
    real(kind=dp)                :: Evdw,csq,ctheta,exp1,exp2
    real(kind=dp)                :: tmpE,rcsq

    real(kind=dp),dimension(3) :: ilj,jlj,klj

    integer :: imol,jmol ! loop counters
    integer :: ji,ki,kmol,ln,ln2

    Evdw  = 0.0_dp

    rcsq = mw_sigma*sw_a*mw_sigma*sw_a


    !------------------------------------------------!
    !         mW water model (Stillinger Weber)      !
    !------------------------------------------------!
    do imol = 1,nwater  ! loop over central molecule imol

       ilj(:) = ljr(:,1,imol,ils)

       do ln = 1,nn(imol,ils) ! loop over other molecule jmol

          jmol = jn(ln,imol,ils)   ! molecule
          ji   = vn(ln,imol,ils)   ! image

          jlj(:) = ljr(:,1,jmol,ils) + ivect(:,ji,ils)

          ! compute separation vector
          tmpvect = jlj(:) - ilj(:)
          r2_ij   = tmpvect(1)**2 + tmpvect(2)**2 + tmpvect(3)**2
             
          ! compute interactions if in range with the current image
          if ( r2_ij < rcsq ) then
                
             r1_ij = sqrt(r2_ij)
                
             ! Pair interaction
             exp2 = exp(mw_sigma/(r1_ij-mw_sigma*sw_a))
             tmpE = sw_bigA*mw_epsilon*(sw_B*(mw_sigma*mw_sigma/r2_ij)**2-1.0_dp)
             tmpE = tmpE*exp2
             exp2 = exp(sw_gamma*mw_sigma/(r1_ij-mw_sigma*sw_a))
             
             Evdw  = Evdw + 0.5_dp*tmpE
                
             ! Three body interactions
             do ln2 = ln+1,nn(imol,ils) ! third body kmol

                kmol = jn(ln2,imol,ils)
                ki   = vn(ln2,imol,ils)
                   
                klj(:) = ljr(:,1,kmol,ils)  + ivect(:,ki,ils)    
                      
                tmpvect2(:) = klj(:) - ilj(:)
                r2_ik      = tmpvect2(1)**2 + tmpvect2(2)**2 + tmpvect2(3)**2
                      
                if ( r2_ik < rcsq ) then
                   r1_ik = sqrt(r2_ik)
                         
                   ctheta = dot_product(tmpvect,tmpvect2)/(r1_ik*r1_ij)
                   csq  = (ctheta - cos0)**2
                   exp1 = exp(sw_gamma*mw_sigma/(r1_ik-mw_sigma*sw_a))
                   Evdw = Evdw + mw_lambda*mw_epsilon*exp1*exp2*csq
                         
                end if ! if kmol within range inside image ki
                
             end do ! 3rd body
                   
          end if ! jmol within range of jmol inside image ji
                       
       end do  ! end loop over neighbours of imol

    end do ! end loop over imol
             
    model_energy(ils) = Evdw 

    return 

  end subroutine compute_model_energy

  subroutine compute_neighbours(ils)


    use userparams, only : nwater
    use model, only      : ljr
    implicit none
    integer,intent(in) :: ils

    integer :: imol,jmol,k,ni
    real(kind=dp) :: r2_ij,rn
    real(kind=dp),dimension(3) :: v_ij,ilj,jlj,tmpvect
    
    rn = sw_a*mw_sigma*1.25_dp

    call compute_ivects(ils)

    do imol = 1,nwater

       ilj(:) = ljr(:,1,imol,ils)

       nn(imol,ils) = 0
       do jmol = 1,nwater

          jlj(:) = ljr(:,1,jmol,ils)

          v_ij(:)   = jlj(:) - ilj(:)

          do k = 1,nivect(ils)
             if ( (k==1).and.(jmol==imol) )cycle

             tmpvect(:) = v_ij(:) + ivect(:,k,ils) ! apply image             
             r2_ij = dot_product(tmpvect,tmpvect)

             if (r2_ij<rn*rn) then
                nn(imol,ils) = nn(imol,ils) + 1
                ni = nn(imol,ils)        ! number of neighbours of imol
                jn(ni,imol,ils) = jmol   ! jmol is a neighbour of imol
                vn(ni,imol,ils) = k      ! in image k
             end if

          end do
       end do
       
       !write(0,'("Molecule ",I5," has ",I5," neighbours")')imol,nn(imol,ils)
       
    end do


  end subroutine compute_neighbours


end module energy
