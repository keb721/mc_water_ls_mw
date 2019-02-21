! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: Git -*-
!=============================================================================!
!                                 I N I T                                     !
!=============================================================================!
! Allocates initial positions either randomly or by reading an Xmol file.     !
! note that no checking is made to determine if the co-ordinates read from    !
! file are consistent with the water model specified.                         !
!=============================================================================!

module init


  use constants,  only : dp                     !Minimal useage where practical

  implicit none                                 !Impose strong typing
  private                                       !Everything is private ...

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  public :: read_xmol
  !public :: continuation

  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

contains

  subroutine read_xmol()
    !------------------------------------------------------------------------------!
    ! Reads initial co-ordinates into the water data structures. Note that only    !
    ! postions of the oxygen and hydrogen atoms are read. Any massless charge sites!
    ! are assigned in this routine. No checking of O-H distances or HOH angles is  !
    ! implemented.                                                                 !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------! 
    use constants,  only : ang_to_bohr
    use comms,      only : myrank,comms_bcastint,comms_bcastreal
    use userparams, only : nwater,model_type,num_lattices
    use model,      only : ljr,ref_ljr,hmatrix, &
                           recip_matrix,ref_hmatrix,ref_recipmatrix
    use util,       only : util_recipmatrix
    implicit none
    
    integer,parameter :: xml = 40 ! Unit number for xmol files

    character(1) :: elemid        ! chemical symbol for the current atom
    character(3) :: lattstr

    integer :: ierr,imol,tmpN,ils

    real(kind=dp),dimension(3) :: o_vec

    do ils = 1,num_lattices

       ! open the correct xmol file for the current lattice number
       write(lattstr,'(I3.3)')ils
       
       if (myrank==0) open(unit=xml,file='input'//lattstr//'.xmol',status='old',iostat=ierr)
       call comms_bcastint(ierr,1)
       if (ierr/=0) stop 'Error opening input coords for reading'

       ! check the number of atoms in the xmol file is consistent with 
       ! that specified in the input file
       if (myrank==0) read(xml,*)tmpN
       call comms_bcastint(tmpN,1)
       if (tmpN/=nwater) stop 'Error wrong number of atoms in input.xmol'

       ! read the matrix of cell vectors
       if (myrank==0) read(xml,*)hmatrix(:,:,ils)
       call comms_bcastreal(hmatrix(1,1,ils),9)

       !write(*,'(3F15.6)')hmatrix(:,1,ils)
       !write(*,'(3F15.6)')hmatrix(:,2,ils)
       !write(*,'(3F15.6)')hmatrix(:,3,ils)

       hmatrix(:,:,ils) = hmatrix(:,:,ils)*ang_to_bohr
       call util_recipmatrix(hmatrix(:,:,ils),recip_matrix(:,:,ils))

       ref_hmatrix(:,:,ils)     = hmatrix(:,:,ils)
       ref_recipmatrix(:,:,ils) = recip_matrix(:,:,ils)

       do imol = 1,Nwater

          ! assume arrangement of OHH
          if (myrank==0) read(xml,*)elemid,o_vec
       
          call comms_bcastreal(o_vec(1),3)

          select case(trim(model_type))
          case ('mW')

             ljr(:,1,imol,ils) = o_vec(:)*ang_to_bohr

             ! .. relative to this configuration
             ref_ljr(:,:,imol,ils) = ljr(:,:,imol,ils)

          case default
             
             write(*,*)'input from file not yet supported for selected'
             write(*,*)'water model'
             stop
             
          end select

       end do

    end do


    if (myrank==0) close(xml)

 

  end subroutine read_xmol


!!$  subroutine renumber()
!!$
!!$    !------------------------------------------------------------------------------!
!!$    ! Attempts to renumber the molecules in lattice 2 such that the orientation    !
!!$    ! of a particular molecule matches the corresponding molecule in lattice 1 as  !
!!$    ! closely as possible. i.e. we minimise the difference between two vectors by  !
!!$    ! re-ordering the second one only.                                             !
!!$    !------------------------------------------------------------------------------!
!!$    ! D.Quigley October 2009                                                       !
!!$    !------------------------------------------------------------------------------! 
!!$    use constants, only  : pi,bohr_to_ang,Pi
!!$    use model, only      : ref_ljr,ref_csr,ljr,csr,ljspm,cspm,hmatrix
!!$    use userparams, only : nwater
!!$    use random, only     : random_uniform_random
!!$    use io,      only    : io_hist_append
!!$    implicit none
!!$    integer :: iwater,jwater,iloop,ilj,ics,ilev=1
!!$
!!$    real(kind=dp) :: sumdiff,sumdiff_old,gsumdiff_old,mix
!!$    real(kind=dp),allocatable,dimension(:,:) :: tmp_csr,tmp_ljr
!!$    integer,dimension(2) :: ierr
!!$
!!$
!!$    real(kind=dp),dimension(3,3) :: a
!!$    real(kind=dp) :: tq0,tq1,tq2,tq3,u1,u2,u3,norm,sq1,sq2,sq3
!!$
!!$
!!$    !return
!!$
!!$    !call io_hist_append(0)
!!$
!!$    ierr = 0
!!$    allocate(tmp_ljr(1:3,1:ljspm),stat=ierr(1))
!!$    allocate(tmp_csr(1:3,1:cspm),stat=ierr(2))
!!$    if (any(ierr/=0)) stop 'Error allocating memory in renumber'
!!$
!!$    sumdiff_old  = difference()
!!$    gsumdiff_old = sumdiff_old
!!$
!!$    !do
!!$
!!$    ! generate random quaternion
!!$    !u1 = random_uniform_random()
!!$    !u2 = random_uniform_random()
!!$    !u3 = random_uniform_random()
!!$
!!$    !tq0 = sqrt(1.0_dp-u1)*sin(2.0_dp*pi*u2)
!!$    !tq1 = sqrt(1.0_dp-u1)*cos(2.0_dp*pi*u2)
!!$    !tq2 = sqrt(u1)*sin(2.0_dp*pi*u3)
!!$    !tq3 = sqrt(u1)*cos(2.0_dp*pi*u3)
!!$
!!$    ! normalise
!!$    !norm = 1.0_dp/sqrt(tq0*tq0+tq1*tq1+tq2*tq2+tq3*tq3)
!!$    !tq0  = tq0*norm
!!$    !tq1  = tq1*norm
!!$    !tq2  = tq2*norm
!!$    !tq3  = tq3*norm
!!$
!!$    !mix = 0.05_dp/sqrt(real(ilev,kind=dp))
!!$
!!$    !tq0 = (1.0_dp-mix)*1.0 + mix*tq0
!!$    !tq1 = (1.0_dp-mix)*0.0 + mix*tq1
!!$    !tq2 = (1.0_dp-mix)*0.0 + mix*tq2
!!$    !tq3 = (1.0_dp-mix)*0.0 + mix*tq3
!!$
!!$
!!$    tq0 = cos(0.125*Pi) 
!!$    tq1 = 1.0_dp*sin(0.125*Pi) 
!!$    tq2 = 0.0_dp
!!$    tq3 = 0.0_dp
!!$
!!$    ! normalise
!!$    norm = 1.0_dp/sqrt(tq0*tq0+tq1*tq1+tq2*tq2+tq3*tq3)
!!$    tq0  = tq0*norm
!!$    tq1  = tq1*norm
!!$    tq2  = tq2*norm
!!$    tq3  = tq3*norm
!!$
!!$    ! build rotation matrix
!!$    a(1,1) = tq0**2 + tq1**2 - tq2**2 - tq3**2
!!$    a(1,2) = 2.0_dp * ( tq1*tq2 - tq0*tq3 )
!!$    a(1,3) = 2.0_dp * ( tq1*tq3 + tq0*tq2 )
!!$    a(2,1) = 2.0_dp * ( tq1*tq2 + tq0*tq3 )
!!$    a(2,2) = tq0**2 - tq1**2 + tq2**2 - tq3**2
!!$    a(2,3) = 2.0_dp * ( tq2*tq3 - tq0*tq1 )
!!$    a(3,1) = 2.0_dp * ( tq1*tq3 - tq0*tq2 )
!!$    a(3,2) = 2.0_dp * ( tq2*tq3 + tq0*tq1 )
!!$    a(3,3) = tq0**2 - tq1**2 - tq2**2 + tq3**2
!!$
!!$    !write(*,'(3F15.6)')a(1,:)
!!$    !write(*,'(3F15.6)')a(2,:)
!!$    !write(*,'(3F15.6)')a(3,:)
!!$
!!$    ! apply to everything in cell 2, including the cell itself
!!$    hmatrix(:,1,2) = matmul(a,hmatrix(:,1,2))
!!$    hmatrix(:,2,2) = matmul(a,hmatrix(:,2,2))
!!$    hmatrix(:,3,2) = matmul(a,hmatrix(:,3,2))
!!$
!!$    !write(*,'(3F15.6)')hmatrix(1,:,2)*bohr_to_ang
!!$    !write(*,'(3F15.6)')hmatrix(2,:,2)*bohr_to_ang
!!$    !write(*,'(3F15.6)')hmatrix(3,:,2)*bohr_to_ang
!!$
!!$
!!$    do iwater = 1,nwater
!!$       do ilj = 1,ljspm
!!$          ljr(:,ilj,iwater,2) = matmul(a,ljr(:,ilj,iwater,2))
!!$          ref_ljr(:,ilj,iwater,2) = ljr(:,ilj,iwater,2) !- ljr(:,1,iwater,2)
!!$       end do
!!$       do ics = 1,cspm
!!$          csr(:,ics,iwater,2) = matmul(a,csr(:,ics,iwater,2))
!!$          ref_csr(:,ics,iwater,2) = csr(:,ics,iwater,2) !- ljr(:,1,iwater,2)
!!$       end do
!!$    end do
!!$
!!$    !call io_hist_append(0)
!!$
!!$    do iloop = 1,100000
!!$
!!$       ! pick two molecules at random
!!$       iwater = int(random_uniform_random()*real(nwater,kind=dp))+1
!!$       jwater = iwater
!!$       do while (jwater==iwater)
!!$          jwater = int(random_uniform_random()*real(nwater,kind=dp))+1
!!$       end do
!!$
!!$       ! swap reference positions
!!$       tmp_ljr(:,:) = ref_ljr(:,:,iwater,2)
!!$       ref_ljr(:,:,iwater,2) = ref_ljr(:,:,jwater,2)
!!$       ref_ljr(:,:,jwater,2) = tmp_ljr(:,:)
!!$
!!$       tmp_csr(:,:) = ref_csr(:,:,iwater,2)
!!$       ref_csr(:,:,iwater,2) = ref_csr(:,:,jwater,2)
!!$       ref_csr(:,:,jwater,2) = tmp_csr(:,:)
!!$
!!$       sumdiff = difference()
!!$
!!$       if ( sumdiff > sumdiff_old ) then
!!$
!!$          ! reject
!!$          ref_ljr(:,:,jwater,2) = ref_ljr(:,:,iwater,2)
!!$          ref_ljr(:,:,iwater,2) = tmp_ljr(:,:)
!!$          ref_csr(:,:,jwater,2) = ref_csr(:,:,iwater,2)
!!$          ref_csr(:,:,iwater,2) = tmp_csr(:,:)
!!$
!!$       else
!!$
!!$          !write(*,'("Swapped ",I5," with ",I5)')iwater,jwater
!!$          sumdiff_old = sumdiff
!!$          print *,iloop,sumdiff
!!$
!!$          tmp_ljr(:,:) = ljr(:,:,iwater,2)
!!$          ljr(:,:,iwater,2) = ljr(:,:,jwater,2)
!!$          ljr(:,:,jwater,2) = tmp_ljr(:,:)
!!$
!!$          tmp_csr(:,:) = csr(:,:,iwater,2)
!!$          csr(:,:,iwater,2) = csr(:,:,jwater,2)
!!$          csr(:,:,jwater,2) = tmp_csr(:,:)
!!$
!!$
!!$       end if
!!$
!!$    end do
!!$
!!$    !do iwater = 1,nwater
!!$
!!$       !do ilj = 1,ljspm
!!$       !   ref_ljr(:,ilj,iwater,2) =  ljr(:,ilj,iwater,2)
!!$       !end do
!!$       !do ics = 1,cspm
!!$       !   ref_csr(:,ics,iwater,2) =  csr(:,ics,iwater,2)
!!$       !end do
!!$
!!$    !end do
!!$
!!$
!!$!    write(*,'("Best permutation at current cell rotation gives sumdiff = ",F15.6)')sumdiff
!!$    !if (sumdiff<gsumdiff_old) then
!!$
!!$     !  gsumdiff_old = sumdiff
!!$     !  write(*,'("Best permutation at current cell rotation gives sumdiff = ",F15.6)')sumdiff
!!$     !call io_hist_append(ilev)
!!$
!!$      ! ilev = ilev + 1
!!$
!!$    !else
!!$
!!$     !  tq1 = -tq1
!!$      ! tq2 = -tq2
!!$       !tq3 = -tq3
!!$
!!$       ! build rotation matrix
!!$       !a(1,1) = tq0**2 + tq1**2 - tq2**2 - tq3**2
!!$       !a(1,2) = 2.0_dp * ( tq1*tq2 - tq0*tq3 )
!!$       !a(1,3) = 2.0_dp * ( tq1*tq3 + tq0*tq2 )
!!$       !a(2,1) = 2.0_dp * ( tq1*tq2 + tq0*tq3 )
!!$       !a(2,2) = tq0**2 - tq1**2 + tq2**2 - tq3**2
!!$       !a(2,3) = 2.0_dp * ( tq2*tq3 - tq0*tq1 )
!!$       !a(3,1) = 2.0_dp * ( tq1*tq3 - tq0*tq2 )
!!$       !a(3,2) = 2.0_dp * ( tq2*tq3 + tq0*tq1 )
!!$       !a(3,3) = tq0**2 - tq1**2 - tq2**2 + tq3**2
!!$
!!$       !write(*,'(3F15.6)')a(1,:)
!!$       !write(*,'(3F15.6)')a(2,:)
!!$       !write(*,'(3F15.6)')a(3,:)
!!$
!!$       ! apply to everything in cell 2, including the cell itself
!!$       !hmatrix(:,1,2) = matmul(a,hmatrix(:,1,2))
!!$       !hmatrix(:,2,2) = matmul(a,hmatrix(:,2,2))
!!$       !hmatrix(:,3,2) = matmul(a,hmatrix(:,3,2))
!!$
!!$       !write(*,'(3F15.6)')hmatrix(1,:,2)*bohr_to_ang
!!$       !write(*,'(3F15.6)')hmatrix(2,:,2)*bohr_to_ang
!!$       !write(*,'(3F15.6)')hmatrix(3,:,2)*bohr_to_ang
!!$
!!$
!!$       !do iwater = 1,nwater
!!$       !   do ilj = 1,ljspm
!!$       !      ljr(:,ilj,iwater,2) = matmul(a,ljr(:,ilj,iwater,2))
!!$       !      ref_ljr(:,ilj,iwater,2) = ljr(:,ilj,iwater,2) - ljr(:,1,iwater,2)
!!$       !   end do
!!$       !   do ics = 1,cspm
!!$       !      csr(:,ics,iwater,2) = matmul(a,csr(:,ics,iwater,2))
!!$       !      ref_csr(:,ics,iwater,2) = csr(:,ics,iwater,2) - ljr(:,1,iwater,2)
!!$       !   end do
!!$       !end do
!!$
!!$
!!$    !end if
!!$
!!$
!!$    ! end do
!!$
!!$    deallocate(tmp_ljr,tmp_csr,stat=ierr(1))
!!$    if(ierr(1)/=0) stop 'Error releasing memory in renumber'
!!$
!!$
!!$  contains 
!!$
!!$    real(kind=dp) function difference()
!!$
!!$      implicit none
!!$      real(kind=dp),dimension(3) :: vec,vec1
!!$
!!$      integer :: ilj,ics,imol
!!$
!!$      difference = 0.0_dp
!!$
!!$      do imol = 1,nwater
!!$
!!$         !do ilj = 1,ljspm
!!$         !   vec = ref_ljr(:,ilj,imol,1) - ref_ljr(:,ilj,imol,2)
!!$         !   difference = difference + dot_product(vec,vec)
!!$         !end do
!!$
!!$         !do ics = 1,cspm
!!$         !   vec = ref_csr(:,ics,imol,1) - ref_csr(:,ics,imol,2)
!!$         !   difference = difference + dot_product(vec,vec)             
!!$         !end do
!!$
!!$         vec  = ref_ljr(:,1,imol,1) - ref_csr(:,1,imol,1)
!!$         vec1 = ref_ljr(:,1,imol,2) - ref_csr(:,1,imol,2)
!!$
!!$         difference = difference + abs(dot_product(vec,vec1))
!!$
!!$      end do
!!$
!!$    end function difference
!!$
!!$
!!$  end subroutine renumber

end module init
