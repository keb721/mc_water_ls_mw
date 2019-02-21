! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: Git -*-
!=============================================================================!
!                                 U T I L                                     !
!=============================================================================!
! Contains ancilliary utility routines used throughout the code. Many of      !
! these will be inlined manually when used for efficiency.                    !
!=============================================================================!

module util

  use constants
  implicit none

  contains

  real(kind=dp) function util_determinant(matrix)
    !------------------------------------------------------------------------------!
    ! Computes the determinant of a 3x3 matrix.                                    !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------!     

    implicit none
    real(kind=dp),dimension(3,3) :: matrix

    real(kind=dp) :: det

    det = matrix(1,1)*(matrix(2,2)*matrix(3,3) - &
          matrix(2,3)*matrix(3,2))
    det = det - matrix(1,2)*(matrix(2,1)*matrix(3,3) - &
          matrix(2,3)*matrix(3,1))
    det = det + matrix(1,3)*(matrix(2,1)*matrix(3,2) - &
          matrix(2,2)*matrix(3,1))

    util_determinant = det

    !if ( det < tiny(0.0_dp) ) stop 'Error in util_determinant'

    return

  end function util_determinant

  subroutine util_recipmatrix(hmatrix,recip_matrix)
    !------------------------------------------------------------------------------!
    ! Calculates the matrix of reciprocal lattive vectors from the h_matrix        !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------! 
    implicit none
    real(kind=dp),dimension(3,3),intent(in)  :: hmatrix
    real(kind=dp),dimension(3,3),intent(out) :: recip_matrix
    real(kind=dp) :: vol

    ! invert hmatrix to get recip_matrix
    recip_matrix(1,1)=hmatrix(2,2)*hmatrix(3,3)-hmatrix(2,3)*hmatrix(3,2)
    recip_matrix(1,2)=hmatrix(2,3)*hmatrix(3,1)-hmatrix(2,1)*hmatrix(3,3)
    recip_matrix(1,3)=hmatrix(2,1)*hmatrix(3,2)-hmatrix(2,2)*hmatrix(3,1)
    
    recip_matrix(2,1)=hmatrix(1,3)*hmatrix(3,2)-hmatrix(1,2)*hmatrix(3,3)
    recip_matrix(2,2)=hmatrix(1,1)*hmatrix(3,3)-hmatrix(1,3)*hmatrix(3,1)
    recip_matrix(2,3)=hmatrix(1,2)*hmatrix(3,1)-hmatrix(1,1)*hmatrix(3,2)
    
    recip_matrix(3,1)=hmatrix(1,2)*hmatrix(2,3)-hmatrix(1,3)*hmatrix(2,2)
    recip_matrix(3,2)=hmatrix(1,3)*hmatrix(2,1)-hmatrix(1,1)*hmatrix(2,3)
    recip_matrix(3,3)=hmatrix(1,1)*hmatrix(2,2)-hmatrix(1,2)*hmatrix(2,1)
    
    ! Calculte cell volume
    vol =hmatrix(1,1)*recip_matrix(1,1) + &
         hmatrix(1,2)*recip_matrix(1,2) + &
         hmatrix(1,3)*recip_matrix(1,3)

    ! Scale reciprocal lattice by 2*pi/volume
    recip_matrix(:,:)=recip_matrix(:,:)*2.0_dp*Pi/vol

    return

  end subroutine util_recipmatrix

  subroutine util_hmatrix_to_abc(hmatrix,alength,blength,clength,alpha,beta,gamma)
    !------------------------------------------------------------------------------!
    ! Calculates A, B, C, alpha, beta, gamma cell lengths and angles from a matrix !
    ! of cell vectors.                                                             !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2015                                                     !
    !------------------------------------------------------------------------------!     
    implicit none
    
    real(kind=dp),dimension(3,3),intent(in) :: hmatrix
    real(kind=dp),intent(out) :: alength,blength,clength
    real(kind=dp),intent(out) :: alpha,beta,gamma
    
    alength= sqrt(dot_product(hmatrix(:,1),hmatrix(:,1)))
    blength= sqrt(dot_product(hmatrix(:,2),hmatrix(:,2)))
    clength= sqrt(dot_product(hmatrix(:,3),hmatrix(:,3)))
    
    alpha   = acos(dot_product(hmatrix(:,1),hmatrix(:,3))/(alength*clength))
    beta    = acos(dot_product(hmatrix(:,2),hmatrix(:,3))/(blength*clength))
    gamma   = acos(dot_product(hmatrix(:,1),hmatrix(:,2))/(alength*blength))
    
    alpha = alpha*180.0_dp/Pi
    beta  = beta*180.0_dp/Pi
    gamma = gamma*180.0_dp/Pi
   
    return
        
  end subroutine util_hmatrix_to_abc

  subroutine util_images(v_ss,hmatrix,recip_matrix)
    !------------------------------------------------------------------------------!
    ! Computes the minimum image of a vector v_ss - Note that this is inlined for  !
    ! efficiency in most cases.                                                    !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------! 
    implicit none
    real(kind=dp),dimension(3),intent(inout) :: v_ss
    real(kind=dp),dimension(3,3),intent(in)  :: hmatrix,recip_matrix

    real(kind=dp) :: sx,sy,sz
    real(kind=dp) :: ssx,ssy,ssz

    integer :: idim

    ! compute fractional co-ordinates
    sx = recip_matrix(1,1)*v_ss(1) + &
         recip_matrix(2,1)*v_ss(2) + &
         recip_matrix(3,1)*v_ss(3)
    sy = recip_matrix(1,2)*v_ss(1) + &
         recip_matrix(2,2)*v_ss(2) + &
         recip_matrix(3,2)*v_ss(3)  
    sz = recip_matrix(1,3)*v_ss(1) + &
         recip_matrix(2,3)*v_ss(2) + &
         recip_matrix(3,3)*v_ss(3) 


    sx = sx*0.5_dp*invPi 
    sy = sy*0.5_dp*invPi
    sz = sz*0.5_dp*invPi 

    ! apply boundary conditions
    ssx = sx - floor(sx+0.5_dp,kind=dp)
    ssy = sy - floor(sy+0.5_dp,kind=dp)
    ssz = sz - floor(sz+0.5_dp,kind=dp)
    

    ! scale back up
    do idim=1,3
       v_ss(idim) = hmatrix(idim,1)*ssx + &
                    hmatrix(idim,2)*ssy + &
                    hmatrix(idim,3)*ssz
    end do   


    return

  end subroutine util_images

end module util


