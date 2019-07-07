!
! Copyright (C) 2016 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE ktetra
  !
  ! ... Variables used by the tetrahedron method
  ! ... Three versions are implemented: Linear, Optimized, Bloechl
  ! ... Linear and Optimized tetrahedra contributed by Mitsuaki Kawamura,
  ! ... University of Tokyo
  !
  USE kinds, ONLY: dp
  !
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  INTEGER:: &
       tetra_type = 0  ! 0 for Bloechl's correction
                       ! 1 for Linear tetrahedron method
  !                    ! 2 for Optimized tetrahedron method
  INTEGER :: &
       ntetra, &         ! number of tetrahedra
       nntetra           ! k-points per tetrahedron used to compute weights
                         ! 4 for linear / 20 for optimized tetrahedron method
  INTEGER, ALLOCATABLE :: &
       tetra(:,:)        ! index of k-points in a given tetrahedron
                         ! shape (nntetra,ntetra)
  !
  REAL(dp), ALLOCATABLE :: &
       wlsm(:,:)         ! Weights for the optimized tetrahedron method
  !
  PUBLIC :: tetra, ntetra, nntetra
  PUBLIC :: tetra_init, tetra_weights, tetra_weights_only, tetra_dos_t
  PUBLIC :: opt_tetra_init, opt_tetra_weights, opt_tetra_weights_only, &
  &         opt_tetra_dos_t, opt_tetra_partialdos, tetra_type, wlsm
  PUBLIC :: deallocate_tetra
  !
CONTAINS
  !
  !-----------------------------------------------------------------------------
  SUBROUTINE tetra_init ( nsym, s, time_reversal, t_rev, at, bg, npk, &
     k1,k2,k3, nk1,nk2,nk3, nks, xk )
  !-----------------------------------------------------------------------
  !
  ! Tetrahedron method according to P. E. Bloechl et al, PRB49, 16223 (1994)
  !
  USE kinds, ONLY: DP
  IMPLICIT NONE
  ! 
  INTEGER, INTENT(IN):: nks, nsym, t_rev(48), s(3,3,48), npk, &
                        k1, k2, k3, nk1, nk2, nk3
  LOGICAL, INTENT (IN) :: time_reversal
  real(DP), INTENT(IN) :: at(3,3), bg(3,3)
  real(DP), INTENT(INOUT) :: xk(3,npk)
  !
  real(DP) :: xkr(3), deltap(3), deltam(3)
  real(DP), PARAMETER:: eps=1.0d-5
  real(DP), ALLOCATABLE :: xkg(:,:)
  INTEGER :: nkr, i,j,k, ns, n, nk, ip1,jp1,kp1, &
       n1,n2,n3,n4,n5,n6,n7,n8
  INTEGER, ALLOCATABLE:: equiv(:)
  !
  ntetra =6*nk1*nk2*nk3
  nntetra=4
  IF(.NOT. ALLOCATED(tetra)) ALLOCATE ( tetra (nntetra, ntetra) )
  !
  ! Re-generate a uniform grid of k-points xkg
  !
  nkr=nk1*nk2*nk3
  ALLOCATE (xkg( 3,nkr))
  ALLOCATE (equiv( nkr))
!
  DO i=1,nk1
     DO j=1,nk2
        DO k=1,nk3
           !  this is nothing but consecutive ordering
           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
           !  xkg are the components of the complete grid in crystal axis
           xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
           xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
           xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
        ENDDO
     ENDDO
  ENDDO

  !  locate k-points of the uniform grid in the list of irreducible k-points
  !  that was previously calculated

  !  bring irreducible k-points to crystal axis
  CALL cryst_to_cart (nks,xk,at,-1)
  !
  DO nk=1,nkr
     DO n=1,nks
        DO ns=1,nsym
           DO i=1,3
              xkr(i) = s(i,1,ns) * xk(1,n) + &
                       s(i,2,ns) * xk(2,n) + &
                       s(i,3,ns) * xk(3,n)
           ENDDO
           IF(t_rev(ns)==1) xkr = -xkr
           !  xkr is the n-th irreducible k-point rotated wrt the ns-th symmetry
           DO i=1,3
              deltap(i) = xkr(i)-xkg(i,nk) - nint (xkr(i)-xkg(i,nk) )
              deltam(i) = xkr(i)+xkg(i,nk) - nint (xkr(i)+xkg(i,nk) )
           ENDDO
           !  deltap is the difference vector, brought back in the first BZ
           !  deltam is the same but with k => -k (for time reversal)
           IF ( sqrt ( deltap(1)**2 + &
                       deltap(2)**2 + &
                       deltap(3)**2 ) < eps .or. ( time_reversal .and. &
                sqrt ( deltam(1)**2 +  &
                       deltam(2)**2 +  &
                       deltam(3)**2 ) < eps ) ) THEN
              !  equivalent irreducible k-point found
              equiv(nk) = n
              GOTO 15
           ENDIF
        ENDDO
     ENDDO
     !  equivalent irreducible k-point found - something wrong
     CALL errore('tetra_init','cannot locate  k point',nk)
15   CONTINUE
  ENDDO

  DO n=1,nks
     DO nk=1,nkr
        IF (equiv(nk)==n) GOTO 20
     ENDDO
     !  this failure of the algorithm may indicate that the displaced grid
     !  (with k1,k2,k3.ne.0) does not have the full symmetry of the lattice
     CALL errore('tetra_init','cannot remap grid on k-point list',n)
20   CONTINUE
  ENDDO

  !  bring irreducible k-points back to cartesian axis
  CALL cryst_to_cart (nks,xk,bg, 1)

  !  construct tetrahedra

  DO i=1,nk1
     DO j=1,nk2
        DO k=1,nk3
           !  n1-n8 are the indices of k-point 1-8 forming a cube
           ip1 = mod(i,nk1)+1
           jp1 = mod(j,nk2)+1
           kp1 = mod(k,nk3)+1
           n1 = (  k-1) + (  j-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n2 = (  k-1) + (  j-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           n3 = (  k-1) + (jp1-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n4 = (  k-1) + (jp1-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           n5 = (kp1-1) + (  j-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n6 = (kp1-1) + (  j-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           n7 = (kp1-1) + (jp1-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n8 = (kp1-1) + (jp1-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           !  there are 6 tetrahedra per cube (and nk1*nk2*nk3 cubes)
           n  = 6 * ( (k-1) + (j-1)*nk3 + (i-1)*nk3*nk2 )

           tetra (1,n+1) = equiv(n1)
           tetra (2,n+1) = equiv(n2)
           tetra (3,n+1) = equiv(n3)
           tetra (4,n+1) = equiv(n6)

           tetra (1,n+2) = equiv(n2)
           tetra (2,n+2) = equiv(n3)
           tetra (3,n+2) = equiv(n4)
           tetra (4,n+2) = equiv(n6)

           tetra (1,n+3) = equiv(n1)
           tetra (2,n+3) = equiv(n3)
           tetra (3,n+3) = equiv(n5)
           tetra (4,n+3) = equiv(n6)

           tetra (1,n+4) = equiv(n3)
           tetra (2,n+4) = equiv(n4)
           tetra (3,n+4) = equiv(n6)
           tetra (4,n+4) = equiv(n8)

           tetra (1,n+5) = equiv(n3)
           tetra (2,n+5) = equiv(n6)
           tetra (3,n+5) = equiv(n7)
           tetra (4,n+5) = equiv(n8)

           tetra (1,n+6) = equiv(n3)
           tetra (2,n+6) = equiv(n5)
           tetra (3,n+6) = equiv(n6)
           tetra (4,n+6) = equiv(n7)
        ENDDO
     ENDDO
  ENDDO

  !  check

  DO n=1,ntetra
     DO i=1,nntetra
        IF ( tetra(i,n)<1 .or. tetra(i,n)>nks ) &
             CALL errore ('tetra_init','something wrong',n)
     ENDDO
  ENDDO

  DEALLOCATE(equiv)
  DEALLOCATE(xkg)

  RETURN
  END SUBROUTINE tetra_init
  !
!
!--------------------------------------------------------------------
subroutine tetra_weights (nks, nspin, nbnd, nelec, ntetra, tetra, et, &
     ef, wg, is, isk )
  !--------------------------------------------------------------------
  !
  ! ... calculates Ef and weights with the tetrahedron method (P.E.Bloechl)
  ! ... Wrapper routine: computes first Ef, then the weights
  !
  USE kinds
  implicit none
  ! I/O variables
  integer, intent(in) :: nks, nspin, is, isk(nks), nbnd, ntetra, &
       tetra (4, ntetra)
  real(DP), intent(in) :: et (nbnd, nks), nelec
  ! wg must be (inout) and not (out) because if is/=0 only terms for
  ! spin=is are initialized; the remaining terms should be kept, not lost
  real(DP), intent(inout) :: wg (nbnd, nks)
  real(DP), intent(out) :: ef
  ! local variables
  real(DP), external :: efermit

  ! Calculate the Fermi energy ef

  ef = efermit (et, nbnd, nks, nelec, nspin, ntetra, tetra, is, isk)
  !
  ! if efermit cannot find a sensible value for Ef it returns Ef=1d10
  !
  if (abs(ef) > 1.0d8) call errore ('tetra_weights', 'bad Fermi energy ',1)
  !
  CALL tetra_weights_only (nks, nspin, is, isk, nbnd, nelec, ntetra, &
       tetra, et, ef, wg)
  !
  return
end subroutine tetra_weights

!--------------------------------------------------------------------
subroutine tetra_weights_only (nks, nspin, is, isk, nbnd, nelec, ntetra, &
     tetra, et, ef, wg)
  !--------------------------------------------------------------------
  !
  ! ... calculates weights with the tetrahedron method (P.E.Bloechl)
  ! ... Fermi energy has to be calculated in previous step
  ! ... Generalization to noncollinear case courtesy of Iurii Timrov

  USE kinds
  implicit none
  ! I/O variables
  integer, intent(in) :: nks, nspin, is, isk(nks), nbnd, ntetra, &
       tetra (4, ntetra)
  real(DP), intent(in) :: et (nbnd, nks), nelec, ef
  ! wg must be (inout) and not (out) because if is/=0 only terms for
  ! spin=is are initialized; the remaining terms should be kept, not lost
  real(DP), intent(inout) :: wg (nbnd, nks)
  ! local variables
  real(DP) :: e1, e2, e3, e4, c1, c2, c3, c4, etetra (4), dosef
  integer :: ik, ibnd, nt, nk, ns, i, kp1, kp2, kp3, kp4, itetra (4)
  integer :: nspin_lsda
  !
  do ik = 1, nks
     if (is /= 0) then
        if (isk(ik) .ne. is) cycle
     end if
     do ibnd = 1, nbnd
        wg (ibnd, ik) = 0.d0
     enddo
  enddo

  IF ( nspin == 2 ) THEN
     nspin_lsda = 2
  ELSE
     nspin_lsda = 1
  END IF

  do ns = 1, nspin_lsda
     if (is /= 0) then
        if (ns .ne. is) cycle
     end if
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     if (ns.eq.1) then
        nk = 0
     else
        nk = nks / 2
     endif
     do nt = 1, ntetra
        do ibnd = 1, nbnd
           !
           ! etetra are the energies at the vertexes of the nt-th tetrahedron
           !
           do i = 1, 4
              etetra (i) = et (ibnd, tetra (i, nt) + nk)
           enddo
           itetra (1) = 0
           call hpsort (4, etetra, itetra)
           !
           ! ...sort in ascending order: e1 < e2 < e3 < e4
           !
           e1 = etetra (1)
           e2 = etetra (2)
           e3 = etetra (3)
           e4 = etetra (4)
           !
           ! kp1-kp4 are the irreducible k-points corresponding to e1-e4
           !
           kp1 = tetra (itetra (1), nt) + nk
           kp2 = tetra (itetra (2), nt) + nk
           kp3 = tetra (itetra (3), nt) + nk
           kp4 = tetra (itetra (4), nt) + nk
           !
           ! calculate weights wg
           !
           if (ef.ge.e4) then
              wg (ibnd, kp1) = wg (ibnd, kp1) + 0.25d0 / ntetra
              wg (ibnd, kp2) = wg (ibnd, kp2) + 0.25d0 / ntetra
              wg (ibnd, kp3) = wg (ibnd, kp3) + 0.25d0 / ntetra
              wg (ibnd, kp4) = wg (ibnd, kp4) + 0.25d0 / ntetra
           elseif (ef.lt.e4.and.ef.ge.e3) then
              c4 = 0.25d0 / ntetra * (e4 - ef) **3 / (e4 - e1) / (e4 - e2) &
                   / (e4 - e3)
              dosef = 3.d0 / ntetra * (e4 - ef) **2 / (e4 - e1) / (e4 - e2) &
                   / (e4 - e3)
              wg (ibnd, kp1) = wg (ibnd, kp1) + 0.25d0 / ntetra - c4 * &
                   (e4 - ef) / (e4 - e1) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
                   (ibnd, kp1) ) / 40.d0
              wg (ibnd, kp2) = wg (ibnd, kp2) + 0.25d0 / ntetra - c4 * &
                   (e4 - ef) / (e4 - e2) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
                   (ibnd, kp2) ) / 40.d0
              wg (ibnd, kp3) = wg (ibnd, kp3) + 0.25d0 / ntetra - c4 * &
                   (e4 - ef) / (e4 - e3) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
                   (ibnd, kp3) ) / 40.d0
              wg (ibnd, kp4) = wg (ibnd, kp4) + 0.25d0 / ntetra - c4 * &
                   (4.d0 - (e4 - ef) * (1.d0 / (e4 - e1) + 1.d0 / (e4 - e2) &
                   + 1.d0 / (e4 - e3) ) ) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * &
                   et (ibnd, kp4) ) / 40.d0
           elseif (ef.lt.e3.and.ef.ge.e2) then
              c1 = 0.25d0 / ntetra * (ef - e1) **2 / (e4 - e1) / (e3 - e1)
              c2 = 0.25d0 / ntetra * (ef - e1) * (ef - e2) * (e3 - ef) &
                   / (e4 - e1) / (e3 - e2) / (e3 - e1)
              c3 = 0.25d0 / ntetra * (ef - e2) **2 * (e4 - ef) / (e4 - e2) &
                   / (e3 - e2) / (e4 - e1)
              dosef = 1.d0 / ntetra / (e3 - e1) / (e4 - e1) * (3.d0 * &
                   (e2 - e1) + 6.d0 * (ef - e2) - 3.d0 * (e3 - e1 + e4 - e2) &
                   * (ef - e2) **2 / (e3 - e2) / (e4 - e2) )
              wg (ibnd, kp1) = wg (ibnd, kp1) + c1 + (c1 + c2) * (e3 - ef) &
                   / (e3 - e1) + (c1 + c2 + c3) * (e4 - ef) / (e4 - e1) + dosef * &
                   (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp1) ) / 40.d0
              wg (ibnd, kp2) = wg (ibnd, kp2) + c1 + c2 + c3 + (c2 + c3) &
                   * (e3 - ef) / (e3 - e2) + c3 * (e4 - ef) / (e4 - e2) + dosef * &
                   (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp2) ) / 40.d0
              wg (ibnd, kp3) = wg (ibnd, kp3) + (c1 + c2) * (ef - e1) &
                   / (e3 - e1) + (c2 + c3) * (ef - e2) / (e3 - e2) + dosef * &
                   (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp3) ) / 40.d0
              wg (ibnd, kp4) = wg (ibnd, kp4) + (c1 + c2 + c3) * (ef - e1) &
                   / (e4 - e1) + c3 * (ef - e2) / (e4 - e2) + dosef * (e1 + e2 + &
                   e3 + e4 - 4.d0 * et (ibnd, kp4) ) / 40.d0
           elseif (ef.lt.e2.and.ef.ge.e1) then
              c4 = 0.25d0 / ntetra * (ef - e1) **3 / (e2 - e1) / (e3 - e1) &
                   / (e4 - e1)
              dosef = 3.d0 / ntetra * (ef - e1) **2 / (e2 - e1) / (e3 - e1) &
                   / (e4 - e1)
              wg (ibnd, kp1) = wg (ibnd, kp1) + c4 * (4.d0 - (ef - e1) &
                   * (1.d0 / (e2 - e1) + 1.d0 / (e3 - e1) + 1.d0 / (e4 - e1) ) ) &
                   + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp1) ) / 40.d0
              wg (ibnd, kp2) = wg (ibnd, kp2) + c4 * (ef - e1) / (e2 - e1) &
                   + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp2) ) / 40.d0
              wg (ibnd, kp3) = wg (ibnd, kp3) + c4 * (ef - e1) / (e3 - e1) &
                   + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp3) ) / 40.d0
              wg (ibnd, kp4) = wg (ibnd, kp4) + c4 * (ef - e1) / (e4 - e1) &
                   + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp4) ) / 40.d0
           endif
        enddo
     enddo


  enddo
  ! add correct spin normalization (2 for LDA, 1 for all other cases)
  IF ( nspin == 1 ) wg (:,1:nks) = wg (:,1:nks) * 2.d0
  !
  return
end subroutine tetra_weights_only
!
!--------------------------------------------------------------------
subroutine tetra_dos_t (et, nspin, nbnd, nks, e, dost)
  !------------------------------------------------------------------
  !
  USE kinds, only : DP
  implicit none
  integer, intent(in) :: nspin, nbnd, nks

  real(DP), intent(in) :: et (nbnd, nks), e
  REAL(dp), INTENT(OUT):: dost (2)
  !
  integer :: itetra (4), nk, ns, nt, ibnd, i
  real(DP) :: etetra (4), e1, e2, e3, e4
  integer :: nspin0

  if (nspin==4) then
     nspin0=1
  else 
     nspin0=nspin
  endif
  do ns = 1, nspin0
     dost (ns) = 0.d0
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     if (ns.eq.1) then
        nk = 0
     else
        nk = nks / 2
     endif
     do nt = 1, ntetra
        do ibnd = 1, nbnd
           ! these are the energies at the vertexes of the nt-th tetrahedron
           do i = 1, 4
              etetra (i) = et (ibnd, tetra (i, nt) + nk)
           enddo
           itetra (1) = 0
           call hpsort (4, etetra, itetra)
           e1 = etetra (1)
           e2 = etetra (2)
           e3 = etetra (3)
           e4 = etetra (4)
           if (e.lt.e4.and.e.ge.e3) then
              dost (ns) = dost (ns) + 1.d0 / ntetra * (3.0d0 * (e4 - e) **2 / &
                   (e4 - e1) / (e4 - e2) / (e4 - e3) )
           elseif (e.lt.e3.and.e.ge.e2) then
              dost (ns) = dost (ns) + 1.d0 / ntetra / (e3 - e1) / (e4 - e1) &
                   * (3.0d0 * (e2 - e1) + 6.0d0 * (e-e2) - 3.0d0 * (e3 - e1 + e4 - e2) &
                   / (e3 - e2) / (e4 - e2) * (e-e2) **2)
           elseif (e.lt.e2.and.e.gt.e1) then
              dost (ns) = dost (ns) + 1.d0 / ntetra * 3.0d0 * (e-e1) **2 / &
                   (e2 - e1) / (e3 - e1) / (e4 - e1)
           endif
        enddo

     enddo

     ! add correct spin normalization : 2 for LDA, 1 for LSDA or
     ! noncollinear calculations 

     if ( nspin == 1 ) dost (ns) = dost (ns) * 2.d0

  enddo
  return
end subroutine tetra_dos_t

!
SUBROUTINE deallocate_tetra ( )
   IF ( ALLOCATED(tetra) ) DEALLOCATE (tetra)
   IF ( ALLOCATED(wlsm ) ) DEALLOCATE (wlsm )
END SUBROUTINE deallocate_tetra
!
END MODULE ktetra
