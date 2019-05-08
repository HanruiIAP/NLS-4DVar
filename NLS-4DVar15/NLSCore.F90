subroutine NLSCore(AAA,R,Psi_AB,Psi_BA)
use common,r8=> shr_kind_r8
 
implicit none

real(r8) :: AAA(nrobs,nrens),R(nrobs)
real(r8) :: Psi_A(nrens,nrens),Psi_AB(nrens,nrobs)
real(r8) :: TT(nrens,nrens),TRAN_LHA(nrens,nrobs)

!---------------------------------------------------------------------------
real(r8) :: SS(nrens,nrens),TRAN_AAA(nrens,nrobs)
real(r8) :: Psi_BA(nrens,nrobs),Psi_B(nrens,nrens)

integer  :: i,j,k,L
!---------------------------------------------------------------------------

!--------------------------------------------------------------------------
  Psi_A=0.0
  do i=1,nrens
    Psi_A(i,i)=real(nrens-1)
  enddo
  CALL TRANS_MATRIX(AAA,TRAN_LHA,nrobs,nrens)
  TRAN_AAA=TRAN_LHA
 
  do i=1,nrobs
    TRAN_LHA(:,i)=TRAN_LHA(:,i)*R(i)
  enddo
 
  do i=1,nrens
     do j=1,nrens
        TT(i,j)=0.0
        SS(i,j)=0.0
        do k=1,nrobs
          TT(i,j)=TT(i,j)+TRAN_LHA(i,k)*AAA(k,j)
          SS(i,j)=SS(i,j)+TRAN_AAA(i,k)*AAA(k,j)
        enddo
     enddo
     SS(i,i) = SS(i,i) + 0.00000001d0
  enddo

  Psi_A=Psi_A+TT
  CALL BSSGJ(Psi_A,nrens,L)
  CALL BSSGJ(SS,nrens,L)

  CALL BRMULNew(Psi_A,TRAN_LHA,nrens,nrens,nrobs,Psi_AB)
    
  Psi_A=-real(nrens-1)*Psi_A

  CALL BRMULNew(Psi_A,SS,nrens,nrens,nrens,Psi_B)
  CALL BRMULNew(Psi_B,TRAN_AAA,nrens,nrens,nrobs,Psi_BA)
!-------------------------------------------------------------------------
end subroutine NLSCore
