subroutine RSEnUpdate(R,AAA,V)
!--------------------------------------------------------------------------
use common,r8=> shr_kind_r8
!--------------------------------------------------------------------------
real(r8) :: AAA(nrobs,nrens),R(nrobs)
real(r8) :: Psi_A(nrens,nrens)
real(r8) :: TRAN_LHA(nrens,nrobs)

real(r8) :: V(nrens,nrens)
real(r8) :: TT(nrens,nrens),V_BAK(nrens)
real(r8) :: EPS,PA,sum1
integer  :: i,j,s,k,L
!--------------------------------------------------------------------------
  Psi_A=0.0
  do i=1,nrens
    Psi_A(i,i)=real(nrens-1)
  enddo
  CALL TRANS_MATRIX(AAA,TRAN_LHA,nrobs,nrens)
  do i=1,nrobs
    TRAN_LHA(:,i)=TRAN_LHA(:,i)*R(i)
  enddo
  do i=1,nrens
     do j=1,nrens
        TT(i,j)=0.0
        do k=1,nrobs
          TT(i,j)=TT(i,j)+TRAN_LHA(i,k)*AAA(k,j)
        enddo
     enddo
  enddo
  Psi_A=Psi_A+TT
  CALL BSSGJ(Psi_A,nrens,L)
  TT=Psi_A*real(nrens-1)
!--------------------------------------------------------------------------
  EPS=0.000001D01
  CALL CJCBI(TT,nrens,EPS,V,L)
  do i=1,nrens-1
       s=i
       do j=i+1, nrens
         if(TT(j,j)>TT(s,s))s=j
       end do
       PA=TT(i,i)
       TT(i,i)=TT(s,s)
       TT(s,s)=PA
       V_BAK(:)=V(:,i)
       V(:,i)=V(:,s)
       V(:,s)=V_BAK(:)
 enddo
 do i=1,nrens
   if(TT(i,i)>0.01)then
     V(:,i)=V(:,i)*sqrt(TT(i,i))
   endif
 enddo
!--------------------------------------------------------------------------
end subroutine RSEnUpdate
