subroutine NLSUpdate(ndim,Obs_Y,AA,Psi_AB,Psi_BA,HX_a,X_bb)
use common,r8=> shr_kind_r8
implicit none

integer             :: ndim,nrobs1
real(r8)            :: Obs_Y(nrobs),X_bb(ndim),Psi_BA(nrens,nrobs)
real(r8)            :: Psi_AB(nrens,nrobs),V(nrens,nrens),ST(nrens)
real(r8)            :: AA(ndim,nrens),PP(ndim,nrobs),AA_BAK(ndim,nrens)
real(r8)            :: TEMP1(nrens),PP1(ndim,nrobs),X_bbb(ndim),HX_a(nrobs)
real                :: Rad_xy,Rou_xy
integer             :: i,j,k,L,s,s1,i1,j1
real(r8)            :: Obs_Y_bak(nrobs)
!--------------------------------------------------------------------------
nrobs1=nrobs/2
!--------------------------------------------------------------------------
CALL BRMULNew(AA,Psi_AB,ndim,nrens,nrobs,PP)
CALL BRMULNew(AA,Psi_BA,ndim,nrens,nrobs,PP1)

do i=1,ndim
  do s=1,nrobs1
     Rad_xy=real(abs(i-s*2))/4.0
     CALL LocRou(Rad_xy,Rou_xy)
     PP(i,s)=PP(i,s)*Rou_xy
     PP(i,s+nrobs1)=PP(i,s+nrobs1)*Rou_xy
     PP1(i,s)=PP1(i,s)*Rou_xy
     PP1(i,s+nrobs1)=PP1(i,s+nrobs1)*Rou_xy
  enddo
enddo
Obs_Y_bak=Obs_Y
Obs_Y=Obs_Y_bak-HX_a

CALL BRMULNew(PP,Obs_Y,ndim,nrobs,1,X_bb)
CALL BRMULNew(PP1,HX_a,ndim,nrobs,1,X_bbb)
X_bb=X_bb+X_bbb
!--------------------------------------------------------------------------
end subroutine NLSUpdate
