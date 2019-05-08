subroutine MPs_POD(AA,V,ST,ndim)
!--------------------------------------------------------------------------
use common,r8=> shr_kind_r8
implicit none
!--------------------------------------------------------------------------
real(r8) :: AA(ndim,nrens),V(nrens,nrens)
real(r8) :: ST(nrens),AA_BAK(ndim,nrens)
integer  :: i,ndim
!--------------------------------------------------------------------------
  AA_BAK=AA
  call BRMULNew(AA_BAK,V,ndim,nrens,nrens,AA)
  do i=1,nrens
  if(ST(i)>0.01)then
     AA(:,i)=AA(:,i)/sqrt(ST(i))
  endif
  enddo
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
end subroutine MPs_POD
