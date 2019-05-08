subroutine OPs_POD(AA,V,ST)
!--------------------------------------------------------------------------
use common,r8=> shr_kind_r8
implicit none
!--------------------------------------------------------------------------
real(r8) :: AA(nrobs,nrens),V(nrens,nrens),AA_BAK(nrobs,nrens)
real(r8) :: TT(nrens,nrens),V_BAK(nrens),ST(nrens)
real(r8) :: EPS,PA,sum1
integer  :: i,j,s,k,L
!--------------------------------------------------------------------------
           AA_BAK=AA
           do i=1,nrens
              do j=1,nrens
                 TT(i,j)=0.0
                 do k=1,nrobs
                    TT(i,j)=TT(i,j)+AA(k,i)*AA(k,j)
                 enddo
              enddo
           enddo
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
           call BRMULNew(AA_BAK,V,nrobs,nrens,nrens,AA)
           do i=1,nrens
              if(TT(i,i)>0.01)then
               AA(:,i)=AA(:,i)/sqrt(TT(i,i))
              endif
              ST(i)=TT(i,i)
           enddo
!--------------------------------------------------------------------------
end subroutine OPs_POD
