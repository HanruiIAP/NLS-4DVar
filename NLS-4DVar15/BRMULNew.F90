SUBROUTINE BRMULNew(A,B,M,N,K,C)
integer::M,N,K
double precision:: A(M,N),B(N,K),C(M,K)
double precision:: TEMP
integer::i,j,l

C=0.0


do j=1,K
  do l=1,N
     TEMP=B(l,j)
     do i=1,M
       C(i,j)=C(i,j)+TEMP*A(i,l)
     enddo
  enddo
enddo

END SUBROUTINE BRMULNew



