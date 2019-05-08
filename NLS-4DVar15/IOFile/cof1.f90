program  main
integer,parameter::ndim=40
real::x1(ndim),x2(ndim),x3(ndim),d,c,sum1,sum2
integer::yr,day,i,l,it
open(11,file='true.txt')
open(20,file='Xbinner.txt')
open(12,file='ass.txt')
open(13,file='./rmse.txt')
open(14,file='mean.txt')

print*,'=================================='
do i=1,101
  read(11,*)(x1(l),l=1,ndim)
  read(12,*)(x2(l),l=1,ndim)
  read(20,*)(x3(l),l=1,ndim)
  call rmse(ndim,x1,x2,c)  
  call rmse(ndim,x1,x3,d)  
  write(13,*),c,d
enddo
stop
do yr=1,1


  sum2=0.0
  do day=1,25
    
	sum1=0.0
    do it=1,4
    read(11,*)(x1(l),l=1,ndim)
    read(12,*)(x2(l),l=1,ndim)
    !call corrcoef(ndim,ndim,x1,x2,c)
	!print*,c
	call rmse(ndim,x1,x2,c)
	print*,c
	sum1=sum1+c
	enddo
	sum1=sum1/4
    sum2=sum2+sum1
	write(13,*)sum1
  enddo
  sum2=sum2/365
  write(14,*)sum2
enddo
close(11)
close(12)
close(13)



end program main

subroutine corrcoef(n,cn,a,b,c)
integer cn,i,n
real a(n),b(n),c,ma,mb,va,vb
ma=0.0
mb=0.0
do i=1,cn
   ma=ma+a(i)
   mb=mb+b(i)
end do
ma=ma/cn
mb=mb/cn
va=0.0
vb=0.0
do i=1,cn
   va=va+(a(i)-ma)**2.0
   vb=vb+(b(i)-mb)**2.0
end do
c=0.0
do i=1,cn
   c=c+(a(i)-ma)*(b(i)-mb)
end do
c=c/sqrt(va*vb)
end 

subroutine rmse(n,a,b,c)
integer::n,i
real a(n),b(n),c
c=0.0
do i=1,n
  c=c+(a(i)-b(i))**2.0
enddo
c=sqrt(c/(n-1))
end subroutine rmse
