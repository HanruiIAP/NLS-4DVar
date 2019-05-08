subroutine sample(ndim,xx,AA)
  use common,r8=> shr_kind_r8
  use lorenz96

  implicit none
  integer              :: ndim
  integer              :: i,j,ktoneday,ktcyc
  real(r8)             :: xx(ndim),AA(ndim,nrens),x(ndim),bak(ndim)
  !===============================================
  dt=0.005d0
  oneday=0.2d0 
  ktoneday = int(oneday/dt)
  ktcyc = ktoneday/4
  !================================================
  x=xx
  call tinteg_rk4(ktcyc,x,x)
  bak=xx
  do i=1,nrens
  call tinteg_rk4(ktcyc,x,x)
  AA(:,i)=x
  x=(x+2*bak)/3.0d0
  bak=x
  enddo
end subroutine sample
