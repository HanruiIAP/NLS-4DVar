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
  do i=1,nrens
  call tinteg_rk4(ktcyc,x,x)
  AA(:,i)=x
  enddo
end subroutine sample
