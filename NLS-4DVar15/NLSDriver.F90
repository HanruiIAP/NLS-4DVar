program NLSDriver
  use common,r8=> shr_kind_r8
  use lorenz96

  implicit none

  integer          :: i,j,Imax,nass
  integer,parameter:: ndim=nx
  
  integer          :: UNIT1,UNIT2,UNIT3
  integer          :: ktoneday,ktcyc
  integer          :: it
 
  character(len=120)   :: IOFile,infile
  real(r8)             :: x(ndim),xb(ndim),AA(ndim,nrens),H(nrobs/2,ndim)
  real(r8)             :: R(nrobs),Hxb(nrobs),AA_bak(ndim,nrens)
  real(r8)             :: Hxj(nrobs,nrens),Y(nrobs),obs_y(nrobs)
  real(r8)             :: Hxb_bak(nrobs),obs_y_bak(nrobs),Hxj_bak(nrobs,nrens)
  real(r8)             :: AAA(nrobs,nrens),ST(nrens),V(nrens,nrens)
  real(r8)             :: xa(ndim),Hxa(nrobs),Hxa_bak(nrobs)
  real(r8)             :: Psi_AB(nrens,nrobs),Psi_BA(nrens,nrobs)
  real(r8)             :: X_bb(ndim),AA1(ndim,nrens)
  !==================================================
  dt=0.005d0
  oneday=0.2d0 
  ktoneday = int(oneday/dt)
  ktcyc = ktoneday/4
  
  UNIT1=111
  IOFile='./IOFile/'
  infile=trim(IOFile)//'errini.txt'
  open(UNIT1,file=trim(infile))
  do i=1,ndim
  read(UNIT1,*) x(i)
  enddo
  close(UNIT1)

  UNIT1=111
  IOFile='./IOFile/'
  infile=trim(IOFile)//'Obs.txt'
  open(UNIT1,file=trim(infile))

  UNIT2=222
  IOFile='./IOFile/'
  infile=trim(IOFile)//'ass.txt'
  open(UNIT2,file=trim(infile))

  UNIT3=333
  IOFile='./IOFile/'
  infile=trim(IOFile)//'Xbinner.txt'
  open(UNIT3,file=trim(infile))

!--------------------------------------
  xb=x    
  call sample(ndim,xb,AA)    
  !------------------------------------
  H=0.0
  do j=1,nrobs/2
    H(j,2+(j-1)*2)=1.0
  enddo
  !-------------------------------------
  do i=1,nrobs
   R(i)=1.0/(0.01)**2.0
  enddo  
  !------------------------------------
 do nass=1,25
  xa=xb
  AA_BAK=AA
  !----------------------------------------------
  x=xb    
  if(nass==1)then
    write(UNIT3,*)(x(j),j=1,ndim)
  endif

  do it=1,4
    call tinteg_rk4(ktcyc,x,x)
    write(UNIT3,*)(x(j),j=1,ndim)
    if (it==2) then
    call BRMULNew(H,X,nrobs/2,ndim,1,Hxb(1:nrobs/2))     
    endif
    if (it==4) then
    call BRMULNew(H,X,nrobs/2,ndim,1,Hxb(1+nrobs/2:nrobs))     
    endif
  enddo 
  Hxb_bak=Hxb
  !-------------------------------------
  do it=1,4
    do i=1,nrens
       x(:)=AA(:,i)
       call tinteg_rk4(ktcyc,x,x)  
       AA(:,i)=x(:)
    enddo
    if (it==2) then
    call BRMULNew(H,AA,nrobs/2,ndim,nrens,Hxj(1:nrobs/2,:))     
    endif
    if (it==4) then
    call BRMULNew(H,AA,nrobs/2,ndim,nrens,Hxj(1+nrobs/2:nrobs,:))     
    endif
    AA1=AA
  enddo
  Hxj_bak=Hxj 
  !----------------------------------------
  do it=1,4
  if (it==2) then
    read(UNIT1,*) (Y(i),i=1,nrobs)
    call BRMULNew(H,Y,nrobs/2,ndim,1,obs_y(1:nrobs/2))
  endif 
  if (it==4) then
    read(UNIT1,*) (Y(i),i=1,nrobs)
    call BRMULNew(H,Y,nrobs/2,ndim,1,obs_y(1+nrobs/2:nrobs))
  endif 
  enddo
  obs_y_bak=obs_y
  !-----------------------------------------
  obs_y=obs_y_bak-Hxb_bak

  !------------------------------------------
  do i=1,nrens
    AA(:,i)=AA_bak(:,i)-xb(:)
    AAA(:,i)=Hxj_bak(:,i)-Hxb_bak(:)
  enddo
  !------------------------------------------
  call OPs_POD(AAA,V,ST)
  call MPs_POD(AA,V,ST,ndim)
  !------------------------------------------  
  call NLSCore(AAA,R,Psi_AB,Psi_BA)
  !===================================================
  xa=xb  
  do Imax=1,1
    write(*,*) Imax,"================================"
    x=xa       
    do it=1,4      
    call tinteg_rk4(ktcyc,x,x)
    if (it==2) then
    call BRMULNew(H,x,nrobs/2,ndim,1,Hxa(1:nrobs/2))     
    endif
    if (it==4) then
    call BRMULNew(H,x,nrobs/2,ndim,1,Hxa(1+nrobs/2:nrobs))     
    endif
    enddo  
    Hxa_bak=Hxa    
    !--------------------------------------------------
    Hxa=Hxa_bak-Hxb_bak     
    !--------------------------------------------------
    obs_y=obs_y_bak-Hxb_bak
    call NLSUpdate(ndim,obs_y,AA,Psi_AB,Psi_BA,Hxa,X_bb)    
    xa=xa+X_bb     
  enddo
  !-----------------------------------------------------
  do i=1,nrens
     AAA(:,i)=Hxj_bak(:,i)-Hxa_bak(:)
  enddo  
  call RSEnUpdate(R,AAA,V)
  !------------------------------------------------------
  x=xa
  if(nass==1)then
    write(UNIT2,*)(x(j),j=1,ndim)
  endif
  do it=1,4
    call tinteg_rk4(ktcyc,x,x)
    write(UNIT2,*)(x(j),j=1,ndim)
  enddo
  !------------------------------------------------------
  xb=x
  !------------------------------------------------------
  do i=1,nrens
    AA1(:,i)=AA1(:,i)-xb
  enddo   
  call BRMULNew(AA1,V,ndim,nrens,nrens,AA)
  do i=1,nrens
    AA(:,i)=AA(:,i)+xb(:)
  enddo
  !-------------------------------------------------------     
 enddo
close(UNIT1)
close(UNIT2)
close(UNIT3)
end program NLSDriver
