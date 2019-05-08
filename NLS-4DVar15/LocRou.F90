subroutine LocRou(Rad,Rou)
implicit none
real::Rad,Rou
 
if(Rad<=1)then
   Rou=-1.0/4*Rad**5.0+1.0/2*Rad**4.0+5./8*Rad**3.0-5.0/3*Rad**2.0+1.0
endif
if(Rad>1.and.Rad<=2.0)then
   Rou=1.0/12*Rad**5.0-1.0/2*Rad**4.0+5.0/8*Rad**3.0+5.0/3*Rad**2.0&
     &-5.0*Rad+4.0-2.0/3/Rad
endif
if(Rad>2)then
  Rou=0.0
endif
end subroutine LocRou
