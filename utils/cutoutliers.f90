subroutine cutoutliers(npt,x,y,yerr,xpos,ypos,xnep,ynep)
use precision
implicit none
!import variable
integer :: npt
real(double), dimension(:) :: x,y,yerr,xpos,ypos,xnep,ynep
!local variables
integer :: i,j
integer, allocatable, dimension(:) :: icut
real(double) :: threshold,vp,vm

threshold=0.0005
allocate(icut(npt))
icut=0 !intialize to keep all data

do i=2,npt-1
   vp=y(i)-y(i+1)
   vm=y(i)-y(i-1)
   if((abs(vp).gt.threshold).and.(abs(vm).gt.threshold).and.(vp/vm.gt.0))then
      icut(i)=1 !cut
   endif
enddo

j=0
do i=1,npt
   if(icut(i).eq.0)then
      j=j+1
      x(j)=x(i)
      y(j)=y(i)
      yerr(j)=yerr(i)
      xpos(j)=xpos(i)
      ypos(j)=ypos(i)
      xnep(j)=xnep(i)
      ynep(j)=ynep(i)
   endif
enddo
npt=j

return
end
