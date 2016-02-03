!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine exportdata(npt,x,y,yerr,r,sp,xpos,ypos,xnep,ynep)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!import vars
integer :: npt
real(double), dimension(npt) :: x,y,yerr,r,sp,xpos,ypos,xnep,ynep
!local vars
integer :: i
real(double) :: newdata

do i=1,npt
   newdata=y(i)-r(i)+sp(i) !modeled dataset.
   write(6,500) x(i),newdata,yerr(i),y(i),r(i),sp(i),xpos(i),ypos(i),&
     xnep(i),ynep(i)
   500 format(10(1PE17.10,1X))
enddo

return
end
