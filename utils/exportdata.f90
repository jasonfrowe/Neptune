!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine exportdata(npt,x,y,yerr,r,sp,xpos,ypos,xnep,ynep,minx,mean)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!import vars
integer :: npt
real(double) :: minx,mean
real(double), dimension(npt) :: x,y,yerr,r,sp,xpos,ypos,xnep,ynep
!local vars
integer :: i
real(double) :: newdata

do i=1,npt
   newdata=y(i)-r(i)+sp(i) !modeled dataset.
   write(6,500) x(i)+minx,mean*(newdata+1.0d0),mean*yerr(i),            &
    mean*(y(i)+1.0d0),mean*(r(i)+1.0d0),mean*(sp(i)+1.0d0),xpos(i),     &
    ypos(i),xnep(i),ynep(i)
   500 format(10(1PE17.10,1X))
enddo

return
end
