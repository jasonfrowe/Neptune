subroutine lininterp(x,y,n,xin,yout)
!linear interpolation.  X is expected to be in order.
use precision
implicit none
integer :: n,i
real(double) :: x(n),y(n),xin,yout

if(xin.le.x(1))then
  i=1
  yout=y(i)+(xin-x(i))/(x(i+1)-x(i))*(y(i+1)-y(i))
elseif(xin.ge.x(n))then
  i=n-1
  yout=y(i)+(xin-x(i))/(x(i+1)-x(i))*(y(i+1)-y(i))
else
  do i=1,n-1
    if((xin.ge.x(i)).and.(xin.lt.x(i+1)))then
      yout=y(i)+(xin-x(i))/(x(i+1)-x(i))*(y(i+1)-y(i))
    endif
  enddo
endif

return
end
