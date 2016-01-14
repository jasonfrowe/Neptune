subroutine plotdata(npt,x,y,yerr,bb)
use precision
implicit none
integer :: npt,i
real :: diff
real, allocatable, dimension(:) :: xp,yp,yperr
real, dimension(:) :: bb
real(double), dimension(:) :: x,y,yerr

!allocate space for real variables that are used by PGPLOT
allocate(xp(npt),yp(npt),yperr(npt))

!convert from dble to real
xp=real(x(1:npt))
yp=real(y(1:npt))
yperr=real(yerr(1:npt))

!calculate size of plotting window
if(bb(1).eq.bb(2))then
   diff=maxval(xp)-minval(xp)
   bb(1)=minval(xp)-0.1*diff
   bb(2)=maxval(xp)+0.1*diff
endif
if(bb(3).eq.bb(4))then
   diff=maxval(yp)-minval(yp)
   bb(3)=minval(yp)-0.1*diff
   bb(4)=maxval(yp)+0.1*diff
endif

!plotting commands
call pgsci(1)
call pgvport(0.15,0.95,0.15,0.95) !make room around the edges for labels
call pgwindow(bb(1),bb(2),bb(3),bb(4)) !plot scale
call pgbox("BCNTS1",0.0,0,"BCNTS",0.0,0)
call pglabel("X","Y","")  !labels
call pgline(npt,xp,yp)    !plot a line
call pgpt(npt,xp,yp,17)   !plot points
call PGERRB(6,npt,xp,yp,yperr,1.0) !plot errors

return
end subroutine plotdata
