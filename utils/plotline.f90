subroutine plotline(bb,ans,eans)
use precision
implicit none
integer nplot,i
real, dimension(:) :: bb
real, allocatable, dimension(:) :: px,py,py2,rans,reans,val
real(double), dimension(:) :: ans,eans

allocate(rans(2),reans(2)) !convert dble to real
rans=real(ans)
reans=real(eans)

nplot=2 !sampling for plotting the line
allocate(px(nplot),py(nplot))
px(1)=bb(1)
px(2)=bb(2)
py(1)=rans(1)+px(1)*rans(2)
py(2)=rans(1)+px(2)*rans(2)
call pgsci(5) !change colour
call pgline(nplot,px,py)
call pgsci(1)
deallocate(px,py)

nplot=1000 !sampling for plotting uncertainty
allocate(px(nplot),py(nplot),py2(nplot),val(4))

do i=1,nplot
   px(i)=bb(1)+(bb(2)-bb(1))*real(i-1)/real(nplot-1)
   val(1)=rans(1)+reans(1)+(rans(2)+reans(2))*px(i)
   val(2)=rans(1)-reans(1)+(rans(2)+reans(2))*px(i)
   val(3)=rans(1)+reans(1)+(rans(2)-reans(2))*px(i)
   val(4)=rans(1)-reans(1)+(rans(2)-reans(2))*px(i)
   py(i)=maxval(val)
   py2(i)=minval(val)
enddo
call pgsci(4)
call pgline(nplot,px,py)
call pgline(nplot,px,py2)
call pgsci(1)
deallocate(px,py,py2)


end subroutine plotline
