subroutine plotspec(nh,nfft,amp,dt,bb)
use precision
implicit none
integer :: nh,nfft,nplot
real, dimension(4) :: bb
real(double) :: dt
real(double), dimension(nh) :: amp
!local vars
integer :: i
real :: minx,maxx,miny,maxy,diff,dumr
real, allocatable, dimension(:) :: x,y
real(double) :: f,p,cd2uhz

!convert from c/d to uHz
cd2uhz=1.0d6/86400.0d0

nplot=nh-1
!PGPlot only accepts real
allocate(x(nplot),y(nplot))

do i=2,nh

   f=cd2uhz*dble(i-1)/(dt*dble(nfft))
   p=1.0d6*amp(i-1)

   x(i-1)=log10(real(f))
   y(i-1)=real(p)

enddo

!compute range of data
minx=minval(x)
maxx=maxval(x)
miny=minval(y)
maxy=maxval(y)

!calculate size of plotting window
if(bb(1).eq.bb(2))then
   diff=maxx-minx
   bb(1)=minx-0.0*diff
   bb(2)=maxx+0.0*diff
else
   !recalculate bounds for Y direction
   dumr=maxy
   maxy=miny
   miny=dumr
   do i=1,nplot
      if((x(i).gt.bb(1)).and.(x(i).lt.bb(2)))then
         if(y(i).lt.miny) miny=y(i)
         if(y(i).gt.maxy) maxy=y(i)
      endif
   enddo
endif
if(bb(3).eq.bb(4))then
   diff=maxy-miny
   bb(3)=miny-0.0*diff
   bb(4)=maxy+0.1*diff
endif

!write(0,*) bb(1),bb(2),bb(3),bb(4)

!plotting commands
call pgsci(1)
call pgvport(0.15,0.95,0.20,0.95) !make room around the edges for labels
call pgwindow(bb(1),bb(2),bb(3),bb(4)) !plot scale
call pgbox("BLCNTS",0.0,0,"BCNTVS",0.0,0)
call pgptxt((bb(1)+bb(2))/2.0,bb(3)-0.20*(bb(4)-bb(3)),0.0,0.5,         &
   "Frequency (\(0638)Hz)")
call pgptxt(bb(1)-0.10*(bb(2)-bb(1)),(bb(4)+bb(3))/2,90.0,0.5,          &
   "Amplitude (ppm)")
call pgline(nplot,x,y)   !plot line

end subroutine plotspec

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine plotpspec(nh,nfft,amp,dt,bb)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!plot power-spectrum
use precision
implicit none
integer :: nh,nfft,nplot
real, dimension(4) :: bb
real(double) :: dt
real(double), dimension(nh) :: amp
!local vars
integer :: i
real :: minx,maxx,miny,maxy,diff,dumr,rcd2uhz
real, allocatable, dimension(:) :: x,y,rbb
real(double) :: f,p,cd2uhz

!convert from c/d to uHz
cd2uhz=1.0d6/86400.0d0
rcd2uhz=real(cd2uhz)

nplot=nh-1
!PGPlot only accepts real
allocate(x(nplot),y(nplot))

do i=2,nh

   f=cd2uhz*dble(i-1)/(dt*dble(nfft))
   p=1.0d12*amp(i)*amp(i)/f

   x(i-1)=log10(real(f))
   y(i-1)=log10(real(p))

enddo

!compute range of data
minx=minval(x)
maxx=maxval(x)
miny=minval(y)
maxy=maxval(y)

!calculate size of plotting window
if(bb(1).eq.bb(2))then
   diff=maxx-minx
   bb(1)=minx-0.0*diff
   bb(2)=maxx+0.0*diff
else
   !recalculate bounds for Y direction
   dumr=maxy
   maxy=miny
   miny=dumr
   do i=1,nplot
      if((x(i).gt.bb(1)).and.(x(i).lt.bb(2)))then
         if(y(i).lt.miny) miny=y(i)
         if(y(i).gt.maxy) maxy=y(i)
      endif
   enddo
endif
if(bb(3).eq.bb(4))then
   diff=maxy-miny
   bb(3)=miny-0.0*diff
   bb(4)=maxy+0.1*diff
endif

!write(0,*) bb(1),bb(2),bb(3),bb(4)

!plotting commands
call pgsci(1)
call pgvport(0.15,0.95,0.20,0.95) !make room around the edges for labels
call pgwindow(bb(1),bb(2),bb(3),bb(4)) !plot scale
call pgbox("BLNTS",0.0,0,"BCLNTVS",0.0,0)
call pgptxt((bb(1)+bb(2))/2.0,bb(3)-0.20*(bb(4)-bb(3)),0.0,0.5,         &
   "Frequency (\(0638)Hz)")
call pgptxt(bb(1)-0.10*(bb(2)-bb(1)),(bb(4)+bb(3))/2,90.0,0.5,          &
   "PD (ppm\u2\d \(0638)Hz\u-1\d)")

call pgline(nplot,x,y)   !plot line

allocate(rbb(4))
rbb(1)=log10(10.0**bb(1)/rcd2uhz)
rbb(2)=log10(10.0**bb(2)/rcd2uhz)
rbb(3)=bb(3)
rbb(4)=bb(4)
call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4)) !plot scale
call pgbox("CLMTS",0.0,0,"",0.0,0)
call pgptxt((rbb(1)+rbb(2))/2.0,rbb(4)+0.15*(rbb(4)-rbb(3)),0.0,0.5,         &
   "Frequency (c/d)")

call pgwindow(bb(1),bb(2),bb(3),bb(4)) !restore plot scale

end subroutine plotpspec
