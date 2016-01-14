subroutine displaykernel(nxmax,nymax,parray,bpix)
!Jason Rowe 2015 - jasonfrowe@gmail.com
use precision
implicit none
integer :: nxmax,nymax,npt,i,j,ncol,dumi,ndiff,k
integer, dimension(4) :: nr
integer, allocatable, dimension(:) :: p
integer, allocatable, dimension(:,:) :: ia
real :: r,g,b,dumr,x2,y2,xr
real, dimension(4) :: rj
real(double) :: bpix,maxp,minp,z1,z2,med,std,stdev2,lmin,lmax, &
   ratio,lmed,lstd,minlp
real(double), dimension(:,:) :: parray
real(double), allocatable, dimension(:) :: a
real(double), allocatable, dimension(:,:) :: lparray

allocate(lparray(nxmax,nymax)) !used for making a log-scale plot

ncol=64 !number of colours for display

nr=1 !range of pixels for plotting
npt=0
do i=1,nxmax
   do j=1,nymax
      if(parray(i,j).lt.bpix)then
         npt=npt+1
         if(npt.eq.1)then
            maxp=parray(i,j)
            minp=parray(i,j)
         else
            maxp=max(maxp,parray(i,j))
            minp=min(minp,parray(i,j))
         endif
      endif
   enddo
enddo
!write(0,*) "Min/Max:",minp,maxp
nr(1)=1
nr(3)=1
nr(2)=nxmax
nr(4)=nymax

if(npt.le.1)then
   return
endif

rj=real(nr)

!write(0,*) "hello",nr

allocate(a(npt),p(npt))
k=0
do i=nr(1),nr(2)
   do j=nr(3),nr(4)
      if(parray(i,j).lt.bpix)then
         k=k+1
         a(k)=parray(i,j)
!         write(0,*) k,a(k)
         if(k.ge.1000) goto 10
      endif
   enddo
enddo
10 continue
!write(0,*) "K:",k
if(k.ge.3)then
   call rqsort(k,a,p)
   med=a(p(k/2))
   std=stdev2(k,a,med)
else
   med=0.0
   std=1.0
endif
deallocate(a,p)
!write(0,*) "Med,std: ",med,std

allocate(ia(nxmax,nymax))

!write(6,*) "yo..",nr
minlp=minp
lmin=1000.0
lmax=-1000.0
do i=nr(1),nr(2)
   do j=nr(3),nr(4)
      if(parray(i,j).lt.bpix)then
         if(parray(i,j)-minlp+1.0.le.0.0d0)then
            lparray(i,j)=0.0
         else
            lparray(i,j)=log10(parray(i,j)-minlp+1.0)
         endif
         lmin=min(lparray(i,j),lmin)
         lmax=max(lparray(i,j),lmax)
      endif
   enddo
enddo
!write(0,*) "lmin,lmax: ",lmin,lmax
!z1=log10(minlp-minlp+1.0)
z1=log10(max(0.0,med-2.0*std-minlp)+1.0)
z2=log10(maxp-minlp+1.0)

!uncomment following 3 lines for a sqrt scale opposed to log
!lparray=sqrt(parray-minp)
!z1=sqrt(med-std-minp)
!z2=sqrt(maxp-minp)

!linear display
lparray(nr(1):nr(2),nr(3):nr(4))=parray(nr(1):nr(2),nr(3):nr(4))
z1=minval(lparray(nr(1):nr(2),nr(3):nr(4)))
z2=maxval(lparray(nr(1):nr(2),nr(3):nr(4)))

do i=nr(1),nr(2)
   do j=nr(3),nr(4)
      if(parray(i,j).lt.bpix)then
         IA(i,j)=int((lparray(i,j)-z1)/(z2-z1)*dble(NCOL-1))+16
         if(lparray(i,j).le.z1) then
            ia(i,j)=16
         endif
      else
         ia(i,j)=15
      endif
      if(lparray(i,j).gt.z2) ia(i,j)=ncol+15
   enddo
enddo

!set up pgplot window
call pgscr(15,0.0,0.3,0.2)

call pgsch(1.5) !make the font a bit bigger
call pgslw(3)  !make the lines a bit thicker

call pgsci(1)
call pgvport(0.0,1.00,0.0,1.0)
call pgwindow(0.0,1.0,0.0,1.0)

call pgvport(0.15,0.95,0.15,0.95) !make room around the edges for labels
call pgsci(1)
call pgwindow(rj(1),rj(2),rj(3),rj(4)) !plot scale
call pgbox("BCNTS1",0.0,0,"BCNTS1",0.0,0)
call pglabel("X1","X2","")
call pgsci(1)

do i=1,ncol
   call heatlut(i*4-3,r,g,b)
   call heatlut(i*4-2,r,g,b)
   call heatlut(i*4-1,r,g,b)
   call heatlut(i*4  ,r,g,b)
   CALL PGSCR(I+15, R, G, B)
enddo

xr=real(max(nr(2)-nr(1),nr(4)-nr(3)))
x2=real(nr(2)-nr(1))/xr
y2=real(nr(4)-nr(3))/xr

call pgpixl(ia,nxmax,nymax,nr(1),nr(2),nr(3),nr(4),rj(1),rj(2),rj(3),rj(4))
call pgsci(1)
call pgbox("BCNTS1",0.0,0,"BCNTS1",0.0,0)
call pgsch(1.0) !reset font size

deallocate(ia,lparray)

return
end subroutine displaykernel
