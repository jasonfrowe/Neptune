subroutine getdata(filename,npt,nmax,x,y,yerr,xpos,ypos,xnep,ynep,minx,mean)
use precision
implicit none
integer :: npt,nmax,nunit,filestatus,i,npix
real(double) :: time,flux,dumr,sky,xp,yp,xn,yn,minx,mean,cencut(4)
real(double), dimension(:) :: x,y,yerr,xpos,ypos,xnep,ynep
character(80) :: filename,dumc

!Centroid Cuts for Neptune
cencut(1)=-0.2
cencut(2)= 0.4
cencut(3)=-0.15
cencut(4)= 0.25 

!Centroid Cuts for Uranus
cencut(1)=-0.2
cencut(2)= 0.3
cencut(3)=-0.15
cencut(4)= 0.25

nunit=10
open(unit=nunit,file=filename,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",filename
   stop
endif

i=0
do
   if(i.gt.nmax)then
      write(0,*) "Increase nmax to match data points"
      write(0,*) "nmax: ",nmax
      stop
   endif
   read(nunit,*,iostat=filestatus) time,flux,sky,dumr,xn,yn,npix,   &
      dumc,xp,yp
   if(filestatus == 0) then
      if((xp.gt.cencut(1)).and.(xp.lt.cencut(2)).and.(yp.gt.cencut(3)).and.      &
       (yp.lt.cencut(4)))then
         i=i+1
         x(i)=time
         y(i)=flux-sky*dble(npix) !sky correction
         yerr(i)=sqrt(abs(flux)) !approx error
         xpos(i)=xp
         ypos(i)=yp
         xnep(i)=xn
         ynep(i)=yn
      endif
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i+1
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit) !close file
npt=i-1

!scale and shift data.
minx=minval(x(1:npt))
x(1:npt)=x(1:npt)-minx
mean=Sum(y(1:npt))/dble(npt)
y(1:npt)=y(1:npt)/mean-1.0
yerr(1:npt)=yerr(1:npt)/mean

return
end subroutine getdata

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine getpdata(filename,npt,nmax,x,y,yerr,oflux,pmod,smod,xpos,    &
 ypos,xnep,ynep,minx,mean)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!Routine to readin data for pixelfit
use precision
implicit none
integer :: nunit,filestatus,i,nmax,npt
real(double) :: minx,mean,time,flux,ferr,origflux,pmodel,spmodel,xn,yn, &
 xp,yp
real(double), dimension(:) :: x,y,yerr,oflux,pmod,smod,xpos,ypos,xnep,  &
   ynep
character(80) :: filename

nunit=10
open(unit=nunit,file=filename,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",filename
   stop
endif

i=0
do
   if(i.gt.nmax)then
      write(0,*) "Increase nmax to match data points"
      write(0,*) "nmax: ",nmax
      stop
   endif
   read(nunit,*,iostat=filestatus) time,flux,ferr,origflux,pmodel,      &
      spmodel,xn,yn,xp,yp
   if(filestatus == 0) then
      i=i+1
      x(i)=time
      y(i)=flux
      yerr(i)=ferr
      oflux(i)=origflux
      pmod(i)=pmodel
      smod(i)=spmodel
      xpos(i)=xn
      ypos(i)=yn
      xnep(i)=xp
      ynep(i)=yp
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i+1
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit) !close file
npt=i-1

!scale and shift data.
minx=minval(x(1:npt))
x(1:npt)=x(1:npt)-minx
mean=Sum(y(1:npt))/dble(npt)
y(1:npt)=y(1:npt)/mean-1.0
yerr(1:npt)=yerr(1:npt)/mean

return
end subroutine getpdata
