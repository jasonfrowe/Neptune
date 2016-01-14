subroutine getdata(filename,npt,nmax,x,y,yerr,xpos,ypos,xnep,ynep)
use precision
implicit none
integer :: npt,nmax,nunit,filestatus,i,npix
real(double) :: minx,mean,time,flux,dumr,sky,xp,yp,xn,yn
real(double), dimension(:) :: x,y,yerr,xpos,ypos,xnep,ynep
character(80) :: filename,dumc

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
      if((xp.gt.-0.2).and.(xp.lt.0.4).and.(yp.gt.-0.15).and.      &
       (yp.lt.0.25))then
         i=i+1
         x(i)=time
         y(i)=flux-sky*dble(npix)
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
