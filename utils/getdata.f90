subroutine getdata(filename,npt,nmax,x,y,yerr)
use precision
implicit none
integer :: npt,nmax,nunit,filestatus,i
real(double) :: minx,mean
real(double), dimension(:) :: x,y,yerr
character(80) :: filename

nunit=10
open(unit=nunit,file=filename,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",filename
   stop
endif

i=1
do
   if(i.gt.nmax)then
      write(0,*) "Increase nmax to match data points"
      write(0,*) "nmax: ",nmax
      stop
   endif
   read(nunit,*,iostat=filestatus) x(i),y(i),yerr(i)
   if(filestatus == 0) then
      i=i+1
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

minx=minval(x(1:npt))
x(1:npt)=x(1:npt)-minx
mean=Sum(y(1:npt))/dble(npt)
y(1:npt)=y(1:npt)-mean
!yerr(1:npt)=yerr(1:npt)/mean

return
end subroutine getdata
