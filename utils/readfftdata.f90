subroutine readfftdata(filename,nmax,npt,time,flux,ferr,minx,mean)
use precision
implicit none
!input vars
integer :: nmax,npt
real(double) :: minx,mean
real(double), dimension(:) :: time,flux,ferr
character(80) :: filename
!local vars
integer :: nunit,filestatus,i
real(double) :: t,f,fe

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
   read(nunit,*,iostat=filestatus) t,f,fe
   if(filestatus == 0) then
      i=i+1
      time(i)=t
      flux(i)=f
      ferr(i)=fe
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
minx=minval(time(1:npt))
!time(1:npt)=time(1:npt)-minx
mean=Sum(flux(1:npt))/dble(npt)
!flux(1:npt)=-(flux(1:npt)-mean)
flux(1:npt)=flux(1:npt)-mean
!flux(1:npt)=flux(1:npt)/mean-1.0
!ferr(1:npt)=ferr(1:npt)/mean

return
end subroutine
