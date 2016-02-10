program gendata
use precision
implicit none
integer nmax,npt,iargc,seed,filestatus,nunit,i,nf,j
integer, dimension(3) :: now
real(double), parameter :: pi = 3.1415926535897932_8
real(double) :: ran2,gasdev,dumr,t,flux,dt
real(double), allocatable, dimension(:) :: time,amp,phase,freq
character(80) :: filename

!check that we have enough information from the commandline
if(iargc().lt.1)then !if not, spit out the Usage info and stop.
   write(0,*) "Usage: gendata filename"
   stop
endif
!read in filename containing data (3 columns)
call getarg(1,filename)

!set up seed for random number generator
call itime(now)
seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
dumr=ran2(-seed)

nmax=100000 !maximum number of data points

!read in time stamps
allocate(time(nmax))
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
   read(nunit,*,iostat=filestatus) t
   if(filestatus == 0) then
      i=i+1
      time(i)=t
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
write(0,*) "Number of points read: ",npt !report how much data was read in

nf=4 !number of frequencies to generate
allocate(amp(nf),phase(nf),freq(nf))

do i=1,nf
!   freq(i)=10.0**dble(i-1)
   freq(i)=0.1*10**dble(i-1)
   amp(i)=1.0!/dble(i)
   phase(i)=2.0d0*Pi*ran2(seed) !random phase
   write(0,*) freq(i),amp(i),phase(i)
enddo

dt=1.362160997800288E-003
do i=1,npt
!   time(i)=dt*dble(i-1)
!   flux=1.0*gasdev(seed)
   flux=0.0d0
   do j=1,nf
      flux=flux+amp(j)*sin(2.0d0*pi*freq(j)*time(i)+phase(j))
   enddo
   write(6,*) time(i),flux,0.1
enddo

end program gendata
