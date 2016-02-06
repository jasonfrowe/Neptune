program fftpow
!uses spline resampling and FFTs to do quick power spectrum analysis.
use precision
implicit none
integer :: nmax,npt,iargc
real :: tstart,tfinish
real, allocatable, dimension(:) :: bb
real(double) :: minx,mean
real(double), allocatable, dimension(:) :: time,flux,ferr
character(80) :: filename

interface
   subroutine readfftdata(filename,nmax,npt,time,flux,ferr,minx,mean)
      use precision
      implicit none
      integer :: nmax,npt
      real(double) :: minx,mean
      real(double), dimension(:) :: time,flux,ferr
      character(80) :: filename
   end subroutine readfftdata
end interface
interface !makes a plot of your data.
   subroutine plotdatascatter(npt,x,y,yerr,bb)
      use precision
      implicit none
      integer, intent(inout) :: npt
      real, dimension(:), intent(inout) :: bb
      real(double), dimension(:), intent(inout) :: x,y,yerr
   end subroutine plotdatascatter
end interface
interface
   subroutine fftspec(npt,time,flux,ferr)
      use precision
      implicit none
      integer, intent(inout) :: npt
      real(double), dimension(:), intent(inout) :: time,flux,ferr
   end subroutine fftspec
end interface

CALL CPU_TIME(tstart) !for timing runtimes

!check that we have enough information from the commandline
if(iargc().lt.1)then !if not, spit out the Usage info and stop.
   write(0,*) "Usage: pixelfit filename"
   stop
endif

!read in filename containing data (3 columns)
call getarg(1,filename)

nmax=80000 !maximum number of data points
nmax=1400000
allocate(time(nmax),flux(nmax),ferr(nmax))

call readfftdata(filename,nmax,npt,time,flux,ferr,minx,mean)
write(0,*) "Number of points read: ",npt !report how much data was read in

!open PGPLOT device
call pgopen('/xserve')  !'?' lets the user choose the device.
call PGPAP (8.0 ,1.0) !use a square 8" across
call pgsubp(1,4)
call pgpage() !create a fresh page
call pgslw(3) !thicker lines
call pgsch(2.9) !bigger text

allocate(bb(4))
bb=0.0
call plotdatascatter(npt,time,flux,ferr,bb)

call fftspec(npt,time,flux,ferr)

call pgclos()

CALL CPU_TIME(tfinish)

write(0,*) tstart,tfinish,tfinish-tstart

end program fftpow
