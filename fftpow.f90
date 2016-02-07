program fftpow
!uses spline resampling and FFTs to do quick power spectrum analysis.
use precision
implicit none
integer :: nmax,npt,iargc,iresampletype,seed,nover,ns,nfft,debug,nh
integer, dimension(3) :: now
real :: tstart,tfinish
real, allocatable, dimension(:) :: bb
real(double) :: minx,mean,ran2,dumr,gap,dt,mintime,maxtime
real(double), allocatable, dimension(:) :: time,flux,ferr,t1,t2,t3,trs, &
   frs,amp
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
   subroutine fftspec(nfft,frs,amp,npt,dt,debug)
      use precision
      implicit none
      integer :: nfft,debug,npt
      real(double) :: dt
      real(double), dimension(:) :: frs,amp
   end subroutine fftspec
end interface

CALL CPU_TIME(tstart) !for timing runtimes

!parameters controling resampling and FFTs
nover=10 !oversampling for FFT
gap=0.0d0 !identifying gaps in data and replace with white noise.
!resampling routine
! 1-linear interpolation
! 2-sinc interpolation (not working)
! 3-Gaussian-process predictive Mean
iresampletype=1

!check that we have enough information from the commandline
if(iargc().lt.1)then !if not, spit out the Usage info and stop.
   write(0,*) "Usage: pixelfit filename"
   stop
endif

!read in filename containing data (3 columns)
call getarg(1,filename)

nmax=1400000 !maximum number of data points
allocate(time(nmax),flux(nmax),ferr(nmax))

call readfftdata(filename,nmax,npt,time,flux,ferr,minx,mean)
write(0,*) "Number of points read: ",npt !report how much data was read in
!compact data arrays
allocate(t1(npt),t2(npt),t3(npt))
t1=time(1:npt)
t2=flux(1:npt)
t3=ferr(1:npt)
deallocate(time,flux,ferr)
allocate(time(npt),flux(npt),ferr(npt))
time=t1
flux=t2
ferr=t3
deallocate(t1,t2,t3)

!We have data, so we can initialize the random-number generator
!-only used for special cases when resampling with large gaps
call itime(now)
seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
dumr=ran2(-seed)

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

!get stats about datasizes
call getdtns(npt,time,nover,dt,ns,nfft,mintime,maxtime)
!write(0,*) "dt, ns, nfft"
!write(0,*) dt,ns,nfft

debug=0 !if =1, then resampled lightcurve is writen to "interpolate.dat"
allocate(trs(ns),frs(nfft)) !allocate space for resampled
frs=0.0d0 !needs to be initiated to zero for zero-padding for oversampling
call resample(npt,time,flux,ferr,ns,trs,frs,iresampletype,seed,dt,      &
   mintime,gap,debug)

!number of frequency/amplitudes from fftspec to be returned
nh=(nfft/2)+1
allocate(amp(nh))

!calculate amplitude spectrum
debug=1
call fftspec(nfft,frs,amp,npt,dt,debug)

!plot Power-spectrum
call pgpage()
bb=0.0 !auto-scale the plot
call plotspec(nh,nfft,amp,dt,bb)

call pgclos()

CALL CPU_TIME(tfinish)

write(0,*) tstart,tfinish,tfinish-tstart

end program fftpow
