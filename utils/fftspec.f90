subroutine fftspec(npt,time,flux)
use precision
!iso_c_binding is for FFTW3 interface
use, intrinsic :: iso_c_binding
implicit none
!add in the FFTW3 modules
include 'fftw3.f03'
!import vars
integer :: npt,nfft
real(double), dimension(:) :: time,flux
!local vars
integer :: ndt,i,ns,nover
integer, allocatable, dimension(:) :: p
real(double) :: dt,mintime,maxtime,ddt,ts,fs,cd2uhz
real(double), allocatable, dimension(:) :: flux2,dts,trs,frs
!FFTW3 vars
type(C_PTR) :: plan
integer ( kind = 4 ) :: nh
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: frsC

!oversampling for FFT
nover=8

!convert from c/d to uHz
cd2uhz=1.0d6/86400.0d0

!resample data
ndt=npt-1
allocate(dts(ndt))
!$OMP PARALLEL DO
do i=1,ndt
   dts(i)=time(i+1)-time(i) !calculate delta-times
enddo
!$OMP END PARALLEL DO
!sort dts to get median dt
allocate(p(ndt))
call rqsort(ndt,dts,p)
dt=dts(p(ndt/2))
ddt=dble(dt) !precompute int -> double
write(0,*) "dt: ",dt
deallocate(dts,p)

!set up spline
allocate(flux2(npt))
call spline(time,flux,npt,1.d30,1.d30,flux2)

mintime=minval(time(1:npt))
maxtime=maxval(time(1:npt))

ns=(maxtime-mintime)/ddt
nfft=2**int(log10(dble(npt*nover))/log10(2.0d0)+1.0d0)
write(0,*) "nfft: ",nfft
allocate(trs(nfft),frs(nfft)) !allocate space for resampled
frs=0.0d0

!calculate interpolated spline values
!$OMP PARALLEL DO PRIVATE(ts,fs)
do i=1,ns
   ts=mintime+ddt*dble(i-1)
   call splint(time,flux,flux2,npt,ts,fs)
   trs(i)=ts
   frs(i)=fs
enddo
!$OMP END PARALLEL DO
deallocate(flux2) !spline is done, deallocate derivative

nh=(nfft/2)+1
allocate(frsC(nh))

!generate FFTW plan
plan=fftw_plan_dft_r2c_1d(nfft,frs,frsC,FFTW_ESTIMATE)
!execute FFW
call fftw_execute_dft_r2c(plan,frs,frsC)
!destroy plans
call fftw_destroy_plan(plan)

open(unit=11,file="fft.dat")
do i=1,nh
!   write(0,*) i
   write(11,*) cd2uhz*dble(i-1)/(dt*dble(nfft)),abs(frsC(i))/dble(nfft)
enddo
close(11)

return
end subroutine fftspec
