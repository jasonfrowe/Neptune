subroutine fftspec(npt,time,flux)
use precision
implicit none
!import vars
integer :: npt
real(double), dimension(:) :: time,flux
!local vars
integer :: ndt,i,ns
integer, allocatable, dimension(:) :: p
real(double) :: dt,mintime,maxtime,ddt,ts,fs
real(double), allocatable, dimension(:) :: flux2,dts,trs,frs

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
allocate(trs(ns),frs(ns)) !allocate space for resampled

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



return
end subroutine fftspec
