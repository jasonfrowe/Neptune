subroutine resample(npt,time,flux,ferr,ns,trs,frs,iresampletype,seed,dt,&
   mintime,gap,debug)
use precision
implicit none
integer :: iresampletype,ns,seed,npt,debug
real(double) :: dt,mintime,gap
real(double), dimension(npt) :: time,flux,ferr
real(double), dimension(ns) :: trs,frs
!local vars
integer :: i1,i,loop
real(double) :: ts,fs,dtime,gasdev

interface !creates a co-variance matrix
   subroutine makekernel(Kernel,npt1,npt2,x1,x2,npt,yerr,npars,pars)
      use precision
      implicit none
      integer, intent(inout) :: npt1,npt2,npars,npt
      real(double), dimension(:), intent(inout) :: x1,x2,yerr,pars
      real(double), dimension(:,:), intent(inout) :: Kernel
   end subroutine makekernel
end interface

if(iresampletype.eq.1)then

   !calculate interpolated values using simple linear model
   i1=2
   loop=0
   do i=1,ns
      ts=mintime+dt*dble(i-1)

      do while (loop.eq.0)
         if(time(i1).ge.ts)then
            loop=1
            dtime=time(i1)-time(i1-1)
            if((dtime.gt.gap*dt).and.(gap.gt.0.0d0))then
               fs=gasdev(seed)*(ferr(i1-1)+ferr(i1))/2.0d0
            else
               fs=flux(i1-1)+(flux(i1)-flux(i1-1))*(ts-time(i1-1))/dtime
            endif
         else
            i1=i1+1
         endif
      enddo

      trs(i)=ts
      frs(i)=fs
      loop=0
   enddo

endif

!if debug=1 then write out the interpolated lightcurve for inspection.
if(debug.eq.1)then
   open(unit=11,file="interpolate.dat")
   do i=1,ns
      write(11,*) trs(i),frs(i)
   enddo
   close(11)
endif

return
end subroutine resample

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine getdtns(npt,time,nover,dt,ns,nfft,mintime,maxtime)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!get sampling interval and calculate number of new sample points
use precision
implicit none
integer :: npt,ns,nfft,nover
real(double) :: dt,mintime,maxtime
real(double), dimension(npt) :: time
!local vars
integer :: ndt,i
integer, allocatable, dimension(:) :: p
real(double), allocatable, dimension(:) :: dts

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
!write(0,*) "dt: ",dt
deallocate(dts,p)

!range for resampled times.
mintime=minval(time(1:npt))
maxtime=maxval(time(1:npt))

ns=(maxtime-mintime)/dt !number of resampled data points
!write(0,*) "npt,ns: ",npt,ns
!get an estimate of array size for FFTW that power is a 2
nfft=2**int(log10(dble(npt*nover))/log10(2.0d0)+1.0d0)

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function sinc(x)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! returns sin(pi x) /(pi x)
use precision
implicit none
real(double), parameter :: pi = 3.1415926535897932_8
real(double) :: sinc,x

if(x.eq.0)then
   sinc=1.0d0
else
   sinc=sin(pi*x)/(pi*x)
endif

return
end
