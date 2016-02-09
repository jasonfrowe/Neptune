subroutine resample(npt,time,flux,ferr,ns,trs,frs,iresampletype,seed,dt,&
   mintime,gap,debug)
use precision
implicit none
integer :: iresampletype,ns,seed,npt,debug
real(double) :: dt,mintime,gap
real(double), dimension(npt) :: time,flux,ferr
real(double), dimension(ns) :: trs,frs
!local vars
integer :: i1,i2,i,j,k,loop,nbuf,nstep,npars,nss,nbufstep,npts,info,    &
   nrhs
real(double) :: ts,fs,dtime,gasdev,dtmax
real(double), allocatable, dimension(:) :: pars,dts,alpha,yerr2,mu
real(double), allocatable, dimension(:,:) :: Kernel,KernelZ
character(80) :: line

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
               fs=flux(i1-1)+(flux(i1)-flux(i1-1))*(ts-time(i1-1))/dtime
!               fs=fs+gasdev(seed)*(ferr(i1-1)+ferr(i1))/2.0d0
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

elseif(iresampletype.eq.2)then

   write(0,*) "Sorry SINC resampling is currently broken"
   stop

elseif(iresampletype.eq.3)then

   nbuf=1000 !break up resampling to managable size
   nstep=(ns/nbuf)+1 !number of steps to cover sample size

   !Here are the parameters that control the co-variance matrix and fitted
   !parameters
   npars=4 !number of parameters used for model of the matrix
   allocate(pars(npars))  !model parameters for Kernel generation.

   do k=1,nstep
      !number of data points in buffer
      nss=nbuf
      if(k.eq.nstep)then !last step will have less points
         nss=ns-nbuf*(nstep-1)
      endif

      !times for resampling
      nbufstep=nbuf*(k-1) !location in resampled array
      do i=1,nss
         trs(i+nbufstep)=mintime+dt*dble(i-1+nbufstep)
      enddo

      !find range of unsampled data
      i1=1
      loop=0
      do while(loop.eq.0)
         if(time(i1).ge.trs(nbufstep+1))then
            i1=i1-1
            loop=1
         elseif(i1.gt.ns)then
            write(0,*) "Gaussian process interpolator failed.."
            write(0,*) "i1 > ns"
            stop
         else
            i1=i1+1
         endif
      enddo

      i2=1
      loop=0
      do while(loop.eq.0)
         if(time(i2).ge.trs(nss+nbufstep))then
            loop=1
         else
            i2=i2+1
         endif
      enddo


      i1=max(1,i1)
      i2=min(ns,i2)

      npts=i2-i1+1 !how much of the original sample are we using.

!      write(0,*) "iX: ",i1,i2
!      write(0,*) "nss,npts: ",nss,npts
!      read(5,*)

      allocate(dts(npts-1))
      j=0
      do i=i1+1,i2
         j=j+1
         dts(j)=time(i)-time(i-1)
      enddo
      dtmax=maxval(dts(1:j))

      pars(1)=1.0d0 !amp scale for exp
      pars(2)=dt !length scale for exp
      pars(3)=dtmax/dt !second amp scale
      pars(4)=dtmax !second length scale
      deallocate(dts)

      !make Kernel
      !lets make a Kernel/co-variance for the Gaussian process
      allocate(Kernel(npts,npts)) !allocate space
      call makekernel(Kernel,npts,npts,time(i1:i2),time(i1:i2),npt,     &
         ferr(i1:i2),npars,pars) !create Kernel
!      write(0,*) "Kernel made.."

      !Cholesky factorization
!      write(0,*) "Beginning Cholesky factorization"
      call dpotrf('U',npts,Kernel,npts,info) !LAPACK routine for Cholesky
      if (info.ne.0) then !check for errors
         write(0,*) "Cholesky factorization failed"
         write(0,*) "dpotrf info: ",info
         stop
      endif

!      write(0,*) "Beginning Cholesky Solver"
      allocate(alpha(npts))
      alpha(1:npts)=flux(i1:i2) !dpotrs takes alpha as input and output.
      nrhs=1 !how man columns does alpha contain - just one
      call dpotrs('U',npts,nrhs,Kernel,npts,alpha,npts,info) !call LAPACK
      if (info.ne.0) then !check for errors
         write(0,*) "Solver failed.."
         write(0,*) "dpotrs info: ",info
         stop
      endif

      !predict samples
      allocate(KernelZ(nss,npts))
      allocate(yerr2(npts),mu(nss))
      yerr2=0.0d0
!      call makekernel(KernelZ,nss,npts,trs(1+nbufstep:nss+nbufstep),    &
!         time(i1:i2),npts,yerr2,npars,pars)
      call makekernel(KernelZ,nss,npts,trs(1+nbufstep:nss+nbufstep),    &
         time(i1:i2),npts,ferr(i1:i2),npars,pars)
      mu=matmul(KernelZ,alpha)

!     fill in big gaps with noise.
      j=max(i1,2)
      do i=1,nss
         ts=trs(i+nbufstep) !time at sample point
         loop=0
         do while (loop.eq.0)
            if(time(j).ge.ts)then
               loop=1
               dtime=time(j)-time(j-1)
!               write(0,*) dtime,gap*ddt
               if((dtime.gt.gap*dt).and.(gap.gt.0.0d0))then
!                  write(0,*) "gap.."
!                  mu(i)=gasdev(seed)*(ferr(j-1)+ferr(j))/2.0d0
                  mu(i)=flux(j-1)+(flux(j)-flux(j-1))*(ts-time(j-1))/dtime
!                  mu(i)=mu(i)+gasdev(seed)*(ferr(j-1)+ferr(j))/2.0d0
               endif
            else
               j=j+1
            endif
         enddo
      enddo

      frs(1+nbufstep:nss+nbufstep)=mu(1:nss)

      deallocate(Kernel,alpha,KernelZ,yerr2,mu)

      write(line,'(A24,I5,A1,I5)') "Gaussian predictor done ",k,"/",nstep
      call ovrwrt(line,2)

   enddo

   deallocate(pars)

   write(0,*) "Gaussian Predictor Finished                   "

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
nfft=2**int(log10(dble(ns*nover))/log10(2.0d0)+1.0d0)

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
