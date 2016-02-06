subroutine fftspec(npt,time,flux,ferr)
use precision
!iso_c_binding is for FFTW3 interface
use, intrinsic :: iso_c_binding
implicit none
!add in the FFTW3 modules
include 'fftw3.f03'
!import vars
integer :: npt,nfft,iresampletype
real(double), dimension(:) :: time,flux,ferr
!local vars
integer :: ndt,i,ns,nover,i1,loop,seed,j,nB,nbuf,nstep,k,i2,nbufstep,   &
   npts,nss,nrhs,rank,lwork,NLVL,SMLSIZ,info,ii,npars
integer, dimension(3) :: now
integer, allocatable, dimension(:) :: p,iwork
real(double) :: dt,mintime,maxtime,ddt,ts,fs,cd2uhz,dtime,ran2,gasdev,  &
   gap,dumr,T,sinc,rcond,dtmax
real(double), allocatable, dimension(:) :: flux2,dts,trs,frs,B,S,work,  &
   pars,alpha,yerr2,mu
real(double), allocatable, dimension(:,:) :: A,Kernel,KernelZ
!FFTW3 vars
type(C_PTR) :: plan
integer ( kind = 4 ) :: nh
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: frsC

interface !creates a co-variance matrix
   subroutine makekernel(Kernel,npt1,npt2,x1,x2,npt,yerr,npars,pars)
      use precision
      implicit none
      integer, intent(inout) :: npt1,npt2,npars,npt
      real(double), dimension(:), intent(inout) :: x1,x2,yerr,pars
      real(double), dimension(:,:), intent(inout) :: Kernel
   end subroutine makekernel
end interface

!resampling routine
! 1-linear interpolation
! 2-sinc interpolation (not working)
! 3-Gaussian-process predictive Mean
iresampletype=1

call itime(now)
seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
dumr=ran2(-seed)

nover=8 !oversampling for FFT
gap=10.0d0 !identifying gaps in data and replace with white noise.

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

!range for resampled times.
mintime=minval(time(1:npt))
maxtime=maxval(time(1:npt))

ns=(maxtime-mintime)/ddt !number of resampled data points
!get an estimate of array size for FFTW that power is a 2
nfft=2**int(log10(dble(npt*nover))/log10(2.0d0)+1.0d0)
write(0,*) "nfft: ",nfft
allocate(trs(ns),frs(nfft)) !allocate space for resampled
frs=0.0d0 !fill with zeros to insure zero-padding

if(iresampletype.eq.1)then

   !calculate interpolated values using simple linear model
   i1=2
   loop=0
   do i=1,ns
      ts=mintime+ddt*dble(i-1)

      do while (loop.eq.0)
         if(time(i1).ge.ts)then
            loop=1
            dtime=time(i1)-time(i1-1)
            if(dtime.gt.gap*ddt)then
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

elseif(iresampletype.eq.2)then

   nbuf=1000 !break up resampling to managable size
   nstep=(ns/nbuf)+1 !number of steps to cover sample size

   do k=1,nstep
      !number of data points in buffer
      nss=nbuf
      if(k.eq.nstep)then !last step will have less points
         nss=ns-nbuf*(nstep-1)
      endif

      !times for resampling
      nbufstep=nbuf*(k-1) !location in resampled array
      do i=1,nss
         trs(i+nbufstep)=mintime+ddt*dble(i-1+nbufstep)
      enddo

      T=ddt !uniform sample rate
      !find range of unsampled data
      i1=1
      loop=0
      do while(loop.eq.0)
         if(time(i1).ge.trs(nbufstep+1))then
            i1=i1-1
            loop=1
         elseif(i1.gt.ns)then
            write(0,*) "sinc interpolator failed.."
            write(0,*) "i1 > ns"
            stop
         else
            i1=i1+1
         endif
      enddo

      i2=1
      loop=0
      do while(loop.eq.0)
         if(time(i2).ge.trs(nbuf+nbufstep))then
            loop=1
!         elseif(i2.gt.ns)then
!            write(0,*) "sinc interpolator failed.."
!            write(0,*) "i2 > ns"
!            stop
         else
            i2=i2+1
         endif
      enddo

      i1=max(1,i1)
      i2=min(ns,i2)

      npts=i2-i1+1 !how much of the original sample are we using.

      write(0,*) "iX: ",i1,i2
      write(0,*) "nss,npts: ",nss,npts
!      read(5,*)

      nB=max(nss,npts)
      !we are going to solve A*x=b with 2-norm(| b - A*x |) being minimized
      !this is a least-squares problem
      allocate(A(npts,nss),B(nB))
      write(0,*) "Setting up Matrices for LSQ"
      B(1:npts)=flux(i1:i2)
      do j=1,nss
!         !$OMP PARALLEL DO PRIVATE(TS)
         do i=i1,i2 !loop for npts
            ii=i-i1+1
            ts=trs(j+nbufstep)
            A(ii,j)=sinc((time(i)-ts)/T)
!            write(0,*) "1:",time(i),ts,T
!            write(0,*) "2:",i,j,A(ii,j)
!            read(5,*)
         enddo
!         !$OMP END PARALLEL DO
      enddo

!      call DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,          &
!         WORK, LWORK, IWORK, INFO )
      nrhs=1
      rcond=-1.0d0
      SMLSIZ=25
      NLVL = MAX(0,INT( LOG10( DBLE(MIN(npts,nss)/(SMLSIZ+1) ))/log10(2.0d0))+1)
!     NLVL = MAX(0,INT( LOG_2(      MIN(M,N)     /(SMLSIZ+1)  )             )+1)
      LWORK=3*MIN(npts,nss)*NLVL + 11*MIN(npts,nss)
      allocate(S(min(npts,nss)),work(1),IWORK(LWORK))
      LWORK=-1
      call DGELSD(npts,nss,NRHS,A,npts,B,nB,S,RCOND,RANK,          &
         WORK, LWORK, IWORK, INFO )
!      write(0,*) "LWORK: ",LWORK,int(WORK(1))
      LWORK=max(1,int(WORK(1)))
      deallocate(WORK)
      allocate(WORK(LWORK))
      call DGELSD(npts,nss,NRHS,A,npts,B,nB,S,RCOND,RANK,          &
         WORK, LWORK, IWORK, INFO )
      write(0,*) "INFO: ",INFO

      !output for debugging
      open(unit=11,file='sinc.dat')
      do i=1,nss
         write(11,*) trs(i),b(i)
      enddo
      close(11)

      deallocate(A,B,S,WORK,IWORK)
      write(6,*) "sinc interpolator done.."
      read(5,*)

   enddo

else
!resampling with Gaussian Process

   nbuf=10000 !break up resampling to managable size
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
         trs(i+nbufstep)=mintime+ddt*dble(i-1+nbufstep)
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
      pars(2)=ddt !length scale for exp
      pars(3)=dtmax/ddt !second amp scale
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
      call makekernel(KernelZ,nss,npts,trs(1+nbufstep:nss+nbufstep),    &
         time(i1:i2),npts,yerr2,npars,pars)
      mu=matmul(KernelZ,alpha)

!     fill in big gaps with noise.
      j=i1
      loop=0
      do i=2,nss
         ts=trs(i+nbufstep) !time at sample point
         do while (loop.eq.0)
            if(time(j).ge.ts)then
               loop=1
               dtime=time(j)-time(j-1)
!               write(0,*) dtime,gap*ddt
               if(dtime.gt.gap*ddt)then
!                  write(0,*) "gap.."
                  mu(i)=gasdev(seed)*(ferr(j-1)+ferr(j))/2.0d0
               endif
            else
               j=j+1
            endif
         enddo
         loop=0
      enddo

      frs(1+nbufstep:nss+nbufstep)=mu(1:nss)

      !output for debugging
!      open(unit=11,file="gp.dat")
!      do i=1,nss
!         write(11,*) trs(i+nbufstep),mu(i)
!      enddo
!      close(11)

      deallocate(Kernel,alpha,KernelZ,yerr2,mu)
      write(0,*) "Gaussian predictor done ",k,"/",nstep
!      read(5,*)
   enddo

endif


open(unit=11,file="interpolate.dat")
do i=1,ns
   write(11,*) trs(i),frs(i)
enddo
close(11)

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
