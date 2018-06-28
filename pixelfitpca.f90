program pixelfitpca
use precision
implicit none
integer :: iargc,nmax,npt,ixo,iyo,npars,info,nrhs,i
real, allocatable, dimension(:) :: bb
real(double) :: minx,mean,pixel,poly
real(double), allocatable, dimension(:) :: x,y,yerr,xpos,ypos,xnep,ynep,&
 oflux,pmod,smod,ax,ay,xpcor,ypcor,pars,alpha,yerr2,mu,std,res,phase
real(double), allocatable, dimension(:,:) :: Kernel,KernelZ
character(80) :: filename

!These are F90 interfaces that allow one to pass assumed sized arrays
!to subroutines.
interface !reads in a three-column ascii space seperated file
   subroutine getpdata(filename,npt,nmax,x,y,yerr,oflux,pmod,smod,xpos, &
    ypos,xnep,ynep,minx,mean)
      use precision
      implicit none
      integer, intent(inout) :: npt,nmax
      real(double), intent(inout) :: minx,mean
      real(double), dimension(:), intent(inout) :: x,y,yerr,oflux,pmod, &
         smod,xpos,ypos,xnep,ynep
      character(80), intent(inout) :: filename
   end subroutine getpdata
   subroutine fitneptunepos(npt,x,xnep,ynep,ixo,ax,iyo,ay)
      use precision
      implicit none
      integer :: npt,ixo,iyo
      real(double), dimension(:) :: x,xnep,ynep,ax,ay
   end subroutine fitneptunepos
   subroutine makekernel(Kernel,npt1,npt2,x1,x2,npt,yerr,npars,pars)
      use precision
      implicit none
      integer, intent(inout) :: npt1,npt2,npars,npt
      real(double), dimension(:), intent(inout) :: x1,x2,yerr,pars
      real(double), dimension(:,:), intent(inout) :: Kernel
   end subroutine makekernel
   subroutine plotdatascatter(npt,x,y,yerr,bb)
      use precision
      implicit none
      integer, intent(inout) :: npt
      real, dimension(:), intent(inout) :: bb
      real(double), dimension(:), intent(inout) :: x,y,yerr
   end subroutine plotdatascatter
   subroutine plotsamples(npt1,x1,mu,std)
      use precision
      implicit none
      integer, intent(in) :: npt1
      real(double), dimension(:), intent(in) :: x1,mu,std
   end subroutine plotsamples
end interface

!check that we have enough information from the commandline
if(iargc().lt.1)then !if not, spit out the Usage info and stop.
   write(0,*) "Usage: pixelfit filename"
   stop
endif

!read in filename containing data (3 columns)
call getarg(1,filename)

nmax=80000 !initial guess for number of datapoints.
allocate(x(nmax),y(nmax),yerr(nmax),oflux(nmax),pmod(nmax),smod(nmax),  &
   xpos(nmax),ypos(nmax),xnep(nmax),ynep(nmax))

call getpdata(filename,npt,nmax,x,y,yerr,oflux,pmod,smod,xpos,ypos,xnep,ynep,minx,mean)
write(0,*) "Number of points read: ",npt !report how much data was read in

!open PGPLOT device
call pgopen('/xserve')  !'?' lets the user choose the device.
call PGPAP (8.0 ,1.0) !use a square 8" across
call pgpage() !create a fresh page
call pgslw(3) !thicker lines
call pgask(.true.)

!need fits to the motion of Neptune - raw data is too noisy.
!also means that pixel-crossings is now a linear function (good!)
ixo=5 !order of polynomial to fit to x-positions
iyo=5 !order of polynomial to fit to y-positions
allocate(ax(ixo),ay(iyo))
allocate(xpcor(npt),ypcor(npt))
xpcor=xnep+xpos !correct for pointing jitter
ypcor=ynep+ypos
call fitNeptunePos(npt,x,xpcor,ypcor,ixo,ax,iyo,ay)
deallocate(xpcor,ypcor)

!plot the data
allocate(bb(4))
bb=0.0e0 !rescale plot
write(0,*) "plotting.."
call plotdatascatter(npt,x,y,yerr,bb) !plot our original dataset

!Here are the parameters that control the co-variance matrix and fitted
!parameters
npars=4 !number of parameters used for model of the matrix
allocate(pars(npars))  !model parameters for Kernel generation.
pars(1)=1.0d0 !amp scale for exp
pars(2)=0.72d0 !length scale for exp
pars(3)=1.0 !second amp scale
pars(4)=2.0 !second length scale

!lets make a Kernel/co-variance for the Gaussian process
allocate(Kernel(npt,npt)) !allocate space
call makekernel(Kernel,npt,npt,x,x,npt,yerr,npars,pars) !create Kernel

!Cholesky factorization
write(0,*) "Beginning Cholesky factorization"
call dpotrf('U',npt,Kernel,npt,info) !LAPACK routine for Cholesky
if (info.ne.0) then !check for errors
   write(0,*) "Cholesky factorization failed"
   write(0,*) "dpotrf info: ",info
   stop
endif

write(0,*) "Beginning Cholesky Solver"
!calculate solution
! We have to solve,
! Kernel*X=Y,
! for X.
! LAPACK routine dpotrs uses Kfac from dpotrf.  We copy y into a new
! array called alpha, which is overridden with the solution on completion.
allocate(alpha(npt))
alpha=y(1:npt) !dpotrs takes alpha as input and output.
nrhs=1 !how man columns does alpha contain - just one
call dpotrs('U',npt,nrhs,Kernel,npt,alpha,npt,info) !call LAPACK
if (info.ne.0) then !check for errors
   write(0,*) "Solver failed.."
   write(0,*) "dpotrs info: ",info
   stop
endif

!predict samples
allocate(KernelZ(npt,npt))
allocate(yerr2(npt),mu(npt),std(npt))
yerr2=0.0d0
call makekernel(KernelZ,npt,npt,x,x,npt,yerr2,npars,pars)
allocate(res(npt))

!$OMP PARALLEL
!$OMP WORKSHARE
mu=matmul(KernelZ,alpha)
std=0.0d0
res(1:npt)=y(1:npt)-mu(1:npt) !contains 'detrended' data
!$OMP END WORKSHARE
!$OMP END PARALLEL

write(0,*) "done plotting"
call plotsamples(npt,x,mu,std) !plot our predicted sample set on top.

allocate(phase(npt))
!write(0,*) "xnep: ",xnep(1),xnep(100)
!phase(1:npt)=xnep(1:npt)-floor(xnep(1:npt))
!$OMP PARALLEL DO
do i=1,npt
   !pixel=poly(x(i),ixo,ax)
   pixel=xnep(i)!+xpos(i)
   phase(i)=pixel-floor(pixel)
enddo
!$OMP END PARALLEL DO

call pgpage()
bb=0.0
call plotdatascatter(npt,phase,res,yerr,bb)

call pgclos()

end program pixelfitpca

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function poly(x,nfit,ans)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: nfit
real(double) :: x,poly
real(double), dimension(nfit) :: ans
integer :: i

poly=ans(1)
do i=2,nfit
   poly=poly+ans(i)*x**dble(i-1)
enddo

return
end function