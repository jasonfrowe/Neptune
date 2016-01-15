program fitdata
use precision
implicit none
integer :: iargc,nmax,npt,npars,info,nrhs,i,ixo,iyo
real(double) :: bpix
real, allocatable, dimension(:) :: bb
real(double), allocatable, dimension(:) :: x,y,yerr,pars,alpha,yerr2,mu,&
   std,res,xpos,ypos,xnep,ynep,ax,ay
real(double), allocatable, dimension(:,:) :: Kernel,Kfac,newKernelT,cov
character(80) :: filename

!These are F90 interfaces that allow one to pass assumed sized arrays
!to subroutines.
interface !reads in a three-column ascii space seperated file
   subroutine getdata(filename,npt,nmax,x,y,yerr,xpos,ypos,xnep,ynep)
      use precision
      implicit none
      integer, intent(inout) :: npt,nmax
      real(double), dimension(:), intent(inout) :: x,y,yerr,xpos,ypos,  &
         xnep,ynep
      character(80), intent(in) :: filename
   end subroutine getdata
end interface
interface !makes a plot of your data.
   subroutine plotdata(npt,x,y,yerr,bb)
      use precision
      implicit none
      integer, intent(in) :: npt
      real, dimension(:), intent(inout) :: bb
      real(double), dimension(:), intent(in) :: x,y,yerr
   end subroutine plotdata
end interface
interface !makes a plot of your data.
   subroutine plotdatascatter(npt,x,y,yerr,bb)
      use precision
      implicit none
      integer, intent(in) :: npt
      real, dimension(:), intent(inout) :: bb
      real(double), dimension(:), intent(in) :: x,y,yerr
   end subroutine plotdatascatter
end interface
interface !fits a straight line to data
   subroutine fitline(npt,x,y,yerr,ans,eans,chisq)
      use precision
      implicit none
      integer, intent(in) :: npt
      real(double), intent(inout) :: chisq
      real(double), dimension(:), intent(in) :: x,y,yerr
      real(double), dimension(:), intent(inout) :: ans,eans
   end subroutine fitline
end interface
interface !plots a line based on output from fitline
   subroutine plotline(bb,ans,eans)
      use precision
      implicit none
      real, dimension(:), intent(in) :: bb
      real(double), dimension(:), intent(in) :: ans,eans
   end subroutine plotline
end interface
interface !creates a co-variance matrix
   subroutine makekernel(Kernel,npt1,npt2,x1,x2,npt,yerr,npars,pars)
      use precision
      implicit none
      integer, intent(in) :: npt1,npt2,npars,npt
      real(double), dimension(:), intent(in) :: x1,x2,yerr,pars
      real(double), dimension(:,:), intent(inout) :: Kernel
   end subroutine makekernel
end interface
interface !displays a 2D array as a picture
   subroutine displaykernel(nx,ny,Kernel,bpix)
      use precision
      implicit none
      integer, intent(inout) :: nx,ny
      real(double), intent(inout) :: bpix
      real(double), dimension(:,:), intent(inout) :: Kernel
   end subroutine displaykernel
end interface
interface !plots samples and uncertainties
   subroutine plotsamples(npt1,x1,mu,std)
      use precision
      implicit none
      integer, intent(in) :: npt1
      real(double), dimension(:), intent(in) :: x1,mu,std
   end subroutine plotsamples
end interface
interface
   subroutine fitter(npt,Kfac,npars,pars,x,y,yerr,xnep,ynep,ixo,ax,iyo, &
    ay)
      use precision
      implicit none
      integer :: npt,npars,ixo,iyo
      real(double), dimension(:) :: pars,x,y,yerr,xnep,ynep,ax,ay
      real(double), dimension(:,:) :: Kfac
   end subroutine fitter
end interface
interface
   subroutine fitneptunepos(npt,x,xnep,ynep,ixo,ax,iyo,ay)
      use precision
      implicit none
      integer :: npt,ixo,iyo
      real(double), dimension(:) :: x,xnep,ynep,ax,ay
   end subroutine fitneptunepos
end interface

!Here are the parameters that control the co-variance matrix and fitted
!parameters
npars=4 !number of parameters used for model of the matrix
allocate(pars(npars))  !model parameters for Kernel generation.
pars(1)=1.0d0 !amp scale for exp
pars(2)=0.146d0 !length scale for exp
pars(3)=100.0 !second amp scale
pars(4)=500.0 !second length scale

!check that we have enough information from the commandline
if(iargc().lt.1)then !if not, spit out the Usage info and stop.
   write(0,*) "Usage: fitdata filename"
   stop
endif

!read in filename containing data (3 columns)
call getarg(1,filename)

nmax=80000 !initial guess for number of datapoints.
allocate(x(nmax),y(nmax),yerr(nmax),xpos(nmax),ypos(nmax),xnep(nmax),   &
   ynep(nmax)) !allocate arrays
!it is assumed that the data is sorted wrt time.
call getdata(filename,npt,nmax,x,y,yerr,xpos,ypos,xnep,ynep) !subroutine to read in data
write(0,*) "Number of points read: ",npt !report how much data was read in

!open PGPLOT device
call pgopen('?')  !'?' lets the user choose the device.
call PGPAP (8.0 ,1.0) !use a square 8" across
call pgpage() !create a fresh page
call pgslw(3) !thicker lines

!need fits to the motion of Neptune - raw data is too noisy.
!also means that pixel-crossings is now a linear function (good!)
ixo=5 !order of polynomial to fit to x-positions
iyo=5 !order of polynomial to fit to y-positions
allocate(ax(ixo),ay(iyo))
call fitNeptunePos(npt,x,xnep,ynep,ixo,ax,iyo,ay)

!plot the data
allocate(bb(4)) !contains plot boundaries
bb=0.0e0 !tell code to generate scale for plot
call plotdatascatter(npt,x,y,yerr,bb) !plot data

!lets make a Kernel/co-variance for the Gaussian process
allocate(Kernel(npt,npt)) !allocate space
call makekernel(Kernel,npt,npt,x,x,npt,yerr,npars,pars) !create Kernel

!call pgpage() !create fresh page for plotting
!bpix=1.0e30 !cut off for bright pixels
!call displaykernel(npt,npt,Kernel,bpix) !show what the Kernel looks like

!Cholesky factorization
write(0,*) "Beginning Cholesky factorization"
allocate(Kfac(npt,npt)) !allocate array to contain Cholesky factorization
Kfac=Kernel  !make a copy of Kernel as dpotrf overwrites input
call dpotrf('U',npt,Kfac,npt,info) !LAPACK routine for Cholesky
if (info.ne.0) then !check for errors
   write(0,*) "Cholesky factorization failed"
   write(0,*) "dpotrf info: ",info
   stop
endif
deallocate(Kernel) !don't need the original Kernel anymore

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
call dpotrs('U',npt,nrhs,Kfac,npt,alpha,npt,info) !call LAPACK
if (info.ne.0) then !check for errors
   write(0,*) "Solver failed.."
   write(0,*) "dpotrs info: ",info
   stop
endif

allocate(Kernel(npt,npt))
allocate(yerr2(npt),mu(npt),std(npt))
yerr2=0.0d0
call makekernel(Kernel,npt,npt,x,x,npt,yerr2,npars,pars)
mu=matmul(Kernel,alpha)
std=0.0d0

!!uncomment to go after standard-deviations (slow....)
!allocate(newKernelT(npt,npt))  !allocate space for transposed Kernel
!newKernelT=transpose(Kernel) !create transpose
!allocate(cov(npt,npt))
!!make Kernel based on predicted dataset
!call makekernel(cov,npt,npt,x,x,npt,yerr,npars,pars)
!nrhs=npt !note, more than one column!
!write(0,*) "Beginning Cholesky Solver"
!call dpotrs('U',npt,nrhs,Kfac,npt,newKernelT,npt,info) !Solve
!if (info.ne.0) then !error check
!   write(0,*) "Solver failed with newKernelT"
!   write(0,*) "dpotrs info: ",info
!   stop
!endif
!cov=cov-matmul(Kernel,newKernelT)
!do i=1,npt
!   std(i)=sqrt(cov(i,i)) !use diagonals as estimate of standard deviation
!enddo
!deallocate(newKernelT)

deallocate(Kernel)

!call pgpage()!fresh page for plotting
!bb=0.0e0 !rescale plot
!call plotdata(npt,x,y,yerr,bb) !plot our original dataset
call plotsamples(npt,x,mu,std) !plot our predicted sample set on top.

!at this point.. everything looks good, so lets call the fitter.
write(0,*) "Calling the fitter"
call fitter(npt,Kfac,npars,pars,x,y,yerr,xnep,ynep,ixo,ax,iyo,ay)

!!lets have a look at X-position vs residuals
!call pgpage()
!allocate(res(npt))
!res(1:npt)=y(1:npt)-mu(1:npt)
!bb=0.0e0
!!xnep=xnep-floor(xnep)
!call plotdatascatter(npt,xnep,res,yerr2,bb)

!open(unit=10,file="newdata.dat")
!do i=1,npt
!   write(10,'(4(F17.11,1X))') x(i),xnep(i),res(i),yerr(i)
!enddo
!close(10)

!close plotter
call pgclos()

end program fitdata
