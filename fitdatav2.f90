program fitdatav2
!Version 2 to Neptune data cleaner.
use precision
implicit none
integer :: iargc,nmax,npt,i,npars,info,nrhs,npord
integer, allocatable, dimension(:) :: npixel
real, allocatable, dimension(:) :: bb
real(double), allocatable, dimension(:) :: x,y,yerr,xpos,ypos,xnep,ynep,&
   pars,alpha,yerr2,mu,std,res,r,sp
real(double), allocatable, dimension(:,:) :: Kernel,Kfac,KernelZ

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
      character(80), intent(inout) :: filename
   end subroutine getdata
end interface
!interface !makes a plot of your data.
!   subroutine plotdatascatter(npt,x,y,yerr,bb)
!      use precision
!      implicit none
!      integer, intent(inout) :: npt
!      real, dimension(:), intent(inout) :: bb
!      real(double), dimension(:), intent(inout) :: x,y,yerr
!   end subroutine plotdatascatter
!end interface
interface
   subroutine findjumps(npt,x,y,yerr,npixel)
      use precision
      implicit none
      integer :: npt
      integer, dimension(:) :: npixel
      real(double), dimension(:) :: x,y,yerr
   end subroutine findjumps
end interface
interface !creates a co-variance matrix
   subroutine makekernel(Kernel,npt1,npt2,x1,x2,npt,yerr,npars,pars)
      use precision
      implicit none
      integer, intent(inout) :: npt1,npt2,npars,npt
      real(double), dimension(:), intent(inout) :: x1,x2,yerr,pars
      real(double), dimension(:,:), intent(inout) :: Kernel
   end subroutine makekernel
end interface
interface
   subroutine cutoutliers(npt,x,y,yerr,xpos,ypos,xnep,ynep)
      use precision
      implicit none
      integer :: npt
      real(double), dimension(:) :: x,y,yerr,xpos,ypos,xnep,ynep
   end subroutine cutoutliers
end interface
interface
   subroutine fitterv2(npt,x,y,yerr,npixel,npars,pars,npord,r,sp)
      use precision
      implicit none
      integer :: npt,npars,npord
      integer, dimension(:) :: npixel
      real(double), dimension(:) :: x,y,yerr,pars,r,sp
   end subroutine fitterv2
end interface

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

!cut out crap.
call cutoutliers(npt,x,y,yerr,xpos,ypos,xnep,ynep)

!open PGPLOT device
!call pgopen('/xserve')  !'?' lets the user choose the device.
!call PGPAP (8.0 ,1.0) !use a square 8" across
!call pgpage() !create a fresh page
!call pgslw(3) !thicker lines

!Here are the parameters that control the co-variance matrix and fitted
!parameters
npars=4 !number of parameters used for model of the matrix
allocate(pars(npars))  !model parameters for Kernel generation.
pars(1)=1.0d0 !amp scale for exp
pars(2)=0.02d0 !length scale for exp
pars(3)=1.0 !second amp scale
pars(4)=500.0 !second length scale

!lets make a Kernel/co-variance for the Gaussian process
allocate(Kernel(npt,npt)) !allocate space
call makekernel(Kernel,npt,npt,x,x,npt,yerr,npars,pars) !create Kernel

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

!predict samples
allocate(KernelZ(npt,npt))
allocate(yerr2(npt),mu(npt),std(npt))
yerr2=0.0d0
call makekernel(KernelZ,npt,npt,x,x,npt,yerr2,npars,pars)
mu=matmul(KernelZ,alpha)
std=0.0d0

allocate(res(npt))
res(1:npt)=y(1:npt)-mu(1:npt)

!!plot the data
!allocate(bb(4)) !contains plot boundaries - type is real (not double!)
!bb=0.0e0 !tell code to generate scale for plot
!call plotdatascatter(npt,x,res,yerr,bb) !plot data

allocate(npixel(npt))
call findjumps(npt,x,res,yerr,npixel)

!writing some data to check calculations
!open(unit=11,file="pixeltest.dat")
!do i=2,npt
!   write(11,'(4(F17.11,1X),I3)') x(i),y(i),yerr(i),res(i),npixel(i)-npixel(i-1)
!enddo
!close(11)

!clean up arrays no longer needed
deallocate(KernelZ,res,mu,std,yerr2)

!now we can go back to Kernel parameters to describe the rotation
!modulation of Neptune

!Here are the parameters that control the co-variance matrix and fitted
!parameters
npars=4 !number of parameters used for model of the matrix
pars(1)=0.02d0 !amp scale for exp
pars(2)=0.02856!0.146d0 !length scale for exp
pars(3)=0.002d0 !second amp scale
pars(4)=500.0d0 !second length scale

write(0,*) "Calling fitter"
!polynomial order for segment fits
npord=2
allocate(r(npt),sp(npt))
call fitterv2(npt,x,y,yerr,npixel,npars,pars,npord,r,sp)

call exportdata(npt,x,y,yerr,r,sp,xpos,ypos,xnep,ynep)

!call pgend()

end program fitdatav2

