program pixelfit
use precision
implicit none
integer :: iargc,nmax,npt,npars
real, allocatable, dimension(:) :: bb
real(double) :: minx,mean
real(double), allocatable, dimension(:) :: x,y,yerr,xpos,ypos,xnep,ynep,&
 pars
real(double), allocatable, dimension(:,:) :: Kernel
character(80) :: filename

!These are F90 interfaces that allow one to pass assumed sized arrays
!to subroutines.
interface !reads in a three-column ascii space seperated file
   subroutine getpdata(filename,npt,nmax,x,y,yerr,xpos,ypos,xnep,ynep,minx,mean)
      use precision
      implicit none
      integer, intent(inout) :: npt,nmax
      real(double), intent(inout) :: minx,mean
      real(double), dimension(:), intent(inout) :: x,y,yerr,xpos,ypos,  &
         xnep,ynep
      character(80), intent(inout) :: filename
   end subroutine getpdata
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
interface !creates a co-variance matrix
   subroutine makekernel(Kernel,npt1,npt2,x1,x2,npt,yerr,npars,pars)
      use precision
      implicit none
      integer, intent(inout) :: npt1,npt2,npars,npt
      real(double), dimension(:), intent(inout) :: x1,x2,yerr,pars
      real(double), dimension(:,:), intent(inout) :: Kernel
   end subroutine makekernel
end interface

!check that we have enough information from the commandline
if(iargc().lt.1)then !if not, spit out the Usage info and stop.
   write(0,*) "Usage: pixelfit filename"
   stop
endif

!read in filename containing data (3 columns)
call getarg(1,filename)

nmax=80000 !initial guess for number of datapoints.
allocate(x(nmax),y(nmax),yerr(nmax),xpos(nmax),ypos(nmax),xnep(nmax),   &
   ynep(nmax))

call getpdata(filename,npt,nmax,x,y,yerr,xpos,ypos,xnep,ynep,minx,mean)
write(0,*) "Number of points read: ",npt !report how much data was read in

!open PGPLOT device
call pgopen('/xserve')  !'?' lets the user choose the device.
call PGPAP (8.0 ,1.0) !use a square 8" across
call pgpage() !create a fresh page
call pgslw(3) !thicker lines

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
pars(2)=0.69d0 !length scale for exp
pars(3)=1.0 !second amp scale
pars(4)=500.0 !second length scale

!lets make a Kernel/co-variance for the Gaussian process
allocate(Kernel(npt,npt)) !allocate space
call makekernel(Kernel,npt,npt,x,x,npt,yerr,npars,pars) !create Kernel

call pgclos()

end program pixelfit
