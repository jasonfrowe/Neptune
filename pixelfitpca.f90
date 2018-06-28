program pixelfitpca
use precision
implicit none
integer :: iargc,nmax,npt
real(double) :: minx,mean
real(double), allocatable, dimension(:) :: x,y,yerr,xpos,ypos,xnep,ynep,&
 oflux,pmod,smod
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

end program pixelfitpca