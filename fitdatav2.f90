program fitdatav2
!Version 2 to Neptune data cleaner.
use precision
implicit none
integer :: iargc,nmax,npt
real(double), allocatable, dimension(:) :: x,y,yerr,xpos,ypos,xnep,ynep
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

end program fitdatav2
