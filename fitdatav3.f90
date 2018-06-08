program fitdatav3
!Version 2 to Neptune data cleaner.
use precision
implicit none
integer :: iargc,nmax,npt
real(double) :: minx,mean
real(double), allocatable, dimension(:) :: x,y,yerr,xpos,ypos,xnep,ynep
character(80) :: filename
!Temp vars
integer :: i

!These are F90 interfaces that allow one to pass assumed sized arrays
!to subroutines.
interface !reads in a three-column ascii space seperated file
   subroutine getdata(filename,npt,nmax,x,y,yerr,xpos,ypos,xnep,ynep,minx,mean)
      use precision
      implicit none
      integer, intent(inout) :: npt,nmax
      real(double), intent(inout) :: minx,mean
      real(double), dimension(:), intent(inout) :: x,y,yerr,xpos,ypos,  &
         xnep,ynep
      character(80), intent(inout) :: filename
   end subroutine getdata
end interface
interface
   subroutine cutoutliers(npt,x,y,yerr,xpos,ypos,xnep,ynep)
      use precision
      implicit none
      integer :: npt
      real(double), dimension(:) :: x,y,yerr,xpos,ypos,xnep,ynep
   end subroutine cutoutliers
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
!subroutine to read in data
call getdata(filename,npt,nmax,x,y,yerr,xpos,ypos,xnep,ynep,minx,mean)
write(0,*) "Number of points read: ",npt !report how much data was read in

!cut out crap.
call cutoutliers(npt,x,y,yerr,xpos,ypos,xnep,ynep)


do i=1,npt
	write(6,*) x(i),y(i)
enddo

end program fitdatav3