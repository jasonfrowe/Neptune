subroutine fcn(npt,nfit,solin,fvec,iflag)
use precision
use fittingmod
implicit none
integer :: npt,nfit,iflag
real(double), dimension(npt) :: fvec
real(double), dimension(nfit) :: solin
!local variables
integer i,j
real(double), allocatable, dimension(:) :: sol3,r

interface
   subroutine pixelmodel(r,npars2,npord2,sol3,npt,x2,ixo2,ax2,iyo2,ay2)
      use precision
      implicit none
      integer, intent(in) :: npord2,npt,ixo2,iyo2,npars2
      real(double), dimension(:), intent(in) :: sol3,ax2,ay2,x2
      real(double), dimension(:), intent(inout) :: r
   end subroutine pixelmodel
end interface

!sol3 will contain updated solution for model.
allocate(sol3(nfitp2))
!copy over solution from module share
sol3(1:nfitp2)=sol2(1:nfitp2)
j=0
do i=1,nfitp2
   if(isol2(i).ne.0)then
      j=j+1
      sol3(i)=solin(j)
   endif
!   write(0,*) i,sol3(i)
!   read(5,*)
enddo

!1. now we calculate y-model
!allocate space for model.
allocate(r(npt))

!initialize model to zero.
r=0.0d0
!get pixel model.  Expects 'r'to be initiated.
call pixelmodel(r,npars2,npord2,sol3,npt,x2,ixo2,ax2,iyo2,ay2)

!1b. if Kernel is changing, then we need to update Cholosky stuff.

!2. calculate Kernel-1 x r

!3. break out terms of r x (Kernel-1 x r) into fvec

!4. profit


!until code is ready we have a read statement for debugging.
read(5,*)

return
end subroutine fcn
