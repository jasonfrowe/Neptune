subroutine fcn(npt,nfit,solin,fvec,iflag)
use precision
use fittingmod
implicit none
integer :: npt,nfit,iflag,nrhs,info
real(double), dimension(npt) :: fvec
real(double), dimension(nfit) :: solin
!local variables
integer i,j
real(double), allocatable, dimension(:) :: sol3,r,ra,zeros
real(double), allocatable, dimension(:,:) :: Kernel

interface
   subroutine pixelmodel(r,npars2,npord2,sol3,npt,x2,ixo2,ax2,iyo2,ay2)
      use precision
      implicit none
      integer, intent(in) :: npord2,npt,ixo2,iyo2,npars2
      real(double), dimension(:), intent(in) :: sol3,ax2,ay2,x2
      real(double), dimension(:), intent(inout) :: r
   end subroutine pixelmodel
end interface

!write(0,*) "Starting pixelmodel"

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
open(unit=11,file="model.dat")
do i=1,npt
   write(11,*) x2(i),y2(i),r(i)
enddo
close(11)
!r contains model

!1b. if Kernel is changing, then we need to update Cholosky stuff.
!to be done later.

!2. calculate ra = Kernel-1 x r
allocate(ra(npt))
!ra(1:npt)=y2(1:npt)
ra(1:npt)=y2(1:npt)-r(1:npt)  !y-model
nrhs=1
call dpotrs('U',npt,nrhs,Kfac2,npt,ra,npt,info) !call LAPACK solver
if (info.ne.0) then !check for errors
   write(0,*) "Solver failed in fcn.."
   write(0,*) "dpotrs info: ",info
   stop
endif



!3. break out terms of r x (Kernel-1 x r) into fvec
do i=1,npt
   fvec(i)=sqrt(abs((y2(i)-r(i))*ra(i)))
!   write(0,*) i,y2(i)-r(i),fvec(i)
!   read(5,*)
enddo
!write(0,*) "ll: ",Sum(fvec(1:npt))

!allocate(zeros(npt),Kernel(npt,npt))
!zeros=0.0d0
!call makekernel(Kernel,npt,npt,x2,x2,npt,zeros,npars2,pars2)

!fvec=y2(1:npt)-matmul(KernelZ2,ra)-r


!4. profit

!deallocate(r,ra,zeros,Kernel)

!until code is ready we have a read statement for debugging.
!write(0,*) "So far, so good.. "
!read(5,*)

return
end subroutine fcn
