subroutine gradient(nfit,f,g,nfitp,isol,sol,npars,pars,npix,npord,npt,  &
   x,y,yerr,npixel,Kernel,logDK)
use precision
implicit none
!import vars
integer :: nfit,nfitp,npars,npix,npord,npt,ikch
integer, dimension(:) :: isol,npixel
real(double) :: logDK,f
real(double), dimension(:) :: g,sol,x,y,yerr,pars
real(double), dimension(:,:) :: Kernel
!local vars
integer i,j,k,info
real(double) :: small,ftest1,ftest2,loglikelihood,logDK1
real(double), allocatable, dimension(:) :: soltest,r,pars1
real(double), allocatable, dimension(:,:) :: Kernel1

interface !pixel/jump model
   subroutine pixelmodelv2(r,npars,npix,npord,sol,npt,x,npixel)
      use precision
      implicit none
      integer :: npars,npix,npord,npt
      integer, dimension(:) :: npixel
      real(double), dimension(:) :: sol,x,r
   end subroutine pixelmodelv2
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

small=1.0d-8

allocate(soltest(nfitp),pars1(npars),Kernel1(npt,npt))
allocate(r(npt))

ikch=0 !flag for updating Kernel
j=0
do i=1,nfitp
   if(isol(i).ne.0)then
      j=j+1 !update counter

      soltest=sol !make copy of original solution
      soltest(i)=soltest(i)+small

      do k=1,npars !copy Kernel parameters
         pars1(k)=pars(k)
      enddo
      logDK1=logDK  !copy determinant
      Kernel1=Kernel !copy Kernel

      if(i.le.npars)then !see if Kernel needs updating
         pars1(i)=soltest(i)
         ikch=1
      endif

      if(ikch.eq.1)then !if Kernel needs updating, then update it.
         write(0,*) "Updating Kernel1"
         call makekernel(Kernel1,npt,npt,x,x,npt,yerr,npars,pars1)
         call dpotrf('U',npt,Kernel1,npt,info) !LAPACK routine for Cholesky
         if (info.ne.0) then !check for errors
            write(0,*) "Cholesky factorization failed for gradient"
            write(0,*) "dpotrf info: ",info
            stop
         endif
         logDK1=0.0d0
         do k=1,npt
            logDK1=logDK1+Kernel1(k,k)  !when factorized, determinant is sum of diagonal
         enddo
         logDK1=2.0*logDK1
      endif

      call pixelmodelv2(r,npars,npix,npord,soltest,npt,x,npixel)
      ftest1=-loglikelihood(npt,x,y,r,Kernel1,logDK1)
!      write(0,*) j,f,ftest1
      g(j)=(ftest1-f)/small

!      soltest(i)=soltest(i)-2.0d0*small
!      call pixelmodelv2(r,npars,npix,npord,soltest,npt,x,npixel)
!      ftest2=-loglikelihood(npt,x,y,r,Kernel1,logDK1)
!      g(j)=(ftest1-ftest2)/(2.0d0*small)
   endif
   ikch=0
enddo

return
end
