subroutine fitterv2(npt,x,y,yerr,npixel,npars,pars,npord)
use precision
!use fittingmodv2
implicit none
!input vars
integer :: npt,npars,npord
integer, dimension(:) :: npixel
real(double), dimension(:) :: x,y
real(double), dimension(:) :: yerr,pars
!shared vars
integer :: nfitp
integer, allocatable, dimension(:) :: isol
real(double) :: logDK
real(double), allocatable, dimension(:) :: alpha
real(double), allocatable, dimension(:) :: sol
real(double), allocatable, dimension(:,:) :: Kernel
!local vars
integer, allocatable, dimension(:) :: nbd,iwa
integer :: nrhs,info,nfit,i,j,k,npix,iprint,isave(44),n,m
real, allocatable, dimension(:) :: bb
real(double) :: loglikelihood,lgll,f,factr,pgtol,dsave(29)
real(double), allocatable, dimension(:) :: yerr2,mu,std,ymodel,r,solin, &
   l,u,g,wa,sol1
logical :: lsave(4)

real(double), allocatable, dimension(:,:) :: KernelZ
character(60) :: task,csave
!external funk

interface !creates a co-variance matrix
   subroutine makekernel(Kernel,npt1,npt2,x1,x2,npt,yerr,npars,pars)
      use precision
      implicit none
      integer, intent(inout) :: npt1,npt2,npars,npt
      real(double), dimension(:), intent(inout) :: x1,x2,yerr,pars
      real(double), dimension(:,:), intent(inout) :: Kernel
   end subroutine makekernel
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
interface !plots samples and uncertainties
   subroutine plotsamples(npt1,x1,mu,std)
      use precision
      implicit none
      integer, intent(in) :: npt1
      real(double), dimension(:), intent(in) :: x1,mu,std
   end subroutine plotsamples
end interface
interface !pixel/jump model
   subroutine pixelmodelv2(r,npars,npix,npord,sol,npt,x,npixel)
      use precision
      implicit none
      integer :: npars,npix,npord,npt
      integer, dimension(:) :: npixel
      real(double), dimension(:) :: sol,x,r
   end subroutine pixelmodelv2
end interface

!set up Kernel
!lets make a Kernel/co-variance for the Gaussian process
allocate(Kernel(npt,npt)) !allocate space
call makekernel(Kernel,npt,npt,x,x,npt,yerr,npars,pars) !create Kernel

!Cholesky factorization
!write(0,*) "Beginning Cholesky factorization in fitter"
call dpotrf('U',npt,Kernel,npt,info) !LAPACK routine for Cholesky
if (info.ne.0) then !check for errors
   write(0,*) "Cholesky factorization failed in fitterv2"
   write(0,*) "dpotrf info: ",info
   stop
endif
!at this point, Kernel has been factorized.

!precompute determinant
logDK=0.0d0
do i=1,npt
   logDK=logDK+Kernel(i,i)  !when factorized, determinant is sum of diagonal
enddo
logDK=2.0*logDK

!pre-compute alpha
allocate(alpha(npt))
alpha=y(1:npt) !dpotrs takes alpha as input and output.
nrhs=1 !how man columns does alpha contain - just one
call dpotrs('U',npt,nrhs,Kernel,npt,alpha,npt,info) !call LAPACK
if (info.ne.0) then !check for errors
   write(0,*) "Solver failed in fitterv2.."
   write(0,*) "dpotrs info: ",info
   stop
endif
!predictive means
allocate(KernelZ(npt,npt),yerr2(npt),mu(npt),std(npt))
yerr2=0.0d0
call makekernel(KernelZ,npt,npt,x,x,npt,yerr2,npars,pars)
mu=matmul(KernelZ,alpha)
std=0.0d0
!plot the prediction
allocate(bb(4))
bb=0.0e0 !rescale plot
write(0,*) "plotting.."
call plotdatascatter(npt,x,y,yerr,bb) !plot our original dataset
write(0,*) "done plotting"
call plotsamples(npt,x,mu,std) !plot our predicted sample set on top.

!number of jumps detected - sets number of parameters for segment fit
npix=maxval(npixel(1:npt))
!total number of parameters that define the model.
nfitp=npars+npix*npord
allocate(sol(nfitp),isol(nfitp))
isol=-1 !initiate isol array.  if isol(i).ne.0 then the variable is fit.
do i=1,npars
   sol(i)=pars(i)
   isol(i)=0 !keep Kernel constant.
enddo
!add in X-pixel model fit.
do i=1+npars,npars+npix*npord
   sol(i)=0.0d0 !start with a straight line for a guess
enddo

!calculate how many variables are fit.
nfit=0
do i=1,nfitp
   if(isol(i).ne.0)then
      nfit=nfit+1
   endif
enddo


n=nfit !number of variables
m=5 !corrections used in limited memory matrix
!set up parameters for fit
allocate(solin(nfit),sol1(nfitp))
j=0
do i=1,nfitp
   if(isol(i).ne.0)then
      j=j+1
      solin(j)=sol(i) !pick off variables that will be fit
   endif
enddo
allocate(l(npt),u(npt),nbd(npt))
nbd=0 !no bounds set
allocate(g(npt))
factr=1.0d+7
pgtol=1.0d-5
allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m) )
allocate ( iwa(3*n) )
iprint=1  !diagonistic info to print to screen (set negative to quiet)

!vars for model
allocate(r(npt))

task = 'START'

do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
               task.eq.'START')

   call setulb ( n, m, solin, l, u, nbd, f, g, factr, pgtol, &
                       wa, iwa, task, iprint,&
                       csave, lsave, isave, dsave )

   write(0,'(A6,A6)') "task: ",task

   if (task(1:2) .eq. 'FG') then
      j=0
      do i=1,nfitp
         if(isol(i).ne.0)then
            j=j+1
            sol1(i)=solin(j)
         else
            sol1(i)=sol(i)
         endif
      enddo
   call pixelmodelv2(r,npars,npix,npord,sol1,npt,x,npixel)
   f=-loglikelihood(npt,x,y,r,Kernel,logDK)


   endif


   read(5,*)
enddo

return
end subroutine fitterv2

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function funk(p1)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
use fittingmodv2
implicit none
!import vars
real(double) :: funk,loglikelihood
real(double), dimension(npt2) :: p1
!local vars
integer i,j
real(double), allocatable, dimension(:) :: sol3,r

interface !pixel/jump model
   subroutine pixelmodelv2(r,npars,npix,npord,sol,npt,x,npixel)
      use precision
      implicit none
      integer :: npars,npix,npord,npt
      integer, dimension(:) :: npixel
      real(double), dimension(:) :: sol,x,r
   end subroutine pixelmodelv2
end interface


!write(0,*) "Starting func.."

allocate(sol3(nfitp2),r(npt2))

j=0 !counter
do i=1,nfitp2
   if(isol2(i).ne.0)then
      j=j+1
      sol3(i)=p1(j)
   else
      sol3(i)=sol2(i)
   endif
enddo

!write(0,*) "Call pixelmodel"
call pixelmodelv2(r,npars2,npix2,npord2,sol3,npt2,x2,npixel2)
!write(0,*) "Estimate lgll"
funk=-loglikelihood(npt2,x2,y2,r,Kernel2,logDK2)

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function loglikelihood(npt,x,y,r,Kernel,logDK)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!Kernel needs to factorized
use precision
implicit none
!import vars
integer :: npt
real(double) :: logDK
real(double), dimension(npt) :: x,y,r
real(double), dimension(npt,npt) :: Kernel
!local vars
integer :: nrhs,info,i
real(double) :: loglikelihood,log2pi
real(double), parameter :: pi = 3.1415926535897932_8
real(double), allocatable, dimension(:) :: ymr

log2pi=log(2.0d0*pi)

allocate(ymr(npt)) !observations minus model (y-r)
ymr(1:npt) = y(1:npt)-r(1:npt)

!solve for K**-1 * ymr
nrhs=1
!input ymr gets over written with K**-1 * ymr
call dpotrs('U',npt,nrhs,Kernel,npt,ymr,npt,info) !call LAPACK solver
if (info.ne.0) then !check for errors
   write(0,*) "Solver failed in fcn.."
   write(0,*) "dpotrs info: ",info
   stop
endif

!calculate ymr * K**-1 * ymr
loglikelihood=0.0
do i=1,npt
   loglikelihood=loglikelihood+(y(i)-r(i))*ymr(i)
enddo
loglikelihood=-0.5d0*(loglikelihood+logDK+dble(npt)*log2pi)


return
end

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
