subroutine pfitter(npt,x,res,yerr,ixo,ax,nfitp,sol,isol)
use precision
implicit none
integer :: npt,ixo,nfitp
integer, dimension(:) :: isol
real(double), dimension(:) :: x,res,yerr,ax,sol
!local vars
integer nfit,i,j,n,m,iprint,isave(44)
integer, allocatable, dimension(:) :: nbd,iwa
real :: twork
real(double), parameter :: pi = 3.1415926535897932_8
real(double) :: twoPi,factr,pgtol,f,dsave(29),loglikelihood
real(double), allocatable, dimension(:) :: solin,sol1,l,u,g,wa,psmod
character(60) :: task,csave
logical :: lsave(4)

interface
   subroutine pshapemodel(npt,x,ixo,ax,nfit,sol,psmod)
      use precision
      implicit none
      integer :: npt,ixo,nfit
      real(double), dimension(:) :: x,ax,sol,psmod
   end subroutine pshapemodel
end interface

twoPi=2.0d0*Pi

nfit=0 !number of variables that are fit
do i=1,nfitp
   if(isol(i).ne.0)then
      nfit=nfit+1
   endif
enddo
write(0,*) "nfit: ",nfit

n=nfit !number of variables
m=5 !corrections used in limited memory matrix
!set up parameters for fit
allocate(solin(nfit),sol1(nfitp))
allocate(l(n),u(n),nbd(n))

write(0,*) "Hello.."
j=0
do i=1,nfitp
   if(isol(i).ne.0)then
      j=j+1
      solin(j)=sol(i) !pick off variables that will be fit
      nbd(j)=0 !no bounds for line-segments
      if(i.eq.2)then  !set bounds for phase [0:2*Pi]
         l(j)=0.0
         u(j)=twoPi
         nbd(j)=2
      endif
      if(i.ge.2)then !set bounds for Kernel parameters
         l(j)=1.0e-8 !scales must be greater than zero.
         nbd(j)=1 !lower bounds for Kernel parameters
      endif
   endif
enddo

allocate(g(nfit))
factr=1.0d+7
pgtol=1.0d-5
allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m) )
allocate ( iwa(3*n) )
iprint=-1  !diagonistic info to print to screen (set negative to quiet)

allocate(psmod(npt))

task = 'START'

do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
               task.eq.'START')

   call setulb ( n, m, solin, l, u, nbd, f, g, factr, pgtol, &
                       wa, iwa, task, iprint,&
                       csave, lsave, isave, dsave )

   write(0,'(A6,A6)') "task: ",task

   if (task(1:2) .eq. 'FG') then
      j=0
      do i=1,nfitp !loop over all parameters to find fitted parameters
         if(isol(i).ne.0)then
            j=j+1
            sol1(i)=solin(j)  !update fitted parameters
         else
            sol1(i)=sol(i) !fill in from beginning solution.
         endif
      enddo

      !get model
      call pshapemodel(npt,x,ixo,ax,nfitp,sol1,psmod)
!      write(0,*) "pshapemodel done.. get f"
      f=-loglikelihood(npt,res,yerr,psmod)
      CALL CPU_TIME(twork)
      write(0,*) "F: ",f,twork
      call gradient(nfitp,sol1,isol,f,nfit,g,npt,x,res,yerr,ixo,ax)
      write(0,*) "G(1,2): ",g(1),g(2)

   endif
   write(0,*) solin(1),solin(2)

enddo

sol=sol1 !update solution.

open(unit=11,file="soltest.dat")
do i=1,nfitp
   write(11,*) i,sol(i)
enddo
close(11)

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine gradient(nfitp,sol1,isol,f,nfit,g,npt,x,res,yerr,ixo,ax)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: nfitp,isol(nfitp),nfit,npt,ixo
real(double) :: f,g(nfit),sol1(nfitp),x(npt),res(npt),yerr(npt),ax(ixo)
!local vars
integer :: i,j
real(double) :: ftest1,small,loglikelihood
real(double), allocatable, dimension(:) :: soltest,psmod

interface
   subroutine pshapemodel(npt,x,ixo,ax,nfit,sol,psmod)
      use precision
      implicit none
      integer :: npt,ixo,nfit
      real(double), dimension(:) :: x,ax,sol,psmod
   end subroutine pshapemodel
end interface

allocate(soltest(nfitp),psmod(npt))

small=1.0d-8

j=0
do i=1,nfitp
   soltest(1:nfitp)=sol1(1:nfitp)
   if(isol(i).ne.0)then
      j=j+1
      soltest(i)=sol1(i)+small
      call pshapemodel(npt,x,ixo,ax,nfitp,soltest,psmod)
      ftest1=-loglikelihood(npt,res,yerr,psmod)
      g(j)=(ftest1-f)/small
   endif
enddo

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function loglikelihood(npt,res,yerr,psmod)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: npt,i
real(double) :: loglikelihood,logdet,log2pi
real(double), parameter :: pi = 3.1415926535897932_8
real(double), dimension(npt) :: res,yerr,psmod

log2pi=log(2.0d0*pi)
logdet=log(sum(yerr(1:npt)))
!write(0,*) "logdet: ",logdet

loglikelihood=0.0d0
do i=1,npt
   loglikelihood=loglikelihood+(res(i)-psmod(i))*(res(i)-psmod(i))/     &
     (yerr(i)*yerr(i))
!   write(0,*) i,loglikelihood,(res(i)-psmod(i)),yerr(i)
!   read(5,*)
enddo

loglikelihood=-0.5*(loglikelihood+logdet+dble(npt)*log2pi)

return
end
