subroutine fitneptunepos(npt,x,xnep,ynep,ixo,ax,iyo,ay)
!jasonfrowe@gmail.com (2016)
!npt - number of data points
!x - time  vector
!xnep,ynep - position of Neptune (X,Y) pixels
!ixo,iyo - order of polynomial to fit.
use precision
implicit none
integer :: npt,ixo,iyo
real(double), dimension(:) :: x,xnep,ynep,ax,ay

integer :: i,j
integer, allocatable, dimension(:) :: iax,iay
real(double) :: chisq,f
real, allocatable, dimension(:) :: bb
real(double), allocatable, dimension(:) :: sig,res
real(double), allocatable, dimension(:,:) :: covar

interface !makes a plot of your data.
   subroutine plotdatascatter(npt,x,y,yerr,bb)
      use precision
      implicit none
      integer, intent(in) :: npt
      real, dimension(:), intent(inout) :: bb
      real(double), dimension(:), intent(in) :: x,y,yerr
   end subroutine plotdatascatter
end interface

!Fit X-positions
allocate(sig(npt),iax(ixo),covar(ixo,ixo))
iax=1 !fit for all co-efficients.
sig=0.02 !approximate error on position
call lfit(x,xnep,sig,npt,ax,iax,ixo,covar,ixo,chisq)
deallocate(covar)

allocate(res(npt))
do i=1,npt
   f=ax(1)
   do j=2,ixo
      f=f+ax(j)*x(i)**dble(j-1)
   enddo
   res(i)=xnep(i)-f
enddo
bb=0
allocate(bb(4))
bb=0.0d0
call plotdatascatter(npt,x,res,sig,bb)
call pgpage()

!Fit Y-positions
allocate(iay(iyo),covar(iyo,iyo))
iay=1
sig=0.5 !approximate error on position
call lfit(x,ynep,sig,npt,ay,iay,iyo,covar,iyo,chisq)
deallocate(covar)

do i=1,npt
   f=ay(1)
   do j=2,iyo
      f=f+ay(j)*x(i)**dble(j-1)
   enddo
   res(i)=ynep(i)-f
enddo
bb=0
call plotdatascatter(npt,x,res,sig,bb)
call pgpage()

return
end subroutine fitneptunepos
