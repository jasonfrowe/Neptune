subroutine fitter(npt,Kfac,npars,pars,x,y,yerr,xnep,ynep)
use precision
implicit none
integer :: npt,npars,nfitp,i,npord,iplast,nxpix,ipcurrent
real(double), dimension(:) :: pars,x,y,yerr,xnep,ynep
real(double), dimension(:,:) :: Kfac
integer, allocatable, dimension(:) :: isol
real(double), allocatable, dimension(:) :: sol

nfitp=npars !number of parameters that can be fit.
allocate(sol(nfitp),isol(nfitp))
isol=-1 !initiate isol array.  if isol(i).ne.0 then the variable is fit.
do i=1,npars
   sol(i)=pars(i)
enddo

npord=2 !order of polynomial to fit across pixel boundaries.

!figure out number of pixels.
!X-direction pixels - this needs to be done as a function of time
!thus input data should be sorted wrt time.
iplast=int(xnep(1))
nxpix=1
do i=2,npt
   ipcurrent=int(xnep(i))
!   write(0,*) i,iplast,ipcurrent
   if(iplast.ne.ipcurrent)then !did we move a pixel
      nxpix=nxpix+1 !update counter
      iplast=ipcurrent !update pixel location tracker
   endif
enddo
write(0,*) "Number of pixels: ",nxpix

return
end subroutine fitter
