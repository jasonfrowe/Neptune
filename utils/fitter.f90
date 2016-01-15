subroutine fitter(npt,Kfac,npars,pars,x,y,yerr,xnep,ynep,ixo,ax,iyo,ay, &
   npord)
use precision
use fittingmod
implicit none
integer :: npt,i,j,iplast,nxpix,ipcurrent,nypix
real(double) :: poly,tollm
real(double), dimension(:) :: ax,ay
real(double), dimension(:,:) :: Kfac
!targets for pointer reference to feed fcn
integer, target :: npars,npord,ixo,iyo,nfitp
real(double), dimension(:), target :: pars,x,y,yerr,xnep,ynep
integer, allocatable, dimension(:), target :: isol
real(double), allocatable, dimension(:), target :: sol
!lmdif1 variables
integer :: nfit,info,lwa
integer, allocatable, dimension(:) :: iwa
real(double), allocatable, dimension(:) :: solin,fvec,wa
external fcn

!figure out number of pixels.
!X-direction pixels - this needs to be done as a function of time
!thus input data should be sorted wrt time.
iplast=int(poly(x(1),ixo,ax))
nxpix=1
do i=2,npt
   ipcurrent=int(poly(x(i),ixo,ax))
!   write(0,*) i,iplast,ipcurrent
   if(iplast.ne.ipcurrent)then !did we move a pixel
      nxpix=nxpix+1 !update counter
      iplast=ipcurrent !update pixel location tracker
   endif
enddo
write(0,*) "Number of pixels: ",nxpix

!Y-direction pixels - this needs to be done as a function of time
!thus input data should be sorted wrt time.
iplast=int(poly(x(1),iyo,ay))
nypix=1
do i=2,npt
   ipcurrent=int(poly(x(i),iyo,ay))
!   write(0,*) i,iplast,ipcurrent
   if(iplast.ne.ipcurrent)then !did we move a pixel
      nypix=nypix+1 !update counter
      iplast=ipcurrent !update pixel location tracker
   endif
enddo
write(0,*) "Number of pixels: ",nypix

!assemble model parameters
nfitp=npars+npord*(nxpix+nypix) !number of parameters that can be fit.
allocate(sol(nfitp),isol(nfitp))
isol=-1 !initiate isol array.  if isol(i).ne.0 then the variable is fit.
do i=1,npars
   sol(i)=pars(i)
   isol(i)=0 !keep Kernel constant.
enddo
!add in X-pixel model fit.
do i=1+npars,npars+nxpix*npord
   sol(i)=0.0d0 !start with a straight line for a guess
enddo
!add in Y-pixel model fit
do i=1+npars+nxpix*npord,npars+nxpix*npord+nypix*npord
   sol(i)=0.0d0 !start with a straight line for a guess
enddo

!calculate how many variables are fit.
nfit=0
do i=1,nfitp
   if(isol(i).ne.0)then
      nfit=nfit+1
   endif
enddo
!allocate space for fitted variables (solin) to feed to lmdif1 and fill
!in values fo sol.
allocate(solin(nfit))
j=0
do i=1,nfitp
   if(isol(i).ne.0)then
      j=j+1
      solin(j)=sol(i)
   endif
enddo

!update pointers to share data with fcn
npars2 => npars
npord2 => npord
ixo2 => ixo
iyo2 => iyo
pars2 => pars
x2 => x
y2 => y
yerr2 => yerr
xnep2 => xnep
ynep2 => ynep
isol2 => isol
sol2 => sol
nfitp2 => nfitp

!fvec contains model calculated with solin
allocate(fvec(npt))
!work arrays for lmdif1
lwa=npt*nfit+5*nfit+npt
allocate(wa(lwa))

tollm=1.0d-12 !tolerence for fitting
call lmdif1(fcn,npt,nfit,solin,fvec,tollm,info,iwa,wa,lwa)
write(0,*) "lmdif1 info: ",info

return
end subroutine fitter

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
