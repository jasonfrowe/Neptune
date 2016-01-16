subroutine fcn(npt,nfit,solin,fvec,iflag)
use precision
use fittingmod
implicit none
integer :: npt,nfit,iflag
real(double), dimension(npt) :: fvec
real(double), dimension(nfit) :: solin
!local variables
integer i,j,iplast,ipcurrent,nxpix
real(double) :: poly,pixel,fpix
real(double), allocatable, dimension(:) :: sol3,r,pmodelpars

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
!allocate space for polynomial parameters
allocate(pmodelpars(npord2))

!walk through X-positions
iplast=int(poly(x2(1),ixo2,ax2)) !use pointers to reference data
nxpix=1
do j=1,npord2
!   write(0,*) "jj: ",sol3(npars2+j+2*(nxpix-1)),npars2+j+2*(nxpix-1)
   pmodelpars(j)=sol3(npars2+j+2*(nxpix-1))
enddo
do i=1,npt
   pixel=poly(x2(i),ixo2,ax2)
   ipcurrent=int(pixel)
!   write(0,*) i,iplast,ipcurrent
   if(iplast.ne.ipcurrent)then !did we move a pixel
      nxpix=nxpix+1 !update counter
      iplast=ipcurrent !update pixel location tracker
      do j=1,npord2
!         write(0,*) "jj: ",sol3(npars2+j+2*(nxpix-1)),npars2+j+2*(nxpix-1)
         pmodelpars(j)=sol3(npars2+j+2*(nxpix-1))
      enddo
   endif
   fpix=pixel-floor(pixel)
   r(i)=poly(fpix,npord2,pmodelpars)
!   write(0,*) fpix,r(i),pmodelpars(1),pmodelpars(2)
!   read(5,*)
enddo

!Walk through y-positions.

!1b. if Kernel is changing, then we need to update Cholosky stuff.

!2. calculate Kernel-1 x r

!3. break out terms of r x (Kernel-1 x r) into fvec

!4. profit


!until code is ready we have a read statement for debugging.
read(5,*)

return
end subroutine fcn
