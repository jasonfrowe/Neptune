subroutine pixelmodel(r,npars2,npord2,sol3,npt,x2,ixo2,ax2,iyo2,ay2)
!jasonfrowe@gmail.com 2016
use precision
implicit none
!passed parameters
integer :: npord2,npt,ixo2,iyo2,npars2
real(double), dimension(:) :: r,sol3,ax2,ay2,x2
!local parameters
integer i,j,iplast,ipcurrent,nxpix
real(double) :: poly,pixel,fpix
real(double), allocatable, dimension(:) :: pmodelpars

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
   r(i)=r(i)+poly(fpix,npord2,pmodelpars)
!   write(0,*) fpix,r(i),pmodelpars(1),pmodelpars(2)
!   read(5,*)
enddo

!Walk through y-positions.


end subroutine pixelmodel
