subroutine pixelmodelv2(r,npars,npix,npord,sol,npt,x,npixel)
use precision
implicit none
!input vars
integer :: npars,npix,npord,npt
integer, dimension(:) :: npixel
real(double), dimension(:) :: sol,x,r
!local vars
integer :: i,j,iplast,ipcurrent
real(double) :: zpt,fpix,poly
real(double), allocatable, dimension(:) :: pmodelpars

allocate(pmodelpars(npord))

iplast=0!npixel(1)
zpt=x(1) !start of segment
do i=1,npt
   ipcurrent=npixel(i)
   if(ipcurrent.ne.iplast)then !did we move to a new segment?
      zpt=x(i)  !update zero point
      do j=1,npord !update polynomial coefficients
!         write(0,*) npars+j+2*(ipcurrent-1),sol(npars+j+2*(ipcurrent-1))
         pmodelpars(j)=sol(npars+j+2*(ipcurrent-1))
      enddo
      iplast=ipcurrent
   endif
   fpix=x(i)-zpt
   r(i)=poly(fpix,npord,pmodelpars)
!   write(0,*) fpix,r(i),(pmodelpars(j),j=1,npord)
!   read(5,*)
enddo

return
end
