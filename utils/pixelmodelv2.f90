subroutine pixelmodelv2(r,npars,npix,npord,sol,npt,x,npixel)
use precision
implicit none
!input vars
integer :: npars,npix,npord,npt
integer, dimension(:) :: npixel
real(double), dimension(:) :: sol,x,r
!local vars
integer :: i,j,iplast,ipcurrent,poly
real(double) :: zpt
real(double), allocatable, dimension(:) :: pmodelpars

allocate(pmodelpars(npord))

iplast=npixel(1)
zpt=x(1) !start of segment
do i=1,npt
   ipcurrent=npixel(i)
   if(ipcurrent.ne.iplast)then !did we move to a new segment?
      zpt=x(i)  !update zero point
      do j=1,npord !update polynomial coefficients
         pmodelpars(j)=sol(npars+npix*npord+j+2*(ipcurrent-1))
      enddo
      iplast=ipcurrent
   endif
   fpix=x(i)-zpt
   r(i)=poly(fpix,npord,pmodelpars)
enddo

return
end
