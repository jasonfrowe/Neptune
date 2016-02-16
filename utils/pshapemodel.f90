subroutine pshapemodel(npt,x,ixo,ax,nfit,sol,psmod)
use precision
implicit none
integer :: npt,ixo,nfit
real(double), dimension(:) :: x,ax,sol,psmod
!local vars
integer :: i,nxpixel,np1,np2
real(double) :: pixel,poly,Pi,twoPi

Pi=acos(-1.d0)
twoPi=2.0d0*Pi

pixel=poly(x(1),ixo,ax)
nxpixel=1
np1=int(pixel)
psmod(1)=sol(1)*sol(2+1)*sin(twoPi*pixel+sol(2))
do i=2,npt
   pixel=poly(x(i),ixo,ax)
   np2=int(pixel)
   if(np2.ne.np1)then
      nxpixel=nxpixel+1
      np1=np2
      if(nxpixel.gt.nfit-2)then
         write(0,*) "Warning: More pixels than parameters.."
         write(0,*) "nfit-2 : ",nfit-2
         write(0,*) "nxpixel: ",nxpixel
      endif
   endif
   psmod(i)=sol(1)*sol(2+np2)*sin(twoPi*pixel+sol(2))
enddo


return
end
