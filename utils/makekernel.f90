subroutine makekernel(Kernel,npt1,npt2,x1,x2,npt,yerr,npars,pars)
use precision
implicit none
integer :: npt1,npt2,npars,i,j,npt
real(double) :: twol2,a2,a2l,twol2l
real(double), dimension(:) :: x1,x2,yerr,pars
real(double), dimension(:,:) :: Kernel

!initialize Kernel to zero
Kernel=0.0d0

!pre-compute
a2=pars(1)*pars(1)
twol2=2.0d0*pars(2)*pars(2)
if(twol2.le.0.0d0)then
   write(0,*) "Critial Error: pars(2) must be greater than zero!!"
   stop
endif
a2l=pars(3)*pars(3)
twol2l=pars(4)*pars(4)

do i=1,npt1
   do j=1,npt2
      !exponential part
      Kernel(i,j)=Kernel(i,j)+a2*exp(-(x1(i)-x2(j))*(x1(i)-x2(j))/twol2)
      !linear part
      Kernel(i,j)=Kernel(i,j)+a2l*exp(-(x1(i)-x2(j))*(x1(i)-x2(j))/twol2l)
      !add in co-variance
      if(abs(x1(i)-x2(j)).le.1.0d-15)then
         if((i.gt.npt).or.(j.gt.npt))then
            Kernel(i,j)=Kernel(i,j)+1.0e-15
         else
            Kernel(i,j)=Kernel(i,j)+yerr(i)*yerr(j)
         endif
      endif
   enddo
enddo

write(0,*) "min max",minval(Kernel(1:npt1,1:npt2)),maxval(Kernel(1:npt1,1:npt2))

return
end subroutine makekernel
