module fittingmodv2
   use precision
   implicit none
   integer, pointer :: npt2,nfitp2,npars2,npix2,npord2
   real(double), pointer :: logDK2
   real(double), dimension(:), pointer :: pars2,x2,y2
   integer, dimension(:), pointer :: isol2, npixel2
   real(double), dimension(:), pointer :: sol2
   real(double), dimension(:,:), pointer :: Kernel2
end module fittingmodv2
