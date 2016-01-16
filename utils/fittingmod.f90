module fittingmod
   use precision
   implicit none
   integer, pointer :: npars2,npord2,ixo2,iyo2,nfitp2
   real(double), dimension(:), pointer :: pars2,x2,y2,yerr2,ax2,ay2
   integer, dimension(:), pointer :: isol2
   real(double), dimension(:), pointer :: sol2
end module fittingmod
