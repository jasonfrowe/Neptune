subroutine fitter(npt,Kfac,npars,pars)
use precision
implicit none
integer :: npt,npars,nfitp,i
real(double), dimension(:) :: pars
real(double), dimension(:,:) :: Kfac
integer, allocatable, dimension(:) :: isol
real(double), allocatable, dimension(:) :: sol

nfitp=npars !number of parameters that can be fit.
allocate(sol(nfitp),isol(nfitp))
isol=-1 !initiate isol array.  if isol(i).ne.0 then the variable is fit.
do i=1,npars
   sol(i)=pars(i)
enddo


return
end subroutine fitter
