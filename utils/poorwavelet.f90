subroutine poorwavelet(ns,trs,frs,nsamp,nsamprate)
use precision
implicit none
integer ns,nsamp,nsamprate
real(double), dimension(:) :: trs,frs
!local vars
integer :: i,j,nwave

interface
   subroutine fftspec(nfft,frs,amp,npt,dt,debug)
      use precision
      implicit none
      integer :: nfft,debug,npt
      real(double) :: dt
      real(double), dimension(:) :: frs,amp
   end subroutine fftspec
end interface

!nwave is number of FFTs computed to make wavelet
nwave=(ns-nsamprate/2)/nsamprate+1


j=0
do i=nsamprate/2,ns,nsamprate
   j=j+1
   write(0,*) i,j

enddo
write(0,*) nsamprate,ns


return
end subroutine poorwavelet
