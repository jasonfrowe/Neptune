subroutine fitfft(nh,nfft,amp,dt)
! nh - storage size of amp
! nfft - size of FFT
! amp - amplitudes from FFT
! dt - time scaling
use precision
implicit none
!import vars
integer :: nh,nfft
real(double) :: dt
real(double), dimension(:) :: amp
!local vars
real(double) :: cd2uhz

!convert from c/d to uHz
cd2uhz=1.0d6/86400.0d0

end subroutine fitfft
