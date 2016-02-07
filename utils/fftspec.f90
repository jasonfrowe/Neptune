subroutine fftspec(nfft,frs,amp,npt,dt,debug)
use precision
!iso_c_binding is for FFTW3 interface
use, intrinsic :: iso_c_binding
implicit none
!add in the FFTW3 modules
include 'fftw3.f03'
!import vars
integer :: nfft,debug,npt
real(double) :: dt
real(double), dimension(:) :: frs,amp
!local vars
integer :: i
real(double) :: cd2uhz,amp1

!FFTW3 vars
type(C_PTR) :: plan
integer ( kind = 4 ) :: nh
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: frsC

!convert from c/d to uHz
cd2uhz=1.0d6/86400.0d0

!space needed for FFT output
nh=(nfft/2)+1
allocate(frsC(nh))

!generate FFTW plan
plan=fftw_plan_dft_r2c_1d(nfft,frs,frsC,FFTW_ESTIMATE)
!execute FFW
call fftw_execute_dft_r2c(plan,frs,frsC)
!destroy plan
call fftw_destroy_plan(plan)

!calculate amplitudes
do i=1,nh
   amp(i)=abs(frsC(i))/dble(npt/2)
enddo

!if debug=1, then write out the FFT to 'fft.dat'
if(debug.eq.1)then
   open(unit=11,file="fft.dat")
   do i=1,nh
      amp1=abs(frsC(i))/dble(npt/2)
      write(11,*) cd2uhz*dble(i-1)/(dt*dble(nfft)),amp1
   enddo
   close(11)
endif

return
end subroutine fftspec
