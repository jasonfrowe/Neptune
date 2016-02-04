subroutine spcor(npt,sp,x,r,npixel)
use precision
implicit none
integer :: npt
integer, dimension(:) :: npixel
real(double), dimension(:) :: x,r,sp
!local vars
integer :: i,j,np,ipf,ipl
real(double), allocatable, dimension(:) :: xp,rp,rp2

allocate(xp(npt),rp(npt))

j=0 !counter
ipf=1
do i=1,npt-1
   ipl=i
   if(npixel(i+1).ne.npixel(i))then !jump found
      j=j+1
      xp(j)=(x(ipl)+x(ipf))/2.0d0 !average position
      rp(j)=(r(ipl)+r(ipf))/2.0d0 !average amplitude
      ipf=i+1
      write(0,*) ipf,ipl,xp(j),rp(j)
   endif
enddo
np=j !number of pixel positions
write(0,*) "np: ",np

!set up spline
allocate(rp2(np))
call spline(xp,rp,np,1.d30,1.d30,rp2)

!calculate interpolated spline values
!$OMP PARALLEL DO
do i=1,npt
   call splint(xp,rp,rp2,np,x(i),sp(i))
enddo
!$OMP END PARALLEL DO

return
end
