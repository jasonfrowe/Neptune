subroutine fitline(npt,x,y,yerr,ans,eans,chisq)
use precision
implicit none
integer :: npt,i
real(double) :: S,Sx,Sy,Sxx,Sxy,d,chisq
real(double), dimension(:) :: x,y,yerr,ans,eans
real(double), allocatable, dimension(:) :: yerr2

!pre-compute
allocate(yerr2(npt))
yerr2=yerr(1:npt)*yerr(1:npt) !yerr^2

!15.2.4 from Numerical recipes
S=Sum(1.0d0/yerr2(1:npt))
Sx=Sum(x(1:npt)/yerr2(1:npt))
Sy=Sum(y(1:npt)/yerr2(1:npt))
Sxx=Sum(x(1:npt)*x(1:npt)/yerr2(1:npt))
Sxy=Sum(x(1:npt)*y(1:npt)/yerr2(1:npt))

!15.2.6 from Numerical Recipes
d=S*Sxx-Sx*Sx
ans(1)=(Sxx*Sy-Sx*Sxy)/d
ans(2)=(S*Sxy-Sx*Sy)/d

!15.2.9 from Numerical Recipes
eans(1)=sqrt(Sxx/d)  !uncertainties in the fit
eans(2)=sqrt(S/d)

write(0,*) "ans1: ",ans(1),"+/-",eans(1)
write(0,*) "ans2: ",ans(2),"+/-",eans(2)

!estimate chi-squared
chisq=Sum( ((y(1:npt)-ans(1)-ans(2)*x(1:npt))/yerr(1:npt))**2.0d0 )
write(0,*) "Chisq: ",chisq,chisq/dble(npt)

return
end subroutine
