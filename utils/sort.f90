subroutine sort(n,x,y)
!simple routine that uses rqsort to sort data by x in-order 
!x,y and returned as sorted.
use precision
implicit none
!import vars
integer :: n
real(double), dimension(n) :: x,y
!local vars
integer :: i
integer, allocatable, dimension(:) :: p
real(double), allocatable, dimension(:) :: px,py

allocate(px(n),py(n),p(n)) !arrays for temp storage and sorting
px=x
py=y
call rqsort(n,x,p)
do i=1,n
	x(i)=px(p(i))
	y(i)=py(p(i))
enddo

return
end