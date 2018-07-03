subroutine pca(nspl,npca,xpca,nsplmax,npcamax,y)
!nspl = number of samples for each data array
!npca = nimber of data arrays
!xpca = holds data dim (nsplmax,npcamax)
use precision
implicit none
!import vars
integer :: nspl,npca,nsplmax,npcamax
real(double), dimension(nsplmax,npcamax) :: xpca,y
!local vars
integer :: i,j
real(double) :: dnspl,var
real(double), allocatable, dimension(:) :: m
real(double), allocatable, dimension(:,:) :: C,xpcac
!LAPACK vars
integer :: lwork,info,ldvl,ldvr,n,lda
real(double), allocatable, dimension(:) :: wr,wi,work
real(double), allocatable, dimension(:,:) :: vl,vr

dnspl=dble(nspl) !precalculate int->dble

allocate(m(npca))
!calculate expected values (mean)
do i=1,npca
	m(i)=Sum(xpca(1:nspl,i))/dnspl
	!write(0,*) i,m(i)
enddo

!calculate mean-centered data set
allocate(xpcac(nspl,npca))
do i=1,npca
	xpcac(:,i)=xpca(:,i)-m(i)
enddo

!calculate co-variance matrix
! ** this can be made faster by a factor of ~2 since the matrix is symmetric
allocate(C(npca,npca))
do i=1,npca
	do j=1,npca
		var=Sum(xpcac(1:nspl,i)*xpcac(1:nspl,j))/dnspl
		C(i,j)=var
		!write(0,*) i,j,var
	enddo
enddo

n=npca
lda=n
allocate(wr(n),wi(n))
ldvl=n
ldvr=n
allocate(vl(ldvl,n),vr(ldvr,n))
lwork=4*n
allocate(work(lwork))

lwork=-1
!write(0,*) "Calling DGEEV"
call dgeev(  'V',   'V',   n, C, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
!call DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
lwork=int(work(1)) !get optimum workspace
!write(0,*) "lwork",lwork,4*npca
deallocate(work)    !allocate optimum workspace
allocate(work(lwork)) 
call dgeev(  'V',   'V',   n, C, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)

!calculate PCAs
y(1:nspl,1:npca)=matmul(xpcac,vr)

return
end



