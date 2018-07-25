program pixelfitpca
use precision
implicit none
integer :: iargc,nmax,npt,ixo,iyo,npars,info,nrhs,i,j,k,imaxp,iminp,npix,&
 npixsampmax,nspl,minsamp,npca,ngpflag,npixtype,nplots,ii
integer, allocatable, dimension(:) :: npixsamp
integer, allocatable, dimension(:,:) :: npixmap
real, allocatable, dimension(:) :: bb
real(double) :: minx,mean,pixel,poly,maxp,minp,timediv,px,py,dnsplm1,pxold
real(double), allocatable, dimension(:) :: x,y,yerr,xpos,ypos,xnep,ynep,&
 oflux,pmod,smod,ax,ay,xpcor,ypcor,pars,alpha,yerr2,mu,std,res,phase,xspl,&
 yspl,pixels,xp,yp,zp
real(double), allocatable, dimension(:,:) :: Kernel,KernelZ,xpca,pcavec
character(80) :: filename

!These are F90 interfaces that allow one to pass assumed sized arrays
!to subroutines.
interface !reads in a three-column ascii space seperated file
   subroutine getpdata(filename,npt,nmax,x,y,yerr,oflux,pmod,smod,xpos, &
    ypos,xnep,ynep,minx,mean)
      use precision
      implicit none
      integer, intent(inout) :: npt,nmax
      real(double), intent(inout) :: minx,mean
      real(double), dimension(:), intent(inout) :: x,y,yerr,oflux,pmod, &
         smod,xpos,ypos,xnep,ynep
      character(80), intent(inout) :: filename
   end subroutine getpdata
   subroutine fitneptunepos(npt,x,xnep,ynep,ixo,ax,iyo,ay,nplots)
      use precision
      implicit none
      integer :: npt,ixo,iyo,nplots
      real(double), dimension(:) :: x,xnep,ynep,ax,ay
   end subroutine fitneptunepos
   subroutine makekernel(Kernel,npt1,npt2,x1,x2,npt,yerr,npars,pars)
      use precision
      implicit none
      integer, intent(inout) :: npt1,npt2,npars,npt
      real(double), dimension(:), intent(inout) :: x1,x2,yerr,pars
      real(double), dimension(:,:), intent(inout) :: Kernel
   end subroutine makekernel
   subroutine plotdatascatter(npt,x,y,yerr,bb)
      use precision
      implicit none
      integer, intent(inout) :: npt
      real, dimension(:), intent(inout) :: bb
      real(double), dimension(:), intent(inout) :: x,y,yerr
   end subroutine plotdatascatter
   subroutine plotsamples(npt1,x1,mu,std)
      use precision
      implicit none
      integer, intent(in) :: npt1
      real(double), dimension(:), intent(in) :: x1,mu,std
   end subroutine plotsamples
end interface

!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!Parameters!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!
nplots=0 !0 - make plots, 1 - do not make plots
nmax=80000 !initial guess for number of datapoints.
ixo=5 !order of polynomial to fit to x-positions
iyo=5 !order of polynomial to fit to y-positions
npixtype=1 !0-use native measurements, 1-use polynomial+pointing jitter.
minsamp=50 !minimum number of samples to use a pixel lightcurve
nspl=200 !number of samples for intrapixel model.
!Here are the parameters that control the co-variance matrix and fitted parameters
ngpflag=1 !0 - use GP, 1 do not use GP
npars=4 !number of parameters used for model of the matrix
allocate(pars(npars))  !model parameters for Kernel generation.
pars(1)=1.0d0 !amp scale for exp
pars(2)=0.72d0 !length scale for exp
pars(3)=1.0 !second amp scale
pars(4)=2.0 !second length scale
!!!!!!!!!!!!!!!!!!!!!!!

!check that we have enough information from the commandline
if(iargc().lt.1)then !if not, spit out the Usage info and stop.
   write(0,*) "Usage: pixelfit filename"
   stop
endif

!read in filename containing data (3 columns)
call getarg(1,filename)

allocate(x(nmax),y(nmax),yerr(nmax),oflux(nmax),pmod(nmax),smod(nmax),  &
   xpos(nmax),ypos(nmax),xnep(nmax),ynep(nmax))

call getpdata(filename,npt,nmax,x,y,yerr,oflux,pmod,smod,xpos,ypos,xnep,ynep,minx,mean)
write(0,*) "Number of points read: ",npt !report how much data was read in
write(0,*) "Min x: ", minx
write(0,*) "Mean : ", mean

!open PGPLOT device
if(nplots.eq.0)then
	call pgopen('/xserve')  !'?' lets the user choose the device.
	call PGPAP (6.0 ,1.0) !use a square 8" across
	call pgpage() !create a fresh page
	call pgslw(3) !thicker lines
	call pgask(.true.)
endif

!need fits to the motion of Neptune - raw data is too noisy.
!also means that pixel-crossings is now a linear function (good!)
allocate(ax(ixo),ay(iyo))
allocate(xpcor(npt),ypcor(npt))
xpcor=xnep+xpos !correct for pointing jitter
ypcor=ynep+ypos
call fitNeptunePos(npt,x,xpcor,ypcor,ixo,ax,iyo,ay,nplots)
deallocate(xpcor,ypcor)

if(nplots.eq.0)then
	!plot the data
	allocate(bb(4))
	bb=0.0e0 !rescale plot
	write(0,*) "plotting.."
	call plotdatascatter(npt,x,y,yerr,bb) !plot our original dataset
endif

!GP Smoothing
if(ngpflag.eq.0)then

	!lets make a Kernel/co-variance for the Gaussian process
	allocate(Kernel(npt,npt)) !allocate space
	call makekernel(Kernel,npt,npt,x,x,npt,yerr,npars,pars) !create Kernel

	!Cholesky factorization
	write(0,*) "Beginning Cholesky factorization"
	call dpotrf('U',npt,Kernel,npt,info) !LAPACK routine for Cholesky
	if (info.ne.0) then !check for errors
   		write(0,*) "Cholesky factorization failed"
   		write(0,*) "dpotrf info: ",info
   		stop
	endif

	write(0,*) "Beginning Cholesky Solver"
	!calculate solution
	! We have to solve,
	! Kernel*X=Y,
	! for X.
	! LAPACK routine dpotrs uses Kfac from dpotrf.  We copy y into a new
	! array called alpha, which is overridden with the solution on completion.
	allocate(alpha(npt))
	alpha=y(1:npt) !dpotrs takes alpha as input and output.
	nrhs=1 !how man columns does alpha contain - just one
	call dpotrs('U',npt,nrhs,Kernel,npt,alpha,npt,info) !call LAPACK
	if (info.ne.0) then !check for errors
   	write(0,*) "Solver failed.."
   	write(0,*) "dpotrs info: ",info
   	stop
	endif

	!predict samples
	allocate(KernelZ(npt,npt))
	allocate(yerr2(npt),mu(npt),std(npt))
	yerr2=0.0d0
	call makekernel(KernelZ,npt,npt,x,x,npt,yerr2,npars,pars)
	allocate(res(npt))

!$OMP PARALLEL
!$OMP WORKSHARE
	mu=matmul(KernelZ,alpha)
	std=0.0d0
	res(1:npt)=y(1:npt)-mu(1:npt) !contains 'detrended' data
!$OMP END WORKSHARE
!$OMP END PARALLEL

	if(nplots.eq.0)then
    	write(0,*) "done plotting"
		call plotsamples(npt,x,mu,std) !plot our predicted sample set on top.
	endif

!deallocate Kernel resources
	deallocate(alpha,Kernel,KernelZ,yerr2,mu,std)

!code for skipping GP
else
	allocate(res(npt))
	res=y !residuals are equal to input
endif


allocate(phase(npt),pixels(npt))

if(npixtype.eq.1)then
	!$OMP PARALLEL DO PRIVATE(pixel)
	do i=1,npt
		pixel=poly(x(i),ixo,ax)
   		pixels(i)=pixel-xpos(i)
   	enddo
   	!$OMP END PARALLEL DO
else
	pixels(1:npt)=xnep(1:npt)
endif


maxp=-1000.0
minp=1000.0
do i=1,npt
   pixel=pixels(i)
   if(pixel.gt.maxp) maxp=pixel
   if(pixel.lt.minp) minp=pixel
   !pixel=xnep(i)
   phase(i)=pixel-floor(pixel)
enddo

if(nplots.eq.0)then
	call pgpage()
	bb=0.0
	call plotdatascatter(npt,phase,res,yerr,bb)
endif

write(0,*) "Min/max Pixel: ",minp,maxp
iminp=int(minp+0.5) !round up to first pixel with full coverage
imaxp=int(maxp-0.5) !round down to last pixel with full coverage

npix=(imaxp-iminp)+1
allocate(npixsamp(npix))
npixsamp=0 !initiate counters to zero
timediv=39.1762-minx

!start by getting the maximum number of values for any individual pixel
do i=1,npt
	pixel=pixels(i)
	j=int(pixel)-iminp+1
	if((j.gt.0).and.(j.le.npix))then
		npixsamp(j)=npixsamp(j)+1
	endif
enddo

!now we know the largest array size.
!add 2 for edges to avoid extrapolation
npixsampmax=maxval(npixsamp)+2 
!write(0,*) npixsampmax

deallocate(npixsamp)
allocate(npixsamp(npix*2)) !now we need twice as much
allocate(npixmap(npixsampmax,npix*2))
npixsamp=0 !initiate counters to zero

!fill in using pre-quadrature lightcurve
do i=1,npt
	if(x(i).lt.timediv)then
		pixel=pixels(i)
		j=int(pixel)-iminp+1
		if((j.gt.0).and.(j.le.npix))then
			npixsamp(j)=npixsamp(j)+1
			npixmap(npixsamp(j),j)=i !store location
		endif
	endif
enddo

!fill in using post-quadrature lightcurve
do i=1,npt
	if(x(i).gt.timediv)then
		pixel=pixels(i)
		j=int(pixel)-iminp+1+npix 
		if((j.gt.npix).and.(j.le.npix*2))then
			npixsamp(j)=npixsamp(j)+1
			npixmap(npixsamp(j),j)=i !store location
		endif
	endif
enddo
dnsplm1=dble(nspl-1) !pre-calculate int->dble
allocate(xpca(nspl,npix*2)) !saves resampled pixel lightcurves
allocate(xspl(npixsampmax),yspl(npixsampmax)) !storage for individual pixel LCs

k=0
do j=1,npix*2
	if(npixsamp(j).gt.minsamp)then !only use well sampled pixels.

		ii=0
		if(j.gt.1)then
			if(npixsamp(j-1).gt.0)then
				ii=ii+1
				pxold=pixels(npixmap(1,j-1))
				do i=1,npixsamp(j-1)
					px=pixels(npixmap(i,j-1))
					if(px.ge.pxold)then
						xspl(ii)=px-floor(px)-1.0d0 !x-position
						yspl(ii)=res(npixmap(i,j-1)) !flux
						pxold=px
					endif
				enddo
				!write(0,*) ii,xspl(ii),yspl(ii)
			endif
		endif

		k=k+1
		do i=1,npixsamp(j)
			px=pixels(npixmap(i,j))
			xspl(i+ii)=px-floor(px) !x-position
			yspl(i+ii)=res(npixmap(i,j)) !flux
			!write(0,*) i+ii,xspl(i+ii),yspl(i+ii)
		enddo

		if(j.lt.npix*2)then
			if(npixsamp(j+1).gt.0)then
				ii=ii+1
				pxold=pixels(npixmap(1,j+1))
				do i=1,npixsamp(j+1)
					px=pixels(npixmap(i,j+1))
					if(px.le.pxold)then
						xspl(npixsamp(j)+ii)=px-floor(px)+1.0d0 !x-position
						yspl(npixsamp(j)+ii)=res(npixmap(i,j+1)) !flux
						pxold=px
					endif
				enddo
				!write(0,*) "last"
				!write(0,*) npixsamp(j)+ii,xspl(npixsamp(j)+ii),yspl(npixsamp(j)+ii)
			endif
		endif

		!read(5,*)

		!need to sort xspl,yspl 
		call sort(npixsamp(j)+ii,xspl,yspl)
		!now we can resample.
		do i=1,nspl !resample onto common grid
			px=dble(i-1)/dnsplm1
			call lininterp(xspl,yspl,npixsamp(j)+ii,px,py)
			!write(0,*) px,py
			xpca(i,k)=py !store resampled value into grid.
			!write(0,*) i,k,py
		enddo
	endif
enddo
npca=k !number of lightcurves for PCA
!write(0,*) "k: ",k,npix*2

deallocate(xspl,yspl) !deallocate spline resources.

allocate(pcavec(nspl,npix*2)) !contains PCAs
call pca(nspl,npca,xpca,nspl,npix*2,pcavec) !call PCA routine

!make a plot of the PCA
if(nplots.eq.0)then
	call pgpage()
	bb=0.0
	allocate(xp(nspl),yp(nspl),zp(nspl))
	zp=0.0

	do i=1,nspl
		xp(i)=dble(i-1)/dnsplm1
	enddo

	do i=1,5
		yp=pcavec(:,i)
		call plotdatascatter(nspl,xp,yp,yerr,bb)
	enddo

	do i=1,nspl
		write(6,500) xp(i),(pcavec(i,j),j=1,5)
		500 format(6(F9.6,1X))
	enddo
endif

!fit PCA to data
!call pcafit(nspl,npix*2,pcavec,npixsampmax,npix*2,npixsamp,npixmap,npt,pixels,res)

if(nplots.eq.0) call pgclos()

end program pixelfitpca

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function poly(x,nfit,ans)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: nfit
real(double) :: x,poly
real(double), dimension(nfit) :: ans
integer :: i

poly=ans(1)
do i=2,nfit
   poly=poly+ans(i)*x**dble(i-1)
enddo

return
end function
