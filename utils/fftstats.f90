subroutine fftstats(nh,amp,meanamp,stdamp,nbin)
use precision
implicit none
integer :: nh,nbin
real(double), dimension(nh) :: amp,meanamp,stdamp
!local vars
integer :: i,j,k,i1,i2,nbind2,nsamp,ncut,iter,itermax,ilast
integer, allocatable, dimension(:) :: cut
real(double) :: stdev,sigcut,std,mean,lmean,lstd
real(double), allocatable, dimension(:) :: samples

nbind2=nbin/2
sigcut=4.0d0

allocate(samples(2*nbind2+1),cut(2*nbind2+1))
ilast=0
do i=1,nh+nbind2,nbin
   i1=max(1,i-nbind2)
   i2=min(nh,i+nbind2)
   nsamp=i2-i1+1
   samples(1:nsamp)=amp(i1:i2)
   stdamp(i)=stdev(nsamp,samples,meanamp(i))

   ncut=1
   iter=0
   itermax=5
   do while ((ncut.gt.0).and.(iter.lt.itermax))

      !cut outliers
      do k=1,nsamp
         if(abs(samples(k)-meanamp(i)).lt.sigcut*stdamp(i))then
            cut(k)=0
         else
            cut(k)=1
         endif
      enddo
      j=0
      do k=1,nsamp
         if(cut(k).eq.0)then
            j=j+1
            samples(j)=samples(k)
         endif
      enddo
      ncut=nsamp-j
!      write(0,*) "ncut: ",ncut
      nsamp=j
      std=stdev(nsamp,samples,mean)

      meanamp(i)=mean
      stdamp(i)=std

      iter=iter+1
   enddo

   if(ilast.gt.0)then
      do k=ilast+1,i
         meanamp(k)=lmean+(mean-lmean)*dble(k-ilast)/dble(i-ilast)
         stdamp(k) = lstd+( std- lstd)*dble(k-ilast)/dble(i-ilast)
      enddo
   else
      meanamp(i1:i)=mean
      stdamp(i1:i)=std
   endif
   lmean=mean
   lstd=std
   ilast=i
enddo
meanamp(ilast:nh)=mean !fill in end
stdamp(ilast:nh)=std
deallocate(samples)


return
end subroutine

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine plotpstats(nh,meanamp,stdamp,nfft,dt)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: nh,nfft
real(double) :: dt
real(double), dimension(nh) :: meanamp,stdamp
!local vars
integer nplot,i
real, allocatable, dimension(:) :: x,y,y2
real(double) :: cd2uhz,f,p,p2,sig

!convert from c/d to uHz
cd2uhz=1.0d6/86400.0d0

!noise threshold
sig=1.0

nplot=nh-1
!PGPlot only accepts real
allocate(x(nplot),y(nplot),y2(nplot))

do i=2,nh

   f=cd2uhz*dble(i-1)/(dt*dble(nfft))
   p=(1.0d6*meanamp(i-1))**2.0/f
   p2=(1.0d6*(meanamp(i-1)+sig*stdamp(i-1)))**2.0/f

   x(i-1) =log10(real(f))
   y(i-1) =log10(real(p))
   y2(i-1)=log10(real(p2))

enddo

call pgsci(2)
call pgline(nplot,x,y)
call pgsci(3)
call pgline(nplot,x,y2)
call pgsci(1)

return
end subroutine

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine plotstats(nh,meanamp,stdamp,nfft,dt)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: nh,nfft
real(double) :: dt
real(double), dimension(nh) :: meanamp,stdamp
!local vars
integer nplot,i
real, allocatable, dimension(:) :: x,y,y2
real(double) :: cd2uhz,f,p,p2,sig

!convert from c/d to uHz
cd2uhz=1.0d6/86400.0d0

!noise threshold
sig=1.0

nplot=nh-1
!PGPlot only accepts real
allocate(x(nplot),y(nplot),y2(nplot))

do i=2,nh

   f=cd2uhz*dble(i-1)/(dt*dble(nfft))
   p=1.0d6*meanamp(i-1)
   p2=1.0d6*(meanamp(i-1)+sig*stdamp(i-1))

   x(i-1) =log10(real(f))
   y(i-1) =real(p)
   y2(i-1)=real(p2)

enddo

call pgsci(2)
call pgline(nplot,x,y)
call pgsci(3)
call pgline(nplot,x,y2)
call pgsci(1)

return
end subroutine
