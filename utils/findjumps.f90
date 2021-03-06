subroutine findjumps(npt,x,y,yerr,npixel)
use precision
implicit none
!import variables
integer :: npt
integer, dimension(:) :: npixel
real(double), dimension(:) :: x,y,yerr
!local variables
integer i,nline,ndata,ns,nf,j,imaxth,imaxc,currentpixel,nsampmax,i1,i2, &
   nsamp
real(double) :: threshold,th,ap,bp,af,bf,abdevp,abdevf,sig,stdev,mean,  &
   chisq,maxth,th2,meddiff,sig_glob
real(double), allocatable, dimension(:) :: var,xdfit,ydfit,yderr,ans,   &
   eans,tpixels,samps

interface !fits a straight line to data
   subroutine fitline(npt,x,y,yerr,ans,eans,chisq)
      use precision
      implicit none
      integer, intent(in) :: npt
      real(double), intent(inout) :: chisq
      real(double), dimension(:), intent(in) :: x,y,yerr
      real(double), dimension(:), intent(inout) :: ans,eans
   end subroutine fitline
end interface

!global sigma scatter
sig_glob=stdev(npt,y,mean)
!write(0,*) "sig: ",sig
!sig=1.1e-4

!new thresholding vars. 
nsampmax=25 !number of +/- nearby samples to use for stats
threshold=6.0
allocate(samps(nsampmax*2+1))

nline=10 !number of points to use to estimate slope.
allocate(xdfit(nline),ydfit(nline),yderr(nline),ans(2),eans(2),         &
   tpixels(npt))

allocate(var(3))

!do i=1,5
!   npixel(i)=0
!enddo

tpixels=0.0d0 !initalize

do i=5,npt-5

   !estimate scatter by taking a local sample, then estimating the standard deviation  
   i1=max(1,i-nsampmax)
   i2=min(npt,i+nsampmax)
   nsamp=i2-i1+1
   samps(1:nsamp)=y(i1:i2)
   sig=meddiff(nsamp,samps) !calculate the median point-to-point change.
   !write(0,*) i,sig_glob,sig

   var(1)=y(i-1)-y(i-2)
   var(2)=y(i)-y(i-1)
   var(3)=y(i)-y(i+1)

   ns=max(1,i-nline)
   nf=i-1
   ndata=nf-ns+1
!   write(0,*) i,ns,nf,ndata
   xdfit(1:ndata)=x(ns:nf)
   ydfit(1:ndata)=y(ns:nf)
   yderr(1:ndata)=sig
!   call medfit(xdfit,ydfit,ndata,ap,bp,abdevp)
   call fitline(ndata,xdfit,ydfit,yderr,ans,eans,chisq)
   ap=ans(1)
   bp=ans(2)
!   write(0,*) ap,bp,abdevp

   ns=i
   nf=min(npt,i+nline-1)
   ndata=nf-ns+1
!   write(0,*) i,ns,nf,ndata
   xdfit(1:ndata)=x(ns:nf)
   ydfit(1:ndata)=y(ns:nf)
   yderr(1:ndata)=sig
!   call medfit(xdfit,ydfit,ndata,af,bf,abdevf)
   call fitline(ndata,xdfit,ydfit,yderr,ans,eans,chisq)
   af=ans(1)
   bf=ans(2)

   th2=var(2)/sqrt(var(1)*var(1)+var(3)*var(3)+sig*sig)
   th=((ap+bp*x(i))-(af+bf*x(i)))/sig!max(sig,abdevp,abdevf)
!   write(6,'(6(F17.11,1X))') x(i),ap+bp*x(i),af+bf*x(i),th,th2
   if((abs(th).gt.threshold).and.(abs(th2).lt.2.0))then !get rid of more crap.
      tpixels(i)=abs(th2)
   else
      tpixels(i)=abs(th)
   endif



   !if(abs(th).gt.threshold)then
   !   write(0,*) "Jump.."
   !endif


enddo


currentpixel=1
i=1  !counter
imaxc=i !current 'flux jump' position
do while(i.lt.npt)
!   write(0,*) "i :",i
   ns=max(1,i-nline)
   nf=min(npt,i+nline)
   maxth=tpixels(ns)
   imaxth=ns
   do j=ns+1,nf
      if(tpixels(j).gt.maxth)then
         maxth=tpixels(j)
         imaxth=j
      endif
   enddo
   if(maxth.gt.threshold)then
      if(imaxc.eq.imaxth)then !this is the best within 10-pixels
         npixel(ns:imaxth-1)=currentpixel
         currentpixel=currentpixel+1 !update jump count
         npixel(imaxth:nf)=currentpixel
         i=nf+nline+1 !jump ahead
      else
         npixel(ns:imaxth-1)=currentpixel
         imaxc=imaxth
         i=imaxc !jump ahead to flagged pixel and check.
      endif
   else
      npixel(ns:nf)=currentpixel
      i=nf+nline+1 !jump ahead..
   endif
!   write(0,*) "Current pixel: ",currentpixel
enddo
npixel(imaxc+1:npt)=currentpixel !make sure we fill in ends

return
end subroutine findjumps
