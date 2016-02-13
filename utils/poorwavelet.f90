subroutine poorwavelet(ns,trs,frs,nover,dt,nsamp,nsamprate,minamp,      &
   maxamp,bb,scaletype)
use precision
implicit none
integer ns,nsamp,nsamprate,nover,scaletype
real, dimension(:) :: bb
real(double) :: dt,minamp,maxamp
real(double), dimension(:) :: trs,frs
!local vars
integer :: i,j1,j2,j,k,nwave,nfftl,nhl,nsampd2,nss,debug,i1,i2,nplot,ncol,&
   ia1,nbin
integer, allocatable, dimension(:,:) :: ia
real:: r,g,b,pt1,pt2
real(double) :: cd2uhz,f,lminamp,lmaxamp,lamp,f1,f2,dtd2,lfmin,lfmax,   &
   dlf,dnplot,lf1,lf2
real(double), allocatable, dimension(:) :: amp,frsl,meanamp,stdamp
!real(double), allocatable, dimension(:) :: wavelet

interface
   subroutine fftspec(nfft,frs,amp,ns,dt,debug)
      use precision
      implicit none
      integer :: nfft,debug,ns
      real(double) :: dt
      real(double), dimension(:) :: frs,amp
   end subroutine fftspec
end interface

!scaletype -> 0-linearscale,1-logscale

!convert from c/d to uHz
cd2uhz=1.0d6/86400.0d0

!precompute half the sample size
nsampd2=nsamp/2
!precompute half the frequency step size
dtd2=dt/2.0d0

!compute logscale for wavelet plot
lminamp=log10(minamp)
lmaxamp=log10(maxamp)

!nwave is number of FFTs computed to make wavelet
nwave=(ns-nsamprate/2)/nsamprate+1

!oversampled datasize for FFT that is power of 2.
nfftl=2**int(log10(dble(nsamp*nover))/log10(2.0d0)+1.0d0)
allocate(frsl(nfftl))

!number of frequency/amplitudes from fftspec to be returned
nhl=(nfftl/2)+1
allocate(amp(nhl))

write(0,*) "nwave,nhl: ",nwave,nhl

!number of point for plotting
!nplot=nhl-1
nplot=400
allocate(ia(nplot,1))
dnplot=dble(nplot)  !precompute int->dble

!set up plot
ncol=64 !number of colours for display
call pgscr(15,0.0,0.3,0.2)
call pgsci(1)

!calculate frequency range.
if(bb(1).eq.bb(2))then
   f=cd2uhz*dble(2)/(dt*dble(nfftl))
   lfmin=log10(f)
   bb(1)=log10(real(f))
   f=cd2uhz*dble(nhl)/(dt*dble(nfftl))
   lfmax=log10(f)
   bb(2)=log10(real(f))
   dlf=lfmax-lfmin !plot range.
else
   lfmin=dble(bb(1))
   lfmax=dble(bb(2))
   dlf=lfmax-lfmin
endif

bb(3)=trs(1)
bb(4)=trs(ns)
!write(0,*) bb(1),bb(2),bb(3),bb(4)
!read(5,*)

call PGSCLP(0) !turn off clipping
call pgvport(0.15,0.95,-1.8,0.95) !make room around the edges for labels
call pgsci(1)
call pgwindow(bb(1),bb(2),bb(3),bb(4)) !plot scale
call pgbox("BCLNTS",0.0,0,"BCNTS",0.0,0)
call pgptxt((bb(1)+bb(2))/2.0,bb(3)-0.06*(bb(4)-bb(3)),0.0,0.5,         &
   "Frequency (\(0638)Hz)")
call pgptxt(bb(1)-0.10*(bb(2)-bb(1)),(bb(4)+bb(3))/2,90.0,0.5,          &
   "Time (days)")

do i=1,ncol
   call heatlut(i*4-3,r,g,b)
   call heatlut(i*4-2,r,g,b)
   call heatlut(i*4-1,r,g,b)
   call heatlut(i*4  ,r,g,b)
   CALL PGSCR(I+15, R, G, B)
enddo

do i=nsamprate/2,ns+nsamprate/2,nsamprate
   pt1=real(trs(max(1,i-nsamprate/2)))
   pt2=real(trs(min(ns,i+nsamprate/2)))
!   write(0,*) i-nsamprate/2,i+nsamprate/2,pt1,pt2
!   read(5,*)
   ia=0 !initialize ia.
!   write(0,*) "trs: ",trs(i),j
   i1=max(1,i-nsampd2)
   i2=min(ns,i+nsampd2)
   nss=i2-i1+1
   frsl=0.0d0 !zero out array for zero-padding
   frsl(1:nss)=frs(i1:i2)
   !calculate amplitude spectrum
   debug=0
   call fftspec(nfftl,frsl,amp,nss,dt,debug)

!   nbin=100
!   allocate(meanamp(nhl),stdamp(nhl))
!   call fftstats(nhl,amp,meanamp,stdamp,nbin)
!   amp=amp-meanamp
!   deallocate(meanamp,stdamp)

   do k=2,nhl
      lamp=log10(amp(k))
      f=cd2uhz*dble(k)/(dt*dble(nfftl))
      f1=cd2uhz*dble(k-1)/(dt*dble(nfftl))
      f2=cd2uhz*dble(k+1)/(dt*dble(nfftl))
      lf1=log10((f+f1)/2.0d0)
      lf2=log10((f+f2)/2.0d0)
!      write(0,*) "f: ",f,(f1-lfmin)/dlf
      j1=max(1,int((lf1-lfmin)/dlf*dnplot)+1)
      j2=min(nplot,int((lf2-lfmin)/dlf*dnplot)+1)
      if(scaletype.eq.0)then
         ia1=int((amp(k)-minamp)/(maxamp-minamp)*dble(NCOL-1))+16
         if(amp(k).lt.minamp) ia1=16
         if(amp(k).gt.maxamp) ia1=ncol+15
      else
         ia1=int((lamp-lminamp)/(lmaxamp-lminamp)*dble(NCOL-1))+16
         if(lamp.lt.lminamp) ia1=16
         if(lamp.gt.lmaxamp) ia1=ncol+15
      endif
!      write(0,*) j1,j2,ia1
      do j=j1,j2
         ia(j,1)=max(ia(j,1),ia1)
      enddo
!      read(5,*)
   enddo
!   write(0,*) "f: ",f,f1,f2
!   write(0,*) "j: ",j1,j2,nhl

!   call pgpixl(ia,nplot,1,1,nplot,1,1,bb(1),bb(2),real(trs(i1)),real(trs(i2)))
   call pgpixl(ia,nplot,1,1,nplot,1,1,bb(1),bb(2),pt1,pt2)

!   write(0,*) "Wave # ",i,"done"
!   read(5,*)
enddo

call pgbox("BCLNTS",0.0,0,"BCNTS",0.0,0) !replot axes

call PGSCLP(1) !enable clipping

return
end subroutine poorwavelet
