subroutine fitfft(nh,nfft,amp,dt,nover,meanamp,stdamp)
! nh - storage size of amp
! nfft - size of FFT
! amp - amplitudes from FFT
! dt - time scaling
use precision
use fittingfft
implicit none
!import vars
integer,target :: nh,nfft,nover
real(double),target :: dt
real(double), dimension(:), target :: amp,meanamp,stdamp
!local vars
integer :: npt,lwa,info,i,nbin,j,nfitin,nbins,bin,bsum
integer, target :: nfit
integer, allocatable, dimension(:) :: iwa,pbin
integer, allocatable, dimension(:), target :: isol
real(double) :: cd2uhz,tollm,dnfft,pmodel,ftemp,ptemp,bmin,bmax,thres,  &
   dN,Pm,opn
real(double), allocatable, dimension(:), target :: p,f,sol
real(double), allocatable, dimension(:) :: fvec,wa,solt1,solt2,     &
   solt3,solin
external fcn

!convert from c/d to uHz
cd2uhz=1.0d6/86400.0d0

!real frequency threshold
thres=9.89
dn=dble(nh)/dble(nover)
Pm=1.0d0/dn
opn=1.0d0/(dn*2.8d0)
thres=-log10(1.0d0-(1.0d0-Pm)**opn)
write(0,*) "thres: ",thres

!number of paramters to fit
nfit=13

!initial guess at best fit parameters
allocate(sol(nfit),isol(nfit))
sol=0.0d0
 sol(1)=1.96e-6
isol(1)=1
 sol(2)=9.01e6  !A
isol(2)=1
 sol(3)=1.69e1  !th
isol(3)=1
 sol(4)=2.27  !alpha
isol(4)=1

 sol(5)=9.97e4  !A
isol(5)=1
 sol(6)=8.30e-2 !th
isol(6)=1
 sol(7)=4.38  !alpha
isol(7)=1

 sol(8)=2.0d4 !amplitude
isol(8)=1
 sol(9)=1.0 !vmax
isol(9)=1
 sol(10)=25.0 !sigma (width)
isol(10)=1

 sol(11)=1.0d-3 !amplitude
isol(11)=1
 sol(12)=3000.0 !vmax
isol(12)=1
 sol(13)=100.0 !sigma (width)
isol(13)=1


allocate(solin(nfit))
nfitin=0
do i=1,nfit
   if(isol(i).ne.0)then
      nfitin=nfitin+1
      solin(nfitin)=sol(i)
   endif
enddo


!maximum number of datapoints
npt=nh-1
dnfft=dble(nfft)
allocate(f(nh-1),p(nh-1))
nbin=nh/int(cd2uhz*dble(nh-1)/(dt*dnfft)/4.0)
write(0,*) "mf:",cd2uhz*dble(nh-1)/(dt*dnfft),nbin
npt=0
do i=1+nbin,nh,nbin
   npt=npt+1
   f(npt)=0.0d0
   p(npt)=0.0d0
!   write(0,*) i-nbin+1,i
!   read(5,*)
   do j=i-nbin+1,i
      !calculate frequency
      ftemp=cd2uhz*dble(j-1)/(dt*dnfft)
      !calculate power
      ptemp=1.0d12*amp(j)*amp(j)/ftemp
      f(npt)=f(npt)+ftemp
      p(npt)=p(npt)+ptemp
   enddo
   f(npt)=f(npt)/dble(nbin)
   p(npt)=p(npt)/dble(nbin)
enddo
!do i=2,nh
!   !calculate frequency
!   f(i-1)=cd2uhz*dble(i-1)/(dt*dnfft)
!   !calculate power
!   p(i-1)=1.0d12*amp(i)*amp(i)/f(i-1)
!enddo

!update pointers to pass info to fcn subroutine
p2 => p
f2 => f
sol2 => sol
isol2 => isol
nfit2 => nfit

!fvec contains model calculated with sol
allocate(fvec(npt))
!work arrays for lmdif1
lwa=npt*nfitin+5*nfitin+npt
allocate(wa(lwa),iwa(nfitin))

write(0,*) "Starting lmdif1"
tollm=1.0d-12 !tolerence for fitting
!    lmdif1(fcn,m,  n,   x,    fvec,tol,  info,iwa,wa,lwa)
call lmdif1(fcn,npt,nfitin,solin,fvec,tollm,info,iwa,wa,lwa)
write(0,*) "lmdif1 info: ",info

j=0
do i=1,nfit
   if(isol(i).ne.0)then
      j=j+1
      sol(i)=solin(j)
   endif
enddo

allocate(solt1(nfit),solt2(nfit),solt3(nfit))
solt1=sol
solt2=sol
solt3=sol
solt1(11)=0.0d0
solt2(1)=0.0d0
solt2(2)=0.0d0
solt2(8)=0.0d0
solt2(11)=0.0d0
solt3(1)=0.0d0
solt3(5)=0.0d0
solt3(8)=0.0d0
solt3(11)=0.0d0
write(0,*) "sol:"
write(0,'(15(1PE17.10,1X))') sol
open(unit=11,file="fittest.dat")
do i=2,nh
   ftemp=cd2uhz*dble(i-1)/(dt*dnfft)
   ptemp=1.0d12*amp(i)*amp(i)/ftemp
   write(11,'(10(1PE17.10,1X))') ftemp,ptemp,pmodel(nfit,sol,ftemp),    &
    pmodel(nfit,solt1,ftemp),pmodel(nfit,solt2,ftemp),pmodel(nfit,solt3,ftemp)
   meanamp(i)=sqrt(ftemp*pmodel(nfit,sol,ftemp)/1.0d12)
   stdamp(i)=sqrt(ftemp*pmodel(nfit,sol,ftemp)*thres/1.0d12)-meanamp(i)
enddo
close(11)
meanamp(1)=meanamp(2)
stdamp(1)=stdamp(2)


!binning
nbins=100  !number of bins
bmin=0.0  !minimum range for binning
bmax=10.0 !maximum range for binning
allocate(pbin(nbins))
pbin=0 !initialize bin counts to zero
do i=2,nh
   ftemp=cd2uhz*dble(i-1)/(dt*dnfft)
   if(ftemp.gt.100.0)then
      ptemp=1.0d12*amp(i)*amp(i)/ftemp/pmodel(nfit,sol,ftemp) !scale p by model
      bin=int(dble(nbins)*(ptemp-bmin)/(bmax-bmin))+1
      if((bin.gt.0).and.(bin.le.nbins))then
         pbin(bin)=pbin(bin)+1
      endif
   endif
enddo
bsum=sum(pbin(1:nbins))
open(unit=11,file="pbin.dat")
do i=1,nbins
!   write(0,*) bmin+(dble(i)+0.5)/dble(nbins)*(bmax-bmin),pbin(i)
!   write(11,'(F7.3,1X,F9.6)') bmin+(dble(i)+0.5)/dble(nbins)*(bmax-bmin),dble(pbin(i))/dble(bsum)
   write(11,'(F7.3,1X,I6)') bmin+(dble(i)+0.5)/dble(nbins)*(bmax-bmin),pbin(i)
enddo
close(11)


return
end subroutine fitfft

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine fcn(npt,nfitin,solin,fvec,iflag)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
use fittingfft
!use fittingmod
implicit none
integer :: npt,nfitin,iflag
real(double), dimension(npt) :: fvec
real(double), dimension(nfitin) :: solin
!local vars
integer :: i,j
real(double) :: cd2uhz,pm,pmodel
real(double), allocatable, dimension(:) :: sol

!convert from c/d to uHz
cd2uhz=1.0d6/86400.0d0

allocate(sol(nfit2))
sol=sol2

j=0
do i=1,nfit2
   if(isol2(i).ne.0)then
      j=j+1
      sol(i)=solin(j)
   endif
enddo

!write(0,*) "solin:"
!write(0,'(10(1PE17.10,1X))') sol

do i=1,npt
   !calculate power models
   pm=pmodel(nfit2,sol,f2(i))
!   fvec(i)=(pm-p2(i))!/pm
   fvec(i)=log10(pm)-log10(p2(i))
!   write(0,*) i,p2(i),pm
!   read(5,*)
enddo

return
end subroutine

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function pmodel(nfit,sol,f)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
integer :: nfit
real(double) :: pmodel,f
real(double), dimension(nfit) :: sol
!local vars
integer :: i,nc,inc
real(double) :: A,th,al,D,vmax,sig

pmodel=abs(sol(1))+1.0e-12 !white noise floor

nc=2! (nfit-1)/3 !number of components

do i=1,nc
   inc=3*(i-1)+1
   A=abs(sol(inc+1))
   th=abs(sol(inc+2))
   al=abs(sol(inc+3))
!   write(0,*) A,th,al
   pmodel=pmodel+A/(1.0d0+(th*f)**al)
enddo

!first Gaussian
D=abs(sol(8))
vmax=abs(sol(9))
sig=abs(sol(10))
pmodel=pmodel+D*exp(-(vmax-f)**2.0/(2.0d0*sig*sig))

!second Gaussian
D=abs(sol(11))
vmax=abs(sol(12))
sig=abs(sol(13))
pmodel=pmodel+D*exp(-(vmax-f)**2.0/(2.0d0*sig*sig))


return
end
