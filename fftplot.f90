program fftplot
implicit none
integer :: nmax,npt,nptn
real, allocatable, dimension(:) :: freq,amp,bb,px,py,freqn,ampn
character(80) :: fitfile,normfile

nmax=300000
allocate(freq(nmax),amp(nmax))
fitfile="fittest.20160428.dat"
call readpow(fitfile,nmax,npt,freq,amp)

!open PGPLOT device
call pgopen('?')!('/xserve')  !'?' lets the user choose the device.
call PGPAP (8.0 ,0.5) !use a square 8" across
call pgsubp(1,1)
call pgpage() !create a fresh page
call pgslw(3) !thicker lines
call pgsch(1.5) !bigger text

allocate(bb(4))
bb(1)=minval(freq(1:npt))
bb(2)=maxval(freq(1:npt))
bb(3)=0.0
bb(4)=8.0

call pgvport(0.15,0.95,0.15,0.95) !make room around the edges for labels
call pgwindow(bb(1),bb(2),bb(3),bb(4)) !plot scale
call pgbox("BCNTS1",0.0,0,"BCNTS",0.0,0)
call pglabel("Frequency (\(0638)Hz)","Normalized Amplitude","")

call pgline(npt,freq,amp)

allocate(px(2),py(2))

px(1)=bb(1)
px(2)=bb(2)
py(1)=sqrt(9.89)
py(2)=sqrt(9.89)
call pgsci(3)
call pgline(2,px,py)
call pgsci(1)

normfile="norm2_avg.dat"
allocate(freqn(nmax),ampn(nmax))
call readnorm(normfile,nmax,nptn,freqn,ampn)

call pgsci(2)
call pgline(nptn,freqn,ampn)
call pgsci(1)

call pgclos()



end program fftplot

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine readnorm(filename,nmax,npt,freqn,ampn)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
implicit none
!import vars
integer :: nmax,npt
real, dimension(nmax) :: freqn,ampn
character(80) :: filename
!local vars
integer :: nunit,filestatus,i

nunit=10
open(unit=nunit,file=filename,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",filename
   stop
endif

i=1
do
   if(i.gt.nmax)then
      write(0,*) "Increase nmax to match data points"
      write(0,*) "nmax: ",nmax
      stop
   endif
   read(nunit,*,iostat=filestatus) freqn(i),ampn(i)
   ampn(i)=sqrt(ampn(i))
   if(filestatus == 0) then
      i=i+1
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit) !close file
npt=i-1

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine readpow(filename,nmax,npt,freq,amp)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
implicit none
!import vars
integer :: nmax,npt
real, dimension(nmax) :: freq,amp
character(80) :: filename
!local vars
integer :: nunit,filestatus,i
real :: pd1,pd2,dumr

nunit=10
open(unit=nunit,file=filename,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",filename
   stop
endif

i=1
do
   if(i.gt.nmax)then
      write(0,*) "Increase nmax to match data points"
      write(0,*) "nmax: ",nmax
      stop
   endif
   read(nunit,*,iostat=filestatus) freq(i),pd1,pd2
   amp(i)=sqrt(pd1/pd2)
   if(filestatus == 0) then
      i=i+1
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit) !close file
npt=i-1

return
end
