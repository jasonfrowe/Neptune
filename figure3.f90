program figure3
implicit none
integer :: nmax,npt1,npt2,npt3
real :: medval
real, allocatable, dimension(:) :: time1,phot1,time2,phot2,time3,phot3, &
   bb
character(80) :: filename

nmax=100000
medval=7.885898e6

allocate(time1(nmax),phot1(nmax))
filename="phot.m2.20150723.c.dat"
call readphot(filename,nmax,npt1,time1,phot1)
write(0,*) "npt1: ",npt1
!scale data
phot1=phot1/medval-1.0

allocate(time2(nmax),phot2(nmax))
filename="phot.m2.20160217.cjpm.dat"
call readphot(filename,nmax,npt2,time2,phot2)
write(0,*) "npt2: ",npt2

allocate(time3(nmax),phot3(nmax))
filename="phot.m2.20160217.cjpmd50.dat"
call readphot(filename,nmax,npt3,time3,phot3)
write(0,*) "npt3: ",npt3


!open PGPLOT device
call pgopen('?')!('/xserve')  !'?' lets the user choose the device.
call PGPAP (8.0 ,1.0) !use a square 8" across
call pgsubp(1,2)
call pgpage() !create a fresh page
call pgslw(3) !thicker lines
call pgsch(1.5) !bigger text

allocate(bb(4))
bb(1)=minval(time1(1:npt1))
bb(2)=maxval(time1(1:npt1))
bb(3)=-0.04
bb(4)= 0.16

call pgvport(0.15,0.95,0.15,0.95) !make room around the edges for labels
call pgwindow(bb(1),bb(2),bb(3),bb(4)) !plot scale
call pgbox("BCNTS1",0.0,0,"BCNTS",0.0,0)
call pglabel("Time (days)","Normalized Flux","")

phot1=phot1+0.11
call pgpt(npt1,time1,phot1,-1)

phot2=phot2+0.06
call pgpt(npt2,time2,phot2,-1)

phot3=phot3+0.00
call pgpt(npt3,time3,phot3,-1)

call pgpage()
bb(1)=minval(time1(1:npt1))
bb(2)=22.0
bb(3)=-0.04
bb(4)= 0.16

call pgvport(0.15,0.95,0.15,0.95) !make room around the edges for labels
call pgwindow(bb(1),bb(2),bb(3),bb(4)) !plot scale
call pgbox("BCNTS1",0.0,0,"BCNTS",0.0,0)
call pglabel("Time (days)","Normalized Flux","")

call pgpt(npt1,time1,phot1,-1)

call pgpt(npt2,time2,phot2,-1)

call pgpt(npt3,time3,phot3,-1)

call pgclos()

end program figure3

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine readphot(filename,nmax,npt,time,phot)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
implicit none
!import vars
integer :: nmax,npt
real, dimension(nmax) :: time,phot
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
   read(nunit,*,iostat=filestatus) time(i),phot(i)
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
