program joinparts
use precision
implicit none
integer :: npt1,npt2,nmax,i,j,k,npt
real(double) :: t1,t2,f1,f2
real(double), allocatable, dimension(:) :: x1,y1,yerr1,xpos1,ypos1,     &
   xnep1,ynep1,x2,y2,yerr2,xpos2,ypos2,xnep2,ynep2,x,y,yerr,xpos,ypos,  &
   xnep,ynep,oflux1,pmod1,smod1,oflux2,pmod2,smod2,oflux,pmod,smod
character(80) :: filename1,filename2

!These are F90 interfaces that allow one to pass assumed sized arrays
!to subroutines.
interface !reads in a three-column ascii space seperated file
   subroutine getpdata(filename,npt,nmax,x,y,yerr,oflux,pmod,smod,xpos, &
    ypos,xnep,ynep)
      use precision
      implicit none
      integer, intent(inout) :: npt,nmax
      real(double), dimension(:), intent(inout) :: x,y,yerr,oflux,pmod, &
        smod,xpos,ypos,xnep,ynep
      character(80), intent(inout) :: filename
   end subroutine getpdata
end interface

nmax=80000

allocate(x(nmax),y(nmax),yerr(nmax),xpos(nmax),ypos(nmax),xnep(nmax),   &
   ynep(nmax),oflux(nmax),pmod(nmax),smod(nmax))
!filename="phot.m2.20150723.c.dat"
!call getpdata(filename,npt,nmax,x,y,yerr,xpos,ypos,xnep,ynep,minx,mean)
!write(0,*) "Number of points read: ",npt !report how much data was read in

filename1="phot.p1.j.dat"
allocate(x1(nmax),y1(nmax),yerr1(nmax),xpos1(nmax),ypos1(nmax),         &
  xnep1(nmax),ynep1(nmax),oflux1(nmax),pmod1(nmax),smod1(nmax))
call getpdata(filename1,npt1,nmax,x1,y1,yerr1,oflux1,pmod1,smod1,xpos1, &
  ypos1,xnep1,ynep1)
write(0,*) "Number of points read: ",npt1

filename2="phot.p2.j.dat"
allocate(x2(nmax),y2(nmax),yerr2(nmax),xpos2(nmax),ypos2(nmax),         &
  xnep2(nmax),ynep2(nmax),oflux2(nmax),pmod2(nmax),smod2(nmax))
call getpdata(filename2,npt2,nmax,x2,y2,yerr2,oflux2,pmod2,smod2,xpos2, &
  ypos2,xnep2,ynep2)
write(0,*) "Number of points read: ",npt2

t1=x2(1)  !start of first segment
i=1
npt=0
do while(x1(i).lt.t1)
   npt=npt+1
   x(npt)=x1(i)
   y(npt)=y1(i)
   yerr(npt)=yerr1(i)
   oflux(npt)=oflux1(i)
   pmod(npt)=pmod1(i)
   smod(npt)=smod1(i)
   xpos(npt)=xpos1(i)
   ypos(npt)=ypos1(i)
   xnep(npt)=xnep1(i)
   ynep(npt)=ynep1(i)
   i=i+1
enddo
!write(0,*) "npt: ",npt

k=2
do while(k.le.6)

   t2=x1(npt1) !go till end of next segment
   j=1
   do while(x1(i).lt.t2)
      if(abs(x1(i)-x2(j)).gt.1.0e-10)then
         write(0,*) "T mismatch: ",x1(i),x2(j)
         write(0,*) i,j
         do while(abs(x1(i)-x2(j)).gt.1.0e-10)
            if(x1(i).gt.x2(j))then
               j=j+1
               if(j.gt.npt2) write(0,*) "j.gt.npt2"
            else
               if(j.eq.1)then
                  npt=npt+1
                  x(npt)=x1(i)
                  y(npt)=y1(i)
                  yerr(npt)=yerr1(i)
                  oflux(npt)=oflux1(i)
                  pmod(npt)=pmod1(i)
                  smod(npt)=smod1(i)
                  xpos(npt)=xpos1(i)
                  ypos(npt)=ypos1(i)
                  xnep(npt)=xnep1(i)
                  ynep(npt)=ynep1(i)
                  t1=x1(i)
               endif
               i=i+1
               if(i.gt.npt1) write(0,*) "i.gt.npt1"
            endif
         enddo
         write(0,*) i,j,x1(i),x2(j)
!         read(5,*)
      endif
      npt=npt+1
      x(npt)=x1(i)
      f1=(t2-x(npt))/(t2-t1)
      f2=(x(npt)-t1)/(t2-t1)
      y(npt)=y1(i)*f1+y2(j)*f2
      yerr(npt)=yerr1(i)

      oflux(npt)=oflux1(i)*f1+oflux2(j)*f2
      pmod(npt)=pmod1(i)*f1+pmod2(j)*f2
      smod(npt)=smod1(i)*f1+smod2(j)*f2

      xpos(npt)=xpos1(i)
      ypos(npt)=ypos1(i)
      xnep(npt)=xnep1(i)
      ynep(npt)=ynep1(i)

      i=i+1
      j=j+1
   enddo
   write(0,*) "i,j: ",i,j

   !move dataset 2->dataset 1
   npt1=npt2
   t1=t2
   x1=x2
   y1=y2
   yerr1=yerr2
   oflux1=oflux2
   pmod1=pmod2
   smod1=smod2
   xpos1=xpos2
   ypos1=ypos2
   xnep1=xnep2
   i=j

   k=k+1
   if(k.le.6)then
      write(filename2,'(A6,I1,A6)') "phot.p",k,".j.dat"
      write(0,*) "K: ",k
      write(0,'(A80)') filename2
      call getpdata(filename2,npt2,nmax,x2,y2,yerr2,oflux2,pmod2,smod2, &
         xpos2,ypos2,xnep2,ynep2)
      write(0,*) "Number of points read: ",npt2
   endif
enddo

do k=i,npt1
   npt=npt+1
   x(npt)=x1(k)
   y(npt)=y1(k)
   yerr(npt)=yerr1(k)
   oflux(npt)=oflux1(k)
   pmod(npt)=pmod1(k)
   smod(npt)=smod1(k)
   xpos(npt)=xpos1(k)
   ypos(npt)=ypos1(k)
   xnep(npt)=xnep1(k)
   ynep(npt)=ynep1(k)
enddo

!scale and shift data.
call exportdata(npt,x,y,yerr,oflux,pmod,smod,xpos,ypos,xnep,ynep)


end program joinparts

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine exportdata(npt,x,y,yerr,oflux,pmod,smod,xpos,ypos,xnep,ynep)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!import vars
integer :: npt
real(double), dimension(npt) :: x,y,yerr,oflux,pmod,smod,xpos,ypos,xnep,&
  ynep
!local vars
integer :: i

do i=1,npt
   write(6,500) x(i),y(i),yerr(i),oflux(i),pmod(i),smod(i),xpos(i),     &
     ypos(i),xnep(i),ynep(i)
   500 format(10(1PE17.10,1X))
enddo

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine getpdata(filename,npt,nmax,x,y,yerr,oflux,pmod,smod,xpos,    &
  ypos,xnep,ynep)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!Routine to readin data for pixelfit
use precision
implicit none
integer :: nunit,filestatus,i,nmax,npt
real(double) :: time,flux,ferr,origflux,pmodel,spmodel,xn,yn, &
 xp,yp
real(double), dimension(:) :: x,y,yerr,oflux,pmod,smod,xpos,ypos,xnep,ynep
character(80) :: filename

nunit=10
open(unit=nunit,file=filename,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",filename
   stop
endif

i=0
do
   if(i.gt.nmax)then
      write(0,*) "Increase nmax to match data points"
      write(0,*) "nmax: ",nmax
      stop
   endif
   read(nunit,*,iostat=filestatus) time,flux,ferr,origflux,pmodel,      &
      spmodel,xn,yn,xp,yp
   if(filestatus == 0) then
      i=i+1
      x(i)=time
      y(i)=flux
      yerr(i)=ferr
      oflux(i)=origflux
      pmod(i)=pmodel
      smod(i)=spmodel
      xpos(i)=xn
      ypos(i)=yn
      xnep(i)=xp
      ynep(i)=yp
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i+1
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit) !close file
npt=i-1

return
end subroutine getpdata
