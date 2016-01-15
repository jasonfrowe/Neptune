CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE lfit(x,y,sig,ndat,a,ia,ma,covar,npc,chisq)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER ma,ia(ma),npc,ndat,MMAX
      REAL*8 chisq,a(ma),covar(npc,npc),sig(ndat),x(ndat),y(ndat)
      EXTERNAL funcs2
      PARAMETER (MMAX=50) !Set to the maximum number of coefficients ma.
C     USES covsrt,gaussj
      INTEGER i,j,k,l,m,mfit
      REAL*8 sig2i,sum,wt,ym,afunc(MMAX),beta(MMAX)

      mfit=0
      do 11 j=1,ma
        if(ia(j).ne.0) mfit=mfit+1
 11   continue
      if(mfit.eq.0) pause "lfit: no parameters to be fitted"
      do 13 j=1,mfit !Initialize the (symmetric) matrix.
        do 12 k=1,mfit
            covar(j,k)=0.
 12     continue
        beta(j)=0.
 13   continue
      do 17 i=1,ndat !Loop over data to accumulate coefficients of the normal
        call funcs2(x(i),afunc,ma) !equations.
        ym=y(i)
        if(mfit.lt.ma) then !Subtract off dependences on known pieces of the fitting
            do 14 j=1,ma !function.
                if(ia(j).eq.0) ym=ym-a(j)*afunc(j)
 14         continue
        endif
        sig2i=1./sig(i)**2
        j=0
        do 16 l=1,ma
            if (ia(l).ne.0) then
                j=j+1
                wt=afunc(l)*sig2i
                k=0
                do 15 m=1,l
                    if (ia(m).ne.0) then
                        k=k+1
                        covar(j,k)=covar(j,k)+wt*afunc(m)
                    endif
 15             continue
                beta(j)=beta(j)+ym*wt
            endif
 16     continue
 17   continue
      do 19 j=2,mfit !Fill in above the diagonal from symmetry.
        do 18 k=1,j-1
            covar(k,j)=covar(j,k)
 18     continue
 19   continue
      call gaussj(covar,mfit,npc,beta,1,1) !Matrix solution.
      j=0
      do 21 l=1,ma
        if(ia(l).ne.0) then
            j=j+1
            a(l)=beta(j) !Partition solution to appropriate coefficients a.
        endif
 21   continue
      chisq=0.
      do 23 i=1,ndat
        call funcs2(x(i),afunc,ma)
        sum=0.
        do 22 j=1,ma
            sum=sum+a(j)*afunc(j)
 22     continue
        chisq=chisq+((y(i)-sum)/sig(i))**2
 23   continue
      call covsrt(covar,npc,ma,ia,mfit) !Sort covariance matrix to true order of fitting
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FUNCS2(X,P,NP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
C
C     Taken from NUMERICAL RECIPES: The Art of Scientific Computing by
C     Press, Flannery, Teukolsky, Vetterling (1986).
C
      integer j,np
      real*8 P(NP),X
      P(1)=1.
      DO 11 J=2,NP
         P(J)=P(J-1)*X
 11   CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE gaussj(a,n,np,b,m,mp)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER m,mp,n,np,NMAX
      REAL*8 a(np,np),b(np,mp)
      PARAMETER (NMAX=1000)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),
     *     ipiv(NMAX)
      REAL*8 big,dum,pivinv
      do 11 j=1,n
         ipiv(j)=0
 11   enddo
      do 22 i=1,n
         big=0.
         do 13 j=1,n
            if(ipiv(j).ne.1)then
               do 12 k=1,n
                  if (ipiv(k).eq.0) then
                     if (abs(a(j,k)).ge.big)then
                        big=abs(a(j,k))
                        irow=j
                        icol=k
                     endif
                  else if (ipiv(k).gt.1) then
                     pause 'singular matrix in gaussj'
                  endif
 12            enddo
            endif
 13      enddo
         ipiv(icol)=ipiv(icol)+1
         if (irow.ne.icol) then
            do 14 l=1,n
               dum=a(irow,l)
               a(irow,l)=a(icol,l)
               a(icol,l)=dum
 14         enddo
            do 15 l=1,m
               dum=b(irow,l)
               b(irow,l)=b(icol,l)
               b(icol,l)=dum
 15         enddo
         endif
         indxr(i)=irow
         indxc(i)=icol
         if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
         pivinv=1./a(icol,icol)
         a(icol,icol)=1.
         do 16 l=1,n
            a(icol,l)=a(icol,l)*pivinv
 16      enddo
         do 17 l=1,m
            b(icol,l)=b(icol,l)*pivinv
 17      enddo
         do 21 ll=1,n
            if(ll.ne.icol)then
               dum=a(ll,icol)
               a(ll,icol)=0.
               do 18 l=1,n
                  a(ll,l)=a(ll,l)-a(icol,l)*dum
 18            enddo
               do 19 l=1,m
                  b(ll,l)=b(ll,l)-b(icol,l)*dum
 19            enddo
            endif
 21      enddo
 22   enddo
      do 24 l=n,1,-1
         if(indxr(l).ne.indxc(l))then
            do 23 k=1,n
               dum=a(k,indxr(l))
               a(k,indxr(l))=a(k,indxc(l))
               a(k,indxc(l))=dum
 23         enddo
         endif
 24   enddo
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER ma,mfit,npc,ia(ma)
      REAL*8 covar(npc,npc)
      INTEGER i,j,k
      REAL*8 swap
      do 12 i=mfit+1,ma
         do 11 j=1,i
            covar(i,j)=0.
            covar(j,i)=0.
 11      enddo
 12   enddo
      k=mfit
      do 15 j=ma,1,-1
         if(ia(j).ne.0)then
            do 13 i=1,ma
               swap=covar(i,k)
               covar(i,k)=covar(i,j)
               covar(i,j)=swap
 13         enddo
            do 14 i=1,ma
               swap=covar(k,i)
               covar(k,i)=covar(j,i)
               covar(j,i)=swap
 14         enddo
            k=k-1
         endif
 15   enddo
      return
      END
