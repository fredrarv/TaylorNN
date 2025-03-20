module linesearch
use globals,        only: P, Nslips, Nslipsmax
use SCYLglobals,    only: SCYLexp
contains
      subroutine findstress(sig, Dp, CRSS, errsig)
!     finds the stress, sig, with normal Dp
      implicit none
      integer i,j,nexp,iter,iflag,Nitermax
      double precision errsig,errsig0,CRSS(:)
      double precision Dp(3,3),errx,erry,k1,k2,tauc(Nslipsmax)
      double precision Dpnorm, sig(3,3), df(3,3),Dpdf(3,3)
      double precision dot1,dot2,fac
      double precision sig0(3,3),dsig0(3,3),dpdir(3,3),dfnorm,tmp
      common/srcline/sig0,dsig0,dpdir
      common/slip/tauc
      common/str/nexp,Nitermax

      Nitermax = 200
      nexp = SCYLexp
      tauc(:Nslips) = CRSS
      errsig0=errsig
      iter=1
      Dpnorm=0.
      call mat3dot(Dpnorm,Dp,Dp)
      Dpnorm=Sqrt(Dpnorm)

      do i =1,3
         do j =1,3
            dpdir(i,j)=Dp(i,j)/Dpnorm
         enddo
      enddo
 333  continue
      do i =1,3
         do j =1,3
            sig0(i,j)=sig(i,j)
         enddo
      enddo
      call ysproj(sig0)
      call derivate(sig0,df,nexp)
!
 222  continue
!      
      call mat3dot(dfnorm,df,df)
      dfnorm=sqrt(dfnorm)
      do i=1,3
         do j=1,3
            df(i,j)=df(i,j)/dfnorm
         enddo
      enddo
      call mat3dot(fac,Dpdir,df)
      tmp = 0.
      do i =1,3
         do j=1,3
            dsig0(i,j)=Dpdir(i,j)-fac*df(i,j)
            tmp=tmp+dsig0(i,j)**2
         enddo
      enddo
      tmp=sqrt(tmp)
      do i=1,3
         do j=1,3
            dsig0(i,j)=dsig0(i,j)/tmp
         enddo
      enddo
!     line search sig= sig0 + k1*dsig0
      errx=1.E-6
      erry=1.E-6
      k1=0.d0
      k2=2.*tauc(1)
      iflag=1
      call Nmod(linesrc,k1,k2,erry,errx,iflag)
      if(iflag.eq.-2)then
         iflag=-2
         write(*,*) "stress: no convergency after Nitermax iterations"
      endif
      if(iflag.eq.-1)then
         write(*,*) "stress: both guesses gives same sign of function"
         k1=1.
         do i =1,3
            do j =1,3
               Sig(i,j)=Sig0(i,j)+k1*dsig0(i,j)
               sig0(i,j)=sig(i,j)
            enddo
         enddo
         stop
      endif
!
      do i =1,3
         do j =1,3
            Sig(i,j)=Sig0(i,j)+k1*dsig0(i,j)
            sig0(i,j)=sig(i,j)
         enddo
      enddo
      call ysproj(sig0)
      call ysproj(sig)
!
      call derivate(sig,df,nexp)
      call mat3dot(dot1,df,df)
      call mat3dot(dot2,Dpdir,df)
      errsig= dot2/sqrt(dot1)-1.d0
      iter=iter+1
      if(iter.gt.Nitermax) goto 223
      if(abs(errsig).gt.errsig0) goto 222
 223  continue
!      write(*,*)"iter on paths= ",iter," err= ", errsig
      return
      end
!
!__________________________________________      
      subroutine linesrc(y,x)
!     y is targeted to zero by changing x
      implicit none
      integer nexp, i, j, Nitermax
      double precision sigma(3,3)
      double precision x, y,ff,dfnorm
      double precision sig0(3,3),dsig0(3,3),df(3,3),dpdir(3,3)
      common /srcline/sig0,dsig0,Dpdir
      common/str/nexp,Nitermax
!
      do i=1,3
         do j=1,3
            sigma(i,j)=sig0(i,j)+x*dsig0(i,j)
         enddo
      enddo
      call ysproj(sigma)
      call derivate(sigma,df,nexp)
      
      call mat3dot(dfnorm,df,df)
      dfnorm = sqrt(dfnorm)
      y=0.
      do i =1,3
         do j=1,3
            y = y+dsig0(i,j)*(dpdir(i,j)-df(i,j)/dfnorm)
         enddo
      enddo
      return
      end
!
!
!_______________________________________
      subroutine ysproj(sig)
      implicit none
      integer i,j,nexp,Nitermax
      double precision sig(3,3),tmp
      common/str/nexp,Nitermax
      call yieldsurf(tmp, sig, nexp)
!     find solution, utilizing first order homologeous
      do i=1,3
         do j=1,3
            sig(i,j)=sig(i,j)/(1.+tmp)
         enddo
      enddo
      return
      end
!
!____________________________________________    
      subroutine Nmod(fun,x1,x2,erry,errx,iflag)
!     modified Newton
!     the solution has to be between x1 and x2 
      implicit none
      integer i,iflag,isection,nexp,Nitermax
      double precision errx,erry,x1,x2, x3, dx3,y1,y2,y3
      double precision  xmax, xmin, signmin,signmax, ymin,ymax
      common/str/nexp,Nitermax
      external fun
      
      interface
         subroutine fun(y, x)
            double precision y, x
         end subroutine fun
      end interface
      
      iflag=1
      i=0
      isection=0
      xmin=min(x1,x2)
      xmax=max(x1,x2)
      call fun(y1,x1)
      x2=x1+0.01*(x2-x1)
      call fun(y2,x2)
      if(abs(y1).lt.erry) then
         x3=x1
         goto 111
      endif
      if(abs(y2).lt.erry) then
         x3=x2
         goto 111
      endif
!      
      call fun(ymin,xmin)
      if(abs(ymin).lt.erry) then
         x3=xmin
         goto 111
      endif
      signmin=ymin/abs(ymin)
      do i =1,5
         call fun(ymax,xmax)
         if(abs(ymax).lt.erry) then
            x3=xmax
            goto 111
         endif
         signmax=ymax/abs(ymax)
         if(signmin*signmax.lt.0) goto 444
         xmax=xmax+abs(xmax-xmin)*0.2
!     xmin=xmin-abs(xmax-xmin)*0.2
      enddo
 444  continue
!     
      do i=1,20      
         if(abs(y2-y1).gt.1E-8) then
            dx3=-y2*(x2-x1)/(y2-y1)
            x3=x2+dx3
            call fun(y3,x3)
         else
            x3=0.5*(xmin+xmax)
            call fun(y3,x3)
            isection=isection+1
         endif
         if((abs(y3).lt.erry).and.(abs(x3-x2).lt.errx)) goto 111
         if((x3.lt.xmin).or.(x3.gt.xmax)) then
            x3=0.5*(xmin+xmax)
            call fun(y3,x3)
            isection=isection+1
         endif
         x1=x2
         y1=y2
         x2=x3
         y2=y3
         if(signmin*y2.gt.0.) xmin=max(xmin,x2)
         if(signmax*y2.gt.0.) xmax=min(xmax,x2)
      enddo
      iflag=-2
 111  continue
      errx=abs(x3-x2)
      erry=abs(y3)
      x1=x3
      i=i-isection
!      write(*,*) i,isection
      return
      end
      
!_________________________________________

      double precision function brak(x)
      implicit none
      double precision x
      brak = max(0.,x)
      return
      end      
!_________________________________________
      double precision function dotprod(a,b,s)
      integer i,j,s
      double precision a(3,3), b(:,:,:)
      dotprod = 0.
      do i =1,3
         do j=1,3
            dotprod = dotprod+a(i,j)*b(i,j,s)
         enddo
      enddo
      return
      end
!_________________________________________
      subroutine yieldsurf(f, sigma, nexp)
      implicit none
      integer nexp, s,i,j
      double precision f, sigma(3,3)
      double precision tmp(Nslipsmax),signorm
      double precision dot
      double precision sighat(3,3), tauc(Nslipsmax)
      common/slip/tauc
!     normalise sig into s, to avoid numerical overflow
      ! scaling by the maximum element
      tmp = 0.d0
      do s=1, Nslips
         tmp(s) = dotprod(sigma,P,s)/tauc(s)
         tmp(s) = brak(tmp(s))
      end do
      signorm = maxval(tmp, dim=1)
      tmp = tmp/signorm
      f=0.
      do s=1, Nslips
            f = f + tmp(s)**nexp
      end do  
      f = signorm*f**(1./nexp)-1.
      return
      end
!     
!_________________________________________
      subroutine derivate(sigma,df,nexp)
      implicit none
      integer nexp
      integer s,i,j
      double precision sigma(3,3), df(3,3)
      double precision bracket(Nslipsmax), tmp
      double precision f, pre, tauc(Nslipsmax)
      common/slip/tauc
!      
      call yieldsurf(f, sigma, nexp)
      pre=(f+1.)**(1.-nexp)
      do s=1,Nslips
         tmp = dotprod(sigma,P,s)/tauc(s)
         bracket(s)=brak(tmp)
      enddo
      do i=1,3
         do j=1,3
            df(i,j)=0.
            do s=1,Nslips
               df(i,j)=df(i,j)+bracket(s)**(nexp-1)*P(i,j,s)/tauc(s)
            enddo
         enddo
      enddo
!
      do i=1,3
         do j=1,3
            df(i,j)=pre*df(i,j)
         enddo
      enddo
!
      return
      end
!
!___________________________________________
      subroutine multmat(C,A,B)
      implicit none
      integer i,j,k
      double precision C(3,3),A(3,3), B(3,3)
      do i = 1,3
         do j = 1,3
            C(i,j)=0.
            do k =1,3
               C(i,j)=C(i,j)+A(i,k)*B(k,j)
            enddo
         enddo
      enddo
      return
      end
!
!_______________________________________
      subroutine mat3dot(dot,A,B)
      implicit none
      integer i,j
      double precision A(3,3),B(3,3),dot
      dot=0.
      do i=1,3
         do j=1,3
            dot=dot+A(i,j)*B(i,j)
         enddo
      enddo
      return
      end
      
end module linesearch