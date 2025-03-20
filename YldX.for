      subroutine Hill48(svec, params, mode, func, grad, hessian)
c**********************************************************************
c      Calculates Hill48 yield function and/or its gradient
c      and/or its Hessian expressed in new notation. 
c
c**********************************************************************
c (IN)    SVEC is REAL*8 array, dimension (6)
c         It is stress tensor given in Voigt notation
c (IN)    PARAMS is REAL*8 array, dimension (6)
c         It is 6 anisotropic parameters of Hill48 building the Hill
c         matrix
c
c         | G+H  -H -G   0   0   0 |         
c         | -H  F+H -F   0   0   0 |
c         | -G   -F F+G  0   0   0 |
c	    |  0   0   0   2L  0   0 |
c         |  0   0   0   0  2M   0 |
c         |  0   0   0   0   0  2N |
c
c (IN)    MODE is INTEGER
c         Controls which of F, GRAD, HESSIAN is to be calculated
c         = 0  :  only F is calculated
c         = 1  :  only F and GRAD is calculated
c         else :  F, GRAD and HESSIAN are calculated
c (OUT)   F is REAL*8
c         It is scalar value of Yld2004 yield function
c (OUT)   GRAD is REAL*8 array, dimension (6)
c         It is gradient of Yld2004 yield function, computed in 
c         Voigt notation
c (OUT)   HESSIAN is REAL*8 array, dimension (6,6)
c         It is the Hessian of Yld2004 yield function computed in 
c         Voigt notation
c         
c**********************************************************************
      implicit none
c     ! in/out
      integer     mode
      real*8      svec(6), params(6), func, grad(6), hessian(6,6)
c     ! internal
      real*8      F, G, H, L, M, N, ovfunc
c
      F = params(1)
      G = params(2)
      H = params(3)
      L = params(4)
      M = params(5)
      N = params(6)
c
      func = sqrt(F*(svec(2)-svec(3))**2
     +     +      G*(svec(3)-svec(1))**2
     +     +      H*(svec(1)-svec(2))**2
     +     + 2.d0*(L*svec(4)**2+M*svec(5)**2+N*svec(6)**2))
c
      if (mode .EQ. 0) return
c
      ovfunc = 1.d0/func
      grad(1) = (G*(svec(1)-svec(3)) + H*(svec(1)-svec(2)))*ovfunc
      grad(2) = (F*(svec(2)-svec(3)) + H*(svec(2)-svec(1)))*ovfunc
      grad(3) = (F*(svec(3)-svec(2)) + G*(svec(3)-svec(1)))*ovfunc
      grad(4) = 2.d0*L*svec(4)*ovfunc
      grad(5) = 2.d0*M*svec(5)*ovfunc
      grad(6) = 2.d0*N*svec(6)*ovfunc
c
      if (mode .EQ. 1) return
c
      hessian(1,1) = (G+H - grad(1)*grad(1))*ovfunc
      hessian(1,2) = ( -H - grad(1)*grad(2))*ovfunc
      hessian(1,3) = ( -G - grad(1)*grad(3))*ovfunc
      hessian(1,4) = -grad(1)*grad(4)*ovfunc  
      hessian(1,5) = -grad(1)*grad(5)*ovfunc  
      hessian(1,6) = -grad(1)*grad(6)*ovfunc  
      hessian(2,2) = (F+H - grad(2)*grad(2))*ovfunc
      hessian(2,3) = ( -F - grad(2)*grad(3))*ovfunc
      hessian(2,4) = -grad(2)*grad(4)*ovfunc  
      hessian(2,5) = -grad(2)*grad(5)*ovfunc  
      hessian(2,6) = -grad(2)*grad(6)*ovfunc  
      hessian(3,3) = (F+G - grad(3)*grad(3))*ovfunc
      hessian(3,4) = -grad(3)*grad(4)*ovfunc  
      hessian(3,5) = -grad(3)*grad(5)*ovfunc  
      hessian(3,6) = -grad(3)*grad(6)*ovfunc  
      hessian(4,4) = (2.d0*L - grad(4)*grad(4))*ovfunc 
      hessian(4,5) = -grad(4)*grad(5)*ovfunc           
      hessian(4,6) = -grad(4)*grad(6)*ovfunc           
      hessian(5,5) = (2.d0*M - grad(5)*grad(5))*ovfunc 
      hessian(5,6) = -grad(5)*grad(6)*ovfunc           
      hessian(6,6) = (2.d0*N - grad(6)*grad(6))*ovfunc 
c
      return
      end subroutine Hill48
c
c
      subroutine KB(svec, cparams, a, f)
c**********************************************************************
c      Calculates Karafilis and Boyce yield function
c**********************************************************************
c (IN)    SVEC is REAL*8 array, dimension (6)
c         It is stress tensor given in Voigt notation
c (IN)    CPARAMS is REAL*8 array, dimension (6)
c         It is 18 anisotropic parameters of Yld2004 expected in Voigt
c         notation, i.e. the first L anisotropic transformation matrix
c         read
c           |  1    beta1   beta2    0     0      0   |         
c           |beta1 alpha1   beta3    0     0      0   |
c         C*|beta2  beta3  alpha2    0     0      0   |
c           |  0      0       0   gamma1   0      0   |
c           |  0      0       0      0   gamma2   0   |
c           |  0      0       0      0     0    gamma3|
c
c         Coefficients in LPARAMS are ordered as
c         C, alpha1, alpha2, gamma1, gamma2, gamma3 
c (IN)    A is REAL*8
c         It is the exponent in Yld2004 (a must be >= 2)
c (OUT)   F is REAL*8
c         It is scalar value of Yld2004 yield function
c         
c**********************************************************************
      implicit none
c     ! in/out
      real*8      svec(6), cparams(6), a, f
c     ! internal
      integer     i
      real*8      s1vec(6), scale, asms(3), ova, xi
c
      real*8      c, alpha1, alpha2, gamma1, gamma2, gamma3,
     +            beta1, beta2, beta3,
     +            s11, s12, s13, Smat(3,3),
     +            eigvec11(3), eigvec12(3), eigvec13(3)
c
Cf2py intent(out) f
c
	  ova     = 1.d0/a
c      
      c  = cparams(1)
      alpha1  = cparams(2)
      alpha2  = cparams(3)
      gamma1  = cparams(4)
      gamma2  = cparams(5)
      gamma3  = cparams(6)
c
      beta1 = (alpha2-alpha1-1.d0)/2.d0
      beta2 = (alpha1-alpha2-1.d0)/2.d0
      beta3 = (1.d0-alpha1-alpha2)/2.d0
      !
      s1vec(1) = c*(svec(1)+beta1*svec(2)+beta2*svec(3))
      s1vec(2) = c*(beta1*svec(1)+alpha1*svec(2)+beta3*svec(3))
      s1vec(3) = c*(beta2*svec(1)+beta3*svec(2)+alpha2*svec(3))
      s1vec(4) = c*gamma1*svec(4)
      s1vec(5) = c*gamma2*svec(5)
      s1vec(6) = c*gamma3*svec(6)
      
c     Calculate eigenvalues and eigenvectors of transformed stresses
      Smat = reshape(
     +     (/ s1vec(1), s1vec(6), s1vec(5),
     +        s1vec(6), s1vec(2), s1vec(4),
     +        s1vec(5), s1vec(4), s1vec(3) /), (/3,3/))
      call eig2(Smat,0,s11,s12,s13,eigvec11,eigvec12,eigvec13)
c      
c     
      asms(1) = abs(s11)
      asms(2) = abs(s12)
      asms(3) = abs(s13)
c
      scale = asms(1)
      do i=2, 3
          if (asms(i) .GT. scale) scale = asms(i)
      enddo
c     
      xi = (2.d0**a+2.d0)/(3.d0**a)
      f = 0.d0
      if (scale .GT. 1.d-16) then
          do i=1, 3
              f = f + (asms(i)/scale)**a
          enddo
c
c         Compute yield function F
          f = scale*(f/xi)**ova
      end if
c
      return    
      end subroutine KB
c
c
      subroutine Yld2004_9p(svec, cparams, a, f)
c**********************************************************************
c      Calculates Barlat's Yld2004-18p yield function
c**********************************************************************
c (IN)    SVEC is REAL*8 array, dimension (6)
c         It is stress tensor given in Voigt notation
c (IN)    CPARAMS is REAL*8 array, dimension (9)
c         It is 18 anisotropic parameters of Yld2004 expected in Voigt
c         notation, i.e. the first L anisotropic transformation matrix
c         read
c         |  0  -C12 -C13  0   0   0 |         
c         |-C21   0  -C23  0   0   0 |
c         |-C31 -C32   0   0   0   0 |
c	      |  0    0    0  C44  0   0 |
c         |  0    0    0   0  C55  0 |
c         |  0    0    0   0   0  C66|
c
c         Coefficients in LPARAMS are ordered as
c         C'12, C'13, C'21, C'23, C'31, C'32, C'44, C'55, C'66, 
c (IN)    A is REAL*8
c         It is the exponent in Yld2004 (a must be >= 2)
c (OUT)   F is REAL*8
c         It is scalar value of Yld2004 yield function
c         
c**********************************************************************
      implicit none
c     ! in/out
      real*8      svec(6), cparams(9), a, f
c     ! internal
      integer     i, j, emod
      real*8      s1vec(6), scale, asms(3), ov3, ova
c
      real*8      c12, c13, c21, c23, c31, c32, c44, c55, c66,          
     +            a11, a12, a13, a21, a22, a23, a31, a32, a33,
     +            s11, s12, s13, Smat(3,3),
     +            eigvec11(3), eigvec12(3), eigvec13(3)
c
Cf2py intent(out) f
c
      ov3     = 1.d0/3.d0
	ova     = 1.d0/a
c      
      c12  = cparams(1)
      c13  = cparams(2)
      c21  = cparams(3)
      c23  = cparams(4)
      c31  = cparams(5)
      c32  = cparams(6)
      c44  = cparams(7)
      c55  = cparams(8)
      c66  = cparams(9)
c
      a11 = ov3*(c12+c13)
      a12 = ov3*(c13-2.0d0*c12)
      a13 = ov3*(c12-2.0d0*c13)
      a21 = ov3*(c23-2.0d0*c21)
      a22 = ov3*(c21+c23)
      a23 = ov3*(c21-2.0d0*c23)
      a31 = ov3*(c32-2.0d0*c31)
      a32 = ov3*(c31-2.0d0*c32)
      a33 = ov3*(c31+c32)
      !
      s1vec(1) = a11*svec(1)+a12*svec(2)+a13*svec(3)
      s1vec(2) = a21*svec(1)+a22*svec(2)+a23*svec(3)
      s1vec(3) = a31*svec(1)+a32*svec(2)+a33*svec(3)
      s1vec(4) = c44*svec(4)
      s1vec(5) = c55*svec(5)
      s1vec(6) = c66*svec(6)
      
c     Calculate eigenvalues and eigenvectors of transformed stresses
      Smat = reshape(
     +     (/ s1vec(1), s1vec(6), s1vec(5),
     +        s1vec(6), s1vec(2), s1vec(4),
     +        s1vec(5), s1vec(4), s1vec(3) /), (/3,3/))
      call eig2(Smat,0,s11,s12,s13,eigvec11,eigvec12,eigvec13)
c      
c     
      asms(1) = abs(s11-s12)
      asms(2) = abs(s11-s13)
      asms(3) = abs(s12-s13)
c
      scale = asms(1)
      do i=2, 3
          if (asms(i) .GT. scale) scale = asms(i)
      enddo
c     
      f = 0.d0
      if (scale .GT. 1.d-16) then
          do i=1, 3
              f = f + (asms(i)/scale)**a
          enddo
c
c         Compute yield function F
          f = scale*(f/2.d0)**ova
      end if
c
      return    
      end subroutine Yld2004_9p
c
c
      subroutine Yld2004_18p(svec, cparams, a, f)
c**********************************************************************
c      Calculates Barlat's Yld2004-18p yield function
c**********************************************************************
c (IN)    SVEC is REAL*8 array, dimension (6)
c         It is stress tensor given in Voigt notation
c (IN)    CPARAMS is REAL*8 array, dimension (18)
c         It is 18 anisotropic parameters of Yld2004 expected in Voigt
c         notation, i.e. the first C anisotropic transformation matrix
c         read
c         |  0  -C12 -C13  0   0   0 |         
c         |-C21   0  -C23  0   0   0 |
c         |-C31 -C32   0   0   0   0 |
c	    |  0    0    0  C44  0   0 |
c         |  0    0    0   0  C55  0 |
c         |  0    0    0   0   0  C66|
c
c         Coefficients in CPARAMS are ordered as
c         C'12, C'13, C'21, C'23, C'31, C'32, C'44, C'55, C'66, 
c         followed by coefficients associated with the second anistropic 
c         transformation matrix C''
c (IN)    A is REAL*8
c         It is the exponent in Yld2004 (a must be >= 1)
c (OUT)   F is REAL*8
c         It is scalar value of Yld2004 yield function
c         
c**********************************************************************
      implicit none
c     ! in/out
      real*8      svec(6), cparams(18), a, f
c     ! internal
      integer     i, j, emod
      real*8      s1vec(6), s2vec(6), scale, asms(9), ov3, ova
c
      real*8      c12, c13, c21, c23, c31, c32, c44, c55, c66,          
     +            d12, d13, d21, d23, d31, d32, d44, d55, d66, 
     +            a11, a12, a13, a21, a22, a23, a31, a32, a33,
     +            b11, b12, b13, b21, b22, b23, b31, b32, b33, 
     +            s11, s12, s13, s21, s22, s23, Smat(3,3),
     +            eigvec11(3), eigvec12(3), eigvec13(3),                
     +            eigvec21(3), eigvec22(3), eigvec23(3)
c
Cf2py intent(out) f
c
      ov3     = 1.d0/3.d0
	ova     = 1.d0/a
c      
      c12  = cparams(1)
      c13  = cparams(2)
      c21  = cparams(3)
      c23  = cparams(4)
      c31  = cparams(5)
      c32  = cparams(6)
      c44  = cparams(7)
      c55  = cparams(8)
      c66  = cparams(9)
      d12  = cparams(10)
      d13  = cparams(11)
      d21  = cparams(12)
      d23  = cparams(13)
      d31  = cparams(14)
      d32  = cparams(15)
      d44  = cparams(16)
      d55  = cparams(17)
      d66  = cparams(18)
c
      a11 = ov3*(c12+c13)
      a12 = ov3*(c13-2.0d0*c12)
      a13 = ov3*(c12-2.0d0*c13)
      a21 = ov3*(c23-2.0d0*c21)
      a22 = ov3*(c21+c23)
      a23 = ov3*(c21-2.0d0*c23)
      a31 = ov3*(c32-2.0d0*c31)
      a32 = ov3*(c31-2.0d0*c32)
      a33 = ov3*(c31+c32)
      !
      b11 = ov3*(d12+d13)
      b12 = ov3*(d13-2.0d0*d12)
      b13 = ov3*(d12-2.0d0*d13)
      b21 = ov3*(d23-2.0d0*d21)
      b22 = ov3*(d21+d23)
      b23 = ov3*(d21-2.0d0*d23)
      b31 = ov3*(d32-2.0d0*d31)
      b32 = ov3*(d31-2.0d0*d32)
      b33 = ov3*(d31+d32)
      !
      s1vec(1) = a11*svec(1)+a12*svec(2)+a13*svec(3)
      s1vec(2) = a21*svec(1)+a22*svec(2)+a23*svec(3)
      s1vec(3) = a31*svec(1)+a32*svec(2)+a33*svec(3)
      s1vec(4) = c44*svec(4)
      s1vec(5) = c55*svec(5)
      s1vec(6) = c66*svec(6)
      !
      s2vec(1) = b11*svec(1)+b12*svec(2)+b13*svec(3)
      s2vec(2) = b21*svec(1)+b22*svec(2)+b23*svec(3)
      s2vec(3) = b31*svec(1)+b32*svec(2)+b33*svec(3)
      s2vec(4) = d44*svec(4)
      s2vec(5) = d55*svec(5)
      s2vec(6) = d66*svec(6)
      
c     Calculate eigenvalues and eigenvectors of transformed stresses
      Smat = reshape(
     +     (/ s1vec(1), s1vec(6), s1vec(5),
     +        s1vec(6), s1vec(2), s1vec(4),
     +        s1vec(5), s1vec(4), s1vec(3) /), (/3,3/))
      call eig2(Smat,0,s11,s12,s13,eigvec11,eigvec12,eigvec13)
c
      Smat = reshape(
     +     (/ s2vec(1), s2vec(6), s2vec(5),
     +        s2vec(6), s2vec(2), s2vec(4),
     +        s2vec(5), s2vec(4), s2vec(3) /), (/3,3/))
      call eig2(Smat,0,s21,s22,s23,eigvec21,eigvec22,eigvec23)
c      
c     
      asms(1) = abs(s11-s21)
      asms(2) = abs(s11-s22)
      asms(3) = abs(s11-s23)
c              
      asms(4) = abs(s12-s21)
      asms(5) = abs(s12-s22)
      asms(6) = abs(s12-s23)
c              
      asms(7) = abs(s13-s21)
      asms(8) = abs(s13-s22)
      asms(9) = abs(s13-s23)
c
      scale = asms(1)
      do i=2, 9
          if (asms(i) .GT. scale) scale = asms(i)
      enddo
c     
      f = 0.d0
      if (scale .GT. 1.d-16) then
          do i=1, 9
              f = f + (asms(i)/scale)**a
          enddo
c
c         Compute yield function F
          f = scale*(f/4.d0)**ova
      end if
c
      return    
      end subroutine Yld2004_18p
c
c
      subroutine Yld2004nat_18p(svec, cparams, a, f)
c**********************************************************************
c      Calculates Barlat's Yld2004-18p yield function
c**********************************************************************
c (IN)    SVEC is REAL*8 array, dimension (6)
c         It is stress tensor given in Voigt notation
c (IN)    CPARAMS is REAL*8 array, dimension (18)
c         It is 18 anisotropic parameters of Yld2004 expected in new
c         notation, i.e. the first L anisotropic transformation matrix
c         read
c         | 0  L12 L13  0   0   0 |         
c         | 0  L22 L23  0   0   0 |
c         | 0  L32 L33  0   0   0 |
c         | 0   0   0  L44  0   0 |
c         | 0   0   0   0  L55  0 |
c         | 0   0   0   0   0  L66|
c
c         Coefficients in CPARAMS are ordered as
c         L'12, L'13, L'22, L'23, L'32, L'33, L'44, L'55, L'66, 
c         followed by coeffiients associated with the second anistropic 
c         transformation matrix L'' as
c         L''12, L''13, L''22, L''23, L''32, L''33, L''44, L''55, L''66
c (IN)    A is REAL*8
c         It is the exponent in Yld2004 (a must be >= 1)
c (OUT)   F is REAL*8
c         It is scalar value of Yld2004 yield function
c         
c**********************************************************************
      implicit none
c     ! in/out
      real*8      svec(6), cparams(18), a, f
c     ! internal
      integer     i, j, emod
      real*8      s1vec(6), s2vec(6), tmp(6), scale, asms(9), ova
c
      real*8      L12, L13, L22, L23, L32, L33, L44, L55, L66,          
     +            M12, M13, M22, M23, M32, M33, M44, M55, M66, 
     +            s11, s12, s13, s21, s22, s23, Smat(3,3),
     +            eigvec11(3), eigvec12(3), eigvec13(3),                
     +            eigvec21(3), eigvec22(3), eigvec23(3)
c
Cf2py intent(out) f
c
	ova     = 1.d0/a
c      
      L12  = cparams(1)
      L13  = cparams(2)
      L22  = cparams(3)
      L23  = cparams(4)
      L32  = cparams(5)
      L33  = cparams(6)
      L44  = cparams(7)
      L55  = cparams(8)
      L66  = cparams(9)
      M12  = cparams(10)
      M13  = cparams(11)
      M22  = cparams(12)
      M23  = cparams(13)
      M32  = cparams(14)
      M33  = cparams(15)
      M44  = cparams(16)
      M55  = cparams(17)
      M66  = cparams(18)
c
c     convert svec from Voigt to natural notation
c     hydrostatic part is not stored
	  tmp(1) = (2.d0*svec(3)-svec(1)-svec(2))/sqrt(6.d0)
	  tmp(2) = (svec(2)-svec(1))/sqrt(2.d0)
	  tmp(3) = sqrt(2.d0)*svec(4)
	  tmp(4) = sqrt(2.d0)*svec(5)
	  tmp(5) = sqrt(2.d0)*svec(6)
c
      s1vec(1) = L12*tmp(1) + L13*tmp(2)
	s1vec(2) = L22*tmp(1) + L23*tmp(2)
	s1vec(3) = L32*tmp(1) + L33*tmp(2)
      s1vec(4) = L44*tmp(3)
      s1vec(5) = L55*tmp(4)
      s1vec(6) = L66*tmp(5)
c
      s2vec(1) = M12*tmp(1) + M13*tmp(2)
	s2vec(2) = M22*tmp(1) + M23*tmp(2)
	s2vec(3) = M32*tmp(1) + M33*tmp(2)
      s2vec(4) = M44*tmp(3)
      s2vec(5) = M55*tmp(4)
      s2vec(6) = M66*tmp(5)
      
c     Calculate eigenvalues and eigenvectors of transformed stresses
      Smat = reshape(
     +     (/ s1vec(1), s1vec(6), s1vec(5),
     +        s1vec(6), s1vec(2), s1vec(4),
     +        s1vec(5), s1vec(4), s1vec(3) /), (/3,3/))
      call eig2(Smat,0,s11,s12,s13,eigvec11,eigvec12,eigvec13)
c
      Smat = reshape(
     +     (/ s2vec(1), s2vec(6), s2vec(5),
     +        s2vec(6), s2vec(2), s2vec(4),
     +        s2vec(5), s2vec(4), s2vec(3) /), (/3,3/))
      call eig2(Smat,0,s21,s22,s23,eigvec21,eigvec22,eigvec23)
c      
c     
      asms(1) = abs(s11-s21)
      asms(2) = abs(s11-s22)
      asms(3) = abs(s11-s23)
c              
      asms(4) = abs(s12-s21)
      asms(5) = abs(s12-s22)
      asms(6) = abs(s12-s23)
c              
      asms(7) = abs(s13-s21)
      asms(8) = abs(s13-s22)
      asms(9) = abs(s13-s23)
c
      scale = asms(1)
      do i=2, 9
          if (asms(i) .GT. scale) scale = asms(i)
      enddo
c     
      f = 0.d0
      if (scale .GT. 1.d-16) then
          do i=1, 9
              f = f + (asms(i)/scale)**a
          enddo
c
c         Compute yield function F
          f = scale*(f/4.d0)**ova
      end if
c
      return    
      end subroutine Yld2004nat_18p
c
c
      subroutine Yld2004_27p(svec, cparams, a, f)
c**********************************************************************
c      Calculates Barlat's Yld2004-27p yield function
c**********************************************************************
c (IN)    SVEC is REAL*8 array, dimension (6)
c         It is stress tensor given in Voigt notation
c (IN)    CPARAMS is REAL*8 array, dimension (27)
c         It is 27 anisotropic parameters of Yld2004 expected in Voigt
c         notation, i.e. the first L anisotropic transformation matrix
c         read
c         |  0  -C12 -C13  0   0   0 |         
c         |-C21   0  -C23  0   0   0 |
c         |-C31 -C32   0   0   0   0 |
c	    |  0    0    0  C44  0   0 |
c         |  0    0    0   0  C55  0 |
c         |  0    0    0   0   0  C66|
c
c         Coefficients in LPARAMS are ordered as
c         C'12, C'13, C'21, C'23, C'31, C'32, C'44, C'55, C'66, 
c         followed by coefficients associated with the second and third
c         anistropic transformation matrices L'' and L'''
c (IN)    A is REAL*8
c         It is the exponent in Yld2004 (a must be >= 2)
c (OUT)   F is REAL*8
c         It is scalar value of Yld2004 yield function
c         
c**********************************************************************
      implicit none
c     ! in/out
      real*8      svec(6), cparams(27), a, f
c     ! internal
      integer     i, j, emod
      real*8      s1vec(6), s2vec(6), s3vec(6), scale, asms(12), 
     +            Smat(3,3), ov3, ova
c
      real*8      c12, c13, c21, c23, c31, c32, c44, c55, c66,          
     +            d12, d13, d21, d23, d31, d32, d44, d55, d66, 
     +            e12, e13, e21, e23, e31, e32, e44, e55, e66, 
     +            a11, a12, a13, a21, a22, a23, a31, a32, a33,
     +            b11, b12, b13, b21, b22, b23, b31, b32, b33, 
     +            z11, z12, z13, z21, z22, z23, z31, z32, z33, 
     +            s11, s12, s13, s21, s22, s23, s31, s32, s33,
     +            eigvec11(3), eigvec12(3), eigvec13(3),                
     +            eigvec21(3), eigvec22(3), eigvec23(3)
c
Cf2py intent(out) f
c
      ov3     = 1.d0/3.d0
	ova     = 1.d0/a
c      
      c12  = cparams(1)
      c13  = cparams(2)
      c21  = cparams(3)
      c23  = cparams(4)
      c31  = cparams(5)
      c32  = cparams(6)
      c44  = cparams(7)
      c55  = cparams(8)
      c66  = cparams(9)
      d12  = cparams(10)
      d13  = cparams(11)
      d21  = cparams(12)
      d23  = cparams(13)
      d31  = cparams(14)
      d32  = cparams(15)
      d44  = cparams(16)
      d55  = cparams(17)
      d66  = cparams(18)
      e12  = cparams(19)
      e13  = cparams(20)
      e21  = cparams(21)
      e23  = cparams(22)
      e31  = cparams(23)
      e32  = cparams(24)
      e44  = cparams(25)
      e55  = cparams(26)
      e66  = cparams(27)
c
      a11 = ov3*(c12+c13)
      a12 = ov3*(c13-2.0d0*c12)
      a13 = ov3*(c12-2.0d0*c13)
      a21 = ov3*(c23-2.0d0*c21)
      a22 = ov3*(c21+c23)
      a23 = ov3*(c21-2.0d0*c23)
      a31 = ov3*(c32-2.0d0*c31)
      a32 = ov3*(c31-2.0d0*c32)
      a33 = ov3*(c31+c32)
      !
      b11 = ov3*(d12+d13)
      b12 = ov3*(d13-2.0d0*d12)
      b13 = ov3*(d12-2.0d0*d13)
      b21 = ov3*(d23-2.0d0*d21)
      b22 = ov3*(d21+d23)
      b23 = ov3*(d21-2.0d0*d23)
      b31 = ov3*(d32-2.0d0*d31)
      b32 = ov3*(d31-2.0d0*d32)
      b33 = ov3*(d31+d32)
      !
      z11 = ov3*(e12+e13)
      z12 = ov3*(e13-2.0d0*e12)
      z13 = ov3*(e12-2.0d0*e13)
      z21 = ov3*(e23-2.0d0*e21)
      z22 = ov3*(e21+e23)
      z23 = ov3*(e21-2.0d0*e23)
      z31 = ov3*(e32-2.0d0*e31)
      z32 = ov3*(e31-2.0d0*e32)
      z33 = ov3*(e31+e32)
      !
      s1vec(1) = a11*svec(1)+a12*svec(2)+a13*svec(3)
      s1vec(2) = a21*svec(1)+a22*svec(2)+a23*svec(3)
      s1vec(3) = a31*svec(1)+a32*svec(2)+a33*svec(3)
      s1vec(4) = c44*svec(4)
      s1vec(5) = c55*svec(5)
      s1vec(6) = c66*svec(6)
      !
      s2vec(1) = b11*svec(1)+b12*svec(2)+b13*svec(3)
      s2vec(2) = b21*svec(1)+b22*svec(2)+b23*svec(3)
      s2vec(3) = b31*svec(1)+b32*svec(2)+b33*svec(3)
      s2vec(4) = d44*svec(4)
      s2vec(5) = d55*svec(5)
      s2vec(6) = d66*svec(6)
      !
      s3vec(1) = z11*svec(1)+z12*svec(2)+z13*svec(3)
      s3vec(2) = z21*svec(1)+z22*svec(2)+z23*svec(3)
      s3vec(3) = z31*svec(1)+z32*svec(2)+z33*svec(3)
      s3vec(4) = e44*svec(4)
      s3vec(5) = e55*svec(5)
      s3vec(6) = e66*svec(6)
      
c     Calculate eigenvalues and eigenvectors of transformed stresses
      Smat = reshape(
     +     (/ s1vec(1), s1vec(6), s1vec(5),
     +        s1vec(6), s1vec(2), s1vec(4),
     +        s1vec(5), s1vec(4), s1vec(3) /), (/3,3/))
      call eig2(Smat,0,s11,s12,s13,eigvec11,eigvec12,eigvec13)
c
      Smat = reshape(
     +     (/ s2vec(1), s2vec(6), s2vec(5),
     +        s2vec(6), s2vec(2), s2vec(4),
     +        s2vec(5), s2vec(4), s2vec(3) /), (/3,3/))
      call eig2(Smat,0,s21,s22,s23,eigvec21,eigvec22,eigvec23)
c
      Smat = reshape(
     +     (/ s3vec(1), s3vec(6), s3vec(5),
     +        s3vec(6), s3vec(2), s3vec(4),
     +        s3vec(5), s3vec(4), s3vec(3) /), (/3,3/))
      call eig2(Smat,0,s31,s32,s33,eigvec21,eigvec22,eigvec23)
c      
c     
      asms(1) = abs(s11-s21)
      asms(2) = abs(s11-s22)
      asms(3) = abs(s11-s23)
c              
      asms(4) = abs(s12-s21)
      asms(5) = abs(s12-s22)
      asms(6) = abs(s12-s23)
c              
      asms(7) = abs(s13-s21)
      asms(8) = abs(s13-s22)
      asms(9) = abs(s13-s23)
c
      asms(10) = abs(s31-s32)
      asms(11) = abs(s31-s33)
      asms(12) = abs(s32-s33)
c
      scale = asms(1)
      do i=2, 12
          if (asms(i) .GT. scale) scale = asms(i)
      enddo
c     
      f = 0.d0
      if (scale .GT. 1.d-16) then
          do i=1, 12
              f = f + (asms(i)/scale)**a
          enddo
c
c         Compute yield function F
          f = scale*(f/6.d0)**ova
      end if
c
      return    
      end subroutine Yld2004_27p
c
c
      subroutine Yld2011_18p(svec, cparams, a, f)
c**********************************************************************
c      Calculates Aretz's Yld2011-18p yield function
c**********************************************************************
c (IN)    SVEC is REAL*8 array, dimension (6)
c         It is stress tensor given in Voigt notation
c (IN)    CPARAMS is REAL*8 array, dimension (18)
c         It is 18 anisotropic parameters of Yld2011 expected in Voigt
c         notation, i.e. the first L anisotropic transformation matrix
c         read
c         |  0  -C12 -C13  0   0   0 |         
c         |-C21   0  -C23  0   0   0 |
c         |-C31 -C32   0   0   0   0 |
c	    |  0    0    0  C44  0   0 |
c         |  0    0    0   0  C55  0 |
c         |  0    0    0   0   0  C66|
c
c         Coefficients in LPARAMS are ordered as
c         C'12, C'13, C'21, C'23, C'31, C'32, C'44, C'55, C'66, 
c         followed by coefficients associated with the second anistropic 
c         transformation matrix L''
c (IN)    A is REAL*8
c         It is the exponent in Yld2011 (a must be >= 2)
c (OUT)   F is REAL*8
c         It is scalar value of Yld2011 yield function
c         
c**********************************************************************
      implicit none
c     ! in/out
      real*8      svec(6), cparams(18), a, f
c     ! internal
      integer     i, j, emod
      real*8      s1vec(6), s2vec(6), scale, asms(9), xi, ov3, ova, ov3a
c
      real*8      c12, c13, c21, c23, c31, c32, c44, c55, c66,          
     +            d12, d13, d21, d23, d31, d32, d44, d55, d66, 
     +            a11, a12, a13, a21, a22, a23, a31, a32, a33,
     +            b11, b12, b13, b21, b22, b23, b31, b32, b33, 
     +            s11, s12, s13, s21, s22, s23, Smat(3,3),
     +            eigvec11(3), eigvec12(3), eigvec13(3),                
     +            eigvec21(3), eigvec22(3), eigvec23(3)
c
Cf2py intent(out) f
c
      ov3     = 1.d0/3.d0
	ova     = 1.d0/a
      ov3a    = ov3**a
c      
      c12  = cparams(1)
      c13  = cparams(2)
      c21  = cparams(3)
      c23  = cparams(4)
      c31  = cparams(5)
      c32  = cparams(6)
      c44  = cparams(7)
      c55  = cparams(8)
      c66  = cparams(9)
      d12  = cparams(10)
      d13  = cparams(11)
      d21  = cparams(12)
      d23  = cparams(13)
      d31  = cparams(14)
      d32  = cparams(15)
      d44  = cparams(16)
      d55  = cparams(17)
      d66  = cparams(18)
c
      a11 = ov3*(c12+c13)
      a12 = ov3*(c13-2.0d0*c12)
      a13 = ov3*(c12-2.0d0*c13)
      a21 = ov3*(c23-2.0d0*c21)
      a22 = ov3*(c21+c23)
      a23 = ov3*(c21-2.0d0*c23)
      a31 = ov3*(c32-2.0d0*c31)
      a32 = ov3*(c31-2.0d0*c32)
      a33 = ov3*(c31+c32)
      !
      b11 = ov3*(d12+d13)
      b12 = ov3*(d13-2.0d0*d12)
      b13 = ov3*(d12-2.0d0*d13)
      b21 = ov3*(d23-2.0d0*d21)
      b22 = ov3*(d21+d23)
      b23 = ov3*(d21-2.0d0*d23)
      b31 = ov3*(d32-2.0d0*d31)
      b32 = ov3*(d31-2.0d0*d32)
      b33 = ov3*(d31+d32)
      !
      s1vec(1) = a11*svec(1)+a12*svec(2)+a13*svec(3)
      s1vec(2) = a21*svec(1)+a22*svec(2)+a23*svec(3)
      s1vec(3) = a31*svec(1)+a32*svec(2)+a33*svec(3)
      s1vec(4) = c44*svec(4)
      s1vec(5) = c55*svec(5)
      s1vec(6) = c66*svec(6)
      !
      s2vec(1) = b11*svec(1)+b12*svec(2)+b13*svec(3)
      s2vec(2) = b21*svec(1)+b22*svec(2)+b23*svec(3)
      s2vec(3) = b31*svec(1)+b32*svec(2)+b33*svec(3)
      s2vec(4) = d44*svec(4)
      s2vec(5) = d55*svec(5)
      s2vec(6) = d66*svec(6)
      
c     Calculate eigenvalues and eigenvectors of transformed stresses
      Smat = reshape(
     +     (/ s1vec(1), s1vec(6), s1vec(5),
     +        s1vec(6), s1vec(2), s1vec(4),
     +        s1vec(5), s1vec(4), s1vec(3) /), (/3,3/))
      call eig2(Smat,0,s11,s12,s13,eigvec11,eigvec12,eigvec13)
c
      Smat = reshape(
     +     (/ s2vec(1), s2vec(6), s2vec(5),
     +        s2vec(6), s2vec(2), s2vec(4),
     +        s2vec(5), s2vec(4), s2vec(3) /), (/3,3/))
      call eig2(Smat,0,s21,s22,s23,eigvec21,eigvec22,eigvec23)
c      
c     
      asms(1) = abs(s11+s21)
      asms(2) = abs(s11+s22)
      asms(3) = abs(s11+s23)
c              
      asms(4) = abs(s12+s21)
      asms(5) = abs(s12+s22)
      asms(6) = abs(s12+s23)
c              
      asms(7) = abs(s13+s21)
      asms(8) = abs(s13+s22)
      asms(9) = abs(s13+s23)
c
      scale = asms(1)
      do i=2, 9
          if (asms(i) .GT. scale) scale = asms(i)
      enddo
c     
      f = 0.d0
      if (scale .GT. 1.d-16) then
          do i=1, 9
              f = f + (asms(i)/scale)**a
          enddo
c
c         Compute yield function F
          xi = ov3a*(4.d0**a + 4.d0*(2.d0**a + 1.d0))
          f = scale*(f/xi)**ova
      end if
c
      return    
      end subroutine Yld2011_18p
c
c
      subroutine Yld2011_27p(svec, cparams, a, f)
c**********************************************************************
c      Calculates Barlat's Yld2011-27p yield function
c**********************************************************************
c (IN)    SVEC is REAL*8 array, dimension (6)
c         It is stress tensor given in Voigt notation
c (IN)    CPARAMS is REAL*8 array, dimension (27)
c         It is 27 anisotropic parameters of Yld2011 expected in Voigt
c         notation, i.e. the first L anisotropic transformation matrix
c         read
c         |  0  -C12 -C13  0   0   0 |         
c         |-C21   0  -C23  0   0   0 |
c         |-C31 -C32   0   0   0   0 |
c	    |  0    0    0  C44  0   0 |
c         |  0    0    0   0  C55  0 |
c         |  0    0    0   0   0  C66|
c
c         Coefficients in LPARAMS are ordered as
c         C'12, C'13, C'21, C'23, C'31, C'32, C'44, C'55, C'66, 
c         followed by coefficients associated with the second and third
c         anistropic transformation matrices L'' and L'''
c (IN)    A is REAL*8
c         It is the exponent in Yld2011 (a must be >= 2)
c (OUT)   F is REAL*8
c         It is scalar value of Yld2011 yield function
c         
c**********************************************************************
      implicit none
c     ! in/out
      real*8      svec(6), cparams(27), a, f
c     ! internal
      integer     i, j, emod
      real*8      s1vec(6), s2vec(6), s3vec(6), scale, asms(12), 
     +            Smat(3,3), xi, ov3, ova, ov3a
c
      real*8      c12, c13, c21, c23, c31, c32, c44, c55, c66,          
     +            d12, d13, d21, d23, d31, d32, d44, d55, d66, 
     +            e12, e13, e21, e23, e31, e32, e44, e55, e66, 
     +            a11, a12, a13, a21, a22, a23, a31, a32, a33,
     +            b11, b12, b13, b21, b22, b23, b31, b32, b33, 
     +            z11, z12, z13, z21, z22, z23, z31, z32, z33, 
     +            s11, s12, s13, s21, s22, s23, s31, s32, s33,
     +            eigvec11(3), eigvec12(3), eigvec13(3),                
     +            eigvec21(3), eigvec22(3), eigvec23(3)
c
Cf2py intent(out) f
c
      ov3     = 1.d0/3.d0
	ova     = 1.d0/a
      ov3a    = ov3**a
c      
      c12  = cparams(1)
      c13  = cparams(2)
      c21  = cparams(3)
      c23  = cparams(4)
      c31  = cparams(5)
      c32  = cparams(6)
      c44  = cparams(7)
      c55  = cparams(8)
      c66  = cparams(9)
      d12  = cparams(10)
      d13  = cparams(11)
      d21  = cparams(12)
      d23  = cparams(13)
      d31  = cparams(14)
      d32  = cparams(15)
      d44  = cparams(16)
      d55  = cparams(17)
      d66  = cparams(18)
      e12  = cparams(19)
      e13  = cparams(20)
      e21  = cparams(21)
      e23  = cparams(22)
      e31  = cparams(23)
      e32  = cparams(24)
      e44  = cparams(25)
      e55  = cparams(26)
      e66  = cparams(27)
c
      a11 = ov3*(c12+c13)
      a12 = ov3*(c13-2.0d0*c12)
      a13 = ov3*(c12-2.0d0*c13)
      a21 = ov3*(c23-2.0d0*c21)
      a22 = ov3*(c21+c23)
      a23 = ov3*(c21-2.0d0*c23)
      a31 = ov3*(c32-2.0d0*c31)
      a32 = ov3*(c31-2.0d0*c32)
      a33 = ov3*(c31+c32)
      !
      b11 = ov3*(d12+d13)
      b12 = ov3*(d13-2.0d0*d12)
      b13 = ov3*(d12-2.0d0*d13)
      b21 = ov3*(d23-2.0d0*d21)
      b22 = ov3*(d21+d23)
      b23 = ov3*(d21-2.0d0*d23)
      b31 = ov3*(d32-2.0d0*d31)
      b32 = ov3*(d31-2.0d0*d32)
      b33 = ov3*(d31+d32)
      !
      z11 = ov3*(e12+e13)
      z12 = ov3*(e13-2.0d0*e12)
      z13 = ov3*(e12-2.0d0*e13)
      z21 = ov3*(e23-2.0d0*e21)
      z22 = ov3*(e21+e23)
      z23 = ov3*(e21-2.0d0*e23)
      z31 = ov3*(e32-2.0d0*e31)
      z32 = ov3*(e31-2.0d0*e32)
      z33 = ov3*(e31+e32)
      !
      s1vec(1) = a11*svec(1)+a12*svec(2)+a13*svec(3)
      s1vec(2) = a21*svec(1)+a22*svec(2)+a23*svec(3)
      s1vec(3) = a31*svec(1)+a32*svec(2)+a33*svec(3)
      s1vec(4) = c44*svec(4)
      s1vec(5) = c55*svec(5)
      s1vec(6) = c66*svec(6)
      !
      s2vec(1) = b11*svec(1)+b12*svec(2)+b13*svec(3)
      s2vec(2) = b21*svec(1)+b22*svec(2)+b23*svec(3)
      s2vec(3) = b31*svec(1)+b32*svec(2)+b33*svec(3)
      s2vec(4) = d44*svec(4)
      s2vec(5) = d55*svec(5)
      s2vec(6) = d66*svec(6)
      !
      s3vec(1) = z11*svec(1)+z12*svec(2)+z13*svec(3)
      s3vec(2) = z21*svec(1)+z22*svec(2)+z23*svec(3)
      s3vec(3) = z31*svec(1)+z32*svec(2)+z33*svec(3)
      s3vec(4) = e44*svec(4)
      s3vec(5) = e55*svec(5)
      s3vec(6) = e66*svec(6)
      
c     Calculate eigenvalues and eigenvectors of transformed stresses
      Smat = reshape(
     +     (/ s1vec(1), s1vec(6), s1vec(5),
     +        s1vec(6), s1vec(2), s1vec(4),
     +        s1vec(5), s1vec(4), s1vec(3) /), (/3,3/))
      call eig2(Smat,0,s11,s12,s13,eigvec11,eigvec12,eigvec13)
c
      Smat = reshape(
     +     (/ s2vec(1), s2vec(6), s2vec(5),
     +        s2vec(6), s2vec(2), s2vec(4),
     +        s2vec(5), s2vec(4), s2vec(3) /), (/3,3/))
      call eig2(Smat,0,s21,s22,s23,eigvec21,eigvec22,eigvec23)
c
      Smat = reshape(
     +     (/ s3vec(1), s3vec(6), s3vec(5),
     +        s3vec(6), s3vec(2), s3vec(4),
     +        s3vec(5), s3vec(4), s3vec(3) /), (/3,3/))
      call eig2(Smat,0,s31,s32,s33,eigvec21,eigvec22,eigvec23)
c      
c     
      asms(1) = abs(s11+s21)
      asms(2) = abs(s11+s22)
      asms(3) = abs(s11+s23)
c              
      asms(4) = abs(s12+s21)
      asms(5) = abs(s12+s22)
      asms(6) = abs(s12+s23)
c              
      asms(7) = abs(s13+s21)
      asms(8) = abs(s13+s22)
      asms(9) = abs(s13+s23)
c
      asms(10) = abs(s31)
      asms(11) = abs(s32)
      asms(12) = abs(s33)
c
      scale = asms(1)
      do i=2, 12
          if (asms(i) .GT. scale) scale = asms(i)
      enddo
c     
      f = 0.d0
      if (scale .GT. 1.d-16) then
          do i=1, 12
              f = f + (asms(i)/scale)**a
          enddo
c
c         Compute yield function F
          xi = ov3a*(4.d0**a + 5.d0*2.d0**a + 6.d0)
          f = scale*(f/xi)**ova
      end if
c
      return    
      end subroutine Yld2011_27p
c
c
      subroutine Yld2013_27p(svec, cparams, a, f)
c**********************************************************************
c      Calculates Aretz's alternative Yld2011-27p yield function
c**********************************************************************
c (IN)    SVEC is REAL*8 array, dimension (6)
c         It is stress tensor given in Voigt notation
c (IN)    CPARAMS is REAL*8 array, dimension (27)
c         It is 27 anisotropic parameters of Yld2011 expected in Voigt
c         notation, i.e. the first L anisotropic transformation matrix
c         read
c         |  0  -C12 -C13  0   0   0 |         
c         |-C21   0  -C23  0   0   0 |
c         |-C31 -C32   0   0   0   0 |
c	    |  0    0    0  C44  0   0 |
c         |  0    0    0   0  C55  0 |
c         |  0    0    0   0   0  C66|
c
c         Coefficients in LPARAMS are ordered as
c         C'12, C'13, C'21, C'23, C'31, C'32, C'44, C'55, C'66, 
c         followed by coefficients associated with the second and third
c         anistropic transformation matrices L'' and L'''
c (IN)    A is REAL*8
c         It is the exponent in Yld2011 (a must be >= 2)
c (OUT)   F is REAL*8
c         It is scalar value of Yld2011 yield function
c         
c**********************************************************************
      implicit none
c     ! in/out
      real*8      svec(6), cparams(27), a, f
c     ! internal
      integer     i, j, emod
      real*8      s1vec(6), s2vec(6), s3vec(6), scale, asms(27), 
     +            Smat(3,3), xi, ov3, ova
c
      real*8      c12, c13, c21, c23, c31, c32, c44, c55, c66,          
     +            d12, d13, d21, d23, d31, d32, d44, d55, d66, 
     +            e12, e13, e21, e23, e31, e32, e44, e55, e66, 
     +            a11, a12, a13, a21, a22, a23, a31, a32, a33,
     +            b11, b12, b13, b21, b22, b23, b31, b32, b33, 
     +            z11, z12, z13, z21, z22, z23, z31, z32, z33, 
     +            s11, s12, s13, s21, s22, s23, s31, s32, s33,
     +            eigvec11(3), eigvec12(3), eigvec13(3),                
     +            eigvec21(3), eigvec22(3), eigvec23(3)
c
Cf2py intent(out) f
c
      ov3     = 1.d0/3.d0
	ova     = 1.d0/a
c      
      c12  = cparams(1)
      c13  = cparams(2)
      c21  = cparams(3)
      c23  = cparams(4)
      c31  = cparams(5)
      c32  = cparams(6)
      c44  = cparams(7)
      c55  = cparams(8)
      c66  = cparams(9)
      d12  = cparams(10)
      d13  = cparams(11)
      d21  = cparams(12)
      d23  = cparams(13)
      d31  = cparams(14)
      d32  = cparams(15)
      d44  = cparams(16)
      d55  = cparams(17)
      d66  = cparams(18)
      e12  = cparams(19)
      e13  = cparams(20)
      e21  = cparams(21)
      e23  = cparams(22)
      e31  = cparams(23)
      e32  = cparams(24)
      e44  = cparams(25)
      e55  = cparams(26)
      e66  = cparams(27)
c
      a11 = ov3*(c12+c13)
      a12 = ov3*(c13-2.0d0*c12)
      a13 = ov3*(c12-2.0d0*c13)
      a21 = ov3*(c23-2.0d0*c21)
      a22 = ov3*(c21+c23)
      a23 = ov3*(c21-2.0d0*c23)
      a31 = ov3*(c32-2.0d0*c31)
      a32 = ov3*(c31-2.0d0*c32)
      a33 = ov3*(c31+c32)
      !
      b11 = ov3*(d12+d13)
      b12 = ov3*(d13-2.0d0*d12)
      b13 = ov3*(d12-2.0d0*d13)
      b21 = ov3*(d23-2.0d0*d21)
      b22 = ov3*(d21+d23)
      b23 = ov3*(d21-2.0d0*d23)
      b31 = ov3*(d32-2.0d0*d31)
      b32 = ov3*(d31-2.0d0*d32)
      b33 = ov3*(d31+d32)
      !
      z11 = ov3*(e12+e13)
      z12 = ov3*(e13-2.0d0*e12)
      z13 = ov3*(e12-2.0d0*e13)
      z21 = ov3*(e23-2.0d0*e21)
      z22 = ov3*(e21+e23)
      z23 = ov3*(e21-2.0d0*e23)
      z31 = ov3*(e32-2.0d0*e31)
      z32 = ov3*(e31-2.0d0*e32)
      z33 = ov3*(e31+e32)
      !
      s1vec(1) = a11*svec(1)+a12*svec(2)+a13*svec(3)
      s1vec(2) = a21*svec(1)+a22*svec(2)+a23*svec(3)
      s1vec(3) = a31*svec(1)+a32*svec(2)+a33*svec(3)
      s1vec(4) = c44*svec(4)
      s1vec(5) = c55*svec(5)
      s1vec(6) = c66*svec(6)
      !
      s2vec(1) = b11*svec(1)+b12*svec(2)+b13*svec(3)
      s2vec(2) = b21*svec(1)+b22*svec(2)+b23*svec(3)
      s2vec(3) = b31*svec(1)+b32*svec(2)+b33*svec(3)
      s2vec(4) = d44*svec(4)
      s2vec(5) = d55*svec(5)
      s2vec(6) = d66*svec(6)
      !
      s3vec(1) = z11*svec(1)+z12*svec(2)+z13*svec(3)
      s3vec(2) = z21*svec(1)+z22*svec(2)+z23*svec(3)
      s3vec(3) = z31*svec(1)+z32*svec(2)+z33*svec(3)
      s3vec(4) = e44*svec(4)
      s3vec(5) = e55*svec(5)
      s3vec(6) = e66*svec(6)
      
c     Calculate eigenvalues and eigenvectors of transformed stresses
      Smat = reshape(
     +     (/ s1vec(1), s1vec(6), s1vec(5),
     +        s1vec(6), s1vec(2), s1vec(4),
     +        s1vec(5), s1vec(4), s1vec(3) /), (/3,3/))
      call eig2(Smat,0,s11,s12,s13,eigvec11,eigvec12,eigvec13)
c
      Smat = reshape(
     +     (/ s2vec(1), s2vec(6), s2vec(5),
     +        s2vec(6), s2vec(2), s2vec(4),
     +        s2vec(5), s2vec(4), s2vec(3) /), (/3,3/))
      call eig2(Smat,0,s21,s22,s23,eigvec21,eigvec22,eigvec23)
c
      Smat = reshape(
     +     (/ s3vec(1), s3vec(6), s3vec(5),
     +        s3vec(6), s3vec(2), s3vec(4),
     +        s3vec(5), s3vec(4), s3vec(3) /), (/3,3/))
      call eig2(Smat,0,s31,s32,s33,eigvec21,eigvec22,eigvec23)
c      
c     
      asms(1)  = abs(s11+s21+s31)
      asms(2)  = abs(s11+s21+s32)
      asms(3)  = abs(s11+s21+s33)    
      asms(4)  = abs(s11+s22+s31)
      asms(5)  = abs(s11+s22+s32)
      asms(6)  = abs(s11+s22+s33)     
      asms(7)  = abs(s11+s23+s31)
      asms(8)  = abs(s11+s23+s32)
      asms(9)  = abs(s11+s23+s33)
      asms(10) = abs(s12+s21+s31)
      asms(11) = abs(s12+s21+s32)
      asms(12) = abs(s12+s21+s33)
      asms(13) = abs(s12+s22+s31)
      asms(14) = abs(s12+s22+s32)
      asms(15) = abs(s12+s22+s33)    
      asms(16) = abs(s12+s23+s31)
      asms(17) = abs(s12+s23+s32)
      asms(18) = abs(s12+s23+s33)     
      asms(19) = abs(s13+s21+s31)
      asms(20) = abs(s13+s21+s32)
      asms(21) = abs(s13+s21+s33)
      asms(22) = abs(s13+s22+s31)
      asms(23) = abs(s13+s22+s32)
      asms(24) = abs(s13+s22+s33)
      asms(25) = abs(s13+s23+s31)
      asms(26) = abs(s13+s23+s32)
      asms(27) = abs(s13+s23+s33)
c
      scale = asms(1)
      do i=2, 27
          if (asms(i) .GT. scale) scale = asms(i)
      enddo
c     
      f = 0.d0
      if (scale .GT. 1.d-16) then
          do i=1, 27
              f = f + (asms(i)/scale)**a
          enddo
c
c         Compute yield function F
          xi = 2.d0**a + 14.d0
          f = scale*(f/xi)**ova
      end if
c
      return    
      end subroutine Yld2013_27p
c
c
      subroutine Yld2013_Xp(svec, cparams, a, NT, f)
c**********************************************************************
c      Calculates Aretz's alternative Yld2013-Xp yield function
c**********************************************************************
c (IN)    SVEC is REAL*8 array, dimension (6)
c         It is stress tensor given in Voigt notation
c (IN)    CPARAMS is REAL*8 array, dimension (200)
c         It is up to 200 anisotropic parameters of Yld2013-Xp with up to
c         20 linear transformations, expected in Voigt
c         notation, i.e. the first L anisotropic transformation matrix
c         read
c         |  0  -C12 -C13  0   0   0 |         
c         |-C21   0  -C23  0   0   0 |
c         |-C31 -C32   0   0   0   0 |
c	    |  0    0    0  C44  0   0 |
c         |  0    0    0   0  C55  0 |
c         |  0    0    0   0   0  C66|
c
c         Coefficients in LPARAMS are ordered as
c         C'12, C'13, C'21, C'23, C'31, C'32, C'44, C'55, C'66, 
c         followed by coefficients associated with the second and third
c         anistropic transformation matrices L'' and L''' etc. 
c (IN)    A is REAL*8
c         It is the exponent in Yld2013 (a must be >= 2)
c (IN)    NT is INTEGER
c         It is the number of linear tranformations to be used
c (OUT)   F is REAL*8
c         It is scalar value of Yld2013 yield function
c         
c**********************************************************************
      implicit none
c     ! in/out
      integer     NT
      real*8      svec(6), cparams(200), a, f
c     ! internal
      integer     i, k, n, j, Nsum, emod, m(16)
      real*8      sTvec(6), Smat(3,3), ov3, ova, ov3a, xi,
     +            c12, c13, c21, c23, c31, c32, c44, c55, c66,
     +            a11, a12, a13, a21, a22, a23, a31, a32, a33,
     +            ST(3,16), eigvec1(3), eigvec2(3), eigvec3(3),
     +            summ
c
Cf2py intent(out) f
c
      ov3     = 1.d0/3.d0
      ova     = 1.d0/a
      ov3a    = ov3**a
      Nsum    = 3**NT
c
      do n = 1, NT
          c12  = cparams(9*n-8)
          c13  = cparams(9*n-7)
          c21  = cparams(9*n-6)
          c23  = cparams(9*n-5)
          c31  = cparams(9*n-4)
          c32  = cparams(9*n-3)
          c44  = cparams(9*n-2)
          c55  = cparams(9*n-1)
          c66  = cparams(9*n)
c
          a11 = ov3*(c12+c13)
          a12 = ov3*(c13-2.0d0*c12)
          a13 = ov3*(c12-2.0d0*c13)
          a21 = ov3*(c23-2.0d0*c21)
          a22 = ov3*(c21+c23)
          a23 = ov3*(c21-2.0d0*c23)
          a31 = ov3*(c32-2.0d0*c31)
          a32 = ov3*(c31-2.0d0*c32)
          a33 = ov3*(c31+c32)
c
          sTvec(1) = a11*svec(1)+a12*svec(2)+a13*svec(3)
          sTvec(2) = a21*svec(1)+a22*svec(2)+a23*svec(3)
          sTvec(3) = a31*svec(1)+a32*svec(2)+a33*svec(3)
          sTvec(4) = c44*svec(4)
          sTvec(5) = c55*svec(5)
          sTvec(6) = c66*svec(6)
c     
c         Calculate eigenvalues and eigenvectors of transformed stresses
          Smat = reshape(
     +         (/ sTvec(1), sTvec(6), sTvec(5),
     +            sTvec(6), sTvec(2), sTvec(4),
     +            sTvec(5), sTvec(4), sTvec(3) /), (/3,3/))
          call eig2(Smat, 0, ST(1,n), ST(2,n), ST(3,n), 
     +              eigvec1,eigvec2,eigvec3)
      end do
c      
c     
      m = 1
      f = 0.d0
      do j=1, Nsum
          summ = 0.d0
          do n = 1, NT
              summ = summ + ST(m(n),n)
          end do
          f = f + (abs(summ))**a 
          if (minval(m(1:NT)) .EQ. 3) exit

          if (m(1) .LT. 3) then
              m(1) = m(1) + 1
          else
              m(1) = 1
              k = 2
              do
                  if (m(k) .LT. 3) then
                      m(k) = m(k) + 1
                      exit
                  else
                      m(k) = 1
                      k = k + 1
                  end if
              end do
          end if
      end do
c
c     Compute yield function F
      select case(NT)
      case(1)
          xi = ov3a*(2.d0**a + 2.d0)
      case(2)
          xi = ov3a*(4.d0**a + 4.d0*2.d0**a + 4.d0)
      case(3)  
          xi = 2.d0**a + 14.d0
      case(4)  
          xi = ov3a*(8.d0**a + 16.d0*4.d0**a + 8.d0*5.d0**a 
     +       + 24.d0*2.d0**a + 32.d0)
      case(5)
          xi = ov3a*(10.d0**a + 40.d0*4.d0**a + 80.d0*2.d0**a
     +       + 32.d0*5.d0**a + 10.d0*7.d0**a + 80.d0)
      case(6)
          xi = 4.d0**a + 12.d0*3.d0**a + 352.d0 + 124.d0*2.d0**a
      case(7)
          xi = ov3a*(14.d0*11.d0**a + 14.d0**a + 84.d0*8.d0**a 
     +       + 448.d0*4.d0**a + 560.d0*2.d0**a + 280.d0*5.d0**a 
     +       + 128.d0*7.d0**a + 672.d0)
      case(8)
          xi = ov3a*(112.d0*10.d0**a + 16.d0*13.d0**a + 16.d0**a
     +       + 256.d0*8.d0**a + 1120.d0*4.d0**a + 1792.d0*2.d0**a
     +       + 1024.d0*5.d0**a + 448.d0*7.d0**a + 1792.d0)
      case(9)
          xi = 144.d0*4.d0**a + 18.d0*5.d0**a + 6.d0**a + 8640.d0 
     +       + 4320.d0*2.d0**a + 1184.d0*3.d0**a
      case(10)
          xi = ov3a*(1024.d0*10.d0**a + 960.d0*11.d0**a 
     +       + 180.d0*14.d0**a + 20.d0*17.d0**a + 3360.d0*8.d0**a 
     +       + 11520.d0*4.d0**a + 13440.d0*2.d0**a + 20.d0**a 
     +       + 8064.d0*5.d0**a + 5120.d0*7.d0**a + 15360.d0)
      case(11) 
          xi = ov3a*(5280.d0*10.d0**a + 2048.d0*11.d0**a 
     +       + 1320.d0*13.d0**a + 22.d0*19.d0**a + 220.d0*16.d0**a 
     +       + 11264.d0*8.d0**a + 29568.d0*4.d0**a + 42240.d0*2.d0**a
     +       + 22.d0**a + 28160.d0*5.d0**a + 14784.d0*7.d0**a
     +       + 42240.d0)   
      case(12)
          xi = 8.d0**a + 24.d0*7.d0**a + 264.d0*6.d0**a 
     +       + 1760.d0*5.d0**a + 12016.d0*4.d0**a + 49920.d0*3.d0**a
     +       + 126720.d0*2.d0**a + 214016.d0
      case(16)
          xi = ov3a*(1966080.d0*10.d0**a + 1464320.d0*11.d0**a 
     +       + 524288.d0*13.d0**a + 512512.d0*14.d0**a 
     +       + 139776.d0*17.d0**a + 65536.d0*16.d0**a 
     +       + 3294720.d0*8.d0**a + 7454720.d0*4.d0**a 
     +       + 8200192.d0*2.d0**a + 29120.d0*20.d0**a 
     +       + 4480.d0*23.d0**a + 480.d0*26.d0**a + 32.d0*29.d0**a
     +       + 32.d0**a + 5857280.d0*5.d0**a + 4587520.d0*7.d0**a 
     +       + 8945664.d0)
      end select
         
      f = (f/xi)**ova
c
      return    
      end subroutine Yld2013_Xp
c
c
      subroutine Yld2004nort_21p(svec, cparams, a, f)
c**********************************************************************
c      Calculates a non-orthotropic Yld2004-type of yield function
c      with 1-linear transformation and 21 parameters
c**********************************************************************
c (IN)    SVEC is REAL*8 array, dimension (6)
c         It is stress tensor given in Voigt notation
c (IN)    CPARAMS is REAL*8 array, dimension (21)
c         It is 21 anisotropic parameters of Yld2004 expected in Voigt
c         notation, i.e. the C anisotropic transformation matrix
c         read
c         |  0  -C12 -C13 C14 C15 C16|         
c         |-C21   0  -C23 C24 C25 C26|
c         |-C31 -C32   0  C34 C35 C36|
c         | C14  C24  C34 C44 C45 C46|
c         | C15  C25  C35 C45 C55 C56|
c         | C16  C26  C36 C46 C56 C66|
c
c (IN)    A is REAL*8
c         It is the exponent in Yld2004 (a must be >= 2)
c (OUT)   F is REAL*8
c         It is scalar value of Yld2004 yield function
c         
c**********************************************************************
      implicit none
c     ! in/out
      real*8      svec(6), cparams(21), a, f
c     ! internal
      integer     i, j, emod
      real*8      s1vec(6), scale, asms(3), ov3, ova
c
      real*8      c12, c13, c21, c23, c31, c32, c44, c55, c66,
     +            c14, c15, c16, c24, c25, c26, c34, c35, c36,
     +            c45, c46, c56,
     +            a11, a12, a13, a21, a22, a23, a31, a32, a33,
     +            a41, a42, a43, a51, a52, a53, a61, a62, a63,
     +            s11, s12, s13, Smat(3,3),
     +            eigvec11(3), eigvec12(3), eigvec13(3)
c
Cf2py intent(out) f
c
      ov3     = 1.d0/3.d0
      ova     = 1.d0/a
c      
      c12  = cparams(1)
      c13  = cparams(2)
      c21  = cparams(3)
      c23  = cparams(4)
      c31  = cparams(5)
      c32  = cparams(6)
      c44  = cparams(7)
      c55  = cparams(8)
      c66  = cparams(9)
      c14  = cparams(10)
      c15  = cparams(11)
      c16  = cparams(12)
      c24  = cparams(13)
      c25  = cparams(14)
      c26  = cparams(15)
      c34  = cparams(16)
      c35  = cparams(17)
      c36  = cparams(18)
      c45  = cparams(19)
      c46  = cparams(20)
      c56  = cparams(21)
c
      a11 = ov3*(c12+c13)
      a12 = ov3*(c13-2.0d0*c12)
      a13 = ov3*(c12-2.0d0*c13)
      a21 = ov3*(c23-2.0d0*c21)
      a22 = ov3*(c21+c23)
      a23 = ov3*(c21-2.0d0*c23)
      a31 = ov3*(c32-2.0d0*c31)
      a32 = ov3*(c31-2.0d0*c32)
      a33 = ov3*(c31+c32)
c
      a41 = ov3*(2.0d0*c14-c24-c34)
      a51 = ov3*(2.0d0*c15-c25-c35)
      a61 = ov3*(2.0d0*c16-c26-c36)
      a42 = ov3*(-c14+2.0d0*c24-c34)
      a52 = ov3*(-c15+2.0d0*c25-c35)
      a62 = ov3*(-c16+2.0d0*c26-c36)
      a43 = ov3*(-c14-c24+2.0d0*c34)
      a53 = ov3*(-c15-c25+2.0d0*c35)
      a63 = ov3*(-c16-c26+2.0d0*c36)
c
      s1vec(1) = a11*svec(1)+a12*svec(2)+a13*svec(3)
     +         + c14*svec(4)+c15*svec(5)+c16*svec(6)
      s1vec(2) = a21*svec(1)+a22*svec(2)+a23*svec(3)
     +         + c24*svec(4)+c25*svec(5)+c26*svec(6)
      s1vec(3) = a31*svec(1)+a32*svec(2)+a33*svec(3)
     +         + c34*svec(4)+c35*svec(5)+c36*svec(6)
      s1vec(4) = a41*svec(1)+a42*svec(2)+a43*svec(3)
     +         + c44*svec(4)+c45*svec(5)+c46*svec(6)
      s1vec(5) = a51*svec(1)+a52*svec(2)+a53*svec(3)
     +         + c45*svec(4)+c55*svec(5)+c56*svec(6)
      s1vec(6) = a61*svec(1)+a62*svec(2)+a63*svec(3)
     +         + c46*svec(4)+c56*svec(5)+c66*svec(6)
      
c     Calculate eigenvalues and eigenvectors of transformed stresses
      Smat = reshape(
     +     (/ s1vec(1), s1vec(6), s1vec(5),
     +        s1vec(6), s1vec(2), s1vec(4),
     +        s1vec(5), s1vec(4), s1vec(3) /), (/3,3/))
      call eig2(Smat,0,s11,s12,s13,eigvec11,eigvec12,eigvec13)
c      
c     
      asms(1) = abs(s11-s12)
      asms(2) = abs(s11-s13)
      asms(3) = abs(s12-s13)
c
      scale = asms(1)
      do i=2, 3
          if (asms(i) .GT. scale) scale = asms(i)
      enddo
c     
      f = 0.d0
      if (scale .GT. 1.d-16) then
          do i=1, 3
              f = f + (asms(i)/scale)**a
          enddo
c
c         Compute yield function F
          f = scale*(f/2.d0)**ova
      end if
c
      return    
      end subroutine Yld2004nort_21p
c
c
      subroutine Yld2004_18p_grad(svec, cparams, a, mode, 
     +                                              f, grad, hessian)
c**********************************************************************
c      Calculates Barlat's Yld2004 yield function and/or its gradient
c      and/or its Hessian expressed in Voigt notation. Algorithm is 
c      based on the paper by Scherzinger, W. M.: A return mapping 
c      algorithm for isotropic and anisotropic plasticity models using
c      a line search method, Computer Methods in Applied Mechanics
c      and Engineering, 2017, DOI: 10.1016/j.cma.2016.11.026
c
c**********************************************************************
c (IN)    SVEC is REAL*8 array, dimension (6)
c         It is stress tensor given in Voigt notation
c (IN)    CPARAMS is REAL*8 array, dimension (18)
c         It is 18 anisotropic parameters of Yld2004 expected in Voigt
c         notation, i.e. the first L anisotropic transformation matrix
c         read
c         |  0  -C12 -C13  0   0   0 |         
c         |-C21   0  -C23  0   0   0 |
c         |-C31 -C32   0   0   0   0 |
c	    |  0    0    0  C44  0   0 |
c         |  0    0    0   0  C55  0 |
c         |  0    0    0   0   0  C66|
c
c         Coefficients in CPARAMS are ordered as
c         C'12, C'13, C'21, C'23, C'31, C'32, C'44, C'55, C'66, 
c         followed by coeffiients associated with the second anistropic 
c         transformation matrix L''
c (IN)    A is REAL*8
c         It is the exponent in Yld2004 (a must be >= 2)
c (IN)    MODE is INTEGER
c         Controls which of F, GRAD, HESSIAN is to be calculated
c         = 0  :  only F is calculated
c         = 1  :  only F and GRAD is calculated
c         else :  F, GRAD and HESSIAN are calculated
c (OUT)   F is REAL*8
c         It is scalar value of Yld2004 yield function
c (OUT)   GRAD is REAL*8 array, dimension (6)
c         It is gradient of Yld2004 yield function, computed in 
c         Voigt notation
c (OUT)   HESSIAN is REAL*8 array, dimension (6,6)
c         It is the Hessian of Yld2004 yield function computed in 
c         Voigt notation
c         
c**********************************************************************
      implicit none
c     ! in/out
      integer     mode
      real*8      svec(6), cparams(18), a, f, grad(6), hessian(6,6)
c     ! internal
      integer     i, j, emod
      real*8      s1vec(6), s2vec(6), scale, tol, asms(9), tmp(6),
     +            tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
c
      real*8      c12, c13, c21, c23, c31, c32, c44, c55, c66,          
     +            d12, d13, d21, d23, d31, d32, d44, d55, d66, 
     +            a11, a12, a13, a21, a22, a23, a31, a32, a33,
     +            b11, b12, b13, b21, b22, b23, b31, b32, b33, 
     +            s11, s12, s13, s21, s22, s23,                         
     +            eigvec11(3), eigvec12(3), eigvec13(3),                
     +            eigvec21(3), eigvec22(3), eigvec23(3),                
     +            s11s21, s11s22, s11s23,                               
     +            s12s21, s12s22, s12s23,                               
     +            s13s21, s13s22, s13s23,                               
     +            as11s21pa2, as11s22pa2, as11s23pa2,                   
     +            as12s21pa2, as12s22pa2, as12s23pa2,                   
     +            as13s21pa2, as13s22pa2, as13s23pa2,                   
     +            dfds11, dfds12, dfds13, dfds21, dfds22, dfds23,       
     +            ddfds11ds11, ddfds12ds12, ddfds13ds13,                
     +            ddfds21ds21, ddfds22ds22, ddfds23ds23,                
     +            ddfds11ds12, ddfds11ds13,                             
     +            ddfds12ds11, ddfds12ds13,                             
     +            ddfds13ds11, ddfds13ds12,                             
     +            ddfds21ds22, ddfds21ds23,                             
     +            ddfds22ds21, ddfds22ds23,                             
     +            ddfds23ds21, ddfds23ds22,                             
     +            ddfds11ds21, ddfds12ds22, ddfds13ds23,                
     +            ddfds11ds22, ddfds11ds23,                             
     +            ddfds12ds21, ddfds12ds23,                             
     +            ddfds13ds21, ddfds13ds22,                             
     +            ddfds21ds11, ddfds21ds12,                             
     +            ddfds21ds13, ddfds22ds11,                             
     +            ddfds22ds12, ddfds22ds13,                             
     +            ddfds23ds11, ddfds23ds12, ddfds23ds13,                
     +            e11xe11v(6), e11xe12v(6), e11xe13v(6),                
     +            e12xe11v(6), e12xe12v(6), e12xe13v(6),                
     +            e13xe11v(6), e13xe12v(6), e13xe13v(6),                
     +            e21xe21v(6), e21xe22v(6), e21xe23v(6),                
     +            e22xe21v(6), e22xe22v(6), e22xe23v(6),                
     +            e23xe21v(6), e23xe22v(6), e23xe23v(6),                
     +            ds11xds11(6,6), ds11xds12(6,6), ds11xds13(6,6),       
     +            ds12xds11(6,6), ds12xds12(6,6), ds12xds13(6,6),       
     +            ds13xds11(6,6), ds13xds12(6,6), ds13xds13(6,6),       
     +            ds21xds21(6,6), ds21xds22(6,6), ds21xds23(6,6),       
     +            ds22xds21(6,6), ds22xds22(6,6), ds22xds23(6,6),       
     +            ds23xds21(6,6), ds23xds22(6,6), ds23xds23(6,6),       
     +            ds11xds21(6,6), ds11xds22(6,6), ds11xds23(6,6),       
     +            ds12xds21(6,6), ds12xds22(6,6), ds12xds23(6,6),       
     +            ds13xds21(6,6), ds13xds22(6,6), ds13xds23(6,6),       
     +            ds21xds11(6,6), ds21xds12(6,6), ds21xds13(6,6),       
     +            ds22xds11(6,6), ds22xds12(6,6), ds22xds13(6,6),       
     +            ds23xds11(6,6), ds23xds12(6,6), ds23xds13(6,6),       
     +            Etmp1212(6,6), Etmp1221(6,6), Etmp2112(6,6),          
     +            Etmp2121(6,6), Etmp2323(6,6), Etmp2332(6,6),          
     +            Etmp3223(6,6), Etmp3232(6,6), Etmp3131(6,6),          
     +            Etmp3113(6,6), Etmp1331(6,6), Etmp1313(6,6),          
     +            E1212p(6,6), E2323p(6,6), E3131p(6,6),                
     +            E1212pp(6,6), E2323pp(6,6), E3131pp(6,6),             
     +            coefE1212p, coefE2323p, coefE3131p,                   
     +            coefE1212pp, coefE2323pp, coefE3131pp,                
     +            H1(6,6), H2(6,6), H3(6,6), H4(6,6), H5(6,6), H6(6,6), 
     +            H7(6,6), H8(6,6), dsymLTxH3xM(6,6), Smat(3,3),       
     +            LTxH7xL(6,6), MTxH8xM(6,6)
c
      real*8      ov3, ov4, mov4, a2, a1, a1ovf, ova, ov2
c
      ov2     = 0.5d0
      ov3     = 1.d0/3.d0
      ov4     = 0.25d0
      mov4    = -ov4
      a2      = a-2.d0
      a1      = a-1.d0 
	ova     = 1.d0/a
	tol     = 1.d-5
c      
      c12  = cparams(1)
      c13  = cparams(2)
      c21  = cparams(3)
      c23  = cparams(4)
      c31  = cparams(5)
      c32  = cparams(6)
      c44  = cparams(7)
      c55  = cparams(8)
      c66  = cparams(9)
      d12  = cparams(10)
      d13  = cparams(11)
      d21  = cparams(12)
      d23  = cparams(13)
      d31  = cparams(14)
      d32  = cparams(15)
      d44  = cparams(16)
      d55  = cparams(17)
      d66  = cparams(18)
c
      a11 = ov3*(c12+c13)
      a12 = ov3*(c13-2.0d0*c12)
      a13 = ov3*(c12-2.0d0*c13)
      a21 = ov3*(c23-2.0d0*c21)
      a22 = ov3*(c21+c23)
      a23 = ov3*(c21-2.0d0*c23)
      a31 = ov3*(c32-2.0d0*c31)
      a32 = ov3*(c31-2.0d0*c32)
      a33 = ov3*(c31+c32)
      !
      b11 = ov3*(d12+d13)
      b12 = ov3*(d13-2.0d0*d12)
      b13 = ov3*(d12-2.0d0*d13)
      b21 = ov3*(d23-2.0d0*d21)
      b22 = ov3*(d21+d23)
      b23 = ov3*(d21-2.0d0*d23)
      b31 = ov3*(d32-2.0d0*d31)
      b32 = ov3*(d31-2.0d0*d32)
      b33 = ov3*(d31+d32)
      !
      s1vec(1) = a11*svec(1)+a12*svec(2)+a13*svec(3)
      s1vec(2) = a21*svec(1)+a22*svec(2)+a23*svec(3)
      s1vec(3) = a31*svec(1)+a32*svec(2)+a33*svec(3)
      s1vec(4) = c44*svec(4)
      s1vec(5) = c55*svec(5)
      s1vec(6) = c66*svec(6)
      !
      s2vec(1) = b11*svec(1)+b12*svec(2)+b13*svec(3)
      s2vec(2) = b21*svec(1)+b22*svec(2)+b23*svec(3)
      s2vec(3) = b31*svec(1)+b32*svec(2)+b33*svec(3)
      s2vec(4) = d44*svec(4)
      s2vec(5) = d55*svec(5)
      s2vec(6) = d66*svec(6)
c
      if (mode .EQ. 0) then
c         Calculate only eigenvalues if only F is required
          emod = 0
      else
c         Calculate eigenvalues and eigenvectors for GRAD and HESSIAN
          emod = 1
      end if
      
c     Calculate eigenvalues and eigenvectors of transformed stresses
      Smat = reshape(
     +     (/ s1vec(1), s1vec(6), s1vec(5),
     +        s1vec(6), s1vec(2), s1vec(4),
     +        s1vec(5), s1vec(4), s1vec(3) /), (/3,3/))
      call eig2(Smat,emod,s11,s12,s13,eigvec11,eigvec12,eigvec13)
c
      Smat = reshape(
     +     (/ s2vec(1), s2vec(6), s2vec(5),
     +        s2vec(6), s2vec(2), s2vec(4),
     +        s2vec(5), s2vec(4), s2vec(3) /), (/3,3/))
      call eig2(Smat,emod,s21,s22,s23,eigvec21,eigvec22,eigvec23)
c      
c     
      asms(1) = abs(s11-s21)
      asms(2) = abs(s11-s22)
      asms(3) = abs(s11-s23)
c              
      asms(4) = abs(s12-s21)
      asms(5) = abs(s12-s22)
      asms(6) = abs(s12-s23)
c              
      asms(7) = abs(s13-s21)
      asms(8) = abs(s13-s22)
      asms(9) = abs(s13-s23)
c
      scale = asms(1)
      do i=2, 9
          if (asms(i) .GT. scale) scale = asms(i)
      enddo
c     
      f = 0.d0
      if (scale .GT. 1.d-16) then
          do i=1, 9
              f = f + (asms(i)/scale)**a
          enddo
c
c         Compute yield function F
          f = scale*(ov4*f)**ova
      end if
c
c
      if (mode .EQ. 0) return    
c
	a1ovf   = a1/f
c     rescaling eigenstresses by F
      s11 = s11/f
	s12 = s12/f
	s13 = s13/f
c
      s21 = s21/f
	s22 = s22/f
	s23 = s23/f
c
      s11s21 = s11-s21
      s11s22 = s11-s22
      s11s23 = s11-s23
c      
      s12s21 = s12-s21
      s12s22 = s12-s22
      s12s23 = s12-s23
c
      s13s21 = s13-s21
      s13s22 = s13-s22
      s13s23 = s13-s23
c
      as11s21pa2 = abs(s11s21)**a2
      as11s22pa2 = abs(s11s22)**a2
      as11s23pa2 = abs(s11s23)**a2
c      
      as12s21pa2 = abs(s12s21)**a2
      as12s22pa2 = abs(s12s22)**a2
      as12s23pa2 = abs(s12s23)**a2
c      
      as13s21pa2 = abs(s13s21)**a2
      as13s22pa2 = abs(s13s22)**a2
      as13s23pa2 = abs(s13s23)**a2
cc
      dfds11 = ov4*( s11s21*as11s21pa2 + s11s22*as11s22pa2            
     +        + s11s23*as11s23pa2 )
      dfds12 = ov4*( s12s21*as12s21pa2 + s12s22*as12s22pa2            
     +        + s12s23*as12s23pa2 )                                   
      dfds13 = ov4*( s13s21*as13s21pa2 + s13s22*as13s22pa2            
     +        + s13s23*as13s23pa2 )                                   
cc                                                                    
      dfds21 = ov4*(-s11s21*as11s21pa2 - s12s21*as12s21pa2            
     +        - s13s21*as13s21pa2 )                                   
      dfds22 = ov4*(-s11s22*as11s22pa2 - s12s22*as12s22pa2            
     +        - s13s22*as13s22pa2 )                                   
      dfds23 = ov4*(-s11s23*as11s23pa2 - s12s23*as12s23pa2            
     +        - s13s23*as13s23pa2 )
cc    
      ddfds11ds11 = a1ovf*( ov4*(as11s21pa2 + as11s22pa2 + as11s23pa2)
     +            - dfds11*dfds11 )                                   
      ddfds12ds12 = a1ovf*( ov4*(as12s21pa2 + as12s22pa2 + as12s23pa2)
     +            - dfds12*dfds12 )                                   
      ddfds13ds13 = a1ovf*( ov4*(as13s21pa2 + as13s22pa2 + as13s23pa2)
     +            - dfds13*dfds13 )                                   
      !                                                               
      ddfds21ds21 = a1ovf*( ov4*(as11s21pa2 + as12s21pa2 + as13s21pa2)
     +            - dfds21*dfds21 )                                   
      ddfds22ds22 = a1ovf*( ov4*(as11s22pa2 + as12s22pa2 + as13s22pa2)
     +            - dfds22*dfds22 )                                   
      ddfds23ds23 = a1ovf*( ov4*(as11s23pa2 + as12s23pa2 + as13s23pa2)
     +            - dfds23*dfds23 )
cc
      ddfds11ds12 = -a1ovf*dfds11*dfds12
      ddfds11ds13 = -a1ovf*dfds11*dfds13
      ddfds12ds11 = ddfds11ds12
      ddfds12ds13 = -a1ovf*dfds12*dfds13
      ddfds13ds11 = ddfds11ds13
      ddfds13ds12 = ddfds12ds13
c
      ddfds21ds22 = -a1ovf*dfds21*dfds22
      ddfds21ds23 = -a1ovf*dfds21*dfds23
      ddfds22ds21 = ddfds21ds22
      ddfds22ds23 = -a1ovf*dfds22*dfds23
      ddfds23ds21 = ddfds21ds23
      ddfds23ds22 = ddfds22ds23
cc
      ddfds11ds21 = a1ovf*(mov4*as11s21pa2 - dfds11*dfds21)
      ddfds12ds22 = a1ovf*(mov4*as12s22pa2 - dfds12*dfds22)
      ddfds13ds23 = a1ovf*(mov4*as13s23pa2 - dfds13*dfds23)
      ddfds11ds22 = a1ovf*(mov4*as11s22pa2 - dfds11*dfds22)
      ddfds11ds23 = a1ovf*(mov4*as11s23pa2 - dfds11*dfds23)
      ddfds12ds21 = a1ovf*(mov4*as12s21pa2 - dfds12*dfds21)
      ddfds12ds23 = a1ovf*(mov4*as12s23pa2 - dfds12*dfds23)
      ddfds13ds21 = a1ovf*(mov4*as13s21pa2 - dfds13*dfds21)
      ddfds13ds22 = a1ovf*(mov4*as13s22pa2 - dfds13*dfds22)
c
      ddfds21ds11 = ddfds11ds21
      ddfds21ds12 = ddfds12ds21
      ddfds21ds13 = ddfds13ds21
      ddfds22ds11 = ddfds11ds22
      ddfds22ds12 = ddfds12ds22
      ddfds22ds13 = ddfds13ds22
      ddfds23ds11 = ddfds11ds23
      ddfds23ds12 = ddfds12ds23
      ddfds23ds13 = ddfds13ds23
ccc
c
      tmp(1) = dfds11*eigvec11(1)*eigvec11(1)
     +       + dfds12*eigvec12(1)*eigvec12(1)
     +       + dfds13*eigvec13(1)*eigvec13(1)
      tmp(2) = dfds11*eigvec11(2)*eigvec11(2)
     +       + dfds12*eigvec12(2)*eigvec12(2)
     +       + dfds13*eigvec13(2)*eigvec13(2)
      tmp(3) = dfds11*eigvec11(3)*eigvec11(3)
     +       + dfds12*eigvec12(3)*eigvec12(3)
     +       + dfds13*eigvec13(3)*eigvec13(3)
	tmp(4) = dfds11*eigvec11(2)*eigvec11(3)
     +       + dfds12*eigvec12(2)*eigvec12(3)
     +       + dfds13*eigvec13(2)*eigvec13(3)
	tmp(5) = dfds11*eigvec11(1)*eigvec11(3)
     +       + dfds12*eigvec12(1)*eigvec12(3)
     +       + dfds13*eigvec13(1)*eigvec13(3)
	tmp(6) = dfds11*eigvec11(1)*eigvec11(2)
     +       + dfds12*eigvec12(1)*eigvec12(2)
     +       + dfds13*eigvec13(1)*eigvec13(2)
c
c	Compute GRAD
c
	grad(1) = a11*tmp(1) + a21*tmp(2) + a31*tmp(3)
	grad(2) = a12*tmp(1) + a22*tmp(2) + a32*tmp(3)
      grad(3) = a13*tmp(1) + a23*tmp(2) + a33*tmp(3)
      grad(4) = c44*tmp(4)
      grad(5) = c55*tmp(5)
      grad(6) = c66*tmp(6)
c	
      tmp(1) = dfds21*eigvec21(1)*eigvec21(1)
     +       + dfds22*eigvec22(1)*eigvec22(1)
     +       + dfds23*eigvec23(1)*eigvec23(1)
      tmp(2) = dfds21*eigvec21(2)*eigvec21(2)
     +       + dfds22*eigvec22(2)*eigvec22(2)
     +       + dfds23*eigvec23(2)*eigvec23(2)
	tmp(3) = dfds21*eigvec21(3)*eigvec21(3)
     +       + dfds22*eigvec22(3)*eigvec22(3)
     +       + dfds23*eigvec23(3)*eigvec23(3)
	tmp(4) = dfds21*eigvec21(2)*eigvec21(3)
     +       + dfds22*eigvec22(2)*eigvec22(3)
     +       + dfds23*eigvec23(2)*eigvec23(3)
	tmp(5) = dfds21*eigvec21(1)*eigvec21(3)
     +       + dfds22*eigvec22(1)*eigvec22(3)
     +       + dfds23*eigvec23(1)*eigvec23(3)
	tmp(6) = dfds21*eigvec21(1)*eigvec21(2)
     +       + dfds22*eigvec22(1)*eigvec22(2)
     +       + dfds23*eigvec23(1)*eigvec23(2)
c
	grad(1) = grad(1) + b11*tmp(1) + b21*tmp(2) + b31*tmp(3)
	grad(2) = grad(2) + b12*tmp(1) + b22*tmp(2) + b32*tmp(3)
      grad(3) = grad(3) + b13*tmp(1) + b23*tmp(2) + b33*tmp(3)
      grad(4) = grad(4) + d44*tmp(4) 
      grad(5) = grad(5) + d55*tmp(5) 
      grad(6) = grad(6) + d66*tmp(6) 
c
      grad(4) = 2.d0*grad(4) ! factor 2 due to Voigt
      grad(5) = 2.d0*grad(5) ! factor 2 due to Voigt
      grad(6) = 2.d0*grad(6) ! factor 2 due to Voigt
c
c     Stop here if Hessian is not required 
c     and only F and GRAD are returned
      if (mode .EQ. 1) return
c
c
      call outer2vec_Voigt(eigvec11, eigvec11, e11xe11v)
      call outer2vec_Voigt(eigvec11, eigvec12, e11xe12v)
      call outer2vec_Voigt(eigvec11, eigvec13, e11xe13v)
      call outer2vec_Voigt(eigvec12, eigvec11, e12xe11v)
      call outer2vec_Voigt(eigvec12, eigvec12, e12xe12v)
      call outer2vec_Voigt(eigvec12, eigvec13, e12xe13v)
      call outer2vec_Voigt(eigvec13, eigvec11, e13xe11v)
      call outer2vec_Voigt(eigvec13, eigvec12, e13xe12v)
      call outer2vec_Voigt(eigvec13, eigvec13, e13xe13v)
c
      call outer2vec_Voigt(eigvec21, eigvec21, e21xe21v)
      call outer2vec_Voigt(eigvec21, eigvec22, e21xe22v)
      call outer2vec_Voigt(eigvec21, eigvec23, e21xe23v)
      call outer2vec_Voigt(eigvec22, eigvec21, e22xe21v)
      call outer2vec_Voigt(eigvec22, eigvec22, e22xe22v)
      call outer2vec_Voigt(eigvec22, eigvec23, e22xe23v)
      call outer2vec_Voigt(eigvec23, eigvec21, e23xe21v)
      call outer2vec_Voigt(eigvec23, eigvec22, e23xe22v)
      call outer2vec_Voigt(eigvec23, eigvec23, e23xe23v)
c    
c
	call symouter6(e11xe11v, e11xe11v, ds11xds11)
	call symouter6(e11xe11v, e12xe12v, ds11xds12)
	call symouter6(e11xe11v, e13xe13v, ds11xds13)
      call symouter6(e12xe12v, e11xe11v, ds12xds11)
	call symouter6(e12xe12v, e12xe12v, ds12xds12)
	call symouter6(e12xe12v, e13xe13v, ds12xds13)
	call symouter6(e13xe13v, e11xe11v, ds13xds11)
	call symouter6(e13xe13v, e12xe12v, ds13xds12)
	call symouter6(e13xe13v, e13xe13v, ds13xds13)
c             					
      call symouter6(e21xe21v, e21xe21v, ds21xds21)
	call symouter6(e21xe21v, e22xe22v, ds21xds22)
	call symouter6(e21xe21v, e23xe23v, ds21xds23)
      call symouter6(e22xe22v, e21xe21v, ds22xds21)
	call symouter6(e22xe22v, e22xe22v, ds22xds22)
	call symouter6(e22xe22v, e23xe23v, ds22xds23)
	call symouter6(e23xe23v, e21xe21v, ds23xds21)
	call symouter6(e23xe23v, e22xe22v, ds23xds22)
	call symouter6(e23xe23v, e23xe23v, ds23xds23)
c           					 		  		   
	call symouter6(e11xe11v, e21xe21v, ds11xds21)
	call symouter6(e11xe11v, e22xe22v, ds11xds22)
	call symouter6(e11xe11v, e23xe23v, ds11xds23)
      call symouter6(e12xe12v, e21xe21v, ds12xds21)
	call symouter6(e12xe12v, e22xe22v, ds12xds22)
	call symouter6(e12xe12v, e23xe23v, ds12xds23)
	call symouter6(e13xe13v, e21xe21v, ds13xds21)
	call symouter6(e13xe13v, e22xe22v, ds13xds22)
	call symouter6(e13xe13v, e23xe23v, ds13xds23)
c             					 		  		   
	call symouter6(e21xe21v, e11xe11v, ds21xds11)
	call symouter6(e21xe21v, e12xe12v, ds21xds12)
	call symouter6(e21xe21v, e13xe13v, ds21xds13)
      call symouter6(e22xe22v, e11xe11v, ds22xds11)
	call symouter6(e22xe22v, e12xe12v, ds22xds12)
	call symouter6(e22xe22v, e13xe13v, ds22xds13)
	call symouter6(e23xe23v, e11xe11v, ds23xds11)
	call symouter6(e23xe23v, e12xe12v, ds23xds12)
	call symouter6(e23xe23v, e13xe13v, ds23xds13)
c
c     note, only upper triangle of 6x6 sym matrices
c     H1 and H2 is calculated. H3 and H4 are not
c     symmetric, but H3 = H4T, so the upper
c     triangle is enough to calculate
	do i=1, 6
		do j=i, 6
              H1(i,j) = ddfds11ds11*ds11xds11(i,j)
     +				+ ddfds11ds12*ds11xds12(i,j)
     +				+ ddfds11ds13*ds11xds13(i,j)
     +				+ ddfds12ds11*ds12xds11(i,j)
     +				+ ddfds12ds12*ds12xds12(i,j)
     +				+ ddfds12ds13*ds12xds13(i,j)
     +				+ ddfds13ds11*ds13xds11(i,j)
     +				+ ddfds13ds12*ds13xds12(i,j)
     +				+ ddfds13ds13*ds13xds13(i,j)
c                                                 
	 		H2(i,j) = ddfds21ds21*ds21xds21(i,j)
     +				+ ddfds21ds22*ds21xds22(i,j)
     +				+ ddfds21ds23*ds21xds23(i,j)
     +				+ ddfds22ds21*ds22xds21(i,j)
     +				+ ddfds22ds22*ds22xds22(i,j)
     +				+ ddfds22ds23*ds22xds23(i,j)
     +				+ ddfds23ds21*ds23xds21(i,j)
     +				+ ddfds23ds22*ds23xds22(i,j)
     +				+ ddfds23ds23*ds23xds23(i,j)
c                                                 
	 		H3(i,j) = ddfds11ds21*ds11xds21(i,j)
     +				+ ddfds11ds22*ds11xds22(i,j)
     +				+ ddfds11ds23*ds11xds23(i,j)
     +				+ ddfds12ds21*ds12xds21(i,j)
     +				+ ddfds12ds22*ds12xds22(i,j)
     +				+ ddfds12ds23*ds12xds23(i,j)
     +				+ ddfds13ds21*ds13xds21(i,j)
     +				+ ddfds13ds22*ds13xds22(i,j)
     +				+ ddfds13ds23*ds13xds23(i,j)
c                                                 
	 		H4(i,j) = ddfds21ds11*ds21xds11(i,j)
     +				+ ddfds21ds12*ds21xds12(i,j)
     +				+ ddfds21ds13*ds21xds13(i,j)
     +				+ ddfds22ds11*ds22xds11(i,j)
     +				+ ddfds22ds12*ds22xds12(i,j)
     +				+ ddfds22ds13*ds22xds13(i,j)
     +				+ ddfds23ds11*ds23xds11(i,j)
     +				+ ddfds23ds12*ds23xds12(i,j)
     +				+ ddfds23ds13*ds23xds13(i,j)
c
          end do
      end do   
c
c     Eq (32) single prime
	call symouter6(e11xe12v, e11xe12v, Etmp1212)
	call symouter6(e11xe12v, e12xe11v, Etmp1221)
	call symouter6(e12xe11v, e11xe12v, Etmp2112)
	call symouter6(e12xe11v, e12xe11v, Etmp2121)
c
      do i=1, 6
          do j=i, 6
	        E1212p(i,j) = Etmp1212(i,j) + Etmp1221(i,j)
     +                    + Etmp2112(i,j) + Etmp2121(i,j)
          end do
      end do
c
      call symouter6(e12xe13v, e12xe13v, Etmp2323)
	call symouter6(e12xe13v, e13xe12v, Etmp2332)
	call symouter6(e13xe12v, e12xe13v, Etmp3223)
	call symouter6(e13xe12v, e13xe12v, Etmp3232)
c
      do i=1, 6
          do j=i, 6
      	    E2323p(i,j) = Etmp2323(i,j) + Etmp2332(i,j)
     +                    + Etmp3223(i,j) + Etmp3232(i,j)
          end do
      end do
c
	call symouter6(e13xe11v, e13xe11v, Etmp3131)
	call symouter6(e13xe11v, e11xe13v, Etmp3113)
	call symouter6(e11xe13v, e13xe11v, Etmp1331)
	call symouter6(e11xe13v, e11xe13v, Etmp1313)
c
      do i=1, 6
          do j=i, 6
	        E3131p(i,j) = Etmp3131(i,j) + Etmp3113(i,j)
     +                    + Etmp1331(i,j) + Etmp1313(i,j)
          end do
      end do
cc
c     Eq (32) double prime
	call symouter6(e21xe22v, e21xe22v, Etmp1212)
	call symouter6(e21xe22v, e22xe21v, Etmp1221)
	call symouter6(e22xe21v, e21xe22v, Etmp2112)
	call symouter6(e22xe21v, e22xe21v, Etmp2121)
c
      do i=1, 6
          do j=i, 6
	        E1212pp(i,j) = Etmp1212(i,j) + Etmp1221(i,j)
     +                     + Etmp2112(i,j) + Etmp2121(i,j)
          end do
      end do
c
      call symouter6(e22xe23v, e22xe23v, Etmp2323)
	call symouter6(e22xe23v, e23xe22v, Etmp2332)
	call symouter6(e23xe22v, e22xe23v, Etmp3223)
	call symouter6(e23xe22v, e23xe22v, Etmp3232)
c
      do i=1, 6
          do j=i, 6
	        E2323pp(i,j) = Etmp2323(i,j) + Etmp2332(i,j)
     +                     + Etmp3223(i,j) + Etmp3232(i,j)
          end do
      end do
c
	call symouter6(e23xe21v, e23xe21v, Etmp3131)
	call symouter6(e23xe21v, e21xe23v, Etmp3113)
	call symouter6(e21xe23v, e23xe21v, Etmp1331)
	call symouter6(e21xe23v, e21xe23v, Etmp1313)
c
      do i=1, 6
          do j=i, 6
	        E3131pp(i,j) = Etmp3131(i,j) + Etmp3113(i,j)
     +                     + Etmp1331(i,j) + Etmp1313(i,j)
          end do
      end do
c
c     single primed
	if (abs(s11-s12) .LT. tol) then
		coefE1212p = ddfds11ds11 - ddfds11ds12
	else
		coefE1212p = (dfds11 - dfds12)/(s11-s12)
	end if
c
	if (abs(s12-s13) .LT. tol) then
		coefE2323p = ddfds12ds12 - ddfds12ds13
	else
		coefE2323p = (dfds12 - dfds13)/(s12-s13)
	end if
c
	if (abs(s13-s11) .LT. tol) then
		coefE3131p = ddfds13ds13 - ddfds13ds11
	else
		coefE3131p = (dfds13 - dfds11)/(s13-s11)
      end if
c
cc    double primed
      if (abs(s21-s22) .LT. tol) then
		coefE1212pp = ddfds21ds21 - ddfds21ds22
	else
		coefE1212pp = (dfds21 - dfds22)/(s21-s22)
	end if
c
	if (abs(s22-s23) .LT. tol) then
		coefE2323pp = ddfds22ds22 - ddfds22ds23
	else
		coefE2323pp = (dfds22 - dfds23)/(s22-s23)
	end if
c
	if (abs(s23-s21) .LT. tol) then
		coefE3131pp = ddfds23ds23 - ddfds23ds21
	else
		coefE3131pp = (dfds23 - dfds21)/(s23-s21)
      end if
c	
c     note, only upper triangle of 6x6 sym matrices
c     H5 and H6 is calculated
	do i=1, 6
		do j=i, 6
			H5(i,j) = (coefE1212p*E1212p(i,j)
     +                +  coefE2323p*E2323p(i,j)         
     +                +  coefE3131p*E3131p(i,j))*ov2/f
c
			H6(i,j) = (coefE1212pp*E1212pp(i,j)
     +                +  coefE2323pp*E2323pp(i,j)       
     +                +  coefE3131pp*E3131pp(i,j))*ov2/f
c
          end do
      end do
c
c     summing up H1 with H5, and H2 with H6, before tarnsforming them
c     together by L and M, respectively
      do i=1, 6
          do j=i, 6
              H7(i,j) = H1(i,j) + H5(i,j)
              H8(i,j) = H2(i,j) + H6(i,j)
          end do
      end do
c
c     copmuting upper part of sym matrix (LT*H3*M + MT*H4*L) (Eq. 25)
c     note that H3 = H4T, then 
c     dsymLTxH3xM = LT*H3*M + (LT*H3*M)T = 2*sym(LT*H3*M)
      tmp1 = (H3(1,1)*b11 + H3(1,2)*b21 + H3(1,3)*b31)
      tmp2 = (H4(1,2)*b11 + H3(2,2)*b21 + H3(2,3)*b31)
      tmp3 = (H4(1,3)*b11 + H4(2,3)*b21 + H3(3,3)*b31)
      tmp4 = (H3(1,1)*a11 + H4(1,2)*a21 + H4(1,3)*a31)
      tmp5 = (H3(1,2)*a11 + H3(2,2)*a21 + H4(2,3)*a31)
      tmp6 = (H3(1,3)*a11 + H3(2,3)*a21 + H3(3,3)*a31)
      dsymLTxH3xM(1,1) = (a11*tmp1 + a21*tmp2 + a31*tmp3)*2.d0
      dsymLTxH3xM(1,2) = a12*tmp1 + a22*tmp2 + a32*tmp3 
     +                 + b12*tmp4 + b22*tmp5 + b32*tmp6
      dsymLTxH3xM(1,3) = a13*tmp1 + a23*tmp2 + a33*tmp3 
     +                 + b13*tmp4 + b23*tmp5 + b33*tmp6
      dsymLTxH3xM(1,4) = c44*(H4(1,4)*b11 + H4(2,4)*b21 + H4(3,4)*b31) 
     +                 + d44*(H3(1,4)*a11 + H3(2,4)*a21 + H3(3,4)*a31)
      dsymLTxH3xM(1,5) = c55*(H4(1,5)*b11 + H4(2,5)*b21 + H4(3,5)*b31) 
     +                 + d55*(H3(1,5)*a11 + H3(2,5)*a21 + H3(3,5)*a31)
      dsymLTxH3xM(1,6) = c66*(H4(1,6)*b11 + H4(2,6)*b21 + H4(3,6)*b31) 
     +                 + d66*(H3(1,6)*a11 + H3(2,6)*a21 + H3(3,6)*a31)
      tmp1 = (H3(1,1)*b12 + H3(1,2)*b22 + H3(1,3)*b32)
      tmp2 = (H4(1,2)*b12 + H3(2,2)*b22 + H3(2,3)*b32)
      tmp3 = (H4(1,3)*b12 + H4(2,3)*b22 + H3(3,3)*b32)
      dsymLTxH3xM(2,2) = (a12*tmp1 + a22*tmp2 + a32*tmp3)*2.d0
      dsymLTxH3xM(2,3) = a13*tmp1 + a23*tmp2 + a33*tmp3 
     +                 + b13*(H3(1,1)*a12 + H4(1,2)*a22 + H4(1,3)*a32) 
     +                 + b23*(H3(1,2)*a12 + H3(2,2)*a22 + H4(2,3)*a32) 
     +                 + b33*(H3(1,3)*a12 + H3(2,3)*a22 + H3(3,3)*a32)
      dsymLTxH3xM(2,4) = c44*(H4(1,4)*b12 + H4(2,4)*b22 + H4(3,4)*b32)
     +                 + d44*(H3(1,4)*a12 + H3(2,4)*a22 + H3(3,4)*a32)
      dsymLTxH3xM(2,5) = c55*(H4(1,5)*b12 + H4(2,5)*b22 + H4(3,5)*b32)
     +                 + d55*(H3(1,5)*a12 + H3(2,5)*a22 + H3(3,5)*a32)
      dsymLTxH3xM(2,6) = c66*(H4(1,6)*b12 + H4(2,6)*b22 + H4(3,6)*b32)
     +                 + d66*(H3(1,6)*a12 + H3(2,6)*a22 + H3(3,6)*a32)
      dsymLTxH3xM(3,3) = (a13*(H3(1,1)*b13 + H3(1,2)*b23 + H3(1,3)*b33)
     +                 +  a23*(H4(1,2)*b13 + H3(2,2)*b23 + H3(2,3)*b33)
     +                 +  a33*(H4(1,3)*b13 + H4(2,3)*b23 + H3(3,3)*b33))
     +                 * 2.d0
      dsymLTxH3xM(3,4) = c44*(H4(1,4)*b13 + H4(2,4)*b23 + H4(3,4)*b33)
     +                 + d44*(H3(1,4)*a13 + H3(2,4)*a23 + H3(3,4)*a33)
      dsymLTxH3xM(3,5) = c55*(H4(1,5)*b13 + H4(2,5)*b23 + H4(3,5)*b33)
     +                 + d55*(H3(1,5)*a13 + H3(2,5)*a23 + H3(3,5)*a33)
      dsymLTxH3xM(3,6) = c66*(H4(1,6)*b13 + H4(2,6)*b23 + H4(3,6)*b33)
     +                 + d66*(H3(1,6)*a13 + H3(2,6)*a23 + H3(3,6)*a33)
      dsymLTxH3xM(4,4) = 2.d0*H3(4,4)*c44*d44
      dsymLTxH3xM(4,5) = H3(4,5)*c44*d55 + H4(4,5)*c55*d44
      dsymLTxH3xM(4,6) = H3(4,6)*c44*d66 + H4(4,6)*c66*d44
      dsymLTxH3xM(5,5) = 2.d0*H3(5,5)*c55*d55
      dsymLTxH3xM(5,6) = H3(5,6)*c55*d66 + H4(5,6)*c66*d55
      dsymLTxH3xM(6,6) = 2.d0*H3(6,6)*c66*d66
      
c     copmuting upper part of sym matrix LT*H7*L = LT*(H1+H5)*L (Eq. 25)
      tmp1 = (H7(1,1)*a11 + H7(1,2)*a21 + H7(1,3)*a31)
      tmp2 = (H7(1,2)*a11 + H7(2,2)*a21 + H7(2,3)*a31)
      tmp3 = (H7(1,3)*a11 + H7(2,3)*a21 + H7(3,3)*a31)
      LTxH7xL(1,1) = a11*tmp1 + a21*tmp2 + a31*tmp3
      LTxH7xL(1,2) = a12*tmp1 + a22*tmp2 + a32*tmp3
      LTxH7xL(1,3) = a13*tmp1 + a23*tmp2 + a33*tmp3
      LTxH7xL(1,4) = c44*(H7(1,4)*a11 + H7(2,4)*a21 + H7(3,4)*a31)
      LTxH7xL(1,5) = c55*(H7(1,5)*a11 + H7(2,5)*a21 + H7(3,5)*a31)
      LTxH7xL(1,6) = c66*(H7(1,6)*a11 + H7(2,6)*a21 + H7(3,6)*a31)
      tmp1 = (H7(1,1)*a12 + H7(1,2)*a22 + H7(1,3)*a32)
      tmp2 = (H7(1,2)*a12 + H7(2,2)*a22 + H7(2,3)*a32)
      tmp3 = (H7(1,3)*a12 + H7(2,3)*a22 + H7(3,3)*a32)  
      LTxH7xL(2,2) = a12*tmp1 + a22*tmp2 + a32*tmp3
      LTxH7xL(2,3) = a13*tmp1 + a23*tmp2 + a33*tmp3
      LTxH7xL(2,4) = c44*(H7(1,4)*a12 + H7(2,4)*a22 + H7(3,4)*a32)
      LTxH7xL(2,5) = c55*(H7(1,5)*a12 + H7(2,5)*a22 + H7(3,5)*a32)
      LTxH7xL(2,6) = c66*(H7(1,6)*a12 + H7(2,6)*a22 + H7(3,6)*a32)
      LTxH7xL(3,3) = a13*(H7(1,1)*a13 + H7(1,2)*a23 + H7(1,3)*a33) 
     +             + a23*(H7(1,2)*a13 + H7(2,2)*a23 + H7(2,3)*a33)
     +             + a33*(H7(1,3)*a13 + H7(2,3)*a23 + H7(3,3)*a33)
      LTxH7xL(3,4) = c44*(H7(1,4)*a13 + H7(2,4)*a23 + H7(3,4)*a33)
      LTxH7xL(3,5) = c55*(H7(1,5)*a13 + H7(2,5)*a23 + H7(3,5)*a33)
      LTxH7xL(3,6) = c66*(H7(1,6)*a13 + H7(2,6)*a23 + H7(3,6)*a33)
      LTxH7xL(4,4) = H7(4,4)*c44**2
      LTxH7xL(4,5) = H7(4,5)*c44*c55
      LTxH7xL(4,6) = H7(4,6)*c44*c66
      LTxH7xL(5,5) = H7(5,5)*c55*c55
      LTxH7xL(5,6) = H7(5,6)*c55*c66
      LTxH7xL(6,6) = H7(6,6)*c66**2
c
c     copmuting upper part of sym matrix MT*H8*M = MT*(H2+H6)*M (Eq. 25)
      tmp1 = (H8(1,1)*b11 + H8(1,2)*b21 + H8(1,3)*b31)
      tmp2 = (H8(1,2)*b11 + H8(2,2)*b21 + H8(2,3)*b31)
      tmp3 = (H8(1,3)*b11 + H8(2,3)*b21 + H8(3,3)*b31)
      MTxH8xM(1,1) = b11*tmp1 + b21*tmp2 + b31*tmp3
      MTxH8xM(1,2) = b12*tmp1 + b22*tmp2 + b32*tmp3
      MTxH8xM(1,3) = b13*tmp1 + b23*tmp2 + b33*tmp3
      MTxH8xM(1,4) = d44*(H8(1,4)*b11 + H8(2,4)*b21 + H8(3,4)*b31)
      MTxH8xM(1,5) = d55*(H8(1,5)*b11 + H8(2,5)*b21 + H8(3,5)*b31)
      MTxH8xM(1,6) = d66*(H8(1,6)*b11 + H8(2,6)*b21 + H8(3,6)*b31)
      tmp1 = (H8(1,1)*b12 + H8(1,2)*b22 + H8(1,3)*b32)
      tmp2 = (H8(1,2)*b12 + H8(2,2)*b22 + H8(2,3)*b32)
      tmp3 = (H8(1,3)*b12 + H8(2,3)*b22 + H8(3,3)*b32)  
      MTxH8xM(2,2) = b12*tmp1 + b22*tmp2 + b32*tmp3
      MTxH8xM(2,3) = b13*tmp1 + b23*tmp2 + b33*tmp3
      MTxH8xM(2,4) = d44*(H8(1,4)*b12 + H8(2,4)*b22 + H8(3,4)*b32)
      MTxH8xM(2,5) = d55*(H8(1,5)*b12 + H8(2,5)*b22 + H8(3,5)*b32)
      MTxH8xM(2,6) = d66*(H8(1,6)*b12 + H8(2,6)*b22 + H8(3,6)*b32)
      MTxH8xM(3,3) = b13*(H8(1,1)*b13 + H8(1,2)*b23 + H8(1,3)*b33) 
     +             + b23*(H8(1,2)*b13 + H8(2,2)*b23 + H8(2,3)*b33)
     +             + b33*(H8(1,3)*b13 + H8(2,3)*b23 + H8(3,3)*b33)
      MTxH8xM(3,4) = d44*(H8(1,4)*b13 + H8(2,4)*b23 + H8(3,4)*b33)
      MTxH8xM(3,5) = d55*(H8(1,5)*b13 + H8(2,5)*b23 + H8(3,5)*b33)
      MTxH8xM(3,6) = d66*(H8(1,6)*b13 + H8(2,6)*b23 + H8(3,6)*b33)
      MTxH8xM(4,4) = H8(4,4)*d44**2
      MTxH8xM(4,5) = H8(4,5)*d44*d55
      MTxH8xM(4,6) = H8(4,6)*d44*d66
      MTxH8xM(5,5) = H8(5,5)*d55*d55
      MTxH8xM(5,6) = H8(5,6)*d55*d66
      MTxH8xM(6,6) = H8(6,6)*d66**2
c
c     suming up the final hessian
      do i=1, 6
          do j=i, 6
              hessian(i,j) = LTxH7xL(i,j) 
     +                     + MTxH8xM(i,j) 
     +                     + dsymLTxH3xM(i,j)              
              ! factor 2 due to Voigt
              if (i .GT. 3) hessian(i,j) = hessian(i,j)*2.d0
              if (j .GT. 3) hessian(i,j) = hessian(i,j)*2.d0
              if (i .NE. j) then
                  hessian(j,i) = hessian(i,j)
              end if
          end do
      end do
c
      return
      end subroutine Yld2004_18p_grad
c
c
      subroutine eig2(Amat, mode, e1, e2, e3, evec1, evec2, evec3)
c**********************************************************************
c      Calculates eigenvalues and eigenvectors of a 3x3 symmetrix 
c      matrix. Based on paper by
c      Scherzinger, W. M. and Dohrmann, C. R., A robust algorithm for 
c      finding the eigenvalues and eigenvectors of 3 x 3 symmetric 
c      matrices, Computer Methods in Applied Mechanics and Engineering,
c      2008, DOI: 10.1016/j.cma.2008.03.031
c
c**********************************************************************
c (IN)    AMAT is REAL*8 array, dimension (3,3)
c         It is a symmetric 3x3 matrix
c (IN)    MODE is integer
c         MODE = 0 returns only eigenvalues, eigenvectors will not be 
c         calculated
c         MODE = else returns both eigenvalues and eigenvectors
c (OUT)   E1, E2, E3 are REAL*8 
c         three eigenvalues
c (OUT)   EVEC1, EVEC2, EVEC3 are REAL*8 array, dimension (3)
c         three eigenvectors
c         
c**********************************************************************
      implicit none  
c     ! in/out
      integer     mode
      real*8      e1, e2, e3, evec1(3), evec2(3), evec3(3)
c     ! internal
      integer     i, j, flag, imax, k, m
      real*8      J2, J3, PI, r(3), nr, sdott,                 
     +            t(3,2), As1(3), As2(3), AR11, AR22, AR23, AR32,   
     +            AR33, s1(3), s2(3), a1, a2, a3, cos3a, ns, w1(3),
     +            Amat(3,3), u1(3), u2(3), p
c
      PI = 4.d0*atan(1.d0)
c
      p = (Amat(1,1)+Amat(2,2)+Amat(3,3))/3.d0
      do i=1, 3
          Amat(i,i) = Amat(i,i) - p
      end do
c
      J2 = Amat(1,1)**2+Amat(2,2)**2+Amat(1,2)**2+Amat(1,3)**2
     +   + Amat(2,3)**2+Amat(1,1)*Amat(2,2)
c
      J3 = Amat(1,1)*(Amat(2,2)*Amat(3,3)-Amat(3,2)**2)
     +   + Amat(1,3)*(2.d0*Amat(2,1)*Amat(3,2) - Amat(2,2)*Amat(1,3))
     +   - Amat(3,3)*Amat(2,1)**2
c     
      if (J2 .LT. 1.d-30) then
          e1 = 0.d0
          e2 = 0.d0
          e3 = 0.d0
          evec1(1) = 1.d0
          evec1(2) = 0.d0
          evec1(3) = 0.d0
          evec2(1) = 0.d0
          evec2(2) = 1.d0
          evec2(3) = 0.d0
          evec3(1) = 0.d0
          evec3(2) = 0.d0
          evec3(3) = 1.d0
      goto 100
      end if
      cos3a = 0.5d0*J3*(3.d0/J2)**(3.d0/2.d0)
c     to make cos3a within [-1, 1] interval
      cos3a = max(-1.d0, min(1.d0, cos3a))
c
      a1 = acos(cos3a)/3.d0
      a3 = a1 + 2.d0/3.d0*PI
      a2 = a1 + 4.d0/3.d0*PI
c
      if (a1 .LT. PI/6.d0) then
          e1 = 2.d0*sqrt(J2/3.d0)*cos(a1)
      else
          e1 = 2.d0*sqrt(J2/3.d0)*cos(a3)
      end if
c      
      do i=1, 3
          Amat(i,i) = Amat(i,i) - e1
      end do
c
c     Find the largest column of Amat and store as s1
      ns = 0.d0
      do j=1, 3
          nr = Amat(1,j)**2+Amat(2,j)**2+Amat(3,j)**2
          if (nr .GT. ns) then
              ns = nr
              imax = j
              do i=1, 3
                  s1(i) = Amat(i,j)
              end do
          end if
      end do
c
      do i=1, 3
          s1(i) = s1(i)/sqrt(ns)
      end do
c
      m = 1
      do j=1, 3
          if (j .NE. imax) then
              sdott = s1(1)*Amat(1,j)+s1(2)*Amat(2,j)+s1(3)*Amat(3,j)
              do i=1, 3
                  t(i,m) = Amat(i,j) - sdott*s1(i)
              end do
              m = m+1
          end if
      end do
c
c     Find the largest t column and store as s2
      ns = 0.d0
      do j=1, 2
          nr = t(1,j)**2+t(2,j)**2+t(3,j)**2
          if (nr .GT. ns) then
              ns = nr
              do i=1, 3
                  s2(i) = t(i,j)
              end do
          end if
      end do
c
      do i=1, 3
          s2(i) = s2(i)/sqrt(ns)
      end do
c
c     First eigenvector v1
      evec1(1) = s1(2)*s2(3)-s2(2)*s1(3)
      evec1(2) = s1(3)*s2(1)-s2(3)*s1(1)
      evec1(3) = s1(1)*s2(2)-s2(1)*s1(2)
c
c     Build reduced form of A' matrix (Eq. 22)
      do i=1, 3
          Amat(i,i) = Amat(i,i) + e1
      end do
c
      AR11 = e1
      do i=1, 3
          As1(i) = Amat(i,1)*s1(1) + Amat(i,2)*s1(2) + Amat(i,3)*s1(3)
          As2(i) = Amat(i,1)*s2(1) + Amat(i,2)*s2(2) + Amat(i,3)*s2(3)
      end do
c
      AR22 = s1(1)*As1(1) + s1(2)*As1(2) + s1(3)*As1(3)
      AR23 = s1(1)*As2(1) + s1(2)*As2(2) + s1(3)*As2(3)
      AR32 = s2(1)*As1(1) + s2(2)*As1(2) + s2(3)*As1(3)
      AR33 = s2(1)*As2(1) + s2(2)*As2(2) + s2(3)*As2(3)
c
c     Find the remaining eigenvalues e2, e3 by the Wilkinsons shift
      e2 = 0.5d0*(AR22+AR33) - 0.5d0*sign(1.d0, AR22-AR33)
     +   * sqrt((AR22-AR33)**2 + 4.d0*AR23*AR32)
      e3 = AR22 + AR33 - e2
c
c     returns here if only eigenvalues are required
      if (mode .EQ. 0) goto 100
c
c     Find eigenvectors evec2 and evec3
      do i=1, 3
          Amat(i,i) = Amat(i,i) - e2
      end do
c
      do i=1, 3
          u1(i) = Amat(i,1)*s1(1) + Amat(i,2)*s1(2) + Amat(i,3)*s1(3)
          u2(i) = Amat(i,1)*s2(1) + Amat(i,2)*s2(2) + Amat(i,3)*s2(3)
      end do
c
      nr = u1(1)**2 + u1(2)**2 + u1(3)**2
      ns = u2(1)**2 + u2(2)**2 + u2(3)**2
c     if s1 and s2 are already second and third eigenvectors, then
c     both u1 and u2 and their norms equal zero (ELSE branch)
      if ((nr .GT. 1.d-30) .or. (ns .GT. 1.d-30)) then
          if (nr .GT. ns) then
              do i=1, 3
                  w1(i) = u1(i)/sqrt(nr)
              end do
          else
              do i=1, 3
                  w1(i) = u2(i)/sqrt(ns)
              end do
          end if
          evec2(1) = w1(2)*evec1(3)-evec1(2)*w1(3)
          evec2(2) = w1(3)*evec1(1)-evec1(3)*w1(1)
          evec2(3) = w1(1)*evec1(2)-evec1(1)*w1(2)
c
          evec3(1) = evec1(2)*evec2(3)-evec2(2)*evec1(3)
          evec3(2) = evec1(3)*evec2(1)-evec2(3)*evec1(1)
          evec3(3) = evec1(1)*evec2(2)-evec2(1)*evec1(2)
      else
          do i=1, 3
              evec2(i) = s1(i)
              evec3(i) = s2(i)
          end do
      end if
c
c     adding pressure to get final eigenvalues
100   e1 = e1 + p
      e2 = e2 + p
      e3 = e3 + p
c      
c     calculate the original Amat
      do i=1, 3
          Amat(i,i) = Amat(i,i) + e2
      end do
      return
      end subroutine eig2
c
c
      subroutine outer2vec(e1, e2, v)
c
c     Uses the natural notation
c
	implicit none
	integer :: i, j
	real(8) :: e1(3), e2(3), v(6), tmp(6)
      real(8), parameter :: ovsqrt3 = 1.d0/sqrt(3.d0), 
     +                      ovsqrt2 = 1.d0/sqrt(2.d0),
     +                      ovsqrt6 = ovsqrt3*ovsqrt2, 
     +                      sqrt2   = sqrt(2.d0)
c	
	tmp(1) = e1(1)*e2(1)
	tmp(2) = e1(2)*e2(2)
	tmp(3) = e1(3)*e2(3)
	tmp(4) = e1(2)*e2(3)
	tmp(5) = e1(1)*e2(3)
	tmp(6) = e1(1)*e2(2)
c	
	v(1) = ovsqrt3*(tmp(1)+tmp(2)+tmp(3))
	v(2) = ovsqrt6*(2.d0*tmp(3)-tmp(1)-tmp(2))
	v(3) = ovsqrt2*(tmp(2)-tmp(1))
	v(4) = sqrt2*tmp(4)
	v(5) = sqrt2*tmp(5)
	v(6) = sqrt2*tmp(6)
c
      return
      end subroutine outer2vec
c
c
      subroutine symouter6(e1, e2, M)
c
	implicit none
	integer :: i, j
	real(8) :: e1(6), e2(6), M(6,6)
c		
      do i=1, 6
		do j=i, 6
			M(i,j) = e1(i)*e2(j)
		end do
      end do
c
      return
      end subroutine symouter6
c
c
      real*8 function vdot(x, y, n)
c          
      implicit none
      integer :: n, i
      real(8) :: x(n), y(n)
c          
      vdot = 0.d0
      do i = 1, n
          vdot = vdot + x(i)*y(i)
      end do
c          
      return
      end function vdot
c
c
      subroutine outer2vec_Voigt(e1, e2, v)
c
c     Uses Voigt notation
c
	implicit none
	integer     i, j
	real*8      e1(3), e2(3), v(6)
c	
	v(1) = e1(1)*e2(1)
	v(2) = e1(2)*e2(2)
	v(3) = e1(3)*e2(3)
	v(4) = e1(2)*e2(3)
	v(5) = e1(1)*e2(3)
	v(6) = e1(1)*e2(2)
c
      return
      end subroutine outer2vec_Voigt
c
c
      SUBROUTINE MMULT(A,B,N,AB)
C
      IMPLICIT NONE
      REAL*8 A(N,N),B(N,N),AB(N,N),P
      INTEGER I,J,K,N
C
      DO I=1,N
         DO J=1,N
            P=0.D0
            DO K=1,N
               P=P+A(I,K)*B(K,J)
            ENDDO
            AB(I,J)=P
         ENDDO
      ENDDO
C
      RETURN
      END SUBROUTINE MMULT