c		This file contains subroutines hand-picked from the file
c		labeled 'FROMCEM09', which is used in LAGSM calculations.
c		This file contains only a few subroutines which are no
c		longer in the CEM program.
c
c		SUBROUTINES needed from the 'FROMCEM09' file include:
c				1) GAMAGU2
c				2) GAMBETN
c				3) GEMBETP
c				4) TKINM3
c				5) VHELP
c
c		The majority of this code is written in F66/F77
c
c		Code gathered by CMJ, 07/2016
c
c

      function gamagu2 (j)

c ======================================================================
c
c   Finds the emission width for a preequilbrium particle of type "j"
c   by integrating Eq. (14) of Gudima, et al. over excitation energy.
c
c   This version incorporates the entire calculation into the subroutine
c   to enhance execution efficiency relative to the calls to GU8 and
c   WBE.
c
c   Called by: PRECOF
c
c   CEM95 written by S. G. Mashnik
c
c    Edited by A. J. Sierk   LANL  T-2,  February, 1996.
c    Edited by AJS, August, 1997.
c    Modified by AJS, January, 1999.
c    Converted from GAMAGU, GU6, GU6LA, WBE by AJS, January, 1999.
c    Small efficiency improvements, AJS, March, 1999.
c   "Last" change: 12-AUG-2003 by NVMokhov
c    Edited by A. J. Sierk, LANL T-16, October, 2003.
c    Modified by NVM, July, 2005.
c
c ======================================================================
 
      implicit real*8 (a-h, o-z), integer (i-n)

c ======================================================================

      dimension abla(6), absc(6), w(6), wtla(6)
      dimension emmp1(4), emmp2(4), mmpr1(4), mmpr2(4)
      dimension ddenom(4)

      common /bl1005/  aj(66), ajthr(66)
      common /bl1011/ vj(6)
      common /bl1015/ rj(7), uej(6), pevapj(6)
      common /bl1018/ afjthr(7), athrd
      common /blac/   ac
      common /blalj/  alj(6), ami(6)
      common /blbj/   bj(7)
      common /blexn/  exn
      common /blr0/   r0
      common /br0j/   r0j(6), bn
      common /gambet/ gb(6)

c  Gauss-Legendre Weights
      data w   / 0.171324492379170d0, 0.360761573048139d0, 
     &           0.467913934572691d0, 0.467913934572691d0,
     &           0.360761573048139d0, 0.171324492379170d0/
c  Gauss-Legendre abscissas
      data absc/-0.932469514203152d0,-0.661209386466265d0,
     &          -0.238619186083197d0, 0.238619186083197d0, 
     &           0.661209386466265d0, 0.932469514203152d0/
c  Gauss-Laguerre Weights [Real weights * exp(abla(k))]
      data wtla/ 0.573535507423d0, 1.36925259071d0, 2.26068459338d0, 
     &           3.35052458236d0,  4.88682680021d0, 7.84901594560d0/
c  Gauss-Laguerre abscissas
      data abla/ 0.222846604179d0, 1.188932101673d0, 2.992736326059d0,
     &           5.775143569105d0, 9.837467418383d0,15.982873980602d0/

      data zro, hlf, one /0.d0, 0.5d0, 1.d0/
c   Define auxiliary arrays to speed up calculation by not calculating
c   them each time routine is called!
c                     d       t      He3     He4
      data ddenom /3.75d0, 8.75d0, 8.75d0, 15.75d0/
      data emmp1  / 1.5d0,  2.5d0,  2.5d0,  3.5d0/
      data emmp2  / 2.5d0,  3.5d0,  3.5d0,  4.5d0/
      data mmpr1  /   1,      2,      2,      3/
      data mmpr2  /   2,      3,      3,      4/

      save ifirst,cst0,cst1
      data ifirst /1/
c ======================================================================

      if (ifirst.eq.1) then
        cst0 = 0.000234d0*r0*r0
c   0.104 is an approximation to 3*c/(2*sqrt(2*940)) = 0.103713. (AJS)
        cst1 = 0.104d0/r0
        ifirst = 0
      endif
      dq = alj(j)*r0j(j)/afjthr(j)
      a1 = vj(j)
      b1 = uej(j) - bj(j)
      b2 = b1 - a1
      nn = nint(exn)
      nm1 = nn - 1
      if (j.le.2) then
c   Use analytic forms for protons & neutrons:
        c1 = b2/uej(j)
        wbe0 = cst0*dq/ac
        dm = exn*(exn - one)
        if (j.eq.2) then
          denom = b2
        elseif (j.eq.1) then
          denom = b1 + exn*bn
        endif
        y = denom*wbe0*c1**nm1/dm
      elseif (j.gt.2) then
        c2 = b1/uej(j)
        n = nint(exn - aj(j) - one)
        f3 = sqrt(b1*aj(j))
        cst2 = cst1*gb(j)*dq/f3
        spnrm = zro
        if (n.eq.0) then
c   Use analytic form for integral when n = 0:
          q2 = bj(j)/b1
          q3 = vj(j)/b1
          q1 = q2 + q3
          sq1 = sqrt(q1)
          sq2 = sqrt(c2)
          spnrm = b1*(sq1*q1**mmpr2(j-2) + (emmp1(j-2) - q2 -
     &                 emmp2(j-2)*q3)/(sq2*c2**mmpr1(j-2)))/ddenom(j-2)
        elseif (nn.le.15) then
c   Use quadrature for composite particles with n > 0:
c   Gauss-Legendre:
          d1 = hlf*b2
          d2 = hlf*(b1 + a1)
            do 20 k = 1,6
            e = d1*absc(k) + d2
            x = e/b1
            fct = max((e + bj(j))/b1, zro)
            sqr = sqrt(fct)
            f2 = (e - a1)/b1
            f4 = f2*(one - x)**n
            wbe = f4*sqr
            if (j.gt.3) then
              wbe = fct*wbe
              if (j.eq.6) wbe = fct*wbe
            endif
            spnrm = spnrm + w(k)*wbe
   20       continue
          spnrm = d1*spnrm
        elseif (nn.gt.15) then
c   Gauss-Laguerre:
          gam1 = 1.5d0*dble(n)
          gam2 = (abla(6) + one)*uej(j)/b2
          gam = max(gam1, gam2)
          d2 = uej(j)/gam
            do 40 k = 1,6
            ela = a1 + d2*abla(k)
            x = ela/b1
            fct = max((ela + bj(j))/b1, zro)
            sqr = sqrt(fct)
            f2 = (ela - a1)/b1
            f4 = f2*(one - x)**n
            wbe = f4*sqr
            if (j.gt.3) then
              wbe = fct*wbe
              if (j.eq.6) wbe = fct*wbe
            endif
            spnrm = spnrm + wbe*wtla(k)
   40       continue
          spnrm = d2*spnrm
        endif
        y = cst2*spnrm*c2**nm1
      endif
      gamagu2 = y
      return
      end

c ======================================================================
			
      subroutine gambetn (gbm)

c ======================================================================
c
c    Linear interpolation of empirical gamma_beta multipliers = 
c    f(A0, t0) for neutron-induced preequilibrium reactions.
c
c    Written by K. G. Gudima using the Baznat's fits, Jul.-Oct. 2004.
c    Edited by A. J. Sierk, LANL T-16,  January, 2005.
c
c ======================================================================

      implicit real*8 (a-h, o-z), integer (i-n)

c ======================================================================

      common /targ0/ atar0, ztar0, ener0

      dimension gbm(4)

      dimension tfit(11), afit(5), xq(11), yq(4,11), y(4), y1(4), y2(4),  
     &      al27d(11) , si28d(11) , fe56d(11) , bi209d(11) , U238d(11) ,
     &      al27t(11) , si28t(11) , fe56t(11) , bi209t(11) , U238t(11) ,
     &      al27h3(11), si28h3(11), fe56h3(11), bi209h3(11), U238h3(11),
     &      al27h4(11), si28h4(11), fe56h4(11), bi209h4(11), U238h4(11)
      data
     &     afit   /27.d0, 28.d0, 56.d0, 209.d0, 238.d0/,
     &     tfit   /0.028d0, 0.031d0, 0.034d0, 0.037d0, 0.041d0, 0.045d0,
     &             0.049d0, 0.053d0, 0.063d0, 0.096d0, 0.300d0/
      data
     &     al27d  /4*1.00d0, 2.00d0, 3.00d0, 4*4.00d0, 1.00d0/,
     &     al27t  /4*1.00d0, 1.00d0, 2.00d0, 4*4.00d0, 1.00d0/,
     &     al27h3 /4*1.00d0, 1.00d0, 2.00d0, 4*4.00d0, 1.00d0/,
     &     al27h4 /9*100.0d0, 30.d0, 1.00d0/
      data 
     &     si28d  /4*1.00d0, 2.00d0, 3.00d0, 4*4.00d0, 1.00d0/,
     &     si28t  /4*1.00d0, 2.00d0, 4.00d0, 3*8.00d0, 4.0d0, 1.00d0/,
     &     si28h3 /4*1.00d0, 2.00d0, 4.00d0, 3*8.00d0, 4.0d0, 1.00d0/,
     &     si28h4 /9*100.0d0, 30.d0, 1.00d0/
      data 
     &     fe56d  /3*2.00d0, 2.50d0, 3.50d0,  3.50d0, 4*4.00d0, 1.00d0/,
     &     fe56t  /3*1.00d0, 2.00d0, 4.00d0,  4.00d0, 4*4.00d0, 1.00d0/,
     &     fe56h3 /3*1.00d0, 2.00d0, 5*4.0d0, 8.00d0, 1.00d0/,
     &     fe56h4 /9*100.0d0, 30.d0, 1.00d0/
      data 
     &     bi209d /9*2.00d0,  4.0d0, 1.00d0/,
     &     bi209t /9*8.00d0, 16.0d0, 1.00d0/,
     &     bi209h3/9*8.00d0, 16.0d0, 1.00d0/,
     &     bi209h4/9*100.d0, 100.d0, 1.0d0/
      data 
     &     u238d  /9* 3.0d0,  4.0d0, 1.00d0/,
     &     u238t  /9*12.0d0, 12.0d0, 1.00d0/,
     &     u238h3 /9*12.0d0, 12.0d0, 1.00d0/,
     &     u238h4 /9*120.d0, 100.d0, 1.0d0/

      data mq /11/

c ======================================================================

      a0 = atar0
      t0 = ener0
      if (a0.le.afit(1))     then
        ia = 1
      elseif (a0.ge.afit(5)) then
        ia = 5
      else
        do i = 1,5
        if (a0.le.afit(i)) then
          if ((a0-afit(i-1)).lt.(afit(i)-a0)) then
            ia = i - 1
          else
            ia = i
          endif
          go to 10
        endif
        end do
      endif
   10 continue
        do i = 1,mq
        xq(i) = tfit(i)
        if (ia.eq.1)     then
          yq(1,i) = al27d(i)
          yq(2,i) = al27t(i)
          yq(3,i) = al27h3(i)
          yq(4,i) = al27h4(i)
        elseif (ia.eq.2) then
          yq(1,i) = si28d(i)
          yq(2,i) = si28t(i)
          yq(3,i) = si28h3(i)
          yq(4,i) = si28h4(i)
        elseif (ia.eq.3) then
          yq(1,i) = fe56d(i)
          yq(2,i) = fe56t(i)
          yq(3,i) = fe56h3(i)
          yq(4,i) = fe56h4(i)
        elseif (ia.eq.4) then
          yq(1,i) = bi209d(i)
          yq(2,i) = bi209t(i)
          yq(3,i) = bi209h3(i)
          yq(4,i) = bi209h4(i)
        elseif (ia.eq.5) then
          yq(1,i) = u238d(i)
          yq(2,i) = u238t(i)
          yq(3,i) = u238h3(i)
          yq(4,i) = u238h4(i)
        endif
        end do
      x = t0       
c   Exclude extrapolation!     
      if (x.lt.xq(1)) then
        y(1) = yq(1,1)
        y(2) = yq(2,1)
        y(3) = yq(3,1)
        y(4) = yq(4,1)
        go to 30
      elseif (x.gt.xq(mq)) then 
        y(1) = yq(1,mq)
        y(2) = yq(2,mq)
        y(3) = yq(3,mq)
        y(4) = yq(4,mq)
        go to 30
      endif
        do i = 1,mq 
        k = i 
        if (abs(x-xq(i)).lt.1.d-10) then
         y(1) = yq(1,k)
         y(2) = yq(2,k)
         y(3) = yq(3,k)
         y(4) = yq(4,k)
         go to 30
        endif
        end do
        do k=2,mq
        if (x.ge.xq(k-1) .and. x.le.xq(k)) then
          k1 = k - 1
          k2 = k
          go to 20
        endif   
        end do 
  20  x1 = xq(k1)
      x2 = xq(k2)
      den1 = x2 - x1
        do iy=1,4
        y1(iy) = yq(iy,k1)
        y2(iy) = yq(iy,k2)
        a = (x2*y1(iy) - x1*y2(iy))/den1
        b = (y2(iy) - y1(iy))/den1
        y(iy) = a + b*x  
        end do
  30    do iy=1,4
        gbm(iy) = y(iy)
        end do
      return 
      end

c ======================================================================
			
      subroutine gambetp (gbm)

c ======================================================================
c
c    Linear interpolation of empirical gamma_beta multipliers = 
c    f(A0, t0) for proton-induced preequilibrium reactions.
c
c    Written by K. G. Gudima using the Baznat's fits, Jul.-Oct. 2004.
c    Edited by A. J. Sierk, LANL T-16,  January, 2005.
c
c ======================================================================

      implicit real*8 (a-h, o-z), integer (i-n)

      real*8 ni60d, ni60t, ni60h3, ni60h4

c ======================================================================

      common /targ0/ atar0, ztar0, ener0

      dimension gbm(4)

      dimension tfit(8), afit(7), xq(8), yq(4,8), y(4), y1(4), y2(4),  
     &         al27d(8),  al27t(8),  al27h3(8),  al27h4(8),
     &         fe54d(8),  fe54t(8),  fe54h3(8),  fe54h4(8),
     &         ni60d(8),  ni60t(8),  ni60h3(8),  ni60h4(8),
     &          y89d(8),   y89t(8),   y89h3(8),   y89h4(8), 
     &        sn120d(8), sn120t(8), sn120h3(8), sn120h4(8),
     &        au197d(8), au197t(8), au197h3(8), au197h4(8),
     &        bi209d(8), bi209t(8), bi209h3(8), bi209h4(8)

      data
     & afit /27.d0, 54.d0, 60.d0, 89.d0, 120.d0, 197.d0, 209.d0/,
     & tfit /0.028d0, 0.038d0, 0.0617d0, 0.0629d0,
     &       0.120d0, 0.160d0, 0.2000d0, 0.3000d0/
      data
     &  al27d  /1.0d0, 1.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 1.0d0/,
     &  al27t  /8.0d0, 6.0d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, 1.0d0/,
     &  al27h3 /8.0d0, 6.0d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, 1.0d0/,
     &  al27h4 /20.d0, 30.d0, 64.d0, 64.d0, 30.d0, 30.d0, 30.d0, 1.0d0/
      data 
     &  fe54d  /1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/,
     &  fe54t  /8.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 1.0d0/,
     &  fe54h3 /8.0d0, 4.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 1.0d0/,
     &  fe54h4 /24.d0, 10.d0, 10.d0, 10.d0, 5.0d0, 5.0d0, 5.0d0, 1.0d0/
      data 
     &  ni60d  /1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/,
     &  ni60t  /8.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 1.0d0/,
     &  ni60h3 /8.0d0, 4.0d0, 2.0d0, 2.0d0, 2.0d0, 4.0d0, 4.0d0, 1.0d0/,
     &  ni60h4 /24.d0, 20.d0, 16.d0, 16.d0, 8.0d0, 8.0d0, 8.0d0, 1.0d0/
      data 
     &  y89d   /1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/,
     &  y89t   /8.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 1.0d0/,
     &  y89h3  /8.0d0, 4.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/,
     &  y89h4  /24.d0, 10.d0, 8.0d0, 8.0d0, 4.0d0, 4.0d0, 4.0d0, 1.0d0/
      data 
     &  sn120d /1.5d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/,
     &  sn120t /8.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 1.0d0/,
     &  sn120h3/8.0d0, 4.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 1.0d0/,
     &  sn120h4/30.d0, 20.d0, 16.d0, 16.d0, 8.0d0, 8.0d0, 8.0d0, 1.0d0/
      data 
     &  au197d /2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 1.0d0/,
     &  au197t /12.d0, 8.0d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, 1.0d0/,
     &  au197h3/12.d0, 8.0d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, 1.0d0/,
     &  au197h4/64.d0, 30.d0, 20.d0, 20.d0, 15.d0, 15.d0, 10.d0, 1.0d0/
      data 
     &  bi209d /2.0d0, 2.0d0, 1.5d0, 1.3d0, 1.2d0, 1.2d0, 1.2d0, 1.0d0/,
     &  bi209t /12.d0, 8.0d0, 8.0d0, 6.0d0, 4.0d0, 4.0d0, 4.0d0, 1.0d0/,
     &  bi209h3/12.d0, 8.0d0, 4.0d0, 3.0d0, 3.0d0, 3.0d0, 3.0d0, 1.0d0/,
     &  bi209h4/64.d0, 64.d0, 40.d0, 30.d0, 20.d0, 20.d0, 15.d0, 1.0d0/

      data mq /8/

c ======================================================================

      a0 = atar0
      t0 = ener0
      if (a0.le.afit(1))     then
        ia = 1
      elseif (a0.ge.afit(7)) then
        ia = 7
      else
        do i = 1,7
        if (a0.le.afit(i)) then
          if ((a0-afit(i-1)).lt.(afit(i)-a0)) then
            ia = i - 1
          else
            ia = i
          endif
          go to 10
        endif
        end do
      endif
   10 continue
        do i = 1,mq
        xq(i) = tfit(i)
        if (ia.eq.1)     then
          yq(1,i) = al27d(i)
          yq(2,i) = al27t(i)
          yq(3,i) = al27h3(i)
          yq(4,i) = al27h4(i)
        elseif (ia.eq.2) then
          yq(1,i) = fe54d(i)
          yq(2,i) = fe54t(i)
          yq(3,i) = fe54h3(i)
          yq(4,i) = fe54h4(i)
        elseif (ia.eq.3) then
          yq(1,i) = ni60d(i)
          yq(2,i) = ni60t(i)
          yq(3,i) = ni60h3(i)
          yq(4,i) = ni60h4(i)
        elseif (ia.eq.4) then
          yq(1,i) = y89d(i)
          yq(2,i) = y89t(i)
          yq(3,i) = y89h3(i)
          yq(4,i) = y89h4(i)
        elseif (ia.eq.5) then
          yq(1,i) = sn120d(i)
          yq(2,i) = sn120t(i)
          yq(3,i) = sn120h3(i)
          yq(4,i) = sn120h4(i)
        elseif (ia.eq.6) then
          yq(1,i) = au197d(i)
          yq(2,i) = au197t(i)
          yq(3,i) = au197h3(i)
          yq(4,i) = au197h4(i)
        elseif (ia.eq.7) then
          yq(1,i) = bi209d(i)
          yq(2,i) = bi209t(i)
          yq(3,i) = bi209h3(i)
          yq(4,i) = bi209h4(i)
        endif
        end do
      x = t0       
c   Exclude extrapolation!     
      if (x.lt.xq(1)) then
        y(1) = yq(1,1)
        y(2) = yq(2,1)
        y(3) = yq(3,1)
        y(4) = yq(4,1)
        go to 30
      elseif (x.gt.xq(mq)) then 
        y(1) = yq(1,mq)
        y(2) = yq(2,mq)
        y(3) = yq(3,mq)
        y(4) = yq(4,mq)
        go to 30
      endif
        do i = 1,mq 
        k = i 
        if (abs(x-xq(i)).lt.1.d-10) then
         y(1) = yq(1,k)
         y(2) = yq(2,k)
         y(3) = yq(3,k)
         y(4) = yq(4,k)
         go to 30
        endif
        end do
        do k=2,mq
        if (x.ge.xq(k-1) .and. x.le.xq(k)) then
          k1 = k-1
          k2 = k
          go to 20
        endif   
        end do 
  20  x1 = xq(k1)
      x2 = xq(k2)
      den1 = x2 - x1
        do iy=1,4
        y1(iy) = yq(iy,k1)
        y2(iy) = yq(iy,k2)
        a = (x2*y1(iy) - x1*y2(iy))/den1
        b = (y2(iy) - y1(iy))/den1
        y(iy) = a + b*x  
        end do
  30  continue 
        do iy=1,4
        gbm(iy) = y(iy)
        end do
      	return
      end
		
c ======================================================================			
			
      function tkinm3 (j, exn)

c ======================================================================
c
c  Kinetic energy for particles in pre-equilibrium decay.
c  For nucleons: 
c    The particles have an energy distribution proportional to:
c    (beta + T)*(Ejmax - T)**n, where T is the kinetic energy
c    of the preequilbrium nucleon, and n = exn - 2.  The maximum of the
c    spectrum occurs for y = 1/(n+1); y = (T + beta)/(Ejmax + beta). 
c    For exn > 2, TKINM3 uses the rejection technique with the rejection
c    function = constant for n <= 5, and C*z*exp(-z) for n > 5,
c    where z = y/ybar, and beta = bn for neutrons, and -Vj(j) for 
c    protons.
c    For exn = 2, the distribution is triangular, and an exact
c    spectrum from the inversion technique is possible.
c  For composite particles: 
c    The particles have an energy distribution proportional to:
c    (T - Vj)*(T + Bj)**(aj - 1.5)*(Ejmax - T)**n, where n = 
c    exn - Aj - 1; Ejmax = uej(j) - bj(j).
c    For n = 0, the distribution is close to triangular, and a 
c    triangular function is used as the comparison function in the   
c    rejection technique.  
c    For n > 0, TKINM3 uses the rejection technique with the rejection
c    function constant for n <= 5, and z*exp(-z) for n > 5; where z = 
c    y/ybar, where Tm is the value of T at the maximum of the spectrum,
c    which is found by solving a quadratic equation for the derivative
c    of the spectrum vanishing. ybar is the corresponding value of y for 
c    the maximum of the spectrum. y = (T - Vj)/(Ejmax - Vj)
c
c   CEM95 written by S. G. Mashnik
c
c   New version using rejection written by A. J. Sierk, LANL T-2,
c   February, 1999.
c   Added new rejection function for exn small; AJS, March, 1999.
c   Fixed problem with negative energy for exn = 2 neutrons,
c   April, 2000.
c   Corrected by KKG, 11/29/01. 
c
c   "Last" change: 13-AUG-2003 by NVM
c   Modified by A. J. Sierk, October, 2003.
c
c ======================================================================

      implicit real*8 (a-h, o-z), integer (i-n)

c ======================================================================

      common /bl1005/  aj(66), ajthr(66)
      common /bl1011/  vj(6)
      common /bl1015/  rj(7), uej(6), pevapj(6)
      common /blbj/    bj(7)
      common /br0j/    r0j(6), bn
      common /c2vals/  c2(300)
c	NVM      common /counter/ icntr
      common /tkmcnt/  itkmgd, itkmtot, ifermi

      data zro, one /0.d0, 1.d0/

c ======================================================================

      itkmgd = itkmgd + 1
      if (j.eq.1) then
        dj = bn
      else
        dj = -vj(j)
      endif
      n = nint(exn - aj(j) - one)
      en = dble(n)
      enp1 = en + one
      ejmax = rj(j) + vj(j)
      if (n.lt.0) then
        tkinm3 = ejmax
        return
      endif
      if (ejmax.le.1.d-3) then
        tkinm3 = vj(j)
        return
      endif
      c1 = ejmax + dj
      if (j.le.2) then
c  Nucleon:
        if (j.eq.1) then
          c3 = dj/c1
        elseif (j.eq.2) then
          c3 = zro
        endif
        if (n.eq.0) then
    5     b1 = sqrt(rndm(-1.))
          itkmtot = itkmtot + 1
          if (b1.lt.c3) go to 5
          tkinm3 = b1*c1 - dj 
          if (tkinm3.lt.zro) go to 5
        else
          if (n.le.5) then
c   Use a constant comparison function for small n; faster than
c   gamma distribution.
   10       b1 = rndm(-1.)
            itkmtot = itkmtot + 1
            if (b1.lt.c3) go to 10
            b1c = one - b1
            phi = enp1*c2(n)*b1*b1c**n
            b2 = rndm(-1.)
            itkmtot = itkmtot + 1
            if (b2.gt.phi) go to 10
            tkinm3 = c1*b1 - dj
          else
c   Use a gamma distribution comparison function:
   20       b1 = rndm(-1.)
            b2 = rndm(-1.)
            itkmtot = itkmtot + 1
            itkmtot = itkmtot + 1
            z = -log(b1*b2)
            y = z/enp1
            if (y.lt.c3 .or. y.gt.one) go to 20
            if (n.gt.300) continue             
            phi = c2(n)*((one - y)**n)*exp(z - one)
            b3 = rndm(-1.)
            itkmtot = itkmtot + 1
            if (b3.gt.phi) go to 20
            tkinm3 = c1*y - dj
          endif
        endif
      else
c   Composite particle:
        em = aj(j) - 1.5d0
        emp1 = em + one
        mpr = nint(em - 0.5d0)
        enpm = en + em
        c4 = (bj(j) - dj)/c1
        if (n.eq.0) then
c   Use a triangular comparison function; N(y) <= Ay for
c   0 <= y <= 1. ; A selected so phi(1) = 1. y = (T - Vj)/(Ejmax - Vj)
          aux2 = one/(one + c4)
   40     b1 = sqrt(rndm(-1.))
          itkmtot = itkmtot + 1
          fct1 = max(zro, (b1 + c4)*aux2)
          phi = sqrt(fct1)
          if (mpr.gt.0) then
            do im = 1,mpr
            phi = phi*fct1
            end do
          endif
          b2 = rndm(-1.)
          itkmtot = itkmtot + 1
          if (b2.gt.phi) go to 40
          tkinm3 = c1*b1 - dj
        else
c   N > 0:
          c5 = max(zro, (bj(j) - dj)/uej(j))
          fct = uej(j)/c1
          aa = em + enp1
          bb = 0.5d0*(emp1 + enpm*c5)/aa
          cc =  em*c5/aa
          fact = bb**2 - cc
          xbar = bb + sqrt(fact)
          ybar = (xbar - c5)*fct
          ybarc = one - ybar
          if (n.le.5) then
c   Use a constant comparison function for small n; faster than
c   gamma distribution.
   60       b1 = rndm(-1.)
            b1c = one - b1
            itkmtot = itkmtot + 1
            fc2 = max(zro, (b1 + c4)/(ybar + c4))
            sqr = sqrt(fc2)
            phi = 0.99999d0*(b1/ybar)*sqr*(b1c/ybarc)**n
            if (mpr.gt.0) then
              do im = 1,mpr
              phi = phi*fc2
              end do
            endif
            b2 = rndm(-1.)
            if (b2.gt.phi) go to 60
            y = b1
          else
c   Use a gamma distribution comparison function:
            zmax = one/ybar
   80       b1 = rndm(-1.)
            b2 = rndm(-1.)
            itkmtot = itkmtot + 1
            itkmtot = itkmtot + 1
            z = -log(b1*b2)
            if (z.gt.zmax) go to 80
            y = z*ybar
            b1c = one - y
            fc2 = max(zro, (y + c4)/(ybar + c4))
            sqr = sqrt(fc2)
            phi = 0.99999d0*sqr*((b1c/ybarc)**n)*exp(z - one)
            if (mpr.gt.0) then
              do im = 1,mpr
              phi = phi*fc2
              end do
            endif
            b3 = rndm(-1.)
            itkmtot = itkmtot + 1
            if (b3.gt.phi) go to 80
          endif
          tkinm3 = c1*y - dj
        endif
      endif
      return
      end
		
c ======================================================================
      
      subroutine vhelp

c ======================================================================
c
c   Auxiliary subroutine for inverse cross section extraction by
c   interpolating the tables of Dostrovsky, et al.
c
c    CEM95 written by S. G. Mashnik
c
c    Edited by A. J. Sierk   LANL  T-2,  February, 1996.
c    Corrected for Z of daughter nucleus instead of compound
c    per SGM, March, 1999.
c    Modified by A. J. Sierk, October, 2003.
c
c ======================================================================
 
      implicit real*8 (a-h, o-z)

c ======================================================================
 
      common /bl1003/ u, a, z
      common /bl1016/ cc(6) 
      common /bl1017/ vk(6) 
      common /stmass/ z1(5), a1(5), c1(5), a2(5), c2(5)

      data zro, one, two, thr, for /0.d0, 1.d0, 2.d0, 3.d0, 4.d0/

c ======================================================================

      cc(1) = zro
      vk(1) = zro
      zz = z - one
      cc(2) = subev (zz, z1, c1, 5)
      vk(2) = subev (zz, z1, a1, 5)
      zz = zz - one
      cc(6) = subev (zz, z1, c2, 5)
      vk(6) = subev (zz, z1, a2, 5)
      cc(3) = cc(2)/two
      vk(3) = vk(2) + 0.06d0
      cc(4) = cc(2)/thr
      vk(4) = vk(2) + 0.12d0
      cc(5) = cc(6)*for/thr
      vk(5) = vk(6) - 0.06d0
      return
      end

! ======================================================================
			
      subroutine egsfill

c ======================================================================
c
c   This subroutine sets up the array of ground-state energies as a
!   function of Z, A, L from the Sierk finite-range LDM.  See the
c   routine BARFIT for references and explanations.  This array is
c   created once so cem03 can execute more efficiently in a parallel
c   implementation than when using the standalone version of the
c   code CEM03.01 as released to RSICC in 2005. It should be executed
c   once in an initialization routine and left unchanged subsequently.
c   The array egs3 is defined on the mesh L = 0,5,...,100 to reduce
c   the size of the array, and quadratic interpolation is used for
c   intermediate L values.
c   
c   Written by A. J. Sierk, LANL T-16, August, 2005.
c
c ======================================================================

      real*8 amin2, amax2, aa, bf2, bfdum, egs3, egs30, egscof, el,
     &       elmax, pa, pl, pz, xa, xl, xz, zz, one, zro
      integer i, ia, iamax, iamin, il, il5, iz, j, k, l, m, na

c ======================================================================

      common /ajsbar3/ egs3(100,106,21), iamin(100), iamax(100)

      dimension egscof(5,7,5), pa(5), pl(9), pz(7)

c  iamin and iamax defined by Audi-Wapstra mass tables for Z,N < 8;
c  by Moller-Nix mass table for Z,N > 8:

      data iamin /  1,  3,  4,  5,  7,  8, 10, 12, 14, 16, 18, 19, 21,
     &             22, 23, 24, 25, 27, 29, 30, 32, 34, 36, 38, 40, 42,
     &             44, 46, 48, 51, 53, 55, 57, 59, 61, 63, 66, 68, 70,
     &             72, 74, 77, 79, 81, 83, 86, 88, 90, 92, 94, 97, 99,
     &            101,103,106,108,110,113,115,118,120,123,125,128,130,
     &            133,136,138,141,143,146,149,151,154,156,159,162,165,
     &            167,170,173,175,178,181,184,186,189,192,195,198,200,
     &            203,206,209,212,215,218,221,224,226/

      data iamax /  6, 10, 12, 14, 19, 22, 24, 34, 38, 41, 44, 47, 51,
     &             54, 57, 60, 63, 67, 70, 73, 76, 80, 83, 86, 89, 92,
     &             96, 99,102,105,108,112,115,118,121,124,128,131,134,
     &            137,140,144,147,150,153,156,160,163,166,169,172,176,
     &            179,182,185,189,192,195,198,201,205,208,211,214,218,
     &            221,224,227,230,234,237,240,243,247,250,253,256,260,
     &            263,266,269,273,276,279,282,286,289,292,295,299,302,
     &            305,308,312,315,318,321,325,328,331/

      data ((egscof(i,j,1),i=1,5),j=1,7) 
     &  /-1.781665232d6,-2.849020290d6, 9.546305856d5, 2.453904278d5,
     &    3.656148926d5,
     &    4.358113622d6, 6.960182192d6,-2.381941132d6,-6.262569370d5,
     &   -9.026606463d5,
     &   -4.804291019d6,-7.666333374d6, 2.699742775d6, 7.415602390d5,
     &    1.006008724d6,
     &    3.505397297d6, 5.586825123d6,-2.024820713d6,-5.818008462d5,
     &   -7.353683218d5,
     &   -1.740990985d6,-2.759325148d6, 1.036253535d6, 3.035749715d5,
     &    3.606919356d5,
     &    5.492532874d5, 8.598827288d5,-3.399809581d5,-9.852362945d4,
     &   -1.108872347d5,
     &   -9.229576432d4,-1.431344258d5, 5.896521547d4, 1.772385043d4,
     &    1.845424227d4/

      data ((egscof(i,j,2),i=1,5),j=1,7) 
     & /4.679351387d6, 7.707630513d6,-2.718115276d6,-9.845252314d5,
     & -1.107173456d6,
     & -1.137635233d7,-1.870617878d7, 6.669154225d6, 2.413451470d6,
     &  2.691480439d6,
     &  1.237627138d7, 2.030222826d7,-7.334289876d6,-2.656357635d6,
     & -2.912593917d6,
     & -8.854155353d6,-1.446966194d7, 5.295832834d6, 1.909275233d6,
     &  2.048899787d6,
     &  4.290642787d6, 6.951223648d6,-2.601557110d6,-9.129731614d5,
     & -9.627344865d5,
     & -1.314924218d6,-2.095971932d6, 8.193066795d5, 2.716279969d5,
     &  2.823297853d5,
     &  2.131536582d5, 3.342907992d5,-1.365390745d5,-4.417841315d4,
     & -4.427025540d4/

      data ((egscof(i,j,3),i=1,5),j=1,7) 
     &/-3.600471364d6,-5.805932202d6, 1.773029253d6, 4.064280430d5,
     &  7.419581557d5,
     &  8.829126250d6, 1.422377198d7,-4.473342834d6,-1.073350611d6,
     & -1.845960521d6,
     & -9.781712604d6,-1.575666314d7, 5.161226883d6, 1.341287330d6,
     &  2.083994843d6,
     &  7.182555931d6, 1.156915972d7,-3.941330542d6,-1.108259560d6,
     & -1.543982755d6,
     & -3.579820035d6,-5.740079339d6, 2.041827680d6, 5.981648181d5,
     &  7.629263278d5,
     &  1.122573403d6, 1.777161418d6,-6.714631146d5,-1.952833263d5,
     & -2.328129775d5,
     & -1.839672155d5,-2.871137706d5, 1.153532734d5, 3.423868607d4,
     &  3.738902942d4/

      data ((egscof(i,j,4),i=1,5),j=1,7) 
     &/ 2.421750735d6, 4.107929841d6,-1.302310290d6,-5.267906237d5,
     & -6.197966854d5,
     & -5.883394376d6,-9.964568970d6, 3.198405768d6, 1.293156541d6,
     &  1.506909314d6,
     &  6.387411818d6, 1.079547152d7,-3.517981421d6,-1.424705631d6,
     & -1.629099740d6,
     & -4.550695232d6,-7.665548805d6, 2.530844204d6, 1.021187317d6,
     &  1.141553709d6,
     &  2.182540324d6, 3.646532772d6,-1.228378318d6,-4.813626449d5,
     & -5.299974544d5,
     & -6.518758807d5,-1.070414288d6, 3.772592079d5, 1.372024952d5,
     &  1.505359294d5,
     &  9.952777968d4, 1.594230613d5,-6.029082719d4,-2.023689807d4,
     & -2.176008230d4/

      data ((egscof(i,j,5),i=1,5),j=1,7) 
     &/-4.902668827d5,-8.089034293d5, 1.282510910d5,-1.704435174d4,
     &  8.876109934d4,
     &  1.231673941d6, 2.035989814d6,-3.727491110d5, 4.071377327d3,
     & -2.375344759d5,
     & -1.429330809d6,-2.376692769d6, 5.216954243d5, 7.268703575d4,
     &  3.008350125d5,
     &  1.114306796d6, 1.868800148d6,-4.718718351d5,-1.215904582d5,
     & -2.510379590d5,
     & -5.873353309d5,-9.903614817d5, 2.742543392d5, 9.055579135d4,
     &  1.364869036d5,
     &  1.895325584d5, 3.184776808d5,-9.500485442d4,-3.406036086d4,
     & -4.380685984d4,
     & -2.969272274d4,-4.916872669d4, 1.596305804d4, 5.741228836d3,
     &  6.669912421d3/

      data zro, one /0.d0, 1.d0/

c ======================================================================

        do iz = 1,100
        zz = dble(iz)
        na = iamax(iz) - iamin(iz) + 1
        amin2 = 1.2d0*zz + 0.01d0*zz*zz
        amax2 = 5.8d0*zz - 0.024d0*zz*zz
        xz = 0.01d0*zz
        if (iz.ge.20) then
          do ia = 1,na
          aa = dble(iamin(iz) + ia - 1)
          xa = 2.5d-3*aa
          if (aa.le.250.d0) then
            call lpoly2 (xz, pz, 7)
            call lpoly2 (xa, pa, 5)
          endif
          if (aa.gt.amin2 - 5.d0 .and. aa.lt.amax2 + 10.d0) then
            call elmaxc (zz, aa, elmax)
            egs3(iz,ia,1) = zro
              do il=2,21   
              el = 5.d0*dble(il - 1)
              xl = el/elmax
              if (xl.le.one .and. aa.le.250.d0) then
c  Use the BARFIT values for the rotating ground state energy:
                call lpoly2 (xl, pl, 9)
                egs30 = zro
                  do k = 1,5
                    do l = 1,7
                     do m = 1,5
                     egs30 = egs30 + egscof(m,l,k)*pz(l)*pa(k)*pl(2*m-1)
                     end do
                    end do
                  end do
              elseif (xl.gt.one .and. aa.le.250.d0) then
c   L > Lmax;  outside the range of the fit to GS energy:
c   Find the value of BARFIT egs(Lmax); then scale by L^2:
                egs30 = zro               
                  do k = 1,5
                    do l = 1,7
                      do m = 1,5
                      egs30 = egs30 + egscof(m,l,k)*pz(l)*pa(k)
                      end do
                    end do  
                  end do
                egs30 = egs30*(el/elmax)**2
              else
c   A > 250; outside the range of the fit to GS energy:
                il5 = 5*(il - 1)
                bfdum = bf2 (aa, zz, il5, egs30)
              endif
              il5 = 5*(il - 1)
              if (egs30.le.zro) bfdum = bf2 (aa, zz, il5, egs30) 
              egs30 = min(aa+one,egs30)
              egs3(iz,ia,il) = egs30
              end do
          else
c   A is outside the range of the fit to GS energy for this Z:
            egs3(iz,ia,1) = zro
              do il = 2,21    
              il5 = 5*(il - 1)
              bfdum = bf2 (aa, zz, il5, egs30)
              egs30 = min(aa+one,egs30)
              egs3(iz,ia,il) = egs30
              end do
          endif
          end do
          if (na.lt.106) then
            do ia = na+1,106
              do il = 1,21    
              egs3(iz,ia,il) = egs3(iz,na,il)
              end do
            end do
          endif
        else
c  Z < 20:
          na = iamax(iz) - iamin(iz) + 1
          zz = dble(iz)
            do ia = 1,na
            aa = dble(iamin(iz)) + dble(ia - 1)
            egs3(iz,ia,1) = zro
              do il = 2,21    
              il5 = 5*(il - 1)
              bfdum = bf2 (aa, zz, il5, egs30)
              egs30 = min(aa+one,egs30)
              egs3(iz,ia,il) = egs30
              end do
            end do
            do ia = na+1,106
              do il = 1,21    
              egs3(iz,ia,il) = egs3(iz,na,il)
              end do
            end do
        endif
        end do
      return
      end
! ======================================================================
			
		