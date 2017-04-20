
      function rndm (rdummy)

c ======================================================================
c                                                                      *
c     This routine merely acts as an interface to F.James rm48  gen.   *
c                                                                      *
c     Created on 03  April  1992   by    Alfredo Ferrari & Paola Sala  *
c            INFN - Milan                                              *
c                                                                      *
c     Modified    on 16-Sep-93     by    Alfredo Ferrari               *
c                                                                      *
c     Revision: 13-FEB-1997    BY NVMokhov                             *
c     Edited by A. J. Sierk, LANL T-16, September-October, 2003.       *
c                                                                      *
c ======================================================================

      implicit real*8 (a-h, o-z)
      integer cntr, idummy, ione, iseed1, iseed2

      real rdummy

c ======================================================================

      parameter (pipipi = 3.1415926535897932270d0)

      dimension rndnum (2)
      common /countrn/ cntr

      data ione /1/

c ======================================================================

      call rm48 (rndnum, ione)
      rndm = rndnum (1)
      cntr = cntr + 1
      return

      entry rd2in (iseed1, iseed2)
c  The following card just to avoid warning messages on the HP compiler
      rd2in = pipipi
      call rm48in (54217137, iseed1, iseed2)
      return

      entry rd2out (iseed1, iseed2)
c  The following card just to avoid warning messages on the HP compiler
      rd2out = pipipi
      call rm48ut (idummy, iseed1, iseed2)
      return

c ======================================================================
      end

      subroutine rm48 (rvec, lenv)

c ======================================================================
c
c     Double-precision version of
c    universal random number generator proposed by Marsaglia and Zaman
c    in report FSU-SCRI-87-50
c      Based on RANMAR, modified by F. James, to generate vectors
c      of pseudorandom numbers rvec of length lenv, where the numbers
c      in rvec are numbers with at least 48-bit mantissas.
c                                                                      *
c     Revision: 21-JAN-1997    BY NVMokhov                             *
c     Modified by A. J. Sierk, LANL T-16, September-October, 2003.
c                                                                      *
c   Input and output entry points: rm48in, rm48ut.
c!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c!!!  Calling sequences for rm48 :                                    ++
c!!!      call rm48  (rvec, len)    returns a vector rvec of len      ++
c!!!      64-bit random floating point numbers between                ++
c!!!      zero and one.                                               ++
c!!!      call rm48in (i1,n1,n2)  initializes the generator from one  ++
c!!!      64-bit integer i1, and number counts n1,n2                  ++
c!!!     (for initializing, set n1=n2=0, but to restart               ++
c!!!       a previously generated sequence, use values                ++
c!!!       output by rm48ut)                                          ++
c!!!      call rm48ut (i1,n1,n2)  outputs the value of the original   ++
c!!!     seed and the two number counts, to be used                   ++
c!!!     for restarting by initializing to I1 and                     ++
c!!!     skipping n2*100000000+n1 numbers.                            ++
c!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      implicit real*8 (a-h, o-z), integer (i-n)

      character*15 chars1
      character*18 chars2
      character*22 chars3

c ======================================================================

      parameter (modcns=1000000000)

      dimension rvec(*)

      common /r48st1/ u(97), c, i97, j97

      save cd, cm, twom24,  zero, one, half, ntot, ntot2, ijkl

      data ntot, ntot2, ijkl /-1, 0, 0/
      data lunout /16/
      data one, half, zero /1.d0, 0.5d0, 0.d0/

c ======================================================================

      if (ntot.ge.0) go to 20
c
c   Default initialization. User has called rm48  without rm48in.
      ijkl = 54217137
      ntot = 0
      ntot2 = 0
      lenv = max (lenv, 1)
      kalled = 0
      go to 10
c
      entry rm48in (ijklin, ntotin, ntot2n)
c    Initializing routine for rm48 , may be called before
c    generating pseudorandom numbers with rm48 .   The input
c    values should be in the ranges:  0<=ijklin<=900 OOO OOO
c           0<=ntotin<=999 999 999
c           0<=ntot2n<<999 999 999!
c To get the standard values in Marsaglia's paper, ijklin=54217137
c     ntotin,ntot2n=0
      ijkl = ijklin
      ntot = max(ntotin, 0)
      ntot2= max(ntot2n, 0)
      lenv = 1
      kalled = 1
c   Always come here to initialize:
   10 continue
      ij = ijkl/30082
      kl = ijkl - 30082*ij
      i = mod(ij/177, 177) + 2
      j = mod(ij, 177)     + 2
      k = mod(kl/169, 178) + 1
      l = mod(kl, 169)
      chars2 = ' RM48 initialized:'
cc    write (lunout, '(a,i10,2x,2i10)')
cc   &       ' RM48  initialized:', ijkl, ntot, ntot2
      write (lunout, 1000) chars2, ijkl, ntot, ntot2
        do ii = 1, 97
        x = zero
        t = half
          do jj = 1, 48
          m = mod(mod(i*j,179)*k, 179)
          i = j
          j = k
          k = m
          l = mod(53*l + 1, 169)
          if (mod(l*m,64).ge.32)  x = x + t
          t = half*t
          end do
        u(ii) = x
        end do
      twom24 = one
        do i24 = 1, 24
        twom24 = half*twom24
        end do
      c  =   362436.d0*twom24
      cd =  7654321.d0*twom24
      cm = 16777213.d0*twom24
      i97 = 97
      j97 = 33
c   Complete initialization by skipping
c            (ntot2*modcns + ntot) random numbers
        do loop2 = 1, ntot2+1
        now = modcns
        if (loop2.eq.ntot2 + 1) now = ntot
        if (now.gt.0) then
          chars3 = ' Rm48in skipping over '
cc        write (lunout,'(a,i15)') ' RM48IN skipping over ', now
          write (lunout, 1400) chars3, now
            do idum = 1, ntot
            uni = u(i97) - u(j97)
            if (uni.lt.zero) uni = uni + one
            u(i97) = uni
            i97 = i97 - 1
            if (i97.eq.0) i97 = 97
            j97 = j97 - 1
            if (j97.eq.0) j97 = 97
            c = c - cd
            if (c.lt.zero) c = c + cm
            end do
        endif
        end do
      if (kalled.eq.1)  return
c
c    Normal entry to generate lenv random numbers:
   20 continue
        do ivec = 1, lenv
        uni = u(i97) - u(j97)
        if (uni.lt.zero) uni = uni + one
        u(i97) = uni
        i97 = i97 - 1
        if (i97.eq.0) i97 = 97
        j97 = j97 - 1
        if (j97.eq.0) j97 = 97
        c = c - cd
        if (c.lt.zero) c = c + cm
        uni = uni - c
        if (uni.lt.zero) uni = uni + one
        rvec(ivec) = uni
        end do
      ntot = ntot + lenv
      if (ntot.ge.modcns) then
        ntot2 = ntot2 + 1
        ntot = ntot - modcns
      endif
      return

c  Entry to output current status
      entry rm48ut (ijklut, ntotut, ntot2t)
      ijklut = ijkl
      ntotut = ntot
      ntot2t = ntot2
      chars1 = ' RM48ut output:'
cc    write (lunout, '(//a,i10,2x,2i10)')
cc   & ' RM48ut output:', ijkl, ntot, ntot2
c     write (lunout, 1000) chars1, ijkl, ntot, ntot2
      return

c   Output routine for RM48 , without skipping numbers:
      entry rm48wr (ioseed)
cc    write (ioseed, '(2i10)') ntot, ntot2
cc    write (ioseed, '(2i10,f18.16)') i97, j97, c
cc    write (ioseed, '(24(4z16,/),z16)') u
      write (ioseed, 1100) ntot, ntot2
      write (ioseed, 1200) i97, j97, c
      write (ioseed, 1300) u
      return

c   Initializing routine for RM48 , without skipping numbers:
      entry  rm48rd (ioseed)
cc    read (ioseed, '(2i10)') ntot, ntot2
cc    read (ioseed, '(2i10,f18.16)') i97, j97, c
cc    read (ioseed, '(24(4z16,/),z16)') u
      write (ioseed, 1100) ntot, ntot2
      write (ioseed, 1200) i97, j97, c
      write (ioseed, 1300) u
      close (unit = ioseed)
      ijkl = 54217137
      ij = ijkl/30082
      kl = ijkl - 30082*ij
      i = mod(ij/177, 177) + 2
      j = mod(ij, 177)     + 2
      k = mod(kl/169, 178) + 1
      l = mod(kl, 169)
      chars2 = ' RM48 initialized:'
cc    write (lunout,'(a,i10,2x,2i10)')
cc   &  ' RM48  initialized:', ijkl, ntot, ntot2
      write (lunout,1000) chars2, ijkl, ntot, ntot2
      twom24 = one
        do i24 = 1, 24
        twom24 = half*twom24
        end do
      cd =  7654321.d0*twom24
      cm = 16777213.d0*twom24
      return

c ======================================================================
 
 1000 format (//a18,i10,2x,2i10)
 1100 format (2i10)
 1200 format (2i10,f18.16)
 1300 format (24(4z16,/),z16)
 1400 format (a22,i15)

c ======================================================================
      end
*****************************************************************

      SUBROUTINE RDMINI
*   Last change: 17-DEC-2003 by NVM
      IMPLICIT INTEGER (I-N)

      ISEED1=54217137
      ISEED2=0
      ISEED3=0
      CALL RM48IN (ISEED1,ISEED2,ISEED3)
      RETURN
      END
C-----------------------------------------------------------------
c ******************************************
      SUBROUTINE RDMIN (ISEED1,ISEED2,ISEED3)
      IMPLICIT INTEGER (I-N)
      if(ISEED1.lt.0.or.ISEED1.gt.900000000)  ISEED1=54217137
      if(ISEED2.lt.0.or.ISEED2.gt.999999999)  ISEED1=0
      if(ISEED3.lt.0.or.ISEED3.gt.900000000)  ISEED1=0
      CALL RM48IN (ISEED1,ISEED2,ISEED3)
c     write(*,*) ' ISEED1,2,3=',ISEED1,ISEED2,ISEED3
      RETURN
c ******************************************
      ENTRY RDMOUT(ISEED1,ISEED2,ISEED3)
      CALL RM48UT (ISEED1,ISEED2,ISEED3)
      RETURN
      END
*************************************************************

