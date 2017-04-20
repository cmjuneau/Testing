C  *************************************************************
C      LAST VERSION   10.01.2005 (rijmin=2.2)
      SUBROUTINE  COALES(MV,P0D,P0T,P0A)
      implicit real*8 (a-h, o-z), integer (i-n)
      COMMON /MEMORY/ PME(9,5999),IME(5,5999)
      COMMON /IDPME/ IDPME(5999)
      COMMON /PORIG/ IORI(3,5999)
      COMMON /RIJMIN/ rijmin
      
      DIMENSION ID1(200),ID2(200),IT1(50),IT2(50),IT3(50),
     *IA1(50),IA2(50),IA3(50),IA4(50),V(3),PJS(3),PJL(3)
      DIMENSION  P1L(3),P1S(3),P2L(3),P2S(3),P3L(3),P3S(3)
      DIMENSION  IND(4) 
	data cmn/0.940d0/, rijmin/2.2d0/              
c                    ! rijmin=1/p0d=5.06/(0.090GeV/c) 
      IF(MV.LE.1)   RETURN
c
c
      ND=0
             I=1
   10 CONTINUE
c      WI=PME(7,I)
       WI=1.         
      IF(IME(4,I).NE.1.OR.IME(2,I).NE.0.OR.IME(3,I).NE.0.OR.WI.LT.1.)
     &              GO TO 12
      J=I+1
   11 IF(J.GT.MV)   GO  TO  12
c     WJ=PME(7,J)
      WJ=1.      
      IF(IME(4,J).NE.1.OR.IME(2,J).NE.0.OR.IME(3,J).NE.0.OR.WJ.LT.1.)
     &              GO TO 13
      I1=I
      P1=SQRT(PME(8,I1)*(PME(8,I1)+2.0*cmn))
                                            E1=PME(8,I1)+cmn
      PJ=SQRT(PME(8,J)*(PME(8,J)+2.0*cmn))
                                          EJ=PME(8,J)+cmn
      PJL(1)=PME(4,J)
      PJL(2)=PME(5,J)
      PJL(3)=PME(6,J)
      ES = EJ+E1
      P1L(1)=PME(4,I1)
      P1L(2)=PME(5,I1)
      P1L(3)=PME(6,I1)
      V(1)=-(P1L(1)+PJL(1))/ES
      V(2)=-(P1L(2)+PJL(2))/ES
      V(3)=-(P1L(3)+PJL(3))/ES
      CALL  CINEMQ(PJL,V,PJS,CTJ,STJ,CFJ,SFJ,TJS,cmn)
c
ccc, SGM, 07/13/05
ccc      CALL  RIJMQ(I,J,RMIN)
c
c      write(16,2005) I,J,TJS,RMIN
 2005 format(1x,2I5,2(1PE11.4))
c
ccc      if(RMIN.gt.rijmin)                   GO  TO  13 ! 2.2=0.2/0.09
c
      IF(SQRT(TJS*(TJS+2.0*cmn)).GT.P0D)   GO  TO  13
      IF(ND.LT.200)   GO  TO  50
      write(16,49)  ND
   49 FORMAT(/20X,'ND=',I5)
      GO TO 14
   50 CONTINUE
      ND=ND+1
                ID1(ND)=I
                            ID2(ND)=J
      IME(2,I)=ND
                    IME(2,J)=ND
c      CALL  RIJC(I,J,P1L,PJL,2)
   12 IF(I.GT.(MV-1))   GO  TO  14
      I=I+1
              GO  TO  10
   13 J=J+1
              GO  TO  11
   14 NA=0
      IF(ND.LT.2)  GO TO 53
      ND1=ND-1
      DO 52 ID=1,ND1
      I1=ID1(ID)
                   I2=ID2(ID)
                                IF(I1.EQ.0)  GO TO 52
      IB1=ID+1
      DO 51 IB=IB1,ND
      I3=ID1(IB)
      J=ID2(IB)
      IF(I3.EQ.0)   GO  TO  51
      IQ=IME(1,I1)+IME(1,I2)+IME(1,I3)+IME(1,J)
      IF(IQ.NE.2)  GO TO 51
      P1=SQRT(PME(8,I1)*(PME(8,I1)+2.0*cmn))
      E1=PME(8,I1)+cmn
      P2=SQRT(PME(8,I2)*(PME(8,I2)+2.0*cmn))
      E2=PME(8,I2)+cmn
      P3=SQRT(PME(8,I3)*(PME(8,I3)+2.0*cmn))
      E3=PME(8,I3)+cmn
      PJ=SQRT(PME(8,J)*(PME(8,J)+2.0*cmn))
      EJ=PME(8,J)+cmn
      P1L(1)=PME(4,I1)
      P1L(2)=PME(5,I1)
      P1L(3)=PME(6,I1)
      P2L(1)=PME(4,I2)
      P2L(2)=PME(5,I2)
      P2L(3)=PME(6,I2)
      P3L(1)=PME(4,I3)
      P3L(2)=PME(5,I3)
      P3L(3)=PME(6,I3)
      PJL(1)=PME(4,J)
      PJL(2)=PME(5,J)
      PJL(3)=PME(6,J)
      ES=EJ+E1+E2+E3
      V(1)=-(P1L(1)+P2L(1)+P3L(1)+PJL(1))/ES
      V(2)=-(P1L(2)+P2L(2)+P3L(2)+PJL(2))/ES
      V(3)=-(P1L(3)+P2L(3)+P3L(3)+PJL(3))/ES
      CALL  CINEMQ(P1L,V,P1S,CT1,ST1,CF1,SF1,T1S,cmn)
      CALL  CINEMQ(P2L,V,P2S,CT2,ST2,CF2,SF2,T2S,cmn)
      CALL  CINEMQ(P3L,V,P3S,CT3,ST3,CF3,SF3,T3S,cmn)
      CALL  CINEMQ(PJL,V,PJS,CTJ,STJ,CFJ,SFJ,TJS,cmn)
c
ccc      CALL  RIJMQ(I1,J,RMIN)
ccc      if(RMIN.gt.rijmin)                   GO  TO  51 ! 2.2=0.2/0.09
ccc      CALL  RIJMQ(I2,J,RMIN)
ccc      if(RMIN.gt.rijmin)                   GO  TO  51 ! 2.2=0.2/0.09
ccc      CALL  RIJMQ(I1,I3,RMIN)
ccc      if(RMIN.gt.rijmin)                   GO  TO  51 ! 2.2=0.2/0.09
ccc      CALL  RIJMQ(I2,I3,RMIN)
ccc      if(RMIN.gt.rijmin)                   GO  TO  51 ! 2.2=0.2/0.09
c
       IF(SQRT(T1S*(T1S+2.0*cmn)).GT.P0A) GO TO 51
      IF(SQRT(T2S*(T2S+2.0*cmn)).GT.P0A)   GO  TO  51
      IF(SQRT(T3S*(T3S+2.0*cmn)).GT.P0A)   GO  TO  51
      IF(SQRT(TJS*(TJS+2.0*cmn)).GT.P0A)   GO  TO  51
      IF(NA.LT.50)  GO TO 48
      write(16,'(20X,''NA='',I5)')  NA
      GO TO 51
   48 CONTINUE
      NA=NA+1
      IA1(NA)=I1
      IA2(NA)=I2
      IA3(NA)=I3
      IA4(NA)=J
      IME(2,I1)=NA
      IME(2,I2)=NA
      IME(2,I3)=NA
      IME(2,J)=NA
      ID1(ID)=0
      ID2(ID)=0
      ID1(IB)=0
      ID2(IB)=0
c      CALL  RIJC(I1,I2,P1L,P2L,4)
c      CALL  RIJC(I1,I3,P1L,P3L,4)
c      CALL  RIJC(I1, J,P1L,PJL,4)
c      CALL  RIJC(I2,I3,P2L,P3L,4)
c      CALL  RIJC(I2, J,P2L,PJL,4)
c      CALL  RIJC(I3, J,P3L,PJL,4)
      GO  TO  52
   51 CONTINUE
   52 CONTINUE
   53 NT=0
      IF(ND.EQ.0)  RETURN
      DO  16  ID=1,ND
      I1=ID1(ID)
      I2=ID2(ID)
      IF(I1.EQ.0)  GO TO16
      DO  15   J=1,MV
c     WJ=PME(7,J)
      WJ=1.         
      IF(IME(4,J).NE.1.OR.IME(2,J).NE.0.OR.IME(3,J).NE.0.OR.WJ.LT.1.)
     &    GO TO 15
      IQ=IME(1,I1)+IME(1,I2)+IME(1,J)
      IF(IQ.EQ.0.OR.IQ.EQ.3)   GO  TO   15
      P1=SQRT(PME(8,I1)*(PME(8,I1)+2.0*cmn))
      E1=PME(8,I1)+cmn
      P2=SQRT(PME(8,I2)*(PME(8,I2)+2.0*cmn))
      E2=PME(8,I2)+cmn
      PJ=SQRT(PME(8,J)*(PME(8,J)+2.0*cmn))
      EJ=PME(8,J)+cmn
      P1L(1)=PME(4,I1)
      P1L(2)=PME(5,I1)
      P1L(3)=PME(6,I1)
      P2L(1)=PME(4,I2)
      P2L(2)=PME(5,I2)
      P2L(3)=PME(6,I2)
      PJL(1)=PME(4,J)
      PJL(2)=PME(5,J)
      PJL(3)=PME(6,J)
      ES = EJ+E1+E2
      V(1)=-(P1L(1)+P2L(1)+PJL(1))/ES
      V(2)=-(P1L(2)+P2L(2)+PJL(2))/ES
      V(3)=-(P1L(3)+P2L(3)+PJL(3))/ES
      CALL  CINEMQ(P1L,V,P1S,CT1,ST1,CF1,SF1,T1S,cmn)
      CALL  CINEMQ(P2L,V,P2S,CT2,ST2,CF2,SF2,T2S,cmn)
      CALL  CINEMQ(PJL,V,PJS,CTJ,STJ,CFJ,SFJ,TJS,cmn)
c
ccc      CALL  RIJMQ(I1,J,RMIN)
ccc      if(RMIN.gt.rijmin)                   GO  TO  15 ! 2.2=0.2/0.09
ccc      CALL  RIJMQ(I2,J,RMIN)
ccc      if(RMIN.gt.rijmin)                   GO  TO  15 ! 2.2=0.2/0.09
c
      IF(SQRT(T1S*(T1S+2.0*cmn)).GT.P0T)   GO  TO  15
      IF(SQRT(T2S*(T2S+2.0*cmn)).GT.P0T)   GO  TO  15
      IF(SQRT(TJS*(TJS+2.0*cmn)).GT.P0T)   GO  TO  15
      IF(NT.LT.50)  GO  TO  55
      write(16,'(20X,''NT='',I5)')  NT
      GO  TO  15
   55 CONTINUE
      NT=NT+1
      IT1(NT)=I1
      IT2(NT)=I2
      IT3(NT)=J
      IME(2,I1)=NT
      IME(2,I2)=NT
      IME(2,J)=NT
      ID1(ID)=0
      ID2(ID)=0
c      CALL  RIJC(I1,I2,P1L,P2L,3)
c      CALL  RIJC(I1, J,P1L,PJL,3)
c      CALL  RIJC(I2, J,P2L,PJL,3)
      GO  TO  16
   15 CONTINUE
   16 CONTINUE
   17 IF(NT.EQ.0)   GO  TO  119
      DO  19  IT=1,NT
      DO  18   J=1,MV
c     WJ=PME(7,J)
      WJ=1.       
      IF(IME(4,J).NE.1.OR.IME(2,J).NE.0.OR.IME(3,J).NE.0.OR.WJ.LT.1.)
     &      GO TO 18
      I1=IT1(IT)
      I2=IT2(IT)
      I3=IT3(IT)
      IQ=IME(1,I1)+IME(1,I2)+IME(1,I3)+IME(1,J)
      IF(IQ.NE.2)   GO  TO  18
      P1=SQRT(PME(8,I1)*(PME(8,I1)+2.0*cmn))
      E1=PME(8,I1)+cmn
      P2=SQRT(PME(8,I2)*(PME(8,I2)+2.0*cmn))
      E2=PME(8,I2)+cmn
      P3=SQRT(PME(8,I3)*(PME(8,I3)+2.0*cmn))
      E3=PME(8,I3)+cmn
      PJ=SQRT(PME(8,J)*(PME(8,J)+2.0*cmn))
      EJ=PME(8,J)+cmn
      P1L(1)=PME(4,I1)
      P1L(2)=PME(5,I1)
      P1L(3)=PME(6,I1)
      P2L(1)=PME(4,I2)
      P2L(2)=PME(5,I2)
      P2L(3)=PME(6,I2)
      P3L(1)=PME(4,I3)
      P3L(2)=PME(5,I3)
      P3L(3)=PME(6,I3)
      PJL(1)=PME(4,J)
      PJL(2)=PME(5,J)
      PJL(3)=PME(6,J)
      ES = EJ+E1+E2+E3
      V(1)=-(P1L(1)+P2L(1)+P3L(1)+PJL(1))/ES
      V(2)=-(P1L(2)+P2L(2)+P3L(2)+PJL(2))/ES
      V(3)=-(P1L(3)+P2L(3)+P3L(3)+PJL(3))/ES
      CALL  CINEMQ(P1L,V,P1S,CT1,ST1,CF1,SF1,T1S,cmn)
      CALL  CINEMQ(P2L,V,P2S,CT2,ST2,CF2,SF2,T2S,cmn)
      CALL  CINEMQ(P3L,V,P3S,CT3,ST3,CF3,SF3,T3S,cmn)
      CALL  CINEMQ(PJL,V,PJS,CTJ,STJ,CFJ,SFJ,TJS,cmn)
c
ccc      CALL  RIJMQ(I1,J,RMIN)
ccc      if(RMIN.gt.rijmin)                   GO  TO  18 ! 2.2=0.2/0.09
ccc      CALL  RIJMQ(I2,J,RMIN)
ccc      if(RMIN.gt.rijmin)                   GO  TO  18 ! 2.2=0.2/0.09
ccc      CALL  RIJMQ(I3,J,RMIN)
ccc      if(RMIN.gt.rijmin)                   GO  TO  18 ! 2.2=0.2/0.09
c
      IF(SQRT(T1S*(T1S+2.0*cmn)).GT.P0A)   GO  TO  18
      IF(SQRT(T2S*(T2S+2.0*cmn)).GT.P0A)   GO  TO  18
      IF(SQRT(T3S*(T3S+2.0*cmn)).GT.P0A)   GO  TO  18
      IF(SQRT(TJS*(TJS+2.0*cmn)).GT.P0A)   GO  TO  18
      IF(NA.LT.50)  GO TO 57
      write(16,'(20X,''NA*='',I5)')  NA
      GO TO 18
   57 CONTINUE
      NA=NA+1
      IA1(NA)=I1
      IA2(NA)=I2
      IA3(NA)=I3
      IA4(NA)=J
      IME(2,I1)=NA
      IME(2,I2)=NA
      IME(2,I3)=NA
      IME(2,J)=NA
      IT1(IT)=0
      IT2(IT)=0
      IT3(IT)=0
c      CALL  RIJC(I1,I2,P1L,P2L,4)
c      CALL  RIJC(I1,I3,P1L,P3L,4)
c      CALL  RIJC(I1, J,P1L,PJL,4)
c      CALL  RIJC(I2,I3,P2L,P3L,4)
c      CALL  RIJC(I2, J,P2L,PJL,4)
c      CALL  RIJC(I3, J,P3L,PJL,4)
      GO  TO  19
   18 CONTINUE
   19 CONTINUE
  119 IF(NA.EQ.0)   GO  TO  21
C--------------------------------------
      DO  20  IA=1,NA
      I1=IA1(IA)
      I2=IA2(IA)
      I3=IA3(IA)
      I4=IA4(IA)
      IND(1)=I1
      IND(2)=I2
      IND(3)=I3
      IND(4)=I4
      CALL  CODIRQ(IND,4)
      IME(1,I1)=2
      IME(2,I1)=0
      IME(4,I1)=4
      IDPME(I1)=0
      IORI(1,I1)=0
      IORI(2,I1)=0
      IORI(3,I1)=4
   20 CONTINUE
   21 IF(NT.EQ.0)   GO  TO  23
C--------------------------------------
      DO  22  IT=1,NT
      IF(IT1(IT).EQ.0)   GO  TO  22
      I1=IT1(IT)
      I2=IT2(IT)
      I3=IT3(IT)
      IND(1)=I1
      IND(2)=I2
      IND(3)=I3
      CALL  CODIRQ(IND,3)
      IME(1,I1)=IME(1,I1)+IME(1,I2)+IME(1,I3)
      IME(2,I1)=0
      IME(4,I1)=3
      IDPME(I1)=0
      IORI(1,I1)=0
      IORI(2,I1)=0
      IORI(3,I1)=4
   22 CONTINUE
C--------------------------------------
   23 DO  25  ID=1,ND
      IF(ID1(ID).EQ.0)   GO  TO  25
      I1=ID1(ID)
      I2=ID2(ID)
      IF((IME(1,I1)+IME(1,I2)).NE.1)   GO  TO  24
      IND(1)=I1
      IND(2)=I2
      CALL  CODIRQ(IND,2)
      IME(1,I1)=1
      IME(2,I1)=0
      IME(4,I1)=2
      IDPME(I1)=0
      IORI(1,I1)=0
      IORI(2,I1)=0
      IORI(3,I1)=4
C--------------------------------------
      GO  TO  25
   24 IME(2,I1)=0
      IME(2,I2)=0
   25 CONTINUE
      I=1
   26 IF(IME(2,I).EQ.0)   GO  TO  29
      DO  27  K=1,9
   27 PME(K,I)=PME(K,MV)
      DO  28  K=1,5
      IF(K.LE.3)  IORI(K,I)=IORI(K,MV)
   28 IME(K,I)=IME(K,MV)
      IDPME(I)=IDPME(MV)
      MV=MV-1
      IF(I.GT.MV)   GO  TO  30
      GO  TO  26
   29 IF(I.GE.MV)   GO  TO  30
      I=I+1
                GO  TO  26
   30 CONTINUE
      RETURN
      END
C
C  *************************************************************


SUBROUTINE  CODIRQ(IND,NN)
implicit real*8 (a-h, o-z), integer (i-n)
COMMON /MEMORY/ PME(9,5999),IME(5,5999)
DIMENSION  PS(3),IND(4)
data cmn/0.940d0/
PS(1)=0.
PS(2)=0.
PS(3)=0.
DO  K=1,NN
  I=IND(K)
  PI=SQRT(PME(8,I)*(PME(8,I)+2.0*cmn))
  PS(1)=PS(1)+PME(4,I)
  PS(2)=PS(2)+PME(5,I)
  PS(3)=PS(3)+PME(6,I)
ENDDO
PSM=SQRT(PS(1)**2+PS(2)**2+PS(3)**2)
I1=IND(1)
PMM=FLOAT(NN)*cmn
PME(8,I1)=SQRT(PSM*PSM+PMM*PMM)-PMM
PME(9,I1)=PMM
PME(4,I1)=PS(1)
PME(5,I1)=PS(2)
PME(6,I1)=PS(3)
PME(7,I1)=1.0d0      !!! WEIGHT
RETURN
END

C**********************************************************************
      SUBROUTINE  RIJMQ(J1,J2,RMIN)
      implicit real*8 (a-h, o-z), integer (i-n)
      COMMON /MEMORY/ PME(9,5999),IME(5,5999)
      DIMENSION  R12(3),V12(3)
      DO  K=1,3
        R12(K)=PME(K,J1)-PME(K,J2)
        V12(K)=PME(K+3,J1)/(PME(8,J1)+0.94)-
     -         PME(K+3,J2)/(PME(8,J2)+0.94)
      ENDDO
      S=R12(1)*V12(1)+R12(2)*V12(2)+R12(3)*V12(3)
      R2=R12(1)**2+R12(2)**2+R12(3)**2
      V2=V12(1)**2+V12(2)**2+V12(3)**2
c
c      write(16,2005) J1,J2,R12,V12,R2,V2,S
 2005 format(1x,2I3,9(1PE11.4))
c
      RMIN=100.
      if(S.gt.0.0001d0)  RETURN
c  kkg 12/17/04
      RMIN=SQRT(ABS(R2-(S**2)/V2))
c      RMIN=SQRT(R2)
      RETURN
      END
      
      
      
      subroutine peqemtl(n, p, h, pz, lm, ipflg, rm, pnx,
     &                   pny, pnz, bf0, ln, erotev)

c ======================================================================
c
c    Preequilibrium emission calculation; extracted from old PRECOF
c    routine, 10/08/03.
c
c   Calls: AUXL GAMAGU2 KINEMQ MOLNIX ROTORQ TKINM3 TRANS8 
c
c    Written by A. J. Sierk, LANL T-16, October, 2003.
c    Modified by K. K. Gudima, June, 2004.  
c    Edited by AJS, January, 2005.
c    Edited by S. G. Mashnik, June, 2005
c    Corrected prequilibrium angular distribution, S. G. Mashnik
c    July, 2005.
c    Modified : 02-Aug-2005 by NVM
c    Updated for MARS by K.K. Gudima, March,2007
c ======================================================================

      implicit real*8 (a-h, o-z), integer (i-n)
      real*8 molnix
      logical emitok, first

c ======================================================================

      dimension gb1(2), gj(6), p12(3), pl(3), ps(3)
c      dimension gb1(2), gj(6), pl(3), ps(3)

      common /aljz/    alj0(60,90,4)
      common /azi/     iz, in
      common /bl1003/  u, a, z
      common /bl1005/  aj(66), ajthr(66)
      common /bl1006/  zj(6)
      common /bl1009/  afj(7)
      common /bl1010/  zfj(6)
      common /bl1011/  vj(6)
      common /bl1015/  rj(7), uej(6), pevapj(6)
      common /bl1018/  afjthr(7), athrd
      common /blac/    ac
      common /blalj/   alj(6), ami(6)
      common /blbj/    bj(7)
      common /blexn/   exn
      common /constan/ thrd, twthrd
      common /gambet/  gb(6)
      common /gbzro/   gb0(300,6)
      common /kktot/   ktot
      common /mnuc/    emnucm, emnucg, emnucb, emnuct, emneut, emprot
      common /pi/      pi, twpi
      common /redmas/  emured(300,28), emuredf(300,28)
      common /redmas2/ redcom(7), redpre(6) 
      common /resid/   angmom(3), v(3), remn
      common /trec/    trec
c     kkg  06/21/04
      common /gbmul/   gbm(4)
      COMMON /MEMORY/ PMEMO(9,5999),IMEMO(5,5999)
      COMMON /IDPME/ IDPME(5999)
      COMMON/PORIG/ IORI(3,5999)
	common /pmultip/ pmult(5,9)
ccc, SGM, 07/13/05
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,t0gevl,EPS1,EPS2,
     *VPI,A1,A2,C1L,C2L,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
ccc
      data (gb1(i),i=1,2) /2*1.d0/
      data zro, one, two, thr, for /0.d0, 1.d0, 2.d0, 3.d0, 4.d0/
      data thsn /1.d3/
      data first /.true./

      save first, gb1           ! 08/02/05 by NVM

c ======================================================================

      if (first) then
         first = .false.
         gb(1) = gb1(1)
         gb(2) = gb1(2)
      endif

 10   if (p.lt.one) then
        p = p + one
        h = h + one
        n = n + 2
        go to 10
      endif
      exn = p + h
      ac0 = ac
      ip = nint(p)
      iaa = nint(a)
      if (n.le.90 .and. ip.le.60) then
        alj(1) = alj0(ip,n,1)
        alj(2) = alj(1)
        alj(3) = alj0(ip,n,2)
        alj(4) = alj0(ip,n,3)
        alj(5) = alj(4)
        alj(6) = alj0(ip,n,4)
      else
        alj(1) = p*(exn - one)
        alj(2) = alj(1)
        alj(3) = alj(2)*(p - one)*(exn - two)/two
        alj(4) = alj(3)*(p - two)*(exn - thr)/6.d0
        alj(5) = alj(4)
        alj(6) = alj(5)*(p - thr)*(exn - for)/12.d0
      endif
c   KKG 03/04/04
c   Empirical multipliers for complex particle emission:
c   MIB, KKG 06/21/04
      gb(3) = gb0(iaa,1)*gbm(1)         ! d
      gb(4) = gb0(iaa,2)*gbm(2)         ! t
      gb(5) = gb0(iaa,2)*gbm(3)         ! He3
      gb(6) = gb0(iaa,3)*gbm(4)         ! He4
      g = zro
        do l = 1,6
c   gj is the emission rate of particles of type j (n,p,d,t,3He,4He)
c   into the continuum.
        emitok = p.ge.aj(l) .and. exn.ge.(aj(l)+one) .and. pz.ge.zj(l)
     &           .and. (l.ne.1 .or. p-pz.ge.one) .and. rj(l).gt.zro
        gj(l) = zro
        if (emitok) then
          ac = 0.595d0*ami(l)
          x1234 = gamagu2(l)
***          write(*,*)' l,gamagu2 =',l,x1234
          gj(l) = redpre(l)*x1234
c   g is total preequilibrium particle emission rate:
          g = g + gj(l)
        endif
        end do
      if (n.le.0) return
      ac = ac0
c   Find rates for changing exciton number with no emission:
      call trans8 (p, h, ac, c1, c2, c3)
        c = c1 + c2 + c3
        if (c.le.zro) then
          ipflg = 1
          return
        endif
        b2 = c + g
        c11 = c1/c
        c21 = (c2 + c1)/c
   20   b1 = rndm(-1.)
        if (b1.gt.g/b2) then
c  Transition in exciton number (-2, 0, or +2)
          b3 = rndm(-1.)
          if (b3.le.c21) then
            bz = rndm(-1.)
            bz1 = z/a
            if (b3.le.c11) then
              n = n + 2
              p = p + one
              h = h + one
              fc = one
            elseif (b3.gt.c11) then
              n = n - 2
              p = p - one
              h = h - one
              fc = -one
              if (pz.le.zro) then
                ipflg = 1
                return
              endif
            endif
            if (bz.le.bz1) pz = pz + fc
            ipflg = 1
            return
          else
c  delta-n = 0 option; do not recalculate c's and g's!
            go to 20
          endif
        else
c   Preequilibrium particle is emitted:
              ppp=p             
c   gj is converted to the sum up to j of all the rates for
c   particles with index less than or equal to j.
            do j = 2,6
            gj(j) = gj(j-1) + gj(j)
            end do
c   Randomly select particle to be emitted, according to partial
c   emission rates.
          b = rndm(-1.)*g
            do j = 1,6
            if (b.le.gj(j)) then
              lm = j
              go to 30
            endif
            end do
   30     continue
          ep1 = tkinm3 (lm, exn)
          izj = nint(zj(lm))
          inj = nint(aj(lm) - zj(lm))
          if (izj.ne.0) then
            emxp = molnix (izj, inj, 2)
          else
            emxp = 8.071d0
          endif
          ep3 = thsn*emnuct*aj(lm) + emxp - zj(lm)*0.511004d0
          ktot = ktot + 1
          k = ktot
          p = p - aj(lm)
          n1 = nint(aj(lm))
          n = n - n1
          pz = pz - zj(lm)
          a = afj(lm)
          athrd = afjthr(lm)
          z = zfj(lm)
          un = a - z
          iz = nint(z)
          in = nint(un)
          iz = max(1,iz)
          in = max(1,in)
          emx = molnix (iz, in, 2)
          v2 = v(1)**2 + v(2)**2 + v(3)**2
          pm = sqrt(ep1*(ep1 + two*ep3))
c  kkg 06/21/04
          if (t0gevl.le.0.21d0) then
            ct = ctkalb(ep1,bj(lm))     ! Kalbach distribution 
          else 
ccc         ct = one - two*rndm(-1.)    ! isotropical emission
            ct = rndm(-1.)   
c   Momentum of 250 MeV/c corresponds to about 34 MeV Fermi energy.
          p2 = 250.d0*rndm(-1.)**thrd
          ct2 = one - two*rndm(-1.)
          st2 = sqrt(one - ct2*ct2)
          fi2 = twpi*rndm(-1.)
c   Factor of 1000. to convert pnx from GeV/c to MeV/c.
          f1 = thsn/ppp
          f2 = p2*st2
          p12(1) = f2*cos(fi2) + pnx*f1
          p12(2) = f2*sin(fi2) + pny*f1
          p12(3) = p2*ct2 + pnz*f1
ccc
          endif
c
          st = sqrt(one - ct*ct)
          fi = twpi*rndm(-1.)
          cf = cos(fi)
          sf = sin(fi)
c   pl is the momentum of the particle with respect to the source:
          pl(1) = pm*st*cf
          pl(2) = pm*st*sf
          pl(3) = pm*ct
          if (v2.gt.1.0d-7) then
c  Nucleus has non-zero velocity; do appropriate kinematic
c  transformation on preequlibrium particle:
c   p12 is the momentum of the emitting source in the lab frame:
cc KKG 06/21/04;  SGM 07/12/05:
            if (t0gevl.le.0.21d0) then
              ps(1) = pl(1) 
              ps(2) = pl(2) 
              ps(3) = pl(3) 
            else
              call rotorq (p12, v, pl, ps)
            endif
            call kinemq (ps, v, pl, ct, st, cf, sf, tl, ep3)
c  After kinemq, pl is the lab momentum of the emitted particle:
          else
c  Nuclear velocity essentially zero;
c  Isotropic emission of preequilibrium particle:
            tl = ep1
          endif

          PMEMO(4,k)=pl(1)/1000.
          PMEMO(5,k)=pl(2)/1000.
          PMEMO(6,k)=pl(3)/1000.
c   kkg  06/26/06
          PMEMO(7,k)=zro           !!! excit. energy of fragment
          PMEMO(8,k)=tl/1000.
          PMEMO(9,k)=ep3/1000.
          IMEMO(1,k)=int(zj(lm) + 1.e-6)
          IMEMO(2,k)=0
          IMEMO(3,k)=0
          IMEMO(4,k)=int(aj(lm) + 1.e-6)
          IMEMO(5,k)=0
          IDPME(k)=0
          IORI(1,k)=NINT(a)
          IORI(2,k)=NINT(z)
          IORI(3,k)=1
          pmult(1,LM)=pmult(1,LM)+1.  ! particle LM from preco
c
c         write(*,*) 'lm,tl,ep3=',lm,tl,ep3
c
          pnx = pnx - pl(1)/thsn 
          pny = pny - pl(2)/thsn 
          pnz = pnz - pl(3)/thsn 
          if (iz.gt.7 .and. in.gt.7) then
            remn = a*emnucb + emx/thsn
          else
            remn = a*emnuct + emx/thsn
          endif
          e = sqrt(pnx**2 + pny**2 + pnz**2 + remn**2)
          trec = (e - remn)*thsn 
          v(1) = pnx/e
          v(2) = pny/e
          v(3) = pnz/e
          fi1 = atan2 (sf, cf)
          if (fi1.lt.zro) fi1 = twpi + fi1
          almax = 0.219327d0*rm*(athrd + ajthr(lm))*
     &            sqrt(emured(iaa,n1)*(ep1 - vj(lm)))
          alp = almax*sqrt(rndm(-1.))
   40     continue
          ctr = one - two*rndm(-1.)
          str = sqrt(one - ctr*ctr)
          fir = twpi*rndm(-1.)
          cfr = cos(fir)
          sfr = sin(fir)
          ctrp = ctr*ct + str*st*(cfr*cf + sfr*sf)
c   Exclude emission inside of parent nucleus:
          if (ctrp.lt.zro) go to 40
          alpx = alp*(str*sfr*ct - st*sf*ctr)
          alpy = alp*(st*cf*ctr  - str*cfr*ct)
          alpz = alp*str*st*(cfr*sf - cf*sfr) 
          angmom(1) = angmom(1) - alpx
          angmom(2) = angmom(2) - alpy
          angmom(3) = angmom(3) - alpz
          call auxl (angmom, bf0, ln, erotev, delu)
          u = uej(lm) + delu - bj(lm) + pevapj(lm) - ep1 + erotev
        endif

c ======================================================================
      end
      

c     *****************************************************************
c
      DOUBLE PRECISION function ctkalb(e,b) 
      implicit real*8 (a-h, o-z), integer (i-n)
c     this function use the Kalbach's systematics to simmulate
c     the angular distribution for preequilibrium emission in the form:
c     f(x)dx=N*exp(a*x), where x=cos(theta), theta is angle
c     in the source system. C.Kalbach ,PRC 37(1988)2350 
c           
      eb=e+b
      a = 0.04*eb+1.8d-6*eb**3
      r1 = rndm(-1.)
      ctkalb = -1.0+log(1.0+r1*(exp(2.0*a)-1.0))/a
      return
      end

C
C  *************
C**********************************************************************
      SUBROUTINE CINEMQ(PSTAR,V,P,CT,ST,CFI,SFI,T,CM)
      implicit real*8 (a-h, o-z), integer (i-n)
C     KINEMATIC BLOCK.
      DIMENSION PSTAR(3),V(3),P(3)
      SPV = PSTAR(1)*V(1)+PSTAR(2)*V(2)+PSTAR(3)*V(3)
      V2=V(1)**2+V(2)**2+V(3)**2
      G=1.
      V2S=SQRT(ABS(1.-V2))
      IF(V2S.GT.0.D0)  G=1./SQRT(ABS(1.-V2))
      ES2=PSTAR(1)**2+PSTAR(2)**2+PSTAR(3)**2+CM**2
      ES=SQRT(ES2)
      do   K=1,3
        P(K)=PSTAR(K)+G*V(K)*(SPV*G/(G+1.)+ES)
      enddo
      E=G*(ES+SPV)
      T=E-CM
      IF(T.LE.0.)  T=1.E-6
      PM2=P(1)**2+P(2)**2+P(3)**2
      PM=SQRT(PM2)
      CT=P(3)/PM
      ST2=1.-CT*CT
      IF(ST2.GT.0.)  then
        ST=SQRT(ST2)
        CFI=P(1)/PM/ST
        SFI=P(2)/PM/ST
      ELSE
        ST=0.
        CFI=1.
        SFI=0.
      ENDIF
      RETURN
         END