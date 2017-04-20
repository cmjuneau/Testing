C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ELEXQ(V,U,TIN1,PARTIN,IPATIN,IPATNE,
     *MV,NP,L,MS,MQ,KSI,ME)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C     BLOCK OF CALCULATION OF PARTICLE CHARACTERISTICS IN ELASTIC AND
C     CHARGE EXCHANGE SCATTERING.
c     Calls: PINDEL,cduarteq, COSELQ, STREX, ABELQ  
c
      DIMENSION PMEMO(9,5999),IMEMO(5,5999),V(3),PARTIN(9),IPATIN(5),
     1IPATNE(5),PIST(3),PNST(3)
      COMMON/MEMORY/PMEMO,IMEMO
      COMMON/ISOB2/ISOB2
      COMMON/NCASCA/NCAS,NCPRI
      IS=IPATIN(3)
      NS=IPATNE(3)
      IF (IPATIN(2)) 10,11,10
   10 CMI = 0.14
      CMN = 0.94
                                GO TO 12
   11 CMI = PARTIN(9)
      CMN = 0.94
                                GO TO 12
   12 SEX=CROSEG(L,MS,MQ,KSI,2,TIN1,PARTIN(9),IPATIN(5))
      SEL=CROSEG(L,MS,MQ,KSI,1,TIN1,PARTIN(9),IPATIN(5))
      BETAEX=SEX/(SEX+SEL)
      DRND=RNDM(-1.)
      IF(DRND-BETAEX)14,13,13
   13 IE = IPATIN(1)
      NE = IPATNE(1)
c  kkg 10/28/03
      IF(ISOB2.EQ.0.or.L.ne.0)  GO  TO  113
      ISO=0
      IF(IPATIN(4).EQ.0.AND.IPATIN(3).EQ.0.AND.IPATIN(5).EQ.0)
     *CALL  PINDEL(TIN1,U,MV,IE,NE,NP,ISO)
      IF(ISO.EQ.1)  RETURN
c  kkg 10/22/03
  113 CONTINUE
      if(MQ.eq.2.and.TIN1.lt.2.0)  then
         CTSTI=cduarteq(TIN1, IPATIN(1), IPATNE(1))
      else
         CTSTI = COSELQ(L,MQ,KSI,TIN1,PARTIN(9))
      endif
      GO TO 18
c  kkg 10/22/03
   14 IF(IPATIN(5).NE.0.AND.IPATIN(4).EQ.1)    GO  TO  117
      IF(IPATIN(3).EQ.0)  GO  TO  114
      IF(IPATIN(3).NE.0)
     *CALL STREX(IPATIN(3),IPATIN(1),IPATNE(1),CMI,CMN,IS,NS,IE,NE)
      GO  TO  17
  114 IF (IPATIN(1)) 15,16,15
   15 IE = 0
      NE = ME
                         GO TO 17
   16 NE = 1-IPATNE(1)
      IE = ME-NE
                         GO TO 17
   17 CTSTI = COSEXQ (L,TIN1,PARTIN(9))
                         GO TO 18
  117 NE=1-IPATNE(1)
      IE=ME-NE
      GO  TO  113
   18 FISTI = 0.
      CALL ABELQ(PARTIN,V,U,PIST,PNST,CTSTI,FISTI,CMI,CMN)
      IF(MV-5997) 19,19,20
   20 NP = 0
      write(16,21)
   21 FORMAT (45X,29H MEMORY IS EXCEEDED IN CASCAD)
      RETURN
   19 CONTINUE
      PMEMO(1,MV+3)=0.
      PMEMO(2,MV+3)=0.
      PMEMO(3,MV+3)=0.
      PMEMO(4,MV+3) = PIST(1)
      PMEMO(5,MV+3) = PIST(2)
      PMEMO(6,MV+3) = PIST(3)
      PMEMO(7,MV+3)=0.
      PMEMO(9,MV+3)=CMI
      IMEMO(1,MV+3)=IE
      IMEMO(2,MV+3)=0
      IMEMO(3,MV+3) = IS
      IMEMO(4,MV+3) = IPATIN(4)
      IMEMO(5,MV+3) = IPATIN(5)
      PMEMO(1,MV+1)=0.
      PMEMO(2,MV+1)=0.
      PMEMO(3,MV+1)=0.
      PMEMO(4,MV+1) = PNST(1)
      PMEMO(5,MV+1) = PNST(2)
      PMEMO(6,MV+1) = PNST(3)
      PMEMO(7,MV+1)=0.
      PMEMO(9,MV+1)=CMN
      IMEMO(1,MV+1)=NE
      IMEMO(2,MV+1)=0
      IMEMO(3,MV+1) = NS
      IMEMO(4,MV+1) = 1
      IMEMO(5,MV+1) = IPATNE(5)
      NP = 2
c
      if(NCAS.gt.NCPRI)  then
        write(16,100) MV+1, TIN1, PIST,CMN,(IMEMO(k,MV+1),k=1,5),
     &                MV+3, CTSTI,PNST,CMI,(IMEMO(k,MV+3),k=1,5)
        write( *,100) MV+1, TIN1, PIST,CMN,(IMEMO(k,MV+1),k=1,5),
     &                MV+3, CTSTI,PNST,CMI,(IMEMO(k,MV+3),k=1,5)
  100   format(' ELEXQ:',I5,5E11.4,4I3,I10/
     &               6x,I5,5E11.4,4I3,I10)
      endif 
c
      RETURN
         END
c   
c*********************************************************************
c
      double precision function cduarteq(tin1,ie1,ie2)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c     angular distribution simulation for n+p and p+p
c     for energies < 2GeV using Duarte's approximations:
c     n+p : ds/dom ~ exp(B*t)+a*exp(B*u)+c*exp(alc*u)
c     P+p:  ds/dom ~ exp(B*t)+a*exp(B1*t)
c     with Mandelstam's variables t,u
c     H.Duarte http://www.fjfi.cvut.cz/con_adtt99/papers/Mo-o-c17.pdf
      real*8 mn
      data mn/0.939/
      data  pi /3.1415926536d0/
c
      elab=tin1*1000.
      ps2=tin1*mn/2.
      if(ie1.eq.ie2)   then      ! p +p or n + n
        if(elab.le.300.)      then
          B = 0.
          a = 0.2
          B1= 0.
        elseif(elab.le.670.)  then
          B = 9.87E-8*(elab-300.)**3
          a = 0.2
          B1= 0.
        elseif(elab.le.1100.)  then
          B = 4.56E-3*(elab-670.)+4.76
          a = 97.02E+3*EXP(-2.0E-2*elab)+0.053
          B1= 9.72E-8*EXP(-5.0E-3*(elab-670.))*(elab-670.)**3  
        else
          B = 7.4/(1.+3.0E+5/(elab-300.)**2.23)
          a = 0.28*EXP(-1.5E-3*elab)
          B1= 1.94*EXP(-7.0E-4*elab)   
        endif 
        if(B.le.0.0001)  then
          CN1 = 4.*pi 
        else
          CN1 = 2.*pi*  (1.-EXP(-2.*ps2*B ))/(ps2*B) 
        endif
        if(B1.le.0.0001) then
          CN2 = 4.*pi*a
        else
          CN2=  2.*pi*a*(1.-EXP(-2.*ps2*B1))/(ps2*B1)
        endif
        CN = CN1+CN2
        C1 = CN1/CN
        C2 = CN2/CN 
        r1=RNDM(-1.)
        tm=-2.0*ps2  
        if(r1.le.C1)   then
          r2=RNDM(-1.)
          if(B.le.0.0001)  then
            t= tm*r2
          else
            t=LOG(1.0-r2*(1.0-EXP(tm*B)))/B
          endif
        else            
          r2=RNDM(-1.)
          if(B1.le.0.0001)  then
            t= tm*r2
          else
            t=LOG(1.0-r2*(1.0-EXP(tm*B1)))/B1
          endif
        endif 
        cts=1.0-t/tm
        if(RNDM(-1.).le.0.5) then
          cduarteq = cts   
        else
          cduarteq =-cts
        endif 
c
      elseif((ie1.eq.0.and.ie2.eq.1).            ! n + p
     &    or.(ie1.eq.1.and.ie2.eq.0))   then     ! p + n
        B = 25.0*EXP(-6.0E-3*elab)+2.0E-3*elab+2.8
        a = (3.0E-3*elab+0.1)*EXP(1.55-4.9E-3*elab)+1.0E-4*elab
        c = 6.0E-5*elab*elab*EXP(-6.0E-3*elab)+15.0E-5*elab
        if(elab.lt.100.)  alc=330.-elab 
        if(elab.ge.100.)  alc=80.+15.0E+3/elab
        CN1 = pi/(ps2*B)*(1.-EXP(-4.*ps2*B))  
        CN2 = CN1*a
        CN3 = pi/(ps2*alc)*c*(1.-EXP(-4.*ps2*alc))
        CN = CN1+CN2+CN3
        C1=CN1/CN
        C2=CN2/CN
        C3=CN3/CN 
        tm=-4.0*ps2
        um=tm  
        r1 = RNDM(-1.)
        if(r1.le.C1)         then
          r2=RNDM(-1.)
          if(B.le.0.0001)     then
            t= tm*r2
          else
            t=LOG(1.0-r2*(1.0-EXP(tm*B)))/B
          endif
          cts = 1.0-2.0*t/tm 
        elseif(r1.le.(C1+C2)) then   
          r2=RNDM(-1.)
          if(B.le.0.0001)     then
            u= um*r2
          else
            u=LOG(1.0-r2*(1.0-EXP(um*B)))/B
          endif
          cts =-1.0+2.0*u/um
        else  
          r2=RNDM(-1.)
          if(alc.le.0.0001)     then
            u= um*r2
          else
            u=LOG(1.0-r2*(1.0-EXP(um*alc)))/alc
          endif
          cts =-1.0+2.0*u/um
        endif
        if(ie1.eq.1.and.ie2.eq.0) then
c  parametrization is for n + p case        
          cduarteq=-cts
        else
          cduarteq= cts
        endif  

      else
        write(*,*) ' cduarteq: ie1 and ie2 are wrong, ie1,ie2=',ie1,ie2
        stop
      endif
      if(ABS(cduarteq).gt.1.0)         then
        cduarteq=sign(1.d0,cduarteq)
      elseif(ABS(cduarteq).lt.1.0d-10) then
        cduarteq=0.0
      endif
      return
      end   
c   
c*********************************************************************
c
      
      SUBROUTINE STREX(IS1,IE1,NE1,CMI,CMN,IS,NS,IE,NE)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Calculation of  channels a) charge exchange K + N ==> K' + N'   
c     b) strange exchange AKA + N ==> pi + Y,
c     where K=(K+,K0), AKA=(K-,AK0)
c 
      IF(IS1.LT.0)  GO  TO  10
      IE=1-IE1
    9 IS=IS1
      NS=0
      CMI=0.492
      CMN=0.940
      GO  TO  23
   10 IF((IE1+NE1).EQ.0)  GO  TO  17
      RND4=4.*RNDM(-1.)
      NBR=INTG(RND4)+1
      GO  TO  (11,12,13,14,14),NBR
   11 IE=1
      GO  TO  15
   12 IE=0
      GO  TO  16
   13 IE=1
      GO  TO  16
   14 IE=-1
      GO  TO   9
   15 IS=0
      NS=-1
      CMI=0.140
      CMN=1.116
      GO  TO  23
   16 IS=0
      NS=-1
      CMI=0.140
      CMN=1.189
      GO  TO  23
   17 RND5=5.*RNDM(-1.)
      NBR=INTG(RND5)+1
      GO  TO  (18,19,20,21,22,22),NBR
   18 IE=0
      GO  TO  15
   19 IE=0
      GO  TO  16
   20 IE=-1
      GO  TO  16
   21 IE=1
      GO  TO  16
   22 IF(IE1.EQ. 0)  IE=-1
      IF(IE1.EQ.-1)  IE= 0
      GO  TO  9
   23 NE=IE1+NE1-IE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PINDEL(T1,DM,MV,IE,NE,NP,ISO)
      IMPLICIT REAL*8(A-H,O-Z)
c
c   forms the Delta from pi+N system
c
      COMMON/MEMORY/PME(9,5999),IME(5,5999)
      ISO=0
      IF(T1.GT.0.5)                  RETURN
      IF(DM.LE.1.082)                RETURN
      CALL  WMD(DM,TAU0,FMD)
      DRND=RNDM(-1.)
      IF(DRND.GT.FMD) RETURN
      PME(1,MV+1)=0.
      PME(2,MV+1)=0.
      PME(3,MV+1)=0.
      PME(4,MV+1)=0.
      PME(5,MV+1)=0.
      PME(6,MV+1)=0.
      PME(7,MV+1)=0.
      PME(9,MV+1)=DM
      IME(1,MV+1)=IE+NE
      IME(2,MV+1)=0
      IME(3,MV+1)=0
      IME(4,MV+1)=1
      IME(5,MV+1)=INTG(1000.*TAU0)
      IF(IME(5,MV+1).EQ.0)  IME(5,MV+1)=1
      NP=1
      ISO=1
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  TMATUR(NP,V,M,P1,P2,MV,PS,TL,TMAT)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c  calculation of formation (maturity) time(fm/c) of produced particle
c   
      REAL*8  MMAT,MMES,MBAR,MLID
      COMMON/MMATUR/MMES,MBAR
      COMMON/INTTYP/ ITYP
      COMMON/ITHEA/ITHEA(11)
       COMMON/MEMORY/PMEMO(9,5999),IMEMO(5,5999)
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      COMMON/IACT/ IACT
      COMMON/CSLID/CLIDER(5999)
      COMMON/CMALI/CMALI
      COMMON/ACTIM/TINT
      DIMENSION V(3),P1(9),P2(9),PS(3),RL(3)
C
C     ITYP=1(ABSORP),2(ELEX),3(BINEL),4(HEINEN)
C
      TMAT=0.
      IF(NP.LE.2.OR.ITYP.LE.3.or.
     &(ITYP.EQ.4.AND.ITHEA(4).EQ.1))  THEN
c   formation time TMAT=0 for  2-body channel
c                              absorption
c                              elastic scattering
c                              low energy inelastic collision
c                              high energy elastic collision
        CLIDER(M)=1.
      ELSE
        IF(IACT.LE.2) then
          MLID=MBAR
          IF(IMEMO(4,M).EQ.0)              MMAT=MMES
          IF(IMEMO(4,M).NE.0)              MMAT=MBAR
C  FOR ANTI-BARYON
          IF(IMEMO(4,M).LT.0)              MMAT=MMES
C
cc          IF(CLIDER(M).GT.0.)              MMAT=MLID !!12.10.97/12.03.98
C
          TAU0=MMAT/(5.06*PMEMO(9,M))  !! 05.06.96,  04/11/07
c         TAU0=MMAT/(5.06*0.14)        !! 14.09.97
          GL=1.+PMEMO(8,M)/PMEMO(9,M)
          V2=V(1)**2+V(2)**2+V(3)**2
          U=(P1(8)+P1(9)+P2(8)+P2(9))*SQRT(1.-V2)
          X=ABS(PS(3))/(U/2.)
c         FX=FXTMAT(X)
          FX=1.0d0
          TAU0L=GL*TAU0*FX
          BRND=RNDM(-1.)
          AF=-LOG(BRND)               !! 05.06.96
cc          AF=1.                        !! 27.02.97
          TMAT=TAU0L*AF
        else
c  For IACT=3
          RTAU=TL-TINT
          TMAT=TL
          IF(CLIDER(M).GT.CMALI)   THEN
            TMAT=0.
            CLIDER(M)=1.
          ENDIF
        endif
        if(TMAT.LT.0.) then
          write(16,100) ITYP,M,IMEMO(5,M),(PMEMO(K,M),K=8,9),X,TMAT
  100     FORMAT(2X,'TMATUR:ITYP=',I2,2X,'M=',I5,2X,
     &    'IM(5)=',I15,2X,'T=',F8.3,2X,'MASS=',F8.3,2X,
     &    'X=',E13.6,2X,'TMAT=',E13.6)
          TMAT=0.
        endif
      ENDIF
c     CALL  HTFORL(M,TMAT)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE KINEMR(V,M,RL,TL)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c   calculation of proper and observer's system time of
c   produced particle M
c
      COMMON/MEMORY/PME(9,5999),IME(5,5999)
      DIMENSION V(3),RL(3)
C PROPER TIME
      TP2=PME(7,M)**2-PME(1,M)**2-PME(2,M)**2-PME(3,M)**2
      TP=SQRT(ABS(TP2))
C TIME IN OBSERVER SYSTEM
      G=(PME(8,M)+PME(9,M))/PME(9,M)
      TL=G*TP
      RL(1)=0.
      RL(2)=0.
      RL(3)=0.
c     CALL  HTFORP(M)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE VMNSPQ (PARTIN,IPATIN,U,MV,NP,ITH,MQ,TIN1,LP)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Calculate secondary particle number and determine absolute 
c     values of momenta in inelastic (pion production) interaction.
c
c   Called by: BINELQ
c
c   Calls: JTYPBQ PMOMQ
c
c
      DIMENSION PMEMO(9,5999),IMEMO(5,5999),PARTIN(9),IPATIN(5)
      COMMON/MEMORY/PMEMO,IMEMO
      KH = 0
      LP = 0
      AM3=PARTIN(9)
      IF(IPATIN(5).NE.0.AND.IPATIN(4).EQ.1)  AM3=0.940
   10 U1 = U
      LAMBDA = 1
   11 LTEMP = MV+LAMBDA
      IF(LTEMP-5999) 37,37,38
   38 NP = 0
      write(16,39)
   39 FORMAT (25X,'MEMORY IS EXCEEDED IN CASCAN')
      RETURN
   37 CONTINUE
      IF (LAMBDA-1) 12,12,13
   12 PMEMO(9,MV+1)=0.94
      IMEMO(2,MV+1)=0
      IMEMO(3,MV+1)=0
      IMEMO(4,MV+1) = 1
      IMEMO(5,MV+1) = 0
                          GO TO 16
   13 IF (LAMBDA-3) 15,14,15
   14 PMEMO(9,MV+3) = AM3
      IMEMO(2,MV+3) = 0
      IMEMO(3,MV+3) = IPATIN(3)
      IMEMO(4,MV+3) = IPATIN(4)
      IMEMO(5,MV+3) = 0
                                  GO TO 16
   15 PMEMO(9,LTEMP) = 0.14
      IMEMO(2,LTEMP) = 0
      IMEMO(3,LTEMP) = 0
      IMEMO(4,LTEMP) = 0
      IMEMO(5,LTEMP) = 0
                                                GO TO 16
   16 JB = JTYPBQ(ITH,MQ,LAMBDA)
      PMEMO(8,LTEMP) = PMOMQ(JB,TIN1)
      EL = SQRT (PMEMO(8,LTEMP)**2+PMEMO(9,LTEMP)**2)
      DELTU = U1-EL
      IF (LAMBDA-2) 23,17,23
   17 IF (DELTU-AM3)       35,35,18
   18 IF (ITH) 19,22,19
   19 PMEMO(8,MV+3) = SQRT (DELTU**2-AM3**2)
      PMEMO(9,MV+3) = AM3
      IMEMO(2,MV+3) = 0
      IMEMO(3,MV+3) = IPATIN(3)
      IMEMO(4,MV+3) = IPATIN(4)
      IMEMO(5,MV+3) = 0
      IF (PMEMO(8,MV+1)-PMEMO(8,MV+2)-PMEMO(8,MV+3)) 20,20,35
   20 IF (PMEMO(8,MV+1)-ABS(PMEMO(8,MV+2)-PMEMO(8,MV+3))) 35,35,21
   21 NP = 3
      RETURN
   22 U1 = DELTU
      LAMBDA = LAMBDA+1
                                       GO TO 11
   23 IF (DELTU-0.14) 24,24,22
   24 IF (LAMBDA-1) 35,35,25
   25 IF (LAMBDA-3) 26,35,26
   26 EL=DELTU+EL
      PMEMO(8,LTEMP) = SQRT (EL**2-PMEMO(9,LTEMP)**2)
      NP = LAMBDA
      I = 1
   28 ITEMP = MV+I
      C = PMEMO(8,ITEMP)
   29 IF (NP-I) 34,34,30
   30 IF (C-PMEMO(8,ITEMP+1)) 32,31,31
   31 I=I+1
      ITEMP=ITEMP+1
                GO  TO  29
   32 I = I+1
                GO TO 28
   34 PMAX = C
      SIGMA = 0.
      DO 33 I=1,NP
      ITEMP = MV+I
      SIGMA = SIGMA+PMEMO(8,ITEMP)
   33 CONTINUE
      IF (2.*PMAX-SIGMA) 27,35,35
   27 RETURN
   35 KH = KH+1
      IF (KH-100) 10,36,36
   36 LP = 2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE CHINELQ (IPATIN,L,MS,MQ,KSI,NP,MV,TIN1,ME,IPATNE,AMIN)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Determinine secondary particles' charges in inelastic (pion prod.)
c     reaction.
c
c   Called by: BINELQ
c
c   Calls: CROSEG
c
      DIMENSION IPATIN(5),IPATNE(5),PMEMO(9,5999),IMEMO(5,5999)
      COMMON/MEMORY/PMEMO,IMEMO
      IF(IPATIN(5).NE.0.AND.IPATIN(4).EQ.1)  GO  TO  41
      IF (NP-3) 21,10,21
   10 SPI0 = CROSEG (L,MS,MQ,KSI,4,TIN1,AMIN,IPATIN(5))
      STH = CROSEG (L,MS,MQ,KSI,7,TIN1,AMIN,IPATIN(5))
      BPI0 = SPI0/STH
      BPIEX = (SPI0+CROSEG(L,MS,MQ,KSI,5,TIN1,AMIN,IPATIN(5)))/STH
      TEMP1 = RNDM(-1.)
      IF (TEMP1-BPI0) 19,11,11
   11 IF (TEMP1-BPIEX) 20,12,12
   12 IMEMO(1,MV+1)=IPATNE(1)-(IPATIN(4)-1)*IPATIN(1)
      IMEMO(1,MV+3)=(IPATIN(4)-1)*IPATIN(1)**2-IPATIN(4)*IPATIN(1)+1
      GO TO 18
   18 IMEMO(1,MV+2)=ME-IMEMO(1,MV+1)-IMEMO(1,MV+3)
                                                     RETURN
   19 IMEMO(1,MV+1) = IPATNE(1)
                                  IMEMO(1,MV+2) = 0
      IMEMO(1,MV+3) = IPATIN(1)
                                  RETURN
   20 IMEMO(1,MV+1) = 1-IPATNE(1)
                                    IMEMO(1,MV+3) = IPATIN(1)
      IMEMO(1,MV+2)=ME-IMEMO(1,MV+1)-IMEMO(1,MV+3)
                                                     RETURN
   21 IF (RNDM(-1.) -0.5) 22,22,23
   22 IMEMO(1,MV+1) = 1
                          GO TO 24
   23 IMEMO(1,MV+1) = 0
                          GO TO 24
   24 IF (MQ-1) 28,28, 25
   25 IF (RNDM(-1.) -0.5) 26,26,27
   26 IMEMO(1,MV+3) = 1
                          GO TO 28
   27 IMEMO(1,MV+3) = 0
                          GO TO 28
   28 LAMBDA = 2
   29 IF (MQ-1) 31,31,30
   30 IF (LAMBDA-3) 31,36,31
   31 TEMP2 = RNDM(-1.)
      MTEMP = MV+LAMBDA
      IF (TEMP2-1./3.) 33,32,32
   32 IF (TEMP2-2./3.) 34,35,35
   33 IMEMO(1,MTEMP) = 1
                           GO TO 36
   34 IMEMO(1,MTEMP) = 0
                           GO TO 36
   35 IMEMO(1,MTEMP) = -1
                           GO TO 36
   36 IF (LAMBDA-NP) 37,38,38
   37 LAMBDA = LAMBDA+1
                          GO TO 29
   38 SIGQ = 0.
      DO 39 I=1,NP
      ITEMP = MV+I
      SIGQ = SIGQ+IMEMO(1,ITEMP)
   39 CONTINUE
      IF (ME-SIGQ) 21,40,21
   40 RETURN
   41 IE1=1
      IF(RNDM(-1.).GT.0.5)  IE1=0
      BR=RNDM(-1.)
      IF(BR.LE.0.33333333)                      IE2=-1
      IF(BR.GT.0.33333333.AND.BR.LE.0.66666667) IE2= 0
      IF(BR.GT.0.66666667)                      IE2= 1
      IE3=ME-IE1-IE2
      IF(IE3.LT.0.OR.IE3.GT.1)  GO  TO  41
      IMEMO(1,MV+1)=IE1
      IMEMO(1,MV+2)=IE2
      IMEMO(1,MV+3)=IE3
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION  POTBAR(NU,X,Y,Z)
c
*  calculated as U=U0*rho(r)/rho(0) with U0=200 MeV according to the paper
*  by Ye.S.Golubeva, A.S.Iljinov at al. Nucl.Phys. A483 (1988) 539.
c
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      DATA PBAR0 /0.200/
      POTBAR=0.
      if(NU.eq.1) then
        A=A1
        C=C1
        D=D1
        AN=ANUCL1
        RMAX=RM1
      else
        A=A2
        C=C2
        D=D2
        AN=ANUCL2
        RMAX=RM2
      endif
      R=SQRT(X**2+Y**2+Z**2)
      if(R.gt.1.5*RMAX)  RETURN
      if(AN.le.10.) then
        RORO0=EXP(-(R/A)**2)
      else
        RORO0=(1.+EXP(-A/C)) / (1.+EXP((R-A)/C))
      endif
      POTBAR=PBAR0*RORO0
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION  POTENQ(P,IP,A,C,D,TF0,VPION,EPS)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Calculation of particle nuclear potential.
c
c
      DIMENSION   IP(5),P(9)
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      IF(IP(5).NE.0)  GO  TO  10
      IF(IP(3))  10,11,10
   10 POTENQ=0.
                     RETURN
   11 IF(IP(4))   13,12,13
   12 POTENQ = VPION
                     RETURN
   13 IF(ABS(A-A1)-.001) 100,101,101
  100 AN=ANUCL1
      RMAX=RM1
      NU=1
                     GO TO 105
  101 AN=ANUCL2
      RMAX=RM2
      NU=2
  105 CONTINUE
      R=SQRT(P(1)**2+P(2)**2+P(3)**2)/RMAX
c     R=RPOTEN(R)
      IF(R-1.5)  14,10,10
   14 R=R*RMAX
      IF(AN-10.) 106,106,107
  106 TF=TF0*EXP(-(2./3.)*(R/A)**2)
                               GO TO 108
  107 TF=TF0*(((1.+EXP(-A/C))/(1.+EXP((R-A)/C)))**0.6666667)
  108 CONTINUE
c     !!! 20.06.1995
      IF(NU.EQ.1) TF=TF*(AN1/ANUCL1)**(2./3.)
      IF(NU.EQ.2) TF=TF*(AN2/ANUCL2)**(2./3.)
c
      POTENQ=TF+EPS
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE RECUL(NU,PX,PY,PZ,X,Y,Z)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c   Change of linear and angular momenta of projectile (NU=1) 
c   or target (NU=2) by recul px,py,pz
c
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      IF(NU-1) 10,10,11
   10 PNUCL1(1)=PNUCL1(1)+PX
      PNUCL1(2)=PNUCL1(2)+PY
      PNUCL1(3)=PNUCL1(3)+PZ
      AMNUC1(1)=AMNUC1(1)+Z*PY-Y*PZ
      AMNUC1(2)=AMNUC1(2)+X*PZ-Z*PX
      AMNUC1(3)=AMNUC1(3)+Y*PX-X*PY
                                      GO TO 12
   11 PNUCL2(1)=PNUCL2(1)+PX
      PNUCL2(2)=PNUCL2(2)+PY
      PNUCL2(3)=PNUCL2(3)+PZ
      AMNUC2(1)=AMNUC2(1)+Z*PY-Y*PZ
      AMNUC2(2)=AMNUC2(2)+X*PZ-Z*PX
      AMNUC2(3)=AMNUC2(3)+Y*PX-X*PY
   12 CONTINUE
      RETURN
               END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c      DOUBLE PRECISION FUNCTION  RPOTEN(R)
c      REAL*8 R
c
c      RPOTEN=R
c      RETURN
c      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  COLDEV(R0X,R0Y,R0Z,T0,ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,
     *P01,R01,IPRO)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Changes of impact parameter and momentum of projectile 
c     acording to Coulomb trajectory
c
      DIMENSION  P01(3),R01(3)
      IF(ANUCL1-0.5)  7,7,8
    7 P0=SQRT(T0*(T0+0.28))
      E0=T0+0.14
      VZ=P0/(E0+ANUCL2*0.940)
      G=1./SQRT(1.-VZ*VZ)
      E0S=G*(E0-P0*VZ)
      P0S=SQRT(E0S**2-0.0196)
      ANEFF=0.148936
                       GO TO 9
    8 P0=ANUCL1*SQRT(T0*(T0+1.88))
      E0=ANUCL1*(T0+0.940)
      VZ=P0/(E0+ANUCL2*0.940)
      G=1./SQRT(1.-VZ*VZ)
      E0S=G*(E0-P0*VZ)
      P0S=SQRT(E0S**2-(ANUCL1*0.940)**2)
      ANEFF=ANUCL1
    9 B0=0.940*(ANEFF *ANUCL2/(ANEFF +ANUCL2))*ZNUCL1*ZNUCL2*
     *0.001440/P0S**2
      B=SQRT(R0X**2+R0Y**2)
      R12=SQRT(B*B+R0Z**2)
      G1=SQRT(B*B+B0*B0)
      SB=B0/G1
      CB=B/G1
      CF=R0X/B
      SF=R0Y/B
      PX=P0S*SB*CF
      PY=P0S*SB*SF
      PZ=G*(P0S*CB+B0*VZ)
      B1=B0+G1
      IF(B1-R12)  10,10,11
   11 IPRO=1
               RETURN
   10 IPRO=0
      R0X=B1*CF
      R0Y=B1*SF
      R0Z=-SQRT(R12**2-B1**2)
      PM=SQRT(PX**2+PY**2+PZ**2)
      BS=(PX*R0X+PY*R0Y+PZ*R0Z)/PM
      AS=SQRT(R12**2-BS**2)
      P01(1)=(R0Z-BS*PZ/PM)/AS
      P01(2)=(PX*R0Y-PY*R0X)/(AS*PM)
      P01(3)=PZ/PM
      R01(1)=(R0X-BS*PX/PM)/AS
      R01(2)=(PY*R0Z-PZ*R0Y)/(AS*PM)
      R01(3)=PX/PM
      RETURN
               END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE HELPQ(R0N,ANUCL,A,C,D,TF0,RM)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c    Determines nucleon density, Fermi momentum, Fermi energy
c
      DIMENSION   W(8),FIKS(8)
      DATA W/0.1012285363,0.2223810345,0.3137066459,0.3626837834,
     *0.3626837834,0.3137066459,0.2223810345,0.1012285363/,
     * FIKS/0.9602898565,0.7966664774,0.5255324099,0.1834346425,
     *-0.1834346425,-0.5255324099,-0.7966664774,-0.9602898565/
      PI=3.141593
      IF(ANUCL.le.10.) then
        A=R0N
        RM=A*SQRT(-LOG(D))
        RO0=ANUCL/((SQRT(PI)*A)**3)
      ELSE
        S=0.
        A=R0N*ANUCL**0.333333333
        RM=A+C*LOG((1.-D)/D)
        DO  K=1,8
          SK=(((FIKS(K)+1.)**2)/(EXP(RM*(FIKS(K)+1.)/(2.*C))+
     &    EXP(A/C)))*W(K)
          S=S+SK
        ENDDO
        S=(RM**3)*S/8.
        RO0=ANUCL/(12.566370*(EXP(A/C)+1.)*S)
      ENDIF
      TF0=0.1985*(RO0/2.)**0.666666667
      RH=A
      IF(ANUCL.LE.10.)  RH=A*SQRT(LOG(2.D0))
c      write(16,100) ANUCL,RH,RO0,TF0
c      write( *,100) ANUCL,RH,RO0,TF0
c  100 FORMAT(1X,'ANUCL = ',F6.0,'  R(1/2)=  ',F6.3,
c     &'  RO0(NUCL)=',F6.4,'  TF0=',F6.4)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION COSEXQ (L,T,CM)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C     BLOCK OF COSINUS CALCULATION FOR CHARGE EXCHANGE SCATTERING.
c  kkg 10/28/03        includes gamma +N --> pi+- + N1
c
      if(L.ne.0)              then
c  gamma + N => pi+-  + N1
        if(T.le.0.51)     then
c          COSEXQ = COSTAQ(24,T)
          COSEXQ = cosgaml(24)
        elseif(T.le.1.0)  then 
c          COSEXQ = COSTAQ(25,T)
          COSEXQ = cosgaml(25)
        else
c          COSEXQ = COSTAQ(26,T)
          COSEXQ = cosgaml(26)
        endif
        return
      endif 
   11 IF (T-0.08) 15,15,16
   15 COSEXQ = COSTAQ(12,T)
                            RETURN
   16 IF (T-0.3) 17,17,18
   17 COSEXQ = COSTAQ(13,T)
                            RETURN
   18 IF (T-1.0) 19,19,20
   19 COSEXQ = COSTAQ(10,T)
                            RETURN
   20 IF (T-2.4) 21,21,22
   21 COSEXQ = COSTAQ(11,T)
                            RETURN
   22 TMAX=3.5344*T*(T+CM)/(1.88*(T+CM)+CM*CM+0.8836)
      BET=BHN(T,CM)
      EXM=BET*TMAX
      IF(EXM.GT.30.)    EXM=30.
      COSEXQ=1.+(2.*LOG(1.+RNDM(-1.)*(EXP(-EXM)-1.)))/EXM
        RETURN
         END
c
c*********************************************************************
c
      double precision function cosgaml (j0)

c ======================================================================
c
c     Cosine calculation for elastic and 
c     charge-exchange gamma+ N reactions using linear interpolation
c   energy of gamma is fixed in inigam
c   Called by: COSELQ COSEXQ 
c
c   written by K. K. Gudima, Mar. 2004
c ======================================================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      common /ixsgpn/ thetai(181),cthetai(181),si(22,182),ri(22,181)
      data DEGRAD/0.0174532925199432958D0/
      if(j0.eq.22.or.j0.eq.23)  then
        jg=2
      elseif(j0.eq.24.or.j0.eq.25.or.j0.eq.26)  then
        jg=1
      else
        write(*,*)  ' cosgaml: j0=',j0
        stop
      endif 
      rr = RNDM(-1.)
      do  ir=1,181
        if(ABS(rr-ri(jg,ir)).lt.0.00001)  then
          cth=COS(thetai(ir)*DEGRAD)
          go  to  2
        endif
      enddo   
      do  ir=2,181
        if(rr.lt.ri(jg,ir))  then
          ir1=ir-1
          ir2=ir
          go  to  1
        endif
      enddo
      ir1=180
      ir2=181
    1 continue
      x1=ri(jg,ir1)  
      x2=ri(jg,ir2)  
      y1=thetai(ir1)
      y2=thetai(ir2)
      th=y1+(y2-y1)*((rr-x1)/(x2-x1))
      cth=COS(th*DEGRAD)
    2 temp1 = ABS(cth)
      if (temp1.le.1.0) then
        cosgaml= cth
      else
        cosgaml= SIGN(1.d0, cth)
      endif
      return

c ======================================================================
      end
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION BHN(T,CM)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c  Determines slope parameters of distribution exp(-b*t) for
c  high energy elastic scattering
c 
      IF(CM-0.9)  10,10,11
   10 BHN=11.
      IF(T.LT.10.)  BHN=7.5
      RETURN
   11 BHN=8.3+0.56*LOG(1.88*(T+1.88))
      RETURN
              END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION COSELQ (L,MQ,KSI,T,CM)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C     BLOCK OF COSINUS CALCULATION FOR ELASTIC SCATTERING.
c  kkg 10/28/03        includes gamma +N --> pi0 + N
c
   13 IF(MQ-2) 24,14,78
c  N-He4 scattering
   78 COSELQ=COSAL(T)
                       RETURN
c  N-N scattering
   14 GO  TO  (17,15,17,15,15,17),KSI
c   n + p:
   15 IF (T-0.97) 16,16,19
   16 COSELQ = COSTAQ(3,T)
                           RETURN
c   n + n or p + p:
   17 IF (T-0.46) 18,18,19
   18 COSELQ=1.-2.*RNDM(-1.)
                             RETURN
c   All types:
   19 IF (T-2.8) 20,20,21
   20 COSELQ = (1.+COSTAQ(1,T))/2.
                                   RETURN
   21 IF (T-10.) 22,22,23
   22 COSELQ = (3.+COSTAQ(2,T))/4.
                                   RETURN
   23 TMAX=3.5344*T*(T+CM)/(1.88*(T+CM)+CM*CM+0.8836)
      BET=BHN(T,CM)
      EXM=BET*TMAX
      IF(EXM.GT.30.)    EXM=30.
      COSELQ=1.+(2.*LOG(1.+RNDM(-1.)*(EXP(-EXM)-1.)))/EXM
      RETURN
   24 if(L.ne.0)  go  to  42
c  pi-N scattering
      IF (KSI-2) 25,34,33
c  pi+ p or pi- n scattering:
   25 IF (T-0.08) 26,26,27
   26 COSELQ = COSTAQ(4,T)
                           RETURN
   27 IF (T-0.3) 28,28,29
   28 COSELQ = COSTAQ(5,T)
                           RETURN
   29 IF (T-1.) 30,30,31
   30 COSELQ = COSTAQ(6,T)
                           RETURN
   31 IF (T-2.4) 32,32,23
   32 COSELQ = COSTAQ(7,T)
                           RETURN
c  pi0 p or pi0 n scattering:
   33 IF(RNDM(-1.)-0.5) 25,25,34
c  pi+ n or pi- p scattering:
   34 IF (T-0.08) 35,35,36
   35 COSELQ = COSTAQ(8,T)
                           RETURN
   36 IF (T-0.3) 37,37,38
   37 COSELQ = COSTAQ(9,T)
                           RETURN
   38 IF (T-1.0) 39,39,40
   39 COSELQ = COSTAQ(10,T)
                            RETURN
   40 IF (T-2.4) 41,41,23
   41 COSELQ = COSTAQ(11,T)
                            RETURN
   42 continue
c  gamma + N => pi0 + N
      if(T.le.0.45)  then
c        COSELQ = COSTAQ (22,T)
        COSELQ = cosgaml(22)
      else
c        COSELQ = COSTAQ (23,T)
        COSELQ = cosgaml(23)
      endif
      RETURN
         END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DHIST1(X,A,B,H,RX,N,W)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C     BLOCK OF BILDING OF HISTOGRAMMS.
c
      DIMENSION RX(N)
      RX(N) = RX(N)+X*W
      IF (X-A) 1,2,2
    1 RX(N-4) = RX(N-4)+W
                            RETURN
    2 IF (X-B) 4,3,3
    3 RX(N-2) = RX(N-2)+W
                            RETURN
    4 L=(X-A)/H
                  NL=N-5
                           IF(L.GT.NL) GO TO 5
      RX(L+1)=RX(L+1)+W
      RX(N-1)=RX(N-1)+X*W
                            RX(N-3)=RX(N-3)+W
                                                RETURN
    5 write(16,6)  X,A,B,H,W
    6 FORMAT(25X,'MISTAKE IN DIMENSION OF DHIST '/25X,5(F10.3,5X))
      RETURN
      END
      SUBROUTINE DHIST2(X,A,B,H,RX,N,W,J)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C       BLOCK OF BUILDING 2-DIMENSIONAL HISTOGRAMS
c
      DIMENSION  RX(N,2)
      RX(N,J) = RX(N,J)+X
      IF(X-A) 1,2,2
    1 RX(N-4,J) = RX(N-4,J)+W
                                  RETURN
    2 IF(X-B) 4,3,3
    3 RX(N-2,J) = RX(N-2,J)+W
                                  RETURN
    4 L = (X-A)/H
                      RX(L+1,J) = RX(L+1,J)+W
      RX(N-1,J) = RX(N-1,J)+X
                                  RX(N-3,J) = RX(N-3,J)+W
                                                              RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DHIST3(X,A,B,H,RX,N,M,W,J)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C       BLOCK OF BUILDING M-DIMENSIONAL HISTOGRAMS
c
      DIMENSION  RX(N,M)
      RX(N,J) = RX(N,J)+X*W
      IF(X-A) 1,2,2
    1 RX(N-4,J) = RX(N-4,J)+W
                                  RETURN
    2 IF(X-B) 4,3,3
    3 RX(N-2,J) = RX(N-2,J)+W
                                  RETURN
    4 L=(X-A)/H
                  NL=N-5
                           IF(L.GT.NL) GO TO 5
      RX(L+1,J)=RX(L+1,J)+W
      RX(N-1,J)=RX(N-1,J)+X*W
                                  RX(N-3,J) = RX(N-3,J)+W
                                                              RETURN
    5 write(16,6)  X,A,B,H,J
    6 FORMAT(25X,'MISTAKE IN DIMENSION OF DHIST3'/25X,5(F10.3,5X),I3)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DHIST4(X,A,B,H,RX,N,M,L,W,I,J)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C       BLOCK OF BUILDING (N,M)-DIMENSIONAL HISTOGRAMS
c
      DIMENSION RX(N,M,L)
      RX(I,J,L)=RX(I,J,L)+X*W
      IF(X-A) 1,2,2
    1 RX(I,J,L-4)=RX(I,J,L-4)+W
                                  RETURN
    2 IF(X-B) 4,3,3
    3 RX(I,J,L-2)=RX(I,J,L-2)+W
                                  RETURN
    4 L1=(X-A)/H
                   NL=L-5
                            IF(L1.GT.NL) GO TO 5
      RX(I,J,L1+1) = RX(I,J,L1+1)+W
      RX(I,J,L-1)=RX(I,J,L-1)+X*W
                                  RX(I,J,L-3)=RX(I,J,L-3)+W
      RETURN
    5 write(16,6)  X,A,B,H,W
    6 FORMAT(25X,'MISTAKE IN DIMENSION OF DHIST4'/25X,5(F10.3,5X))
      RETURN
               END
C     * * * * * * * * * * * * * * * * *
      BLOCK DATA SIGAR
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     This subprogram puts the tables of cross sections
c     and energies into the appropriate common blocks.
c     EDITED by KKG, October, 2003, April 2007 
c
      COMMON /TABELE/SIGMA,ARGUS/TYPECSQ/ICST,NSICST
      DIMENSION SIGMA(30,46) , ARGUS(30,9) , ICST(46) , NSICST(37)
     *,A1(30,3),A2(30,3),A3(30,3),A4(30,3),
     * A5(30,3),A6(30,3),A7(30,3),A8(30,3),A9(30,2),
     * B1(30,2),B2(30,2),B3(30,2),B4(30,2),EX(30,3),
     *AR1(30,3),AR2(30,2),AR3(30,1),
     &ARG(30,3),AG1(30,2),AG2(30,3),AG3(30,3),AG4(30,1)
      EQUIVALENCE (SIGMA(1, 1),A1(1,1)),(SIGMA(1, 4),A2(1,1)),
     *            (SIGMA(1, 7),A3(1,1)),(SIGMA(1,10),A4(1,1)),
     *            (SIGMA(1,13),A5(1,1)),(SIGMA(1,16),A6(1,1)),
     *            (SIGMA(1,19),A7(1,1)),(SIGMA(1,22),A8(1,1)),
     *            (SIGMA(1,25),A9(1,1)),(SIGMA(1,27),B1(1,1)),
     *            (SIGMA(1,29),B2(1,1)),(SIGMA(1,31),B3(1,1)),
     *            (SIGMA(1,33),B4(1,1)),(SIGMA(1,35),EX(1,1)),
     *            (ARGUS(1, 1),AR1(1,1)),(ARGUS(1,4),AR2(1,1)),
     *            (ARGUS(1, 6),AR3(1,1)),
     &            (ARGUS(1,7),ARG(1,1)),(SIGMA(1,38),AG1(1,1)),
     &            (SIGMA(1,40),AG2(1,1)),(SIGMA(1,43),AG3(1,1)),
     &            (SIGMA(1,46),AG4(1,1))
      DATA  A1/
c   (1) N+N , P+P total
     &        17613.0,   302.911,   150.083,    95.861,    69.160, 
     &         54.056,    38.395,    29.639,    23.563,    22.800, 
     &         22.600,    23.100,    24.045,    25.957,    30.747,
     &         42.160,    48.114,    47.430,    47.465,    47.307,  
     &         46.681,    45.082,    42.476,    41.236,    40.866, 
     &         40.143,    40.184,    38.352,    38.375,    38.406,
c   (2) N+N,  P+P elastic
     &        17613.0,   302.911,   150.083,    95.861,    69.160, 
     &         54.056,    38.395,    29.639,    23.563,    22.800, 
     &         22.600,    23.097,    23.399,    23.824,    24.623, 
     &         25.060,    24.306,    23.388,    22.543,    21.101, 
     &         19.701,    17.264,    14.477,    12.848,    11.745, 
     &         10.691,     9.640,     8.583,     8.296,     7.956, 
c   (3) N+P  total
     &        20357.0,   912.566,   499.490,   288.179,   208.121, 
     &        161.621,   106.227,    71.137,    50.588,    42.483, 
     &         37.982,    35.381,    33.820,    33.299,    34.219, 
     &         37.566,    38.801,    39.255,    40.368,    41.570, 
     &         41.624,    42.144,    42.584,    42.091,    41.847, 
     &         41.26029,  39.48050,  39.40114,  39.37191,  39.33528/
      DATA  A2/
c   (4) N+P  elastic
     &        20357.0,   912.566,   499.490,   288.179,   208.121, 
     &        161.621,   106.227,    71.137,    50.588,    42.483, 
     &         37.982,    35.375,    33.000,    30.700,    27.800, 
     &         24.700,    21.100,    20.000,    18.750,    17.700, 
     &         17.000,    15.750,    14.200,    13.300,    12.700, 
     &         11.950,    11.100,     9.97127,   9.60346,   9.16243, 
c   (5) pi+ + N, pi- + P total
     &          5.200,     6.400,     9.300,    17.000,    37.150, 
     &         58.940,    72.000,    59.050,    39.880,    26.460, 
     &         27.430,    32.010,    46.730,    35.900,    45.000, 
     &         59.550,    43.570,    36.440,    36.430,    36.530, 
     &         34.780,    34.370,    36.020,    36.300,    34.600, 
     &         32.950,    32.200,    28.900,    26.000,    25.47546, 
c   (6) pi+ + N, pi- + P elastic
     &          1.593,     1.784,     2.281,     4.117,    10.467, 
     &         19.657,    23.245,    23.589,    15.349,     9.931, 
     &         10.087,    14.277,    20.290,    14.807,    18.139, 
     &         25.873,    18.879,    14.746,    13.019,    12.345, 
     &          9.897,    10.019,     9.553,     8.942,     8.448, 
     &          8.034,     7.188,     5.325,     4.258,     3.991/
      DATA  A3/
c   (7) pi(+0-) + (N,P)  single charge exchange
     &          3.607,     5.303,     8.281,    13.170,    26.134, 
     &         42.727,    48.927,    38.530,    25.550,    14.713, 
     &         10.465,     8.442,     7.448,     4.610,     6.693, 
     &          6.584,     2.828,     2.465,     2.438,     2.360, 
     &          2.026,     1.459,     1.052,     0.792,     0.610, 
     &          0.380,     0.255,     0.0903,    0.0381,    0.0252, 
c   (8) pi- + N,  pi+ + P total
     &          2.015,     5.866,    12.637,    32.909,    97.023,
     &        172.340,   209.160,   167.843,   103.625,    52.363,
     &         30.821,    18.862,    14.257,    15.979,    22.619, 
     &         26.113,    26.892,    31.413,    37.517,    41.354,
     &         37.945,    31.563,    29.571,    29.083,    30.211,
     &         31.127,    28.865,    26.390,    23.357,    23.348,
c   (9) pi- + N,  pi+ + P elastic
     &          2.015,     5.866,    12.637,    32.909,    97.023, 
     &        172.340,   209.160,   167.843,   103.625,    52.363, 
     &         30.566,    17.281,    11.473,     8.214,     9.272, 
     &         12.004,    13.186,    14.052,    15.522,    18.586, 
     &         15.715,    12.785,    10.808,     9.548,     8.649,
     &          7.936,     7.051,     5.668,     4.216,     3.951/  
      DATA  A4/
c   (10) pion + np absorption cross section
c   (Must be multiplied by a function [currently = 4.0] for use)
     &          3.900,     4.350,     6.090,     8.540,    11.050, 
     &         11.750,    10.190,     6.540,     3.400,     1.700, 
     &          0.825,     0.380,     0.200,     0.130,     0.089, 
     &          0.063,     0.049,     0.040,     0.031,     0.026, 
     &          0.018,     0.0125,    0.0076,    0.0040,    0.0021, 
     &          0.00085,   0.00030,   0.00002,   0.000003,  0.0000003, 
c   (11) p + p --> p + p + pi0
     &          0.000,     0.000,     0.00295,   0.01969,   0.07956, 
     &          0.24419,   0.5399,    1.086,     1.953,     3.029, 
     &          3.522,     3.730,     3.845,     4.181,     4.526, 
     &          4.619,     4.571,     4.483,     4.282,     4.164, 
     &          4.099,     4.023,     3.960,     3.800,     3.374, 
     &          3.089,     2.603,     2.889,     1.532,     1.320, 
c   (12) p + p --> p + n + pi+
     &          0.000,     0.000,     0.000,     0.35747,   0.79723, 
     &          1.582,     2.053,     4.773,     7.542,    10.113, 
     &         13.214,    15.564,    16.767,    17.219,    17.589, 
     &         17.851,    18.190,    18.279,    18.330,    17.926, 
     &         17.546,    17.068,    16.420,    16.067,    13.130, 
     &         10.534,     8.688,     7.278,     4.480,     3.686/
      DATA  A5/
c  (13) p + n --> p + n + pi0
     &          0.000,     0.000,     0.05749,   0.28462,   0.5271, 
     &          0.9429,    1.8173,    3.142,     4.761,     6.508, 
     &          6.778,     7.169,     7.400,     7.550,     7.700, 
     &          7.806,     7.982,     8.195,     8.424,     8.640, 
     &          8.920,     9.201,     9.055,     8.600,     7.700, 
     &          6.550,     5.129,     3.985,     2.448,     1.300, 
c  (14) p + n --> p + p + pi-
     &          0.000,     0.000,     0.000,     0.000,     0.06818,
     &          0.15784,   0.29075,   0.6300,    1.0910,    1.3824, 
     &          1.6719,    1.9102,    2.1401,    2.350,     2.526, 
     &          2.675,     2.806,     2.918,     3.009,     3.097, 
     &          3.168,     3.227,     3.199,     3.121,     2.875, 
     &          2.51422,   1.85326,   1.46821,   1.22737,   1.07373, 
c  (15) pi+ + p --> pi+ + p + pi0
     &          0.00227,   0.02183,   0.08394,   0.22867,   0.38106, 
     &          0.63783,   1.0998,    2.0452,    3.1995,    4.574, 
     &          6.276,     9.381,    10.310,    10.343,     9.484, 
     &          8.784,     8.277,    10.559,    12.108,    11.563, 
     &         10.076,     7.845,     6.545,     5.452,     3.589, 
     &          2.6616,    1.77222,   1.15771,   0.79248,   0.60762/ 
      DATA  A6/
c  (16) pi+ + p --> pi+ + n + pi+
     &          0.00147,   0.03056,   0.05862,   0.09654,   0.13958, 
     &          0.21361,   0.34633,   0.51922,   0.68048,   0.93144, 
     &          1.1971,    1.480,     1.781,     1.977,     2.083, 
     &          2.097,     2.117,     2.499,     3.107,     3.474, 
     &          3.648,     3.631,     3.099,     2.309,     2.467, 
     &          2.102,     1.2741,    0.85565,   0.70307,   0.53185, 
c  (17) pi- + p --> pi0 + n + pi0
     &          0.01433,   0.21686,   0.61527,   1.25405,   1.52076, 
     &          1.7479,    1.9039,    2.0563,    2.1526,    2.250, 
     &          2.402,     2.807,     3.294,     3.623,     3.085, 
     &          2.251,     1.841,     1.445,     1.183,     0.9801, 
     &          0.80905,   0.59528,   0.62035,   0.70180,   0.40717, 
     &          0.25280,   0.12009,   0.08365,   0.015573,  0.013499, 
c  (18) pi- + p --> pi- + p + pi0
     &          0.00193,   0.02508,   0.0919,    0.2654,    0.6346, 
     &          1.180,     2.020,     3.237,     5.008,     4.861,
     &          4.529,     5.011,     5.687,     6.188,     6.426, 
     &          5.896,     5.286,     4.383,     4.632,     4.700, 
     &          4.737,     5.115,     4.128,     3.363,     3.395, 
     &          1.876,     2.0066,    1.1843,    0.68066,   0.60043/ 
      DATA  A7/
c  (19) pi- + p --> pi- + n + pi+
     &          0.00908,   0.12116,   0.64364,   1.5597,    3.153, 
     &          3.971,     5.004,     5.860,     6.306,     6.474, 
     &          7.194,     8.300,     9.978,    11.775,    11.789, 
     &          9.938,     7.522,     7.512,     7.930,     7.781, 
     &          7.089,     6.944,     6.127,     5.332,     4.000, 
     &          3.385,     2.967,     2.10879,   1.20557,   0.87941, 
     *   310.00,   270.00,   240.00,   213.00,   172.00,
     *   156.00,   132.00,   116.00,   108.00,   102.00,
     *    98.00,    96.00,    96.00,    98.00,   102.00,
     *   109.00,   116.00,   122.00,   128.00,   132.00,
     *   136.00,   138.00,   140.00,   140.00,   140.00,
     *   141.00,   140.00,   138.00,   134.00,   132.00,
     *   204.00,   165.00,   138.00,   117.00,    87.00,
     *    75.00,    55.00,    40.00,    30.00,    24.00,
     *    22.00,    20.00,    20.00,    20.00,    21.00,
     *    23.00,    26.00,    28.00,    31.00,    33.00,
     *    34.00,    36.00,    36.00,    36.00,    36.00,
     *    36.00,     36.00,     36.00,     34.00,     32.00/
      DATA  A8/
     *   266.00,   228.00,   200.00,   178.00,   148.00,
     *   136.00,   115.00,   102.00,    94.00,    86.00,
     *    79.00,    76.00,    76.00,    78.00,    81.00,
     *    85.00,    91.00,   100.00,   107.00,   112.00,
     *   115.00,   118.00,   118.00,   118.00,   118.00,
     *   118.00,   116.00,   112.00,   110.00,   106.00,
     *   133.00,   105.00,    84.00,    68.00,    46.00,
     *    38.00,    28.00,    22.00,   19.00,    15.00,
     *    11.00,   10.00,     9.50,    10.00,   11.00,
     *    12.00,    15.00,    17.00,    19.00,    21.00,
     *    23.00,    24.00,    24.00,    24.00,    23.00,
     *    22.00,    22.00,    21.00,   20.00,    20.00,
     *     0.00,     0.00,     0.00,     0.20,     0.80,
     *     1.60,     3.20,     6.20,    33.00,    40.30,
     *    41.60,    42.00,    42.00,    42.00,    41.90,
     *    41.50,    40.70,    40.00,    39.30,    38.40,
     *    37.60,    36.60,    35.40,    34.20,    32.80,
     *    31.20,     29.40,     27.70,     26.20,     16.00/
      DATA  A9/
     *000000.00,     0.00,     0.00,     0.30,     1.20,
     *     2.20,     3.40,     4.90,     9.00,    11.20,
     *    13.40,    15.60,    18.10,    20.30,    22.30,
     *    25.40,    27.10,    27.90,    28.30,    28.60,
     *    28.80,    29.00,    29.10,    29.20,    29.30,
     *    29.30,    29.30,    29.30,    29.30,    29.30,
     *     0.00,     1.00,    30.00,    41.00,    43.20,
     *    43.60,    43.40,    43.00,    41.00,    39.00,
     *    32.00,    23.50,    16.00,    12.00,     9.40,
     *     5.60,     3.80,     2.90,     2.20,    1.80,
     *     1.60,     1.30,     1.00,     0.80,     0.70,
     *     0.50,      0.40,      0.40,      0.30,      0.00/
      DATA  B1/
     * 10.0, 11.0, 11.5, 11.9, 12.0, 12.1, 12.2, 12.2, 12.8, 14.5,
     * 16.3, 18.9, 18.0, 17.6, 17.2, 17.1, 17.0, 17.0, 17.1, 17.2,
     * 17.5, 17.8, 18.0, 19.0, 19.4, 20.0, 20.7, 21.4, 22.0, 24.0,
     * 10.0, 11.0, 11.5, 11.9, 12.0, 12.1, 12.2, 12.1, 12.0, 11.9,
     * 11.6, 10.5,  9.0,  6.5,  5.0,  4.1,  3.5,  3.2,  3.2,  3.1,
     *  2.9,  2.5,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4/
      DATA  B2/
     *  6.0,  7.0,  7.5, 10.0, 12.0, 13.0, 13.0, 13.5, 15.0, 16.5,
     * 18.5, 20.5, 19.0, 18.5, 17.8, 17.3, 17.2, 17.2, 17.2, 17.5,
     * 17.7, 18.0, 18.2, 18.5, 19.0, 19.7, 20.0, 20.3, 20.6, 21.0,
     *  3.0,  4.0,  4.5,  6.0,  8.0,  5.5,  5.5,  5.7,  5.7,  5.7,
     *  5.7,  5.5,  5.2,  5.0,  4.8,  4.5,  4.2,  3.8,  3.5,  3.1,
     *  2.3,  2.5,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4/
      DATA  B3/
     * 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 20.0, 26.0, 29.0, 31.0,
     * 31.0, 29.0, 28.0, 23.0, 22.0, 21.0, 20.7, 20.5, 20.3, 20.2,
     * 20.0, 20.0, 20.0, 20.0, 20.1, 20.5, 21.0, 22.0, 23.0, 24.0,
     *  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  7.0,  9.0, 15.0, 18.0,
     * 17.0, 10.0,  7.5,  5.0,  3.4,  3.1,  3.0,  2.9,  2.8,  2.7,
     * 2.65, 2.60, 2.50, 2.51, 2.52, 2.55, 2.56, 2.57, 2.58, 2.58/
      DATA  B4/
     *450.0,323.0,125.0, 82.0, 65.0, 38.0,30.75, 29.0, 32.5, 34.0,
     * 50.0, 32.0, 33.0, 30.0, 28.0, 27.0, 25.0, 24.0, 23.0, 22.0,
     * 21.0, 20.5, 20.0, 20.0, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5,
     *150.0, 98.0, 60.0, 42.0, 33.0, 23.0, 18.0, 16.0, 17.0, 20.0,
     * 22.0, 15.0,  9.0,  8.0,  5.5,  4.5,  4.0,  3.6,  3.2,  2.9,
     *  2.8,  2.6,  2.5, 2.51, 2.52, 2.55, 2.56, 2.57, 2.58, 2.58/
      DATA  EX/
     *  3.0,  3.0,  3.0,  4.0,  4.0,  7.5,  7.5,  7.8,  7.8,  7.8,
     *  7.8,  4.5,  3.0,  1.8,  0.8, 0.45, 0.30, 0.15,0.075,0.035,
     *  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     * 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,  9.0,  8.0,  7.0,
     *  6.5,  4.0, 2.75, 1.50, 0.65, 0.30,  0.0,  0.0,  0.0,  0.0,
     *  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     *300.0,225.0, 65.0, 40.0, 32.0, 15.0,12.75, 12.0, 15.5, 14.0,
     * 9.75,  5.0,  4.0,  2.0, 1.25,  0.5,  0.0,  0.0,  0.0,  0.0,
     *  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
c
c
      DATA AG1/
c   (38) gamma + p --> p + pi0 
c   New constants (August, 1998) (Valid to 10 GeV [use log for low E])
     &          0.00021,   0.00437,   0.02075,   0.09167,   0.2423,
     &          0.2796,    0.1371,    0.08599,   0.05499,   0.04352,
     &          0.03260,   0.03310,   0.03990,   0.04336,   0.03950,
     &          0.03307,   0.02750,   0.02847,   0.02950,   0.02337,
     &          0.01960,   0.01637,   0.01320,   0.00863,   0.00574,
     &          0.00346,   0.00200,   0.00076,   0.000212,  0.000065,
c   (39) gamma + p --> n + pi+:
c   New constants (August, 1998) (Valid to 10 GeV [use log for low E])
     &          0.000001,  0.08169,   0.10202,   0.15975,   0.2343,
     &          0.2031,    0.1230,    0.10029,   0.09024,   0.08760,
     &          0.08789,   0.09241,   0.10347,   0.09100,   0.06329,
     &          0.04951,   0.04944,   0.05184,   0.05612,   0.05019,
     &          0.03622,   0.02501,   0.01922,   0.01350,   0.00774,
     &          0.00516,   0.00250,   0.00075,   0.000232,  0.000083/
      data AG2/
c   (40) gamma + [NN] absorption:
c   New constants (August, 1998) (Valid to 5 GeV)
c   (Must be multiplied by a function [currently = 5.0] for use)
     &          0.000,     1.14301,   1.8797,    2.5348,    1.3234,
     &          0.57795,   0.3391,    0.2289,    0.1259,    0.08875,
     &          0.06946,   0.06128,   0.05130,   0.05410,   0.05460,
     &          0.05805,   0.06317,   0.05787,   0.04589,   0.03161,
     &          0.02021,   0.01355,   0.01090,   0.00830,   0.00600,
     &          0.00450,   0.00340,   0.00270,   0.00100,   0.00050,
c   (41) gamma + p --> pi- + p + pi+
c   New constants (August, 1998) (Valid to 20 GeV [use log for low E])
     &          1.82e-6,   0.000163,  0.00312,   0.00762,   0.02099,
     &          0.04421,   0.06208,   0.07277,   0.07452,   0.07485,
     &          0.07468,   0.07448,   0.07434,   0.07355,   0.07194,
     &          0.07058,   0.06931,   0.06781,   0.06471,   0.06091,
     &          0.05760,   0.05302,   0.04648,   0.03943,   0.03271,
     &          0.02806,   0.02357,   0.01731,   0.01368,   0.01219,
c   (42) gamma + p --> pi0 + p + pi0
c   New constants (August, 1998) (Valid to 20 GeV [use log for low E])
c   These coefficients are found by multiplying those of # 41 by
c   0.142, found by optimally fitting the limited available data
c   in the energy range 464 to 779 MeV gamma lab energy to the #23
c   cross section linearly scaled.
     &          2.58e-7,   0.0000231, 0.000443,  0.00108,   0.00298,
     &          0.00628,   0.00882,   0.01033,   0.01058,   0.01063,
     &          0.01060,   0.01058,   0.01056,   0.01044,   0.01022,
     &          0.01022,   0.00984,   0.00963,   0.00919,   0.00865,
     &          0.00818,   0.00753,   0.00660,   0.00560,   0.00464,
     &          0.00398,   0.00335,   0.00246,   0.00194,   0.00173/
      data AG3/
c   (43) gamma + p --> pi+ + n + pi0
c   New constants (August, 1998) (Valid to 20 GeV [use log for low E])
c   These coefficients are found by fitting the 12 data points from
c   464 to 779 MeV gamma lab energy, and fixing the coefficients
c   from 800 to 20000 MeV to 2/3 * (#41 coefficients).  There is
c   an unphysical bump (with negligible cross section) at energy below 
c   400 MeV due to limitations in the parabolic interpolation method.
     &          2.1e-9,    1.1e-6,    0.00061,   0.00124,   0.00610,
     &          0.01491,   0.02359,   0.03218,   0.04855,   0.05341,
     &          0.04979,   0.04965,   0.04956,   0.04903,   0.04796,
     &          0.04705,   0.04621,   0.04520,   0.04314,   0.04061,
     &          0.03849,   0.03532,   0.03099,   0.02629,   0.02181,
     &          0.01871,   0.01570,   0.011557,  0.009120,  0.008126,
c   (44) gamma + p --> delta + pi
c   New constants (August, 1998) (Valid to 10 GeV [use log for low E])
c   There is an unphysical zero at energy below 400 MeV due to 
c   limitations in the parabolic interpolation method (the cross
c   section is negligible anyway).
     &          1.8e-11,   0.00002,   0.00082,   0.00762,   0.02099,
     &          0.04421,   0.06208,   0.06935,   0.06387,   0.05537,
     &          0.04570,   0.04594,   0.05130,   0.05416,   0.05146,
     &          0.04503,   0.03727,   0.02564,   0.02613,   0.01998,
     &          0.01725,   0.01511,   0.01071,   0.00634,   0.00414,
     &          0.00289,   0.00174,   0.000821,  0.000257,  0.000124,
c   (45)  gamma + p --> 2 pi + N (#41 + #42 + #43)
c   New constants (August, 1998) (Valid to 20 GeV [use log for low E])
     &          2.1e-6,    0.000179,  0.00417,   0.00994,   0.03007,
     &          0.06540,   0.09449,   0.11528,   0.13365,   0.13889,
     &          0.13507,   0.13471,   0.13446,   0.13302,   0.13012,
     &          0.12785,   0.12536,   0.12264,   0.11704,   0.11017,
     &          0.10427,   0.09587,   0.08407,   0.07132,   0.05916,
     &          0.05075,   0.04262,   0.03133,   0.02474,   0.02205/
      data AG4/
c   (46)  gamma + p --> X (TOTAL)
c   New constants (August, 1998) (Valid to 100 GeV [use log for low E])
     &          0.00010,   0.12102,   0.24837,   0.4994,    0.4709, 
     &          0.2727,    0.1877,    0.1699,    0.1880,    0.2146,
     &          0.2339,    0.2659,    0.2772,    0.2316,    0.2093,
     &          0.2102,    0.2013,    0.2142,    0.2120,    0.1842,
     &          0.1693,    0.1539,    0.1509,    0.1524,    0.1446,
     &          0.1394,    0.1283,    0.1289,    0.11739,   0.11714/
      DATA  ICST/
     *  210 ,   211 ,   220 ,   221 ,   120 ,   121 ,   122 ,
     *  110 ,   111 ,   123 ,   214 ,   215 ,   224 ,   225 ,
     *  114 ,   115 ,   126 ,   124 ,   125 ,
     *00510 , 00511 , 00410 , 00411 , 00514 , 00515 , 00516,
     * 1110 ,  1111 ,  1120 ,  1121 , -1110 , -1111 , -1120,
     *-1121 ,  1122 , -1112 , -1122 ,
     & 10111,  10112,  10113,  10115,  10114,  10116,  10118,
     & 10117,  10110/
      DATA  NSICST/
     *  112 ,   113 ,   116 ,   117 ,   127 ,   130 ,   131 ,
     *  132 ,   133 ,   134 ,   135 ,   136 ,   137 ,   212 ,
     *  213 ,   216 ,   217 ,   222 ,   223 ,   226 ,   227 ,
     *  230 ,   231 ,   232 ,   233 ,
     *  240 ,   241 ,   242 ,   243 ,
     *  250 ,   251 ,   252 ,   253 ,
     *  260 ,   261 ,   262 ,   263 /
      DATA  AR1/
c   Lab kinetic energies for nucleon-nucleon reactions (#1-#4):
     &          0.000,     0.010,     0.020,     0.030,     0.040, 
     &          0.050,     0.070,     0.100,     0.150,     0.200, 
     &          0.250,     0.300,     0.350,     0.400,     0.500, 
     &          0.650,     0.850,     0.950,     1.100,     1.300, 
     &          1.500,     2.000,     3.000,     4.000,     5.000, 
     &          7.000,    10.000,    16.000,    22.000,    30.000, 
c  Lab kinetic energies for pion induced reactions(2)  (#5-#10):
     &          0.000,     0.030,     0.050,     0.080,     0.120, 
     &          0.150,     0.175,     0.210,     0.250,     0.320, 
     &          0.400,     0.500,     0.600,     0.700,     0.800, 
     &          0.900,     1.000,     1.100,     1.200,     1.300, 
     &          1.450,     1.600,     1.800,     2.000,     2.250, 
     &          2.500,     3.000,     5.000,    12.000,    20.000, 
c   Lab kinetic energies for pion production reactions (3) (#11-19):
     &          0.200,     0.250,     0.300,     0.350,     0.400, 
     &          0.450,     0.500,     0.550,     0.600,     0.650, 
     &          0.700,     0.750,     0.800,     0.850,     0.900, 
     &          0.950,     1.000,     1.100,     1.200,     1.300, 
     &          1.400,     1.600,     1.800,     2.000,     2.500, 
     &          3.000,     4.000,     5.000,     7.000,    10.000/
      DATA  AR2/
c   Lab kinetic energies for N + He-4 and N + He-3(H-3) reactions
     *   0.05  ,   0.06  ,   0.07  ,   0.08  ,   0.10  ,
     *   0.11  ,   0.13  ,   0.15  ,   0.17  ,   0.20  ,
     *   0.25  ,   0.30  ,   0.35  ,   0.40  ,   0.45  ,
     *   0.50  ,   0.55  ,   0.60  ,   0.65  ,   0.70  ,
     *   0.75  ,   0.80  ,   0.85  ,   0.90  ,    1.00  ,
     *   2.00  ,   3.00  ,   5.00  ,   10.0  ,   30.0  ,
c   Lab kinetic energies for N + He-4 ==> N + N + (A=3)
     *   0.0185,   0.0200,   0.0224,   0.0250,   0.0282,
     *   0.0316,   0.0355,   0.0400,   0.0500,   0.0560,
     *   0.0630,   0.0710,   0.0800,   0.0900,   0.1000,
     *   0.1260,   0.1570,   0.2000,   0.2500,   0.3160,
     *   0.4000,   0.5000,   0.6300,   0.8000,   1.0000,
     *   1.2600,   1.5700,   2.0000,   2.5000,   30.000/
      DATA  AR3/
c   Lab kinetic energies for K + N reactions
     *   0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
     *   1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0,15.0,
     *  20.0,30.0, 50.0, 100.0, 150.0, 200.0, 300.0, 400.0,
     * 500.0, 1000.0/
      data arg/
c   Lab Energies for gamma-induced pion production cross sections 
c   (new version, #'s 38,39,46)
     &          0.150,     0.165,     0.200,     0.250,     0.300,
     &          0.350,     0.400,     0.450,     0.500,     0.550,
     &          0.600,     0.650,     0.700,     0.750,     0.800,
     &          0.850,     0.900,     0.950,     1.000,     1.050,
     &          1.100,     1.150,     1.200,     1.350,     1.500,
     &          1.750,     2.250,     3.500,     6.000,    10.000, 
c   Lab energies for gamma absorption on 2 nucleons (5) (#40):
     &          0.002226,  0.0025,    0.003,     0.005,     0.010, 
     &          0.020,     0.030,     0.040,     0.060,     0.080, 
     &          0.100,     0.120,     0.150,     0.170,     0.200, 
     &          0.220,     0.250,     0.300,     0.330,     0.360, 
     &          0.400,     0.440,     0.500,     0.600,     0.700, 
     &          0.800,     0.900,     1.000,     2.000,     3.000, 
c   Lab energies for gamma induced 2-pion production (6) (#41-45):
     &          0.320,     0.350,     0.400,     0.450,     0.500, 
     &          0.550,     0.600,     0.650,     0.700,     0.750, 
     &          0.800,     0.850,     0.900,     0.950,     1.000, 
     &          1.050,     1.100,     1.150,     1.200,     1.250, 
     &          1.300,     1.500,     1.700,     2.000,     2.500, 
     &          3.000,     4.000,     6.000,    10.000,    14.000/
c
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      BLOCK DATA COEFABQ
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C     SUBROUTINE WHICH PUT THE COEFFICIENTS ANKJ,BNKJ AND CKJ
C     IN THE MAIN PROGRAM
c  kkg  10/28/03
      COMMON /COEFAQ/ ANKJ  /COEFBCQ/ BNKJ,CKJ
      DIMENSION ANKJ(4,4,29),BNKJ(4,4,8),CKJ(3,8)
     *,A1(4,4,4),A2(4,4,4),A3(4,4,4),A4(4,4,4),
     * A5(4,4,4),A6(4,4,1),A7(4,4,7),A8(4,4,1),
     * B1(4,4,4),B2(4,4,4)
      EQUIVALENCE (ANKJ(1,1, 1),A1(1,1,1)),
     *            (ANKJ(1,1, 5),A2(1,1,1)),
     *            (ANKJ(1,1, 9),A3(1,1,1)),
     *            (ANKJ(1,1,13),A4(1,1,1)),
     *            (ANKJ(1,1,17),A5(1,1,1)),
     *            (ANKJ(1,1,21),A6(1,1,1)),
     *            (ANKJ(1,1,22),A7(1,1,1)),      
     *            (ANKJ(1,1,29),A8(1,1,1)),      
     *            (BNKJ(1,1, 1),B1(1,1,1)),
     *            (BNKJ(1,1, 5),B2(1,1,1))
      DATA  A1/
c  ankj Angular distribution coefficients:
c  j = 1;  N + N elastic scattering; Tlab <= 2.8 GeV:
c          (n + n & p + p isotropic below 0.46 GeV.)
     * 2.7404E 00 , -9.6998E 00 ,  1.0400E 01 ,  2.3882E 00 ,
     *-7.5137E 00 ,  4.4096E 01 , -7.4379E 01 ,  4.6038E 01 ,
     * 7.5479E 00 , -3.9274E 01 ,  6.4835E 01 , -4.1609E 01 ,
     *-1.8369E 00 ,  8.6911E 00 , -1.3060E 01 ,  7.1880E 00 ,
c  j = 2;  N + N elastic scattering; 2.8 < Tlab <= 10. GeV:
     *-3.0853E 01 ,  1.0624E 02 , -1.2939E 02 ,  5.4339E 01 ,
     * 1.9465E 01 , -6.8102E 01 ,  9.6358E 01 , -5.6827E 01 ,
     *-3.4831E 00 ,  1.2341E 01 , -1.8592E 01 ,  1.2024E 01 ,
     * 1.8941E-01 , -6.7880E-01 ,  1.0665E 00 , -7.2910E-01 ,
c  j = 3;  n + p elastic scattering; Tlab <= 0.97 GeV:
     * 1.0258E-01 , -1.0542E 00 ,  1.1389E 01 , -1.6638E 01 ,
     *-4.9607E-01 ,  1.1800E 01 , -9.0857E 01 ,  1.6476E 02 ,
     * 1.5437E 00 , -3.3769E 01 ,  2.5192E 02 , -4.5071E 02 ,
     *-1.2021E 00 ,  2.5336E 01 , -1.8658E 02 ,  3.3254E 02 ,
c  j = 4; pi+ p or pi- n elastic scattering; Tlab <= 0.080 GeV:
     * 1.5789E-01 ,  2.9671E 00 , -5.5251E 00 ,  6.8925E 00 ,
     *-7.0218E 00 , -2.0534E 02 ,  5.6951E 02 , -8.9858E 02 ,
     * 1.3496E 02 ,  4.8722E 03 , -1.4674E 04 ,  2.3924E 04 ,
     *-8.2116E 02 , -3.2586E 04 ,  1.0098E 05 , -1.6553E 05 /
      DATA  A2/
c  j = 5; pi+ p or pi- n elastic scattering; 0.08 < Tlab <= 0.3 GeV:
     * 3.1531E-01 , -7.4981E 00 ,  4.3295E 01 , -7.6360E 01 ,
     *-6.5373E 00 ,  1.9307E 02 , -1.0181E 03 ,  1.7426E 03 ,
     * 4.6864E 01 , -1.3030E 03 ,  6.7291E 03 , -1.1075E 04 ,
     *-9.5192E 01 ,  2.6373E 03 , -1.2857E 04 ,  2.0294E 04 ,
c  j = 6; pi+ p or pi- n elastic scattering; 0.30 < Tlab <= 1.0 GeV:
     *-1.7953E 01 ,  1.0972E 02 , -2.3954E 02 ,  2.2826E 02 ,
     * 9.1968E 01 , -5.1963E 02 ,  1.1266E 03 , -1.0740E 03 ,
     *-1.3270E 02 ,  7.4112E 02 , -1.6000E 03 ,  1.5249E 03 ,
     * 5.8598E 01 , -3.1874E 02 ,  6.7751E 02 , -6.4011E 02 ,
c  j = 7; pi+ p or pi- n elastic scattering; 1.0 < Tlab <= 2.4 GeV:
     * 4.2169E-01 ,  1.4705E 02 , -6.5335E 02 ,  9.1507E 02 ,
     *-3.5198E 00 , -2.6019E 02 ,  1.2250E 03 , -1.7481E 03 ,
     * 3.6373E 00 ,  1.5592E 02 , -7.5201E 02 ,  1.0796E 03 ,
     *-7.8041E-01 , -3.0563E 01 ,  1.4795E 02 , -2.1250E 02 ,
c  j = 8; pi+ n or pi- p elastic scattering; Tlab <= 0.080 GeV:
     *-3.8288E-01 ,  3.7587E 00 , -6.5144E 00 ,  6.7740E 00 ,
     * 1.0381E 02 , -2.7282E 02 ,  4.7759E 02 , -5.1222E 02 ,
     *-1.7882E 03 ,  4.3052E 03 , -7.9314E 03 ,  9.3471E 03 ,
     * 7.1475E 03 , -3.3395E 03 , -4.1392E 03 , -4.4364E 03 /
      DATA  A3/
c  j = 9; pi- p or pi+ n elastic scattering; 0.08 < Tlab <= 0.3 GeV:
     * 2.4991E-01 ,  3.2028E 01 , -1.1882E 02 ,  1.5099E 02 ,
     *-2.6994E 00 , -4.6045E 02 ,  1.8959E 03 , -2.5190E 03 ,
     * 1.6268E 01 ,  2.1384E 03 , -9.1262E 03 ,  1.2431E 04 ,
     *-2.9654E 01 , -3.1823E 03 ,  1.3944E 04 , -1.9342E 04 ,
c  j = 10; pi- p or pi+ n elastic or CX scattering;
c          0.30 < Tlab <= 1.0 GeV:
     * 3.9025E 00 , -9.1126E 01 ,  3.2373E 02 , -4.0048E 02 ,
     *-2.0619E 01 ,  4.9170E 02 , -1.7155E 03 ,  2.1143E 03 ,
     * 3.3004E 01 , -7.6684E 02 ,  2.7003E 03 , -3.3525E 03 ,
     *-1.6367E 01 ,  3.7394E 02 , -1.3202E 03 ,  1.6423E 03 ,
c  j = 11; pi- p or pi+ n elastic or CX scattering; 
c          1.0 < Tlab <= 2.4 GeV:
     * 1.9402E 01 , -2.2446E 02 ,  7.4733E 02 , -9.3570E 02 ,
     *-4.4180E 01 ,  4.7194E 02 , -1.4856E 03 ,  1.8055E 03 ,
     * 3.1567E 01 , -3.0176E 02 ,  9.0763E 02 , -1.0773E 03 ,
     *-6.8648E 00 ,  6.0476E 01 , -1.7520E 02 ,  2.0381E 02 ,
c  j = 12; pi- + p --> pi0 n of pi+ + n --> pi0 + p Charge exchange
c          scattering; Tlab <= 0.08 GeV:
     * 1.4988E-01 ,  2.8753E 00 , -5.3078E 00 ,  6.2233E 00 ,
     *-5.9558E 00 , -1.6203E 02 ,  4.3079E 02 , -6.2548E 02 ,
     * 1.2875E 02 ,  3.1402E 03 , -7.9189E 03 ,  1.0983E 04 ,
     *-8.5161E 02 , -1.8780E 04 ,  4.4607E 04 , -5.8790E 04 /
      DATA  A4/
c  j = 13; pi- + p --> pi0 n of pi+ + n --> pi0 + p Charge exchange
c          scattering; 0.08 < Tlab <= 0.30 GeV:
     * 5.3689E-01 , -1.3216E 01 ,  8.1011E 01 , -1.4285E 02 ,
     *-1.0550E 01 ,  2.9629E 02 , -1.6957E 03 ,  2.8935E 03 ,
     * 6.9621E 01 , -1.9245E 03 ,  1.0620E 04 , -1.7468E 04 ,
     *-1.3865E 02 ,  3.9281E 03 , -2.0293E 04 ,  3.2058E 04 ,
c  j = 14;  N + N --> N + N + pi; nucleon distributions:
     * 8.5591E-02 ,  5.0390E 00 , -1.3782E 01 ,  1.4661E 01 ,
     * 5.4284E-02 , -9.2324E 00 ,  3.6397E 01 , -4.2962E 01 ,
     *-5.1111E-02 ,  4.6003E 00 , -2.0534E 01 ,  2.7731E 01 ,
     * 7.4514E-03 , -6.2529E-01 ,  2.9159E 00 , -4.1101E 00 ,
c  j = 15;  N + N --> N + N + pi; pion distributions:
     * 7.1622E-02 ,  3.0960E 00 , -1.1125E 01 ,  1.8130E 01 ,
     * 9.2581E-02 , -3.2186E 00 ,  2.0273E 01 , -3.3245E 01 ,
     *-5.1531E-02 ,  8.9886E-01 , -7.5084E 00 ,  1.3188E 01 ,
     * 5.8258E-03 , -1.7288E-03 ,  7.0224E-01 , -1.4854E 00 ,
c  j = 16;  N + N --> N + N + n*pi, n > 1; nucleon distributions:
     * 8.2300E-02 ,  1.5854E-01 ,  3.7716E 00 , -4.0562E 00 ,
     * 1.0802E-02 , -3.3688E-01 ,  1.1727E 00 , -6.7476E-01 ,
     *-2.1798E-03 ,  5.2166E-02 , -2.5816E-01 ,  3.2048E-01 ,
     * 6.5764E-05 , -1.4711E-03 ,  7.8209E-03 , -1.0580E-02 /
      DATA  A5/
c  j = 17;  N + N --> N + N + n*pi, n > 1; pion distributions:
     * 1.1138E-01 ,  6.0396E-01 ,  3.0174E 00 , -4.4190E 00 ,
     *-1.7709E-02 ,  2.3015E-01 , -1.8187E 00 ,  3.4518E 00 ,
     * 2.0977E-03 , -2.5458E-02 ,  2.1626E-01 , -4.0692E-01 ,
     *-5.4799E-05 ,  5.9111E-04 , -5.5552E-03 ,  1.0647E-02 ,
c  j = 18;  pi + N --> pi + N + pi; nucleon distributions:
     * 1.7288E-01 ,  7.1080E 00 , -1.7961E 01 ,  1.6403E 01 ,
     *-1.4504E-01 , -1.3032E 01 ,  4.1781E 01 , -4.0799E 01 ,
     * 4.5390E-02 ,  8.3515E 00 , -3.0260E 01 ,  3.2882E 01 ,
     *-4.7961E-03 , -1.4095E 00 ,  5.3505E 00 , -6.0946E 00 ,
c  j = 19;  pi + N --> pi + N + pi; pion distributions:
     * 3.7596E-02 ,  1.4331E 00 , -3.1350E 00 ,  6.4864E 00 ,
     * 2.3827E-01 ,  1.8253E 00 ,  1.7648E 00 , -1.6735E 01 ,
     *-1.5410E-01 , -1.5201E 00 , -1.5692E 00 ,  1.7185E 01 ,
     * 2.5037E-02 ,  3.0588E-01 ,  3.2520E-01 , -3.5277E 00 ,
c  j = 20;  pi + N --> pi + N + n*pi, n > 1; nucleon distributions:
     * 1.2489E-01 ,  1.3573E 00 ,  8.2338E-01 , -1.4595E 00 ,
     *-5.1577E-02 , -3.5778E-01 , -1.1690E 00 ,  1.8078E 00 ,
     * 7.4864E-03 ,  3.2888E-02 ,  2.3744E-01 , -3.9802E-01 ,
     *-2.9880E-04 , -7.5117E-04 , -1.1402E-02 ,  1.9505E-02 /
      DATA  A6/
c  j = 21;  pi + N --> pi + N + n*pi, n > 1; pion distributions:
     * 1.8470E-01 ,  1.9269E 00 , -3.2979E 00 ,  3.6843E 00 ,
     *-7.3932E-02 ,  2.7213E-01 ,  1.0600E 00 , -2.3354E 00 ,
     * 1.8907E-02 , -5.6473E-02 , -1.6487E-01 ,  3.8426E-01 ,
     *-9.2984E-04 ,  2.5506E-03 ,  7.3052E-03 , -1.7220E-02 /
      DATA  A7/
c  j = 22;  gamma + N --> pi0 + N, E-g <= 0.45 GeV. 
     & 4.0693D-01 , -4.1404D 00 ,  1.4044D 01 , -1.7265D 01 ,
     &-3.6799D 00 ,  5.9610D 01 , -1.6269D 02 ,  1.8873D 02 ,
     & 1.4556D 01 , -1.7550D 02 ,  4.5839D 02 , -5.3390D 02 ,
     &-1.2621D 01 ,  1.4964D 02 , -3.8118D 02 ,  4.5141D 02 ,
c  j = 23;  gamma + N --> pi0 + N, E-g > 0.45 GeV. 
     &-4.7554D-01 ,  2.2641D 00 , -1.2528D 01 ,  2.4647D 01 ,
     & 5.1620D 00 , -9.9236D 00 ,  5.5623D 01 , -1.0462D 02 ,
     &-8.1117D 00 ,  1.9315D 01 , -8.4255D 01 ,  1.3908D 02 ,
     & 3.5187D 00 , -9.1783D 00 ,  3.4950D 01 , -5.1243D 01 ,
c  j = 24; gamma + p --> n + pi+; E-g <= 0.51 GeV:
     & 4.8173D-01 ,  5.7726D 00 , -1.3745D 01 ,  2.7125D 01 ,
     &-4.4804D 00 , -3.8582D 01 ,  1.1159D 02 , -2.4305D 02 ,
     & 1.6306D 01 ,  1.1046D 02 , -3.3045D 02 ,  7.2270D 02 ,
     &-1.5968D 01 , -8.0140D 01 ,  2.4616D 02 , -6.0753D 02 ,
c  j = 25; gamma + p --> n + pi+; 0.51 < E-g <= 1.0 GeV:
     &-5.1646D 00 , -6.0776D 00 ,  7.8989D 01 , -1.0705D 02 ,
     & 2.1871D 01 ,  5.6915D 01 , -4.0159D 02 ,  5.1215D 02 ,
     &-2.7993D 01 , -9.4670D 01 ,  5.6928D 02 , -6.9621D 02 ,
     & 1.1587D 01 ,  4.5998D 01 , -2.4566D 02 ,  2.8452D 02 ,
c  j = 26; gamma + p --> n + pi+; 1.0 < E-g <= 10 GeV:
     &-5.3067D 01 ,  5.7612D 02 , -1.5438D 03 ,  1.6455D 05 ,
     & 1.4750D 02 , -1.6638D 03 ,  4.5923D 03 , -4.9949D 03 ,
     &-1.3436D 02 ,  1.5780D 03 , -4.4463D 03 ,  4.9022D 03 ,
     & 4.0253D 01 , -4.8860D 02 ,  1.4001D 03 , -1.5606D 03 ,
c  j = 27; gamma + N --> delta + pi; pion distribution; T < 1.0
     &-1.0306D 00 ,  3.2849D 01 , -7.5052D 01 ,  6.0255D 01 ,
     & 7.9586D 00 , -1.2572D 02 ,  2.5604D 02 , -1.6547D 02 ,
     &-1.4797D 01 ,  1.6590D 02 , -2.7991D 02 ,  1.1333D 02 ,
     & 8.2309D 00 , -6.7871D 01 ,  8.5762D 01 ,  5.9727D 00 ,
c  j = 28; gamma + N --> delta + pi; pion distribution; T > 1.0
     &-2.3722D 02 ,  9.6890D 02 , -1.6219D 03 ,  1.3637D 03 ,
     & 6.5800D 02 , -2.6941D 03 ,  4.5480D 03 , -3.8460D 03 ,
     &-6.0653D 02 ,  2.4983D 03 , -4.2498D 03 ,  3.6136D 03 ,
     & 1.8604D 02 , -7.6933D 02 ,  1.3166D 03 , -1.1242D 03 /
      DATA  A8/
c  j = 29; coefficients for K absorption ang. dist., Tlab <= 0.455
     & 6.5288D-01 ,  3.8977D-01 ,  8.4078D-01 ,  1.8893D-01 ,
     &-4.3964D 00 ,  3.4309D 01 , -7.3692D 01 ,  8.4308D 01 ,
     & 1.4889D 01 , -1.4380D 02 ,  3.1227D 02 , -3.5014D 02 ,
     &-1.5658D 01 ,  1.7160D 02 , -3.7212D 02 ,  4.1299D 02 /
      DATA  B1/
     * 5.0278E-01 ,  3.1442E 00 , -7.8172E 00 ,  8.1667E 00 ,
     * 9.3482E-01 , -1.0590E 01 ,  2.9227E 01 , -3.4550E 01 ,
     *-9.6685E-02 ,  4.7335E 00 , -1.4298E 01 ,  1.7685E 01 ,
     *-2.5041E-02 , -6.2478E-01 ,  2.0282E 00 , -2.5895E 00 ,
     * 1.1965E 00 , -8.2889E-01 ,  1.0426E 00 , -1.9090E 00 ,
     * 2.8703E-01 , -4.9065E 00 ,  1.6264E 01 , -1.9904E 01 ,
     *-2.4492E-01 ,  2.9191E 00 , -9.5776E 00 ,  1.1938E 01 ,
     * 3.7297E-02 , -4.2200E-01 ,  1.3883E 00 , -1.7476E 00 ,
     * 1.3508E 00 , -4.3139E 00 ,  1.2291E 01 , -1.5288E 01 ,
     *-2.0086E-01 ,  1.3641E 00 , -3.4030E 00 ,  3.8559E 00 ,
     * 1.2583E-02 , -8.3492E-02 ,  1.8600E-01 , -2.0043E-01 ,
     *-2.3628E-04 ,  1.3514E-03 , -2.4324E-03 ,  2.1906E-03 ,
     * 1.2419E 00 , -4.3633E 00 ,  1.3743E 01 , -1.8592E 01 ,
     *-2.4404E-01 ,  1.3158E 00 , -3.5691E 00 ,  4.3867E 00 ,
     * 1.5693E-02 , -8.2579E-02 ,  2.1427E-01 , -2.5846E-01 ,
     *-2.9386E-04 ,  1.4060E-03 , -3.3835E-03 ,  3.8664E-03 /
      DATA  B2/
     * 6.3054E-01 , -3.7333E 00 ,  1.3464E 01 , -1.8594E 01 ,
     *2.1801E 00  ,  1.5163E 00 , -1.6380E 01 ,  2.7944E 01 ,
     *-1.2886E 00 , -2.4570E 00 ,  1.5129E 01 , -2.3295E 01 ,
     * 2.0915E-01 ,  5.2279E-01 , -2.8687E 00 ,  4.2688E 00 ,
     * 9.3363E-01 , -1.8181E 00 ,  5.5157E 00 , -8.5216E 00 ,
     * 1.7811E 00 , -8.2927E 00 ,  2.0607E 01 , -2.0827E 01 ,
     *-1.5264E 00 ,  6.8433E 00 , -1.6067E 01 ,  1.6845E 01 ,
     * 2.7128E-01 , -1.1944E 00 ,  2.7495E 00 , -2.9045E 00 ,
     * 1.9439E 00 , -4.6268E 00 ,  9.7879E 00 , -9.6074E 00 ,
     *-3.4640E-01 ,  1.1093E 00 , -1.9313E 00 ,  1.7064E 00 ,
     * 2.7054E-02 , -1.1638E-01 ,  2.6969E-01 , -3.1853E-01 ,
     *-6.6092E-04 ,  5.0728E-03 , -1.4995E-02 ,  1.9605E-02 ,
     * 1.8693E 00 , -5.5678E 00 ,  1.4795E 01 , -1.6903E 01 ,
     *-4.9965E-01 ,  1.7874E 00 , -4.1330E 00 ,  3.8393E 00 ,
     * 4.6194E-02 , -1.8536E-01 ,  4.5315E-01 , -4.6273E-01 ,
     *-1.3341E-03 ,  5.7710E-03 , -1.4554E-02 ,  1.5554E-02 /
      DATA  CKJ/
     * 1.4509E-01 ,  4.6520E-01 , -3.3005E-02 ,  1.5376E-01 ,
     * 2.7436E-01 , -1.4604E-02 ,  6.2959E-01 ,  1.7866E-01 ,
     *-2.6216E-03 ,  8.3810E-01 ,  8.6137E-03 ,  3.2946E-03 ,
     * 9.2852E-02 ,  5.3886E-01 , -5.4493E-02 ,  1.3032E-01 ,
     * 4.0709E-01 , -2.8782E-02 ,  1.4909E-01 ,  3.8502E-01 ,
     *-1.2775E-02 ,  1.8024E-01 ,  3.3022E-01 , -9.4491E-03 /
      END
C     * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION SIGMAG (L,MS,MQ,KSI,IKS,T)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C     SUBPROGRAM OF CHOOSING CROSS SECTION TYPE AND
C     CALCULATION CROSS SECTION VALUE FOR GIVEN ENERGY.
c     EDITED by KKG, October, 2003 
c
      COMMON /TYPECSQ/ICST,NSICST
      DIMENSION ICST(46),NSICST(37)
c     data ICST/
c  Translation:total:|elas:   |total: |elas:  |total:  |elas:   |pi+ n or
c              p+p or|p+p or  |       |       |pi- p or|pi- p or|pi- p  |
c              n+n   |n+n     |p+n    |n+p    |pi+ n   |pi+ n   |SCX    |
c    &         210,    211,    220,    221,     120,     121,     122, 
c                           | Absorp: |
c  Translation:total:|elas: |pi+ np:pp|p + p->|p + p->|p + n->|p + n->|
c            pi+ p or|pi+ p |pi- pp:np|p p pi0|p n pi+|p + n +|p + p +|
c            pi- n   |pi- n |pi- pn:nn|n + n->|n + n->|pi0    |pi-    |
c                           |pi+ nn:np|n n pi0|p n pi-
c    &         110,    111,    123,    214,     215,     224,     225, 
c
c  Transl:   pi+ p->|pi+ p->|pi- p ->|pi- + p->|pi- + p->|gam + p|gam +p|
c         pi0 pi+ p |n +2pi+|2pi0 + n|p + pi- +|n pi- pi+|-> p + |-> n +|
c            pi- n->|pi- n->|pi+ n ->|pi0      |pi+ + n->| pi0   |pi+
c         pi0 pi- n |p +2pi-|2pi0 + p|         |p pi- pi+| (?)   | (?)
c    &         114,    115,    126,     124,      125,    10111,  10112, 
c            Absorp:
c  Transl:   gam +  |gam + p |gam + p |gam + p |gam + p |gam + p|gam + p|
c            2N ->  |-> pi+ p|->pi0 p |->n pi+ |-> delta|-> N + |total  |
c            2N     |+ pi-   |+ pi0   |+ pi0   |++ + pi-|2 pi   |       |
c    &       10113,  10115,  10114,   10116,    10118,   10117,   10110/
c
c  Transl:   total(elastic) for K + N
c            K+  + p or| K0  + p or| K-  + p or   |AK0 + p or  |
c            K0  + n   | K+  + n   |AK0  + n      | K- + n     | 
c            1110(1111)| 1120(1121)| -1120(-1121) |-1110(-1111)|
c
c  Transl:   charge exchange for K + N
c            K+  + p or| K0  + p<=>| K-  + p <==> |AK0 + p or  |
c            K0  + n   | K+  + n   |AK0  + n      | K- + n     | 
c              sig=0   |           |              |   sig=0.   |
c              1112    |    1122   |     -1122    |  -1112     |
c
c     data nsicst/   Absorp:  |  =0! |
c   Key:    pi+ p or|pi+ pn:pp|pi+ p:|pi+ p:|pi+ n:  |pi0 n |pi0 n  |
c            pi- n  |   or    |n 2pi0|N+ 2pi|N + 2pi | or   | or    |
c            SCX=0! |pi- pn:nn|pi- n:|pi- n:|pi- p:  |pi0 p |pi0 p  |
c                             |p 2pi0|N+ 2pi|N + 2pi |total |elastic|
c    &         112,    113,    116,    117,    127,    130,    131, 
c                    Absorp:
c   Key:     pi0+p: |pi0+pn:pn|pi0+p:|pi0+p:  |pi0+p:  |pi0+p: |p + p: |
c            pi+ + n|pi0+pp:pp|p 2pi0|p pi+pi-|n pi+pi0|N + 2pi| SCX   |
c            pi0+n: |pi0+np:np|pi0+n:|pi0+n:  |pi0+n:  |pi0+n: |n + n: |
c            pi- + p|pi0+nn:nn|n 2pi0|n pi+pi-|p pi-pi0|N + 2pi| = 0!  |
c    &         132,    133,    134,    135,     136,     137,    212, 
c
c   Key:     p + p: |p+p;n+n |p + p->|n + p: |n + p:  |n + p: |n + p->|
c      pi absorption|chg exch|2N + pi| SCX   | Pion   |n p pi0|N N pi |
c            n + n: |+ pi0   |n + n->|       |absorp. |Same as|       |
c             = 0!  |  = 0!  |2N + pi| = 0!  | = 0!   |  224! |       |
c    &         213,    216,    217,    222,     223,     226,    227/
c
c
c  Transl:      total(elastic) for Delta + N
c           D-  + p or| D0  + p or| D+  + p or| D++ + p or|
c           D++ + n   | D+  + n   | D0  + n   | D-  + n   |
c            230(231) |  240(241) |  250(251) |  260(261) |
c  
c  Transl:     charge exchange for Delta + N
c           D-  + p =>| D0  + p <=>| D+  + p <=>| D++ + p or|
c           D0  + n or| D+  + n    | D++ + n or | D-  + n   |
c           D++ + n =>|            | D0  + n => |  sig=0.   |
c           D+  + p   |            | D-  +p     |           | 
c            232      |  242       |     252    |  262      |
c  
c  Transl:     absorption  (Delta + N ==> N + N)
c           D-  + p =>| D0  + p  =>| D+  + p =>| D++ + p or|
c           n   + n or| p   + n or | p   + p or| D-  + n   |
c           D++ + n =>| D+  + n  =>| D0  + n =>|  sig=0.   |
c           p   + p   | p   + n    | n   + n   | |         | 
c            233      |  243       |     253   |  263      |
c  
c ======================================================================
      data absgam/5.0/,abspi/4.0/
c   Start subroutine:

      ICS = 10000*L+1000*IABS(MS)+100*MQ+10*KSI+IKS
      IF(MS.LT.0)  ICS=-ICS
      JS = 1
   17 IF (ICS-ICST(JS)) 11,10,11
   10 SIGMAG = QINTG(T,JS)
c     if(JS.eq.10)  SIGMAG = abspi*SIGMAG 
      if(JS.eq.40)  SIGMAG = absgam*SIGMAG 
      RETURN
   11 IF(JS-46)12,13,13
   12 JS = JS + 1
      GO TO 17
   13 NSJS = 1
   16 IF (ICS-NSICST(NSJS)) 15,14,15
   15 IF(NSJS-21)100,39,39
  100 NSJS=NSJS+1
      GO TO 16
   14 KNS = NSJS
      GO TO (18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     137,38),KNS
   18 SIGMAG = 0.
      RETURN
   19 SIGMAG = QINTG(T,10)
      RETURN
   20 SIGMAG = 0.
      RETURN
   21 SIGMAG = QINTG(T,15)+QINTG(T,16)
      RETURN
   22 SIGMAG = QINTG(T,18)+QINTG(T,19)+QINTG(T,17)
      RETURN
   23 SIGMAG = (QINTG(T,8)+QINTG(T,5))/2.
      RETURN
   24 SIGMAG = (QINTG(T,9)+QINTG(T,6)-QINTG(T,7))/2.
      RETURN
   25 SIGMAG = QINTG(T,7)
      RETURN
   26 SIGMAG = QINTG(T,10)*0.5
      RETURN
   27 SIGMAG = (QINTG(T,15)+QINTG(T,18))/2.
      RETURN
   28 SIGMAG = (QINTG(T,16)+QINTG(T,19))/2.
      RETURN
   29 SIGMAG = (QINTG(T,17))/2.
      RETURN
   30 SIGMAG = (QINTG(T,15)+QINTG(T,16)+QINTG(T,18)+
     1QINTG(T,19)+QINTG(T,17))/2.
      RETURN
   31 SIGMAG = 0.
      RETURN
   32 SIGMAG = 0.
      RETURN
   33 SIGMAG = 0.
      RETURN
   34 SIGMAG = QINTG(T,11)+QINTG(T,12)
      RETURN
   35 SIGMAG = 0.
      RETURN
   36 SIGMAG = 0.
      RETURN
   37 SIGMAG = QINTG(T,14)
      RETURN
   38 SIGMAG = QINTG(T,13)+2.*QINTG(T,14)
      RETURN
   39 DO  40  ISD=1,16
      KSD=ISD
      IF(ICS.EQ.NSICST(21+ISD))  GO  TO  41
   40 CONTINUE
      GO  TO  35
   41 CONTINUE
      W2=1.88*T+(1.232+0.940)**2
      TN=((SQRT(W2)+0.94-1.232)**2)/1.88-1.88
      PNN2=W2/4.-0.940**2
      PDN2=(W2+0.940**2-1.232**2)**2/(4.*W2)-0.940**2
      R=0.5*PNN2/PDN2
      HEX=0.5
C  * * * * * * *  PRELIMINARY VERSION  * * * * * * * * * *
      GO  TO  (230,231,232,233,
     *         240,241,242,243,
     *         250,251,252,253,
     *         260,261,262,263),KSD
  230 SIGMAG=QINTG(TN,3)
      RETURN
  240 GO  TO  230
  250 SIGMAG=QINTG(TN,1)
      RETURN
  260 GO  TO  250
  231 SIGMAG=QINTG(TN,4)*(1.-HEX)
      RETURN
  241 GO  TO  231
  251 SIGMAG=QINTG(TN,2)*(1.-HEX)
      RETURN
  261 GO  TO  251
  232 SIGMAG=QINTG(TN,4)*HEX
      RETURN
  242 GO  TO  232
  252 SIGMAG=QINTG(TN,2)*HEX
      RETURN
  262 GO  TO  263
  233 SIGMAG=QINTG(TN,12)/2.*R
      RETURN
  243 SIGMAG=(QINTG(TN,12)/2.+QINTG(TN,11))*R
      RETURN
  253 SIGMAG=(QINTG(TN,14)/2.+QINTG(TN,13)/2.)*R
      RETURN
  263 SIGMAG=0.
      RETURN
         END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION  CROSEG(L,MS,MQ,KSI,IKS,T,AM,IR)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Choose cross section type and calculate cross section for 
c     a given struck-nucleon rest frame kinetic energy.
c     l always 1 for photons; ms always 0; mb = baryon number.
c     mb may be 1 or 2.
c     ksi = 1 for n - n, p - p, pi+ - p & pi- - n.
c     ksi = 2 for n - p, pi+ - n & pi- - p.
c     ksi = 3 for pi0 - n, pi0 - p
c     iks = 0:   total cross section
c     iks = 1:   elastic cross section
c     iks = 2:   pion charge exchange cross section
c     iks = 3:   pion or gamma absorption cross section
c     iks = 4:   neutral pion production cross section
c     iks = 5:   charged pion production cross section
c     iks = 6:   target & projectile isospin change & neutral pion
c                production (9/11/97)
c     iks = 7:   total one-pion production cross section (3/15/99)
c     iks = 8:   delta production (by gamma) cross section
c
      CROSEG=0.
      X=T
      W=SQRT(1.88*T+(1.88+AM)**2)
      IF(IR.EQ.0.AND.MS.EQ.0)      GO  TO  10
      IF(IR.NE.0.AND.MS.EQ.0)      GO  TO  11
      IF(IR.EQ.0.AND.MS.NE.0)      GO  TO  12
      IF(IR.NE.0.AND.MS.NE.0)      GO  TO  14
      RETURN
   10 CROSEG=SIGMAG(L,MS,MQ,KSI,IKS,X)
      GO  TO  15
   11 IF(MQ.EQ.2)                  GO  TO  10
      IF(IKS.GE.3)                 RETURN
      T0=((W+0.140-AM)**2-(0.940+0.140)**2)/1.88
      X=T0
      GO  TO  10
   12 IF(MQ.NE.1)                  GO  TO  13
      P0=SQRT(T*(T+2.*AM))
      X=P0
      GO  TO  10
   13 IF(IKS.GT.1)                 RETURN
      TN=((W+0.940-AM)**2-(0.940+0.940)**2)/1.88
      TP=((W+0.140-AM)**2-(0.940+0.140)**2)/1.88
      TK=((W+0.492-AM)**2-(0.940+0.492)**2)/1.88
      PK=SQRT(TK*(TK+2.*0.492))
      SPP  =QINTG(TN, 1+IKS)
      SPN  =QINTG(TN, 3+IKS)
      SPIPP=QINTG(TP, 8+IKS)
      SPIMP=QINTG(TP, 5+IKS)
      SKMP =QINTG(PK,33+IKS)
      SKMN =QINTG(PK,31+IKS)
      SLP=SPN+SKMP-SPIMP
      SLN=SPP+SKMN-SPIPP
      IF(MS.EQ.-1.AND.KSI.EQ.1)  CROSEG=   SLP+SPP-SPN
      IF(MS.EQ.-1.AND.KSI.EQ.2)  CROSEG=   SLP+SPN-SPP
      IF(MS.EQ.-1.AND.KSI.EQ.3)  CROSEG=   SLP
      IF(MS.EQ.-1.AND.KSI.EQ.4)  CROSEG=   SLN+SPN-SPP
      IF(MS.EQ.-1.AND.KSI.EQ.5)  CROSEG=   SLN+SPP-SPN
      IF(MS.EQ.-1.AND.KSI.EQ.6)  CROSEG=   SLN
      IF(MS.EQ.-2.AND.KSI.EQ.1)  CROSEG=2.*SLP-SPP
      IF(MS.EQ.-2.AND.KSI.EQ.2)  CROSEG=2.*SLP-SPN
      IF(MS.EQ.-2.AND.KSI.EQ.3)  CROSEG=2.*SLN-SPN
      IF(MS.EQ.-2.AND.KSI.EQ.4)  CROSEG=2.*SLN-SPP
      IF(MS.EQ.-3.AND.KSI.EQ.1)  CROSEG=3.*SLP-SPP-SPN
      IF(MS.EQ.-3.AND.KSI.EQ.2)  CROSEG=3.*SLN-SPP-SPN
      GO  TO  15
   14 IF(MQ.EQ.2)           GO  TO  13
      T0=((W+0.492-AM)**2-(0.940+0.492)**2)/1.88
      P0=SQRT(T0*(T0+2.*0.492))
      X=P0
      GO  TO  10
   15 IF(CROSEG.LT.0.)  CROSEG=0.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ROTORQ (AR,BR,PSTAR,PR)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C    BLOCK OF ROTATION.
      DIMENSION AR(3),BR(3),PSTAR(3),PR(3),AN(3)
      SP = 0.
      DO 31 IR=1,3
      SP = SP+AR(IR)*BR(IR)
   31 CONTINUE
      AMOD = SQRT (AR(1)**2+AR(2)**2+AR(3)**2)
      IF(AMOD.LT.1.E-30)  GO  TO  10
      ALPHA1 = SP/AMOD
      BMOD2 = BR(1)**2+BR(2)**2+BR(3)**2
	TEMP=BMOD2-ALPHA1**2
      IF(TEMP.LE.0.D0)  GO  TO  10
      ALPHA2 = SQRT (TEMP)
      IF(ALPHA2.LT.1.E-30)  GO  TO  10
      AN(1) = AR(2)*BR(3)-AR(3)*BR(2)
      AN(2) = AR(3)*BR(1)-AR(1)*BR(3)
      AN(3) = AR(1)*BR(2)-AR(2)*BR(1)
      PR(1)=PSTAR(1)*BR(1)/ALPHA2+(PSTAR(3)-ALPHA1*PSTAR(1)/ALPHA2)
     1*AR(1)/AMOD+(PSTAR(2)*AN(1))/(ALPHA2*AMOD)
      PR(2)=PSTAR(1)*BR(2)/ALPHA2+(PSTAR(3)-ALPHA1*PSTAR(1)/ALPHA2)
     1*AR(2)/AMOD+(PSTAR(2)*AN(2))/(ALPHA2*AMOD)
      PR(3)=PSTAR(1)*BR(3)/ALPHA2+(PSTAR(3)-ALPHA1*PSTAR(1)/ALPHA2)
     1*AR(3)/AMOD+(PSTAR(2)*AN(3))/(ALPHA2*AMOD)
      RETURN
   10 DO  11  K=1,3
   11 PR(K)=PSTAR(K)
      RETURN
         END
c
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      SUBROUTINE BINELQ(PARTIN,IPATIN,IPATNE,L,MS,MQ,KSI,ME,V,U,
     *TIN1,MV,NP,NIN)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C     BLOCK OF INELASTIC SCATTERING at low energies TIN1 < 5-20 GeV
c
      COMMON /I3ACT/ SIG3,SIGIN,IND,ITH
      COMMON/NCASCA/NCAS,NCPRI
      DIMENSION PARTIN(9),IPATIN(5),IPATNE(5),V(3)
      data uthr/1.220/
      NIN = 0
      IK = 0
c  kkg 03/15/04
      if(IPATIN(2).ne.0)  then
        if(u.lt.uthr)  then
          NIN = 2
          return
        endif 
c                    gamma + N => Delta+pi or  => N+2pi 
        betais = QINTG(TIN1,44)/QINTG(TIN1,45)
        if(RNDM(-1.).le.betais)  then
c                    gamma + N => Delta+pi 
          call isobal(U,V,TIN1,PARTIN,IPATNE,MV,NP)
          return
        else
c                    gamma + N  => N+2pi 
          call statl(U,V,PARTIN,IPATNE,MV,NP)
          return
        endif
      endif
			 
c  kkg 03/15/04
      AM3=PARTIN(9)
      IF(IPATIN(5).EQ.0)  GO  TO  11
      IF(IPATIN(4).EQ.1)  AM3=0.940
      GO  TO  14
   11 IF(TIN1-4.)12,24,24
   24 BETATH=0.
                  GO TO 13
   12 SIG3=  CROSEG(L,MS,MQ,KSI,7,TIN1,PARTIN(9),IPATIN(5))
      SIGIN= CROSEG(L,MS,MQ,KSI,0,TIN1,PARTIN(9),IPATIN(5))-
     -       CROSEG(L,MS,MQ,KSI,1,TIN1,PARTIN(9),IPATIN(5))-
     -       CROSEG(L,MS,MQ,KSI,2,TIN1,PARTIN(9),IPATIN(5))
      BETATH=SIG3/SIGIN
   13 DRND=RNDM(-1.)
      IF(DRND-BETATH)14,15,15
   14 ITH=1
      TH=1.
                      GO TO 16
   15 ITH=0
      TH=2.
   16 IF(U- AM3     -0.14*TH-0.96)19,19,17
   17 CONTINUE
      IF(NCAS.GE.NCPRI) write( *,101)
      CALL VMNSPQ (PARTIN,IPATIN,U,MV,NP,ITH,MQ,TIN1,LP)
      IF (NP) 26,26,27
   26 RETURN
   27 CONTINUE
      IF (LP) 19,22,19
   22 CONTINUE
      IF(NCAS.GE.NCPRI) write( *,102)
      CALL DIRECTQ (V,TIN1,MQ,MV,NP,PARTIN,KP,ITH)
      IF (KP) 23,18,23
   18 CONTINUE
      IF(NCAS.GE.NCPRI) write( *,103)
      CALL CHINELQ (IPATIN,L,MS,MQ,KSI,NP,MV,TIN1,ME,IPATNE,PARTIN(9))
      IF(IPATIN(5).EQ.0.AND.ITH.EQ.1.AND.NCAS.GE.NCPRI) write( *,104)
      IF(IPATIN(5).EQ.0.AND.ITH.EQ.1)
     *CALL  DISOB(MV,U,IND,NP,MQ)
      RETURN
   23 IK = IK+1
      IF(IK.LT.50)   GO  TO  17
   19 NIN=2
  101 FORMAT(' ====> VMNSPQ')
  102 FORMAT(' ====> DIRECTQ')
  103 FORMAT(' ====> CHINEL')
  104 FORMAT(' ====> DISOB')
      RETURN
      END
c
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      subroutine STATL(u,v,partin,ipatne,mv,np)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     determining of secondary particles characteristics for
c     gamma-n interaction with statistical model.
c last modification: 15 Mar. 2004 by KKG
c
c
      dimension v(3),pv1(3),pv2(3),pv3(3),ps1(3),ps2(3),ps3(3),pin(3),
     1 pinst(3),pmemo(9,5999),imemo(5,5999),partin(9),ipatne(5)
      common /memory/ pmemo,imemo
      data emnucm, emnucg, emnucb / 938.919, 0.938919, 0.9315014/
      data emneut, emprot, empich, empi0 / 0.9395656, 0.9382723,
     &                                     0.139568, 0.134973/
      data  pi /3.1415926536d0/
c
c  determine randomly the resulting particle types 
      thrd = 1.d0/3.d0
      twthrd = 2.d0/3.d0
      temp4 = RNDM(-1.)
      imemo(2,mv+1) = 0
      imemo(3,mv+1) = 0
      imemo(4,mv+1) = 0
      imemo(5,mv+1) = 0
      imemo(2,mv+2) = 0
      imemo(3,mv+2) = 0
      imemo(4,mv+2) = 0
      imemo(5,mv+2) = 0
      imemo(2,mv+3) = 0
      imemo(3,mv+3) = 0
      imemo(4,mv+3) = 1
      imemo(5,mv+3) = 0
      if (ipatne(1).le.0) then
         if (temp4.le.thrd) then
            imemo(1,mv+1) = 0
            pmemo(9,mv+1) = empi0
            imemo(1,mv+2) = 0
            pmemo(9,mv+2) = empi0
            imemo(1,mv+3) = 0
            pmemo(9,mv+3) = emneut
         elseif (temp4.lt.twthrd) then
            imemo(1,mv+1) = 0
            pmemo(9,mv+1) = empi0
            imemo(1,mv+2) = -1
            pmemo(9,mv+2) = empich
            imemo(1,mv+3) = 1
            pmemo(9,mv+3) = emprot
         else
            imemo(1,mv+1) = -1
            pmemo(9,mv+1) = empich
            imemo(1,mv+2) = 1
            pmemo(9,mv+2) = empich
            imemo(1,mv+3) = 0
            pmemo(9,mv+3) = emneut
        endif
      else
         if (temp4.le.thrd) then
            imemo(1,mv+1) = 0
            pmemo(9,mv+1) = empi0
            imemo(1,mv+2) = 0
            pmemo(9,mv+2) = empi0
            imemo(1,mv+3) = 1
            pmemo(9,mv+3) = emprot
         elseif (temp4.le.twthrd) then
            imemo(1,mv+1) = 0
            pmemo(9,mv+1) = empi0
            imemo(1,mv+2) = 1
            pmemo(9,mv+2) = empich
            imemo(1,mv+3) = 0
            pmemo(9,mv+3) = emneut
         else
            imemo(1,mv+1) = 1
            pmemo(9,mv+1) = empich
            imemo(1,mv+2) = -1
            pmemo(9,mv+2) = empich
            imemo(1,mv+3) = 1
            pmemo(9,mv+3) = emprot
         endif
      endif
      np = 3
c
c  other combinations might be possible!?
      twopi = 2.0d0*pi
      empi = max(pmemo(9,mv+1),pmemo(9,mv+2))
      empi2 = empi**2
      emnu = pmemo(9,mv+3)
      emnupi = emnu + min(pmemo(9,mv+1),pmemo(9,mv+2))
c
      epim = (u**2+empi2-emnupi**2)/(2.d0*u)
      tpim = epim-empi
   10 t1 = RNDM(-1.)*tpim
      t2 = RNDM(-1.)*tpim
      e1 = t1+pmemo(9,mv+1)
      e2 = t2+pmemo(9,mv+2)
      f1 = 27.d0*e1*e2*(u-e1-e2)/(u**3)
      b1 = RNDM(-1.)
      if (b1.lt.f1)then
         e3 = u-e1-e2
         t3 = e3-pmemo(9,mv+3)
      else
         goto 10
      endif
      if (t3.le.0.0) goto 10
      p1 = SQRT(t1*(t1+2.0d0*pmemo(9,mv+1)))
      p2 = SQRT(t2*(t2+2.0d0*pmemo(9,mv+2)))
      p3 = SQRT(t3*(t3+2.0*pmemo(9,mv+3)))
      temp1 = (p1+p2-p3)*(p1-p2+p3)*(p2+p3-p1)
      if (temp1.le.0.0) goto 10
      ct3 = 1.d0-2.d0*RNDM(-1.)
      fi3 = twopi*RNDM(-1.)
      temp2 = SQRT(1.d0-ct3**2)
      pv3(1) = p3*temp2*COS(fi3)
      pv3(2) = p3*temp2*SIN(fi3)
      pv3(3) = p3*ct3
      temp3=SQRT(partin(8)*(partin(8)+2.d0*partin(9)))
      pin(1)=temp3*partin(4)*partin(7)
      pin(2)=temp3*partin(4)*partin(6)
      pin(3) = temp3*partin(5)
      call cmsq(pin,v,pinst,partin(8),partin(9))
      call ROTORQ (pinst,v,pv3,ps3)
      ct1 = (p2**2-p1**2-p3**2)/(2.d0*p3*p1)
      ct2 = (p1**2-p2**2-p3**2)/(2.d0*p3*p2)
      fi1 = twopi*RNDM(-1.)
      fi2 = pi+fi1
      st1 = SQRT(1.d0-ct1**2)
      st2 = SQRT(1.d0-ct2**2)
      pv1(1) = p1*st1*COS(fi1)
      pv1(2) = p1*st1*SIN(fi1)
      pv1(3) = p1*ct1
      call ROTORQ (ps3,v,pv1,ps1)
      pv2(1) = p2*st2*COS(fi2)
      pv2(2) = p2*st2*SIN(fi2)
      pv2(3) = p2*ct2
      call ROTORQ (ps3,v,pv2,ps2)
      pmemo(4,mv+1) = ps1(1)
      pmemo(5,mv+1) = ps1(2)
      pmemo(6,mv+1) = ps1(3)
      pmemo(7,mv+1) = 0.
      pmemo(4,mv+2) = ps2(1)
      pmemo(5,mv+2) = ps2(2)
      pmemo(6,mv+2) = ps2(3)
      pmemo(7,mv+2) = 0.
      pmemo(4,mv+3) = ps3(1)
      pmemo(5,mv+3) = ps3(2)
      pmemo(6,mv+3) = ps3(3)
      pmemo(7,mv+3) = 0.
      return
      end
c
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      subroutine isobal (u,v,tin1,partin,ipatne,mv,np)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     determining of secondary particles characteristics for gamma-n
c     interaction with (3/2,3/2) isobar production.
c
C   Last change: 15-MAR-2004 by KKG
c
      dimension vt(3),ppim(3),ppt(3),ppit(3),ppi(3),pp(3),pin(3),
     & pinst(3),ppimst(3),ppist(3),ppst(3),
     & pmemo(9,5999),imemo(5,5999),
     & partin(9),ipatne(5),v(3)
      common /memory/ pmemo,imemo
      data emnucm, emnucg, emnucb / 938.919, 0.938919, 0.9315014/
      data emneut, emprot, empich, empi0 / 0.9395656, 0.9382723,
     &                                     0.139568, 0.134973/
      data  pi /3.1415926536d0/
c
      imemo(2,mv+1) = 0
      imemo(3,mv+1) = 0
      imemo(4,mv+1) = 0
      imemo(5,mv+1) = 0
      imemo(2,mv+2) = 0
      imemo(3,mv+2) = 0
      imemo(4,mv+2) = 0
      imemo(5,mv+2) = 0
      imemo(2,mv+3) = 0
      imemo(3,mv+3) = 0
      imemo(4,mv+3) = 1
      imemo(5,mv+3) = 0
      if (ipatne(1).gt.0) then
         imemo(1,mv+1) = -1
         pmemo(9,mv+1) = empich
         imemo(1,mv+2) = 1
         pmemo(9,mv+2) = empich
         imemo(1,mv+3) = 1
         pmemo(9,mv+3) = emprot
      else
         imemo(1,mv+1) = 1
         pmemo(9,mv+1) = empich
         imemo(1,mv+2) = -1
         pmemo(9,mv+2) = empich
         imemo(1,mv+3) = 0
         pmemo(9,mv+3) = emneut
      endif
      np = 3
c
      twopi = 2.0d0*pi
      empi2 = empich**2
      emnu = pmemo(9,mv+3)
      emnupi = emnu + empich
      emnu2pi = emnupi + empich
      a1 = (u**2+empi2-emnupi**2)/(2.d0*u)
      a2 = SQRT(a1**2-empi2)
      a3 = u-a1
      f1 = a1*a2*a3/u
      alpha = 200.d0*f1
   10 bms = RNDM(-1.)*(u-emnu2pi)+emnupi
      epim = (u**2+empi2-bms**2)/(2.d0*u)
      edn = u-epim
      pim = SQRT(epim**2-empi2)
      f = (pim*epim*edn)/u
      ts = (bms**2-emnupi**2)/emnu/2.d0
      p = f*QINTG(ts,9)/alpha
      b1 = RNDM(-1.)
      if (p.gt.b1) then
         tpi = epim-empich
      else
         goto 10
      endif
      if (tin1.lt.1.d0) then
         ctpi = COSTAQ(27,tin1)
      else
         ctpi = COSTAQ(28,tin1)
      endif
      fipi = twopi*RNDM(-1.)
      epit = (bms**2+empi2-emnu**2)/(2.d0*bms)
      ent = bms-epit
      temp1 = SQRT(1.d0-ctpi**2)
      ppim(1) = pim*temp1*COS(fipi)
      ppim(2) = pim*temp1*SIN(fipi)
      ppim(3) = pim*ctpi
      temp2 = tpi+empich-u
      vt(1) = ppim(1)/temp2
      vt(2) = ppim(2)/temp2
      vt(3) = ppim(3)/temp2
      ctilpi=1.d0-2.d0*RNDM(-1.)
      ftilpi=twopi*RNDM(-1.)
      temp1 = SQRT(1.d0-ctilpi**2)
      pmt = SQRT(epit**2-empi2)
      ppit(1) = pmt*temp1*COS(ftilpi)
      ppit(2) = pmt*temp1*SIN(ftilpi)
      ppit(3) = pmt*ctilpi
      tpi = epit-empich
      vt(1)=-vt(1)
      vt(2)=-vt(2)
      vt(3)=-vt(3)
      call cmsq (ppit,vt,ppi,tpi,empich)
      ppt(1) = -ppit(1)
      ppt(2) = -ppit(2)
      ppt(3) = -ppit(3)
      tn = ent-emnu
      call cmsq (ppt,vt,pp,tn,emnu)
      temp2 = SQRT(partin(8)*(partin(8)+2.d0*partin(9)))
      pin(1)=temp2*partin(4)*partin(7)
      pin(2)=temp2*partin(4)*partin(6)
      pin(3) = temp2*partin(5)
      call cmsq(pin,v,pinst,partin(8),partin(9))
      call ROTORQ(pinst,v,ppim,ppimst)
      pmemo(4,mv+1) = ppimst(1)
      pmemo(5,mv+1) = ppimst(2)
      pmemo(6,mv+1) = ppimst(3)
      pmemo(7,mv+1) = 0.
      call ROTORQ(pinst,v,ppi,ppist)
      pmemo(4,mv+2) = ppist(1)
      pmemo(5,mv+2) = ppist(2)
      pmemo(6,mv+2) = ppist(3)
      pmemo(7,mv+2) = 0.
      call ROTORQ(pinst,v,pp,ppst)
      pmemo(4,mv+3) = ppst(1)
      pmemo(5,mv+3) = ppst(2)
      pmemo(6,mv+3) = ppst(3)
      pmemo(7,mv+3) = 0.
c
      return
      end

c
c     ********************************************************************
c

      subroutine cmsq (p, v, pstar, t, cm)

c ======================================================================
c
c     Momentum calculation in system which has a relative velocity
c     v to given one.
c     Lorentz transformation; see Jackson, 2nd ed., p541.
c
c   Called by: ISOBAR STAT
c
c   Edited by AJS, August, 1997.
c
C   Last change: 15-MAR-2004 by KKG
c ======================================================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      dimension p(3), v(3), pstar(3)

c ======================================================================

      v2 = v(1)**2 + v(2)**2 + v(3)**2
      temp1 = sqrt(1.0 - v2)
      spv = p(1)*v(1) + p(2)*v(2) + p(3)*v(3)
      temp2 = spv/v2*(1.0/temp1 - 1.0)
      sv = (t + cm)/temp1

C should bellow be "+" instead of "-" ?, SGM, 05/25/03

      pstar(1) = p(1) + v(1)*temp2 - v(1)*sv
      pstar(2) = p(2) + v(2)*temp2 - v(2)*sv
      pstar(3) = p(3) + v(3)*temp2 - v(3)*sv
      return

c ======================================================================
      end

C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION  COSAL(T)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c    sampling cos(theta) for elastic N+He-4 scattering
c
      IF(T-0.147) 1,1,2
    1 A1=25.2*T**(-0.843)
                              GO TO 3
    2 A1=130.*T**0.0145
    3 A2=11.3*T**0.432
      IF(T-0.055) 4,4,5
    4 A3=0.22*T**(-1.35)
                             GO TO 6
    5 A3=0.000043*T**(-4.32)
    6 A4=130.*T**1.33
      TEMP1=1.-EXP(-A2*3.141592)
      TEMP3=2.-TEMP1
      TEMP2=1.-EXP(-A4*3.141592)
      TEMP4=2.-TEMP2
      W=A3*TEMP4/((1.+A4**2)*(A1*TEMP3/(1.+A2**2)+A3*TEMP4/(1.+A4**2)))
      DRND=RNDM(-1.)
      IF(DRND-W) 7,7,9
    7 TAU=-LOG(1.-RNDM(-1.)*TEMP2)/A4
      DRND=RNDM(-1.)
      IF(DRND-SIN(TAU))8,8,7
    8 TETA=3.141592-TAU
                            GO TO 10
    9 TETA=-LOG(1.-RNDM(-1.)*TEMP1)/A2
      DRND=RNDM(-1.)
      IF(DRND-SIN(TETA)) 10,10,9
   10 COSAL= COS(TETA)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION COSTAQ(J,T)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Cosine calculation for elastic and charge-exchange reactions.
c
      COMMON/COEFAQ/ANKJ
c   kkg  10/28/03
      DIMENSION ANKJ(4,4,29),ANK(4,4)
      IF(T.LT.0.001)  GO  TO  14
      DO 10 K=1,4
      DO 10 N=1,4
      ANK(N,K) = ANKJ(N,K,J)
   10 CONTINUE
      S1 = 0.
      R1 = RNDM(-1.)
      IF(R1.lt.1.E-10) R1 = RNDM(-1.)
      S2 = 0.
      DO 11 N=1,4
      DO 11 K=1,4
      S1 = S1+ANK(N,K)*(T**(K-1))*(R1**(N-1))
   11 CONTINUE
      DO 12 N=1,4
      DO 12 K=1,4
      S2 = S2+ANK(N,K)*T**(K-1)
   12 CONTINUE
      CTA = 2.*SQRT(R1)*(S1+(1.-S2)*R1**4)-1.
      TEMP1 = ABS(CTA)
      IF (TEMP1-1.) 13,13,14
   13 COSTAQ = CTA
                    RETURN
   14 COSTAQ=1.-2.*RNDM(-1.)
      RETURN
         END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION qintg (x,lq)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      logical ilog
      COMMON /TABELE/sigma,argus
      DIMENSION sigma(30,46),argus(30,9)
c      Parabolic interpolation of tabulated cross sections.
c      For lq = 1-4 and 15-19, use a parabola in log(sigma) vs. log(E)
c      for energies (x) less than the 2nd point.
c   Modified for new cross section tables using AJS's sugestion from 
c   cem2k, October, 2003, KKG 
c ======================================================================

      lpmax = 30
      ilog = .false.
      if (lq.le.4) then
c  1-4
        i = 1
c  kkg 24.03.06
        if(x.gt.argus(30,i))  then 
           qintg = hexsecg(x,lq)
           return
        endif 
        ilog = .true.
      elseif (lq.le.10) then
c  5-10
         i = 2
c  kkg 24.03.06
         if(lq.ge.5.and.lq.le.9.and.x.gt.argus(30,i))  then
            qintg = hexsecg(x,lq)
            return 
         endif
      elseif (lq.le.14) then
c  11-14
        i = 3
      elseif (lq.le.19) then
c  15-19
        i = 3
         ilog = .true.      
      elseif(lq.le.23)  then
c  20-23
        i = 4
      elseif(lq.le.26)  then
c  24-26
        i = 5
      elseif(lq.le.37) then
        i = 6
c  27-37
      elseif(lq.le.39.or.lq.eq.46)  then
c  38,39,46
        i = 7
         ilog = .true.
      elseif(lq.eq.40)  then
c  40
        i = 8
      elseif(lq.le.45)  then
c  41-45
        i = 9
         ilog = .true.
      else
        write(*,1000)  lq
        stop
      endif
c ---------------------------------      
      lpha = 1
   10 continue
      if (x.eq.argus(lpha,i)) then
c  If x = table value, set cross section to tabulated one.
        qintg = sigma(lpha,lq)
        return
      elseif (x.gt.argus(lpha,i)) then
c  Increment lpha until x is <= tabulated x.
        if (lpha.ge.lpmax-1) then
          lpha = lpmax - 1
          go to 20
        else
          lpha = lpha + 1
          go to 10
        endif
      elseif (x.lt.argus(lpha,i)) then
        if (lpha.le.1 .and. .not.ilog .and. i.eq.3) then
c  If x < first table value and not using parabolic log interpolation;
c  below threshold; set cross section to 0.
          qintg = 0.0
          return
        elseif (lpha.le.1 .and. i.eq.2) then
c  For pi+N reactions with energy < 0; use 0-energy cross section:
          qintg = sigma(1,lq)
          return
        else
          lpha = max(2,lpha)
          lpha = min(lpmax-1,lpha)
        endif
      endif
   20 phi1 = sigma(lpha-1,lq)
      psi1 = argus(lpha-1,i)
      phi2 = sigma(lpha,lq)
      psi2 = argus(lpha,i)
      phi3 = sigma(lpha+1,lq)
      psi3 = argus(lpha+1,i)
      if (ilog .and. lpha.eq.2) then
        x = log(x)
        phi1 = log(sigma(lpha-1,lq))
        psi1 = log(max(argus(lpha-1,i), 1.e-6))
        phi2 = log(sigma(lpha,lq))
        psi2 = log(argus(lpha,i))
        phi3 = log(sigma(lpha+1,lq))
        psi3 = log(argus(lpha+1,i))
      endif
      a = psi2 - psi3
      b = psi3 - psi1
      c = psi1 - psi2
      delta = a*psi1**2 + b*psi2**2 + c*psi3**2
      deltaa = phi1*a + phi2*b + phi3*c
      deltab = (phi2 - phi3)*psi1**2 + (phi3 - phi1)*psi2**2 +
     &         (phi1 - phi2)*psi3**2
      deltac = (psi2*phi3 - psi3*phi2)*psi1**2 + 
     &         (psi3*phi1 - psi1*phi3)*psi2**2 +
     &         (psi1*phi2 - psi2*phi1)*psi3**2
      fact = 0.0
      if (delta.ne.0.0) fact = 1.0/delta
      a = deltaa*fact
      b = deltab*fact
      c = deltac*fact
      qintg = a*x**2 + b*x + c
      if (ilog .and. lpha.eq.2) then
        x = exp(x)
        qintg = exp(qintg)
      else
        qintg = max(qintg, 0.0)
      endif
      return

c ======================================================================

 1000 format (1x,'qintg called with improper value of lq = ',i4)

c ======================================================================
      end
c
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      DOUBLE PRECISION function hexsecg(T,I) 
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c  kkg  10.03.06, 04/12/07
c  high energy approximation of total and elastic cross sections
c  pi+p, pi-p, pp,np
      real*8 mpi,mn
      dimension p1(8),p2(8),p3(8),p4(8)
      data mpi/0.139/,mn/0.939/
      data s1/1.00/,ss0/5.380/,et1/0.4580/,et2/0.5450/
      data p1/19.912,-0.5032,21.122,-1.4187,35.530,5.3828,34.425,35.80/
      data p2/5.4529,8.85970,4.7596,11.2230,4.3527,3.9741,9.6702,6.336/
      data p3/3.6838,8.66200,1.2521,11.5740,0.6752,0.4921,9.8362,5.477/
      data p4/0.6230,0.48470,0.5390,0.50790,0.5555,0.2944,0.6545,0.308/
      if(I.le.4)                then
        k=I+4
        s=2.0*T*mn+(mn+mn)**2 
      elseif(I.eq.5.or.I.eq.6)  then
        k=I-2 
        s=2.0*T*mn+(mpi+mn)**2 
      elseif(I.ge.8)            then 
        k=I-7
        s=2.0*T*mn+(mpi+mn)**2 
      else 
      endif
      s0=ss0**2
      Z = p1(k)
      Y1= p2(k)**2
      Y2= p3(k)**2
      B = p4(k)**2  
      hexsecg=Z+B*LOG(s/s0)**2+Y1*(s1/s)**et1-Y2*(s1/s)**et2
      return
      end
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE SLQEKQ(L,MS,MQ,KSI,ME,IIN,IPN)
c
c     Determine cross section type.
C     Edited by KKG, October, 2003
      DIMENSION IIN(5),IPN(5)
c   |  i   |     IIN(i) or IPN(i) translation
c   |  1   |             charge
c   |  2   |           lepton number
c   |  3   |             strangeness 
c   |  4   |           baryon number
c   |  5   |zero for stable particles, INT(1000*t_life) for resonances  
c                  
      MEIN=IIN(1)
c  Total strangenes in system:
      MS = IIN(3)+IPN(3)
c  Total lepton number in system:
      L =  IIN(2)+IPN(2)
c  Total baryon number in system:
      MQ = IIN(4)+IPN(4)
c  Total charge in system:
      ME = IIN(1)+IPN(1)
c    condition for Delta's
      IF(IIN(5).NE.0.AND.IIN(4).EQ.1.AND.IIN(3).EQ.0)  GO  TO  30
    9 CONTINUE
      IF (MS) 10,11,10
   10 IF(IIN(4).NE.0)    GO  TO  110
c    strange mesons + p
      IF(MS.EQ. 1.AND.IIN(1).EQ. 1.AND.IPN(1).EQ.1)    KSI=1
      IF(MS.EQ. 1.AND.IIN(1).EQ. 0.AND.IPN(1).EQ.1)    KSI=2
      IF(MS.EQ.-1.AND.IIN(1).EQ.-1.AND.IPN(1).EQ.1)    KSI=2
      IF(MS.EQ.-1.AND.IIN(1).EQ. 0.AND.IPN(1).EQ.1)    KSI=1
c    strange mesons + n
      IF(MS.EQ. 1.AND.IIN(1).EQ. 1.AND.IPN(1).EQ.0)    KSI=2
      IF(MS.EQ. 1.AND.IIN(1).EQ. 0.AND.IPN(1).EQ.0)    KSI=1
      IF(MS.EQ.-1.AND.IIN(1).EQ.-1.AND.IPN(1).EQ.0)    KSI=1
      IF(MS.EQ.-1.AND.IIN(1).EQ. 0.AND.IPN(1).EQ.0)    KSI=2
      RETURN
c    strange baryon + p or n
  110 IF(MS.EQ.-1.AND.IIN(1).EQ. 1.AND.IPN(1).EQ.1)    KSI=1
      IF(MS.EQ.-1.AND.IIN(1).EQ. 1.AND.IPN(1).EQ.0)    KSI=4
      IF(MS.EQ.-1.AND.IIN(1).EQ.-1.AND.IPN(1).EQ.1)    KSI=2
      IF(MS.EQ.-1.AND.IIN(1).EQ.-1.AND.IPN(1).EQ.0)    KSI=5
      IF(MS.EQ.-1.AND.IIN(1).EQ. 0.AND.IPN(1).EQ.1)    KSI=3
      IF(MS.EQ.-1.AND.IIN(1).EQ. 0.AND.IPN(1).EQ.0)    KSI=6
      IF(MS.EQ.-2.AND.IIN(1).EQ.-1.AND.IPN(1).EQ.1)    KSI=1
      IF(MS.EQ.-2.AND.IIN(1).EQ.-1.AND.IPN(1).EQ.0)    KSI=3
      IF(MS.EQ.-2.AND.IIN(1).EQ. 0.AND.IPN(1).EQ.1)    KSI=2
      IF(MS.EQ.-2.AND.IIN(1).EQ. 0.AND.IPN(1).EQ.0)    KSI=4
      IF(MS.EQ.-3.AND.IIN(1).EQ.-1.AND.IPN(1).EQ.1)    KSI=1
      IF(MS.EQ.-3.AND.IIN(1).EQ.-1.AND.IPN(1).EQ.0)    KSI=2
      RETURN
   11 IF (L) 13,13,12
c   gamma + p or n
   12 KSI = 1
      RETURN
   13 IF(MQ-1) 15,15,100
c   baryon on baryon
  100 IF(MQ-2) 14,14,16
   14 IF (ME-1)16,17,16
c   2 neutrons or 2 protons
   16 KSI = 1
      RETURN
c   1 neutron and 1 proton
   17 KSI = 2
      RETURN
c   meson interacting with baryon
   15 IF (ME-2) 19,18,19
c   pi+ on proton 
   18 KSI = 1
      RETURN
   19 IF (ME+1) 21,20,21
c   pi- on neutron
   20 KSI = 1
      RETURN
   21 IF (ME) 22,23,22
   22 IF (MEIN-1) 27,26,27
c   pi+ incident on neutron
   26 KSI = 2
      RETURN
c   pi0 incident on proton 
   27 KSI = 3
      RETURN
   23 IF (MEIN+1) 25,24,25
c   pi- incident on proton
   24 KSI = 2
      RETURN
c   pi0 incident on neutron
   25 KSI = 3
      RETURN
c   Delta + proton or neutron
c                    p             n
c     D-             3             6
c     D0             4             5
c     D+             5             4
c     D++            6             3
c    
   30 KSI=5+IIN(1)*(IPN(1)-1)+IPN(1)*(IIN(1)-1)
      RETURN
         END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE CINEMA(PSTAR,V,P,CT,ST,CFI,SFI,T,CM)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C     KINEMATIC BLOCK.
      DIMENSION PSTAR(3),V(3),P(3)
      SPV = PSTAR(1)*V(1)+PSTAR(2)*V(2)+PSTAR(3)*V(3)
      G=1./SQRT(1.d0-V(1)**2-V(2)**2-V(3)**2)
      ESTAR=SQRT(PSTAR(1)**2+PSTAR(2)**2+PSTAR(3)**2+CM**2)
      DO  K=1,3
        P(K)=PSTAR(K)+G*V(K)*(SPV*G/(G+1.)+ESTAR)
      ENDDO
      PM = SQRT(P(1)**2+P(2)**2+P(3)**2)
      IF(PM.LT.1.E-30)   THEN
        CT=1.d0
        CFI=1.d0
        SFI=0.d0
        ST=0.d0
      ELSE
        CT = P(3)/PM
        TEMP4 = 1.d0-CT**2
        if(TEMP4.lt.1.E-30) then
          CT=1.d0
          CFI=1.d0
          SFI=0.d0
          ST=0.d0
        else
          ST = SQRT(TEMP4)
          TEMP3 = PM*ST
          CFI=P(1)/TEMP3
          SFI=P(2)/TEMP3
        endif
      ENDIF
      T=SQRT(PM**2+CM**2)-CM
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION PMOMQ (J,T)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      REAL*8 mn,mpi
C     BLOCK OF CALCULATION OF SECONDARY PARTICLES MOMENTUM.
      COMMON/COEFBCQ/BNKJ,CKJ
      DIMENSION BNKJ(4,4,8),CKJ(3,8),BNK(4,4)
      DATA mn/0.939/,mpi/0.139/
c
      DO 10 K=1,4
      DO 10 N=1,4
      BNK(N,K) = BNKJ(N,K,J)
   10 CONTINUE
      S1 = 0.
      R1 = RNDM(-1.)
      S2 = 0.
      Pmax = 0.
      DO 11 N=1,4
      DO 11 K=1,4
      S1 = S1+BNK(N,K)*(T**(K-1))*(R1**(N-1))
   11 CONTINUE
      DO 12 N=1,4
      DO 12 K=1,4
      S2 = S2+BNK(N,K)*T**(K-1)
   12 CONTINUE
c  calculation of maximal momentum
c  back to old version, KKG, Dec. 2004  
        do K=1,3
          Pmax = Pmax+CKJ(K,J)*T**(K-1)
	enddo
c
      PMOMQ = Pmax*SQRT(R1)*(S1+(1.-S2)*R1**4)
      RETURN
         END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      INTEGER FUNCTION JTYPAQ (ITH,MQ,LAMB)
C     DETERMINING OF TYPE OF COEFFICIENTS A(N,K).
      IF (ITH) 10,18,10
   10 IF (MQ-1) 11,11,14
   11 IF (LAMB-1) 13,13,12
   12 JTYPAQ = 19
                   RETURN
   13 JTYPAQ = 18
                   RETURN
   14 IF (LAMB-1) 16,16,15
   15 IF (LAMB-3) 17,16,17
   16 JTYPAQ = 14
                   RETURN
   17 JTYPAQ = 15
                   RETURN
   18 IF (MQ-1) 19,19,22
   19 IF (LAMB-1) 20,20,21
   20 JTYPAQ = 20
                   RETURN
   21 JTYPAQ = 21
                   RETURN
   22 IF (LAMB-1) 24,24,23
   23 IF (LAMB-3) 25,24,25
   24 JTYPAQ = 16
                   RETURN
   25 JTYPAQ = 17
                   RETURN
         END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      INTEGER FUNCTION JTYPBQ (ITH,MQ,LAMB)
C     DETERMINING OF TYPE OF COEFFICIENTS B(N,K).
      IF (ITH) 10,18,10
   10 IF (MQ-1) 11,11,14
   11 IF (LAMB-1) 13,13,12
   12 JTYPBQ = 6
                  RETURN
   13 JTYPBQ = 5
                  RETURN
   14 IF (LAMB-1) 16,16,15
   15 IF (LAMB-3) 17,16,17
   16 JTYPBQ = 1
                  RETURN
   17 JTYPBQ = 2
                  RETURN
   18 IF (MQ-1) 19,19,22
   19 IF (LAMB-1) 20,20,21
   20 JTYPBQ = 7
                  RETURN
   21 JTYPBQ = 8
                  RETURN
   22 IF (LAMB-1) 24,24,23
   23 IF (LAMB-3) 25,24,25
   24 JTYPBQ = 3
                  RETURN
   25 JTYPBQ = 4
                  RETURN
         END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE COTRAN
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c   Changes the projectile an target Coulomb trajectoty 
c
      COMMON/DTINT/DTAU
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      DIMENSION VCM(3),P(3),P1(3),P2(3),R(3),B(3),RN(3)
      IF((AN1+AN2).EQ.0.)   RETURN
      AMR=AN1*AN2/(AN1+AN2)*0.940
      IF(AMR.LE.0.)  RETURN
      IF(DTAU.LT.0.) RETURN
      AMP=AN1*.940
      AMT=AN2*.940
      GP=1./SQRT(1.-VPR(1)**2-VPR(2)**2-VPR(3)**2)
      GT=1./SQRT(1.-VTA(1)**2-VTA(2)**2-VTA(3)**2)
      G1=GP*AN1/(GP*AN1+GT*AN2)
      G2=GT*AN2/(GP*AN1+GT*AN2)
      DO 10  K=1,3
      RN(K)=RADP(K)-RADT(K)
      R(K)=RN(K)-(VPR(K)-VTA(K))*DTAU
      VCM(K)=G1*VPR(K)+G2*VTA(K)
      B(K)=-VCM(K)
   10 P1(K)=GP*AMP*VPR(K)
      CALL CINEMA(P1,B,P,CT,ST,CF,SF,T,AMP)
      P0=SQRT(P(1)**2+P(2)**2+P(3)**2)
      R0=SQRT(R(1)**2+R(2)**2+R(3)**2)
      R1=SQRT(RN(1)**2+RN(2)**2+RN(3)**2)
      PRN=P(1)*RN(1)+P(2)*RN(2)+P(3)*RN(3)
      PL2=(PRN/R1)**2-2.*AMR*(POTNU(R1)-POTNU(R0))
      SPL=1.
      IF(PRN.LT.0.)   SPL=-1.
      IF(PL2.LT.0.)   PL2=0.
      TEMP=(SPL*SQRT(PL2)-PRN/R1)/R1
      DO  11  K=1,3
   11 P(K)=P(K)+TEMP*RN(K)
      CALL  CINEMA(P,VCM,P1,CT1,ST1,CF1,SF1,T1,AMP)
      P(1)=-P(1)
      P(2)=-P(2)
      P(3)=-P(3)
      CALL  CINEMA(P,VCM,P2,CT2,ST2,CF2,SF2,T2,AMT)
      EP=AMP+T1
      ET=AMT+T2
      DO  12  K=1,3
      VPR(K)=P1(K)/EP
   12 VTA(K)=P2(K)/ET
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION  FORCEN (R1,R2)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c    Calculates Coulomb force
c
      IF(ABS(R1-R2).EQ.0.)  GO  TO 10
      FORCEN=-(POTNU(R2)-POTNU(R1))/(R2-R1)
      RETURN
   10 FORCEN=0.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION  POTNU(R)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c    Calculates Coulomb potential acting between projectile and target
c
c
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      R1=0.
      IF(AN1.GT.0.1)  S1=R0N1*AN1**(1./3.)
      IF(AN1.GT.0.1)  R1=S1*(1.+(C1*3.141592/S1)**2)**(1./3.)
      R2=0.
      IF(AN2.GT.0.1)  S2=R0N2*AN2**(1./3.)
      IF(AN2.GT.0.1)   R2=S2*(1.+(C2*3.141592/S2)**2)**(1./3.)
      IF((R1+R2).LE.0.)  GO  TO  14
      RP=R1
      IF(R1.GT.R2)  RP=R2
      RT=R2
      IF(R1.GE.R2)  RT=R1
      IF(RP.LE.0.0)  GO  TO  14
      X=R/(RP+RT)
      B=0.00144*ZN1*ZN2/(RP+RT)
      A=RT/RP
      BM=1.+1./A
      C=3./5./A/A
      D=(A+2.+1./A)/4.
      X0=(A-1.)/(A+1.)
      IF(X-X0)  10,10,11
   10 POTNU=B*(3.-C-(BM*X)**2)*BM/2.
      RETURN
   11 IF(X-1.) 12,12,13
   12 POTNU=B*(1.-3.*D*D*((1.-X)**4)*(1.-2./15.*D*(1.-X)*(5.+X)))/X
      RETURN
   13 POTNU=B/X
      RETURN
   14 POTNU=0.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  KINEMQ(PS,V,P,CT,ST,CF,SF,T,CM)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C   KINEMATIC  BLOCK
c
      DIMENSION  PS(3),V(3),P(3)
      PSV=PS(1)*V(1)+PS(2)*V(2)+PS(3)*V(3)
      ES=SQRT(PS(1)**2+PS(2)**2+PS(3)**2+CM**2)
      G=1./SQRT(1.-V(1)**2-V(2)**2-V(3)**2)
      DO  10  K=1,3
   10 P(K)=PS(K)+G*V(K)*(PSV*G/(G+1.)+ES)
      E=G*(ES+PSV)
      T=E-CM
      PM=SQRT(P(1)**2+P(2)**2+P(3)**2)
      CT=P(3)/PM
      ST2=1.-CT*CT
      IF(ST2.LE.0.)   GO  TO  11
      ST=SQRT(ST2)
      CF=P(1)/PM/ST
      SF=P(2)/PM/ST
      RETURN
   11 ST=0.
      CF=1.
      SF=0.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE RXYZ(R12,R0X,R0Y,R0Z)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c    Samples impact parameter and initial positions of
c    projectile and target 
c
      COMMON /XBMAX/ XBMAX,IFIB0
      COMMON /BIMP/ B00,BX,BY
      COMMON /HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      DIMENSION VPR(3),VTA(3)
      IF(IFIB0.EQ.0)  GO TO 10
      B1=RNDM(-1.)
      B2=6.283185*RNDM(-1.)
      B3=R12*SQRT(B1)*XBMAX
      R0X=B3*COS(B2)
      R0Y=B3*SIN(B2)
      GO  TO  11
   10 R0X=R12*XBMAX
      R0Y=0.
      B3=R0X
   11 R0Z=-SQRT(R12**2-R0X**2-R0Y**2)
      B00=B3
      BX=R0X
      BY=R0Y
        CALL  VINIT(VPR,VTA,ANUCL1,ANUCL2,T0)
        GPR=1./SQRT(1.-VPR(3)**2)
        GTA=1./SQRT(1.-VTA(3)**2)
        R0Z=-RM1/GPR-RM2/GTA
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE VINIT(VPR,VTA,A1,A2,T0)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c    Calculates initial velocities of projectile and target 
c
      DIMENSION VPR(3),VTA(3)
      COMMON /KSYST/KSYST
      VPR(1)=0.
      VPR(2)=0.
       VTA(1)=0.
       VTA(2)=0.
      IF(A1.LT.2.1)  GO TO 10
      GO  TO  (10,11,12), KSYST
   10 VPR(3)=SQRT(T0*(T0+1.88))/(T0+0.94)
       VTA(3)=1.E-6
      RETURN
   11 VPR(3)=SQRT(T0*(T0+1.88))/(T0+1.88)
      VTA(3)=-VPR(3)
      RETURN
   12 VPR(3)= SQRT(T0*(T0+1.88))/(T0+0.94*(1.+A1/A2))
      VTA(3)=-SQRT(T0*(T0+1.88))/(T0+0.94*(1.+A2/A1))
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  DISOB(MV,U,IND,NP,MQ)
c new version 11-16-95 09:59am
c     Formation of Delta from pion and nucleon
c
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      REAL*8  MD,M2
      COMMON /PRIMP/ PP(3)
      COMMON /TAUD3/ TAU0
      COMMON /ISOB3/ ISOB3
      COMMON /MEMORY/ PME(9,5999),IME(5,5999)
      IND=0
      IF(ISOB3.EQ.0)  RETURN
      IATT=0
      RND=RNDM(-1.)
      IF(MQ.EQ.2) THEN
        N1=MV+1
        N2=MV+3
      ELSE
        N1=MV+2
        N2=MV+3
      ENDIF
    1 CONTINUE
      IATT=IATT+1
      IF(IATT.GT.2)        RETURN
      IF((RND.LE.0.5.AND.IATT.EQ.1).OR.IATT.EQ.2)  THEN
	       DO  K=1,9
	         TEMP=PME(K,N1)
	         PME(K,N1)=PME(K,N2)
	         PME(K,N2)=TEMP
	       ENDDO
	       DO  K=1,5
	         ITEMP=IME(K,N1)
	         IME(K,N1)=IME(K,N2)
	         IME(K,N2)=ITEMP
	       ENDDO
      ENDIF
      L2=MV+3-MQ
      LD=MV+MQ
      IF((IME(4,LD)+IME(4,MV+3)).NE.1)  GO  TO  1
      P2=SQRT(PME(4,L2)**2+PME(5,L2)**2+PME(6,L2)**2)
      PP(1)=PME(4,MV+2)
      PP(2)=PME(5,MV+2)
      PP(3)=PME(6,MV+2)
      M2=PME(9,L2)
      E2=SQRT(P2**2+M2**2)
      MD=SQRT(U*(U-2.*E2)+M2**2)
      IF(MD.LT.1.082)                   GO  TO  1
      CALL  WMD(MD,T0,FMD)
			IF(ISOB3.EQ.2)                    GO  TO  2
      DRND=RNDM(-1.)
      IF(DRND.GT.FMD)                   GO  TO  1
    2 CONTINUE
      PME(4,LD)  =-PME(4,L2)
      PME(5,LD)  =-PME(5,L2)
      PME(6,LD)  =-PME(6,L2)
      PME(9,LD)  =MD
      IME(1,LD)  = IME(1,LD)  +IME(1,MV+3)
      IME(2,LD)  =0
      IME(3,LD)  =0
      IME(4,LD)  =1
      IME(5,LD)  =INTG(1000.*T0)
      IF(IME(5,LD).LT.1)  IME(5,LD)=1
      NP=NP-1
      IND=MV+2*MQ-1
      DO  11 K=1,9
      PME(K,MV+3)=PME(K,MV+2)
      IF(K.LE.5) IME(K,MV+3)=IME(K,MV+2)
   11 CONTINUE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  WMD(MD,TAU0,FMD)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     sampling the life time of Delta resonance with mass MD
c
      REAL*8 MD,M0
      DATA   M0/1.232D0/
      EN=(MD**2+0.940**2-0.140**2)/(2.*MD)
      Q=SQRT(EN**2-0.940**2)
      R=Q/0.140
      G=0.47*Q*R*R/(1.+0.6*R*R)
      DRND=RNDM(-1.)
      IF(DRND.LT.1.D-10)  DRND=1.D-10
      TAU0=-1./(5.06*G)*LOG(DRND)
      IF(TAU0.LE.0.001)  TAU0=0.0011
      A=G*G/4.
      FMD=A/((MD-M0)**2+A)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DECREN(LR,NUR,MV,NP)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      CHARACTER*8 PNAR,PNAJ
c
c   decays of instable particle LR 
c   Calls: PANUID, DECAYQ 
c
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON /MEMORY/PME(9,5999),IME(5,5999)
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/IDPME/ IDPME(5999)
      COMMON/TLIMIT/TLIMIT
      COMMON/ACTIM/TINT
      COMMON/PORIG/IORI(3,5999)
cc      COMMON/IDpart/IDp1,IDp2
      DIMENSION PIN(9),IIN(5)
C
      NDECAY=1
      IF(TINT.GT.TLIMIT)  NDECAY=5
C
      NP=0
      DO   K=1,9
        PIN(K)=PME(K,LR)
        IF(K.LE.5)  IIN(K)=IME(K,LR)
      ENDDO
      IDR=IDPME(LR)
      IDp1=IORI(1,LR)
      IDp2=IORI(2,LR)
      CALL  PANUID(IDR,NUR,PNAR)
      IF(NCAS.GE.NCPRI) THEN
        write(16,18) PNAR,(PME(I,LR),I=4,9),(IME(K,LR),K=1,5)
        write( *,18) PNAR,(PME(I,LR),I=4,9),(IME(K,LR),K=1,5)
   18 FORMAT(1X,'DECAY:',A5,6(1X,F7.3),4I2,I10)
      ENDIF
      NPTCL=1
      PPTCL(1,1)=PIN(4)
      PPTCL(2,1)=PIN(5)
      PPTCL(3,1)=PIN(6)
      PPTCL(4,1)=PIN(8)+PIN(9)
      PPTCL(5,1)=PIN(9)
      IDENT(1)=IDR
      IORIG(1)=0
      IDCAY(1)=0
C
      do  ND=1,NDECAY
        IF(NPTCL.LE.0)   GO  TO  116
        NDEC=0
        NPTCL1=NPTCL
        DO 114  I=1,NPTCL1
        IDI=IDENT(I)
C       IF(IDI.EQ.110.OR.IDI.EQ.220.OR.IDI.EQ.230.OR.IDI.EQ.-230)
        IF(IDI.EQ.110.OR.IDI.EQ.230.OR.IDI.EQ.-230)
     *                   GO  TO  114
        IF(IDI.EQ.1230.OR.IDI.EQ.-1230)
     *                   GO  TO  114
        CALL DECAYQ(I,IPOINT)
        IF(IPOINT.LT.0)  GO  TO  114
        NDEC=NDEC+1
        DO  J=1,9
          PPTCL(J,I)=PPTCL(J,NPTCL)
        ENDDO
        IDENT(I)=IDENT(NPTCL)
        IORIG(I)=IORIG(NPTCL)
        IDCAY(I)=IDCAY(NPTCL)
        NPTCL=NPTCL-1
  114   CONTINUE
        IF(NDEC.EQ.0)    GO  TO  116
      enddo
  116 CONTINUE
      IF(NPTCL.LE.1)  THEN
        write(16,21) NPTCL,IDI,IPOINT,PIN(9)
        write( *,21) NPTCL,IDI,IPOINT,PIN(9)
   21 FORMAT(5X,'DECREN: NPTCL=',I2,2X,'IDI=',I6,2X,'IPOINT=',I3,
     &1x,'m=',1PE11.4)
        RETURN
      ENDIF
      IF((MV+NPTCL).GT.5999)  THEN
        write(16,'(5X,''DECREN : MV+NPTCL>5999'')')
        RETURN
      ENDIF
      DO  J=1,NPTCL
        M=MV+J
        IDJ=IDENT(J)
        CALL  PANUID(IDJ,JP,PNAJ)
        IME(1,M)=IDINT(1.001*CHARGE(IDJ))
        IME(2,M)=0
C  !!!
        IF(JP.GE.55.AND.JP.LE.65) IME(2,M)=IDJ
        IME(3,M)=IS(IDJ)
        IME(4,M)=IB(IDJ)
        IME(5,M)=0
c  kkg 01.02.06
        if((JP.ge.8.and.JP.le.18).or.(JP.ge.27.and.JP.le.36).or.
     &     (JP.ge.45.and.JP.le.53)) then
          TAU0= TAUN(JP)
          TAUL=TAU0*PPTCL(4,J)/PPTCL(5,J)
          IME(5,M)=INTG(TAUL*1000.)
          IF(IME(5,M).EQ.0)  IME(5,M)=1
        endif 
c          
        PME(4,M)=PPTCL(1,J)
        PME(5,M)=PPTCL(2,J)
        PME(6,M)=PPTCL(3,J)
        PME(8,M)=PPTCL(4,J)-PPTCL(5,J)
        PME(9,M)=PPTCL(5,J)
        IDPME(M)=IDJ
        IF(NCAS.GE.NCPRI)
     *  write(16,'(1X,''====>:'',A5,6(1X,F11.4),4I3,I10)')
     &  PNAJ,(PME(I,M),I=4,9),(IME(K,M),K=1,5)
      ENDDO
      NP=NPTCL
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DECDVM(L,MV)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      REAL*8  MN,MPI,MD
c
c     Decay of Deltas with variable mass MD
c     Calls: KINEMQ
c
      COMMON /MEMORY/ PME(9,5999),IME(5,5999)
      DIMENSION VD(3),PSN(3),PLN(3),IEN(8),IEP(8),CBR(4)
      DATA CBR/1.0,0.6666667,0.6666667,1.0/,
     *IEN/0,0,0,1,1,0,1,1/,
     *IEP/-1,-1,0,-1,0,1,1,1/,MN/0.940/,MPI/0.140/
      MD=PME(9,L)
      ESN=(MD**2+MN**2-MPI**2)/2./MD
      PS=SQRT(ESN**2-MN**2)
      CTS=1.-2.*RNDM(-1.)
      FIS=6.283185*RNDM(-1.)
      STS=SQRT(1.-CTS**2)
      PSN(1)=PS*STS*COS(FIS)
      PSN(2)=PS*STS*SIN(FIS)
      PSN(3)=PS*CTS
      ED=PME(8,L)+MD
      VD(1)=PME(4,L)/ED
      VD(2)=PME(5,L)/ED
      VD(3)=PME(6,L)/ED
      CALL  KINEMQ(PSN,VD,PLN,CTL,STL,CFL,SFL,TL,MN)
      LN=MV+1
      PME(4,LN)=PLN(1)
      PME(5,LN)=PLN(2)
      PME(6,LN)=PLN(3)
      PME(8,LN)=TL
      PME(9,LN)=MN
      LP=MV+2
      PME(4,LP)=PME(4,L)-PLN(1)
      PME(5,LP)=PME(5,L)-PLN(2)
      PME(6,LP)=PME(6,L)-PLN(3)
      PME(8,LP)=PME(8,L)+MD-TL-MN-MPI
      PME(9,LP)=MPI
      ICD=IME(1,L)+2
      IBR=2*(ICD-1)+1
      DRND=RNDM(-1.)
      IF(DRND.GT.CBR(ICD))  IBR=IBR+1
      IME(1,LN)=IEN(IBR)
      IME(2,LN)=0
      IME(3,LN)=0
      IME(4,LN)=1
      IME(5,LN)=0
      IME(1,LP)=IEP(IBR)
      IME(2,LP)=0
      IME(3,LP)=0
      IME(4,LP)=0
      IME(5,LP)=0
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE  PIIZO(U,PL,IFIZO)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
      REAL*8 MN,MPI,MD,MD1,MD2,MD0
c
c    simulate Delta production in N+N==>Delta + N
c
      COMMON/CPIIZO/ PN1(3),PN2(3)
      DIMENSION  PL(3),VD(3),PS(3),PSN(3)
      DATA MN/0.939/,MPI/0.139/,MD0/1.232/
C      write( *,*) 'in piizo'
      IFIZO=0
      MD1=1.081
      IF(U.LE.(MD1+MN)) GO  TO  999
      MD2=U-MN
      IF(MD2.LE.MD1)    GO  TO  999
      FMD2=1.
      IF(MD2.LE.MD0)  THEN
	CALL  WMD(MD2,TAU2,FMD2)
      ENDIF
  901 MD=MD1+RNDM(-1.)*(MD2-MD1)
      CALL  WMD(MD,TAU,FMD)
      IF(RNDM(-1.).GT.FMD/FMD2)  GO  TO  901
      ED=(U**2+MD**2-MN**2)/2./U
      IF(ED.LE.MD)   THEN
       write( *,*) 'ED<MD', ED,MD
       GO  TO  901
      ENDIF
      PD=SQRT(ED**2-MD**2)
c  Parametrization from PR C52(1995)2037 for angular distribution
c  of the form (1/s)ds/do=B1+3*B3*cos(theta)**2
	IF(U.LE.2.14)     THEN
	   B1= 0.5
	   B3= 0.
	ELSEIF(U.LE.2.4)  THEN
	   B1= 29.03-23.75*U+4.865*U**2
	   B3=-30.33+25.53*U-5.301*U**2
	ELSE
	   B1= 0.06
	   B3= 0.4
	ENDIF
      IF(RNDM(-1.).LE.B1/(B1+B3))  THEN
         CTD=1.-2.*RNDM(-1.)
      ELSE
         RND=RNDM(-1.)
         CTDM=(ABS(2.*RND-1.))**(1./3.)
         IF(RND.LE.0.5)  CTD=-CTDM
         IF(RND.GT.0.5)  CTD= CTDM
	ENDIF
      FID=6.283185*RNDM(-1.)
      STD=SQRT(1.-CTD**2)
      VD(1)=PD/ED*STD*COS(FID)
      VD(2)=PD/ED*STD*SIN(FID)
      VD(3)=PD/ED*CTD
      DO  k=1,3
        PN1(k)=-VD(k)*ED
      ENDDO
      ESPI=(MD**2+MPI**2-MN**2)/2./MD
      IF(ESPI.LE.MPI)   THEN
       write( *,*) 'ESPI<MPI', ESPI,MPI
       GO  TO  901
      ENDIF
      PSM=SQRT(ESPI**2-MPI**2)
      CT=1.-2.*RNDM(-1.)
      FI=6.283185*RNDM(-1.)
      ST=SQRT(1.-CT**2)
      PS(1)=PSM*ST*COS(FI)
      PS(2)=PSM*ST*SIN(FI)
      PS(3)=PSM*CT
      CALL  KINEMQ(PS,VD,PL,CT1,ST1,CF1,SF1,TL1,MPI)
      IFIZO=1
      DO  k=1,3
         PSN(k)=-PS(k)
      ENDDO
      CALL  KINEMQ(PSN,VD,PN2,CTN,STN,CFN,SFN,TLN,MN)
  999 CONTINUE
C      write( *,*) 'out piizo'
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE SARNDT(U,SIJ)
      IMPLICIT REAL*8(A-H,O-Z)
c
c   Calculates Arndt's N+N inelastic cross sections 
c
      REAL*8 MPI,MN,M0,MDEL0,MNUC0,MSDEL,MSN
      DIMENSION SIJ(4),ALFA(4),BETA(4),M0(4),GG(4)
      DATA ALFA /3.772,15.28,146.3,6.030/
      DATA BETA/1.262,0.,0.,1.700/
      DATA M0/1188.,1245.,1472.,1203./
      DATA GG/99.02,137.4,26.49,134.3/
C      write( *,*) 'in sarndt'
      do  i=1,4
        SIJ(i)=0.
      enddo
      S=(U*1000.)**2
      PI=3.1415927
      HC=198.
      PHC=PI*HC**2
      MPI=139.0
      MN=939.0
      MNUC0=1430.
      MDEL0=1220.
      G0D=120.
      G0N=200.
      DS=SQRT(S)
      IF(DS.LE.(2.*MN+MPI))  RETURN
      P2=S/4.-MN**2
      ZPD=(DS-MN-MDEL0)*(2./G0D)
      ZMD=(MN+MPI-MDEL0)*(2./G0D)
      MSDEL=MDEL0+G0D*LOG((1.+ZPD**2)/(1.+ZMD**2))/
     *(4.*(ATAN(ZPD)-ATAN(ZMD)))
      ZPN=(DS-MN-MNUC0)*(2./G0N)
      ZMN=(MN+MPI-MNUC0)*(2./G0N)
      MSN=MNUC0+G0N*LOG((1.+ZPN**2)/(1.+ZMN**2))/
     *(4.*(ATAN(ZPN)-ATAN(ZMN)))
      DO 1 I=1,3
      OM=MSDEL
      IF (I.EQ.3) OM=MSN
C
      IF(DS.LE.(MN+OM))    GO  TO  1
      IF(OM.LE.(MN+MPI))  GO  TO  1
C
      GM2=(M0(I)*GG(I))**2
      PR2=(S-(MN-OM)**2)*(S-(MN+OM)**2)/(4.*S)
      Q2S=(OM**2-(MN-MPI)**2)*(OM**2-(MN+MPI)**2)/(4.*OM**2)
      OP=M0(I)
      Q20=(OP**2-(MN-MPI)**2)*(OP**2-(MN+MPI)**2)/(4.*OP**2)
      S0=(MN+M0(I))**2
      P0=S0/4.-MN**2
      AL=ALFA(I)
      BE=BETA(I)
      SIJ(I)=PHC/(2.*P2)*AL*(PR2/P0)**(BE/2.)*GM2*(Q2S/Q20)**(3./2.)/
     *((OM**2-M0(I)**2)**2+GM2)*10.
   1  CONTINUE
c s(NN=D+pi)
      SPIN=(DS-MN)**2
      PR2=(S-(2.*MN-MPI)**2)*(S-(2.*MN+MPI)**2)/(4.*S)
      S0=(MN+M0(4))**2
      P02=S0/4.-MN**2
      GM2=(M0(4)*GG(4))**2
      AL=ALFA(4)
      BE=BETA(4)
      SIJ(4)=PHC/(2.*P2)*AL*(PR2/P02)**(BE/2.)*GM2/
     *((SPIN-M0(4)**2)**2+GM2)*10.
C      write( *,*) 'out sarndt sij=',SIJ
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      DOUBLE PRECISION FUNCTION  SZWER(JP,U)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
c
C    B+B==>B+Y+K   CROSS SECTION  (ZWERMANN APPROXIMATION)
c
      DIMENSION AM(8)
      DATA  AM/2*0.494,1.115,3*1.192,0.494,0.939/
      AMJ=AM(JP)
      AMK=0.494
      SZWER=0.
      IF(JP.LE.2) AMX=1.115+0.939
      IF(JP.GT.2) AMX=0.494+0.939
      IF(U.LE.(AMX+AMJ))  RETURN
      PMAX=SQRT((U**2-(AMX+AMJ)**2)*(U**2-(AMX-AMJ)**2))/(2.*U)
      SZWER=0.049*(PMAX/AMK)**4
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      DOUBLE PRECISION FUNCTION  SST(IB1,IB2,IE1,IE2,JP,T,IS)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
c
C TOTAL PRODUCTION CROSS SECTION OF strange particle JP
c
      SST=0.
      GO  TO  (101,102,103,104,105,106,107),JP
c   for K+      ***********************************
  101 IF((IB1+IB2).EQ.1)  GO  TO  51
      IF((IE1+IE2).EQ.2)  GO  TO  31
      IF((IE1+IE2).EQ.1)  GO  TO  32
      IF((IE1+IE2).EQ.0)  SST=SIG(T,13)
      GO  TO  200
   31 S1=SIG(T,1)
      S2=SIG(T,5)
      S3=SIG(T,7)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2+S3)
                                      SST=S1
        IF(RS.GT.S1.AND.RS.LE.(S1+S2))  SST=S2
        IF(RS.GT.(S1+S2))               SST=S3
      ELSE
C !!!
        SST=S1+S2+S3
C !!!
      ENDIF
      GO  TO  200
   32 S1=SIG(T,2)
      S2=SIG(T,9)
      S3=SIG(T,11)
      IF(IS.EQ.1) THEN

        RS=RNDM(-1.)*(S1+S2+S3)
                                      SST=S1
        IF(RS.GT.S1.AND.RS.LE.(S1+S2))  SST=S2
        IF(RS.GT.(S1+S2))               SST=S3
      ELSE
C !!!
        SST=S1+S2+S3
C !!!
      ENDIF
      GO  TO  200
   51 IF(IE1.EQ.1.AND.IE2.EQ.1)   SST=SIG(T,15)
      IF(IE1.EQ.-1.AND.IE2.EQ.1)  SST=SIG(T,21)
      IF(IE1.EQ.0.AND.IE2.EQ.1)   GO  TO  21
      IF(IE1.EQ.1.AND.IE2.EQ.0)   GO  TO  22
      IF(IE1.EQ.0.AND.IE2.EQ.0)   SST=SIG(T,27)
      GO  TO  200
   21 S1=SIG(T,23)
      S2=SIG(T,25)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2)
                    SST=S1
        IF(RS.GT.S1)  SST=S2
      ELSE
C !!!
        SST=S1+S2
C !!!
      ENDIF
      GO  TO  200
   22 S1=SIG(T,17)
      S2=SIG(T,18)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2)
                    SST=S1
        IF(RS.GT.S1)  SST=S2
      ELSE
C !!!
        SST=S1+S2
C !!!
      ENDIF
      GO  TO  200
c  for K0   **************************************
  102 IF((IB1+IB2).EQ.1)   GO  TO  52
      IF((IE1+IE2).EQ.2)   SST=SIG(T,6)
      IF((IE1+IE2).EQ.1)   GO  TO  33
      IF((IE1+IE2).EQ.0)   GO  TO  34
      GO  TO  200
   33 S1=SIG(T,3)
      S2=SIG(T,8)
      S3=SIG(T,10)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2+S3)
                                      SST=S1
        IF(RS.GT.S1.AND.RS.LE.(S1+S2))  SST=S2
        IF(RS.GT.(S1+S2))               SST=S3
      ELSE
C !!!
        SST=S1+S2+S3
C !!!
      ENDIF
      GO  TO  200
   34 S1=SIG(T,4)
      S2=SIG(T,14)
      S3=SIG(T,12)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2+S3)
                                      SST=S1
        IF(RS.GT.S1.AND.RS.LE.(S1+S2))  SST=S2
        IF(RS.GT.(S1+S2))               SST=S3
      ELSE
C !!!
        SST=S1+S2+S3
C !!!
      ENDIF
      GO  TO  200
   52 IF(IE1.EQ.-1.AND.IE2.EQ.1)  GO  TO  23
      IF(IE1.EQ.0.AND.IE2.EQ.1)   SST=SIG(T,24)
      IF(IE1.EQ.1.AND.IE2.EQ.0)   SST=SIG(T,16)
      IF(IE1.EQ.-1.AND.IE2.EQ.0)  SST=SIG(T,22)
      IF(IE1.EQ.0.AND.IE2.EQ.0)   GO  TO  24
      GO  TO  200
   23 S1=SIG(T,19)
      S2=SIG(T,20)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2)
                    SST=S1
        IF(RS.GT.S1)  SST=S2
      ELSE
C !!!
        SST=S1+S2
C !!!
      ENDIF
      GO  TO  200
   24 S1=SIG(T,28)
      S2=SIG(T,26)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2)
                    SST=S1
        IF(RS.GT.S1)  SST=S2
      ELSE
C !!!
        SST=S1+S2
C !!!
      ENDIF
      GO  TO  200
c  for Lambda  **********************
  103 IF((IB1+IB2).EQ.1)   GO  TO  53
      IF((IE1+IE2).EQ.2)   SST=SIG(T,1)
      IF((IE1+IE2).EQ.1)   GO  TO  25
      IF((IE1+IE2).EQ.0)   SST=SIG(T,4)
      GO  TO  200
   25 S1=SIG(T,2)
      S2=SIG(T,3)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2)
                    SST=S1
        IF(RS.GT.S1)  SST=S2
      ELSE
C !!!
        SST=S1+S2
C !!!
      ENDIF
      GO  TO  200
   53 IF(IE1.EQ.-1.AND.IE2.EQ.1)  SST=SIG(T,19)
      IF(IE1.EQ.0.AND.IE2.EQ.1)   SST=SIG(T,23)
      IF(IE1.EQ.1.AND.IE2.EQ.0)   SST=SIG(T,18)
      IF(IE1.EQ.0.AND.IE2.EQ.0)   SST=SIG(T,26)
      GO  TO  200
c   for S+   ******************************************
  104 IF((IB1+IB2).EQ.1)   GO  TO  54
      IF((IE1+IE2).EQ.2)   GO  TO  26
      IF((IE1+IE2).EQ.1)   SST=SIG(T,8)
      GO  TO  200
   26 S1=SIG(T,5)
      S2=SIG(T,6)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2)
                    SST=S1
        IF(RS.GT.S1)  SST=S2
      ELSE
C !!!
        SST=S1+S2
C !!!
      ENDIF
      GO  TO  200
   54 IF(IE1.EQ.1.AND.IE2.EQ.1)   SST=SIG(T,15)
      IF(IE1.EQ.0.AND.IE2.EQ.1)   SST=SIG(T,24)
      IF(IE1.EQ.1.AND.IE2.EQ.0)   SST=SIG(T,16)
      GO  TO  200
c  for S-   *****************************************
  105 IF((IB1+IB2).EQ.1)   GO  TO  55
      IF((IE1+IE2).EQ.1)   SST=SIG(T,11)
      IF((IE1+IE2).EQ.0)   GO  TO  27
      GO  TO  200
   27 S1=SIG(T,13)
      S2=SIG(T,14)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2)
                    SST=S1
        IF(RS.GT.S1)  SST=S2
      ELSE
C !!!
        SST=S1+S2
C !!!
      ENDIF
      GO  TO  200
   55 IF(IE1.EQ.-1.AND.IE2.EQ.1)  SST=SIG(T,21)
      IF(IE1.EQ.-1.AND.IE2.EQ.0)  SST=SIG(T,22)
      IF(IE1.EQ.0.AND.IE2.EQ.0)   SST=SIG(T,27)
      GO  TO  200
c    for    S0   *************************************
  106 IF((IB1+IB2).EQ.1)   GO  TO  56
      IF((IE1+IE2).EQ.2)   SST=SIG(T,7)
      IF((IE1+IE2).EQ.1)   GO  TO  28
      IF((IE1+IE2).EQ.0)   SST=SIG(T,12)
      GO  TO  200
   28 S1=SIG(T,9)
      S2=SIG(T,10)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2)
                    SST=S1
        IF(RS.GT.S1)  SST=S2
      ELSE
C !!!
        SST=S1+S2
C !!!
      ENDIF
      GO  TO  200
   56 IF(IE1.EQ.-1.AND.IE2.EQ.1)  SST=SIG(T,20)
      IF(IE1.EQ.0.AND.IE2.EQ.1)   SST=SIG(T,25)
      IF(IE1.EQ.1.AND.IE2.EQ.0)   SST=SIG(T,17)
      IF(IE1.EQ.0.AND.IE2.EQ.0)   SST=SIG(T,28)
      GO  TO  200
  107 CONTINUE
  200 CONTINUE
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      DOUBLE PRECISION FUNCTION SSL(IB1,IB2,IE1,IE2,JP,T)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
c
C TOTAL PRODUCTION CROSS SECTION OF JP+L(AS A PARTENER)
c
      SSL=0.
      GO  TO  (201,202,103,300,300,300,300),JP
  103 SSL=SST(IB1,IB2,IE1,IE2,JP,T,0)
      GO  TO  300
  201 IF((IB1+IB2).EQ.1)  GO  TO  61
      IF((IE1+IE2).EQ.2)   SSL=SIG(T,1)
      IF((IE1+IE2).EQ.1)   SSL=SIG(T,2)
      GO  TO  300
   61 IF(IE1.EQ.0.AND.IE2.EQ.1)   SSL=SIG(T,23)
      IF(IE1.EQ.1.AND.IE2.EQ.0)   SSL=SIG(T,18)
      GO  TO  300
  202 IF((IB1+IB2).EQ.1)   GO  TO  62
      IF((IE1+IE2).EQ.1)   SSL=SIG(T,3)
      IF((IE1+IE2).EQ.0)   SSL=SIG(T,4)
      GO  TO  300
   62 IF(IE1.EQ.-1.AND.IE2.EQ.1)  SSL=SIG(T,19)
      IF(IE1.EQ.0.AND.IE2.EQ.0)   SSL=SIG(T,26)
  300 CONTINUE
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      DOUBLE PRECISION FUNCTION  SIG(T,N)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
c
C                        N                  REACTION
C			 1*             p+p==>p + L  + K+
C			 2              p+n==>n + L  + K+
C			 3              p+n==>p + L  + K0
C			 4              n+n==>n + L  + K0
C			 5*             p+p==>n + S+ + K+
C			 6*             p+p==>p + S+ + K0
C			 7*             p+p==>p + S0 + K+
C			 8              p+n==>n + S+ + K0
C			 9              p+n==>n + S0 + K+
C			10              p+n==>p + S0 + K0
C			11              p+n==>p + S- + K+
C			12              n+n==>n + S0 + K0
C			13              n+n==>n + S- + K+
C			14              n+n==>p + S- + K0
C			15*             (pi+)+p==>S+ + K+
C			16              (pi+)+n==>S+ + K0
C			17              (pi+)+n==>S0 + K+
C			18              (pi+)+n==>L  + K+
C			19*             (pi-)+p==>L  + K0
C			20*             (pi-)+p==>S0 + K0
C			21*             (pi-)+p==>S- + K+
C			22              (pi-)+n==>S- + K0
C			23              (pi0)+p==>L  + K+
C			24              (pi0)+p==>S+ + K0
C			25              (pi0)+p==>S0 + K+
C			26              (pi0)+n==>L  + K0
C			27              (pi0)+n==>S- + K+
C			28              (pi0)+n==>S0 + K0
C                        *    -          basic channel
C   PP   --> L K+ P
      SN1(X)=FAR(0.132D0,0.709D0,3.432D0,-0.564D0,X)
C   PP   --> S+ K+ N
      SN5(X)=FAR(0.144D0,0.709D0,3.432D0,-0.564D0,X)
C   PP   --> S+ K0 P
      SN6(X)=FAR(0.062D0,0.709D0,3.432D0,-0.564D0,X)
C   PP   --> S0 K+ P
      SN7(X)=FAR(0.058D0,0.709D0,3.432D0,-0.564D0,X)
C   PI+P --> K+ S+
      S15(X)=FAR(0.565D0,1.129D0,0.662D0,-1.539D0,X)
C   PI-P --> K0 L
      S19(X)=FAR(0.183D0,1.375D0,0.130D0,-1.213D0,X)
C   PI-P --> K0 S0
      S20(X)=FAR(0.098D0,1.221D0,0.150D0,-1.073D0,X)
C   PI-P --> K+ S-
      S21(X)=FAR(0.112D0,0.873D0,0.457D0,-1.724D0,X)
C
      TL2=T-0.759
      TS2=T-0.888
      TL3=T-1.579
      TS3=T-1.780

c      TL2=T-0.760
c      TS2=T-0.900
c      TL3=T-1.570
c      TS3=T-1.800
      GO  TO  (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     *21,22,23,24,25,26,27,28),N
    1 SIG=SN1(TL3)
      RETURN
    2 SIG=5.*SN1(TL3)
      RETURN
    3 GO  TO  2
    4 SIG=1.*SN1(TL3)
      RETURN
    5 SIG=SN5(TS3)
      RETURN
    6 SIG=SN6(TS3)
      RETURN
    7 SIG=SN7(TS3)
      RETURN
    8 SIG=1./3.*SN5(TS3)-2.*SN7(TS3)+4.*SN6(TS3)
      RETURN
    9 SIG=4./3.*SN5(TS3)-7.*SN7(TS3)+6.*SN6(TS3)
      RETURN
   10 GO  TO  9
   11 SIG=9.*SN5(TS3)-8.*SN7(TS3)+9.*SN6(TS3)
      RETURN
   12 GO  TO  7
   13 SIG=2./3.*SN5(TS3)-4.*SN7(TS3)+3.*SN6(TS3)
   14 GO  TO  5
   15 SIG=S15(TS2)
      RETURN
   16 GO  TO  21
   17 GO  TO  20
   18 GO  TO  19
   19 SIG=S19(TL2)
      RETURN
   20 SIG=S20(TS2)
      RETURN
   21 SIG=S21(TS2)
      RETURN
   22 GO  TO  15
   23 SIG=S19(TL2)/2.
      RETURN
   24 SIG=2./3.*(S21(TS2)+S15(TS2)/3.-S20(TS2)/2.)
      RETURN
   25 SIG=1./2.*(S21(TS2)+S15(TS2)-S20(TS2))
      RETURN
   26 GO  TO  23
   27 GO  TO  20
   28 SIG=1./2.*(S21(TS2)+S15(TS2)-S20(TS2))
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      DOUBLE PRECISION FUNCTION  FAR(A1,A2,A3,A4,X)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
c   Function form of approximated strange particle 
c   production cross sections 
      FAR=0.0
      IF(X.LE.0.)   RETURN
      S=A3**2+X**2
      FAR=A1*(X**A2)*(S**A4)
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      DOUBLE PRECISION FUNCTION SPINET(T)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
      REAL*8  mpi,mn,meta
C      (pi-)+p=>eta+n   CROSS SECTIONS  [mb]
c     T - lab. kinetic energy of pion
c
      DATA mpi/0.140/,mn/0.940/,meta/0.549/
      ss0=mn+meta
      E=T+mpi
      P=SQRT(T*(T+2.*mpi))
      ss=SQRT((E+mn)**2-P**2)
      SPINET=0.
      if(ss.le.ss0) RETURN
      x=ss-ss0
      SPINET=FAR(0.485D0,0.643D0,0.039D0,-0.649D0,x)
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      DOUBLE PRECISION FUNCTION SPPETA0(T,KSI)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
      REAL*8  mn,meta
C      p+p=>eta+p+p   CROSS SECTIONS  [mb]
c  De Paoli et al. Phys.Lett. B219(1989),194 [A=0.167,B=0.253]
c  New parametrization (Gudima,1996)
c     T - lab. kinetic energy of proton
      DATA mn/0.939/,meta/0.549/,A/0.723/,B/0.411/
      ss0=2.*mn+meta
      E=T+mn
      P=SQRT(T*(T+2.*mn))
      ss=SQRT((E+mn)**2-P**2)
      SPPETA0=0.
      if(ss.le.ss0) RETURN
      x=ss-ss0
      SPPETA0=A*x/(B+x**2)
c  !!!   s(p+n=p+n+eta)=2*s(p+p=p+p+eta)
      IF(KSI.eq.2)  SPPETA0=2.*SPPETA0
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      DOUBLE PRECISION FUNCTION SPPETA(T,KSI)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
      REAL*8  mn,meta
C      p+p=>eta+p+p   CROSS SECTIONS  [mb]
c  Parametrization of M. Thomere et al nucl-th/0702004v1
c  sig=a*x**b, sig(mkb), x=ss-ss0 (GeV)
c     T - lab. kinetic energy of proton
      DATA mn/0.939/,meta/0.549/,A/0.723/,B/0.411/
      ss0=2.*mn+meta
      E=T+mn
      P=SQRT(T*(T+2.*mn))
      ss=SQRT((E+mn)**2-P**2)
      SPPETA=0.
      if(ss.le.ss0) RETURN
      x=ss-ss0
      xm=x*1000.
      if(KSI.eq.1)  then 
         if(xm.lt.283.0)  then
            aw=1213.8
            bw=1.50
         elseif(xm.lt.651.)  then
            aw=162.1
            bw=-0.08
         else
            aw=99.6
            bw=-1.24
         endif
         SPPETA=aw*(x**bw)/1000. ! in mb 
      else
         if(xm.lt.200.0)  then
            aw=25623.0
            bw=2.03
         elseif(xm.lt.651.)  then
            aw=324.3
            bw=-0.08
         else
            aw=199.0
            bw=-1.24
         endif
         SPPETA=aw*(x**bw)/1000. ! in mb 
      endif
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      DOUBLE PRECISION FUNCTION SPIDPP(Tpi)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
      REAL*8 mpi,md
c  cross section for pi+d=p+p; Tpi is pion lab.kin.ener.
c  Parametrization B.G.Ritche,Phys.Rev. C28(1983)926
      DATA a/-1.2/,b/3.5/,c/7.4/,d/5600./,ER/2136./
      DATA  mpi/139.0/,md/1878./
      TpiM=Tpi*1000.
      E=SQRT((mpi+md)**2+2.*TpiM*md)
      SPIDPP=a+b/SQRT(TpiM)+c*10000./((E-ER)**2+d)
	IF(SPIDPP.LE.0.)  SPIDPP=0.
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      DOUBLE PRECISION FUNCTION SSDT(IB1,IB2,IE1,IE2,JP,T,IS)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
c
C  PRODUCTION CROSS SECTION N+N=>DELTA+JP+X
c
      SSDT=0.
      IF((IB1+IB2).NE.2)   RETURN
      IF(IE1+IE2-1)  100,101,102
  100 GO  TO  (11,12,13,14,15,16,17),JP
   11 S1=SIGD(T,5)
      S2=SIGD(T,17)
      S3=SIGD(T,20)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2+S3)
                                      SSDT=S1
        IF(RS.GT.S1.AND.RS.LE.(S1+S2))  SSDT=S2
        IF(RS.GT.(S1+S2))               SSDT=S3
      ELSE
C !!!
        SSDT=S1+S2+S3
C !!!
      ENDIF
      RETURN
   12 S1=SIGD(T,6)
      S2=SIGD(T,16)
      S3=SIGD(T,18)
      S4=SIGD(T,19)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2+S3+S4)
                                              SSDT=S1
        IF(RS.GT.S1.AND.RS.LE.(S1+S2))          SSDT=S2
        IF(RS.GT.(S1+S2).AND.RS.LE.(S1+S2+S3))  SSDT=S3
        IF(RS.GT.(S1+S2+S3))                    SSDT=S4
      ELSE
C !!!
        SSDT=S1+S2+S3+S4
C !!!
      ENDIF
      RETURN
   13 S1=SIGD(T,5)
      S2=SIGD(T,6)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2)
                    SSDT=S1
        IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
        SSDT=S1+S2
C !!!
      ENDIF
      RETURN
   14 SSDT=SIGD(T,19)
      RETURN
   15 S1=SIGD(T,16)
      S2=SIGD(T,17)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2)
                    SSDT=S1
        IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
        SSDT=S1+S2
C !!!
      ENDIF
      RETURN
   16 S1=SIGD(T,18)
      S2=SIGD(T,20)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2)
                    SSDT=S1
        IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
        SSDT=S1+S2
C !!!
      ENDIF
   17 CONTINUE
      RETURN
  101 GO  TO  (21,22,23,24,25,26,27),JP
   21 S1=SIGD(T,4)
      S2=SIGD(T,14)
      S3=SIGD(T,15)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2+S3)
                                      SSDT=S1
        IF(RS.GT.S1.AND.RS.LE.(S1+S2))  SSDT=S2
        IF(RS.GT.(S1+S2))               SSDT=S3
      ELSE
C !!!
        SSDT=S1+S2+S3
C !!!
      ENDIF
      RETURN
   22 S1=SIGD(T,11)
      S2=SIGD(T,12)
      S3=SIGD(T,13)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2+S3)
                                      SSDT=S1
        IF(RS.GT.S1.AND.RS.LE.(S1+S2))  SSDT=S2
        IF(RS.GT.(S1+S2))               SSDT=S3
      ELSE
C !!!
        SSDT=S1+S2+S3
C !!!
      ENDIF
      RETURN
   23 S1=SIGD(T,3)
      S2=SIGD(T,4)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2)
                    SSDT=S1
        IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
        SSDT=S1+S2
C !!!
      ENDIF
      RETURN
   24 S1=SIGD(T,13)
      S2=SIGD(T,14)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2)
                    SSDT=S1
        IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
        SSDT=S1+S2
C !!!
      ENDIF
      RETURN
   25 SSDT=SIGD(T,11)
      RETURN
   26 S1=SIGD(T,12)
      S2=SIGD(T,15)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2)
                    SSDT=S1
        IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
        SSDT=S1+S2
C !!!
      ENDIF
   27 RETURN
  102 GO  TO  (31,32,33,34,35,36,37),JP
   31 S1=SIGD(T,2)
      S2=SIGD(T,8)
      S3=SIGD(T,9)
      S4=SIGD(T,10)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2+S3+S4)
                                              SSDT=S1
        IF(RS.GT.S1.AND.RS.LE.(S1+S2))          SSDT=S2
        IF(RS.GT.(S1+S2).AND.RS.LE.(S1+S2+S3))  SSDT=S3
        IF(RS.GT.(S1+S2+S3))                    SSDT=S4
      ELSE
C !!!
        SSDT=S1+S2+S3+S4
C !!!
      ENDIF
      RETURN
   32 S1=SIGD(T,1)
      S2=SIGD(T,7)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2)
                    SSDT=S1
        IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
        SSDT=S1+S2
C !!!
      ENDIF
      RETURN
   33 S1=SIGD(T,1)
      S2=SIGD(T,2)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2)
                    SSDT=S1
        IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
        SSDT=S1+S2
C !!!
      ENDIF
      RETURN
   34 S1=SIGD(T,7)
      S2=SIGD(T,9)
      IF(IS.EQ.1) THEN
        RS=RNDM(-1.)*(S1+S2)
                    SSDT=S1
        IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
        SSDT=S1+S2
C !!!
      ENDIF
      RETURN
   35 SSDT=SIGD(T,10)
      RETURN
   36 SSDT=SIGD(T,8)
   37 RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      DOUBLE PRECISION FUNCTION SSDL(IB1,IB2,IE1,IE2,JP,T)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
c
C  PRODUCTION CROSS SECTION N+N=>L+DELTA+K
c
      SSDL=0.
      IF((IB1+IB2).NE.2)  RETURN
      IF(IE1+IE2-1)  100,101,102
  100 GO  TO  (11,12,13,14,14,14,14),JP
   11 SSDL=SIGD(T,5)
      RETURN
   12 SSDL=SIGD(T,6)
      RETURN
   13 SSDL=SSDT(IB1,IB2,IE1,IE2,JP,T,0)
   14 RETURN
  101 GO  TO  (21,22,13,14,14,14,14),JP
   21 SSDL=SIGD(T,4)
      RETURN
   22 SSDL=SIGD(T,3)
      RETURN
  102 GO  TO  (31,32,13,14,14,14,14),JP
   31 SSDL=SIGD(T,2)
      RETURN
   32 SSDL=SIGD(T,1)
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      DOUBLE PRECISION FUNCTION  SIGD(T,N)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
c
C
C   PP  --->  D+ L K+
      SD2(X)=FAR(0.132D0,0.709D0,3.432D0,-0.564D0,X)
C   PP  --->  D+ S+ K0
      SD7(X)=FAR(0.062D0,0.709D0,3.432D0,-0.564D0,X)
C   PP  --->  D+ S0 K+
      SD8(X)=FAR(0.058D0,0.709D0,3.432D0,-0.564D0,X)
C   PP  --->  D0 S+ K+
      SD9(X)=0.25*FAR(0.144D0,0.709D0,3.432D0,-0.564D0,X)
C
      TDL=T-2.432
      TDS=T-2.870
      GO  TO  (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),N
    1 SIGD=3.*SD2(TDL)
      RETURN
    2 SIGD=SD2(TDL)
      RETURN
    3 SIGD=2.*SD2(TDL)
      RETURN
    4 SIGD=2.*SD2(TDL)
      RETURN
    5 SIGD=3.*SD2(TDL)
      RETURN
    6 SIGD=SD2(TDL)
      RETURN
    7 SIGD=SD7(TDS)
      RETURN
    8 SIGD=SD8(TDS)
      RETURN
    9 SIGD=SD9(TDS)
      RETURN
   10 SIGD=1./2.*(9.*SD7(TDS)-6.*SD8(TDS)+2.*SD9(TDS))
      RETURN
   11 SIGD=3.*SD9(TDS)
      RETURN
   12 SIGD=1./6.*(9.*SD7(TDS)-6.*SD8(TDS)+8.*SD9(TDS))
      RETURN
   13 SIGD=1./6.*(15.*SD7(TDS)-6.*SD8(TDS)+2.*SD9(TDS))
      RETURN
   14 GO  TO  11
   15 GO  TO  12
   16 GO  TO  9
   17 SIGD=1./3.*(9.*SD7(TDS)-12.*SD8(TDS)+8.*SD9(TDS))
      RETURN
   18 GO  TO  8
   19 GO  TO  10
   20 SIGD=1./2.*(9.*SD7(TDS)-12.*SD8(TDS)+8.*SD9(TDS))
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      DOUBLE PRECISION FUNCTION  STHE(U,JP,ST4)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
c
C   S19    =  TOTAL CROSS SECTION AT T0=19.0 GeV  JP+N
C   FOR  JP=     K+  K0  L  S+  S-  S0  K-
c
      DIMENSION  S19(7)
      DATA  S19/1.8,1.8,0.9,0.9,0.9,0.9,0.6/
      U4=3.325
      U19=6.13
      STHE=ST4+(S19(JP)-ST4)*LOG(U/U4)/LOG(U19/U4)
      RETURN
      END
C
C * * * * * * * * * *  G: 20.05.93  * * * * * * * * * *
C
      DOUBLE PRECISION FUNCTION  SIKMI(PMAX)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
c
C   PRODUCTION CROSS SECTION IN  NN ---> N N K+ K-
C     SIKMI=PMAX/40.   (mb)
c
      SIKMI=PMAX/40.
      RETURN
      END
C ********************************************************
      SUBROUTINE  Tpiabs
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
c
c  test subroutine for absorption pi+(NN) cross section
c
      REAL*8 mpi
      DATA mpi/0.139/, A/0.230/,B/7.14/
      JP=10
      Tpi0=0.001
      DTpi=0.010
      DO  it=1,35
       Tpi=(it-1)*Dtpi
       IF(it.eq.1) Tpi=Tpi0
       CALL SIGPI(JP,Tpi,Sold,3)
       Snew=SPIDPP(Tpi)
       ppi2=Tpi*(Tpi+2.*mpi)
       ppi=SQRT(ppi2)
c Metropolis parametrization (see Garpman et al. P.L.B86(1979)133
       Smet=10.*A*(mpi**2+B*ppi2)/(ppi*mpi)
       TpiM=Tpi*1000.
       write(16,100) TpiM,Sold,Smet,Snew
  100  format(1x,F5.1,3(1x,1PE11.4))
      ENDDO
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      DOUBLE PRECISION FUNCTION CROSSS(JP1,IKS,A,Z,T)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
c
c     calculates isospin averaged cross section
c
      DIMENSION K1(7),K2(7)
C     JP    1   2   3   4   5   6  7
C          PI- PI0 PI+  P   N   d eta
      DATA K1/2,2,1,1,2,1,1/
      DATA K2/1,3,3,2,1,2,2/
      JP=JP1-1
      IF(RNDM(-1.).LE.(Z/A))  THEN
        KSI=K1(JP)
      ELSE
        KSI=K2(JP)
      ENDIF
      IF(JP.LE.3) THEN
C pi+N cross section
	ST=SIGMAG(0,0,1,KSI,IKS,T)
      ELSEIF(JP.EQ.4.OR.JP.EQ.5) THEN
C N+N cross section
	ST=SIGMAG(0,0,2,KSI,IKS,T)
      ELSEIF(JP.EQ.6)  THEN
C d+N cross section
	ST=SIGMAD(KSI,IKS,T)
      ELSEIF(JP.EQ.7)  THEN
C eta+N cross section
        ST=SIGETA(KSI,IKS,T)
      ELSE
        ST=0.
      ENDIF
      CROSSS=ST
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE SIGPI(JP,T,SIG,IS)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
c
c     Calculates pion absorption cross section on pair (pn)
c
      COMMON/PIDABS/PIDABS
      DIMENSION IEPI(3),KSIP(3),KSIN(3)
      DATA IEPI/-1,0,1/,KSIP/2,2,1/,KSIN/1,3,3/
      JPI=JP-8
      IEP=IEPI(JPI)
      KSI=KSIN(JPI)
      IF(RNDM(-1.).GT.0.5)  KSI=KSIP(JPI)
      SIG=SIGMAG(0,0,1,KSI,IS,T)
c      IF(IS.EQ.3) SIG=SIG*PIDABS
c as compared to free pi+d=pp cr.sec. abs. probab. is higher
c by factor 3, see J.Hufner, Phys.Rep.C21,1(1975))
c and by weight  Z/A*N/A~1/4 of pn pairs: PIDABS=3*1/4=0.75
      IF(IS.EQ.3) SIG=SPIDPP(T)*PIDABS
      RETURN
      END
C ********************************************************
      DOUBLE PRECISION FUNCTION SIGETA(KSI,IKS,T)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
c
c   eta + N cross sections 
c
      REAL*8 mn,mpi,met,MS,k
      DATA  mn/0.939/,mpi/0.139/,met/0.549/,MS/1.535/,
     &      G0/0.150/,Sm/74./,Bpi/0.55/
      P0=SQRT(T*(T+2.*met))
      E0=T+met
      ss=SQRT((E0+mn)**2-P0**2)
      s=ss**2
      Eet=(s+met**2-mn**2)/(2.*ss)
      Epi=(s+mpi**2-mn**2)/(2.*ss)
      q=SQRT(Eet**2-met**2)
      Eet0=(MS**2+met**2-mn**2)/(2.*MS)
      q0=SQRT(Eet0**2-met**2)
      G=G0*(q/q0)
      k=SQRT(Epi**2-mpi**2)
      T0pi=(s-(mpi+mn)**2)/(2.*mn)
      sinv=SPINET(T0pi)
      sdir=3.*sinv*(k/q)**2
      sres=Sm*((q0/q)**2)*(G/2.)**2/((ss-MS)**2+(G/2.)**2)
      IF(IKS.eq.0)      then
c total cr.sec.
c       SIGETA=sres
        SIGETA=sdir/Bpi
      ELSEIF(IKS.eq.1)  then
c elastic cr.sec.
c       SIGETA=0.5*sres
        SIGETA=(1.-Bpi)/Bpi*sdir
      ELSEIF(IKS.eq.2)  then
c eta+N-->pi+N
        SIGETA=sdir
      ELSEIF(IKS.eq.3)  then
c abs. cr. sec.
c       SIGETA=0.5*sres
        SIGETA=sdir
      ELSE
        SIGETA=sres
      ENDIF
      IF(KSI.eq.2)  SIGETA=SIGETA*2.
      RETURN
      END
C ********************************************************
      DOUBLE PRECISION FUNCTION SIGMAD(KSI,IKS,T)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
c
c     Approximated d+N cross sections
c
      TN=T/2.
      IF(IKS.EQ.0)      THEN
C total Nd cross section
	SIGMAD=QINDCR(TN,0)
	RETURN
      ELSEIF(IKS.EQ.1)  THEN
C elastic Nd cross section
	SIGMAD=QINDCR(TN,0)-QINDCR(TN,1)
	RETURN
      ELSEIF(IKS.EQ.3)  THEN
C absorption Nd cross section=inelastic cross section
	SIGMAD=QINDCR(TN,1)
	RETURN
      ELSE
	SIGMAD=0.
	WRITE(16,*) 'SIGMAD:IKS=',IKS
	WRITE( *,*) 'SIGMAD:IKS=',IKS
	RETURN
      ENDIF
c     RETURN
      END
C ********************************************************
      DOUBLE PRECISION FUNCTION QINDCR(TN,IS)
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
c
c   Interpolation of d+N cross sections
c
      DIMENSION TNTO(30),TNIN(30),SNDTO(30),SNDIN(30),XX(30),YY(30)
      DATA TNTO/
     *    .0000,     .0010,     .0020,     .0030,     .0040,
     *    .0050,     .0060,     .0070,     .0080,     .0090,
     *    .0100,     .0200,     .0300,     .0400,     .0500,
     *    .0600,     .0700,     .0800,     .0900,     .1000,
     *    .2000,     .3000,     .4000,     .5000,     .6000,
     *    .7000,     .8000,     .9000,    1.0000,    2.0000/
      DATA TNIN/
     *    .0003339,  .0034,     .0035,     .00375,    .0040,
     *    .0050,     .0060,     .0070,     .0080,     .0090,
     *    .0100,     .0200,     .0300,     .0400,     .0500,
     *    .0600,     .0700,     .0800,     .0900,     .1000,
     *    .2000,     .3000,     .4000,     .5000,     .6000,
     *    .7000,     .8000,     .9000,    1.0000,    2.0000/
      DATA SNDIN/
     *    .0   ,    1.2   ,    2.9   ,    8.0   ,   13.5   ,
     *  37.0   ,   60.2   ,   82.0   ,  102.7   ,  121.0   ,
     * 137.0   ,  205.3   ,  175.0   ,  140.0   ,  105.0   ,
     *  85.0   ,   72.0   ,   61.0   ,   54.0   ,   50.0   ,
     *  30.0   ,   26.0   ,   24.0   ,   23.0   ,   22.0   ,
     *  21.5   ,   21.2   ,   21.0   ,   20.5   ,   18.0   /
      DATA SNDTO/
     *3154.0000, 2889.6817, 2572.7316, 2246.4534, 1890.2873,
     *1688.6910, 1522.4142, 1357.7684, 1241.8348, 1158.5485,
     *1067.9506,  596.6951,  379.1020,  288.4043,  219.2028,
     * 177.6805,  150.7435,  126.9888,  113.7065,  105.0809,
     *  64.7585,   56.1204,   57.3817,   64.8824,   70.7762,
     *  73.8892,   76.2111,   79.0000,   80.0000,   80.0000/
C     IS=0  for total Nd cross section
C     IS=1  for inelastic(n,2n) Nd cross section
      MQ=30
      X=TN
      IF(IS.EQ.0)     THEN
	DO  K=1,30
	  XX(K)=TNTO(K)
	  YY(K)=SNDTO(K)
	ENDDO
      ELSEIF(IS.EQ.1) THEN
	DO  K=1,30
	  XX(K)=TNIN(K)
	  YY(K)=SNDIN(K)
	ENDDO
      ELSE
	QINDCR=0.
	WRITE(16,*) 'QINDCR: IS=',IS
	WRITE( *,*) 'QINDCR: IS=',IS
	RETURN
      ENDIF
      IF(X.LT.XX(1))   THEN
        QINDCR=YY(1)
	RETURN
      ENDIF
      DO  j=1,MQ
       K=j
       IF(ABS(X-XX(j)).LT.1.D-10)  go  to  16
      ENDDO
      K=1
   10 IF(X-XX(K))  11,16,17
   11 IF(K-1)      12,12,13
   12 Y1=YY(1)
      X1=XX(1)
      Y2=YY(2)
      X2=XX(2)
      Y3=YY(3)
      X3=XX(3)
      GO  TO  18
   13 IF(K-(MQ-1)) 15,14,14
   14 Y1=YY(MQ-2)
      X1=XX(MQ-2)
      Y2=YY(MQ-1)
      X2=XX(MQ-1)
      Y3=YY(MQ)
      X3=XX(MQ)
      GO  TO  18
   15 Y1=YY(K-1)
      X1=XX(K-1)
      Y2=YY(K)
      X2=XX(K)
      Y3=YY(K+1)
      X3=XX(K+1)
      GO  TO  18
   16 QINDCR= YY(K)
      RETURN
   17 K = K + 1
      IF(K.GT.MQ)  GO  TO  14
      GO  TO  10
   18 DT =(X2-X3)*X1**2+(X3-X1)*X2**2+(X1-X2)*X3**2
      DA=(X2-X3)*Y1   +(X3-X1)*Y2+(X1-X2)*Y3
      DB=(Y2-Y3)*X1**2+(Y3-Y1)*X2**2+(Y1-Y2)*X3**2
      DG=(X2*Y3-X3*Y2)*X1**2+(X3*Y1-X1*Y3)*X2**2+(X1*Y2-X2*Y1)*X3**2
      AL=DA/DT
      BE=DB/DT
      GA=DG/DT
      QINDCR = AL*X**2 + BE*X + GA
      IF(QINDCR.LT.0.)  QINDCR=0.
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      BLOCK  DATA  DGTAB
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
      COMMON/SABCDE/EG(8),ST(8),A(8),B(8),C(8),D(8),E(8)
c
c Photodesintegration cr.section and coeff. for diff.cr.sec
c from F.Partovi, Ann.Physics 27,79(1964)
c
      DATA  EG/
     *   10.,   20.,    40.,  60.,   80.,   100.,   120.,   140./
      DATA  ST/
     * 1387., 558.2, 224.2, 126.6,  87.4,  66.25,  53.41,  44.53/
      DATA   A/
     * 4.623, 5.387, 6.236, 6.101, 5.651,  5.114,  4.618,  4.156/
      DATA   B/
     * 162.0, 65.50, 19.14, 7.209, 3.009,  1.025, 0.0653,-0.4096/
      DATA   C/
     *0.3744,0.7314,0.9983,1.0010,0.9420, 0.8808, 0.8230, 0.7737/
      DATA   D/
     * 29.96, 18.25, 7.914, 4.452, 2.774,  1.812,  1.237, 0.8757/
      DATA   E/
     *-4.230,-4.223,-2.164,-1.564,-1.317,-0.9853,-0.7715,-0.6367/
      END
C ********************************************************
      BLOCK  DATA  WXW
      IMPLICIT REAL*8  (A-H,O-Z), INTEGER (I-N)
      COMMON/WX/W(8),XW(8)
c
c     Data for Gauss 8-point integration
c
      DATA W/
     * 0.1012285352, 0.2223809958, 0.3137066364, 0.3626837730,
     * 0.3626837730, 0.3137066364, 0.2223809958, 0.1012285352/
      DATA XW/
     *-0.9602898359,-0.7966664433,-0.5255323648,-0.1834346056,
     * 0.1834346056, 0.5255323648, 0.7966664433, 0.9602898359/
      END

c **********************************************************************
      subroutine gengamn(eg,emn,ichtar,mv,np) 
      implicit real*8 (a-h, o-z), integer (i-n)
c
c  This subroutine generates a event gamma + N ==> hadrons
c  eg - energy of gamma, GeV; 
c  ichtar = 1(0) for proton(neutron) target
c  resuring products are in memory arrays pme, ime 
c 
	common /uinit/ u,v(3)
	dimension pin(9),pn(9),iin(5),ipn(5)
      data emneut, emprot, empich, empi0 / 0.9395656d0, 0.9382723d0,
     &                                     0.139568d0, 0.134973d0/
	data zro/0.0d0/
c 
      do  k=1,9
	  pin(k) = zro
	   pn(k) = zro
	  if(k.le.5)  then
	    iin(k) = 0
		ipn(k) = 0
	  endif
	enddo
	pin(6) = eg
	pin(8) = eg
	pn(6)  = 0.001
	pn(4)  = 0.001
	pn(9)  = emn
	if(ichtar.eq.1) then      ! target is proton
	  ipn(1)= 1
	else                      ! target is neutron
	  ipn(1)=0     	    
      endif
	ipn(4) = 1
	pn(8) = sqrt(pn(4)**2+pn(5)**2+pn(6)**2+pn(9)**2)-pn(9)
	call tinvuq(pin,pn,u,v,tin1)
	call  gntoh(v,u,tin1,pin,iin,pn,ipn,mv,np)
	return
	end      
c  **************************************************************    
      subroutine gnappr(eg,sig)
      implicit real*8 (a-h, o-z), integer (i-n)
c
C    The function presents a calculation of approximation gamma-
C    nucleon partial (gamma + N ==> m(pi) + N, m=2-8) cross-sections 
c    from experiments: S.I.Alekhin et al., CERN-HERA 87-01,
c    received from Igor Pshenichnov, Modified by KKG for LAQGSM 
c
	dimension sig(7),al(7,4),tth(7), wr(7), nal(7)
	dimension pl(0:3)
	data
     &(al(1,j),j=1,4) / 0.3317900,-0.0821800, 0.0997600,0.0001300/,
     &(al(2,j),j=1,4) / 0.0537900, 0.1350400,-0.0317300,0.0000000/,
     &(al(3,j),j=1,4) /-0.0009800, 0.0983250,-0.0475448,0.0000000/,
     &(al(4,j),j=1,4) / 0.2056800,-0.0629900, 0.0000000,0.0000000/,
     &(al(5,j),j=1,4) / 0.0619000,-0.0192100, 0.0000000,0.0000000/,
     &(al(6,j),j=1,4) / 0.1113700,-0.0409400, 0.0000000,0.0000000/,
     &(al(7,j),j=1,4) / 0.0336780,-0.0130250, 0.0000000,0.0000000/               	 
	data tth/ 0.321, 0.506 , 0.727, 0.952, 1.215, 1.481, 1.788/
	data wr/ 0.7, 0.75, 0.2667, 0.4381, 0.125, 0.2755, 0.05614/
	data nal/4,3,3,2,2,2,2/
	data alfa/2.0d0/
	do  nr=1,7
	  if(eg.le.tth(nr))  then
	    sig(nr) = 0.0d0
        else
	    x = log(eg/tth(nr))
	    nl = nal(nr)-1
	    call  plaguer(alfa,x,pl,nl)
	    f = 0.0d0
	    do  m=0,nl
	      f = f + al(nr,m+1)*pl(m)
          enddo
          f=f*x/exp(x/2.0d0)
          sig(nr) = f*f/wr(nr)*1000.0d0   !  [mkb]
        endif 
      enddo
	return
	end 
c  **************************************************************    
      subroutine plaguer(al,x,pl,n)
      implicit real*8 (a-h, o-z), integer (i-n)  
c  
C The recurrent calculation of LAGERR's polinoms up to order=3 
C (See for details: Handbook of mathematical functions. Ed. by 
C M.Abramowitz and I.A.Stegun), received from Igor Pshenichnov
c Modified by KKG for LAQGSM 
c
	dimension pl(0:3)
	data one, two/1.0d0,2.0d0/
	pl(0) = one
	if(n.eq.0)  then
	  return
      else  
	  pl(1) = al + one - x
        if(n.eq.1)  return
	  do  m=1,n-1
	    rm = float(m) 
	    pl(m+1) = ((two*rm + al + one - x)*pl(m) -
     -		      (rm + al)*pl(m-1))/(rm+one)
        enddo
	endif   		 
      return
      end
c  **************************************************************    
      function csgntot(ip,eg,emn)
      implicit real*8 (a-h, o-z), integer (i-n)
c
c     calculate total gamma + N cross section at energy eg
c
      common /csgn3_9/ s2_9(8)
	dimension s3_9(7)
      s2_9(1) = csgn2(ip,eg,emn)
	csgntot = s2_9(1)
	call  gnappr(eg,s3_9)
      do ir=1,7
	  s2_9(ir+1) = s3_9(ir) 
	  csgntot = csgntot + s3_9(ir) 
      enddo
      return
c
      end
c  **************************************************************    
	function csgn2(ip,eg,emn)
      implicit real*8 (a-h, o-z), integer (i-n)
c
c     calculate two body gamma + N cross section at energy eg
c
      common /nchapn/ nchap(13),nchan(13)
c
      csgn2 = 0.0d0
      do ir=1,13
        if(ip.eq.1) then
          ich = nchap(ir)   ! gamma + p
        else
          ich = nchan(ir)   ! gamma + n
        endif
        csgn2 = csgn2 + csgnh(ich,eg,emn)
c      for channel 20 isotopicaly symmetric channal is added:
c    sig0 + K+ = sig+ + K0 (for gamma + p) or 
c    sig0 + K0 = sig- + K+ (for gamma + n)
	  if(ich.eq.20) csgn2 = csgn2 + csgnh(ich,eg,emn)
      enddo
      return
c
      end
c  **************************************************************    
      function csgnh(ich,eg,emn)
      implicit real*8 (a-h, o-z), integer (i-n)
c
c     calculate n body gamma + N cross section at energy eg
c
      COMMON /csgnt/ WR(251),SIGR(22,251)
      data zro, two /0.d0, 2.d0/
      w = sqrt((two*eg+emn)*emn)
	csgnh = 0.
	if(w.lt.WR(1))  then
	  csgnh = zro
        return
	elseif(w.gt.WR(251))  then 
	  w1=WR(250)
        w2=WR(251)
		s1=SIGR(ich,250)
		s2=SIGR(ich,251)
		csgnh = s1 + (s2-s1)/(w2-w1)*(w-w1)
	    if(csgnh.lt.zro) csgnh = zro  
	  return
	endif
	do  iw=1,251
	  if(abs(w-WR(iw)).lt.1.0d-6)  then
	    csgnh = SIGR(ich,iw)
		return
	  elseif(w.lt.WR(iw).or.(w.gt.WR(iw).and.iw.eq.251))  then
	    w1=WR(iw-1)
		w2=WR(iw)
		s1=SIGR(ich,iw-1)
		s2=SIGR(ich,iw)
		csgnh = s1 + (s2-s1)/(w2-w1)*(w-w1)
		return
	  endif
	enddo
	return
	end  		     
c   ****************************************************************
       subroutine gntoh(v,u,tin1,pin,iin,pn,ipn,mv,np)
      implicit real*8 (a-h, o-z), integer (i-n)
c
c      main subroutine to calculate gamma + N ==> hadrons
c
      COMMON/NCASCA/NCAS,NCPRI
      common /csgn3_9/ s2_9(8)
      dimension v(3),pin(9),iin(5),pn(9),ipn(5)
	np = 0
	rnd = rndm(-1.)
	emn = pn(9)
c      separate two body channels
	stot = csgntot(ipn(1),tin1,emn)
      sig2  = s2_9(1)
	beta2 = sig2/stot
	if(rnd.le.beta2)  then
	  if(NCAS.ge.NCPRI) write(*,*) ' to gnto2h: tin1,stot,sig2=',
     &   tin1,stot,sig2
	  call gnto2h(v,u,tin1,sig2,pin,emn,ipn,mv,np)
	  if(NCAS.ge.NCPRI) write(*,*) ' from gnto2h: np=',np
	  if(np.lt.2)  then
c	    write(*,*) ' from gnto2h(np<2): eg,stot,sig2,np=',
c    &     tin1,stot,sig2,np
	  endif 
	  return
	else
c         channels gamma + N ==> n*pi + N, n=2,3,4,5,6,7,8
 	  if(NCAS.ge.NCPRI) write(*,*) ' to gntonh: tin1,stot,sig2=',
     &   tin1,stot,sig2
       call gntonh(u,tin1,ipn,mv,np)
	  if(NCAS.ge.NCPRI)  write(*,*) ' from gntonh: np=',np
	  if(np.lt.3)  then
	    write(*,*) ' from gntonh(np<3): eg,stot,s2_9,np=',
     &        tin1,stot,s2_9,np
	  endif 
	  return
	endif  		    	   	 	        
      return
	end
c   ****************************************************************
       subroutine gnto2h(v,u,tin1,sig2,pin,emn,ipn,mv,np)
      implicit real*8 (a-h, o-z), integer (i-n)
c
c      main subroutine to calculate gamma + N ==> 2 hadrons
c
      COMMON/MEMORY/PME(9,5999),IME(5,5999)
      COMMON /IDPME/ IDPME(5999)
      common /nchapn/ nchap(13),nchan(13)
      common /gbrems/ tgmin, tgmax, teqv, sxabs, ibrems 
      dimension v(3),pin(9),ipn(5),sig(13),pist(3),pnst(3)
      data zro  /0.0d0/, eps /0.0001d0/
	nrep = 0
    9	ssum = zro
	rnd = rndm(-1.)
	do  ic=1,13
	  if(ipn(1).eq.1)  then
	    ich = nchap(ic)
	  else
	    ich = nchan(ic)
	  endif 		 
	  sig(ic) =  csgnh(ich,tin1,emn)
	  ssum = ssum + sig(ic)
	enddo 
	ss = zro
	do  ic=1,13
	  ss = ss + sig(ic)
	  if(rnd.le.ss/ssum)  then
	    nch = ic
	    go  to  10
	  endif	 	     
	enddo    
      write(*,*) ' gnto2h: not found 2-body channel:',
     &	 'tin1,u,emn,sig2,sig,',
     &             'ssum=',tin1,u,emn,sig2,sig,ssum	  
c	stop
c     go  to  9
      return 
   10 continue
      if(ipn(1).eq.1)  then
c      target is proton
        if(nch.eq.1)        then  ! gamma + p ==> pi+ + n
	    id1 = 120
	    id2 = 1220	
	  elseif(nch.eq.2)    then  ! gamma + p ==> pi0 + p
	    id1 = 110
		id2 = 1120
	  elseif(nch.eq.3)    then  ! gamma + p ==> pi- + D++
	    id1 =-120
		id2 = 1111
	  elseif(nch.eq.4)    then  ! gamma + p ==> pi0 + D+
	    id1 = 110
		id2 = 1121
	  elseif(nch.eq.5)    then  ! gamma + p ==> pi+ + D0
	    id1 = 120
		id2 = 1221
	  elseif(nch.eq.6)    then  ! gamma + p ==> rho0 + p
	    id1 = 111
		id2 = 1120
	  elseif(nch.eq.7)    then  ! gamma + p ==> rho+ + n
	    id1 = 121
		id2 = 1220
	  elseif(nch.eq.8)    then  ! gamma + p ==> eta + p
	    id1 = 220
		id2 = 1120
	  elseif(nch.eq.9)    then  ! gamma + p ==> omeg + p
	    id1 = 221
		id2 = 1120
	  elseif(nch.eq.10)   then  ! gamma + p ==> K+ + L
	    id1 = 130
		id2 = 2130
	  elseif(nch.eq.11)   then  
	    if(rndm(-1.).le.0.5d0) then ! gamma + p ==> K+ + S0
	      id1 = 130
		  id2 = 1230
          else                        ! gamma + p ==> K0 + S+
	      id1 = 230
		  id2 = 1130
          endif
	  elseif(nch.eq.12)   then  ! gamma + p ==> etap + p
	    id1 = 330
		id2 = 1120
	  elseif(nch.eq.13)   then  ! gamma + p ==> phi + p
	    id1 = 331
		id2 = 1120
	  endif
	  ich = nchap(nch)
	else  	 	 
c      target is neutron
        if(nch.eq.1)        then  ! gamma + n ==> pi- + p
	    id1 =-120
	    id2 = 1120
	  elseif(nch.eq.2)    then  ! gamma + n ==> pi0 + n
	    id1 = 110
		id2 = 1220
	  elseif(nch.eq.3)    then  ! gamma + n ==> pi- + D+
	    id1 =-120
		id2 = 1121
	  elseif(nch.eq.4)    then  ! gamma + n ==> pi0 + D0
	    id1 = 110
		id2 = 1221
	  elseif(nch.eq.5)    then  ! gamma + n ==> pi+ + D-
	    id1 = 120
		id2 = 2221
	  elseif(nch.eq.6)    then  ! gamma + n ==> rho- + p
	    id1 =-121
		id2 = 1120
	  elseif(nch.eq.7)    then  ! gamma + n ==> rho0 + n
	    id1 = 111
		id2 = 1220
	  elseif(nch.eq.8)    then  ! gamma + n ==> eta + n
	    id1 = 220
		id2 = 1220
	  elseif(nch.eq.9)    then  ! gamma + n ==> omeg + n
	    id1 = 221
		id2 = 1220
	  elseif(nch.eq.10)   then  ! gamma + n ==> K0 + L
	    id1 = 230
		id2 = 2130
	  elseif(nch.eq.11)   then
	    if(rndm(-1.).le.0.5d0) then  ! gamma + n ==> K0 + S0
	      id1 = 230
		  id2 = 1230            
          else                         ! gamma + n ==> K+ + S-  
	      id1 = 130
		  id2 = 2230
		endif              
	  elseif(nch.eq.12)   then  ! gamma + n ==> etap + n
	    id1 = 330
		id2 = 1220
	  elseif(nch.eq.13)   then  ! gamma + n ==> phi + n
	    id1 = 331
		id2 = 1220
	  endif
	  ich = nchan(nch)
	endif
c  compute mass, strangeness, electic and baryon charge of hadrons
	am1 = AMASSF(id1)
	am2 = AMASSF(id2)
	is1 = IS(id1)
	is2 = IS(id2)
	ic1 = CHARGE(id1)
	ic2 = CHARGE(id2)
	ib1 = 0
	ib2 = 1 
c   check the threshold
      if(u.le.(am1+am2))  then
	if(id2.eq.1111.or.id2.eq.1121.or.id2.eq.1221.or.id2.eq.2221)
     &                     	   then  ! change mass of Delta
	    am2 = 1.080d0 + rndm(-1.)*(u-am1-1.080d0)
	    if(u.le.(am1+am2))  go  to  9
          go  to 11
        elseif(id1.eq.121.or.id1.eq.-121.or.id1.eq.111)  
     &                           then  ! change mass of Rho
	    am1 = 0.281d0 + rndm(-1.)*(u-am2-0.281d0)
	    if(u.le.(am1+am2))  go  to  9
	    go  to  11
        elseif(id1.eq.221)       then  ! change mass of omega
	    am1 = 0.660d0 + rndm(-1.)*(u-am2-0.660d0)
	    if(u.le.(am1+am2))  go  to  9
	    go  to  11
	else	 
c        write(*,*) ' gnto2h:  u < am1 + am2, u, am1,am2,id1,id2,',
c    &  'nch,sig(nch),ssum,sig2=',u, am1,am2,id1,id2,nch,sig(nch),sig2
            nrep = nrep + 1
	    if(nrep.lt.100)  then
	      go  to  9
            else
              np = 0
              return
	    endif
        endif
      endif  	    
c  compute the resonance life time
   11 ir1 = 0
	ir2 = 0
      if(id1.eq.121)        then      !  rho+
	  tau = TAUN(10)
	  ir1  = 1
	elseif(id1.eq.-121)   then      !  rho- 
	  tau = TAUN(11)
	  ir1 = 1            
	elseif(id1.eq. 111)   then      !  rho0 
	  tau = TAUN(16)
	  ir1 = 1            
	elseif(id1.eq. 221)   then      !  omeg 
	  tau = TAUN(17)
	  ir1 = 1            
	elseif(id1.eq. 331)   then      !  phi 
	  tau = TAUN(18)
	  ir1 = 1            
	elseif(id1.eq. 220)   then      !  eta 
	  tau = TAUN(8)
	  ir1 = 1            
	elseif(id1.eq. 330)   then      !  etap 
	  tau = TAUN(9)
	  ir1 = 1 
	endif
	if(id2.eq.1111.or.id2.eq.1121.or.id2.eq.1221.or.id2.eq.2221) then
	  call WMD(am2,tau,fmd)
	  ir2 = 1
	endif               
	r1 = rndm(-1.)
c  select the scattering angle in c.m.s.
c       for bremstrahlung gamma call inigamn for given tin1 
      if(ibrems.eq.1) call inigamn (tin1)
c
	cts=  cosgamm(ich,r1)
	phis= zro
c  compute the momenta of produced hadrons 
	call ABELQ(pin,v,u,pist,pnst,cts,phis,am1,am2) 
	e1 = sqrt(pist(1)**2+pist(2)**2+pist(3)**2+am1**2)
	e2 = sqrt(pnst(1)**2+pnst(2)**2+pnst(3)**2+am2**2)
	pme(1,mv+3) = zro   
	pme(2,mv+3) = zro   
	pme(3,mv+3) = zro   
	pme(4,mv+3) = pist(1)   
	pme(5,mv+3) = pist(2)   
	pme(6,mv+3) = pist(3)   
	pme(7,mv+3) = zro   
	pme(9,mv+3) = am1  
	ime(1,mv+3) = ic1 
	ime(2,mv+3) = 0 
	ime(3,mv+3) = is1
	ime(4,mv+3) = ib1 
	if(ir1.ne.0)  then
	  ime(5,mv+3) = nint(1000. *tau) 
      else
	  ime(5,mv+3) = 0
	endif   
	pme(8,mv+3) = e1 - am1
	pme(1,mv+1) = zro   
	pme(2,mv+1) = zro   
	pme(3,mv+1) = zro   
	pme(4,mv+1) = pnst(1)   
	pme(5,mv+1) = pnst(2)   
	pme(6,mv+1) = pnst(3)   
	pme(7,mv+1) = zro   
	pme(9,mv+1) = am2  
	ime(1,mv+1) = ic2 
	ime(2,mv+1) = 0 
	ime(3,mv+1) = is2
	ime(4,mv+1) = ib2 
	if(ir2.ne.0)  then
	  ime(5,mv+1) = nint(1000. *tau) 
        else
	  ime(5,mv+1) = 0
	endif 
	pme(8,mv+1) = e2 - am2
	IDPME(mv+3) = id1
	IDPME(mv+1) = id2
	np = 2
	psx = pme(4,mv+1) + pme(4,mv+3) 
	psy = pme(5,mv+1) + pme(5,mv+3) 
	psz = pme(6,mv+1) + pme(6,mv+3)
	es  = e1 + e2
	ies = ime(1,mv+1) + ime(1,mv+3) 
	iss = ime(3,mv+1) + ime(3,mv+3) 
	ibs = ime(4,mv+1) + ime(4,mv+3) 
c      check the conservation law
      if(abs(psx).gt.eps.or.abs(psy).gt.eps.or.abs(psz).gt.eps.
     &or.abs(es-u).gt.eps.or.ies.ne.ipn(1).or.iss.ne.0.or.ibs.ne.1)
     &write(*,*) ' gnto2h: psx,psy,psz,es,u,am1,am2,ies,iss,ibs=',
     &                     psx,psy,psz,es,u,am1,am2,ies,iss,ibs
	return
	end
c
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      subroutine gntonh(u,tin1,ipatne,mv,np)
      implicit real*8 (a-h, o-z), integer (i-n)
c
c     determining of secondary particles characteristics for
c     gamma + N ==> npi + N, n=2-8, interaction using GENBODL
c    last modification: 6 Dec. 2004 by KKG
c
	common /su_2wc/ su_2w3(3),su_2w4(4),
     &su_2w5(5),su_2w6(6),su_2w7(7),su_2w8(8),su_2w9(9),
     &isu_2p3(3,3),isu_2p4(4,4), 
     &isu_2p5(5,5),isu_2p6(6,6),isu_2p7(7,7),isu_2p8(8,8),isu_2p9(9,9),
     &isu_2n3(3,3),isu_2n4(4,4), 
     &isu_2n5(5,5),isu_2n6(6,6),isu_2n7(7,7),isu_2n8(8,8),isu_2n9(9,9)
      common /memory/ pmemo(9,5999),imemo(5,5999)
      common /csgn3_9/ s2_9(8)
      COMMON /IDPME/ IDPME(5999)
      dimension ipatne(5),
     & ama3(3),ama4(4),ama5(5),ama6(6),ama7(7),ama8(8),ama9(9),
     & ps3(5,3),ps4(5,4),ps5(5,5),ps6(5,6),ps7(5,7),ps8(5,8),ps9(5,9),
c
c SGM, 10/27/2011 correction sugested by KKG 
c     & sig(7),sigk(7),
     & sig(7),sigk(7)
c     & w3(3),w4(4),w5(5),w6(6),w7(7),w8(8),w9(9)
c
      data emneut, emprot, empich, empi0 /0.9395656d0,0.9382723d0,
     &                                    0.139568d0, 0.134973d0/
	data zro /0.0d0/, eps/0.0001d0/  
c
c  determine  the resulting pion multiplicity
	s28 = zro
      do  k=1,7
	  sig(k) = s2_9(k+1)
	  s28 = s28 + sig(k)
	  sigk(k) = s28
      enddo
	if(s28.le.zro)  then
	  np = 0
	  return
	endif
	rnd = rndm(-1.)
	do  k=1,7
	  npi = k + 1  
	  if(rnd. le. sigk(k)/s28) go  to  10
	enddo
	write(*,*) ' gntonh: npi, rnd, sigk=', npi, rnd, sigk
   10 continue
      np = npi + 1        
c  determine  the resulting particle types 
      if(np. eq. 3)  then      !  gamma + N ==> 2pi + N
	  s = zro
	  rnd = rndm(-1.)
	  do  n=1,np
	    nr = n
          s = s + su_2w3(n)
		if(rnd. le. s) go  to  30
	  enddo
   30   continue
        do  k =1,np
	    if(ipatne(1).eq.1)  then    
		  imemo(1,mv+k) = isu_2p3(nr,k)  ! target is proton 
		else   
		  imemo(1,mv+k) = isu_2n3(nr,k)  ! target is neutron
          endif
	    if(k.eq.1) then
		  imemo(4,mv+k) = 1            ! nucleon
		  if(imemo(1,mv+k).eq.1)  then
		    ama3(k) = emprot           ! proton
	  	  else
			ama3(k) = emneut           ! neutron
		  endif
		else
		  imemo(4,mv+k) = 0            ! pion
	      if(imemo(1,mv+k).eq.0)  then
	        ama3(k) = empi0            ! pi0
		  else
			ama3(k) = empich           ! pi+ or pi-
	      endif
		endif
        enddo
	  call  genbodl(np,u,ama3,ps3,w3)
	  do k=1,np
	    pmemo(4,mv+k) = ps3(1,k)
	    pmemo(5,mv+k) = ps3(2,k)
	    pmemo(6,mv+k) = ps3(3,k)
	    pmemo(7,mv+k) = zro
	    pmemo(9,mv+k) = ama3(k)
        enddo
      elseif(np. eq. 4)  then      !  gamma + N ==> 3pi + N
	  s = zro
	  rnd = rndm(-1.)
	  do  n=1,np
	    nr = n
          s = s + su_2w4(n)
		if(rnd. le. s) go  to  40
	  enddo
   40   continue
        do  k =1,np
	    if(ipatne(1).eq.1)  then    
		  imemo(1,mv+k) = isu_2p4(nr,k)  ! target is proton 
		else   
		  imemo(1,mv+k) = isu_2n4(nr,k)  ! target is neutron
          endif
	    if(k.eq.1) then
		  imemo(4,mv+k) = 1            ! nucleon
		  if(imemo(1,mv+k).eq.1)  then
		    ama4(k) = emprot           ! proton
	  	  else
			ama4(k) = emneut           ! neutron
		  endif
		else
		  imemo(4,mv+k) = 0            ! pion
	      if(imemo(1,mv+k).eq.0)  then
	        ama4(k) = empi0            ! pi0
		  else
			ama4(k) = empich           ! pi+ or pi-
	      endif
		endif
        enddo
	  call  genbodl(np,u,ama4,ps4,w4)
	  do k=1,np
	    pmemo(4,mv+k) = ps4(1,k)
	    pmemo(5,mv+k) = ps4(2,k)
	    pmemo(6,mv+k) = ps4(3,k)
	    pmemo(7,mv+k) = zro
	    pmemo(9,mv+k) = ama4(k)
        enddo
      elseif(np. eq. 5)  then      !  gamma + N ==> 4pi + N
	  s = zro
	  rnd = rndm(-1.)
	  do  n=1,np
	    nr = n
          s = s + su_2w5(n)
		if(rnd. le. s) go  to  50
	  enddo
   50   continue
        do  k =1,np
	    if(ipatne(1).eq.1)  then    
		  imemo(1,mv+k) = isu_2p5(nr,k)  ! target is proton 
		else   
		  imemo(1,mv+k) = isu_2n5(nr,k)  ! target is neutron
          endif
	    if(k.eq.1) then
		  imemo(4,mv+k) = 1            ! nucleon
		  if(imemo(1,mv+k).eq.1)  then
		    ama5(k) = emprot           ! proton
	  	  else
			ama5(k) = emneut           ! neutron
		  endif
		else
		  imemo(4,mv+k) = 0            ! pion
	      if(imemo(1,mv+k).eq.0)  then
	        ama5(k) = empi0            ! pi0
		  else
			ama5(k) = empich           ! pi+ or pi-
	      endif
		endif
        enddo
	  call  genbodl(np,u,ama5,ps5,w5)
	  do k=1,np
	    pmemo(4,mv+k) = ps5(1,k)
	    pmemo(5,mv+k) = ps5(2,k)
	    pmemo(6,mv+k) = ps5(3,k)
	    pmemo(7,mv+k) = zro
	    pmemo(9,mv+k) = ama5(k)
        enddo
      elseif(np. eq. 6)  then      !  gamma + N ==> 5pi + N
	  s = zro
	  rnd = rndm(-1.)
	  do  n=1,np
	    nr = n
          s = s + su_2w6(n)
		if(rnd. le. s) go  to  60
	  enddo
   60   continue
        do  k =1,np
	    if(ipatne(1).eq.1)  then    
		  imemo(1,mv+k) = isu_2p6(nr,k)  ! target is proton 
		else   
		  imemo(1,mv+k) = isu_2n6(nr,k)  ! target is neutron
          endif
	    if(k.eq.1) then
		  imemo(4,mv+k) = 1            ! nucleon
		  if(imemo(1,mv+k).eq.1)  then
		    ama6(k) = emprot           ! proton
	  	  else
			ama6(k) = emneut           ! neutron
		  endif
		else
		  imemo(4,mv+k) = 0            ! pion
	      if(imemo(1,mv+k).eq.0)  then
	        ama6(k) = empi0            ! pi0
		  else
			ama6(k) = empich           ! pi+ or pi-
	      endif
		endif
        enddo
	  call  genbodl(np,u,ama6,ps6,w6)
	  do k=1,np
	    pmemo(4,mv+k) = ps6(1,k)
	    pmemo(5,mv+k) = ps6(2,k)
	    pmemo(6,mv+k) = ps6(3,k)
	    pmemo(7,mv+k) = zro
	    pmemo(9,mv+k) = ama6(k)
        enddo
      elseif(np. eq. 7)  then      !  gamma + N ==> 6pi + N
	  s = zro
	  rnd = rndm(-1.)
	  do  n=1,np
	    nr = n
          s = s + su_2w7(n)
		if(rnd. le. s) go  to  70
	  enddo
   70   continue
        do  k =1,np
	    if(ipatne(1).eq.1)  then    
		  imemo(1,mv+k) = isu_2p7(nr,k)  ! target is proton 
		else   
		  imemo(1,mv+k) = isu_2n7(nr,k)  ! target is neutron
          endif
	    if(k.eq.1) then
		  imemo(4,mv+k) = 1            ! nucleon
		  if(imemo(1,mv+k).eq.1)  then
		    ama7(k) = emprot           ! proton
	  	  else
			ama7(k) = emneut           ! neutron
		  endif
		else
		  imemo(4,mv+k) = 0            ! pion
	      if(imemo(1,mv+k).eq.0)  then
	        ama7(k) = empi0            ! pi0
		  else
			ama7(k) = empich           ! pi+ or pi-
	      endif
		endif
        enddo
	  call  genbodl(np,u,ama7,ps7,w7)
	  do k=1,np
	    pmemo(4,mv+k) = ps7(1,k)
	    pmemo(5,mv+k) = ps7(2,k)
	    pmemo(6,mv+k) = ps7(3,k)
	    pmemo(7,mv+k) = zro
	    pmemo(9,mv+k) = ama7(k)
        enddo
      elseif(np. eq. 8)  then      !  gamma + N ==> 7pi + N
	  s = zro
	  rnd = rndm(-1.)
	  do  n=1,np
	    nr = n
          s = s + su_2w8(n)
		if(rnd. le. s) go  to  80
	  enddo
   80   continue
        do  k =1,np
	    if(ipatne(1).eq.1)  then    
		  imemo(1,mv+k) = isu_2p8(nr,k)  ! target is proton 
		else   
		  imemo(1,mv+k) = isu_2n8(nr,k)  ! target is neutron
          endif
	    if(k.eq.1) then
		  imemo(4,mv+k) = 1            ! nucleon
		  if(imemo(1,mv+k).eq.1)  then
		    ama8(k) = emprot           ! proton
	  	  else
			ama8(k) = emneut           ! neutron
		  endif
		else
		  imemo(4,mv+k) = 0            ! pion
	      if(imemo(1,mv+k).eq.0)  then
	        ama8(k) = empi0            ! pi0
		  else
			ama8(k) = empich           ! pi+ or pi-
	      endif
		endif
        enddo
	  call  genbodl(np,u,ama8,ps8,w8)
	  do k=1,np
	    pmemo(4,mv+k) = ps8(1,k)
	    pmemo(5,mv+k) = ps8(2,k)
	    pmemo(6,mv+k) = ps8(3,k)
	    pmemo(7,mv+k) = zro
	    pmemo(9,mv+k) = ama8(k)
        enddo
      elseif(np. eq. 9)  then      !  gamma + N ==> 8pi + N
	  s = zro
	  rnd = rndm(-1.)
	  do  n=1,np
	    nr = n
          s = s + su_2w9(n)
		if(rnd. le. s) go  to  90
	  enddo
   90   continue
        do  k =1,np
	    if(ipatne(1).eq.1)  then    
		  imemo(1,mv+k) = isu_2p9(nr,k)  ! target is proton 
		else   
		  imemo(1,mv+k) = isu_2n9(nr,k)  ! target is neutron
          endif
	    if(k.eq.1) then
		  imemo(4,mv+k) = 1            ! nucleon
		  if(imemo(1,mv+k).eq.1)  then
		    ama9(k) = emprot           ! proton
	  	  else
			ama9(k) = emneut           ! neutron
		  endif
		else
		  imemo(4,mv+k) = 0            ! pion
	      if(imemo(1,mv+k).eq.0)  then
	        ama9(k) = empi0            ! pi0
		  else
			ama9(k) = empich           ! pi+ or pi-
	      endif
		endif
        enddo
	  call  genbodl(np,u,ama9,ps9,w9)
	  do k=1,np
	    pmemo(4,mv+k) = ps9(1,k)
	    pmemo(5,mv+k) = ps9(2,k)
	    pmemo(6,mv+k) = ps9(3,k)
	    pmemo(7,mv+k) = zro
	    pmemo(9,mv+k) = ama9(k)
        enddo
      else
	  write(*,*) ' gntonh: np>9 ?, np=',np
	  stop
	endif 
	do  k=1,np
	  imemo(2,mv+k) = 0   
	  imemo(3,mv+k) = 0   
	  imemo(5,mv+k) = 0
	  pmemo(8,mv+k) = sqrt(pmemo(4,mv+k)**2+pmemo(5,mv+k)**2+
     +  pmemo(6,mv+k)**2+pmemo(9,mv+k)**2)-pmemo(9,mv+k)
	enddo 
	if(imemo(1,mv+1).eq.1)  then
	  IDPME(mv+1) = 1120        ! p
	else
	  IDPME(mv+1) = 1220        ! n 
	endif
	do  n=2,np
	  if(imemo(1,mv+n).eq.0)      then
	    IDPME(mv+n) =  110        ! pi0
	  elseif(imemo(1,mv+n).eq.1)  then
	    IDPME(mv+n) =  120        ! pi+ 
	  else     
	    IDPME(mv+n) = -120        ! pi-
	  endif   
      enddo
	psx = zro
	psy = zro
	psz = zro
	es  = zro
	ies = 0
	iss = 0
	ibs = 0
	do  k=1,np
        psx = psx + pmemo(4,mv+k)
        psy = psy + pmemo(5,mv+k)
        psz = psz + pmemo(6,mv+k)
        es  = es  + pmemo(8,mv+k) + pmemo(9,mv+k)
	  ies = ies + imemo(1,mv+k)
	  iss = iss + imemo(3,mv+k)
	  ibs = ibs + imemo(4,mv+k)
	enddo 
c      check the conservation law
      if(abs(psx).gt.eps.or.abs(psy).gt.eps.or.abs(psz).gt.eps.
     &or.abs(es-u).gt.eps.or.ies.ne.ipatne(1).or.iss.ne.0.or.ibs.ne.1)
     &write(*,*) ' gntonh: np,psx,psy,psz,es,u,ies,iss,ibs=',
     &                     np,psx,psy,psz,es,u,ies,iss,ibs
      return
      end
c
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      SUBROUTINE GENBODL(NP,TECM,AMASS,PCMS,WT)
      implicit real*8 (a-h, o-z), integer (i-n)
* Revision 1.1.1.1  1996/03/22 16:42:49  mclareni
c Received from Igor Pshenichnov, Oct.,2003
c Mofified by K.K. Gudima for LAQGSM, Nov., 2004
C   SUBROUTINE TO GENERATE N-BODY EVENT
C   ACCORDING TO FERMI LORENTZ-INVARIANT PHASE SPACE
C   ADAPTED FROM FOWL (CERN W505) SEPT. 1974 BY F. JAMES
C   EVENTS ARE GENERATED IN THEIR OWN CENTER-OF-MASS,
C   BUT MAY BE TRANSFORMED TO ANY FRAME USING LOREN4
C
C   INPUT TO SUBROUTINE IS THRU COMMON BLOCK GENIN
C             NP=NUMBER OF OUTGOING PARTICLES (.LT. 19)
C             TECM=TOTAL ENERGY IN CENTER-OF-MASS
C             AMASS(I)=MASS OF ITH OUTGOING PARTICLE
C             KGENEV=1 FOR CONSTANT CROSS SECTION
C                      2 FOR FERMI ENERGY-DEPENDANCE
C
C   OUTPUT FROM SUBROUTINE IS 
C             PCMS(1,I)=X-MOMENTUM IF ITH PARTICLE
C             PCMS(2,I)=Y-MOMENTUM IF ITH PARTICLE
C             PCMS(3,I)=Z-MOMENTUM IF ITH PARTICLE
C             PCMS(4,I)=ENERGY OF ITH PARTICLE
C             PCMS(5,I)=MOMENTUM OF ITH PARTICLE
C             WT=WEIGHT OF EVENT
      dimension AMASS(NP),PCMS(5,NP)
      DIMENSION EMM(18)
      DIMENSION RNO(50)
C--       PCM1 IS LINEAR EQUIV. OF PCM TO AVOID DOUBLE INDICES
      DIMENSION EM(18),PD(18),EMS(18),SM(18),FFQ(18),PCM1(90),PCM(5,18)
      EQUIVALENCE (PCM1(1),PCM(1,1))
C FFQ(N)=PI * (TWOPI)**(N-2) / (N-2)FACTORIAL
      DATA FFQ/0.d0,3.141592d0, 19.73921d0, 62.01255d0, 129.8788d0,
     1   204.0131d0,256.3704d0, 268.4705d0, 240.9780d0, 189.2637d0,
     3              132.1308d0,  83.0202d0,  47.4210d0,  24.8295d0,
     4               12.0006d0,   5.3858d0,   2.2560d0,   0.8859d0/
      DATA KNT,TWOPI,zro,one,two/0,6.2831853073d0,0.0d0,1.0d0,2.0d0/
c#if defined(CERNLIB_CDC)
      DATA KGENEV/ 2 /
c#endif
C        INITIALIZATION
      NT = NP
	do  i=1,NP
	  EM(i) = AMASS(i)
      enddo
      KNT=KNT + 1
      IF(KNT.GT.1) GOTO 100
c     WRITE(6,1160)
c     WRITE(6,1200) NP,TECM,(AMASS(JK),JK=1,NP)
  100 CONTINUE
      IF(NT.LT.2) GOTO 1001
      IF(NT.GT.18) GOTO 1002
      NTM1=NT-1
      NTM2=NT-2
      NTP1=NT+1
      NTNM4=3*NT - 4
      EMM(1)=EM(1)
      TM=zro
      DO 200 I=1,NT
        EMS(I)=EM(I)**2
        TM=TM+EM(I)
 200  SM(I)=TM
C        CONSTANTS DEPENDING ON TECM
      TECMTM=TECM-TM
      IF(TECMTM.LE.zro) GOTO 1000
      EMM(NT)=TECM
      IF(KGENEV.GT.1) GOTO 400
C        CONSTANT CROSS SECTION AS FUNCTION OF TECM
      EMMAX=TECMTM+EM(1)
      EMMIN=zro
      WTMAX=one
      DO 350 I=2,NT
      EMMIN=EMMIN+EM(I-1)
      EMMAX=EMMAX+EM(I)
  350 WTMAX=WTMAX*PDK(EMMAX,EMMIN,EM(I))
      WTMAXQ=one/WTMAX
      GOTO 455
C--      FERMI ENERGY DEPENDENCE FOR CROSS SECTION
  400 WTMAXQ=TECMTM**NTM2*FFQ(NT) / TECM
C        CALCULATION OF WT BASED ON EFFECTIVE MASSES EMM
  455 CONTINUE
C--               FILL RNO WITH 3*NT-4 RANDOM NUMBERS,
C--            OF WHICH THE FIRST NT-2 ARE ORDERED.
      DO 457 I= 1, NTNM4
  457 RNO(I)=rndm(-1.)
      IF(NTM2) 900,509,460
  460 CONTINUE
      CALL FLPSOR(RNO,NTM2)
      DO 508 J=2,NTM1
  508 EMM(J)=RNO(J-1)*(TECMTM)+SM(J)
  509 WT=WTMAXQ
      IR=NTM2
      DO 530 I=1,NTM1
      PD(I)=PDK(EMM(I+1),EMM(I),EM(I+1))
  530 WT=WT*PD(I)
C--       COMPLETE SPECIFICATION OF EVENT (RAUBOLD-LYNCH METHOD)
      PCM(1,1)=zro
      PCM(2,1)=PD(1)
      PCM(3,1)=zro
      DO 570 I=2,NT
        PCM(1,I)=zro
        PCM(2,I)=-PD(I-1)
        PCM(3,I)=zro
        IR=IR+1
        BANG=TWOPI*RNO(IR)
        CB=COS(BANG)
        SB=SIN(BANG)
        IR=IR+1
        C=two*RNO(IR)-one
        S=SQRT(one-C*C)
        IF(I.EQ.NT) GOTO 1567
        ESYS=SQRT(PD(I)**2+EMM(I)**2)
        BETA=PD(I)/ESYS
        GAMA=ESYS/EMM(I)
        DO 568 J=1,I
          NDX=5*J - 5
          AA= PCM1(NDX+1)**2 + PCM1(NDX+2)**2 + PCM1(NDX+3)**2
          PCM1(NDX+5)=SQRT(AA)
          PCM1(NDX+4)=SQRT(AA+EMS(J))
          CALL ROTES2(C,S,CB,SB,PCM,J)
          PSAVE=GAMA*(PCM(2,J)+BETA*PCM(4,J))
  568   PCM(2,J)=PSAVE
        GOTO 570
 1567   DO 1568 J=1,I
          AA=PCM(1,J)**2 + PCM(2,J)**2 + PCM(3,J)**2
          PCM(5,J)=SQRT(AA)
          PCM(4,J)=SQRT(AA+EMS(J))
          CALL ROTES2(C,S,CB,SB,PCM,J)
 1568   CONTINUE
  570 CONTINUE
      do  i=1,NT
	  do  j=1,5
	    PCMS(j,i) = PCM(j,i)
        enddo
      enddo
  900 CONTINUE
      RETURN
C          ERROR RETURNS 
 1000 WRITE(*,1100)
      GOTO 1050
 1001 WRITE(*,1101)
      GOTO 1050
 1002 WRITE(*,1102)
 1050 WRITE(*,1150) KNT
      WRITE(*,1200) NP,TECM,(AMASS(JK),JK=1,NP)
      STOP
 1100 FORMAT(' GENBOD:AVAILABLE ENERGY NEGATIVE' )
 1101 FORMAT(' GENBOD:LESS THAN 2 OUTGOING PARTICLES' )
 1102 FORMAT(' GENBOD:MORE THAN 18 OUTGOING PARTICLES' )
 1150 FORMAT(' GENBOD:ABOVE ERROR DETECTED IN GENBOD AT CALL NUMBER',I7)
 1160 FORMAT(' GENBOD: FIRST CALL TO SUBROUTINE GENBOD' )
 1200 FORMAT(' GENBOD: INPUT DATA NP=   ',I6/'   TECM=',E16.7,
     &       '  PARTICLE MASSES=',5E15.5/(42X,5E15.5))
      END
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      FUNCTION PDK(A,B,C)
      implicit real*8 (a-h, o-z), integer (i-n)
c Received from Igor Pshenichnov, Oct.,2003
c Mofified by K.K. Gudima for LAQGSM, Nov., 2004
C--  CALLED FROM -  GENBOD
C     PDK = SQRT(A*A+(B*B-C*C)**2/(A*A) - 2.0*(B*B+C*C))/2.0
      A2 = A*A
      B2 = B*B
      C2 = C*C
      SQR_ARGUMENT=A2 + (B2-C2)**2/A2 - 2.0d0*(B2+C2)
      IF(SQR_ARGUMENT.GT.0.0d0) THEN
        PDK = 0.5d0*SQRT(SQR_ARGUMENT)
      ELSE
        PDK = 0.0d0
      ENDIF
      RETURN
      END
c   ****************************************************************
      SUBROUTINE FLPSOR(RNO,NTM2)
      implicit real*8 (a-h, o-z), integer (i-n)
c Received from Igor Pshenichnov, Oct.,2003
c Mofified by K.K. Gudima for LAQGSM, Nov., 2004
C--       ORDER THE FIRST NTM2 RANDOM NUMBERS
C--         TWO IS A SPECIAL CASE (FASTER)
C--       CALLED FROM GENBOD
      dimension RNO(50) 
      IF(NTM2 - 2) 200,160,110
  110 KM1 = NTM2 - 1
      DO 150 I= 1, KM1
        IQUIT = 0
        NI = NTM2 - I
        DO 140 J= 1, NI
          IF(RNO(J) - RNO(J+1)) 140,140,120
  120     SAV = RNO(J)
          RNO(J) = RNO(J+1)
          RNO(J+1) = SAV
          IQUIT = 1
  140   CONTINUE
        IF(IQUIT) 200,200,150
  150 CONTINUE
      GOTO 200
  160 IF(RNO(1).LE.RNO(2)) GOTO 200
      SAV = RNO(1)
      RNO(1) = RNO(2)
      RNO(2) = SAV
  200 CONTINUE
      RETURN
      END
c   ****************************************************************
      SUBROUTINE ROTES2(C,S,C2,S2,PR,I)
      implicit real*8 (a-h, o-z), integer (i-n)
c Received from Igor Pshenichnov, Oct.,2003
c Mofified by K.K. Gudima for LAQGSM, Nov., 2004
C--  CALLED FROM - GENBOD
C         THIS SUBROUTINE NOW DOES TWO ROTATIONS (XY AND XZ)
      DIMENSION PR(50)
      K1 = 5*I - 4
      K2 = K1 + 1
      SA = PR(K1)
      SB = PR(K2)
      A      = SA*C - SB*S
      PR(K2) = SA*S + SB*C
      K2 = K2 + 1
      B = PR(K2)
      PR(K1) = A*C2 - B*S2
      PR(K2) = A*S2 + B*C2
      RETURN
      END
c   ****************************************************************
         subroutine inigamn (egamma)

c ======================================================================
c
c    Main routine to extract ds/do for channel 1-22:
c
c    Written by K. K. Gudima, Fall 2003?
c    Modified by AJS, May, 2004.
c    Modified by KKG, Nov., 2004
c ======================================================================

      implicit real*8 (a-h, o-z), integer (i-n)

c ======================================================================
 
      common /ixsgpn/ thetai(181),cthetai(181),si(22,182),ri(22,181)
      common /xsecgn0/ xsectd(22,50,0:18),ecm(22,50),elg(22,50)
      common /xsecgpq/theta(19), ctheta(19)

      dimension s(22,19), st(22), r(22,19)

      data dtheti /1.0d0/
      data zro, two /0.d0, 2.d0/, pi/3.141592d0/
c ======================================================================
      twpi =two*pi
      degrad = pi/180.0d0 

      do j = 1,181
        thetai(j) = dble(j-1)*dtheti
        cthetai(j) = cos(thetai(j)*degrad)
      enddo
      do jch = 1,22
        if (egamma.le.elg(jch,2))      then
          ieg1 = 2
          ieg2 = 3
        elseif (egamma.ge.elg(jch,50))  then
          ieg1 = 49
          ieg2 = 50
        else
          do ie = 3,50
            if (egamma.ge.elg(jch,ie-1) .and. egamma.le.elg(jch,ie))
     &        then
              ieg1 = ie-1
              ieg2 = ie
              go to 10
            endif
          enddo
        endif
   10   if (ieg1.lt.2 .or. ieg1.gt.49 .or. ieg2.lt.3 .or. ieg2.gt.50)
     &    then
          write (*,*) '  stop in INIGAMN: ieg1, ieg2 = ', ieg1, ieg2
          stop
        endif
        eg1 = elg(jch,ieg1)                          
        eg2 = elg(jch,ieg2)
        sint = zro
        do j = 1,19                        
		s1 = xsectd(jch,ieg1,j-1)     
          s2 = xsectd(jch,ieg2,j-1)     
          s(jch,j) = s1 + ((s2 - s1)/(eg2 - eg1))*(egamma - eg1)
          if (j.ge.2) then
            dom = twpi*(ctheta(j-1) - ctheta(j))
            sint = sint + dom*(s(jch,j-1) + s(jch,j))/two
          endif
          r(jch,j) = sint
          enddo
	    if(sint.lt.zro)  sint = zro
          st(jch) = sint
          do j = 1,19
            if (sint.gt.zro)  r(jch,j) = r(jch,j)/sint
          enddo
      end do
      do jch = 1,22
        sint = zro
        do j = 1,181
          si(jch,j) = qintxsq(thetai(j), theta, s, jch, 22, 19)
          if (j.ge.2) then
            dom = twpi*(cthetai(j-1) - cthetai(j))
            sint = sint + dom*(si(jch,j-1) + si(jch,j))/two
          endif
          ri(jch,j) = sint
        end do
	  if(sint.lt.zro)  sint = zro
        si(jch,182) = sint
        do j = 1,181
          if (sint.gt.zro)  ri(jch,j) = ri(jch,j)/sint
        enddo
      end do

      return
      end                                           

c ======================================================================

      function qintxsq (x, th, se, l, m, n)

c ======================================================================
c
c    Interpolation of gamma+N differential cross sections.
c
c    Written by K. K. Gudima, Fall 2003?
c    Modified by AJS, May, 2004.
c
c ======================================================================

      implicit real*8 (a-h, o-z), integer (i-n)

c ======================================================================

      dimension th(n), se(m,n)

      data zro, one /0.d0, 1.d0/

c ======================================================================

      do k = 1,n-1
        if (abs(x - th(k)).le.1.0d-5)  then
          qintxsq = se(l,k)
          qintxsq = max(qintxsq, zro)
          return
        endif
      enddo
      do k = 2,n-1
        if (x.lt.th(k))  then
          k1 = k - 1
          k2 = k
          k3 = k + 1
          go to 10
        endif
      enddo
      k1 = n - 2
      k2 = n - 1
      k3 = n
   10 continue
      x1 = th(k1)
      x2 = th(k2)
      x3 = th(k3)
      y1 = se(l,k1)
      y2 = se(l,k2)
      y3 = se(l,k3)   
      d  = (x2 - x3)*x1*x1 + (x3 - x1)*x2*x2 + (x1 - x2)*x3*x3
      da = (x2 - x3)*y1    + (x3 - x1)*y2    + (x1 - x2)*y3
      db = (y2 - y3)*x1*x1 + (y3 - y1)*x2*x2 + (y1 - y2)*x3*x3
      dc = (x2*y3 - x3*y2)*x1*x1 + (x3*y1 - x1*y3)*x2*x2 + 
     &     (x1*y2 - x2*y1)*x3*x3
      d1 = zro
      if (d.ne.zro) d1 = one/d
      a = da*d1
      b = db*d1
      c = dc*d1
      qintxsq = a*x*x + b*x + c
      qintxsq = max(qintxsq,zro)
      return

c ======================================================================
      end
c  ***************************************************************
      function cosgamm (j0, r1)

c ======================================================================
c
c     Cosine calculation for two body gamma +N ==> hadrons reactions 
c     using linear interpolation of function cos(theta)= f(r1), 
c     r1=random number, j0=number of channel (see channel1.tab)
c     Energy of gamma is fixed in INIGAM.
c
c   Called by: COSEL COSEX 
c
c   Written by K. K. Gudima, Nov. 2003
c   Edited by A. J. Sierk, LANL T-16, March, 2004.
c   Removed call to RNDM to calling subprogram, AJS, March, 2004.
c   Extended to channels 1-22 by KKG, Nov., 2004
c ======================================================================

      implicit real*8 (a-h, o-z), integer (i-n)

c ======================================================================

      common /ixsgpn/ thetai(181),cthetai(181),si(22,182),ri(22,181)

      data zro, one, pi /0.d0, 1.d0, 3.141592d0/

c ======================================================================
      jg = j0
      if (j0.lt.1 .or. j0.gt.22)  then
        write (*, *)  ' In COSGAMM, j0 =', j0
        stop
      endif 
	degrad = pi/180.0d0
        do ir = 1,181
        rrri = r1 - ri(jg,ir)
        if (abs(rrri).lt.1.0d-5)  then
          cth = cos(thetai(ir)*degrad)
          go to 20
        endif
        enddo   
        do ir = 2,181
        rrri = r1 - ri(jg,ir)
        if (rrri.lt.zro) then
          ir1 = ir - 1
          ir2 = ir
          go to 10
        endif
        enddo
      ir1 = 180
      ir2 = 181
   10 continue
      x1 = ri(jg,ir1)  
      x2 = ri(jg,ir2)  
      y1 = thetai(ir1)
      y2 = thetai(ir2)
      th = y1 + (y2 - y1)*((r1 - x1)/(x2 - x1))
      cth = cos(th*degrad)
   20 temp1 = abs(cth)
      if (temp1.le.one) then
        cosgamm = cth
      else
        cosgamm = sign(one, cth)
      endif
c   kkg  12/10/04
      if(jg.eq.17.or.jg.eq.18)  cosgamm = -cosgamm
c
      return

c ======================================================================
      end