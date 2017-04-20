      SUBROUTINE POINTN(M,P1,P2,IP1,IP2,V1,U1,TR2,SIG,SAB,TINT,NR,N1,
     *N2,NA1,NA2,MV,DELTA,R0X,R0Y,R0Z)
c ======================================================================
c
c  Determination of interaction hadron pair (N1,P1,IP1),(N2,P2,IP2)
c  and type NR of collision:
c  NR =1,      collision of nucleon N1  from  projetile and
c                           nucleon N2  from  target;
c  NR =2,      collision of cascade particle N1 and 
c                           nucleon N2  from  projetile            
c  NR =3,      collision of cascade particle N1 and 
c                           nucleon N2  from  target;
c  NR =5,      collision of cascade particle N1 and 
c                           cascade particle N2
c  NR =4,      decay of resonance N1. 
c  Interaction time TINT is replaced by TINT + TAU, where 
c  TAU is min(TAU_NR), NR=1,2,3,4,5.
c  Calls:  RAPID1,RAPID2,RAPID3,RAPIDD and RAPID4 - to calculate
c          time of next interaction TAU_NR
c
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      COMMON/INDINT/INDINT
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
       COMMON/MEMORY/PMEMO(9,5999),IMEMO(5,5999)
      COMMON/TAUE/TPTE,TYPE,TYTE
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      COMMON/DTINT/DTINT
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T00,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/INTCEN/IPO,INT1(300),INT2(5999),INT3(5999),
     *              INT4(100000),IJPA,IST,JST,NRST
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/PIDABS/PIDABS
      COMMON/INTCC/ INTCC
      COMMON/TAUIJ/ TAUK,TPTS,TYPS,TYTS,TYYS,TIJ4(100000)
      COMMON/TLIMIT/TLIMIT
      COMMON/CSLID/CLIDER(5999)
      COMMON/IACT/ IACT/CVALON/ IVALON
      COMMON /HOLPT/ PH1(3),PH2(3),RH1(3),RH2(3),EHOL1,EHOL2,TFP,TFT
      COMMON /IDPME/ IDPME(5999)
      COMMON /IDN12/ ID1,ID2
      COMMON /BARPOT/ POT
	COMMON /PARINC/ PINC(3),EINC
      DIMENSION P1(9),P2(9),IP1(5),IP2(5),V1(3)
     *,PT1(9),PT2(9),IPT1(5),IPT2(5),YP1(9),IYP1(5),YT1(9),IYT1(5)
      TFRO1=TINT
			
    8 CONTINUE
      R0=SQRT(R0X**2+R0Y**2+R0Z**2)*1.01
      IPO=0
    9 TPT=-0.1
      TYP=-0.1
      TYT=-0.1
      TAUD=-0.1
      TYY=-0.1
      TAU=-0.1
      TY=-.1
      TFRO2=TINT
c      write( *,1991) NCAS,TINT
1991  FORMAT('+',50X,I6,1X,E10.3)
      RT=SQRT((RADP(1)-RADT(1))**2+(RADP(2)-RADT(2))**2+
     +(RADP(3)-RADT(3))**2)
      IF(RT-R0)10,10,12
   10 IF(NA1.LT.1.OR.NA2.LT.1)   GO  TO  12
      IF(IPO.NE.1)               GO  TO  11
      IF(TPTS.LE.0.)             GO  TO  12
      IF((TPTS-TAUK).LE.0.0001)  GO  TO  11
      TPT=TPTS-TAUK
      GO  TO  12
   11 CONTINUE
      CALL  RAPID1(NA1,NA2,DELTA,PT1,IPT1,PT2,IPT2,NPT1,NPT2,TPT,DLPT)
   12 IF(MV) 124,124,100
  100 IF(NA1.LE.0)                             GO  TO  101
      IF(IPO.NE.1)                             GO  TO  13
      IF(TYPS.LE.0.)                           GO  TO  101
      IF((TYPS-TAUK).LE.0.0001)                GO  TO  13
      TYP=TYPS-TAUK
      IF(TYP.LT.PMEMO(7,NYP1).AND.CLIDER(NYP1).LT.0.3)
     *                                         GO  TO  13
      IF(TYP.LT.PMEMO(7,NYP1).AND.IVALON.EQ.0)
     *                                         GO  TO  13
      GO  TO  101
   13 CONTINUE
      CALL RAPID2(MV,NA1,DELTA,YP1,IYP1,NYP1,NYP2,TYP,DLYP)
  101 IF(NA2.LT.1)   GO  TO  14
  102 CONTINUE
      IF(IPO.NE.1)                             GO  TO  103
      IF(TYTS.LE.0.)                           GO  TO  14
      IF((TYTS-TAUK).LE.0.0001)                GO  TO  103
      TYT=TYTS-TAUK
      IF(TYT.LT.PMEMO(7,NYT1).AND.CLIDER(NYT1).LT.0.3)
     *                                         GO  TO  103
      IF(TYT.LT.PMEMO(7,NYT1).AND.IVALON.EQ.0)
     *                                         GO  TO  103
      GO  TO  14
  103 CONTINUE
      CALL RAPID3(MV,NA2,DELTA,YT1,IYT1,NYT1,NYT2,TYT,DLYT)
   14 CONTINUE
      CALL  RAPIDD(MV,TAUD,ND)
1992  CONTINUE
      IF(MV.LT.2.OR.INTCC.EQ.0.OR.TINT.GT.TLIMIT) GO  TO 124
      IF(IPO.NE.1)                               GO  TO  114
      IF(TYYS.LE.0.)                             GO  TO  124
      IF((TYY-TAUK).LE.0.0001)                   GO  TO  114
      TYY=TYYS-TAUK
      IF(TYY.LT.PMEMO(7,NYY1).AND.CLIDER(NYY1).LT.0.3)
     *                                            GO  TO 114
      IF(TYY.LT.PMEMO(7,NYY2).AND.CLIDER(NYY2).LT.0.3)
     *                                            GO  TO 114
      IF(TYY.LT.PMEMO(7,NYY1).AND.IVALON.EQ.0)
     *                                            GO  TO 114
      IF(TYY.LT.PMEMO(7,NYY2).AND.IVALON.EQ.0)
     *                                            GO  TO 114
      IF(NYY1.EQ.NYY2)                            GO  TO 114
      GO  TO  124
  114 CONTINUE
      CALL  RAPID4(MV,NYY1,NYY2,TYY)
  124 CONTINUE
      IF(TYP)19,19,15
   15 IF(TYT)16,16,17
   16 TY=TYP
      NR=2
      CLID1=CLIDER(NYP1)
      CLID2=1.
C
                                       GO  TO  20
   17 IF(TYP-TYT) 16,16,18
   18 TY=TYT
      NR=3
      CLID1=CLIDER(NYT1)
      CLID2=1.
                                        GO TO 20
   19 IF(TYT)20,20,18
   20 IF(TPT)24,24,21
   21 IF(TY)23,23,22
   22 IF(TPT-TY)  23,23,25
   23 TAU=TPT
      NR=1
      CLID1=1.
      CLID2=1.
C
                                    GO  TO   26
   24 IF(TY)26,26,25
   25 TAU=TY
C
                                    GO  TO   26
   26 IF(TAUD.LT.0.)                GO  TO  126
      IF(TAUD.GT.TAU.AND.TAU.GE.0.) GO  TO  126
      TAU=TAUD
      N1=ND
      NR=4
  126 IF(TYY.LT.0.)                 GO  TO  226
      IF(TYY.GT.TAU.AND.TAU.GE.0.)  GO  TO  226
      TAU=TYY
      N1=NYY1
      N2=NYY2
      NR=5
      CLID1=CLIDER(NYY1)
      CLID2=CLIDER(NYY2)
  226 IF(TAU.LT.0.)                 GO  TO  41
C
      IF(NCAS.GE.NCPRI)  THEN
c      write(16,6000) MV,IPO,NR,TPT,TYP,TYT,TYY,TAUD,TAU
 6000 FORMAT(1X,'MV,IPO,NR,TPT,TYP,TYT,TYY,TAUD,TAU=',I4,2I2,6(1X,F6.3))
      ENDIF
C
      TAUK=TAU
      IF(MV)29,29,27
   27 DO 28 K=1,MV
      EK=PMEMO(8,K)+PMEMO(9,K)
      PMEMO(1,K)=PMEMO(1,K)+PMEMO(4,K)/EK*TAUK
      PMEMO(2,K)=PMEMO(2,K)+PMEMO(5,K)/EK*TAUK
      PMEMO(3,K)=PMEMO(3,K)+PMEMO(6,K)/EK*TAUK
      PMEMO(7,K)=PMEMO(7,K)-TAUK
      IF(PMEMO(7,K).GT.0.)          GO  TO  128
      PMEMO(7,K)=0.
      CLIDER(K)=1.
  128 CONTINUE
      IF(IMEMO(5,K).EQ.0)           GO  TO  28
      IMEMO(5,K)=IMEMO(5,K)-INTG(TAUK*1000.)
      IF(IMEMO(5,K).LE.0)  IMEMO(5,K)=1
   28 CONTINUE
   29 TINT=TINT+TAU
      IF(IJPA.LE.0)                 GO  TO  291
      DO  290  IJ=1,IJPA
  290 TIJ4(IJ)=TIJ4(IJ)-TAUK
  291 CONTINUE
      DO  30  K=1,3
      RADT(K)=RADT(K)+TAU*VTA(K)
      RADP(K)=RADP(K)+TAU*VPR(K)
   30 CONTINUE
      GO  TO (33,35,37,130,230),NR
C  ====> RAPIDD
  130 continue
      ID1=IDPME(N1)
      M=1
      GO  TO  135
C  ====> RAPID4
  230 DO 231  K=1,9
      P1(K)=PMEMO(K,N1)
      P2(K)=PMEMO(K,N2)
      IF(K.LE.5)   IP1(K)=IMEMO(K,N1)
      IF(K.LE.5)   IP2(K)=IMEMO(K,N2)
  231 CONTINUE
      ID1=IDPME(N1)
      ID2=IDPME(N2)
      GO  TO  32
C  ====> RAPID1
   33 DO  34  K=1,9
      P1(K)=PT1(K)
      P2(K)=PT2(K)
      IF(K.LE.5)   IP1(K)=IPT1(K)
      IF(K.LE.5)   IP2(K)=IPT2(K)
   34 CONTINUE
      N1=NPT1
      N2=NPT2
      DLK=DLPT
      ID1=1120
      IF(IP1(1).EQ.0) ID1=1220
      ID2=1120
      IF(IP2(1).EQ.0) ID2=1220
      GO  TO   32
C  ====> RAPID2
   35 DO  36  K=1,9
      P1(K)=YP1(K)
      IF(K.LE.5)   IP1(K)=IYP1(K)
   36 CONTINUE
	PINC(1)=YP1(4)
	PINC(2)=YP1(5)
	PINC(3)=YP1(6)
	EINC=SQRT(PINC(1)**2+PINC(2)**2+PINC(3)**2+YP1(9)**2)
      N1=NYP1
      N2=NYP2
      DLK=DLYP
                     GO  TO  39
C  ====> RAPID3
   37 DO  38  K=1,9
      P1(K)=YT1(K)
      IF(K.LE.5)   IP1(K)=IYT1(K)
   38 CONTINUE
	PINC(1)=YT1(4)
	PINC(2)=YT1(5)
	PINC(3)=YT1(6)
	EINC=SQRT(PINC(1)**2+PINC(2)**2+PINC(3)**2+YT1(9)**2)
      N1=NYT1
      N2=NYT2
      DLK=DLYT
      CALL  PARTNQ(2,N2,P2,IP2)
      TFT=TFERMIQ(P2(1),P2(2),P2(3),2)
      EHOL2=P2(8)
      P2(9)=0.94
      POT=0.
      IF(IP1(3).EQ.0.AND.IP1(5).EQ.0.AND.IP1(4).EQ.1)  POT=TFT+EPS2
c  kkg 28.10.03
c      IF(IP1(3).EQ.0.AND.IP1(5).EQ.0.AND.IP1(4).EQ.0)  POT=VPI
      IF(IP1(3).EQ.0.AND.IP1(5).EQ.0.AND.IP1(4).EQ.0.
     &   and.ABS(P1(9)-0.140).lt.0.1)                  POT=VPI
C ANTI-NUCLEON
      IF(IP1(3).EQ.0.AND.IP1(5).EQ.0.AND.IP1(4).EQ.-1)
     *  POT=POTBAR(2,P2(1),P2(2),P2(3))
      P1(8)=P1(8)+POT
      C=SQRT(P1(8)/(P1(8)-POT)*(P1(8)+2.*P1(9))/(P1(8)-POT+2.*P1(9)))
      DO  339 K=1,3
      P1(K+3)=P1(K+3)*C
      PH2(K)=-P2(3+K)
  339 RH2(K)=P2(K)
      ID1=IDPME(N1)
      ID2=1120
      IF(IP2(1).EQ.0) ID2=1220
      GO  TO  32
   39 CONTINUE
      CALL  PARTNQ(1,N2,P2,IP2)
      TFP=TFERMIQ(P2(1),P2(2),P2(3),1)
      EHOL1=P2(8)
      P2(9)=0.94
      POT=0.
      IF(IP1(3).EQ.0.AND.IP1(5).EQ.0.AND.IP1(4).EQ.1)  POT=TFP+EPS1
      IF(IP1(3).EQ.0.AND.IP1(5).EQ.0.AND.IP1(4).EQ.0)  POT=VPI
C ANTI-NUCLEON
      IF(IP1(3).EQ.0.AND.IP1(5).EQ.0.AND.IP1(4).EQ.-1)
     *  POT=POTBAR(1,P2(1),P2(2),P2(3))
      P1(8)=P1(8)+POT
      C=SQRT(P1(8)/(P1(8)-POT)*(P1(8)+2.*P1(9))/(P1(8)-POT+2.*P1(9)))
      DO  239 K=1,3
      P1(K+3)=P1(K+3)*C
      PH1(K)=-P2(3+K)
  239 RH1(K)=P2(K)
      ID1=IDPME(N1)
      ID2=1120
      IF(IP2(1).EQ.0) ID2=1220
      GO  TO  32
   32 CONTINUE
C
      IF((ID1.EQ.120.OR.ID1.EQ.-120.OR.ID1.EQ.110.OR.ID1.EQ.1120.OR.
c  kkg  28.10.03     ! including gamma(ID1=10)
c     &    ID1.EQ.1220).AND.(ID2.EQ.1120.OR.ID2.EQ.1220)) THEN
     &    ID1.EQ.1220.or.ID1.eq.10).
     &    AND.(ID2.EQ.1120.OR.ID2.EQ.1220)) THEN
c  kkg  28.10.03
        CALL  SLQEKQ(L,MS,MQ,KSIN,ME,IP1,IP2)
        CALL  TINVUQ(P1,P2,U1,V1,TR2)
        if(ID1.eq.10)  then
          SITO= csgntot(IP2(1),TR2,P2(9))/1000.   ! tot. g+N cr.sec.
        else 
          SITO=CROSEG(L,MS,MQ,KSIN,0,TR2,P1(9),IP1(5))
        endif
      ELSE
        PX1=P1(4)
        PY1=P1(5)
        PZ1=P1(6)
        AM1=P1(9)
        PX2=P2(4)
        PY2=P2(5)
        PZ2=P2(6)
        AM2=P2(9)
        CALL CROSEC(1,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
        CALL  SLQEKQ(L,MS,MQ,KSIN,ME,IP1,IP2)
        CALL  TINVUQ(P1,P2,U1,V1,TR2)
      ENDIF
      SIG=SITO
      SAB=0.
      GO  TO  (135,131,132,135,135),NR
  131 IF(AN1.LT.2.1.OR.ZN1.LT.1.1)          GO  TO  135
      GO  TO  133
  132 IF(AN2.LT.2.1.OR.ZN2.LT.1.1)          GO  TO  135
  133 CONTINUE
c  kkg  28.10.03     ! including gamma(ID1=10)
      IF(ID1.EQ.120.OR.ID1.EQ.-120.OR.ID1.EQ.110.or.ID1.eq.10)
     &   GO  TO  134
                                            GO  TO  135
  134 CONTINUE
c  kkg  28.10.03     ! including gamma(ID1=10)
C   !!!  ONLY FOR PION(+-0) and gamma  !!!
      SAB=CROSEG(L,MS,MQ,KSIN,3,P1(8),P1(9),IP1(5))
      if(ID1.ne.10)  SAB=SAB*PIDABS
c  kkg  28.10.03     
  135 CONTINUE
      IPO=1
      IST=N1
      JST=N2
      NRST=NR
      TPTS=TPT
      TYPS=TYP
      TYTS=TYT
      TYYS=TYY
      IF(NR.EQ.4)        GO  TO  42
      SIGV=SIG
c  kkg  28.10.03     ! not include  gamma(ID1=10)
      IF(IVALON.NE.0.and.ID1.ne.10)    SIGV=CLID1*CLID2*SIG
c  kkg  04/02/04
      if(ID1.eq.10)      SIGV=SIG+SAB
      IF(INDINT.EQ.2)    GO  TO  139
      IF(NR.EQ.5)        GO  TO  42
      TEMP1=1.-SIGV/(31.41592*(DLK**2))
      DRND=RNDM(-1.)
      IF(DRND.LT.TEMP1)  GO  TO  9
                         GO  TO  42
  139 TEMP1=SIGV/31.41592
      BI2=BIM2(P1,P2)
      IF(BI2.GT.TEMP1)   GO  TO  9
                         GO  TO  42
   41 M=0
      TSH2=0.
c      IF(NCAS.GE.NCPRI)  write( *,*) ' M=0'
      RETURN
   42 CONTINUE
      M=1
      DTINT=TINT-TFRO1
      IF(NCAS.LT.NCPRI)  RETURN
      write(16,665) IJPA,NR,TAUK,TINT
 665  FORMAT(1X,'IJPA,NR,TAUK,TINT=',I5,I2,2(1X,F7.3))
      IF(NR.GE.2.AND.NR.LE.4)
     *write(16,666) NR,N1,(PMEMO(K,N1),K=1,9),(IMEMO(K,N1),K=1,5),
     *                 N2,(P2(K),K=1,9),(IP2(K),K=1,5)
 666  FORMAT(' NR,N1=',I1,I5,9(1X,F10.3),4I2,I8/
     *       '   +N2=',1X,I5,9(1X,F10.3),4I2,I8)
      IF(NR.EQ.1.OR.NR.EQ.5)
     *write(16,667) NR,N1,(P1(K),K=1,9),(IP1(K),K=1,5),
     *                 N2,(P2(K),K=1,9),(IP2(K),K=1,5)
 667  FORMAT(' NR,N1=',I1,I5,9(1X,F10.3),4I2,I8/
     *       '   +N2=',1X,I5,9(1X,F10.3),4I2,I8)
      RETURN
      END