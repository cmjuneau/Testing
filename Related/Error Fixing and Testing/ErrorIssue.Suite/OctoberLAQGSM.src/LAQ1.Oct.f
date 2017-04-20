
c          For  LAQGSM-MARS interface      KKG  04/09/07
c aa2k7g.f,  gengamn7.f, gadd7.f,    ->     laqgsm2007_1.f
c qgsm7.f                            ->     laqgsm2007_2.f
c coales07.f, gemaa07.f, preco07g.f, spegem7g.f, fermi07.f,
c gambetm7.f                         ->     laqgsm2007_3.f
c
c
c      =================================================================
c          laqgsm2007_1.f   subroutines :
c      =================================================================
c
C  Dimensions of arrays UP1 and UP2 are 66000 and 70000
c  Dimensions of arrays XC,YC,ZC,IZ,MPA,INT1 are 300,kkg 23.11.04    
C Fermi energies are changed via TF=TF*(AN1/ANUCL1)**(2./3.)
c see PANTNQ, POTENQ, TFERMIQ;                              20.06.95
c anti-baryon potential is added                         05.11.95
c introduced /hadr1/                                     16.03.96
c interactions in RIJ<DELTA are included in RAPID4       01.10.97
c pi+Y=>AK+N channel is added                            08.02.99
c scattering width in ELASLE is added for Deltas         27.04.99
c corr. B+B==>B+Y+K (DDNYK,DDDYK,DNNYK,DNDYK,NNNYK,NNDYK 14.05.01
c on unit 10 3 numbers are writen instead of IS9         21.05.01
c last corrections in PAULI2,PAULI3             26.06.02,13.12.04
c new angular distributions for n+p,p+p at T<2.0GeV(ELEX)22.10.03
c gamma as projectile is included                        28.10.03
c extended to Egamma up to 10 GeV, Dec. 2004
c corrections in RAPID4                                  14.02.05  
c last changes by KKG for interface with MARS code       23.03.07  
c  /isecgpn/ ==> /ixsgpn/, KINEMA ==> KINEMQ             28.03.07 
c
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
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE RAPID1(NA1,NA2,DELTA,P1,IP1,P2,IP2,N1,N2,TAU,DL1)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c  Determination of interaction  pair (N1,P1,IP1),(N2,P2,IP2)
c  where nucleon N1 is from  projetile and
c        nucleon N2 is from  target; the time TAU (fm/c) is calculated 
c  in observer's system (lab. or equal velocity system depending of
c  projectile VPR and target VTA nuclei velocities)    
c  Calls: PARTNQ, TFERMIQ, KINEMQ
c  Input: NA1, NA2 - mass numbers of projectile and target nuclei
c         DL1 - interaction parameter
c  Output: P1,IP1,P2,IP2,N1,N2,TAU 
c
      COMMON/RCOR/ RCOR
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/INTCEN/IPO,INT1(300),INT2(5999),INT3(5999),
     *              INT4(100000),IJPA,IST,JST,NRST
      COMMON/ACTIV/MPA(300),MYP(5999),MYT(5999),MYY(5999)
      COMMON/CENTER/XC(2,300),YC(2,300),ZC(2,300),IZ(2,300)
      COMMON /HOLPT/ PH1(3),PH2(3),RH1(3),RH2(3),EHOL1,EHOL2,TFP,TFT
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      DIMENSION P1(9),P2(9),IP1(5),IP2(5),RI(3),RJ(3),RIJ(3),VIJ(3),
     *PI(3),PJ(3),R1(3),R2(3),PLI(3),PLJ(3)
      DL1=1./(5.06*0.940*SQRT(GRE*GRE-1.))+DELTA
      GPR=1./SQRT(ABS(1.-VPR(1)**2-VPR(2)**2-VPR(3)**2))
      GTA=1./SQRT(ABS(1.-VTA(1)**2-VTA(2)**2-VTA(3)**2))
      DO   8  K=1,3
      PI(K)=0.940*GPR*VPR(K)
      PJ(K)=0.940*GTA*VTA(K)
      VIJ(K)=VPR(K)-VTA(K)
    8 CONTINUE
      VIJ2=VIJ(1)**2+VIJ(2)**2+VIJ(3)**2
      EI=0.940*GPR
      EJ=0.940*GTA
      SIJ=(EI+EJ)**2-(PI(1)+PJ(1))**2-(PI(2)+PJ(2))**2-
     -(PI(3)+PJ(3))**2
      TAU=-0.1
      DO  13  I=1,NA1
      IF(MPA(I).EQ.0)   GO  TO  13
      NI=0
      VRI=XC(1,I)*VPR(1)+YC(1,I)*VPR(2)+ZC(1,I)*VPR(3)
      RI(1)=XC(1,I)-VPR(1)*VRI*GPR/(GPR+1.)+RADP(1)
      RI(2)=YC(1,I)-VPR(2)*VRI*GPR/(GPR+1.)+RADP(2)
      RI(3)=ZC(1,I)-VPR(3)*VRI*GPR/(GPR+1.)+RADP(3)
      DO  12  J=1,NA2
      VRJ=XC(2,J)*VTA(1)+YC(2,J)*VTA(2)+ZC(2,J)*VTA(3)
      RJ(1)=XC(2,J)-VTA(1)*VRJ*GTA/(GTA+1.)+RADT(1)
      RJ(2)=YC(2,J)-VTA(2)*VRJ*GTA/(GTA+1.)+RADT(2)
      RJ(3)=ZC(2,J)-VTA(3)*VRJ*GTA/(GTA+1.)+RADT(3)
      DO  9  K=1,3
      RIJ(K)=RI(K)-RJ(K)
    9 CONTINUE
      TIJ=-(RIJ(1)*VIJ(1)+RIJ(2)*VIJ(2)+RIJ(3)*VIJ(3))/VIJ2
      IF(TIJ.LE.0.)      GO  TO  12
      SP=RIJ(1)*PI(1)+RIJ(2)*PI(2)+RIJ(3)*PI(3)
      ST=RIJ(1)*PJ(1)+RIJ(2)*PJ(2)+RIJ(3)*PJ(3)
      RIJ2=RIJ(1)**2+RIJ(2)**2+RIJ(3)**2
      B2=RIJ2+4.*(SIJ*SP*ST-(0.940*(SP+ST))**2)/
     /SIJ/(SIJ-4.*0.940**2)
      IF(B2.GT.(DL1**2))          GO  TO  12
      IF(IPO.EQ.1.AND.NRST.EQ.1.AND.I.EQ.IST.AND.J.EQ.JST)  GO  TO 12
      NI=NI+1
      IF(TAU.LT.0.)        GO  TO  10
      IF(TAU.LT.TIJ)       GO  TO  12
   10 TAU=TIJ
      N1=I
      N2=J
      DO  11  K=1,3
      R1(K)=RI(K)+TIJ*VPR(K)
      R2(K)=RJ(K)+TIJ*VTA(K)
   11 CONTINUE
   12 CONTINUE
       IF(NI.EQ.0)   MPA(I)=0
   13 CONTINUE
      IF(TAU.LE.0)  RETURN
      CALL  PARTNQ(1,N1,P1,IP1)
      TFP=TFERMIQ(P1(1),P1(2),P1(3),1)
      EHOL1=P1(8)
      CALL  PARTNQ(2,N2,P2,IP2)
      TFT=TFERMIQ(P2(1),P2(2),P2(3),2)
      EHOL2=P2(8)
      DO  14  K=1,3
      PI(K)=P1(3+K)
      PJ(K)=P2(3+K)
   14 CONTINUE
      CALL  KINEMQ(PI,VPR,PLI,CTI,STI,CFI,SFI,TLI,P1(9))
      CALL  KINEMQ(PJ,VTA,PLJ,CTJ,STJ,CFJ,SFJ,TLJ,P2(9))
      DO  15  K=1,3
      RH1(K)=P1(K)
      RH2(K)=P2(K)
      P1(K) =R1(K)
      P2(K) =R2(K)
      PH1(K)=-P1(3+K)
      PH2(K)=-P2(3+K)
      P1(3+K)=PLI(K)
      P2(3+K)=PLJ(K)
   15 CONTINUE
      P1(7)=0.
      P2(7)=0.
      P1(8)=SQRT(PLI(1)**2+PLI(2)**2+PLI(3)**2+0.940**2)-0.940
      P2(8)=SQRT(PLJ(1)**2+PLJ(2)**2+PLJ(3)**2+0.940**2)-0.940
      P1(9)=0.940
      P2(9)=0.940
      RETURN
               END
C * * * * * * *  * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE RAPID2(MV,NA1,DELTA,P1,IP1,N1,N2,TAU,DL)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c  Determination of interaction  pair (N1,P1,IP1),N2,
c  where nucleon N1 is number of cascade particle in memory array and
c  N2 is number of spectator nucleon from projectile nucleus; 
c  the time interval TAU (fm/c) is calculated in observer's system 
c  (lab. or equal velocity system in dependence of
c  projectile VPR and target VTA nuclei velocities)    
c  Calls: CENUM1, CINEMA
c  Input: MV - total number of cascade particles, 
c         NA1 -mass number of projectile nucleus
c         DL - interaction parameter
c  Output: P1,IP1,N1,N2,TAU 
c
c
      COMMON/RCOR/ RCOR
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/ACTIV/MPA(300),MYP(5999),MYT(5999),MYY(5999)
       COMMON/MEMORY/PMEMO(9,5999),IMEMO(5,5999)
      COMMON/CENTER/XC(2,300),YC(2,300),ZC(2,300),IZ(2,300)
      COMMON/CENPAR/NUKC(100)
       COMMON/TAUE/TPTE,TYPE,TYTE
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      COMMON/INTCEN/IPO,INT1(300),INT2(5999),INT3(5999),
     *              INT4(100000),IJPA,IST,JST,NRST
      COMMON/CSLID/CLIDER(5999)
      COMMON/IACT/ IACT/CVALON/ IVALON
      DIMENSION V0(3),PIN(9),IIN(5),CS(3),C(3),P1(9),IP1(5)
      K = 1
      TAU = -.1
      TAUL=-.1
      V0(1)=-VPR(1)
      V0(2)=-VPR(2)
      V0(3)=-VPR(3)
      GPR=1./SQRT(1.-VPR(1)**2-VPR(2)**2-VPR(3)**2)
      GG=GPR*GPR/(GPR+1.)
    9 IF(MYP(K)) 26,26,10
   10 VRK=(PMEMO(1,K)-RADP(1))*VPR(1)+(PMEMO(2,K)-RADP(2))*VPR(2)+
     +(PMEMO(3,K)-RADP(3))*VPR(3)
      XK0=PMEMO(1,K)-RADP(1)+VPR(1)*VRK*GG
      YK0=PMEMO(2,K)-RADP(2)+VPR(2)*VRK*GG
      ZK0=PMEMO(3,K)-RADP(3)+VPR(3)*VRK*GG
      CS(1)=PMEMO(4,K)
      CS(2)=PMEMO(5,K)
      CS(3)=PMEMO(6,K)
      CALL CINEMA(CS,V0,C,CT,ST,CF,SF,TL,PMEMO(9,K))
      PIN(4)=C(1)
      PIN(5)=C(2)
      PIN(6)=C(3)
      PIN(7)=PMEMO(7,K)
      PIN(9)=PMEMO(9,K)
      IIN(1)=IMEMO(1,K)
      IIN(2)=IMEMO(2,K)
      IIN(3)=IMEMO(3,K)
      IIN(4)=IMEMO(4,K)
      IIN(5)=IMEMO(5,K)
      CLID=CLIDER(K)
      PIN(1)=XK0
      PIN(2)=YK0
      PIN(3)=ZK0
      VLI=SQRT(TL*(TL+2.*PIN(9)))/(TL+PIN(9))
      PIN(8)=TL
      DLK=1./(5.06*SQRT(PIN(8)*(PIN(8)+2.*PIN(9))))+DELTA
      ICE=0
      IF(IPO.EQ.0)   GO  TO  12
      K1=INT2(K)
      GO  TO  16
   12 IF(ICE.NE.0)   GO  TO  13
      CALL  CENUM1(NA1,PIN,IIN,DLK,NC,K1,K2,0,1)
      ICE=1
   13 IF(NC.NE.0)   GO  TO  15
       MYP(K)=0
                  GO  TO  26
   15 K1=NUKC(ICE)
      ICE=ICE+1
      NC=NC-1
   16 CONTINUE
      DR=((XC(1,K1)-PIN(1))*CF+(YC(1,K1)-PIN(2))*SF)*ST
     *+(ZC(1,K1)-PIN(3))*CT
      IF(DR.LE.0.)    GO TO  12
      IF(IPO.EQ.1.AND.NRST.EQ.2.AND.K.EQ.IST.AND.K1.EQ.JST) GO TO 12
      PIN(1)=PIN(1)+DR*ST*CF
      PIN(2)=PIN(2)+DR*ST*SF
      PIN(3)=PIN(3)+DR*CT
      DX=PIN(1)-XK0
      DY=PIN(2)-YK0
      DZ=PIN(3)-ZK0
      DRK=SQRT(DX**2+DY**2+DZ**2)
      TAUK=DRK/VLI
      TAUK1=TAUK*GEP*(1.-VLI*(VEP(1)*ST*CF+VEP(2)*ST*SF+VEP(3)+CT))
      TAUKL=TAUK*GPR*(1.+VLI*(VPR(1)*ST*CF+VPR(2)*ST*SF+VPR(3)*CT))
c TAUK  is time interval in projectile rest system
c TAUKL is time interval in observer's rest system
c TAUK1 is time interval in equal velocity  system
C    !!!!!
      IF(TAUKL.LT.PIN(7).AND.CLID.LT.0.3)   GO  TO  12
      IF(TAUKL.LT.PIN(7).AND.IVALON.EQ.0)   GO  TO  12
C    !!!!!
      INT2(K)=K1
      IF(TAU)  23,23,22
ccc22 IF(TAU-TAUK)  26,26,23
   22 IF(TAUL-TAUKL)  26,26,23
   23 TAU=TAUK
      TYPE=TAUK1
      TAUL=TAUKL
      DO 24 L=1,9
      P1(L)=PIN(L)
   24 CONTINUE
      DO 25 L=1,5
      IP1(L)=IIN(L)
   25 CONTINUE
      DL=DLK
      N1=K
      N2=K1
   26 IF(K-MV) 27,28,28
   27 K=K+1
              GO TO 9
   28 CONTINUE
      TAU=TAUL
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE RAPID3(MV,NA2,DELTA,P1,IP1,N1,N2,TAU,DL)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c  Determination of interaction  pair (N1,P1,IP1),N2,
c  where nucleon N1 is number of cascade particle in memory array and
c  N2 is number of spectator nucleon from target nucleus; 
c  the time interval TAU (fm/c) is calculated in observer's system 
c  (lab. or equal velocity system in dependence of
c  projectile VPR and target VTA nuclei velocities)    
c  Calls: CENUM1, CINEMA
c  Input: MV - total number of cascade particles, 
c         NA2 -mass number of target nucleus,
c         DL - interaction parameter
c  Output: P1,IP1,N1,N2,TAU 
c
c
      COMMON/RCOR/ RCOR
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      COMMON/ACTIV/MPA(300),MYP(5999),MYT(5999),MYY(5999)
       COMMON/MEMORY/PMEMO(9,5999),IMEMO(5,5999)
      COMMON/CENTER/XC(2,300),YC(2,300),ZC(2,300),IZ(2,300)
      COMMON/CENPAR/NUKC(100)
       COMMON/TAUE/TPTE,TYPE,TYTE
      COMMON/INTCEN/IPO,INT1(300),INT2(5999),INT3(5999),
     *              INT4(100000),IJPA,IST,JST,NRST
      COMMON/CSLID/CLIDER(5999)
      COMMON/IACT/ IACT/CVALON/ IVALON
      DIMENSION V0(3),PIN(9),IIN(5),CS(3),C(3),P1(9),IP1(5)
      V0(1)=-VTA(1)
      V0(2)=-VTA(2)
      V0(3)=-VTA(3)
      GTA=1./SQRT(1.-VTA(1)**2-VTA(2)**2-VTA(3)**2)
      GG=GTA*GTA/(GTA+1.)
      K = 1
      TAU = -.1
      TAUL=-.1
    9 IF(MYT(K)) 26,26,10
   10 VRK=(PMEMO(1,K)-RADT(1))*VTA(1)+(PMEMO(2,K)-RADT(2))*VTA(2)+
     +(PMEMO(3,K)-RADT(3))*VTA(3)
      XK0=PMEMO(1,K)-RADT(1)+VTA(1)*VRK*GG
      YK0=PMEMO(2,K)-RADT(2)+VTA(2)*VRK*GG
      ZK0=PMEMO(3,K)-RADT(3)+VTA(3)*VRK*GG
      CS(1)=PMEMO(4,K)
      CS(2)=PMEMO(5,K)
      CS(3)=PMEMO(6,K)
      CALL  CINEMA(CS,V0,C,CT,ST,CF,SF,TL,PMEMO(9,K))
      PIN(4)=C(1)
      PIN(5)=C(2)
      PIN(6)=C(3)
      PIN(7)=PMEMO(7,K)
      PIN(9)=PMEMO(9,K)
      PIN(8)=TL
      IIN(1)=IMEMO(1,K)
      IIN(2)=IMEMO(2,K)
      IIN(3)=IMEMO(3,K)
      IIN(4)=IMEMO(4,K)
      IIN(5)=IMEMO(5,K)
      CLID=CLIDER(K)
      PIN(1)=XK0
      PIN(2)=YK0
      PIN(3)=ZK0
      VLI=SQRT(TL*(TL+2.*PIN(9)))/(TL+PIN(9))
      PIN(8)=TL
      DLK=1./(5.06*SQRT(PIN(8)*(PIN(8)+2.*PIN(9))))+DELTA
      ICE=0
      IF(IPO.EQ.0)   GO  TO  12
      K1=INT3(K)
      GO  TO  16
   12 IF(ICE.NE.0)   GO  TO  13
      CALL  CENUM1(NA2,PIN,IIN,DLK,NC,K1,K2,0,2)
      ICE=1
   13 IF(NC.NE.0)   GO  TO  15
       MYT(K)=0
                   GO  TO  26
   15 K1=NUKC(ICE)
      ICE=ICE+1
      NC=NC-1
   16 CONTINUE
      DR=((XC(2,K1)-PIN(1))*CF+(YC(2,K1)-PIN(2))*SF)*ST+
     *(ZC(2,K1)-PIN(3))*CT
      IF(DR.LE.0.)    GO TO  12
      IF(IPO.EQ.1.AND.NRST.EQ.3.AND.K.EQ.IST.AND.K1.EQ.JST) GO TO 12
      PIN(1)=PIN(1)+DR*ST*CF
      PIN(2)=PIN(2)+DR*ST*SF
      PIN(3)=PIN(3)+DR*CT
      DX=PIN(1)-XK0
      DY=PIN(2)-YK0
      DZ=PIN(3)-ZK0
      DRK=SQRT(DX**2+DY**2+DZ**2)
      TAUK=DRK/VLI
      TAUK1=TAUK*GET*(1.-VLI*(VET(1)*ST*CF+VET(2)*ST*SF+VET(3)*CT))
      TAUKL=TAUK*GTA*(1.+VLI*(VTA(1)*ST*CF+VTA(2)*ST*SF+VTA(3)*CT))
c TAUK  is time interval in target     rest system
c TAUKL is time interval in observer's rest system
c TAUK1 is time interval in equal velocity  system
C    !!!!!
      IF(TAUKL.LT.PIN(7).AND.CLID.LT.0.3)   GO  TO  12
      IF(TAUKL.LT.PIN(7).AND.IVALON.EQ.0)   GO  TO  12
C    !!!!!
      INT3(K)=K1
      IF(TAU)  23,23,22
ccc22 IF(TAU-TAUK) 26,26,23
   22 IF(TAUL-TAUKL) 26,26,23
   23 TAU=TAUK
      TYTE=TAUK1
      TAUL=TAUKL
      DO 24 L=1,9
      P1(L)=PIN(L)
   24 CONTINUE
      DO 25 L=1,5
      IP1(L)=IIN(L)
   25 CONTINUE
      DL=DLK
      N1=K
      N2=K1
   26 IF(K-MV) 27,28,28
   27 K=K+1
              GO TO 9
   28 CONTINUE
      TAU=TAUL
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  RAPIDD(MV,TAU,MD)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Find number MD of nearest in decay time of resonance particle 
c     from MV cascade particles
 
      COMMON /MEMORY/ PME(9,5999),IME(5,5999)
      COMMON/TLIMIT/TLIMIT
      COMMON/ACTIM/TINT
      COMMON/IDPME/ IDPME(5999)
      TAU=-.1
      IF(MV.LE.0)   RETURN
      DO  11  M=1,MV
      IF(IME(5,M).EQ.0)  GO  TO  11
      TAUM=DBLE(IME(5,M))/1000.
      IF(TAU.LT.0.)  GO  TO  10
      IF(TAUM.GT.TAU)  GO  TO  11
   10 TAU=TAUM
      MD=M
   11 CONTINUE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  RAPID4(MV,N1,N2,TAU)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c  Determination of  interaction pair (N1,N2) from MV cascade particles
c  the time interval TAU (fm/c) is calculated 
c  in observer's system (lab. or equal velocity system depending of
c  projectile VPR and target VTA nuclei velocities)    
c  Calls: B2IJ, SGIJ
c  Input: MV- total number of cascade particles
c  Output: N1,N2,TAU 
c  For first time all possible numbers of interacting pairs are 
c  memorized in INT4 with respective increasing ordered
c  time intervals in TAUIJ4
c
      COMMON/NCASCA/NCAS,NCPRI
      COMMON /MEMORY/ PME(9,5999),IME(5,5999)
      COMMON/INTCEN/IPO,INT1(300),INT2(5999),INT3(5999),
     *              INT4(100000),IJPA,IST,JST,NRST
      COMMON/TAUIJ/ TAUK,TPTS,TYPS,TYTS,TYYS,TIJ4(100000)
      COMMON/NUCOLL/ NUCOLL,MVCOLL(5999)
      COMMON/ACTIV/MPA(300),MYP(5999),MYT(5999),MYY(5999)
      COMMON/TLIMIT/TLIMIT
      COMMON/ACTIM/TINT
      COMMON/IPAUL/IP
      COMMON/INTCC/INTCC
      COMMON/IACT/ IACT/CVALON/ IVALON
      COMMON/CSLID/CLIDER(5999)
      DIMENSION RIJ(3),VIJ(3),NYY(5999)
      TAU=-.1
      IF(MV.LE.1)                                           RETURN
      IF(IP.EQ.0.OR.IPO.EQ.1)                               GO  TO  15
      DO  8  I=1,MV
    8 NYY(I)=MYY(I)
      DO  13 I=2,MV
      IF(IME(2,I).NE.0)                                     GO  TO  13
      CLIDI=CLIDER(I)
      EI=PME(8,I)+PME(9,I)
      INU=0
      IF( IME(2,I).EQ.0.AND.IME(3,I).EQ.0.AND.IME(4,I).EQ.1.
     *AND.IME(5,I).EQ.0)  INU=1
      JM=I-1
      DO  12  J=1,JM
      IF(IME(2,J).NE.0)                                     GO  TO  12
      IF(MVCOLL(I).EQ.MVCOLL(J))                            GO  TO  12
      IF(MYY(I).EQ.1)                                       GO  TO  9
      IF(MYY(J).EQ.1)                                       GO  TO  9
                                                            GO  TO  12
    9 JNU=0
      IF( IME(2,J).EQ.0.AND.IME(3,J).EQ.0.AND.IME(4,J).EQ.1.
     *AND.IME(5,J).EQ.0)  JNU=1
      IF((INU+JNU).LT.1.AND.INTCC.EQ.1)                     GO  TO  12
      CLIDJ=CLIDER(J)
      EJ=PME(8,J)+PME(9,J)
      DO  10  K=1,3
      RIJ(K)=PME(K,I)-PME(K,J)
      VIJ(K)=PME(K+3,I)/EI-PME(K+3,J)/EJ
   10 CONTINUE
      RVIJ=RIJ(1)*VIJ(1)+RIJ(2)*VIJ(2)+RIJ(3)*VIJ(3)
      IF(RVIJ.GE.0.)                                        GO  TO  12
      VIJ2=VIJ(1)**2+VIJ(2)**2+VIJ(3)**2
      TIJ=-RVIJ/VIJ2
c                ! 14.02.05    
      IF(TIJ.LT.0.001.OR.(TINT+TIJ).GT.TLIMIT)              GO  TO  12
      IF(NRST.eq.5.and.((IST.eq.I.and.JST.eq.J).or.
     &                  (JST.eq.I.and.IST.eq.J)).and.
     & ABS(TIJ-TAUK).le.0.001)           GO  TO  12       !!14.02.05
C   !!!!!
      IF(TIJ.LT.PME(7,I).AND.CLIDI.LT.0.3)                  GO  TO  12
      IF(TIJ.LT.PME(7,I).AND.IVALON.EQ.0)                   GO  TO  12
      IF(TIJ.LT.PME(7,J).AND.CLIDJ.LT.0.3)                  GO  TO  12
      IF(TIJ.LT.PME(7,J).AND.IVALON.EQ.0)                   GO  TO  12
C   !!!!!
      CALL  B2IJ(I,J,B2)
      CALL  SGIJ(I,J,SIG)
      SIGV=SIG
      IF(IVALON.NE.0)    SIGV=CLIDI*CLIDJ*SIG
      IF(B2.GT.(SIGV/31.41592))                             GO  TO  12
      IF(IJPA.LT.100000)                                     GO  TO  11
      write(16,100)
      write( *,100)
  100 FORMAT(1X,'IJPA>100000')
      GO  TO  12
   11 IJPA=IJPA+1
      TIJ4(IJPA)=TIJ
      INT4(IJPA)=10000*I+J
CC
C     IF(IME(4,I).EQ.0.AND.IME(4,J).EQ.0)
C    *CALL NAMIJ(IP,IPO,INT4(IJPA),I,J,IJPA,TIJ)
CC
      NYY(J)=NYY(J)+1
      NYY(I)=NYY(I)+1
   12 CONTINUE
   13 CONTINUE
      DO  14  I=1,MV
      IF(NYY(I).EQ.MYY(I).AND.MYY(I).LE.1)    NYY(I)=0
      IF(NYY(I).EQ.1.AND.MYY(I).EQ.0)         NYY(I)=2
      MYY(I)=NYY(I)
   14 CONTINUE
   15 CONTINUE
      IF(IJPA.LE.0)                                          RETURN
      IJM=1
      IJPA0=IJPA
      DO  19  IJP=1,IJPA0
   16 IF(IJP.GT.IJPA)                                      GO  TO  20
      TIJ=TIJ4(IJP)
      IF(TIJ.LT.0.000001)                                  GO  TO  17
      IJ=INT4(IJP)
      I=(IJ+1)/10000
      J=IJ-10000*I
      IF(I.EQ.J)                                           GO  TO  17
      IF(I.LT.1.OR.I.GT.5999.OR.J.LT.1.OR.J.GT.5999)
     *write(16,101) IJ,I,J
  101 FORMAT(' IJ,I,J=',3I10)
      IF(IPO.EQ.1.AND.NRST.EQ.5.AND.(I+J).EQ.(IST+JST).
     *                          AND.(I*J).EQ.(IST*JST))    GO  TO  17
      IF((TINT+TIJ).GT.TLIMIT)                             GO  TO  17
C   !!!!!
      IF(TIJ.LT.PME(7,I).AND.CLIDER(I).LT.0.3)             GO  TO  17
      IF(TIJ.LT.PME(7,I).AND.IVALON.EQ.0)                  GO  TO  17
      IF(TIJ.LT.PME(7,J).AND.CLIDER(J).LT.0.3)             GO  TO  17
      IF(TIJ.LT.PME(7,J).AND.IVALON.EQ.0)                  GO  TO  17
C   !!!!!
      IF(TAU.LT.0.)                                        GO  TO  18
      IF(TIJ.GT.TAU)                                       GO  TO  19
                                                           GO  TO  18
   17 TIJ4(IJP)=TIJ4(IJPA)
      INT4(IJP)=INT4(IJPA)
      IJPA=IJPA-1
      IF(IJPA.LE.0)                                        GO  TO  21
                                                           GO  TO  16
   18 TAU=TIJ
      N1=I
      N2=J
      IJM=IJP
   19 CONTINUE
   20 CONTINUE
C      IJ=INT4(IJM)
C      INT4(IJM)=INT4(IJPA)
C      TIJ4(IJM)=TIJ4(IJPA)
C      IJPA=IJPA-1
   21 IF(TAU.LT.0.)                                          RETURN
CC
C     IF(IME(4,N1).EQ.0.AND.IME(4,N2).EQ.0)
C    *CALL NAMIJ(IP,IPO,IJ,N1,N2,IJPA,TAU)
CC
      IF(INTCC.GE.2)                                         RETURN
      IF((IME(4,N1)+IME(4,N2)).EQ.2.OR.IME(4,N2).EQ.1)       RETURN
      M=N2
      N2=N1
      N1=M
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ERAIJ(N)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Erase the particle N from INT4, TIJ4 arrays
c
      COMMON/INTCEN/IPO,INT1(300),INT2(5999),INT3(5999),
     *              INT4(100000),IJPA,IST,JST,NRST
      COMMON/TAUIJ/ TAUK,TPTS,TYPS,TYTS,TYYS,TIJ4(100000)
      IF(IJPA.LE.0)          RETURN
      IJPA0=IJPA
      DO  12  IJP=1,IJPA0
   10 IF(IJP.GT.IJPA)        RETURN
      IJ=INT4(IJP)
      I=(IJ+1)/10000
      J=IJ-10000*I
      IF(N.EQ.I.OR.N.EQ.J)  GO  TO  11
      GO  TO  12
   11 TIJ4(IJP)=TIJ4(IJPA)
      INT4(IJP)=INT4(IJPA)
      IJPA=IJPA-1
      GO  TO  10
   12 CONTINUE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE REPIJ(N,M)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Replace the particle M to particle in INT4 array
c
      COMMON/INTCEN/IPO,INT1(300),INT2(5999),INT3(5999),
     *              INT4(100000),IJPA,IST,JST,NRST
      COMMON/TAUIJ/ TAUK,TPTS,TYPS,TYTS,TYYS,TIJ4(100000)
      IF(IJPA.LE.0)          RETURN
      DO  11  IJP=1,IJPA
      IJ=INT4(IJP)
      I=(IJ+1)/10000
      J=IJ-10000*I
      IF(M.EQ.I.OR.M.EQ.J)  GO  TO  10
      GO  TO  11
   10 IF(M.EQ.I)  I=N
      IF(M.EQ.J)  J=N
      IJ=10000*I+J
      INT4(IJP)=IJ
   11 CONTINUE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE NAMIJ(IP,IPO,IJ,I,J,IJPA,TAU)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Determines the name of interacting particles I,J
c     from memory array
c
      CHARACTER*8 PNI,PNJ
      COMMON /MEMORY/ PME(9,5999),IME(5,5999)
      COMMON /IDPME/ IDPME(5999)
      DIMENSION PI(9),IPI(5),PJ(9),IPJ(5)
      DO  10  K=1,9
      PI(K)=PME(K,I)
      PJ(K)=PME(K,J)
   10 CONTINUE
      DO  11  K=1,5
      IPI(K)=IME(K,I)
      IPJ(K)=IME(K,J)
   11 CONTINUE
      IDI=IDPME(I)
      CALL PANUID(IDI,IKI,PNI)
      IDJ=IDPME(J)
      CALL PANUID(IDJ,IKJ,PNJ)
C     IF(IPI(4).EQ.0.AND.IPJ(4).EQ.0)
C    *write(16,100) IP,IPO,PNI,PNJ,IJ,I,J,IJPA,TAU
      IF(IPI(4).EQ.0.AND.IPJ(4).EQ.0)
     *write( *,100) IP,IPO,PNI,PNJ,IJ,I,J,IJPA,TAU
100   FORMAT(11X,2I3,1X,2A8,I7,2I5,I6,F8.3)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  B2IJ(I,J,B2)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Calculates the impact parameter B2=B**2 of cascade particles I,J
c 
      COMMON /MEMORY/ PME(9,5999),IME(5,5999)
      DIMENSION R(3),P(3),RS(3),PS(3)
      EI=PME(8,I)+PME(9,I)
      EJ=PME(8,J)+PME(9,J)
      E=EI+EJ
      PR=0.
      PPI=0.
      PPJ=0.
      DO  10  K=1,3
      R(K)=PME(K,I)-PME(K,J)
      P(K)=PME(K+3,I)+PME(K+3,J)
      PR=PR+P(K)*R(K)
      PPI=PPI+P(K)*PME(K+3,I)
      PPJ=PPJ+P(K)*PME(K+3,J)
   10 CONTINUE
      U=SQRT(E**2-P(1)**2-P(2)**2-P(3)**2)
      EIS=(U**2+PME(9,I)**2-PME(9,J)**2)/(2.*U)
      EJS=U-EIS
      PRS=0.
      DO  11  K=1,3
      PS(K)=(PME(K+3,I)+P(K)/U*(PPI/(E+U)-EI)-
     -       PME(K+3,J)-P(K)/U*(PPJ/(E+U)-EJ))/2.
      RS(K)=R(K)+P(K)*PR/(E+U)/U
      PRS=PRS+PS(K)*RS(K)
   11 CONTINUE
      PS2=PS(1)**2+PS(2)**2+PS(3)**2
      RS2=RS(1)**2+RS(2)**2+RS(3)**2
      B2=RS2-(PRS**2)/PS2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION  BIM2(P1,P2)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Calculates the impact parameter BIM2=B**2 of particles P1,P2
c     Definition of P1(P2):
c                          P1(1) - x coordinate in observer's system
c                          P1(2) - y coordinate in observer's system
c                          P1(3) - z coordinate in observer's system
c                          P1(4) - x component of momentum (GeV/c)  
c                          P1(5) - y component of momentum (GeV/c)  
c                          P1(6) - z component of momentum (GeV/c)
c                          P1(7) - the maturity (formation) time (fm/c)  
c                          P1(8) - kinetic energy (GeV)
c                          P1(9) - mass (GeV)
c
      DIMENSION R(3),P(3),RS(3),PS(3),P1(9),P2(9)
      EI=P1(8)+P1(9)
      EJ=P2(8)+P2(9)
      E=EI+EJ
      PR=0.
      PPI=0.
      PPJ=0.
      DO  10  K=1,3
      R(K)=P1(K)-P2(K)
      P(K)=P1(K+3)+P2(K+3)
      PR=PR+P(K)*R(K)
      PPI=PPI+P(K)*P1(K+3)
      PPJ=PPJ+P(K)*P2(K+3)
   10 CONTINUE
      U=SQRT(E**2-P(1)**2-P(2)**2-P(3)**2)
      EIS=(U**2+P1(9)**2-P2(9)**2)/(2.*U)
      EJS=U-EIS
      PRS=0.
      DO  11  K=1,3
      PS(K)=(P1(K+3)+P(K)/U*(PPI/(E+U)-EI)-
     -       P2(K+3)-P(K)/U*(PPJ/(E+U)-EJ))/2.
      PS(K)=(P1(K+3)*(EJ+EJS)-P2(K+3)*(EI+EIS))/(E+U)
      RS(K)=R(K)+P(K)*PR/(E+U)/U
      PRS=PRS+PS(K)*RS(K)
   11 CONTINUE
      PS2=PS(1)**2+PS(2)**2+PS(3)**2
      RS2=RS(1)**2+RS(2)**2+RS(3)**2
      BIM2=RS2-(PRS**2)/PS2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  SGIJ(I,J,SG)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c   Calculates the total interaction cross section SG of particles I,J 
c
      CHARACTER*8 PNA1,PNA2
      COMMON/MEMORY/PME(9,5999),IME(5,5999)
      COMMON /IDPME/ IDPME(5999)
      DIMENSION P1(9),P2(9),IP1(5),IP2(5)
      DO  10  K=1,9
      P1(K)=PME(K,I)
      P2(K)=PME(K,J)
      IF(K.LE.5)  IP1(K)=IME(K,I)
      IF(K.LE.5)  IP2(K)=IME(K,J)
   10 CONTINUE
      ID1=IDPME(I)
      CALL PANUID(ID1,IK1,PNA1)
      ID2=IDPME(J)
      CALL PANUID(ID2,IK2,PNA2)
      PX1=P1(4)
      PY1=P1(5)
      PZ1=P1(6)
      AM1=P1(9)
      PX2=P2(4)
      PY2=P2(5)
      PZ2=P2(6)
      AM2=P2(9)
      CALL CROSEC(1,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      SG=SITO
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  UPACOV(MV,MVU,IU)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Packs in UP1 array all cascade particle which more no interact 
c     after time interval TLIMIT
c     
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      COMMON/ACTIM/TINT/TPROD/TPROD(5999)
      COMMON/PORIG/IORI(3,5999)
      COMMON/NUCOLL/ NUCOLL,MVCOLL(5999)
      COMMON/CSLID/CLIDER(5999)
      COMMON/BGLAB/BLAB,GLAB,KOBR
      COMMON/UPAC/UP1(66000),IPER
      COMMON/MEMORY/PM(9,5999),IM(5,5999)
      COMMON/ACTIV/MPA(300),MYP(5999),MYT(5999),MYY(5999)
      COMMON/NPIMI/NPIMI
      COMMON /NCASCA/ NCAS,NCPRI
      COMMON/INTCC/INTCC
      COMMON /IDPME/ IDPME(5999)
      COMMON/TLIMIT/TLIMIT
      COMMON /SORI/ SORI(5999),SSOR
      COMMON /HADR1/HADR1(4,5999),HADR2(4,5999),HADI1(4),HADI2(4)
      MV1=MV
                 IF(MV.EQ.0)   RETURN
      M=1
    9 CONTINUE
      IF(IU.EQ.1)   GO  TO  10
      IF(AN1.LT.1.) MYP(M)=0
      IF(AN2.LT.1.) MYT(M)=0
      IF(MYP(M).NE.0.OR.MYT(M).NE.0)                 GO TO 20
      IF(IM(5,M).NE.0)                               GO TO 20
      IF(INTCC.NE.0.AND.IM(2,M).EQ.0)                GO TO 20
   10 CONTINUE
      I=MVU+1
      IF((11*I).LT.66000)      GO  TO  201
      IPER=1
      write(16,'(15X,''ARRAY  UP1 IS EXCEEDED IN CASCAN''/)')
      write( *,'(15X,''ARRAY  UP1 IS EXCEEDED IN CASCAN''/)')
                                RETURN
  201 CONTINUE
      ENM=PM(8,M)+PM(9,M)
      UP1(11*I-10)=DBLE(IDPME(M))
      UP1(11*I- 9)=DBLE(IORI(1,M))
         ICODE=IORI(3,M)*10**4 + IABS(IORI(2,M))
         ICODE=ISIGN(ICODE,IORI(2,M))
      UP1(11*I- 8)=DBLE(ICODE)
      UP1(11*I- 7)=PM(9,M)
      IF(IM(2,M).NE.0)  UP1(11*I- 7)=DBLE(IM(2,M))
      UP1(11*I- 6)=PM(1,M)-PM(4,M)/ENM*(TINT-TPROD(M))
      UP1(11*I- 5)=PM(2,M)-PM(5,M)/ENM*(TINT-TPROD(M))
      UP1(11*I- 4)=PM(3,M)-PM(6,M)/ENM*(TINT-TPROD(M))
      UP1(11*I- 3)=PM(4,M)
      UP1(11*I- 2)=PM(5,M)
      UP1(11*I- 1)=PM(6,M)
      IF(KOBR.EQ.1)
     *UP1(11*I- 1)=-GLAB*(PM(6,M)-BLAB*(PM(8,M)+PM(9,M)))
      UP1(11*I   )=TPROD(M)    !!!  production time
*      UP1(11*I   )=SORI(M)
      IF(IM(4,M).EQ.0.AND.IM(1,M).EQ.-1.AND.IM(3,M).EQ.0)
     & 	NPIMI=NPIMI+1
      MVU=MVU+1
      MV1=MV1-1
      IF(NCAS.GE.NCPRI)
     *write(16,202) M,MVU,MV1,(PM(J,M),J=1,9),(IM(JJ,M),JJ=1,5)
  202 FORMAT(1X,'UPACOV: M=',I5,2X,'MVU=',I5,2X,'MV1=',I5/
     *1X,'==>',9(1X,F 8.3),2X,4I2,I15,'<==')
      IF(M.GT.MV1)   GO  TO  21
      DO  18  J=1,9
   18 PM(J,M)=PM(J,MV1+1)
      DO  19  J=1,5
   19 IM(J,M)=IM(J,MV1+1)
      MYP(M)=MYP(MV1+1)
      MYT(M)=MYT(MV1+1)
      MYY(M)=MYY(MV1+1)
      TPROD(M)=TPROD(MV1+1)
      IORI(1,M)=IORI(1,MV1+1)
      IORI(2,M)=IORI(2,MV1+1)
      IORI(3,M)=IORI(3,MV1+1)
      SORI(M)=SORI(MV1+1)
      do j=1,4
        HADR1(j,M)=HADR1(j,MV1+1)
        HADR2(j,M)=HADR2(j,MV1+1)
      enddo
      MVCOLL(M)=MVCOLL(MV1+1)
      CLIDER(M)=CLIDER(MV1+1)
      IDPME(M)=IDPME(MV1+1)
      CALL REPIJ(M,MV1+1)
      GO  TO  9
   20 IF(M.EQ.MV1)   GO  TO  21
      M=M+1
                GO  TO  9
   21 MV=MV1
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE RPTS(X,Y,Z,X1,Y1,Z1,NU)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c    
c     Transforms the observer system's coordinates X,Y,Z into
c     nucleus NU rest system   
c
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      DIMENSION V(3),R(3)
      DO  11  K=1,3
      IF(NU.EQ.2)  GO  TO  10
      V(K)=VPR(K)
      R(K)=RADP(K)
      GO  TO  11
   10 V(K)=VTA(K)
      R(K)=RADT(K)
   11 CONTINUE
      VR=(X-R(1))*V(1)+(Y-R(2))*V(2)+(Z-R(3))*V(3)
      G=1./SQRT(1.-V(1)**2-V(2)**2-V(3)**2)
      GG=G*G/(G+1.)
      X1=X-R(1)+V(1)*VR*GG
      Y1=Y-R(2)+V(2)*VR*GG
      Z1=Z-R(3)+V(3)*VR*GG
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION  TFERMIQ(X,Y,Z,NU)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Calculates the Fermi energy of nucleus NU  
c
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      TFERMIQ=0.
      IF(NU.EQ.1.AND.ANUCL1.LT.2.1)  RETURN
      IF(NU.EQ.2.AND.ANUCL2.LT.2.1)  RETURN
      IF(NU.EQ.1)  THEN
        AN=ANUCL1
        A=A1
        C=C1
        D=D1
        TF0=TF01
        RM=RM1
        Ares=AN1
      ELSE
        AN=ANUCL2
        A=A2
        C=C2
        D=D2
        TF0=TF02
        RM=RM2
        Ares=AN2
      ENDIF
      R=SQRT(X**2+Y**2+Z**2)
      IF((R/RM).GT.1.5)  RETURN
      IF(AN.le.10.)  THEN
        TFERMIQ=TF0*EXP(-(2./3.)*(R/A)**2) *(Ares/AN)**0.6666667
      ELSE
        TFERMIQ=TF0*( (1.+EXP(-A/C))/(1.+EXP((R-A)/C))
     &  *(Ares/AN) )**0.6666667            ! 09.04.95
      ENDIF
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  DISSIP(IP,T0)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Change the projectile and target momenta due of Coulomb interaction
c     Calls: COTRAN, CINEMA 
c
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      DIMENSION PPL(3),PTL(3)
      IF(IP.EQ.-1)  GO TO  14
      IF(T0.LT.5.0.AND.AN1.GT.1.1) CALL  COTRAN
      IF(IP)  14,14,9
    9 IF(AN1.LT.0.5)   GO  TO  11
      AMP=0.940*AN1
      CALL  CINEMA(PNUCL1,VPR,PPL,CT,ST,CF,SF,TLP,AMP)
      EP=TLP+AMP
      VPR(1)=PPL(1)/EP
      VPR(2)=PPL(2)/EP
      VPR(3)=PPL(3)/EP
      DO  10  K=1,3
   10 PNUCL1(K)=0.0
   11 IF(AN2.LT.0.5)  GO  TO  14
      AMT=0.940*AN2
      CALL  CINEMA(PNUCL2,VTA,PTL,CT,ST,CF,SF,TLT,AMT)
      ET=TLT+AMT
      VTA(1)=PTL(1)/ET
      VTA(2)=PTL(2)/ET
      VTA(3)=PTL(3)/ET
      DO  12  K=1,3
   12 PNUCL2(K)=0.0
   14 GP=1./SQRT(1.-VPR(1)**2-VPR(2)**2-VPR(3)**2)
      GT=1./SQRT(1.-VTA(1)**2-VTA(2)**2-VTA(3)**2)
      VTP=VPR(1)*VTA(1)+VPR(2)*VTA(2)+VPR(3)*VTA(3)
      DO  15  K=1,3
      VRE(K)=(VPR(K)+GT*VTA(K)*(VTP*GT/(GT+1.)-1.))/(1.-VTP)/GT
   15 VEV(K)=(VPR(K)*GP+VTA(K)*GT)/(GT+GP)
      GEV=1./SQRT(1.-VEV(1)**2-VEV(2)**2-VEV(3)**2)
      GRE=1./SQRT(1.-VRE(1)**2-VRE(2)**2-VRE(3)**2)
      SET=VEV(1)*VTA(1)+VEV(2)*VTA(2)+VEV(3)*VTA(3)
      SEP=VEV(1)*VPR(1)+VEV(2)*VPR(2)+VEV(3)*VPR(3)
      DO 115 K=1,3
      VET(K)=(VEV(K)+GT*VTA(K)*(SET*GT/(GT+1.)-1.))/(1.-SET)/GT
      VEP(K)=(VEV(K)+GP*VPR(K)*(SEP*GP/(GP+1.)-1.))/(1.-SEP)/GP
  115 CONTINUE
      GET=GEV*GT*(1.-SET)
      GEP=GEV*GP*(1.-SEP)
   16 CONTINUE
  200 FORMAT(5X,'VPR,VTA,GP,GT',2(3(F7.4,2X),5X),2(F8.4,2X))
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  CENUM1(NA,P,IP,DELTA,NC,K1,K2,IK,NU)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Calculates all possible interaction partners in cylinder of
c     radius DELTA along the direction of incoming particle (P,IP) in 
c     projectile (NU=1) or target (NU=2) nucleus; 
c     all partners are ordered in "time" and its numbers are stored in
c     NUKC array
c
      COMMON/CENTER/XC(2,300),YC(2,300),ZC(2,300),IZ(2,300)
      COMMON/CENPAR/NUKC(100)
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON /NCASCA/ NCAS,NCPRI
      DIMENSION   P(9),RIN1(3),RC1(3),IP(5),ZKC(100)
      PM=SQRT(P(4)**2+P(5)**2+P(6)**2)
      IF(IK)  9,9,16
    9 CONTINUE
      NC=0
      K=1
      K2=0
      DELTA2=DELTA**2
      R=RM2
      IF(NU.EQ.1)   R=RM1
      RIN1(3)=(P(1)*P(4)+P(2)*P(5)+P(3)*P(6))/PM
      DEV=SQRT(P(1)**2+P(2)**2+P(3)**2-RIN1(3)**2)
      IF(DEV.GT.(R+DELTA))   GO  TO  26
  104 RC1(3)=(XC(NU,K)*P(4)+YC(NU,K)*P(5)+ZC(NU,K)*P(6))/PM
      DR1=RC1(3)-RIN1(3)
      IF(DR1.LT.0.)   GO  TO  14
      IF(ABS(DR1).LT.0.0001)   GO  TO  14
      DEL2=(XC(NU,K)-P(1))**2+(YC(NU,K)-P(2))**2+(ZC(NU,K)-P(3))**2-
     -DR1**2
      IF(DEL2-DELTA2)   11,11,14
   11 Z1=RC1(3)
      KA=K
      IF(NC.EQ.0)  GO  TO  113
      J=1
   12 IF(Z1.GT.ZKC(J))  GO  TO  13
      Z=ZKC(J)
      ZKC(J)=Z1
      Z1=Z
      KC=NUKC(J)
      NUKC(J)=KA
      KA=KC
   13 IF(J.EQ.NC)   GO  TO  113
      J=J+1
                GO  TO  12
  113 IF(NC.EQ.100)   GO  TO  114
      NC=NC+1
      ZKC(NC)=Z1
      NUKC(NC)=KA
                       GO  TO  14
  114 K2=1
   14 IF(K-NA)  15,26,26
   15 K=K+1
      GO  TO  104
   16 IF(NA-1) 107,107,108
  107 K2=0
             GO TO 26
  108 IF(IP(4).EQ.1)  GO  TO  107
   17 J=1
   18 IF(K1-1)  20,19,20
   19 J=2
   20 R=SQRT((XC(NU,K1)-XC(NU,J))**2+(YC(NU,K1)-YC(NU,J))**2+
     *(ZC(NU,K1)-ZC(NU,J))**2)
      K2=J
   21 IF(J-K1)  22,24,22
   22 RJ=SQRT((XC(NU,K1)-XC(NU,J))**2+(YC(NU,K1)-YC(NU,J))**2+
     *(ZC(NU,K1)-ZC(NU,J))**2)
      IF(RJ-R)  23,24,24
   23 R=RJ
      K2=J
   24 IF(J-NA)  25,26,26
   25 J=J+1
                  GO  TO  21
c   kkg 11/05/03
   26 continue
      if(NCAS.ge.NCPRI)  then
c        write(*,*) ' CENUM1: NU,NC,K1,K2=',NU,NC,K1,K2
      endif
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  ABELQ(PIN,V,U,P1,P2,CT,FI,CM1,CM2)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     calculates the final momenta of elastic interacting particles
c     with masses CM1,CM2 of total cms energy U, velocity V 
c     Calls: CINEMA, ROTORQ;
c     scattering angles CT and FI are as input, P1,P2 - output
c
      DIMENSION  PIN(9),V(3),B(3),P1(3),P2(3),PINL(3),PINS(3)
      B(1)=-V(1)
      B(2)=-V(2)
      B(3)=-V(3)
      PINL(1)=PIN(4)
      PINL(2)=PIN(5)
      PINL(3)=PIN(6)
      CALL  CINEMA(PINL,B,PINS,CTS,STS,CFS,SFS,TS,PIN(9))
      ST=SQRT(1.-CT*CT)
      CF=COS(FI)
      SF=SIN(FI)
      E1=(U*U+CM1**2-CM2**2)/(2.*U)
      T1=ABS(E1-CM1)
      P1M=SQRT(T1*(T1+2.*CM1))
      P2(1)=P1M*ST*CF
      P2(2)=P1M*ST*SF
      P2(3)=P1M*CT
      CALL  ROTORQ(PINS,V,P2,P1)
      P2(1)=-P1(1)
      P2(2)=-P1(2)
      P2(3)=-P1(3)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  TINVUQ(P1,P2,U,V,TIN1)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     calculates the total cms energy, velocity of the system of two
c     particles (P1,P2), and the kinetic energy TIN1 of particle "1" in 
c     rest system of particle "2"
c
      DIMENSION  V(3),P1(9),P2(9)
c kkg 10/14/03
c      E1=P1(8)+P1(9)
c      E2=P2(8)+P2(9)
      E1=SQRT(P1(4)**2+P1(5)**2+P1(6)**2+P1(9)**2)
      E2=SQRT(P2(4)**2+P2(5)**2+P2(6)**2+P2(9)**2)
      V(1)=(P1(4)+P2(4))/(E1+E2)
      V(2)=(P1(5)+P2(5))/(E1+E2)
      V(3)=(P1(6)+P2(6))/(E1+E2)
      U=(E1+E2)*SQRT(1.-V(1)**2-V(2)**2-V(3)**2)
      TIN1=(U**2-(P1(9)+P2(9))**2)/(2.*P2(9))
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ABSORPQ(PARTIN,IPATIN,PARTNE,PAR1,NE,IE,MV,NP,V,U)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C     Block of calculation of outgoing characteristics 
c     in two nucleons (PARTNE,PAR1) absorption of incidente  
c     particle (pion) (PARTIN,IPATIN);
c     Calls: TINVUQ, ABELQ
c  
      REAL*8 MASN
      COMMON/MEMORY/PMEMO,IMEMO
	COMMON/TEFABS/ TEFABS,EHOL3,TFR3
      DIMENSION PARTIN(9),IPATIN(5),PARTNE(9),PAR1(9),
     *PAF(9),V(3),PIST(3),PNST(3),PMEMO(9,5999),IMEMO(5,5999)
      MASN=0.940
      IF(IPATIN(5).NE.0.AND.IPATIN(4).EQ.1)  GO  TO  13
      PAF(4)=PARTNE(4)+PAR1(4)
      PAF(5)=PARTNE(5)+PAR1(5)
      PAF(6)=PARTNE(6)+PAR1(6)
      PAF(9)=2.*PAR1(9)
      PAF(8)=SQRT(PAF(4)**2+PAF(5)**2+PAF(6)**2+PAF(9)**2)-PAF(9)
	TEFABS=PAF(8)
      CALL  TINVUQ(PARTIN,PAF,U,V,TIN1)
c  kkg 03/19/04
   13 CTST = 1.-2.*RNDM(-1.)
      FIST = 6.283185*RNDM(-1.)
   16 CALL ABELQ(PARTIN,V,U,PIST,PNST,CTST,FIST,MASN,MASN)
      IF(MV-5997)17,17,18
   18 NP = 0
      write(16,19)
   19 FORMAT (25X,'MEMORY IS EXCEEDED IN CASCAD')
      RETURN
   17 CONTINUE
      PMEMO(1,MV+3) = 0.
      PMEMO(2,MV+3) = 0.
      PMEMO(3,MV+3) = 0.
      PMEMO(4,MV+3) = PIST(1)
      PMEMO(5,MV+3) = PIST(2)
      PMEMO(6,MV+3) = PIST(3)
      PMEMO(7,MV+3) = 0.
      PMEMO(9,MV+3) = 0.940
      IMEMO(1,MV+3) = IE
      IMEMO(2,MV+3) = 0
      IMEMO(3,MV+3) = 0
      IMEMO(4,MV+3) = 1
      IMEMO(5,MV+3) = 0
      PMEMO(1,MV+1) = 0.
      PMEMO(2,MV+1) = 0.
      PMEMO(3,MV+1) = 0.
      PMEMO(4,MV+1) = -PIST(1)
      PMEMO(5,MV+1) = -PIST(2)
      PMEMO(6,MV+1) = -PIST(3)
      PMEMO(7,MV+1) = 0.
      PMEMO(9,MV+1) = 0.940
      IMEMO(1,MV+1) = NE
      IMEMO(2,MV+1) = 0
      IMEMO(3,MV+1) = 0
      IMEMO(4,MV+1) = 1
      IMEMO(5,MV+1) = 0
      NP = 2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PARTNQ(NU,I,P,IP)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Partner selection in the projectile (NU=1) or target (NU=2) 
c     nucleus.
c     P and IP refer to the partner nuclear particle which is 
c     potentially interacting with the cascade particle.
c
      COMMON/CENTER/XC(2,300),YC(2,300),ZC(2,300),IZ(2,300)
     */HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      DIMENSION P(9),IP(5)

      P(1)=XC(NU,I)
      P(2)=YC(NU,I)
      P(3)=ZC(NU,I)
      R = SQRT(P(1)**2+P(2)**2+P(3)**2)
c     R=RPOTEN(R)
      IF(NU-1) 10,10,11
   10 IF(ANUCL1-10.) 120,120,110
  110 TFR=TF01*(((1.+EXP(-A1/C1))/(1.+EXP((R-A1)/C1)))**0.6666667 )
      GO TO 12
  120 TFR=TF01*EXP(-(2./3.)*(R/A1)**2)
                                         GO TO 12
   11 IF(ANUCL2-10.) 121,121,111
  111 TFR = TF02*(((1.+EXP(-A2/C2))/(1.+EXP((R-A2)/C2)))**.6666667)
      GO TO 12
  121 TFR=TF02*EXP(-(2./3.)*(R/A2)**2)
                                         GO TO 12
   12 CONTINUE
c                    !!! 20.06.1995
      IF(NU.EQ.1) TFR=TFR*(AN1/ANUCL1)**(2./3.)
      IF(NU.EQ.2) TFR=TFR*(AN2/ANUCL2)**(2./3.)
c
      TN = TFR*(RNDM(-1.)**(2./3.))
      PN=SQRT(TN*(TN+1.88))
      CT=1.-2.*RNDM(-1.)
      ST=SQRT(1.-CT*CT)
      FI=6.283185*RNDM(-1.)
      CF=COS(FI)
      SF=SIN(FI)
      P(4)=PN*ST*CF
      P(5)=PN*ST*SF
      P(6)=PN*CT
      P(7)=0.
      P(8)=TN
      P(9)=0.94
      IP(2)=0
      IP(3)=0
      IP(4)=1
      IP(1)=IZ(NU,I)
      IP(5)=0
	IF(IP(1).GT.1)  then
          WRITE(16,*) ' PARTNQ: NU,I,IP(1)=',NU,I,IP(1)
        IP(1)=0
      END IF
      RETURN
               END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DIRECTQ (V,TIN1,MQ,MV,NP,PARTIN,KP,ITH)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C   DETERMINING OF DIRECTION OF SECONDARY PARTICLES MOTION.
c   Calls: COSTAQ JTYPAQ ROTORQ, CINEMA
c
      DIMENSION PMEMO(9,5999),IMEMO(5,5999),V(3),PAKV(3),PARTIN(9),
     1PLST(3),PL(3),PAKST(3),PIN(3)
     *,B(3),PIL(3)
      COMMON/MEMORY/PMEMO,IMEMO
      ND = 0
      KP = 0
      IF (MQ-1) 12,12,10
   10 IF (RNDM(-1.) -0.5) 15,11,11
   11 IF (RNDM(-1.) -0.5) 14,16,16
   12 IF (RNDM(-1.) -0.5) 14,13,13
   13 IF (RNDM(-1.) -0.5) 15,16,16
   14 M1 = 2
      M2 = 3
                        GO TO 17
   15 M1 = 1
      M2 = 3
                        GO TO 17
   16 M1 = 1
      M2 = 2
                        GO TO 17
   17 LAMBDA=1
      PAKV(1)=0.
      PAKV(2)=0.
      PAKV(3)=0.
      M1TEMP = MV+M1
      M2TEMP = MV+M2
   18 LTEMP = MV+LAMBDA
      IF (LAMBDA-M1) 19,21,19
   19 IF (LAMBDA-M2) 20,21,20
   20 JA = JTYPAQ(ITH,MQ,LAMBDA)
      CTL = COSTAQ(JA,TIN1)
      FL = 6.283185*RNDM(-1.)
      STL = SQRT (1.-CTL**2)
      TEMP1=COS(FL)
      TEMP2=SIN(FL)
      TEMP3=PMEMO(8,LTEMP)
      PMEMO(1,LTEMP) = 0.
      PMEMO(2,LTEMP) = 0.
      PMEMO(3,LTEMP) = 0.
      PMEMO(4,LTEMP) = TEMP3*STL*TEMP1
      PMEMO(5,LTEMP) = TEMP3*STL*TEMP2
      PMEMO(6,LTEMP) = TEMP3*CTL
      PMEMO(7,LTEMP) = 0.
      PAKV(1) = PAKV(1)+PMEMO(4,LTEMP)
      PAKV(2) = PAKV(2)+PMEMO(5,LTEMP)
      PAKV(3) = PAKV(3)+PMEMO(6,LTEMP)
   21 IF (LAMBDA-NP) 22,23,23
   22 LAMBDA = LAMBDA+1
                          GO TO 18
   23 PAKVM = SQRT (PAKV(1)**2+PAKV(2)**2+PAKV(3)**2)
      IF (NP-3) 25,24,25
   24 PIL(1)=PARTIN(4)
      PIL(2)=PARTIN(5)
      PIL(3)=PARTIN(6)
      B(1)=-V(1)
      B(2)=-V(2)
      B(3)=-V(3)
      CALL  CINEMA(PIL,B,PIN,CTS,STS,CFS,SFS,TIN,PARTIN(9))
      LAMBDA = 1
                   GO TO 27
   25 IF (PMEMO(8,M1TEMP)-PAKVM-PMEMO(8,M2TEMP)) 26,32,32
   26 IF (PMEMO(8,M1TEMP)-ABS(PAKVM-PMEMO(8,M2TEMP))) 32,32,24
   27 LTEMP = MV+LAMBDA
      IF (LAMBDA-M1) 29,34,29
   28 LAMBDA = LAMBDA+1
                          GO TO 27
   29 IF (LAMBDA-M2) 30,34,30
   30 PL(1) = PMEMO(4,LTEMP)
      PL(2) = PMEMO(5,LTEMP)
      PL(3) = PMEMO(6,LTEMP)
      CALL ROTORQ (PIN,V,PL,PLST)
      PMEMO(1,LTEMP) = 0.
      PMEMO(2,LTEMP) = 0.
      PMEMO(3,LTEMP) = 0.
      PMEMO(4,LTEMP) = PLST(1)
      PMEMO(5,LTEMP) = PLST(2)
      PMEMO(6,LTEMP) = PLST(3)
      PMEMO(7,LTEMP) = 0.
   34 IF (LAMBDA-NP) 28,31,31
   31 CALL ROTORQ (PIN,V,PAKV,PAKST)
      CTM1 = (PMEMO(8,M2TEMP)**2-PMEMO(8,M1TEMP)**2-PAKVM**2)/
     1(2.*PAKVM*PMEMO(8,M1TEMP))
      CTM2 = (PMEMO(8,M1TEMP)**2-PMEMO(8,M2TEMP)**2-PAKVM**2)/
     1(2.*PAKVM*PMEMO(8,M2TEMP))
      FM1 = 6.283185*RNDM(-1.)
      FM2 = 3.141592+FM1
      STM1 = SQRT (1.-CTM1**2)
      STM2 = SQRT (1.-CTM2**2)
      CFM1 = COS(FM1)
      SFM1 = SIN(FM1)
      CFM2 = COS(FM2)
      SFM2 = SIN(FM2)
      PL(1) = PMEMO(8,M1TEMP)*STM1*CFM1
      PL(2) = PMEMO(8,M1TEMP)*STM1*SFM1
      PL(3) = PMEMO(8,M1TEMP)*CTM1
      CALL ROTORQ (PAKST,V,PL,PLST)
      PMEMO(1,M1TEMP) = 0.
      PMEMO(2,M1TEMP) = 0.
      PMEMO(3,M1TEMP) = 0.
      PMEMO(4,M1TEMP) = PLST(1)
      PMEMO(5,M1TEMP) = PLST(2)
      PMEMO(6,M1TEMP) = PLST(3)
      PMEMO(7,M1TEMP) = 0.
      PL(1) = PMEMO(8,M2TEMP)*STM2*CFM2
      PL(2) = PMEMO(8,M2TEMP)*STM2*SFM2
      PL(3) = PMEMO(8,M2TEMP)*CTM2
      CALL ROTORQ (PAKST,V,PL,PLST)
      PMEMO(1,M2TEMP) = 0.
      PMEMO(2,M2TEMP) = 0.
      PMEMO(3,M2TEMP) = 0.
      PMEMO(4,M2TEMP) = PLST(1)
      PMEMO(5,M2TEMP) = PLST(2)
      PMEMO(6,M2TEMP) = PLST(3)
      PMEMO(7,M2TEMP) = 0.
      RETURN
   32 ND = ND+1
      IF (ND-100) 17,33,33
   33 KP = 2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PINPNQ(RN,R0X,R0Y,R0Z,MV,DELTA,NEL)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Calculation of entry point of particle into nucleus and
c     positions of target nucleons;
c     in the case of A+A interactions simulates positions of
c     nucleons in projectile nucleus (see /CENTER/ )
c     in the case of bremsstrahlung gamma + A interaction
c     energy of gamma is samplied from Schiff's spectra (1/E)
c     Calls: RXYZ,PANUN, IDPANUN, ROTNUC 
c     
      CHARACTER*8 PNA1
      character*60 FGAM
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
     */CENTER/XC(2,300),YC(2,300),ZC(2,300),IZ(2,300)
      COMMON/ACTIV/MPA(300),MYP(5999),MYT(5999),MYY(5999)
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
       COMMON/MEMORY/PMEMO(9,5999),IMEMO(5,5999)
C
      COMMON/XBMAX/XBMAX,IFIB0
      COMMON /STIN/ STIN,AMIN
      COMMON/CSLID/CLIDER(5999)
      COMMON /IDPME/ IDPME(5999)
      common /gbrems/ tgmin, tgmax, teqv, sxabs, ibrems 
      common /filegam/ FGAM 
      DIMENSION PARTIN(9),IPATIN(5)
      T=(2.*RN)**2
      IF(ANUCL1-1.1) 98,98,99
   98 RM1=0.
   99 IF(ANUCL1.EQ.2.)  RM1=2.158
      R12=RM1+RM2+DELTA+1./(5.06*SQRT(T0*(T0+1.88)))
      CALL  RXYZ(R12,R0X,R0Y,R0Z)
      IF(ANUCL1-1.) 100,100,101
  100 PMEMO(1,1)=R0X
      PMEMO(2,1)=R0Y
      PMEMO(3,1)=R0Z
      IMEMO(1,1)=ZNUCL1
      IMEMO(2,1)=0
c  kkg 10.12.04
      if(AMIN.le.0.00001)   then
        IMEMO(2,1)=1   ! gamma
        if(ibrems.eq.1)  then
c      sampling the energy of gamma from Schiff's spectra (1/E)
	    rdm = RNDM(-1.)
	    T0 = tgmin*(tgmax/tgmin)**rdm
        endif
      endif 
c  kkg 10.12.03
      IMEMO(3,1)=INTG(STIN)
      IMEMO(4,1)=INTG(ANUCL1)
      IF(ANUCL1.GT.0.5)  IMEMO(4,1)=1
      IMEMO(5,1)=0
      PMEMO(9,1)=AMIN
      PMEMO(4,1)=0.
      PMEMO(5,1)=0.
      PMEMO(6,1)=SQRT(T0*(T0+2.*PMEMO(9,1)))
      PMEMO(7,1)=0.
      PMEMO(8,1)=T0
      MYP(1)=0
      MYT(1)=1
      MYY(1)=1
      CLIDER(1)=1.
      DO  1 K=1,9
      PARTIN(K)=PMEMO(K,1)
      IF(K.LE.5)  IPATIN(K)=IMEMO(K,1)
    1 CONTINUE
c  kkg 29.10.03
      if(IMEMO(2,1).ne.0)  then
	  ID1=10                 ! gamma
      else
        CALL PANUN(PARTIN,IPATIN,IK1)
        CALL IDPANU(ID1,IK1,PNA1)
      endif
c  kkg 29.10.03
      IDPME(1)=ID1
      MV=1
      AN1=0.
      ZN1=0.
                         GO  TO  102
  101 IF(ANUCL1.GT.2.1)   GO  TO  104
      G=1.+T0/0.940
      CT1=1.-2.*RNDM(-1.)
      FI1=6.283185*RNDM(-1.)
      ST1=SQRT(1.-CT1*CT1)
      Z1=RM1*CT1/G
      Y1=RM1*ST1*SIN(FI1)
      X1=RM1*ST1*COS(FI1)
      PMEMO(1,1)=R0X+X1
      PMEMO(2,1)=R0Y+Y1
      PMEMO(3,1)=R0Z+Z1
      PMEMO(1,2)=R0X-X1
      PMEMO(2,2)=R0Y-Y1
      PMEMO(3,2)=R0Z-Z1
  103 AL=10.*RNDM(-1.)
      BR=RNDM(-1.)
      FBR=4.*AL*AL/(1.+AL*AL)**2
      IF(BR.GT.FBR)   GO  TO  103
      EPM=0.940*0.00223
      PM1=SQRT(EPM)*AL
      CT1=1.-RNDM(-1.)
      ST1=SQRT(1.-CT1*CT1)
      FI1=6.283185*RNDM(-1.)
      PX1=PM1*ST1*COS(FI1)
      PY1=PM1*ST1*SIN(FI1)
      PZ1=PM1*CT1
      PMEMO(4,1)=PX1
      PMEMO(4,2)=-PX1
      PMEMO(5,1)=PY1
      PMEMO(5,2)=-PY1
      PMEMO(6,1)=G*PZ1+SQRT((PM1*PM1+0.94*0.94)*(G*G-1.))
      PMEMO(6,2)=PMEMO(6,1)-2.*G*PZ1
      PMEMO(7,1)=0.
      PMEMO(7,2)=0.
      PMEMO(8,1)=G*SQRT(PM1*PM1+0.94*0.94)+PZ1*SQRT(G*G-1.)-0.940
      PMEMO(8,2)=PMEMO(8,1)-2.*PZ1*SQRT(G*G-1.)
      PMEMO(9,1)=0.940
      PMEMO(9,2)=0.940
      IMEMO(1,1)=0
      IDPME(1)=1220
      IMEMO(2,1)=0
      IMEMO(3,1)=0
      IMEMO(4,1)=1
      IMEMO(1,2)=1
      IDPME(2)=1120
      IMEMO(2,2)=0
      IMEMO(3,2)=0
      IMEMO(4,2)=1
      IMEMO(5,1)=0
      IMEMO(5,2)=0
      MYP(1)=0
      MYP(2)=0
      MYT(1)=1
      MYY(1)=1
      MYT(2)=1
      MYY(2)=1
      CLIDER(1)=1.
      CLIDER(2)=1.
      MV=2
      AN1=0
      ZN1=0
                         GO  TO  102
  104 CONTINUE
      IF(NEL.NE.0)  CALL  ROTNUC(1,ANUCL1)
      IF(NEL.NE.0)  GO TO 102
      IM=INTG(ANUCL1)
      NZ1=INTG(ZNUCL1)
      DO 21 I=1,IM
   11 B1=RNDM(-1.)
      RI=RM1*(B1**(1./3.))
      IF(ANUCL1-10.) 12,12,13
   12 FB=EXP(-(RI**2)/(A1**2))
                                 GO TO 14
   13 FB=(1.+EXP(-A1/C1))/(1.+EXP((RI-A1)/C1))
   14 DRND=RNDM(-1.)
      IF(DRND-FB) 15,15,11
   15 CT=1.-2.*RNDM(-1.)
      FI=6.283185*RNDM(-1.)
      ST=SQRT(1.-CT**2)
      RS=RI*ST
      XC(1,I)=RS*COS(FI)
      YC(1,I)=RS*SIN(FI)
      ZC(1,I)=RI*CT
      IF(I-NZ1) 16,16,17
   16 IZ(1,I)=1
                  GO TO 18
   17 IZ(1,I)=0
   18 IF(I-1) 221,221,19
   19 KM=I-1
      DO 20 K=1,KM
      IF((XC(1,I)-XC(1,K))**2+(YC(1,I)-YC(1,K))**2+(ZC(1,I)-ZC(1,K))**2
     *-T) 11,20,20
   20 CONTINUE
  221 CONTINUE
   21 CONTINUE
  102 IF(ANUCL2.GT.1.1)   GO  TO  114
      PMEMO(1,1)=0.
      PMEMO(2,1)=0.
      PMEMO(3,1)=1.E-20
      PMEMO(4,1)=0.
      PMEMO(5,1)=0.
      PMEMO(6,1)=1.E-20
      PMEMO(7,1)=0.
      PMEMO(8,1)=0.
      PMEMO(9,1)=0.940
      IF(ANUCL2.LT.0.5)  PMEMO(9,1)=0.140
      IMEMO(1,1)=ZNUCL2
      IDPME(2)=1120
      IF(ZNUCL2.LE.0.5) IDPME(2)=1220
      IMEMO(2,1)=0
      IMEMO(3,1)=0
      IMEMO(4,1)=ANUCL2
      IMEMO(5,1)=0
      CLIDER(1)=1.
      AN2=0.
      ZN2=0.
      MV=1
      MYP(1)=1
      MYT(1)=0
      MYY(1)=1
      RETURN
  114 CONTINUE
      IF(NEL.NE.0)  GO TO 105
      IM=INTG(ANUCL2)
      NZ2=INTG(ZNUCL2)
      DO 33 I=1,IM
   22 B1=RNDM(-1.)
   		RI=RM2*(B1**(1./3.))
      IF(ANUCL2-10.) 23,23,24
   23 FB=EXP(-(RI**2)/(A2**2))
                                 GO TO 25
   24 FB=(1.+EXP(-A2/C2))/(1.+EXP((RI-A2)/C2))
   25 DRND=RNDM(-1.)
      IF(DRND-FB) 26,26,22
   26 CT=1.-2.*RNDM(-1.)
      FI=6.283185*RNDM(-1.)
      ST=SQRT(1.-CT**2)
      RS=RI*ST
      XC(2,I)=RS*COS(FI)
      YC(2,I)=RS*SIN(FI)
      ZC(2,I)=RI*CT
      IF (I-NZ2) 27,27,28
   27 IZ(2,I)=1
                  GO TO 29
   28 IZ(2,I)=0
   29 IF(I-1) 33,33,30
   30 KM=I-1
      DO 32 K=1,KM
      IF((XC(2,I)-XC(2,K))**2+(YC(2,I)-YC(2,K))**2+(ZC(2,I)-ZC(2,K))**2
     *-T) 22,32,32
   32 CONTINUE
   33 CONTINUE
      RETURN
  105 CALL  ROTNUC(2,ANUCL2)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ROTNUC(N,AN)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     random rotation of projectile (N=1) or target nucleus (N=2)
c     changing the positions of nucleons
c     Calls: ROTORQ  
c
      COMMON/CENTER/XC(2,300),YC(2,300),ZC(2,300),IZ(2,300)
      DIMENSION A(3),B(3),R(3),RR(3)
      IF(AN.LT.2.1)   RETURN
      CT=1.-2.*RNDM(-1.)
      ST=SQRT(1.-CT*CT)
      FI=6.283185*RNDM(-1.)
      A(1)=ST*COS(FI)
      A(2)=ST*SIN(FI)
      A(3)=CT
      B(1)=0.
      B(2)=0.
      B(3)=1.
      IA=AN+0.1
      DO 10 I=1,IA
      R(1)=XC(N,I)
      R(2)=YC(N,I)
      R(3)=ZC(N,I)
      CALL ROTORQ(A,B,R,RR)
      XC(N,I)=RR(1)
      YC(N,I)=RR(2)
      ZC(N,I)=RR(3)
   10 CONTINUE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  KKG 09/05/08 : parameter IRET was added
      SUBROUTINE   TYPNEW(PARTIN,IPATIN,PARTNE,IPATNE,V,U,TIN1,
     *         SABS,MV,NP,NABS,PAR1,IPA1,N3,NU,N2,NA1,NA2,IRET)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c
c     Determine interaction type and calculate
c     secondary particles' characteristics.
c   Called by: CASCAW
c   Calls: CINEMA, PANUID, SLQEKQ, PARTNQ, CENUM1, CROSEG, ABSORPQ,
c          TFERMIQ, ELEXQ, PINETA, PINSEX, BBSEX, PIPIKK,AKANNI,
c          PIYSEX, BINELQ, HEINEN 
c           
      CHARACTER*4 ATYP(5)
      CHARACTER*8 DTYP(11),TDIAG,PNA1,PNA2,PNAK
      LOGICAL YESELA
      COMMON/YESELA/ YESELA
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/INTTYP/ ITYP
      COMMON/ITHEA/ITHEA(11)
      COMMON/LOWMIS/LOWMIS
      COMMON/IACT/ IACT
      COMMON /MEMORY/ PME(9,5999),IME(5,5999)
      COMMON /HELE/ IHELE
      COMMON /IDN12/ ID1,ID2
      COMMON /IDN120/ ID10,ID20
      COMMON /IDPME/ IDPME(5999)
      COMMON/COMELX/ SEL
      COMMON/COMCRO/ STO
      COMMON /DATA2/ PUD,PS1,SIGMA,CX2
      COMMON /NRAPI/ NRAPI
      COMMON /SORI/ SORI(5999),SSOR
      COMMON /HADR1/HADR1(4,5999),HADR2(4,5999),HADI1(4),HADI2(4)
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     &     ,VEP(3),VET(3),GEP,GET
      COMMON/TEFABS/ TEFABS,EHOL3,TFR3

      COMMON /COMENB/ ENBOU
      common/Uold/ Uold
      COMMON/STREXC/STREXC

      DIMENSION PARTIN(9),IPATIN(5),PARTNE(9),IPATNE(5),V(3),
     *PAR1(9),IPA1(5),PA(9),IPA(5),PS(3),PL(3)
      DATA ATYP/'ABS:','ELE:','BIN:','HEI:','SEX:'/
      DATA DTYP/'DIFTRI','PLANAR','UNCYLI','ELAST  ','ANNIH','DIFSMA',
     *          'CHAINS','BINAR ','??????','REGTRI ','DOUBDI'/
C
C
c  KKG 09/05/08 
      IRET=0
      Ntry=0 
      STREXC = 0.0                 !   kkg 19.02.07
      AM1=PARTIN(9)
      SSOR=U
c  for a 'projectile' hadron
      do j=1,3
        PS(j)=PARTIN(3+j)
      enddo
      if(NRAPI.eq.2)         then
        CALL  CINEMA(PS,VPR,PL,CTL,STL,CFL,SFL,TL,PARTIN(9))
      elseif(NRAPI.eq.3)     then
        CALL  CINEMA(PS,VTA,PL,CTL,STL,CFL,SFL,TL,PARTIN(9))
      else
        do j=1,3
          PL(j)=PARTIN(3+j)
        enddo
      endif
      HADI1(1)=PL(1)
      HADI1(2)=PL(2)
      HADI1(3)=PL(3)
      HADI1(4)=PARTIN(9)
c  for a 'target-partner' hadron
      do j=1,3
        PS(j)=PARTNE(3+j)
      enddo
      if(NRAPI.eq.2)         then
        CALL  CINEMA(PS,VPR,PL,CTL,STL,CFL,SFL,TL,PARTNE(9))
      elseif(NRAPI.eq.3)     then
        CALL  CINEMA(PS,VTA,PL,CTL,STL,CFL,SFL,TL,PARTNE(9))
      else
        do j=1,3
          PL(j)=PARTNE(3+j)
        enddo
      endif
      HADI2(1)=PL(1)
      HADI2(2)=PL(2)
      HADI2(3)=PL(3)
      HADI2(4)=PARTNE(9)
 11   CONTINUE
c  KKG 09/05/08 
      Ntry=Ntry + 1
      if(Ntry.gt.100)  then
        IRET=1
        return
      end if   
c
      ENBOUS=ENBOU
      Uolds =Uold
      if((IPATIN(4)+IPATNE(4)).eq.1)  then
        ENBOU=10.0
	Uold =1.3
      else
        ENBOU=2.0
	Uold =2.0
      endif  
c
      ID10=ID1
      ID20=ID2
      LOWMIS=0
      ITYP=0
      NABS=0
      YESELA=.TRUE.
      CALL PANUID(ID1,IK1,PNA1)
      CALL PANUID(ID2,IK2,PNA2)
      IF(NCAS.GE.NCPRI) THEN
        write(16,599) PNA1,ID1,PNA2,ID2,TIN1,SABS
        write( *,599) PNA1,ID1,PNA2,ID2,TIN1,SABS
  599 FORMAT(1X,A4,'(',I5,')','+',A4,'(',I5,')','  TIN1=',F7.3,
     &' SABS=',1pe9.3)
      ENDIF
      IF(TIN1.LE.0) THEN
c  kkg 10/14/03
c       write(16,*) 'PI=', (PARTIN(KK),KK=4,9)
c       write(16,*) 'PN=', (PARTNE(KK),KK=4,9)
      ENDIF
c       for gamma + N interaction, KKG, 12/10/04
   98 continue   
      if(ID1.eq.10)  then
	  STO = csgntot(IPATNE(1),TIN1,PARTNE(9))/1000.
	  betabs = SABS/(STO+SABS)
	  IF(RNDM(-1.).le.betabs) then
            go  to  101
	  else    
	    go  to  13
          endif
       endif
c
  112 CONTINUE
      IF(IK2.LT.37.OR.IK2.GT.38)           GO  TO  18
c  kkg 10/28/03
      IF(IK1.EQ.1.OR.IK1.EQ.2.OR.IK1.EQ.7.or.ID1.eq.10)
     &                                     GO  TO  99
      IF(IK1.LT.37.OR.IK1.GT.38)           GO  TO  18
                                           GO  TO  13
   99 CONTINUE
      CALL SLQEKQ (L,MS,MQ,KSI,ME,IPATIN,IPATNE)
      STO=CROSEG(L,MS,MQ,KSI,0,TIN1,AM1,IPATIN(5))
      BETABS=SABS/(STO+SABS)
      DRND=RNDM(-1.)
      IF(DRND-BETABS) 101,101,13
  101 CONTINUE
      DL=0.
      IF(NU.EQ.1)  CALL  CENUM1(NA1,PARTIN,IPATIN,DL,NC,N2,N3,N2,1)
      IF(NU.EQ.2)  CALL  CENUM1(NA2,PARTIN,IPATIN,DL,NC,N2,N3,N2,2)
      IF(N3.EQ.0)   then
	      write(*,*) ' TYPNEW: N3=0,ID1=',ID1
	      GO  TO  13
      endif  
      CALL  PARTNQ(NU,N3,PAR1,IPA1)
c  kkg 03/19/04   only gamma +(pn) !
c     if(ID1.eq.10.and.(IPATNE(1)+IPA1(1)).ne.1)  go to 13
      TFR3=TFERMIQ(PAR1(1),PAR1(2),PAR1(3),NU)
	EHOL3=PAR1(8)
      PAR1(9) = PARTNE(9)
      IE2=IPATNE(1)
      IE3=IPA1(1)
      IE1=IPATIN(1)
  104 IF(IE1) 108,103,105
  103 NE1=IE2
      NE2=IE3
                      GO TO 14
  105 IF(IE2+IE3-1) 107,106,13
  106 NE1=1
      NE2=1
                      GO TO 14
  107 NE1=1
      NE2=0
                      GO TO 14
  108 IF(IE2+IE3-1) 13,109,110
  109 NE1=0
      NE2=0
                      GO TO 14
  110 NE1=0
      NE2=1
   14 CONTINUE
      CALL ABSORPQ(PARTIN,IPATIN,PARTNE,PAR1,NE1,NE2,MV,NP,V,U)
      NABS = 1
      ITYP=1
      GO  TO  27
   13 CONTINUE
c      for gamma + N ==> hadrons, KKG, 12/10/04   
      if(ID1.eq.10)  then
      IF(NCAS.GE.NCPRI) write(*,*) ' to gntoh: TIN1=',TIN1
	  call  gntoh(V,U,TIN1,PARTIN,IPATIN,PARTNE,IPATNE,MV,NP)
      IF(NCAS.GE.NCPRI) write(*,*) ' from gntoh: NP=',NP
	    if(np.lt.2)  go to 101
          ITYP=5
	    NIN=0
	  go  to  27
      endif
c
C    ****** TIN1=4.5  IS INTRODUCED ******* Kostya
c      IF(IACT.GE.2.AND.TIN1.GT.4.5)    GO  TO  18
       IF(IACT.GE.2.AND.TIN1.GT.0.8)    GO  TO  18   ! KKG 21.11.07
c      IF(IACT.GE.2.AND.U.GT.Uold)        GO  TO  18 ! KKG 04.03.07

      CALL SLQEKQ (L,MS,MQ,KSI,ME,IPATIN,IPATNE)
      STO=CROSEG(L,MS,MQ,KSI,0,TIN1,AM1,IPATIN(5))
      SEL=CROSEG(L,MS,MQ,KSI,1,TIN1,AM1,IPATIN(5))
      SEX=CROSEG(L,MS,MQ,KSI,2,TIN1,AM1,IPATIN(5))
      BETAEL=(SEL+SEX)/STO
      DRND=RNDM(-1.)
      IF (DRND-BETAEL) 16,16,15
   16 CONTINUE
      CALL ELEXQ(V,U,TIN1,PARTIN,IPATIN,IPATNE,MV,NP,L,MS,
     *MQ,KSI,ME)
      ITYP=2
      GO  TO  27
   15 CONTINUE
c      IF(TIN1.LE.5.0d0)  GO  TO  20
      IF(TIN1.LE.0.8d0)  GO  TO  20      ! KKG 21.11.07

      YESELA=.FALSE.
      GO  TO  18
   20 CONTINUE
C
      YESELA=.FALSE.

C     kkg  11/21/07
      IETA1=0
      IETA2=0
      ISEX1=0
      ISEX2=0
      IBBKAK=0
      IPIBKAK=0

      CALL PINETA(ID1,ID2,V,U,PARTIN,IPATIN,PARTNE,IPATNE,
     &     MV,NP,IETA1)
      if((ID1.eq.1120.and.ID2.eq.1120).or.
     &   (ID1.eq.1220.and.ID2.eq.1220).or. 
     &   (ID1.eq.1120.and.ID2.eq.1220).or. 
     &   (ID1.eq.1220.and.ID2.eq.1120).and.U.le.Uold) ! KKG 04.03.07 
     &CALL PNETA(ID1,ID2,V,U,TIN1,PARTIN,PARTNE,MV,NP,IETA2)

      if(IETA1.eq.0) 
     & CALL PINSEX(ID1,ID2,V,U,PARTIN,IPATIN,PARTNE,IPATNE,
     &            MV,NP,ISEX1)
      IF(IPATIN(4).EQ.1.AND.IPATIN(3).EQ.0.AND.
     &   IPATNE(4).EQ.1.AND.IPATNE(3).EQ.0.and.U.le.Uold.and.
     &IETA2.eq.0) CALL BBSEX(ID1,ID2,V,U,PARTIN,PARTNE,MV,NP,ISEX2)
      if(IETA2.eq.0.and.ISEX2.eq.0.and.
     &     IPATIN(4).eq.1.and.IPATNE(4).eq.1)
     &call  BBKAK(ID1,ID2,V,U,PARTIN,PARTNE,IPATIN,IPATNE,MV,NP,
     & IBBKAK)
       if(IBBKAK.ne.0)   then
             ITYP=5
             NIN=0
             go to  27
       endif
       if(ISEX1.eq.0.and.IETA1.eq.0)
     & call  PIBKAK(ID1,ID2,V,U,PARTIN,PARTNE,IPATIN,IPATNE,MV,NP,
     & IPIBKAK)
       if(IPIBKAK.ne.0)   then
         ITYP=5
         NIN=0
         go to  27
       endif
c
      IF(ISEX1.NE.0.OR.ISEX2.NE.0.OR.IETA1.NE.0.or.IETA2.ne.0)
     &                THEN
	ITYP=5
	NIN=0
      ELSE
        CALL BINELQ(PARTIN,IPATIN,IPATNE,L,MS,MQ,KSI,ME,V,U,TIN1,
     *MV,NP,NIN)
        ITYP=3
        IF (NIN.NE.0)   GO  TO  11
      ENDIF
      GO  TO   27
   18 CONTINUE
C
C     kkg  07/25/06
      IKAK=0
      IETA1=0
      IETA2=0
      ISEX1=0
      ISEX2=0
      ISEX3=0
      ISEX4=0
      IETAN=0
      IBBKAK=0
      IPIBKAK=0

      CALL PIPIKK(ID1,ID2,V,U,PARTIN,IPATIN,PARTNE,IPATNE,
     &            MV,NP,IKAK)
      CALL PINETA(ID1,ID2,V,U,PARTIN,IPATIN,PARTNE,IPATNE,
     &            MV,NP,IETA1)
      if(IETA1.eq.0)
     & CALL PINSEX(ID1,ID2,V,U,PARTIN,IPATIN,PARTNE,IPATNE,
     &            MV,NP,ISEX1)
      if((ID1.eq.1120.and.ID2.eq.1120).or.
     &   (ID1.eq.1220.and.ID2.eq.1220).or. 
     &   (ID1.eq.1120.and.ID2.eq.1220).or. 
     &   (ID1.eq.1220.and.ID2.eq.1120).and.U.le.Uold) ! KKG 04.03.07
     &CALL PNETA(ID1,ID2,V,U,TIN1,PARTIN,PARTNE,MV,NP,IETA2)
      if(IETA1.eq.0.and.ISEX1.eq.0)
     &CALL PIYSEX(ID1,ID2,V,U,PARTIN,PARTNE,MV,NP,ISEX2)
	IF(IPATIN(4).EQ.1.AND.IPATIN(3).EQ.0.AND.
c    &   IPATNE(4).EQ.1.AND.IPATNE(3).EQ.0.and.U.le.Uold.
     &   IPATNE(4).EQ.1.AND.IPATNE(3).EQ.0.and.U.le.2.0d0.
     & and.IETA2.eq.0)                     !21.11.07
     &CALL BBSEX(ID1,ID2,V,U,PARTIN,PARTNE,MV,NP,ISEX3)
c
      CALL AKANNI(ID1,ID2,V,U,PARTIN,MV,NP,ISEX4)
c   kkg  20.09.06
      if(ID1.eq.220) CALL ETAN(ID1,ID2,V,U,PARTIN,IPATIN,
     &    PARTNE,IPATNE,MV,NP,IETAN)
      if(IETA2.eq.0.and.ISEX3.eq.0.and.
     &     IPATIN(4).eq.1.and.IPATNE(4).eq.1)
     &call  BBKAK(ID1,ID2,V,U,PARTIN,PARTNE,IPATIN,IPATNE,MV,NP,
     & IBBKAK)
       if(IBBKAK.ne.0)   then
             ITYP=5
             NIN=0
             go to  27
       endif
       if(ISEX1.eq.0.and.IETA1.eq.0)
     & call  PIBKAK(ID1,ID2,V,U,PARTIN,PARTNE,IPATIN,IPATNE,MV,NP,
     & IPIBKAK)
       if(IPIBKAK.ne.0)   then
         ITYP=5
         NIN=0
         go to  27
       endif
      IF(IKAK.NE.0.or.IETA1.ne.0.or.IETA2.ne.0.or.
     &   ISEX1.ne.0.or.ISEX2.ne.0.or.ISEX3.ne.0.or.
     &   ISEX4.ne.0.or.IETAN.ne.0) THEN
      
         ITYP=5
         NIN=0
      ELSE
C                            12.03.94
        PS1=0.750
C
        CALL  HEINEN(PARTIN,IPATIN,PARTNE,IPATNE,MV,NP,NIN)
        ITYP=4
c
        DO  19  I=1,11
   19   IF(ITHEA(I).NE.0) TDIAG=DTYP(I)
        IF(NIN.NE.0)  GO  TO  11
      ENDIF
   27 CONTINUE
C
      ENBOU=ENBOUS
      Uold =Uolds 
c
      IF(ITYP.EQ.0) GO  TO  603
      IF(NCAS.GE.NCPRI.OR.LOWMIS.NE.0)  THEN
c      IF(ITYP.NE.4)
c     *write(16,601) ATYP(ITYP),      PARTIN,IPATIN,PARTNE,IPATNE
c	IF(ITYP.EQ.1) WRITE(16,600)    PAR1,IPA1
c      IF(ITYP.EQ.4.AND.IHELE.EQ.1)
c     *write(16,611) ATYP(ITYP),TDIAG,PARTIN,IPATIN,PARTNE,IPATNE
c      IF(ITYP.EQ.4.AND.IHELE.EQ.2)
c     *write(16,612) ATYP(ITYP),TDIAG,PARTIN,IPATIN,PARTNE,IPATNE
      DO    K=1,NP
      M=MV+K
      IF(NP.EQ.2.AND.K.EQ.2.AND.(IABS(IPATIN(4))+IABS(IPATNE(4))).GT.0)
     &  M=M+1
c      write(16,602) M,(PME(I,M),I=4,6),PME(9,M),(IME(J,M),J=1,5)
      ENDDO
	ENDIF
C
      IF(ITYP.EQ.4)  GO  TO  3
c   kkg 12/13/04
	if(ID1.eq.10.and.ITYP.eq.5)  go  to 3
      DO  2  MM=1,NP
      M=MM
      IF(NP.EQ.2.AND.M.EQ.2.AND.(IABS(IPATIN(4))+IABS(IPATNE(4))).GT.0)
     &  M=3
      DO  1 K=1,9
      PA(K)=PME(K,MV+M)
      IF(K.LE.5)  IPA(K)=IME(K,MV+M)
    1 CONTINUE
      CALL PANUN(PA,IPA,IKK)
      CALL IDPANU(IDK,IKK,PNAK)
      IDPME(MV+M)=IDK
    2 CONTINUE
    3 CONTINUE
      RETURN
  600 FORMAT(1X,4X,9(1X,F10.3),4I2,I10)
  601 FORMAT(1X,A4,9(1X,F10.3),4I2,I10/
     *          5X,9(1X,F10.3),4I2,I10)
  611 FORMAT(1X,A4,1X,A8,'(LE)'/
     *          5X,9(1X,F10.3),4I2,I10/
     *          5X,9(1X,F10.3),4I2,I10)
  612 FORMAT(1X,A4,1X,A8,'(HE)'/
     *          5X,9(1X,F10.3),4I2,I10/
     *          5X,9(1X,F10.3),4I2,I10)
  602 FORMAT(1X,I3,4(1X,F10.3),4I2,I10)
  603 write(16,604)
  604 FORMAT(2X,'TYPNEW: ITYP=0 !!!!!'/)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  AKANNI(ID1,ID2,V,U,PIN,MV,NP,ISEX)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     charge exchange and annihilation of K- or AK0
c     Calls: ABELQ
c
      REAL*8 MPI,MN,MK,ML,MS,M1,M2
      COMMON/MEMORY/PME(9,5999),IME(5,5999)
      DIMENSION V(3),PIN(9),P1(3),P2(3)
      DATA  MPI/0.140/,MN/0.940/,MK/0.494/,ML/1.116/,MS/1.189/
      ISEX=0
      IR=0
      IF((ID1.EQ.-130.AND.ID2.EQ.1120).OR.(ID1.EQ.1120.AND.ID2.EQ.-130))
     &      IR=1  ! K- + p
      IF((ID1.EQ.-130.AND.ID2.EQ.1220).OR.(ID1.EQ.1220.AND.ID2.EQ.-130))
     &      IR=2  ! K- + n
      IF((ID1.EQ.-230.AND.ID2.EQ.1120).OR.(ID1.EQ.1120.AND.ID2.EQ.-230))
     &      IR=3  ! AK0 + p
      IF((ID1.EQ.-230.AND.ID2.EQ.1220).OR.(ID1.EQ.1220.AND.ID2.EQ.-230))
     &      IR=4  ! AK0 + n
      IF(IR.EQ.0)        RETURN
C
      IF(U.LE.(MK+MN))   RETURN
      EK0=(U**2-MN**2-MK**2)/(2.*MN)
      PK0=SQRT(EK0**2-MK**2)
C   Calculation of K- + p=>L + pi0 cross section
c   param. of G.Q.Li et al NP A625(1997)342
      IF(PK0.LE.0.6)        then
        sla=1.205*PK0**(-1.428)
      ELSEIF(PK0.LE.1.0)    then
	sla=3.5*PK0**0.659
      ELSE
	sla=3.5*PK0**(-3.97)
      ENDIF
C   Calculation of K- + p=>S0 + pi0 cross section
c   param. of G.Q.Li et al NP A625(1997)342
      IF(PK0.LE.0.345)      then
        ssi=0.624*PK0**(-1.830)
      ELSEIF(PK0.LE.0.425)  then
	ssi=0.0138/((PK0-0.385)**2+0.0017)
      ELSE
	ssi=0.7*PK0**(-2.09)
      ENDIF
C   Calculation of K- + p total,elastic,
c   and charge exchange cross sections
c   param. of G.Q.Li et al NP A625(1997)342
      IF(PK0.le.0.35)      then
        sto=23.5*PK0**(-1.04)
      ELSEIF(PK0.le.0.46)  then 
        sto=0.504/((PK0-0.39)**2+0.0056)
      ELSEIF(PK0.le.1.05)  then       
        sto=181.9*(PK0-0.75)**2+34.0
      ELSE
        sto=55.2*PK0**(-1.85)
      ENDIF
      IF(PK0.le.0.70)      then
        sel=11.2*PK0**(-0.986)
      ELSE
        sel=5.0/((PK0-0.95)**2+0.25)
      ENDIF
      IF(PK0.le.0.35)      then
        sex=1.813*PK0**(-1.14)
      ELSEIF(PK0.le.0.43)  then 
        sex=0.0192/((PK0-0.39)**2+0.0016)
      ELSE       
        sex=15.9/((PK0-0.9)**2+2.65)
      ENDIF
      san=sla+ssi
      ban=(san+sex)/sto
      IF(RNDM(-1.).gt.ban)  RETURN
      bex=sex/(sex+san)
      IF(RNDM(-1.).le.bex)  then
c  charge exchange K- p=>AK0 n or AK0 n=> K- p
        IS1=-1
        IS2=0
        M1=MK 
        M2=MN
        if(IR.eq.1)      then 
c K- p=> AK0 n            
          IE1=0
          IE2=0
        elseif(IR.eq.4)  then
c AK0 n=> K- p            
          IE1=-1
c  kkg 07/26/06
c         IE2=0
          IE2=1
        else
          return
        endif
      ELSE
c  annihilation: AK+N==>PI+Y
        IS1=0
        IS2=-1
        M1=MPI
        bla=sla/(sla+ssi)
        if(RNDM(-1.).le.bla)  then
          M2=ML
          IE2=0
          if(IR.eq.1)      then
c K- p=> pi0 L            
            IE1=0
          elseif(IR.eq.2)  then
c K- n=> pi- L            
            IE1=-1
          elseif(IR.eq.3)  then
c AK0 p=> pi+ L            
            IE1=1
          else
c AK0 n=> pi0 L            
            IE1=0
          endif
        else
          M2=MS
          br=RNDM(-1.)
          if(IR.eq.1)      then
            if(br.le.0.333333)     then
c K- p=> pi0 S0
              IE1=0
              IE2=0 
            elseif(br.ge.0.666667) then  
c K- p=> pi- S+
              IE1=-1
              IE2= 1
            else
c K- p=> pi+ S-
              IE1= 1
              IE2=-1
            endif
          elseif(IR.eq.2)  then
            if(br.le.0.5)          then
c K- n=> pi- S0
              IE1=-1
              IE2= 0
            else
c K- n=> pi0 S-
              IE1= 0
              IE2=-1
            endif
          elseif(IR.eq.3)  then
            if(br.le.0.5)          then
c AK0 p=> pi0 S+
              IE1= 0
              IE2= 1
            else
c AK0 p=> pi+ S0
              IE1= 0
              IE2= 1
            endif
          else
            if(br.le.0.5)          then
c AK0 n=> pi0 S0
              IE1= 0
              IE2= 0
            else
c AK0 n=> pi- S+
              IE1=-1
              IE2= 1
            endif
          endif
        endif
      ENDIF  
 
c
      CT=2.*RNDM(-1.)-1.
      FI=6.283185*RNDM(-1.)
      CALL  ABELQ(PIN,V,U,P1,P2,CT,FI,M1,M2)
      PME(1,MV+3)=0.
      PME(2,MV+3)=0.
      PME(3,MV+3)=0.
      PME(4,MV+3)=P1(1)
      PME(5,MV+3)=P1(2)
      PME(6,MV+3)=P1(3)
      PME(7,MV+3)=0.
      PME(9,MV+3)=M1
      IME(1,MV+3)=IE1
      IME(2,MV+3)=0
      IME(3,MV+3)=IS1
      IME(4,MV+3)=0
      IME(5,MV+3)=0
      PME(1,MV+1)=0.
      PME(2,MV+1)=0.
      PME(3,MV+1)=0.
      PME(4,MV+1)=-P1(1)
      PME(5,MV+1)=-P1(2)
      PME(6,MV+1)=-P1(3)
      PME(7,MV+1)=0.
      PME(9,MV+1)=M2
      IME(1,MV+1)=IE2
      IME(2,MV+1)=0
      IME(3,MV+1)=IS2
      IME(4,MV+1)=1
      IME(5,MV+1)=0
      ISEX=1
      NP = 2
c       write(*,*)  ' akanni: IR=',IR,' ID1,ID2=',ID1,ID2 
c       write(*,*)  'IE1,IS1,IE2,IS2=',IE1,IS1,IE2,IS2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  ETAN(ID1,ID2,V,U,PIN,IIN,PN,IPN,MV,NP,IETAN)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 V(3),U,PIN(9),PN(9),P1(3),P2(3),MPI,MN,META,M3
      LOGICAL YESELA
      COMMON/YESELA/ YESELA
      COMMON/MEMORY/PME(9,5999),IME(5,5999)
      DIMENSION  IIN(5),IPN(5)
      DATA  MPI/0.140/,MN/0.940/,META/0.549/
      IETAN=0
      IF(ID1.EQ.220.AND.(ID2.EQ.1120.OR.ID2.EQ.1220)) go  to  10
      RETURN
   10 CONTINUE
      KSI=1
      IF(IPN(1).eq.0) KSI=2 
      TETN=(U**2-(PIN(9)+PN(9))**2)/(2.*PN(9))
      stot=SIGETA(KSI,0,TETN)
      sela=SIGETA(KSI,1,TETN)
      RNS=rndm(-1.)
      if (RNS.le.sela/stot)  then
        M3=META              ! elastic scattering
        IE3=0
        IE1=IPN(1)
        IEL=1 
      else
        M3=MPI
        IE1=1 
        if(rndm(-1.).gt.0.5)  IE1=0
        IE3=IPN(1)-IE1
        IEL=0
      endif   
      CT=1.0 - 2.0*rndm(-1.)
      FI=6.283185*rndm(-1.)
      CALL  ABELQ(PIN,V,U,P1,P2,CT,FI,M3,MN)
      PME(1,MV+3)=0.
      PME(2,MV+3)=0.
      PME(3,MV+3)=0.
      PME(4,MV+3)=P1(1)
      PME(5,MV+3)=P1(2)
      PME(6,MV+3)=P1(3)
      PME(7,MV+3)=0.
      PME(9,MV+3)=M3
      IME(1,MV+3)=IE3
      IME(2,MV+3)=0
      IME(3,MV+3)=0
      IME(4,MV+3)=0
      if(IEL.eq.1)  then 
        IME(5,MV+3)=IDINT(1000.0d0 *TAUN(8))
        IF(IME(5,MV+3).LT.1) IME(5,MV+3)=1
      else
        IME(5,MV+3)=0
      endif
      PME(1,MV+1)=0.
      PME(2,MV+1)=0.
      PME(3,MV+1)=0.
      PME(4,MV+1)=-P1(1)
      PME(5,MV+1)=-P1(2)
      PME(6,MV+1)=-P1(3)
      PME(7,MV+1)=0.
      PME(9,MV+1)=MN
      IME(1,MV+1)=IE1
      IME(2,MV+1)=0
      IME(3,MV+1)=0
      IME(4,MV+1)=1
      IME(5,MV+1)=0
      IETAN=2-IEL
      NP = 2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  BBKAK(IDA,IDB,V,U,PIN,PN,IIN,IPN,MV,NP,IBBKAK)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 V(3),U,PIN(9),PN(9),MK,MN
      LOGICAL YESELA
      COMMON/YESELA/ YESELA
      COMMON/MEMORY/PME(9,5999),IME(5,5999)
      COMMON /IDPME/ IDPME(5999)
      COMMON/NCASCA/NCAS,NCPRI
      dimension IIN(5),IPN(5),amas4(4),ps4(5,4)
      DATA  MK/0.494d0/,MN/0.939d0/,zro/0.0d0/,
     &      onth/0.33333333d0/,twth/0.66666667d0/,
     &      onfo/0.25d0/,ontw/0.5d0/,thfo/0.75d0/
      IBBKAK = 0
      Umin=2.0*MN + 2.0*MK
      if(IIN(4).EQ.1.AND.IIN(3).EQ.0.AND.IIN(2).eq.0.and.
     &   IPN(4).EQ.1.AND.IPN(3).EQ.0.AND.IPN(2).eq.0.and.
     &   U.gt.Umin) then
        go  to  1
      else
        return
      endif
   1  continue
      PX1=PIN(4)
      PY1=PIN(5)
      PZ1=PIN(6)
      AM1=PIN(9)
      PX2=PN(4)
      PY2=PN(5)
      PZ2=PN(6)
      AM2=PN(9)
      CALL  CROSEC(1,IDA,IDB,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      SINO=SITO
      IF(YESELA)         GO TO 2
      CALL  CROSEC(0,IDA,IDB,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEL,0)
      CALL  CROSEC(2,IDA,IDB,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEX,0)
      SINO=SITO-SIEL-SIEX
    2 IF(SINO.LE.0.)                     RETURN
      ru = (Umin/U)**2
      skk = 1.5*((1.0d0 - ru)**3.17)*ru**1.96
      if(rndm(-1.).gt.skk/SINO)  return
      if(IIN(5).eq.0.and.IPN(5).eq.0)  then     
c  N + N
        if(IIN(1).eq.1.and.IPN(1).eq.1)  then
c  p + p 
          rnd=rndm(-1.)
          if(rnd.le.onth)  then
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 1    ! p
	    IDPME(MV+2) = 1120        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130        
          elseif(rnd.le.twth)  then
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 1    ! p
	    IDPME(MV+2) = 1120        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
          else
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
          endif
        elseif((IIN(1)+IPN(1)).eq.1)  then
c  p + n 
          rnd=rndm(-1.)
          if(rnd.le.onfo)          then
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130        
          elseif(rnd.le.ontw)      then
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
          elseif(rnd.le.thfo)      then
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 1    ! p
	    IDPME(MV+2) = 1120        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130        
          else
            IME(1,MV+1) = 0    ! n 
	    IDPME(MV+1) = 1220        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
          endif
        else
c  n + n 
          rnd=rndm(-1.)
          if(rnd.le.onth)  then
            IME(1,MV+1) = 0    ! n 
	    IDPME(MV+1) = 1220        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130        
          elseif(rnd.le.twth)  then
            IME(1,MV+1) = 0    ! n 
	    IDPME(MV+1) = 1220        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
          else
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130        
          endif
        endif
      elseif((IIN(5).ne.0.and.IPN(5).eq.0).or.
     &       (IIN(5).eq.0.and.IPN(5).ne.0))  then     
c  D + N
        if((IIN(1).eq.2.and.IPN(1).eq.1).or.
     &     (IIN(1).eq.1.and.IPN(1).eq.2))  then
c D++ + p
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 1    ! p
	    IDPME(MV+2) = 1120        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
        elseif((IIN(1).eq.2.and.IPN(1).eq.0).or.
     &     (IIN(1).eq.0.and.IPN(1).eq.2))  then
c D++ + n
          rnd=rndm(-1.)
          if(rnd.le.onth)  then
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 1    ! p
	    IDPME(MV+2) = 1120        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
          elseif(rnd.le.twth)  then
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 1    ! p
	    IDPME(MV+2) = 1120        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130        
          else
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
          endif
        elseif(IIN(1).eq.1.and.IPN(1).eq.1)  then
c D+ + p
          rnd=rndm(-1.)
          if(rnd.le.onth)  then
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 1    ! p
	    IDPME(MV+2) = 1120        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
          elseif(rnd.le.twth)  then
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 1    ! p
	    IDPME(MV+2) = 1120        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130        
          else
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
          endif
        elseif((IIN(1)+IPN(1)).eq.1)  then
c D+  + n or D0 + p
          rnd=rndm(-1.)
          if(rnd.le.onfo)          then
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130        
          elseif(rnd.le.ontw)      then
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
          elseif(rnd.le.thfo)      then
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 1    ! p
	    IDPME(MV+2) = 1120        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130        
          else
            IME(1,MV+1) = 0    ! n 
	    IDPME(MV+1) = 1220        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
          endif
        elseif((IIN(1)+IPN(1)).eq.0)  then
c D0  + n or D- + p
          rnd=rndm(-1.)
          if(rnd.le.onth)  then
            IME(1,MV+1) = 0    ! n 
	    IDPME(MV+1) = 1220        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130        
          elseif(rnd.le.twth)  then
            IME(1,MV+1) = 0    ! n 
	    IDPME(MV+1) = 1220        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
          else
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130        
          endif
        elseif((IIN(1).eq.-1.and.IPN(1).eq.0).or.
     &     (IIN(1).eq.0.and.IPN(1).eq.-1))  then
c D- + n
            IME(1,MV+1) = 0    ! n 
	    IDPME(MV+1) = 1220        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130      
        else
            write( *,*) ' BBKAK1 D+N : IDA,IDB=',IDA,IDB 
            write(16,*) ' BBKAK1 D+N : IDA,IDB=',IDA,IDB 
            return
        endif
      elseif(IIN(5).ne.0.and.IPN(5).ne.0)  then
c  D + D
        if((IIN(1)+IPN(1)).eq.3)    then
c  D++  + D+
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 1    ! p
	    IDPME(MV+2) = 1120        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
        elseif((IIN(1)+IPN(1)).eq.2)  then
c  D++  + D0  or  D+  +  D+
          rnd=rndm(-1.)
          if(rnd.le.onth)  then
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 1    ! p
	    IDPME(MV+2) = 1120        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130        
          elseif(rnd.le.twth)  then
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 1    ! p
	    IDPME(MV+2) = 1120        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
          else
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
          endif
        elseif((IIN(1)+IPN(1)).eq.1)  then
c  D+  + D0  or  D-  +  D++
          rnd=rndm(-1.)
          if(rnd.le.onfo)          then
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130        
          elseif(rnd.le.ontw)      then
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
          elseif(rnd.le.thfo)      then
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 1    ! p
	    IDPME(MV+2) = 1120        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130        
          else
            IME(1,MV+1) = 0    ! n 
	    IDPME(MV+1) = 1220        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
          endif
        elseif((IIN(1)+IPN(1)).eq.0)  then
c  D-  + D+  or  D0  +  D0
          rnd=rndm(-1.)
          if(rnd.le.onth)  then
            IME(1,MV+1) = 0    ! n 
	    IDPME(MV+1) = 1220        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 1    ! K+
	    IDPME(MV+3) = 130        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130        
          elseif(rnd.le.twth)  then
            IME(1,MV+1) = 0    ! n 
	    IDPME(MV+1) = 1220        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) = 0    ! AK0
	    IDPME(MV+4) =-230        
          else
            IME(1,MV+1) = 1    ! p 
	    IDPME(MV+1) = 1120        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130        
          endif
        elseif((IIN(1)+IPN(1)).eq.-1)  then
c  D-  + D0
            IME(1,MV+1) = 0    ! n 
	    IDPME(MV+1) = 1220        
            IME(1,MV+2) = 0    ! n
	    IDPME(MV+2) = 1220        
            IME(1,MV+3) = 0    ! K0
	    IDPME(MV+3) = 230        
            IME(1,MV+4) =-1    ! K-
	    IDPME(MV+4) =-130        
        else
            write( *,*) ' BBKAK2 D+N : IDA,IDB=',IDA,IDB 
            write(16,*) ' BBKAK2 D+N : IDA,IDB=',IDA,IDB 
            return
        endif
      else
        write( *,*) ' BBKAK3 D+N : IDA,IDB=',IDA,IDB 
        write(16,*) ' BBKAK3 D+N : IDA,IDB=',IDA,IDB 
        return
      endif
      NP=4
      amas4(1) = MN
      amas4(2) = MN
      amas4(3) = MK
      amas4(4) = MK 
      call  genbodl(NP,U,amas4,ps4,w4)
      do  k=1,NP
        PME(4,MV+k) = ps4(1,k)
        PME(5,MV+k) = ps4(2,k)
        PME(6,MV+k) = ps4(3,k)
        PME(7,MV+k) = zro
        PME(9,MV+k) = amas4(k)
        IME(2,MV+k) = 0
        if(k.le.2) IME(3,MV+k) = 0  !  p or n
        if(k.eq.3) IME(3,MV+k) = 1  !  K+ or K0
        if(k.eq.4) IME(3,MV+k) =-1  !  K- or AK0
        if(k.le.2) IME(4,MV+k) = 1  !  p or n
        if(k.gt.2) IME(4,MV+k) = 0  !  K+ or K0,K- or AK0
        IME(5,MV+k) = 0
      enddo
      IBBKAK = 1
      return
      end
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  PIBKAK(IDA,IDB,V,U,PIN,PN,IIN,IPN,MV,NP,IPIBKAK)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 V(3),U,PIN(9),PN(9),MK,MN
      LOGICAL YESELA
      COMMON/YESELA/ YESELA
      COMMON/MEMORY/PME(9,5999),IME(5,5999)
      COMMON /IDPME/ IDPME(5999)
      COMMON/NCASCA/NCAS,NCPRI
      dimension IIN(5),IPN(5),amas3(3),ps3(5,3)
      DATA  MK/0.494d0/,MN/0.939d0/,zro/0.0d0/,
     &      onth/0.33333333d0/,twth/0.66666667d0/,
     &      onfo/0.25d0/,ontw/0.5d0/,thfo/0.75d0/,
     &      onfi/0.20d0/,thfi/0.6d0/,fofi/0.80d0/
      IPIBKAK = 0
      Umin=MN + 2.0*MK
      if((IDA.eq.120.or.IDA.eq.-120.or.IDA.eq.110).and.
     &   (IDB.eq.1120.or.IDB.eq.1220.or.
     &    IDB.eq.1111.or.IDB.eq.1121.or.IDB.eq.2221.or.
     &    IDB.eq.1221).and.U.gt.Umin) then
        go  to  1
      else
        return
      endif
   1  continue
      ru = (Umin/U)**2
      sig0 = 1.121*((1.0d0 - ru)**1.86)*ru**2
      if(IDA.eq.-120.and.IDB.eq.1120)      then
        skk=5.0*sig0                         ! pi-  +  p 
      elseif(IDA.eq.-120.and.IDB.eq.1220)  then
        skk=1.0*sig0                         ! pi-  +  n 
      elseif(IDA.eq. 120.and.IDB.eq.1120)  then
        skk=1.0*sig0                         ! pi+  +  p 
      elseif(IDA.eq. 120.and.IDB.eq.1220)  then
        skk=5.0*sig0                         ! pi+  +  n 
      elseif(IDA.eq. 110.and.IDB.eq.1120)  then
        skk=3.0*sig0                         ! pi0  +  p 
      elseif(IDA.eq. 110.and.IDB.eq.1220)  then
        skk=3.0*sig0                         ! pi0  +  n
      elseif(IDA.eq. 120.and.IDB.eq.1121)  then
        skk=2.0*sig0                         ! pi+  +  D+
      elseif(IDA.eq. 110.and.IDB.eq.1111)  then
        skk=3.0*sig0                         ! pi0  +  D++
      elseif(IDA.eq. 120.and.IDB.eq.2221)  then
        skk=4.0*sig0                         ! pi+  +  D0
      elseif(IDA.eq. 110.and.IDB.eq.1121)  then
        skk=5.0*sig0                         ! pi0  +  D+
      elseif(IDA.eq.-120.and.IDB.eq.1111)  then
        skk=6.0*sig0                         ! pi-  +  D++
      elseif(IDA.eq. 120.and.IDB.eq.1221)  then
        skk=6.0*sig0                         ! pi+  +  D-
      elseif(IDA.eq. 110.and.IDB.eq.2221)  then
        skk=5.0*sig0                         ! pi0  +  D0
      elseif(IDA.eq.-120.and.IDB.eq.1121)  then
        skk=4.0*sig0                         ! pi-  +  D+
      elseif(IDA.eq. 110.and.IDB.eq.1221)  then
        skk=3.0*sig0                         ! pi0  +  D-
      elseif(IDA.eq.-120.and.IDB.eq.2221)  then
        skk=2.0*sig0                         ! pi-  +  D0
      else
        return
      endif 
      PX1=PIN(4)
      PY1=PIN(5)
      PZ1=PIN(6)
      AM1=PIN(9)
      PX2=PN(4)
      PY2=PN(5)
      PZ2=PN(6)
      AM2=PN(9)
      CALL  CROSEC(1,IDA,IDB,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      SINO=SITO
      IF(YESELA)         GO TO 2
      CALL  CROSEC(0,IDA,IDB,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEL,0)
      CALL  CROSEC(2,IDA,IDB,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEX,0)
      SINO=SITO-SIEL-SIEX
    2 IF(SINO.LE.0.)  return       
      if(rndm(-1.).gt.(skk/SINO))  return
      if(IDA.eq.-120.and.IDB.eq.1120)  then
c  pi-   +  p
        rnd=rndm(-1.)
        if(rnd.le.onfi)         then
          IME(1,MV+1) = 1              ! p
          IDPME(MV+1) = 1120 
          IME(1,MV+2) = 0              ! K0
          IDPME(MV+2) = 230 
          IME(1,MV+3) =-1              ! K-
          IDPME(MV+3) =-130 
        elseif(rnd.le.thfi)         then
          IME(1,MV+1) = 0              ! n
          IDPME(MV+1) = 1220 
          IME(1,MV+2) = 1              ! K+
          IDPME(MV+2) = 130 
          IME(1,MV+3) =-1              ! K-
          IDPME(MV+3) =-130 
        else                        
          IME(1,MV+1) = 0              ! n
          IDPME(MV+1) = 1220 
          IME(1,MV+2) = 0              ! K0
          IDPME(MV+2) = 230 
          IME(1,MV+3) = 0              ! AK0
          IDPME(MV+3) =-230 
        endif
      elseif(IDA.eq. 120.and.IDB.eq.1120)  then
c  pi+   +  p
          IME(1,MV+1) = 1              ! p
          IDPME(MV+1) = 1120 
          IME(1,MV+2) = 1              ! K+
          IDPME(MV+2) = 130 
          IME(1,MV+3) = 0              ! AK0
          IDPME(MV+3) =-230 
      elseif(IDA.eq.-120.and.IDB.eq.1220)  then
c  pi-   +  n
          IME(1,MV+1) = 0              ! n
          IDPME(MV+1) = 1220 
          IME(1,MV+2) = 0              ! K0
          IDPME(MV+2) = 230 
          IME(1,MV+3) =-1              ! K-
          IDPME(MV+3) =-130 
      elseif(IDA.eq. 120.and.IDB.eq.1220)  then
c  pi+   +  n
        rnd=rndm(-1.)
        if(rnd.le.onfi)         then
          IME(1,MV+1) = 0              ! n
          IDPME(MV+1) = 1220 
          IME(1,MV+2) = 1              ! K+
          IDPME(MV+2) = 130 
          IME(1,MV+3) = 0              ! AK0
          IDPME(MV+3) =-230 
        elseif(rnd.le.thfi)         then
          IME(1,MV+1) = 1              ! p
          IDPME(MV+1) = 1120 
          IME(1,MV+2) = 1              ! K+
          IDPME(MV+2) = 130 
          IME(1,MV+3) =-1              ! K-
          IDPME(MV+3) =-130 
        else                        
          IME(1,MV+1) = 1              ! p
          IDPME(MV+1) = 1120 
          IME(1,MV+2) = 0              ! K0
          IDPME(MV+2) = 230 
          IME(1,MV+3) = 0              ! AK0
          IDPME(MV+3) =-230 
        endif
      elseif(IDA.eq. 110.and.IDB.eq.1120)  then
c  pi0   +  p
        rnd=rndm(-1.)
        if(rnd.le.twth)         then
          IME(1,MV+1) = 0              ! n
          IDPME(MV+1) = 1220 
          IME(1,MV+2) = 1              ! K+
          IDPME(MV+2) = 130 
          IME(1,MV+3) = 0              ! AK0
          IDPME(MV+3) =-230 
        elseif(rnd.le.5./6.)         then
          IME(1,MV+1) = 1              ! p
          IDPME(MV+1) = 1120 
          IME(1,MV+2) = 1              ! K+
          IDPME(MV+2) = 130 
          IME(1,MV+3) =-1              ! K-
          IDPME(MV+3) =-130 
        else                        
          IME(1,MV+1) = 1              ! p
          IDPME(MV+1) = 1120 
          IME(1,MV+2) = 0              ! K0
          IDPME(MV+2) = 230 
          IME(1,MV+3) = 0              ! AK0
          IDPME(MV+3) =-230 
        endif
      elseif(IDA.eq. 110.and.IDB.eq.1220)  then
c  pi0   +  n
        rnd=rndm(-1.)
        if(rnd.le.twth)         then
          IME(1,MV+1) = 1              ! p
          IDPME(MV+1) = 1120 
          IME(1,MV+2) = 0              ! K0
          IDPME(MV+2) = 230 
          IME(1,MV+3) =-1              ! K-
          IDPME(MV+3) =-130 
        elseif(rnd.le.5./6.)         then
          IME(1,MV+1) = 0              ! n
          IDPME(MV+1) = 1220 
          IME(1,MV+2) = 1              ! K+
          IDPME(MV+2) = 130 
          IME(1,MV+3) =-1              ! K-
          IDPME(MV+3) =-130 
        else                        
          IME(1,MV+1) = 0              ! n
          IDPME(MV+1) = 1220 
          IME(1,MV+2) = 0              ! K0
          IDPME(MV+2) = 230 
          IME(1,MV+3) = 0              ! AK0
          IDPME(MV+3) =-230 
        endif
      elseif(IDA.eq. 120.and.IDB.eq.1121)  then
c  pi+  +  D+
          IME(1,MV+1) = 1              ! p
          IDPME(MV+1) = 1120 
          IME(1,MV+2) = 1              ! K+
          IDPME(MV+2) = 130 
          IME(1,MV+3) = 0              ! AK0
          IDPME(MV+3) =-230 
      elseif(IDA.eq. 110.and.IDB.eq.1111)  then
c  pi0  +  D++
          IME(1,MV+1) = 1              ! p
          IDPME(MV+1) = 1120 
          IME(1,MV+2) = 1              ! K+
          IDPME(MV+2) = 130 
          IME(1,MV+3) = 0              ! AK0
          IDPME(MV+3) =-230 
      elseif(IDA.eq. 120.and.IDB.eq.2221)  then
c  pi+  +  D0
        rnd=rndm(-1.)
        if(rnd.le.onfo)         then
          IME(1,MV+1) = 1              ! p
          IDPME(MV+1) = 1120 
          IME(1,MV+2) = 0              ! K0
          IDPME(MV+2) = 230 
          IME(1,MV+3) = 0              ! AK0
          IDPME(MV+3) =-230 
        elseif(rnd.le.ontw)         then
          IME(1,MV+1) = 1              ! p
          IDPME(MV+1) = 1120 
          IME(1,MV+2) = 1              ! K+
          IDPME(MV+2) = 130 
          IME(1,MV+3) =-1              ! K-
          IDPME(MV+3) =-130 
        else                        
          IME(1,MV+1) = 0              ! n
          IDPME(MV+1) = 1220 
          IME(1,MV+2) = 1              ! K+
          IDPME(MV+2) = 130 
          IME(1,MV+3) = 0              ! AK0
          IDPME(MV+3) =-230 
        endif
      elseif(IDA.eq. 110.and.IDB.eq.1121)  then
c  pi0  +  D+
        rnd=rndm(-1.)
        if(rnd.le.onfi)         then
          IME(1,MV+1) = 1              ! p
          IDPME(MV+1) = 1120 
          IME(1,MV+2) = 0              ! K0
          IDPME(MV+2) = 230 
          IME(1,MV+3) = 0              ! AK0
          IDPME(MV+3) =-230 
        elseif(rnd.le.fofi)         then
          IME(1,MV+1) = 1              ! p
          IDPME(MV+1) = 1120 
          IME(1,MV+2) = 1              ! K+
          IDPME(MV+2) = 130 
          IME(1,MV+3) =-1              ! K-
          IDPME(MV+3) =-130 
        else                        
          IME(1,MV+1) = 0              ! n
          IDPME(MV+1) = 1220 
          IME(1,MV+2) = 1              ! K+
          IDPME(MV+2) = 130 
          IME(1,MV+3) = 0              ! AK0
          IDPME(MV+3) =-230 
        endif
      elseif(IDA.eq.-120.and.IDB.eq.1111)  then
c  pi-  +  D++
        rnd=rndm(-1.)
        if(rnd.le.ontw)         then
          IME(1,MV+1) = 1              ! p
          IDPME(MV+1) = 1120 
          IME(1,MV+2) = 0              ! K0
          IDPME(MV+2) = 230 
          IME(1,MV+3) = 0              ! AK0
          IDPME(MV+3) =-230 
        else                        
          IME(1,MV+1) = 1              ! p
          IDPME(MV+1) = 1120 
          IME(1,MV+2) = 1              ! K+
          IDPME(MV+2) = 130 
          IME(1,MV+3) =-1              ! K-
          IDPME(MV+3) =-130 
        endif
      elseif(IDA.eq. 120.and.IDB.eq.1221)  then
c  pi+  +  D-
        rnd=rndm(-1.)
        if(rnd.le.ontw)         then
          IME(1,MV+1) = 0              ! n
          IDPME(MV+1) = 1220 
          IME(1,MV+2) = 1              ! K+
          IDPME(MV+2) = 130 
          IME(1,MV+3) =-1              ! K-
          IDPME(MV+3) =-130 
        else                        
          IME(1,MV+1) = 0              ! n
          IDPME(MV+1) = 1220 
          IME(1,MV+2) = 0              ! K0
          IDPME(MV+2) = 230 
          IME(1,MV+3) = 0              ! AK0
          IDPME(MV+3) =-230 
        endif
      elseif(IDA.eq. 110.and.IDB.eq.2221)  then
c  pi0  +  D0
        rnd=rndm(-1.)
        if(rnd.le.onfi)         then
          IME(1,MV+1) = 1              ! p
          IDPME(MV+1) = 1120 
          IME(1,MV+2) = 0              ! K0
          IDPME(MV+2) = 230 
          IME(1,MV+3) =-1              ! K-
          IDPME(MV+3) =-130 
        elseif(rnd.le.fofi)         then
          IME(1,MV+1) = 0              ! n
          IDPME(MV+1) = 1220 
          IME(1,MV+2) = 0              ! K0
          IDPME(MV+2) = 230 
          IME(1,MV+3) = 0              ! AK0
          IDPME(MV+3) =-230 
        else
          IME(1,MV+1) = 0              ! n
          IDPME(MV+1) = 1220 
          IME(1,MV+2) = 1              ! K+
          IDPME(MV+2) = 130 
          IME(1,MV+3) =-1              ! K-
          IDPME(MV+3) =-130 
        endif
      elseif(IDA.eq.-120.and.IDB.eq.1121)  then
c  pi-  +  D+
        rnd=rndm(-1.)
        if(rnd.le.onfo)         then
          IME(1,MV+1) = 0              ! n
          IDPME(MV+1) = 1220 
          IME(1,MV+2) = 1              ! K+
          IDPME(MV+2) = 130 
          IME(1,MV+3) =-1              ! K-
          IDPME(MV+3) =-130 
        elseif(rnd.le.ontw)         then
          IME(1,MV+1) = 0              ! n
          IDPME(MV+1) = 1220 
          IME(1,MV+2) = 0              ! K0
          IDPME(MV+2) = 230 
          IME(1,MV+3) = 0              ! AK0
          IDPME(MV+3) =-230 
        else
          IME(1,MV+1) = 1              ! p
          IDPME(MV+1) = 1120 
          IME(1,MV+2) = 0              ! K0
          IDPME(MV+2) = 230 
          IME(1,MV+3) =-1              ! K-
          IDPME(MV+3) =-130 
        endif
      elseif(IDA.eq. 110.and.IDB.eq.1221)  then
c  pi0  +  D-
          IME(1,MV+1) = 0              ! n
          IDPME(MV+1) = 1220 
          IME(1,MV+2) = 0              ! K0
          IDPME(MV+2) = 230 
          IME(1,MV+3) =-1              ! K-
          IDPME(MV+3) =-130 
      elseif(IDA.eq.-120.and.IDB.eq.2221)  then
c  pi-  +  D0
          IME(1,MV+1) = 0              ! n
          IDPME(MV+1) = 1220 
          IME(1,MV+2) = 0              ! K0
          IDPME(MV+2) = 230 
          IME(1,MV+3) =-1              ! K-
          IDPME(MV+3) =-130 
      else
        write( *,*) ' PIBKAK3 pi+B : IDA,IDB=',IDA,IDB 
        return
      endif
      NP=3
      amas3(1) = MN
      amas3(2) = MK
      amas3(3) = MK 
      call  genbodl(NP,U,amas3,ps3,w3)
      do  k=1,NP
        PME(4,MV+k) = ps3(1,k)
        PME(5,MV+k) = ps3(2,k)
        PME(6,MV+k) = ps3(3,k)
        PME(7,MV+k) = zro
        PME(9,MV+k) = amas3(k)
        IME(2,MV+k) = 0
        if(k.eq.1) IME(3,MV+k) = 0  !  p or n
        if(k.eq.2) IME(3,MV+k) = 1  !  K+ or K0
        if(k.eq.3) IME(3,MV+k) =-1  !  K- or AK0
        if(k.eq.1) IME(4,MV+k) = 1  !  p or n
        if(k.gt.1) IME(4,MV+k) = 0  !  K+ or K0,K- or AK0
        IME(5,MV+k) = 0
      enddo
      IPIBKAK = 1
c
      ie0=IIN(1)+IPN(1)
      is0=IIN(3)+IPN(3)
      ib0=IIN(4)+IPN(4)
      ies=0
      iss=0
      ibs=0
      do  k=1,NP
        ies=ies+IME(1,MV+k) 
        iss=iss+IME(3,MV+k) 
        ibs=ibs+IME(4,MV+k)
      enddo
c     if(ies.ne.ie0.or.iss.ne.is0.or.ibs.ne.ib0)  then
c       write( *,100)  ies,iss,ibs,ie0,is0,ib0
c 100   format(6I5)
c     endif          
      return
      end
      SUBROUTINE  PNETA(IDA,IDB,V,U,TIN1,PIN,PN,MV,NP,IETA)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 V(3),U,TIN1,PIN(9),PN(9),MN,META,M1,M2,M3
      LOGICAL YESELA
      COMMON/YESELA/ YESELA
      COMMON/MEMORY/PME(9,5999),IME(5,5999)
      COMMON/NCASCA/NCAS,NCPRI
      common/psigeta/ seta
      DIMENSION  P0(3),P0S(3),P1C(3),P2C(3),P3C(3),
     &           P1S(3),P2S(3),P3S(3),PSUM(3)
      DATA  MN/0.939/,META/0.549/
      IETA=0
      seta = 0.0d0
      if(IDA.eq.1120.and.IDB.eq.1120)  then
        ksi=1
        IC1=1                 ! p
        IC2=1                 ! p 
      elseif(IDA.eq.1220.and.IDB.eq.1220)  then
        ksi=1
        IC1=0                 ! n
        IC2=0                 ! n
      elseif(IDA.eq.1120.and.IDB.eq.1220)  then
        ksi=2
        IC1=1                 ! p
        IC2=0                 ! n
      elseif(IDA.eq.1220.and.IDB.eq.1120)  then
        ksi=2
        IC1=0                 ! n
        IC2=1                 ! p
      else
        return
      endif
      U0=META+MN+MN
      IF(U.LE.U0)                                RETURN
      PX1=PIN(4)
      PY1=PIN(5)
      PZ1=PIN(6)
      AM1=PIN(9)
      PX2=PN(4)
      PY2=PN(5)
      PZ2=PN(6)
      AM2=PN(9)
      CALL  CROSEC(1,IDA,IDB,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      SINO=SITO
      IF(YESELA)         GO TO 1
      CALL  CROSEC(0,IDA,IDB,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEL,0)
      CALL  CROSEC(2,IDA,IDB,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEX,0)
      SINO=SITO-SIEL-SIEX
    1 IF(SINO.LE.0.)                     RETURN
      seta=SPPETA(TIN1,ksi)
      reta=seta/SINO
      if(rndm(-1.).gt.reta)         return
c
c     write(*,*) ' TIN1,SINO,seta,reta=',TIN1,SINO,seta,reta
c
      ID1=IDA
      ID2=IDB
      id3=220
      M1=MN
      M2=MN
      M3=META
      IC3=0                  ! eta     
C       SAMPLING ACCORDING TO 3-BODY PHASE SPACE VOLUME
      IF(U.le.(M1+M2+M3))  RETURN
      EM1=(U**2+M1**2-(M2+M3)**2)/(2.*U)
      EM2=(U**2+M2**2-(M1+M3)**2)/(2.*U)
    2 CONTINUE
      E1=M1+(EM1-M1)*rndm(-1.)
      E2=M2+(EM2-M2)*rndm(-1.)
      E3=U-E1-E2
      IF(E3.LE.M3)                                  GO  TO  2
      FNORM=27.*E1*E2*E3/U**3
      IF(rndm(-1.).GT.FNORM)                        GO  TO  2
      P1=DSQRT(E1**2-M1**2)
      P2=DSQRT(E2**2-M2**2)
      P3=DSQRT(E3**2-M3**2)
      IF(((P1+P2-P3)*(P1-P2+P3)*(P2+P3-P1)).LE.0.)  GO  TO  2
      CT3=1.-2.*rndm(-1.)
      FI3=6.283185*rndm(-1.)
      P0(1)=PX1
      P0(2)=PY1
      P0(3)=PZ1
      CALL  KINEMQ(P0,V,P0S,CT0,ST0,CF0,SF0,T0,AM1)
      ST3=DSQRT(1.-CT3**2)
      P3C(1)=P3*ST3*DCOS(FI3)
      P3C(2)=P3*ST3*DSIN(FI3)
      P3C(3)=P3*CT3
      CALL ROTORQ(P0S,V,P3C,P3S)
      CT1=-(P3**2+P1**2-P2**2)/(2.*P3*P1)
      CT2=-(P3**2+P2**2-P1**2)/(2.*P3*P2)
      ST1=DSQRT(1.-CT1**2)
      ST2=DSQRT(1.-CT2**2)
      FI1=6.283185*rndm(-1.)
      FI2=3.141592+FI1
      P1C(1)=P1*ST1*DCOS(FI1)
      P1C(2)=P1*ST1*DSIN(FI1)
      P1C(3)=P1*CT1
      P2C(1)=P2*ST2*DCOS(FI2)
      P2C(2)=P2*ST2*DSIN(FI2)
      P2C(3)=P2*CT2
      CALL ROTORQ(P3S,V,P1C,P1S)
      CALL ROTORQ(P3S,V,P2C,P2S)
      PME(1,MV+1)=0.
      PME(2,MV+1)=0.
      PME(3,MV+1)=0.
      PME(4,MV+1)=P1S(1)
      PME(5,MV+1)=P1S(2)
      PME(6,MV+1)=P1S(3)
      PME(7,MV+1)=0.
      PME(8,MV+1)=P1
      PME(9,MV+1)=M1
      PME(1,MV+2)=0.
      PME(2,MV+2)=0.
      PME(3,MV+2)=0.
      PME(4,MV+2)=P3S(1)
      PME(5,MV+2)=P3S(2)
      PME(6,MV+2)=P3S(3)
      PME(7,MV+2)=0.
      PME(8,MV+2)=P3
      PME(9,MV+2)=M3
      PME(1,MV+3)=0.
      PME(2,MV+3)=0.
      PME(3,MV+3)=0.
      PME(4,MV+3)=P2S(1)
      PME(5,MV+3)=P2S(2)
      PME(6,MV+3)=P2S(3)
      PME(7,MV+3)=0.
      PME(8,MV+3)=P2
      PME(9,MV+3)=M2
      IME(1,MV+1)=IC1
      IME(2,MV+1)=0
      IME(3,MV+1)=0
      IME(4,MV+1)=1
      IME(5,MV+1)=0
      IME(1,MV+2)=IC3
      IME(2,MV+2)=0
      IME(3,MV+2)=0
      IME(4,MV+2)=0
      IME(5,MV+2)=IDINT(1000.0d0 *TAUN(8))
      IF(IME(5,MV+2).LT.1) IME(5,MV+2)=1
      IME(1,MV+3)=IC2
      IME(2,MV+3)=0
      IME(3,MV+3)=0
      IME(4,MV+3)=1
      IME(5,MV+3)=0
      IETA=1
      NP = 3
      DO  K=1,3
        PSUM(K)=P1S(K)+P2S(K)+P3S(K)
      ENDDO
c     IF(DABS(PSUM(1)).GT.1.D-3.OR.DABS(PSUM(2)).GT.1.D-3.OR.
c    &   DABS(PSUM(3)).GT.1.D-3) THEN
c	   WRITE(16,*) 'PNETA: PSUM=',PSUM
c	   WRITE( *,*) 'PNETA: PSUM=',PSUM
c     ENDIF
      if(NCAS.ge.NCPRI)  then
        pm1=DSQRT(P1S(1)**2+P1S(2)+P1S(3)**2)
        pm2=DSQRT(P2S(1)**2+P2S(2)+P2S(3)**2)
        pm3=DSQRT(P3S(1)**2+P3S(2)+P3S(3)**2)
        en1=DSQRT(pm1**2+M1**2)
        en2=DSQRT(pm2**2+M2**2)
        en3=DSQRT(pm3**2+M3**2)
        write(16,100) U,V
  100   format(28x,'PNETA: U,V=',1PE11.4,3(1PE11.4))
        write(16,101) ID1,P1C,P1S,E1,en1,P1,pm1,M1,IC1
        write(16,101) ID2,P2C,P2S,E2,en2,P2,pm2,M2,IC2
        write(16,101) ID3,P3C,P3S,E3,en3,P3,pm3,M3,IC3
  101   format(1x,I5,11(1PE11.4),I3)
      endif
c
c     write(*,*)  ' PNETA: TIN1,IC1,IC2,IC3',TIN1,IC1,IC2,IC3
c
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  PIYSEX(ID1,ID2,V,U,PIN,PN,MV,NP,ISEX)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     calculation of channel pi + Y ==> AKA + N 
c     Calls: CROSEC, ABELQ
c
      REAL*8 MPI,MN,MK,ML,MS
c
      COMMON/MEMORY/PME(9,5999),IME(5,5999)
      DIMENSION V(3),PIN(9),PN(9),P1(3),P2(3)
      DATA  MPI/0.140/,MN/0.940/,MK/0.494/,ML/1.116/,MS/1.189/
      ISEX=0
	IR=0
      IF((ID1.EQ.110.AND.ID2.EQ.2130).OR.(ID1.EQ.2130.AND.ID2.EQ. 110))
     &      IR=1  ! pi0 + L
      IF((ID1.EQ. 120.AND.ID2.EQ.2130).OR.(ID1.EQ.2130.AND.ID2.EQ. 120))
     &      IR=2  ! pi+ + L
      IF((ID1.EQ.-120.AND.ID2.EQ.2130).OR.(ID1.EQ.2130.AND.ID2.EQ.-120))
     &      IR=3  ! pi- + L
      IF((ID1.EQ. 110.AND.ID2.EQ.1230).OR.(ID1.EQ.1230.AND.ID2.EQ. 110))
     &      IR=4  ! pi0 + S0
      IF((ID1.EQ. 120.AND.ID2.EQ.1230).OR.(ID1.EQ.1230.AND.ID2.EQ. 120))
     &      IR=5  ! pi+ + S0
      IF((ID1.EQ.-120.AND.ID2.EQ.1230).OR.(ID1.EQ.1230.AND.ID2.EQ.-120))
     &      IR=6  ! pi- + S0
      IF((ID1.EQ. 110.AND.ID2.EQ.1130).OR.(ID1.EQ.1130.AND.ID2.EQ. 110))
     &      IR=7  ! pi0 + S+
      IF((ID1.EQ.-120.AND.ID2.EQ.1130).OR.(ID1.EQ.1130.AND.ID2.EQ.-120))
     &      IR=8  ! pi- + S+
      IF((ID1.EQ. 110.AND.ID2.EQ.2230).OR.(ID1.EQ.2230.AND.ID2.EQ. 110))
     &      IR=9  ! pi0 + S-
      IF((ID1.EQ. 120.AND.ID2.EQ.2230).OR.(ID1.EQ.2230.AND.ID2.EQ. 120))
     &      IR=10 ! pi+ + S-
	IF(IR.EQ.0)        RETURN
C
	IF(U.LE.(MN+MK))   RETURN
	EK0=(U**2-MN**2-MK**2)/(2.*MN)
	PK0=SQRT(EK0**2-MK**2)
	EKS=(U**2+MK**2-MN**2)/(2.*U)
	PKS=SQRT(EKS**2-MK**2)
	IF(IR.LE.3)  THEN
C   Calculation of K- + p=>L + pi0 cross section
c   param. of G.Q.Li et al NP A625(1997)342
	   IF(PK0.LE.0.6)        then
            skmp=1.205*PK0**(-1.428)
	   ELSEIF(PK0.LE.1.0)    then
	      skmp=3.5*PK0**0.659
	   ELSE
	      skmp=3.5*PK0**(-3.97)
	   ENDIF
	   EPIS=(U**2+MPI**2-ML**2)/(2.*U)
	   PIS=SQRT(EPIS**2-MPI**2)
c   inverse cross section
	   spiy=2.*((PKS/PIS)**2)*skmp
	ELSE
C   Calculation of K- + p=>S0 + pi0 cross section
c   param. of G.Q.Li et al NP A625(1997)342
	   IF(PK0.LE.0.345)      then
            skmp=0.624*PK0**(-1.830)
	   ELSEIF(PK0.LE.0.425)  then
	      skmp=0.0138/((PK0-0.385)**2+0.0017)
	   ELSE
	      skmp=0.7*PK0**(-2.09)
	   ENDIF
	   EPIS=(U**2+MPI**2-MS**2)/(2.*U)
	   PIS=SQRT(EPIS**2-MPI**2)
c   inverse cross section
	   spiy=4./3.*((PKS/PIS)**2)*skmp
	ENDIF
      PX1=PIN(4)
      PY1=PIN(5)
      PZ1=PIN(6)
      AM1=PIN(9)
      PX2=PN(4)
      PY2=PN(5)
      PZ2=PN(6)
      AM2=PN(9)
      CALL  CROSEC(1,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      IF(RNDM(-1.).GT.(spiy/(SITO+spiy)))  RETURN
	GO  TO  (1,2,3,4,5,6,7,8,9,10),IR
    1 if(RNDM(-1.).le.0.5)  then     ! pi0 + L  ==> K-   + p
        IEK=-1
	  IEN= 1
	else                           ! pi0 + L  ==> AK0  + n
        IEK= 0
	  IEN= 0
	endif
	go  to  11
    2 continue                       ! pi+ + L  ==> AK0  + p
        IEK= 0
	  IEN= 1
	go  to  11
    3 continue                       ! pi- + L  ==> K-   + n
        IEK=-1
	  IEN= 0
	go  to  11
    4 if(RNDM(-1.).le.0.5)  then     ! pi0 + S0 ==> K-   + p
        IEK=-1
	  IEN= 1
      else                           ! pi0 + S0 ==>AK0   + n
        IEK= 0
	  IEN= 0
	endif
	go  to  11
    5 continue                       ! pi+ + S0 ==>AK0   + p
        IEK= 0
	  IEN= 1
	go  to  11
    6 continue                       ! pi- + S0 ==> K-   + n
        IEK=-1
	  IEN= 0
	go  to  11
    7 continue                       ! pi0 + S+ ==>AK0   + p
        IEK= 0
	  IEN= 1
	go  to  11
    8 if(RNDM(-1.).le.0.5)  then     ! pi- + S+ ==> K-   + p
        IEK=-1
	  IEN= 1
      else                           ! pi- + S+ ==>AK0   + n
        IEK= 0
	  IEN= 0
	endif
	go  to  11
    9 continue                       ! pi0 + S- ==> K-   + n
        IEK=-1
	  IEN= 0
	go  to  11
   10 if(RNDM(-1.).le.0.5)  then     ! pi+ + S- ==> K-   + p
        IEK=-1
	  IEN= 1
      else                           ! pi+ + S- ==>AK0   + n
        IEK= 0
	  IEN= 0
	endif
   11 continue
      CT=2.*RNDM(-1.)-1.
      FI=6.283185*RNDM(-1.)
      CALL  ABELQ(PIN,V,U,P1,P2,CT,FI,MK,MN)
      PME(1,MV+3)=0.
      PME(2,MV+3)=0.
      PME(3,MV+3)=0.
      PME(4,MV+3)=P1(1)
      PME(5,MV+3)=P1(2)
      PME(6,MV+3)=P1(3)
      PME(7,MV+3)=0.
      PME(9,MV+3)=MK
      IME(1,MV+3)=IEK
      IME(2,MV+3)=0
      IME(3,MV+3)=-1
      IME(4,MV+3)=0
      IME(5,MV+3)=0
      PME(1,MV+1)=0.
      PME(2,MV+1)=0.
      PME(3,MV+1)=0.
      PME(4,MV+1)=-P1(1)
      PME(5,MV+1)=-P1(2)
      PME(6,MV+1)=-P1(3)
      PME(7,MV+1)=0.
      PME(9,MV+1)=MN
      IME(1,MV+1)=IEN
      IME(2,MV+1)=0
      IME(3,MV+1)=0
      IME(4,MV+1)=1
      IME(5,MV+1)=0
      ISEX=1
      NP = 2
c	write(*,*)  'piysex: IR=',IR
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  PINSEX(ID1,ID2,V,U,PIN,IIN,PN,IPN,MV,NP,ISEX)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     calculation of strange production by pion on nucleon
c     Calls: CROSEC, ABELQ
c
      REAL*8 MPI,MN
      LOGICAL YESELA
      COMMON/YESELA/ YESELA
      COMMON/STREXC/STREXC
      COMMON/MEMORY/PME(9,5999),IME(5,5999)
      DIMENSION V(3),PIN(9),PN(9),P1(3),P2(3)
      DIMENSION  IIN(5),IPN(5)
      DATA  MPI/0.140/,MN/0.940/
      ISEX=0
      STREXC=0.
      IE1=IIN(1)
      IS1=IIN(3)
      IB1=IIN(4)
      IE2=IPN(1)
      IS2=IPN(3)
      IB2=IPN(4)
      IF((IB1+IB2).ne.1.or.IS1.ne.0.or.IS2.ne.0) RETURN
      IF((IE1+IE2).gt.2.or.(IE1+IE2).lt.-1)      RETURN
*      IF(PIN(9).gt.1.d0.or.PN(9).gt.1.d0)        RETURN
c      IF(IB1.EQ.0.AND.IS1.EQ.0.AND.(ID2.EQ.1120.OR.ID2.EQ.1220))
c     *   GO  TO  10
c      RETURN
c   10 CONTINUE
c      IF(ID1.EQ.120.OR.ID1.EQ.-120.OR.ID1.EQ.110) THEN
c        TPIN=(U**2-(PIN(9)+PN(9))**2)/(2.*PN(9))
c      ELSE
        TPIN=(U**2-(MPI+MN)**2)/2./MN
c      ENDIF
      IF(TPIN.LE.0.759)           RETURN
      SSP=SST(0,1,IE1,IE2,4,TPIN,0)
      if((IE1+IE2).eq.2.and.(IE1.eq.2.or.IE2.eq.2))
     &SSP=SST(0,1,  1,  1,4,TPIN,0)
      SSM=SST(0,1,IE1,IE2,5,TPIN,0)
      if((IE1+IE2).eq.-1.and.(IE1.eq.-1.or.IE2.eq.-1))
     &SSM=SST(0,1, -1,  0,5,TPIN,0)
      SS0=SST(0,1,IE1,IE2,6,TPIN,0)
      SL =SST(0,1,IE1,IE2,3,TPIN,0)
      SSEX=SSP+SSM+SS0+SL
      IF(SSEX.LE.1.D-10)            RETURN
      STREXC=SSEX
      PX1=PIN(4)
      PY1=PIN(5)
      PZ1=PIN(6)
      AM1=PIN(9)
      PX2=PN(4)
      PY2=PN(5)
      PZ2=PN(6)
      AM2=PN(9)
      CALL  CROSEC(1,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      SINO=SITO
      IF(YESELA)         GO TO 1
      CALL  CROSEC(0,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEL,0)
      CALL  CROSEC(2,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEX,0)
      SINO=SITO-SIEL-SIEX
    1 IF(SINO.LE.0.)                RETURN
C
c      WRITE(16,100) TPIN,SL,SS0,SSM,SSP,SSEX,SINO
c  100 FORMAT(2X,'TPIN,SL,SS0,SSM,SSP,SSEX,SINO='/1X,7(1PE10.3,1X))
C
      RNS=RNDM(-1.)
c      WRITE(16,'(5x,''U,ID1,ID2,RNS='',5x,F8.3,1x,2I6,1x,F8.3)')
c     &  U,ID1,ID2,RNS
      IF(RNS.GT.(SSEX/SINO))  RETURN
      RN1=RNDM(-1.)
      IF(RN1.LE.(SL/SSEX))    THEN
C     K L
	AMY=1.116
	IEY= 0
      ELSE IF(RN1.LE.((SL+SS0)/SSEX)) THEN
C     K S0
	AMY=1.189
	IEY= 0
      ELSE IF(RN1.LE.((SL+SS0+SSM)/SSEX)) THEN
C     K S-
	AMY=1.189
	IEY=-1
      ELSE
C     K S+
	AMY=1.189
	IEY=1
      ENDIF
      AMK=0.494
      IF(U.LE.(AMK+AMY+0.001))       RETURN
      if(PN(9).gt.1..or.PIN(9).gt.1.)   AMK=0.896
      IF(U.LE.(AMK+AMY+0.001))       AMK=0.494
      IEK=IE1+IE2-IEY
      IF(IEK.LT.0)  THEN
        WRITE(16,*) '  PINSEX: K- ???'
        WRITE( *,*) '  PINSEX: K- ???'
      ENDIF

c      WRITE(16,'(1x,''U,AMK,AMY,ID1,ID2='',3(F8.3,1x),2I6))')
c    & U,AMK,AMY,ID1,ID2
      TMAX=3.5344*TPIN*(TPIN+MPI)/(1.88*(TPIN+MPI)+MPI*MPI+0.8836)
      BET=BHN(TPIN,MPI)
      EXM=BET*TMAX
      IF(EXM.GT.30.)    EXM=30.
      CT=1.+(2.*LOG(1.+RNDM(-1.)*(EXP(-EXM)-1.)))/EXM
      FI=6.283185*RNDM(-1.)
      CALL  ABELQ(PIN,V,U,P1,P2,CT,FI,AMK,AMY)
      PME(1,MV+3)=0.
      PME(2,MV+3)=0.
      PME(3,MV+3)=0.
      PME(4,MV+3)=P1(1)
      PME(5,MV+3)=P1(2)
      PME(6,MV+3)=P1(3)
      PME(7,MV+3)=0.
      PME(9,MV+3)=AMK
      IME(1,MV+3)=IEK
      IME(2,MV+3)=0
      IME(3,MV+3)=1
      IME(4,MV+3)=0
      IME(5,MV+3)=0
      if(AMK.gt.0.50)  IME(5,MV+3)=INTG(1000. *TAUN(12))
*      if(AMK.gt.0.50)  IME(5,MV+3)=INTG(1000. *TAUN(14)) !for K*0
      PME(1,MV+1)=0.
      PME(2,MV+1)=0.
      PME(3,MV+1)=0.
      PME(4,MV+1)=-P1(1)
      PME(5,MV+1)=-P1(2)
      PME(6,MV+1)=-P1(3)
      PME(7,MV+1)=0.
      PME(9,MV+1)=AMY
      IME(1,MV+1)=IEY
      IME(2,MV+1)=0
      IME(3,MV+1)=-1
      IME(4,MV+1)=1
      IME(5,MV+1)=0
      ISEX=1
      NP = 2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  PINETA(ID1,ID2,V,U,PIN,IIN,PN,IPN,MV,NP,IETA)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     calculation of channel pi + N ==> eta + N
c     Calls: CROSEC, ABELQ
c
      REAL*8 MPI,MN,META
      LOGICAL YESELA
c
      COMMON/YESELA/ YESELA
      COMMON/SPINEC/SPINEC
      COMMON/MEMORY/PME(9,5999),IME(5,5999)
      DIMENSION V(3),PIN(9),PN(9),P1(3),P2(3)
      DIMENSION  IIN(5),IPN(5)
      DATA  MPI/0.140/,MN/0.940/,META/0.549/
      IETA=0
      SPINEC=0.
      IF((ID1.EQ.120.OR.ID1.EQ.-120.OR.ID1.EQ.110).AND.
     &   (ID2.EQ.1120.OR.ID2.EQ.1220)) go  to  10
      RETURN
   10 CONTINUE
      IE1=IIN(1)
      IE2=IPN(1)
      IF((IE1+IE2).EQ.2.OR.(IE1+IE2).EQ.-1)     RETURN
      TPIN=(U**2-(PIN(9)+PN(9))**2)/(2.*PN(9))
      SPINEC=SPINET(TPIN)
      IF(SPINEC.LE.1.D-10)                      RETURN
      PX1=PIN(4)
      PY1=PIN(5)
      PZ1=PIN(6)
      AM1=PIN(9)
      PX2=PN(4)
      PY2=PN(5)
      PZ2=PN(6)
      AM2=PN(9)
      CALL  CROSEC(1,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      SINO=SITO
      IF(YESELA)         GO TO 1
      CALL  CROSEC(0,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEL,0)
      CALL  CROSEC(2,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEX,0)
      SINO=SITO-SIEL-SIEX
    1 IF(SINO.LE.0.)                RETURN
C
c      WRITE(16,100) TPIN,SPINEC,SINO
c  100 FORMAT(2X,'TPIN,SPINEC,SINO='/1X,3(1PE10.3,1X))
C
      RNS=RNDM(-1.)
      IF(RNS.GT.(SPINEC/SINO))      RETURN
      TMAX=3.5344*TPIN*(TPIN+MPI)/(1.88*(TPIN+MPI)+MPI*MPI+0.8836)
      BET=BHN(TPIN,MPI)
      EXM=BET*TMAX
      IF(EXM.GT.30.)    EXM=30.
      CT=1.+(2.*LOG(1.+RNDM(-1.)*(EXP(-EXM)-1.)))/EXM
      FI=6.283185*RNDM(-1.)
      CALL  ABELQ(PIN,V,U,P1,P2,CT,FI,META,MN)
      PME(1,MV+3)=0.
      PME(2,MV+3)=0.
      PME(3,MV+3)=0.
      PME(4,MV+3)=P1(1)
      PME(5,MV+3)=P1(2)
      PME(6,MV+3)=P1(3)
      PME(7,MV+3)=0.
      PME(9,MV+3)=META
      IME(1,MV+3)=0
      IME(2,MV+3)=0
      IME(3,MV+3)=0
      IME(4,MV+3)=0
      IME(5,MV+3)=INTG(1000. *TAUN(8))
      PME(1,MV+1)=0.
      PME(2,MV+1)=0.
      PME(3,MV+1)=0.
      PME(4,MV+1)=-P1(1)
      PME(5,MV+1)=-P1(2)
      PME(6,MV+1)=-P1(3)
      PME(7,MV+1)=0.
      PME(9,MV+1)=MN
      IME(1,MV+1)=IE1+IE2
      IME(2,MV+1)=0
      IME(3,MV+1)=0
      IME(4,MV+1)=1
      IME(5,MV+1)=0
      IETA=1
      NP = 2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  PIPIKK(ID1,ID2,V,U,PIN,IIN,PN,IPN,MV,NP,IKAK)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c      calculation of channel pi + pi ==> K + AKA
c     Calls: PIPICS, ABELQ
c
      REAL*8 MK
      COMMON/MEMORY/PME(9,5999),IME(5,5999)
      DIMENSION V(3),PIN(9),PN(9),P1(3),P2(3)
      DIMENSION  IIN(5),IPN(5)
      DATA  MK/0.495/
	DATA a1/8.6335481/,a2/0.554395/,a3/2.750998/,a4/0.893608/
      IKAK=0
      IF((ID1.EQ. 120.AND.ID2.EQ.-120).OR.
     &   (ID1.EQ. 120.AND.ID2.EQ. 110).OR.
     &   (ID1.EQ.-120.AND.ID2.EQ. 120).OR.
     &   (ID1.EQ.-120.AND.ID2.EQ. 110).OR.
     &   (ID1.EQ. 110.AND.ID2.EQ. 120).OR.
     &   (ID1.EQ. 110.AND.ID2.EQ. 110).OR.
     &   (ID1.EQ. 110.AND.ID2.EQ.-120))  GO  TO 10
      RETURN
   10 CONTINUE
      IF(U.LE.(2.*MK))  RETURN

c  kkg  24.03.06
c	x=(U/(2.*MK))**2
c	SKK=a1*((x-1.)**a2)/(a3**2+x)**a4
        STK=1.6    ! B.Tomasik, E.Kolomeitsev,nucl-th/0512088
        IF((ID1.EQ. 120.AND.ID2.EQ.-120).or.
     &     (ID1.EQ.-120.AND.ID2.EQ. 120))   SKK=STK
        IF( ID1.EQ. 110.AND.ID2.EQ. 110)    SKK=2.0/5.0*STK
        IF((ID1.EQ. 120.AND.ID2.EQ. 110).or.
     &     (ID1.EQ.-120.AND.ID2.EQ. 110))   SKK=6.0/5.0*STK

	STOT=PIPICS(ID1,ID2,U,1)
	IF(RNDM(-1.).GT.(SKK/(SKK+STOT)))  RETURN
c	 write(16,*) 'PIPIKK:SKK,STOT=',SKK,STOT
      CT=1.-2.*RNDM(-1.)
      FI=6.283185*RNDM(-1.)
      CALL  ABELQ(PIN,V,U,P1,P2,CT,FI,MK,MK)
      IE1=IIN(1)
      IE2=IPN(1)
	IF((IE1+IE2).EQ.0)         THEN
	  IF(RNDM(-1.).LT.0.5)  THEN
           IME(1,MV+1)= 1
           IME(1,MV+2)=-1
	  ELSE
           IME(1,MV+1)= 0
           IME(1,MV+2)= 0
	  ENDIF
	ELSEIF((IE1+IE2).EQ.1)  THEN
           IME(1,MV+1)= 1
           IME(1,MV+2)= 0
	ELSEIF((IE1+IE2).EQ.-1) THEN
           IME(1,MV+1)= 0
           IME(1,MV+2)=-1
	ELSE
	     WRITE(16,*) 'PIPIKK:ID1,ID2,IE1,IE2=',ID1,ID2,IE1,IE2
	     RETURN
	ENDIF
      PME(1,MV+1)=0.
      PME(2,MV+1)=0.
      PME(3,MV+1)=0.
      PME(4,MV+1)=P1(1)
      PME(5,MV+1)=P1(2)
      PME(6,MV+1)=P1(3)
      PME(7,MV+1)=0.
      PME(9,MV+1)=MK
      IME(2,MV+1)=0
      IME(3,MV+1)=1
      IME(4,MV+1)=0
      IME(5,MV+1)=0
      PME(1,MV+2)=0.
      PME(2,MV+2)=0.
      PME(3,MV+2)=0.
      PME(4,MV+2)=-P1(1)
      PME(5,MV+2)=-P1(2)
      PME(6,MV+2)=-P1(3)
      PME(7,MV+2)=0.
      PME(9,MV+2)=MK
      IME(2,MV+2)=0
      IME(3,MV+2)=-1
      IME(4,MV+2)=0
      IME(5,MV+2)=0
      IKAK=1
      NP = 2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  BBSEX(IDA,IDB,V,U,PIN,PN,MV,NP,ISEX)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     calculation of channels B + B ==> Y + B + K,
c     where B=N, Delta, Y=Lambda,SIGMA, K=K+,K0
c     Calls: CROSEC, BBBYK, ROTORQ
c
      REAL*8 MN,MK,ML,MS,MD,M1,M2,M3
c
      LOGICAL YESELA
      COMMON/YESELA/ YESELA
      COMMON/MEMORY/PME(9,5999),IME(5,5999)
      COMMON/NCASCA/NCAS,NCPRI
      DIMENSION V(3),PIN(9),PN(9)
      DIMENSION  P0(3),P0S(3),P1C(3),P2C(3),P3C(3),
     &           P1S(3),P2S(3),P3S(3),PSUM(3)
      DATA  MN/0.939/,MK/0.494/,ML/1.115/,MS/1.192/,MD/1.236/
      ISEX=0
	U0=MK+ML+MN
      IF(U.LE.U0)                                RETURN
      PX1=PIN(4)
      PY1=PIN(5)
      PZ1=PIN(6)
      AM1=PIN(9)
      PX2=PN(4)
      PY2=PN(5)
      PZ2=PN(6)
      AM2=PN(9)
      CALL  CROSEC(1,IDA,IDB,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      SINO=SITO
      IF(YESELA)         GO TO 1
      CALL  CROSEC(0,IDA,IDB,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEL,0)
      CALL  CROSEC(2,IDA,IDB,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEX,0)
      SINO=SITO-SIEL-SIEX
    1 IF(SINO.LE.0.)                     RETURN
	CALL BBBYK(IDA,IDB,ID1,ID2,ID3,U,SIGSEX)
	IF(SIGSEX.LT.1.D-10)               RETURN
	IF(RNDM(-1.).GT.SIGSEX/SINO)       RETURN
	IF(ID1.EQ.1120.OR.ID1.EQ.1220)   THEN
	  M1=MN
	  IF(ID1.EQ.1120)   IC1= 1
	  IF(ID1.EQ.1220)   IC1= 0
	ELSEIF(ID1.EQ.1111.OR.ID1.EQ.1121.OR.ID1.EQ.2221.OR.ID1.EQ.1221)
     &                                 THEN
	  M1=MD
	  IF(ID1.EQ.1111)   IC1= 2
	  IF(ID1.EQ.1121)   IC1= 1
	  IF(ID1.EQ.2221)   IC1=-1
	  IF(ID1.EQ.1221)   IC1= 0
        CALL  WMD(MD,TAU,FMD)
	ELSE
	  WRITE(16,*) 'BBSEX: ID1=',ID1
	  WRITE( *,*) 'BBSEX: ID1=',ID1
	ENDIF
	IF(ID2.EQ.2130)  THEN
	  M2=ML
	  IC2=0
	ELSEIF(ID2.EQ.1130.OR.ID2.EQ.2230.OR.ID2.EQ.1230)  THEN
	  M2=MS
	  IF(ID2.EQ.1130)   IC2= 1
	  IF(ID2.EQ.2230)   IC2=-1
	  IF(ID2.EQ.1230)   IC2= 0
	ELSE
	  WRITE(16,*) 'BBSEX: ID2=',ID2
	  WRITE( *,*) 'BBSEX: ID2=',ID2
	ENDIF
	IF(ID3.EQ.130.OR.ID3.EQ.230)  THEN
	  M3=MK
	  IF(ID3.EQ.230)   IC3= 0
	  IF(ID3.EQ.130)   IC3= 1
	ELSE
	  WRITE(16,*) 'BBSEX: ID3=',ID3
	  WRITE( *,*) 'BBSEX: ID3=',ID3
	ENDIF
C       SAMPLING ACCORDING TO 3-BODY PHASE SPACE VOLUME
      IF(U.le.(M1+M2+M3))  RETURN
      EM1=(U**2+M1**2-(M2+M3)**2)/(2.*U)
      EM2=(U**2+M2**2-(M1+M3)**2)/(2.*U)
    2 CONTINUE
      E1=M1+(EM1-M1)*RNDM(-1.)
      E2=M2+(EM2-M2)*RNDM(-1.)
      E3=U-E1-E2
      IF(E3.LE.M3)                                  GO  TO  2
      FNORM=27.*E1*E2*E3/U**3
      IF(RNDM(-1.).GT.FNORM)                        GO  TO  2
      P1=SQRT(E1**2-M1**2)
      P2=SQRT(E2**2-M2**2)
      P3=SQRT(E3**2-M3**2)
      IF(((P1+P2-P3)*(P1-P2+P3)*(P2+P3-P1)).LE.0.)  GO  TO  2
      CT3=1.-2.*RNDM(-1.)
	FI3=6.283185*RNDM(-1.)
	P0(1)=PX1
	P0(2)=PY1
	P0(3)=PZ1
      CALL  KINEMQ(P0,V,P0S,CT0,ST0,CF0,SF0,T0,AM1)
	ST3=SQRT(1.-CT3**2)
	P3C(1)=P3*ST3*COS(FI3)
	P3C(2)=P3*ST3*SIN(FI3)
	P3C(3)=P3*CT3
      CALL ROTORQ (P0S,V,P3C,P3S)
	CT1=-(P3**2+P1**2-P2**2)/(2.*P3*P1)
	CT2=-(P3**2+P2**2-P1**2)/(2.*P3*P2)
	ST1=SQRT(1.-CT1**2)
	ST2=SQRT(1.-CT2**2)
	FI1=6.283185*RNDM(-1.)
	FI2=3.141592+FI1
	P1C(1)=P1*ST1*COS(FI1)
	P1C(2)=P1*ST1*SIN(FI1)
	P1C(3)=P1*CT1
	P2C(1)=P2*ST2*COS(FI2)
	P2C(2)=P2*ST2*SIN(FI2)
	P2C(3)=P2*CT2
      CALL ROTORQ (P3S,V,P1C,P1S)
      CALL ROTORQ (P3S,V,P2C,P2S)
      PME(1,MV+1)=0.
      PME(2,MV+1)=0.
      PME(3,MV+1)=0.
      PME(4,MV+1)=P1S(1)
      PME(5,MV+1)=P1S(2)
      PME(6,MV+1)=P1S(3)
      PME(7,MV+1)=0.
      PME(8,MV+1)=P1
      PME(9,MV+1)=M1
      PME(1,MV+2)=0.
      PME(2,MV+2)=0.
      PME(3,MV+2)=0.
      PME(4,MV+2)=P3S(1)
      PME(5,MV+2)=P3S(2)
      PME(6,MV+2)=P3S(3)
      PME(7,MV+2)=0.
      PME(8,MV+2)=P3
      PME(9,MV+2)=M3
      PME(1,MV+3)=0.
      PME(2,MV+3)=0.
      PME(3,MV+3)=0.
      PME(4,MV+3)=P2S(1)
      PME(5,MV+3)=P2S(2)
      PME(6,MV+3)=P2S(3)
      PME(7,MV+3)=0.
      PME(8,MV+3)=P2
      PME(9,MV+3)=M2
	IME(1,MV+1)=IC1
	IME(2,MV+1)=0
	IME(3,MV+1)=0
	IME(4,MV+1)=1
	IF(ID1.EQ.1120.OR.ID1.EQ.1220)   THEN
	  IME(5,MV+1)=0
	ELSE
        IME(5,MV+1)=INTG(1000.*TAU)
	  IF(IME(5,MV+1).LT.1) IME(5,MV+1)=1
	ENDIF
	IME(1,MV+2)=IC3
	IME(2,MV+2)=0
	IME(3,MV+2)=1
	IME(4,MV+2)=0
	IME(5,MV+2)=0
	IME(1,MV+3)=IC2
	IME(2,MV+3)=0
	IME(3,MV+3)=-1
	IME(4,MV+3)=1
	IME(5,MV+3)=0
      ISEX=1
      NP = 3
	DO  K=1,3
	 PSUM(K)=P1S(K)+P2S(K)+P3S(K)
	ENDDO
	IF(ABS(PSUM(1)).GT.1.D-6.OR.ABS(PSUM(2)).GT.1.D-6.OR.
     &   ABS(PSUM(3)).GT.1.D-6) THEN
c	   WRITE(16,*) 'BBSEX: PSUM=',PSUM
c	   WRITE( *,*) 'BBSEX: PSUM=',PSUM
	ENDIF
      if(NCAS.ge.NCPRI)  then
        write(16,100) U,V
  100   format(28x,'BBSEX: U,V=',1PE11.4,3(1PE11.4))
        write(16,101) ID1,P1C,P1S,E1,P1,M1,IC1
        write(16,101) ID2,P2C,P2S,E2,P2,M2,IC2
        write(16,101) ID3,P3C,P3S,E3,P3,M3,IC3
  101   format(1x,I5,9(1PE11.4),I3)
      endif
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE BBBYK(IDA,IDB,ID1,ID2,ID3,SS,SIG)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     calculation probabilities of channels B + B ==> Y + B + K,
c     where B=N, Delta, Y=Lambda,SIGMA, K=K+,K0
c     Calls: DDNYK,DDDYK,DNNYK, DNDYK,NNNYK,NNDYK 
c
      SIG=0.
      IF((IDA.EQ.1111.OR.IDA.EQ.1121.OR.IDA.EQ.2221.OR.IDA.EQ.1221).AND.
     &     (IDB.EQ.1111.OR.IDB.EQ.1121.OR.IDB.EQ.2221.OR.IDB.EQ.1221))
     &     THEN
C     Delta+Delta
	 CALL  DDNYK(IDA,IDB,IDN1,IDN2,IDN3,SS,SIGN)
	 CALL  DDDYK(IDA,IDB,IDD1,IDD2,IDD3,SS,SIGD)

c   kkg  04/12/07
         SIG=SIGN+SIGD
         IF(SIG.LT.1.D-10)     RETURN  
         IF(RNDM(-1.).LE.SIGN/SIG)  THEN
	    ID1=IDN1
	    ID2=IDN2
            ID3=IDN3
         ELSE
	    ID1=IDD1
            ID2=IDD2
	    ID3=IDD3
         ENDIF
      ELSEIF((IDA.EQ.1120.OR.IDA.EQ.1220).AND.
     &(IDB.EQ.1111.OR.IDB.EQ.1121.OR.IDB.EQ.2221.OR.IDB.EQ.1221))  THEN
C  N+Delta
	 CALL  DNNYK(IDB,IDA,IDN1,IDN2,IDN3,SS,SIGN)
	 CALL  DNDYK(IDB,IDA,IDD1,IDD2,IDD3,SS,SIGD)
	 SIG=SIGN+SIGD
         IF(SIG.LT.1.D-10)     RETURN
	 IF(RNDM(-1.).LE.SIGN/SIG)  THEN
	   ID1=IDN1
	   ID2=IDN2
	   ID3=IDN3
	 ELSE
	   ID1=IDD1
	   ID2=IDD2
	   ID3=IDD3
	 ENDIF
      ELSEIF((IDB.EQ.1120.OR.IDB.EQ.1220).AND.
     &(IDA.EQ.1111.OR.IDA.EQ.1121.OR.IDA.EQ.2221.OR.IDA.EQ.1221)) THEN
C  Delta+N
	 CALL  DNNYK(IDA,IDB,IDN1,IDN2,IDN3,SS,SIGN)
	 CALL  DNDYK(IDA,IDB,IDD1,IDD2,IDD3,SS,SIGD)
	 SIG=SIGN+SIGD
         IF(SIG.LT.1.D-10)     RETURN
	 IF(RNDM(-1.).LE.SIGN/SIG)  THEN
	   ID1=IDN1
	   ID2=IDN2
	   ID3=IDN3
	 ELSE
	   ID1=IDD1
	   ID2=IDD2
	   ID3=IDD3
	 ENDIF
      ELSEIF((IDA.EQ.1120.OR.IDA.EQ.1220).AND.
     &       (IDB.EQ.1120.OR.IDB.EQ.1220))         THEN
C  N+N
	 CALL  NNNYK(IDB,IDA,IDN1,IDN2,IDN3,SS,SIGN)
	 CALL  NNDYK(IDB,IDA,IDD1,IDD2,IDD3,SS,SIGD)
	 SIG=SIGN+SIGD
         IF(SIG.LT.1.D-10)     RETURN
	 IF(RNDM(-1.).LE.SIGN/SIG)  THEN
	   ID1=IDN1
	   ID2=IDN2
	   ID3=IDN3
	 ELSE
	   ID1=IDD1
	   ID2=IDD2
	   ID3=IDD3
         ENDIF
      ELSE
         SIG=0.
	 WRITE(16,*) 'BBBYK: IDA,IDB=',IDA,IDB
	 WRITE( *,*) 'BBBYK: IDA,IDB=',IDA,IDB
	 RETURN
      ENDIF
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  DDNYK(IDA,IDB,ID1,ID2,ID3,SS,SIG)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     calculation of channels Delta + Delta ==> Y + N + K,
c     K=K+,K0
c
      REAL*8 MN,MK,MS,ML
      DATA  MN/0.940/,MK/0.494/,ML/1.115/,MS/1.192/
C Parameters for isospin-averaged NN==>NLK cross section
	DATA AL/2.330/,BL/2.140/,CL/5.024/
C Parameters for isospin-averaged NN==>NSK cross section
	DATA AS/15.49/,BS/2.768/,CS/7.222/
	SIG=0.
	if((IDA.eq.1111.and.IDB.eq.1111).or.
     &   (IDA.eq.2221.and.IDB.eq.2221))    RETURN
	SS0L=MN+ML+MK
	IF(SS.LE.SS0L)                  RETURN
	SS0S=MN+MS+MK
	SIG1=1./4.*SIBIR(AL,BL,CL,SS,SS0L)
	SIG2=1./4.*SIBIR(AS,BS,CS,SS,SS0S)
	if((IDA.eq.1111.and.IDB.eq.1121).or.
     &   (IDA.eq.1121.and.IDB.eq.1111))    then
c D++ D+
        SIG=1./4.*SIG2
	  IF(SIG.LT.1.D-10)               RETURN
C p S+ K+
	  ID1=1120           ! p
	  ID2=1130           ! S+
	  ID3=130            ! K+
	elseif((IDA.eq.1111.and.IDB.eq.1221).or.
     &       (IDA.eq.1221.and.IDB.eq.1111).or.
     &       (IDA.eq.1121.and.IDB.eq.1121))    then
c D++ D0 or D+ D+
          SIGL=1./4.*SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
            SIGS=3./4.*SIG2
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)               RETURN
	  IF(RNDM(-1.).LT.SIGL/SIG)  THEN
C p L K+
	    ID1=1120           ! p
	    ID2=2130           ! L
	    ID3=130            ! K+
        ELSE
	    RND=RNDM(-1.)
	    IF(RND.LT.1./3.)  THEN
C p S+ K0
	      ID1=1120           ! p
            ID2=1130           ! S+
	      ID3=230            ! K0
          ELSEIF(RND.LT.2./3.)  THEN
C p S0 K+
	      ID1=1120           ! p
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ELSE
C n S+ K+
	      ID1=1220           ! n
            ID2=1130           ! S+
            ID3=130            ! K+
	    ENDIF
	  ENDIF
	elseif((IDA.eq.1121.and.IDB.eq.1221).or.
     &       (IDA.eq.1221.and.IDB.eq.1121).or.
     &       (IDA.eq.1111.and.IDB.eq.2221).or.
     &       (IDA.eq.2221.and.IDB.eq.1111))    then
c D+ D0 or D++ D-
          SIGL=2./4.*SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
            SIGS=4./4.*SIG2
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)               RETURN
	  IF(RNDM(-1.).LT.SIGL/SIG)  THEN
	    IF(RNDM(-1.).LT.0.5)  THEN
C p L K0
	      ID1=1120           ! p
      	ID2=2130           ! L
	      ID3=230            ! K0
	    ELSE
C n L K+
            ID1=1220           ! n
      	ID2=2130           ! L
            ID3=130            ! K+
	    ENDIF
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LT.0.25)  THEN
C p S0 K0
	      ID1=1120           ! p
            ID2=1230           ! S0
            ID3=230            ! K0
	    ELSEIF(RND.LT.0.50)  THEN
C p S+ K0
	      ID1=1120           ! p
	      ID2=1130           ! S+
	      ID3=230            ! K0
          ELSEIF(RND.LT.0.75)  THEN
C n S0 K+
	      ID1=1220           ! n
            ID2=1230           ! S0
            ID3=130            ! K+
          ELSE
C p S- K+
	      ID1=1220           ! n
	      ID2=2230           ! S-
	      ID3=130            ! K+
          ENDIF
	  ENDIF
	elseif((IDA.eq.1121.and.IDB.eq.2221).or.
     &       (IDA.eq.2221.and.IDB.eq.1121).or.
     &       (IDA.eq.1221.and.IDB.eq.1221))    then
c D+ D- or D0 D0
          SIGL=1./4.*SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
            SIGS=3./4.*SIG2
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)               RETURN
	  IF(RNDM(-1.).LT.SIGL/SIG)  THEN
C n L K0
	    ID1=1220           ! n
          ID2=2130           ! L
	    ID3=230            ! K0
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LT.1./3.)  THEN
C n S0 K0
	      ID1=1220           ! n
	      ID2=1230           ! S0
            ID3=230            ! K0
          ELSEIF(RND.LT.2./3.)  THEN
C n S- K+
	      ID1=1220           ! n
	      ID2=2230           ! S-
	      ID3=130            ! K+
          ELSE
C p S- K0
	      ID1=1120           ! p
	      ID2=2230           ! S-
	      ID3=230            ! K0
	    ENDIF
        ENDIF
	elseif((IDA.eq.1221.and.IDB.eq.2221).or.
     &       (IDA.eq.2221.and.IDB.eq.1221))    then
c D0 D-
        SIG=1./4.*SIG2
	  IF(SIG.LT.1.D-10)               RETURN
C n S- K0
	  ID1=1220           ! n
	  ID2=2230           ! S-
	  ID3=230            ! K0
	else
	  write(16,*) 'DDNYK: IDA,IDB=',IDA,IDB
	  write( *,*) 'DDNYK: IDA,IDB=',IDA,IDB
	endif
	RETURN
	END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  DDDYK(IDA,IDB,ID1,ID2,ID3,SS,SIG)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     calculation of channels Delta + Delta ==> Y + Delta + K,
c     K=K+,K0
c
      REAL*8 MD,MK,ML,MS
      DATA  MD/1.236/,MK/0.494/,ML/1.115/,MS/1.192/
C Parameters for isospin-averaged DD==>DLK cross section
      DATA AL/0.337/,BL/2.149/,CL/7.967/
C Parameters for isospin-averaged DD==>DSK cross section
      DATA AS/5.140/,BS/2.952/,CS/12.05/
	SIG=0.
	SS0L=MD+ML+MK
	IF(SS.LE.SS0L)                  RETURN
	SS0S=MD+MS+MK
	SIG1=1./2.*SIBIR(AL,BL,CL,SS,SS0L)
	SIG2=1./2.*SIBIR(AS,BS,CS,SS,SS0S)
	if(IDA.eq.1111.and.IDB.eq.1111)        then
c D++ D++  ==>  DYK
	  SIG=SIG2
	  IF(SIG.LT.1.D-10)               RETURN
C D++ S+ K+
	  ID1=1111           ! D++
	  ID2=1130           ! S+
	  ID3=130            ! K+
	elseif((IDA.eq.1111.and.IDB.eq.1121).or.
     &       (IDA.eq.1121.and.IDB.eq.1111))  then
c D++ D+  ==>  DYK
	  SIGL=SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=3.*SIG2
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)               RETURN
	  IF(RNDM(-1.).LT.SIGL/SIG)  THEN
C D++ L K+
	    ID1=1111           ! D++
	    ID2=2130           ! L
	    ID3=130            ! K+
        ELSE
	    RND=RNDM(-1.)
	    IF(RND.LT.0.33333)  THEN
C D++ S+ K0
            ID1=1111           ! D++
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ELSEIF(RND.LT.0.66667)  THEN
C D++ S0 K+
	      ID1=1111           ! D++
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ELSE
C D+ S+ K+
	      ID1=1121           ! D+
	      ID2=1130           ! S+
	      ID3=130            ! K+
	    ENDIF
	  ENDIF
	elseif((IDA.eq.1111.and.IDB.eq.1221).or.
     &       (IDA.eq.1221.and.IDB.eq.1111).or.
     &       (IDA.eq.1121.and.IDB.eq.1121))    then
c D++ D+ or D+ D+ ==>  DYK
	  SIGL=2.*SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=5.*SIG2
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)               RETURN
	  IF(RNDM(-1.).LT.SIGL/SIG)  THEN
	    IF(RNDM(-1.).LT.0.5)  THEN
C D++ L K0
	      ID1=1111           ! D++
              ID2=2130           ! L    KKG 01/19/09
	      ID3=230            ! K0
          ELSE
C D+ L K+
            ID1=1121             ! D+
            ID2=2130             ! L   KKG 01/19/09
	      ID3=130            ! K+
	    ENDIF
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LT.0.2)  THEN
C D++ S0 K0
	      ID1=1111           ! D++
              ID2=1230           ! S0
	      ID3=230            ! K0
	    ELSEIF(RND.LT.0.4)  THEN
C D++ S- K+
	      ID1=1111           ! D++
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ELSEIF(RND.LT.0.6)  THEN
C D+ S0 K+
	      ID1=1121           ! D+
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ELSEIF(RND.LT.0.8)  THEN
C D+ S+ K0
	      ID1=1121           ! D+
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ELSE
C D0 S+ K+
	      ID1=1221           ! D0
	      ID2=1130           ! S+
	      ID3=130            ! K+
	    ENDIF
	  ENDIF
	elseif((IDA.eq.1111.and.IDB.eq.2221).or.
     &       (IDA.eq.2221.and.IDB.eq.1111).or.
     &       (IDA.eq.1121.and.IDB.eq.1221).or.
     &       (IDA.eq.1221.and.IDB.eq.1121))    then
c D++ D- or D+ D0
	  SIGL=2.*SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
  	    SIGS=6.*SIG2
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)               RETURN
	  IF(RNDM(-1.).LT.SIGL/SIG)  THEN
	    IF(RNDM(-1.).LT.0.5)  THEN
C D+ L K0
	      ID1=1121           ! D+
	      ID2=2130           ! L
	      ID3=230            ! K0
	    ELSE
C D0 L K+
            ID1=1221           ! D0
	      ID2=2130           ! L
	      ID3=130            ! K+
	    ENDIF
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LT.0.16667)  THEN
C D++ S- K0
	      ID1=1111           ! D++
	      ID2=2230           ! S-
            ID3=230            ! K0
	    ELSEIF(RND.LT.0.33333)  THEN
C D+ S0 K0
	      ID1=1121           ! D+
	      ID2=1230           ! S0
            ID3=230            ! K0
	    ELSEIF(RND.LT.0.5)  THEN
C D+ S- K+
	      ID1=1121           ! D+
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ELSEIF(RND.LT.0.66667)  THEN
C D0 S+ K0
	      ID1=1221           ! D0
	      ID2=1130           ! S+
	      ID3=230            ! K0
          ELSEIF(RND.LT.0.83333)  THEN
C D0 S0 K+
	      ID1=1221           ! D0
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ELSE
C D- S+ K+
	      ID1=2221           ! D-
	      ID2=1130           ! S+
	      ID3=130            ! K+
	    ENDIF
        ENDIF
	elseif((IDA.eq.1121.and.IDB.eq.2221).or.
     &       (IDA.eq.2221.and.IDB.eq.1121).or.
     &       (IDA.eq.1221.and.IDB.eq.1221))    then
c D+ D- or D0 D0 ==>  DYK
	  SIGL=2.*SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=5.*SIG2
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)               RETURN
	  IF(RNDM(-1.).LT.SIGL/SIG)  THEN
	    IF(RNDM(-1.).LT.0.5)  THEN
C D0 L K0
	      ID1=1221           ! D0
	      ID2=2130           ! L
	      ID3=230            ! K0
	    ELSE
C D- L K+
            ID1=2221           ! D-
	      ID2=2130           ! L
            ID3=130            ! K+
	    ENDIF
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LT.0.2)  THEN
C D+ S- K0
	      ID1=1121           ! D+
	      ID2=2230           ! S-
	      ID3=230            ! K0
	    ELSEIF(RND.LT.0.4)  THEN
C D0 S0 K0
	      ID1=1221           ! D0
	      ID2=1230           ! S0
            ID3=230            ! K0
	    ELSEIF(RND.LT.0.6)  THEN
C D0 S- K+
	      ID1=1221           ! D0
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ELSEIF(RND.LT.0.8)  THEN
C D- S+ K0
	      ID1=2221           ! D-
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ELSE
C D- S0 K+
	      ID1=2221           ! D-
	      ID2=1230           ! S0
	      ID3=130            ! K+
          ENDIF
        ENDIF
	elseif((IDA.eq.2221.and.IDB.eq.1221).or.
     &       (IDA.eq.1221.and.IDB.eq.2221))  then
c D- D0  ==>  DYK
	  SIGL=1.*SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=3.*SIG2
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)               RETURN
	  IF(RNDM(-1.).LT.SIGL/SIG)  THEN
C D- L K0
	    ID1=2221           ! D-
	    ID2=2130           ! L
          ID3=230            ! K0
        ELSE
	    RND=RNDM(-1.)
	    IF(RND.LT.0.33333)  THEN
C D- S0 K0
	      ID1=2221           ! D-
	      ID2=1230           ! S0
	      ID3=230            ! K0
	    ELSEIF(RND.LT.0.66667)  THEN
C D- S- K+
	      ID1=2221           ! D-
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ELSE
C D0 S- K0
	      ID1=1221           ! D0
	      ID2=2230           ! S-
	      ID3=230            ! K0
          ENDIF
	  ENDIF
	elseif(IDA.eq.2221.and.IDB.eq.2221)        then
c D- D-  ==>  DYK
	  SIG=1.*SIG2
	  IF(SIG.LT.1.D-10)               RETURN
C D- S- K0
	  ID1=2221           ! D-
	  ID2=2230           ! S-
	  ID3=230            ! K0
	else
	  write(16,*) 'DDDYK: IDA,IDB=',IDA,IDB
	  write( *,*) 'DDDYK: IDA,IDB=',IDA,IDB
	endif
	RETURN
	END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  DNDYK(IDA,IDB,ID1,ID2,ID3,SS,SIG)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     calculation of channels Delta + N ==> Y + Delta + K,
c     K=K+,K0
c
      REAL*8 MD,MK,ML,MS
      DATA  MD/1.236/,MK/0.494/,ML/1.115/,MS/1.192/
C Parameters for D++ p ==>D++ L  K+ cross section
	DATA  A1/2.704/, B1/2.303/,C1/5.551/
C Parameters for D++ p ==>D++ S0 K+ cross section
	DATA  A2/10.30/, B2/2.748/,C2/9.321/
C Parameters for D++ n ==>D++ S- K+ cross section
	DATA  A7/10.33/, B7/2.743/,C7/8.915/
C Parameters for D+  p ==>D+  L  K+ cross section
	DATA A15/2.917/,B15/2.350/,C15/6.557/
C Parameters for D+  p ==>D+  S0 K+ cross section
	DATA A17/10.62/,B17/2.759/,C17/10.20/
C Parameters for D+  p ==>D0  S+ K+ cross section
	DATA A18/0.647/,B18/2.830/,C18/3.862/
C Parameters for D0  p ==>D+  S- K+ cross section
	DATA A21/2.128/,B21/2.843/,C21/5.986/
C Parameters for D+  n ==>D+  S- K+ cross section
	DATA A29/10.57/,B29/2.757/,C29/10.11/
C Parameters for D+  n ==>D0  L  K+ cross section
	DATA A30/0.312/,B30/2.110/,C30/2.165/
C Parameters for D+  n ==>D0  S0 K+ cross section
	DATA A31/1.112/,B31/2.846/,C31/5.943/
	SIG=0.
	SS0L=MD+ML+MK
	IF(SS.LE.SS0L)                  RETURN
	SS0S=MD+MS+MK
	SIG1 =SIBIR(A1,B1,C1,SS,SS0L)
	SIG2 =SIBIR(A2,B2,C2,SS,SS0S)
	SIG7 =SIBIR(A7,B7,C7,SS,SS0S)
	SIG11=0.
	SIG15=SIBIR(A15,B15,C15,SS,SS0L)
	SIG17=SIBIR(A17,B17,C17,SS,SS0S)
	SIG18=SIBIR(A18,B18,C18,SS,SS0S)
	SIG21=SIBIR(A21,B21,C21,SS,SS0S)
	SIG29=SIBIR(A29,B29,C29,SS,SS0S)
	SIG30=SIBIR(A30,B30,C30,SS,SS0L)
	SIG31=SIBIR(A31,B31,C31,SS,SS0S)
	if(IDA.eq.1111.AND.IDB.EQ.1120) THEN
C D++ p ==>DYK
	  SIGL=SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=SIG2+SIG7+3./4.*SIG18
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
C  D++ L K+
	    ID1=1111           ! D++
	    ID2=2130           ! L
	    ID3=130            ! K+
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.SIG2/SIGS)  THEN
C  D++ S0 K+
	      ID1=1111           ! D++
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ELSEIF(RND.LE.(SIG2+SIG7)/SIGS)  THEN
C  D++ S+ K0
	      ID1=1111           ! D++
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ELSE
C  D+ S+ K+
	      ID1=1121           ! D+
	      ID2=1130           ! S+
	      ID3=130            ! K+
	    ENDIF
	  ENDIF
	elseif(IDA.eq.1111.AND.IDB.EQ.1220) THEN
C D++ n ==>DYK
	  SIGL=SIG1+3./4.*SIG30
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=SIG2+SIG7+3./4.*SIG21+3./4.*SIG31+SIG11
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
	    RND=RNDM(-1.)
	    IF(RND.LE.SIG1/SIGL)  THEN
C  D++ L K0
	      ID1=1111           ! D++
	      ID2=2130           ! L
	      ID3=230            ! K0
	    ELSE
C  D+ L K+
	      ID1=1121           ! D+
	      ID2=2130           ! L
	      ID3=130            ! K+
	    ENDIF
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.SIG2/SIGS)  THEN
C  D++ S0 K0
	      ID1=1111           ! D++
	      ID2=1230           ! S0
	      ID3=230            ! K0
	    ELSEIF(RND.LE.(SIG2+SIG7)/SIGS)  THEN
C  D++ S- K+
	      ID1=1111           ! D++
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ELSEIF(RND.LE.(SIG2+SIG7+3./4.*SIG21)/SIGS)  THEN
C  D+ S+ K0
	      ID1=1121           ! D+
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ELSEIF(RND.LE.(SIG2+SIG7+3./4.*SIG21+SIG11)/SIGS)  THEN
C  D0 S+ K+
	      ID1=1221           ! D0
	      ID2=1130           ! S+
	      ID3=130            ! K+
	    ELSE
C  D+ S0 K+
	      ID1=1121           ! D+
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ENDIF
	  ENDIF
	elseif(IDA.eq.1121.AND.IDB.EQ.1120) THEN
C D+ p ==>DYK
	  SIGL=3./4.*SIG30+SIG15
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=3./4.*SIG31+3./4.*SIG21+SIG15+SIG21+SIG17+SIG18
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
	    RND=RNDM(-1.)
	    IF(RND.LE.3./4.*SIG30/SIGL)  THEN
C  D++ L K0
	      ID1=1111           ! D++
	      ID2=2130           ! L
	      ID3=230            ! K0
	    ELSE
C  D+ L K+
	      ID1=1121           ! D+
	      ID2=2130           ! L
	      ID3=130            ! K+
	    ENDIF
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.3./4.*SIG31/SIGS)  THEN
C  D++ S0 K0
	      ID1=1111           ! D++
	      ID2=1230           ! S0
	      ID3=230            ! K0
	    ELSEIF(RND.LE.(3./4.*(SIG31+SIG21))/SIGS)  THEN
C  D++ S- K+
	      ID1=1111           ! D++
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ELSEIF(RND.LE.(3./4.*(SIG31+SIG21)+SIG21)/SIGS)  THEN
C  D+ S+ K0
	      ID1=1121           ! D+
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ELSEIF(RND.LE.(3./4.*(SIG31+SIG21)+SIG21+SIG18)/SIGS) THEN
C  D0 S+ K+
	      ID1=1221           ! D0
	      ID2=1130           ! S+
	      ID3=130            ! K+
	    ELSE
C  D+ S0 K+
	      ID1=1121           ! D+
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ENDIF
	  ENDIF
	elseif(IDA.eq.1221.AND.IDB.EQ.1120) THEN
C D0 p ==>DYK
	  SIGL=SIG30+SIG15
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=SIG31+SIG21+SIG17+SIG29+3./4.*SIG18+3./2.*SIG11
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
	    RND=RNDM(-1.)
	    IF(RND.LE.SIG30/SIGL)  THEN
C  D+ L K0
	      ID1=1121           ! D+
	      ID2=2130           ! L
	      ID3=230            ! K0
	    ELSE
C  D0 L K+
	      ID1=1221           ! D0
	      ID2=2130           ! L
	      ID3=130            ! K+
	    ENDIF
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.SIG31/SIGS)  THEN
C  D+ S0 K0
	      ID1=1121           ! D+
	      ID2=1230           ! S0
	      ID3=230            ! K0
	    ELSEIF(RND.LE.(SIG31+SIG21)/SIGS)  THEN
C  D+ S- K+
	      ID1=1121           ! D+
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ELSEIF(RND.LE.(SIG31+SIG21+SIG17)/SIGS)   THEN
C  D0 S0 K+
	      ID1=1221           ! D0
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ELSEIF(RND.LE.(SIG31+SIG21+SIG17+SIG29))  THEN
C  D0 S+ K0
	      ID1=1221           ! D0
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ELSEIF(RND.LE.(SIG31+SIG21+SIG17+SIG29+3./4.*SIG18))  THEN
C  D- S+ K+
	      ID1=2221           ! D-
	      ID2=1130           ! S+
	      ID3=130            ! K+
	    ELSE
C  D++ S- K0
	      ID1=1111           ! D++
	      ID2=2230           ! S-
	      ID3=230            ! K0
	    ENDIF
	  ENDIF
	elseif(IDA.eq.1121.AND.IDB.EQ.1220) THEN
C D+ n ==>DYK
	  SIGL=SIG15+SIG30
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=SIG17+SIG29+SIG31+SIG21+SIG11+3./4.*SIG18
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
	    RND=RNDM(-1.)
	    IF(RND.LE.SIG15/SIGL)  THEN
C  D+ L K0
	      ID1=1121           ! D+
	      ID2=2130           ! L
	      ID3=230            ! K0
	    ELSE
C  D0 L K+
	      ID1=1221           ! D0
	      ID2=2130           ! L
	      ID3=130            ! K+
	    ENDIF
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.SIG17/SIGS)  THEN
C  D+ S0 K0
	      ID1=1121           ! D+
	      ID2=1230           ! S0
	      ID3=230            ! K0
	    ELSEIF(RND.LE.(SIG17+SIG29)/SIGS)  THEN
C  D+ S- K+
	      ID1=1121           ! D+
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ELSEIF(RND.LE.(SIG17+SIG29+SIG31)/SIGS)  THEN
C  D0 S0 K+
	      ID1=1221           ! D0
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ELSEIF(RND.LE.(SIG17+SIG29+SIG31+SIG21)/SIGS)  THEN
C  D0 S+ K0
	      ID1=1221           ! D0
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ELSEIF(RND.LE.(SIG17+SIG29+SIG31+SIG21+SIG11)/SIGS)  THEN
C  D- S+ K+
	      ID1=2221           ! D-
	      ID2=1130           ! S+
	      ID3=130            ! K+
	    ELSE
C  D++ S- K0
	      ID1=1111           ! D++
	      ID2=2230           ! S-
	      ID3=230            ! K0
	    ENDIF
	  ENDIF
	elseif(IDA.eq.2221.AND.IDB.EQ.1120) THEN
C D- P ==>DYK
	  SIGL=3./4.*SIG30+SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=SIG11+3./4.*SIG31+3./4.*SIG21+SIG7+SIG2
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
	    RND=RNDM(-1.)
	    IF(RND.LE.3./4.*SIG30/SIGL)  THEN
C  D0 L K0
	      ID1=1221           ! D0
	      ID2=2130           ! L
	      ID3=230            ! K0
	    ELSE
C  D- L K+
	      ID1=2221           ! D-
	      ID2=2130           ! L
	      ID3=130            ! K+
	    ENDIF
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.SIG11/SIGS)  THEN
C  D+ S- K0
	      ID1=1121           ! D+
	      ID2=2230           ! S-
	      ID3=230            ! K0
	    ELSEIF(RND.LE.(SIG11+3./4.*SIG31)/SIGS)  THEN
C  D0 S0 K0
	      ID1=1221           ! D0
	      ID2=1230           ! S0
	      ID3=230            ! K0
	    ELSEIF(RND.LE.(SIG11+3./4.*SIG31+3./4.*SIG21)/SIGS)  THEN
C  D0 S- K+
	      ID1=1221           ! D0
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ELSEIF(RND.LE.(SIG11+3./4.*SIG31+3./4.*SIG21+SIG7)/SIGS) THEN
C  D- S+ K0
	      ID1=2221           ! D-
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ELSE
C  D- S0 K+
	      ID1=2221           ! D-
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ENDIF
	  ENDIF
	elseif(IDA.eq.1221.AND.IDB.EQ.1220) THEN
C D0 n ==>DYK
	  SIGL=SIG17+3./4.*SIG30
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=SIG11+SIG15+SIG29+3./4.*SIG21+3./4.*SIG31
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
	    RND=RNDM(-1.)
	    IF(RND.LE.SIG17/SIGL)  THEN
C  D0 L K0
	      ID1=1221           ! D0
	      ID2=2130           ! L
	      ID3=230            ! K0
	    ELSE
C  D- L K+
	      ID1=2221           ! D-
	      ID2=2130           ! L
	      ID3=130            ! K+
	    ENDIF
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.SIG11/SIGS)  THEN
C  D+ S- K0
	      ID1=1121           ! D+
	      ID2=2230           ! S-
	      ID3=230            ! K0
	    ELSEIF(RND.LE.(SIG11+SIG15)/SIGS)  THEN
C  D0 S0 K0
	      ID1=1221           ! D0
	      ID2=1230           ! S0
	      ID3=230            ! K0
	    ELSEIF(RND.LE.(SIG11+SIG15+SIG29)/SIGS)  THEN
C  D0 S- K+
	      ID1=1221           ! D0
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ELSEIF(RND.LE.(SIG11+SIG15+SIG29+3./4.*SIG21)/SIGS)  THEN
C  D- S+ K0
	      ID1=2221           ! D-
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ELSE
C  D- S0 K+
	      ID1=2221           ! D-
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ENDIF
	  ENDIF
	elseif(IDA.eq.2221.AND.IDB.EQ.1220) THEN
C D- N ==>DYK
	  SIGL=SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=SIG18+SIG2+SIG7
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
C  D- L K0
	    ID1=2221           ! D-
	    ID2=2130           ! L
	    ID3=230            ! K0
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.SIG18/SIGS)  THEN
C  D0 S- K0
	      ID1=1221           ! D0
	      ID2=2230           ! S-
	      ID3=230            ! K0
	    ELSEIF(RND.LE.(SIG18+SIG2)/SIGS)  THEN
C  D- S0 K0
	      ID1=2221           ! D-
	      ID2=1230           ! S0
	      ID3=230            ! K0
	    ELSE
C  D- S- K+
	      ID1=2221           ! D-
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ENDIF
	  ENDIF
	else
	  write(16,*) 'DNDYK: IDA,IDB=',IDA,IDB
	  write( *,*) 'DNDYK: IDA,IDB=',IDA,IDB
	endif
	RETURN
	END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  DNNYK(IDA,IDB,ID1,ID2,ID3,SS,SIG)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     calculation of channels Delta + N ==> Y + N + K,
c     K=K+,K0
c
      REAL*8 MN,MK,ML,MS
      DATA  MN/0.939/,MK/0.494/,ML/1.115/,MS/1.192/
C Parameters for D++ n ==>p  L  K+ cross section
	DATA  A1/8.337/, B1/2.227/,C1/2.511/
C Parameters for D-  p ==>n  S- K+ cross section
	DATA A22/52.72/, B22/2.799/,C22/6.303/
	SIG=0.
	SS0L=MN+ML+MK
	IF(SS.LE.SS0L)                  RETURN
	SS0S=MN+MS+MK
	SIG1 =SIBIR(A1,B1,C1,SS,SS0L)
	SIG22=SIBIR(A22,B22,C22,SS,SS0S)
	if(IDA.eq.1111.AND.IDB.EQ.1220) THEN
C D++ n ==>NYK
	  SIGL=SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=(0.5+0.5+1.)*SIG22
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
C  p L K+
	    ID1=1120           ! p
	    ID2=2130           ! L
	    ID3=130            ! K+
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.0.25)  THEN
C  n S+ K+
	      ID1=1220           ! n
	      ID2=1130           ! S+
	      ID3=130            ! K+
	    ELSEIF(RND.LE.0.5)  THEN
C  p S0 K+
	      ID1=1120           ! p
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ELSE
C  p S+ K0
	      ID1=1120           ! p
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ENDIF
	  ENDIF
	elseif(IDA.eq.1111.AND.IDB.EQ.1120) THEN
C D++ p ==>NYK
	  SIG=0.5*SIG22
	  IF(SIG.LT.1.D-10)         RETURN
C  p S+ K+
	      ID1=1120           ! p
	      ID2=1130           ! S+
	      ID3=130            ! K+
	elseif(IDA.eq.1121.AND.IDB.EQ.1220) THEN
C D+ n ==>NYK
	  SIGL=2./3.*SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=SIG22
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
	    RND=RNDM(-1.)
	    IF(RND.LE.0.5)  THEN
C  p L K0
	      ID1=1120           ! p
	      ID2=2130           ! L
	      ID3=230            ! K0
	    ELSE
C  n L K+
	      ID1=1220           ! n
	      ID2=2130           ! L
	      ID3=130            ! K+
	    ENDIF
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.1./6.)  THEN
C  p S0 K0
	      ID1=1120           ! p
	      ID2=1230           ! S0
	      ID3=230            ! K0
	    ELSEIF(RND.LE.2./6.)  THEN
C  n S0 K+
	      ID1=1220           ! n
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ELSEIF(RND.LE.2./3.)  THEN
C  n S+ K0
	      ID1=1220           ! n
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ELSE
C  p S- K+
	      ID1=1120           ! p
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ENDIF
	  ENDIF
	elseif(IDA.eq.1121.AND.IDB.EQ.1120) THEN
C D+ p ==>NYK
	  SIGL=1./3.*SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=1./3.*SIG22
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
C  p L K+
	      ID1=1120           ! p
	      ID2=2130           ! L
	      ID3=130            ! K+
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.0.333333)  THEN
C  p S0 K+
	      ID1=1120           ! p
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ELSE
C  p S+ K0
	      ID1=1120           ! p
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ENDIF
	  ENDIF
	elseif(IDA.eq.1221.AND.IDB.EQ.1220) THEN
C D0 n ==>NYK
	  SIGL=1./3.*SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=SIG22
          endif 
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
C  n L K0
	      ID1=1220           ! n
	      ID2=2130           ! L
	      ID3=230            ! K0
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.1./6.)  THEN
C  n S0 K0
	      ID1=1220           ! n
	      ID2=1230           ! S0
	      ID3=230            ! K0
	    ELSEIF(RND.LE.2./3.)  THEN
C  p S- K0
	      ID1=1120           ! p
	      ID2=2230           ! S-
	      ID3=230            ! K0
	    ELSE
C  n S- K+
	      ID1=1220           ! n
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ENDIF
	  ENDIF
	elseif(IDA.eq.2221.AND.IDB.EQ.1120) THEN
C D- p ==>NYK
	  SIGL=SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=2.*SIG22
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
C  n L K0
	      ID1=1220           ! n
	      ID2=2130           ! L
	      ID3=230            ! K0
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.0.25)  THEN
C  n S0 K0
	      ID1=1220           ! n
	      ID2=1230           ! S0
	      ID3=230            ! K0
	    ELSEIF(RND.LE.0.5)  THEN
C  p S- K0
	      ID1=1120           ! p
	      ID2=2230           ! S-
	      ID3=230            ! K0
	    ELSE
C  n S- K+
	      ID1=1220           ! n
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ENDIF
	  ENDIF
	elseif(IDA.eq.2221.AND.IDB.EQ.1220) THEN
C D- n ==>NYK
	  SIG=0.5*SIG22
	  IF(SIG.LT.1.D-10)         RETURN
C  n S- K0
	      ID1=1220           ! n
	      ID2=2230           ! S-
	      ID3=230            ! K0
	elseif(IDA.eq.1221.AND.IDB.EQ.1120) THEN
C D0 p ==>NYK
	  SIGL=2./3.*SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=SIG22
          endif 
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
	    RND=RNDM(-1.)
	    IF(RND.LE.0.5)  THEN
C  n L K+
	      ID1=1220           ! n
	      ID2=2130           ! L
	      ID3=130            ! K+
	    ELSE
C  p L K0
	      ID1=1120           ! p
	      ID2=2130           ! L
	      ID3=230            ! K0
	    ENDIF
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.1./6.)  THEN
C  n S0 K+
	      ID1=1220           ! n
	      ID2=1230           ! S0
	      ID3=130            ! K+4
	    ELSEIF(RND.LE.1./3.)  THEN
C  p S0 K0
	      ID1=1120           ! p
	      ID2=1230           ! S0
	      ID3=230            ! K0
	    ELSEIF(RND.LE.2./3.)  THEN
C  n S+ K0
	      ID1=1220           ! n
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ELSE
C  p S- K+
	      ID1=1120           ! p
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ENDIF
	  ENDIF
	else
	  write(16,*) 'DNNYK: IDA,IDB=',IDA,IDB
	  write( *,*) 'DNNYK: IDA,IDB=',IDA,IDB
	endif
	RETURN
	END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  NNDYK(IDA,IDB,ID1,ID2,ID3,SS,SIG)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     calculation of channels N + N ==> Y + Delta + K,
c     K=K+,K0
c
      REAL*8 MD,MK,ML,MS
      DATA  MD/1.236/,MK/0.494/,ML/1.115/,MS/1.192/
C Parameters for p p ==>D++ L  K0 cross section
      DATA  A1/6.166/, B1/2.842/,C1/1.960/
C Parameters for p p ==>D++ S- K+ cross section
      DATA A10/10.00/, B10/2.874/,C10/2.543/
	SIG=0.
	SS0L=MD+ML+MK
	IF(SS.LE.SS0L)                  RETURN
	SS0S=MD+MS+MK
	SIG1 =SIBIR(A1,B1,C1,SS,SS0L)
	SIG10=SIBIR(A10,B10,C10,SS,SS0S)
	if(IDA.eq.1120.AND.IDB.EQ.1120) THEN
C p  p ==>DYK
	  SIGL=4./3.*SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=(1./3.+1./6.+1./3.+1.+1./2.)*SIG10
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
	    RND=RNDM(-1.)
	    IF(RND.LE.0.75)  THEN
C  D++ L K0
	      ID1=1111           ! D++
	      ID2=2130           ! L
	      ID3=230            ! K0
	    ELSE
C  D+ L K+
	      ID1=1121           ! D+
	      ID2=2130           ! L
	      ID3=130            ! K+
	    ENDIF
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.2./14.)  THEN
C  D+ S+ K0
	      ID1=1121           ! D+
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ELSEIF(RND.LE.3./14.)  THEN
C  D+ S0 K+
	      ID1=1121           ! D+
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ELSEIF(RND.LE.5./14.)  THEN
C  D0 S+ K+
	      ID1=1221           ! D0
	      ID2=1130           ! S+
	      ID3=130            ! K+
	    ELSEIF(RND.LE.11./14.)  THEN
C  D++ S- K+
	      ID1=1111           ! D++
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ELSE
C  D++ S0 K0
	      ID1=1111           ! D++
	      ID2=1230           ! S0
	      ID3=230            ! K0
	    ENDIF
	  ENDIF
	elseif(IDA.eq.1120.AND.IDB.EQ.1220.OR.
     &       IDA.eq.1220.AND.IDB.EQ.1120) THEN
C n p(p n) ==>DYK
	  SIGL=2./3.*SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=(1./3.+1./6.+1./3.+1./6.+1./2.+1./3.)*SIG10
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
	    RND=RNDM(-1.)
	    IF(RND.LE.0.5)  THEN
C  D+ L K0
	      ID1=1121           ! D+
	      ID2=2130           ! L
	      ID3=230            ! K0
	    ELSE
C  D0 L K+
	      ID1=1221           ! D0
	      ID2=2130           ! L
	      ID3=130            ! K+
	    ENDIF
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.2./11.)  THEN
C  D++ S- K0
	      ID1=1111           ! D++
	      ID2=2230           ! S-
	      ID3=230            ! K0
	    ELSEIF(RND.LE.3./11.)  THEN
C  D+ S0 K0
	      ID1=1121           ! D+
	      ID2=1230           ! S0
	      ID3=230            ! K0
	    ELSEIF(RND.LE.5./11.)  THEN
C  D0 S+ K0
	      ID1=1221           ! D0
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ELSEIF(RND.LE.6./11.)  THEN
C  D0 S0 K+
	      ID1=1221           ! D0
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ELSEIF(RND.LE.9./11.)  THEN
C  D- S+ K+
	      ID1=2221           ! D-
	      ID2=1130           ! S+
	      ID3=130            ! K+
	    ELSE
C  D+ S- K+
	      ID1=1121           ! D+
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ENDIF
	  ENDIF
	elseif(IDA.eq.1220.AND.IDB.EQ.1220) THEN
C n  n ==>DYK
	  SIGL=4./3.*SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=(1./6.+1./3.+1./6.+1.+0.5)*SIG10
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
	    RND=RNDM(-1.)
	    IF(RND.LE.0.75)  THEN
C  D- L K+
	      ID1=2221           ! D-
	      ID2=2130           ! L
	      ID3=130            ! K+
	    ELSE
C  D0 L K0
	      ID1=1221           ! D0
	      ID2=2130           ! L
	      ID3=230            ! K0
	    ENDIF
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.1./13.)  THEN
C  D+ S- K0
	      ID1=1121           ! D+
	      ID2=2230           ! S-
	      ID3=230            ! K0
	    ELSEIF(RND.LE.3./13.)  THEN
C  D0 S- K+
	      ID1=1221           ! D0
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ELSEIF(RND.LE.4./13.)  THEN
C  D0 S0 K0
	      ID1=1221           ! D0
	      ID2=1230           ! S0
	      ID3=230            ! K0
	    ELSEIF(RND.LE.10./13.)  THEN
C  D- S+ K0
	      ID1=2221           ! D-
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ELSE
C  D- S0 K+
	      ID1=2221           ! D-
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ENDIF
	  ENDIF
	else
	  write(16,*) 'NNDYK: IDA,IDB=',IDA,IDB
	  write( *,*) 'NNDYK: IDA,IDB=',IDA,IDB
	endif
	RETURN
	END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  NNNYK(IDA,IDB,ID1,ID2,ID3,SS,SIG)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     calculation of channels N + N ==> Y + N + K,
c     K=K+,K0
c
      REAL*8 MN,MK,ML,MS
      DATA  MN/0.939/,MK/0.494/,ML/1.115/,MS/1.192/
C Parameters for p p ==>p L  K+ cross section
      DATA  A1/1.879/, B1/2.176/,C1/5.264/
C Parameters for p n ==>n L  K+ cross section
      DATA  A2/2.812/, B2/2.121/,C2/4.893/
C Parameters for p p ==>n S+ K+ cross section
      DATA  A5/1.466/, B5/2.743/,C5/3.271/
C Parameters for p p ==>p S+ K0 cross section
      DATA  A6/7.079/, B6/2.760/,C6/8.164/
C Parameters for p p ==>p S0 K+ cross section
      DATA  A7/5.321/, B7/2.753/,C7/8.510/
C Parameters for p n ==>p S- K+ cross section
      DATA  A8/11.02/, B8/2.782/,C8/7.674/
C Parameters for p n ==>n S0 K+ cross section
      DATA  A9/6.310/, B9/2.773/,C9/7.820/
	SIG=0.
	SS0L=MN+ML+MK
	IF(SS.LE.SS0L)                  RETURN
	SS0S=MN+MS+MK
	SIG1 =SIBIR(A1,B1,C1,SS,SS0L)
	SIG2 =SIBIR(A2,B2,C2,SS,SS0L)
	SIG5 =SIBIR(A5,B5,C5,SS,SS0S)
	SIG6 =SIBIR(A6,B6,C6,SS,SS0S)
	SIG7 =SIBIR(A7,B7,C7,SS,SS0S)
	SIG8 =SIBIR(A8,B8,C8,SS,SS0S)
	SIG9 =SIBIR(A9,B9,C9,SS,SS0S)
	if(IDA.eq.1120.AND.IDB.EQ.1120) THEN
C p  p ==>NYK
	  SIGL=SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=SIG5+SIG6+SIG7
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
C  p L K+
	    ID1=1120           ! p
	    ID2=2130           ! L
	    ID3=130            ! K+
        ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.SIG5/SIGS)  THEN
C  n S+ K+
	      ID1=1220           ! n
	      ID2=1130           ! S+
	      ID3=130            ! K+
	    ELSEIF(RND.LE.(SIG5+SIG6)/SIGS)  THEN
C  p S+ K0
	      ID1=1120           ! p
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ELSE
C  p S0 K+
	      ID1=1120           ! p
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ENDIF
	  ENDIF
	elseif(IDA.eq.1120.AND.IDB.EQ.1220.OR.
     &       IDA.eq.1220.AND.IDB.EQ.1120) THEN
C n p(p n) ==>NYK
	  SIGL=2.*SIG2
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=SIG8+SIG9+SIG9+SIG8
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
	    RND=RNDM(-1.)
	    IF(RND.LE.0.5)  THEN
C  p L K0
	      ID1=1120           ! p
	      ID2=2130           ! L
	      ID3=230            ! K0
	    ELSE
C  n L K+
	      ID1=1220           ! n
	      ID2=2130           ! L
	      ID3=130            ! K+
	    ENDIF
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.SIG8/SIGS)  THEN
C  n S+ K0
	      ID1=1220           ! n
	      ID2=1130           ! S+
	      ID3=230            ! K0
	    ELSEIF(RND.LE.(SIG8+SIG9)/SIGS)  THEN
C  n S0 K+
	      ID1=1220           ! n
	      ID2=1230           ! S0
	      ID3=130            ! K+
	    ELSEIF(RND.LE.(SIG8+SIG9+SIG9)/SIGS)  THEN
C  p S0 K0
	      ID1=1120           ! p
	      ID2=1230           ! S0
	      ID3=230            ! K0
	    ELSE
C  p S- K+
	      ID1=1120           ! p
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ENDIF
	  ENDIF
	elseif(IDA.eq.1220.AND.IDB.EQ.1220) THEN
C n  n ==>NYK
	  SIGL=SIG1
          if(SS.le.SS0S)  then
            SIGS=0.
          else
	    SIGS=SIG7+SIG6+SIG5
          endif
	  SIG=SIGL+SIGS
	  IF(SIG.LT.1.D-10)         RETURN
	  IF(RNDM(-1.).LE.SIGL/SIG)  THEN
C  n L K0
	    ID1=1220           ! n
	    ID2=2130           ! L
	    ID3=230            ! K0
	  ELSE
	    RND=RNDM(-1.)
	    IF(RND.LE.SIG7/SIGS)  THEN
C  n S0 K0
	      ID1=1220           ! n
	      ID2=1230           ! S0
	      ID3=230            ! K0
	    ELSEIF(RND.LE.(SIG7+SIG6)/SIGS)  THEN
C  n S- K+
	      ID1=1220           ! n
	      ID2=2230           ! S-
	      ID3=130            ! K+
	    ELSE
C  p S- K0
	      ID1=1120           ! p
	      ID2=2230           ! S-
	      ID3=230            ! K0
	    ENDIF
	  ENDIF
	else
	  write(16,*) 'NNNYK: IDA,IDB=',IDA,IDB
	  write( *,*) 'NNNYK: IDA,IDB=',IDA,IDB
	endif
	RETURN
	END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      DOUBLE PRECISION FUNCTION  SIBIR(A,B,C,SS,SS0)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
c     Calculation of Sibirtsev's functional form of 
c     strange proguction cross sections in B+B==>Y+B+K
c
      SIBIR=0.0
      IF(SS.LE.SS0)   RETURN
      S=SS*SS
      S0=SS0*SS0
      SIBIR=A*((S/S0-1.)**B)*(S0/S)**C
      RETURN
      END
