
C  LAQGSM Initial set module,  K. GUDIMA 04/27/07             *
C                                                                *
C  This module is used for standalone LAQGSM calculations     *
C                                                                *
C
	     SUBROUTINE INITGU
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      REAL*8 MMES,MBAR
	    CHARACTER*4 TSYST(9)
      logical first, fisonly
      
      common /targ0/ atar0,ztar0,ener0
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/KAPPA/XAP/CINSID/INSIDE/CVALON/IVALON
      COMMON/CMALI/CMALI
      COMMON/RCOR/ RCOR/RINT/ RINT
      COMMON/PIDABS/PIDABS
      COMMON/TIMGO/TIMGO,MIN,ITMGO
      COMMON/MMATUR/MMES,MBAR
      COMMON /STIN/ STIN,AMIN
      COMMON/INDINT/INDINT
      COMMON/INTCC/INTCC
      COMMON/TLIMIT/TLIMIT
      COMMON /IACT/ IACT/INDDEC/INDDEC/ISOB3/ISOB3/ISOB2/ISOB2
      COMMON/XBMAX/XBMAX,IFIB0
      COMMON/BGLAB/BLAB,GLAB,KOBR
      COMMON/KSYST/KSYST
      COMMON/CMNN/BNN,GNN
      COMMON/NCASCA/NCAS,NCPRI
      COMMON /IW10/ IW10
      COMMON/NCOLLT/ ARES1,ARES2,COLLT(4)
      COMMON/BRAN/ IRAN
	    COMMON/RNDELT/ RN,DELTA
      COMMON/IFSPE/IFSPE,NCAS1
      COMMON/GEOCRS/ SIG1,SIG2
      COMMON/VGSYS/T0A,VLA,GLA,VEV,GEV,VCM,GCM,ISYS
      common /gbrems/ tgmin, tgmax, teqv, sxabs, ibrems 
      common /nopreco/ sigpre,noprec
c kkg 04/17/05
      COMMON/COMIND/ PUD,SIGMA,ALFA,BETA
      COMMON/COMECB/ ECMB
      COMMON /COMENB/ ENBOU
c sgm 12/04/05
      common /limita/jamin(135), jamax(135)
! CMJ 9/22/166
      ! From iniprec routine
      common /adbf/    amf, r0m, ijsp, nhump
      common /blcemn1/ icemalone
      common /gbzro/   gb0(300,6)
      common /menu1/   ipar1, ipar2, ifam, idel, fisonly
      common /pairc/   cevap,cfis
      common /stopr/   istp
		  common /TNEL/    TNEL
      ! From RTAP10 routine
      COMMON/UPAC2/UP2(70000),LU2
      
      
      DIMENSION RMS(10)
! CMJ 9/22/166
      data first /.true./
      
c sgm 12/04/05
c  limiting the valus of A(iz) only betwen jamin and jamax after INC
      data 
     *jamin/1,	 2,	3,	4,	5,	6,	7,	  
     *8,	9,	11,	13,	15,	17,	19,	   
     *22,	23,	25,	27,	29,	29,	32,
     *31,	36,	35,	40,	37,	43,	41,
     *48,	49,	53,	51,	57,	55,	61,
     *59,	64,	63,	68,	67,	72,	70,
     *76,	74,	80,	78,	84,	82,	88,
     *86,	96,	95,	100,	99,	104,	103,
     *109,	107,	114,	112,	118,	115,	125,
     *120,	130,	126,	134,	131,	140,	136,
     *144,	140,	148,	146,	154,	152,	162,
     *158,	166,	163,	172,	170,	178,	178,
     *184,	182,	189,	188,	195,	192,	200,
     *201,	206,	208,	212,	212,	218,	216,
     *224,	222,	229,	228,	235,	233,	241,
     *240,	247,	240,	253,	256,	259,	262,
     *266,	269,	272,	275,	278,	281,	284,
     *287,	290,	294,	297,	300,	303,	306,
     *310,	313,	316,	319,	323,	326,	329,
     *332,	336/	    						
      data 
     *jamax/12, 14,	17,	22,	25,	28,	31,	
     *34,	39,	44,	47,	54,	55,	60,
     *61,	68,	69,	74,	75,	80,	81,
     *82,	85,	88,	91,	100,	103,	112,
     *114,	118,	119,	120,	121,	126,	127,
     *136,	137,	142,	143,	148,	149,	154,
     *157,	162,	163,	168,	169,	172,	175,
     *178,	181,	183,	185,	188,	191,	192,
     *193,	195,	198,	202,	205,	208,	211,
     *214,	218,	221,	224,	227,	230,	234,
     *237,	240,	243,	247,	250,	253,	256,
     *260,	263,	266,	269,	273,	276,	279,
     *282,	286,	289,	292,	295,	299,	302,
     *305,	308,	312,	315,	318,	312,	325,
     *328,	331,	334,	338,	339,	339,	339,
     *339,	339,	339,	339,	339,	339,	339,
     *339,	339,	339,	339,	339,	339,	339,
     *339,	339,	339,	339,	339,	339,	339,
     *339,	339,	339,	339,	339,	339,	339,
     *339,	339/							
c
      DATA RMS/0.85,2.095,1.976,1.671,2.50,2.57,2.45,2.519,2.45,2.42/
      DATA TSYST/' LAB','ORAT','ORY ','EQUA','L VE','LOC.',
     *'CENT','ER M','ASS '/
C
      CALL  INITAM
      CLOSE(18)
C
      IRAN=12345
c     CALL RDMINI
c
      ARES1=0.
      ARES2=0.
      do k=1,4
       COLLT(k)=0.
      enddo
c
c    kkg  04/17/07
      ECMB  = 10.0 
      ENBOU = 2.000
      PUD=0.415
      SIGMA=0.51           ! 05.05.06
c
	    TLIMIT=100.
	    TIMGO =9000.
	    XBMAX =1.0
	    RCOR  =0.0
	    RINT  =0.0
	    RN    =0.2
	    DELTA =1.3
	    PIDABS=0.01
c	    MMES  =0.1
	    MMES  =0.7
	    MBAR  =0.0
	    IFIB0 =1
c  kkg 07/13/04
      if(ANUCL1.gt.ANUCL2)  then
        tmp=ANUCL2
        ANUCL2=ANUCL1
        ANUCL1=tmp         
        tmp=ZNUCL2
        ZNUCL2=ZNUCL1
        ZNUCL1=tmp         
	      KOBR  = 1        ! in this case all results are boosted 
c                          ! into antilab. system 
      else
        KOBR  = 0 
      endif    
	    INDINT=1
	    INTCC =2
	    KSYST =1       !  Lab.  syst
c	    KSYST =2       !  Eq.V. syst.   
	    IACT  =2
	    INDDEC=0
	    ISOB2 =1
	    ISOB3 =1
c
	    NCPRI =999999999
c	    NCPRI =5760
c       read(15,*)  NCPRI
c	 
	    IW10  =0
	    IFSPE =1
	    XAP   =1.0
	    CMALI =1.1
	    INSIDE=0
	    IVALON=0
	    EPS1  =0.007
	    EPS2  =0.007
	    VPI   =0.025
	    C1    =0.545
	    C2    =0.545
	    D1    =0.05
	    D2    =0.05
	    R0N1  =1.07
	    R0N2  =1.07
      IF(ANUCL1.GT.2.1.AND.ANUCL1.LT.10.1)  THEN
        IA1=INT(ANUCL1+0.1)
        R0N1=RMS(IA1)
      ENDIF
      IF(ANUCL2.GT.2.1.AND.ANUCL2.LT.10.1)  THEN
        IA2=INT(ANUCL2+0.1)
        R0N2=RMS(IA2)
      ENDIF
C for D/p ratio
      IF(ANUCL2.LT.2.1)  D2=0.001
C
      if(AMIN.lt.0.0001.and.ibrems.eq.0)  then
        write(16,275) T0
  275   format(/2x,'gamma + A interaction at Eg=',f7.3,' GeV'/)  
      elseif(AMIN.lt.0.0001.and.ibrems.eq.1)  then
        write(16,276) tgmin,tgmax
  276   format(/2x,'Bremsstrahlung gamma(Egmin=',f7.3,', Egmax=',f7.3,
     &') + A'/)
      endif
      BNN=SQRT(T0/(T0+1.88))
      GNN=SQRT(1.+T0/1.88)
      RM1=0.
      RM2=0.
  100 IF(ANUCL1-2.)  102,102,101
  101 CALL  HELPQ(R0N1,ANUCL1,A1,C1,D1,TF01,RM1)
  102 IF(ANUCL2.GT.1.1)  CALL  HELPQ(R0N2,ANUCL2,A2,C2,D2,TF02,RM2)
      if(ANUCL1.eq.2.0)  RM1=2.158    ! kkg 26.11.01
      if(ANUCL2.eq.2.0)  RM2=2.158    ! kkg 26.11.01
      AM=0.940
      IF(ANUCL1.LE.1.1)  AM=AMIN
      DLAM=1./(5.06*SQRT(T0*(T0+2.*AM)))
      Z0=RM1+RM2+DLAM+DELTA
      SIG1=31.459*(Z0**2)
     **XBMAX**2
      SIG2=31.459*((RM1+RM2)**2)
     **XBMAX**2
      P0=SQRT(T0*(T0+2.*AM))
      BLAB=P0/(T0+AM)
      VLA=BLAB
      if(abs(ANUCL1).lt.0.001)  then
        VCM=P0/(T0+ANUCL2*0.940)
      else
        VCM=P0/(T0+(1.+ANUCL2/abs(ANUCL1))*0.94)
      endif  
c  kkg  12/10/04
      if(AM.le.0.0001) then
        GLAB=0.0
        GLA =0.0 
      else
        GLAB=1./SQRT(1.-BLAB**2)
        GLA =GLAB 
      endif
      T0A=T0
      VEV=P0/(T0+1.88)
      GCM=1./SQRT(1.-VCM**2)
      GEV=1./SQRT(1.-VEV**2)
c
      atar0=ANUCL2
      ztar0=ZNUCL2
      ener0=T0

! CMJ - Variables from iniprec
      if (first) then
	      first = .false.
        en = 300.d0
				do k = 1, 300
c   Can be generalized for more than 6 preeq. ejectiles:
          gb0(k,1) = 16.d0/en
         	gb0(k,2) = 243.d0/en**2
         	gb0(k,3) = gb0(k,1)**3
				enddo
        ifam=12
        fisonly = .false.
        if (idel.eq.2) fisonly = .true.
        nhump=1
        ijsp=0
        cevap = 12.d0
        r0m = 1.2d0
        cfis=14.0d0
        amf=0.1d0
        istp = 0
        TNEL = 0
! CMJ - Variables from RTAP10
        LU2 = 1
c  kkg  03/27/07
        icemalone = 0
        call egsfill
      end if

	    RETURN
	     END
  
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  
      SUBROUTINE INITAM
      IMPLICIT  INTEGER (I-N)
      REAL*8 ENBOU
      
      COMMON/NODCAY/NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOKA0
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/PRINTS/IPRINT
      COMMON/NOTRE/ NOTRE
      COMMON/COMENB/ ENBOU
      COMMON/COMFR/  ICMS
      COMMON/YESELA/ YESELA
      COMMON/KEYHH/ KEYHH
      COMMON/KEYPLA/ KEYPLA
      COMMON/COMQSE/ QSEE,QVSEE
      COMMON/CINSID/ INSIDE
      COMMON/CVALON/ IVALON
      COMMON/COMWTI/ WTIME
      COMMON/IORHEI/ IORHEI
      COMMON/ISOB3/ ISOB3
      LOGICAL WTIME
      LOGICAL NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOKA0
      LOGICAL IPRINT
      LOGICAL LPRNT
      LOGICAL NOTRE
      LOGICAL YESELA
      LOGICAL KEYHH
      LOGICAL KEYPLA
      LOGICAL QSEE,QVSEE
      ITCOM=15
      ITLIS=16
      ITEVT=17
      ITDKY=18
C
      NODCAY=.FALSE.
      NOTRE=.FALSE.
      LPRNT=.FALSE.
c      LPRNT=.TRUE.
      IPRINT=.FALSE.
      YESELA=.FALSE.
      KEYHH=.FALSE.
      KEYPLA=.FALSE.
      QSEE=.FALSE.
      QVSEE=.FALSE.
C     WTIME=.TRUE.
      WTIME=.FALSE.
      ENBOU= 4.4
      INSIDE=0
      IVALON=0
      IORHEI=1
      ISOB3=1
C
C       IF LAB. FRAME
      ICMS=0
C       IF  C.M. FRAME
      ICMS=1
      CALL SETCON
      CALL SETDKY(LPRNT)
      RETURN
      END
      
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      
      SUBROUTINE SETCON
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C          THIS SUBROUTINE SETS THE CONSTANTS IN /CONST/.
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/CONST/PI,SQRT2,ALF,GF,UNITS
      PI=4.*ATAN(1.0D0)
      SQRT2=SQRT(2.0D0)
      ALF=1./137.036
      GF=1.16570D-5
      UNITS=1./2.56815
      RETURN
      END
      
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      
      SUBROUTINE SETDKY(LPRINT)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      CHARACTER*4 IQUIT
      CHARACTER*8 IBLANK,LREAD(10),LMODE(6),LRES
C          THIS SUBROUTINE READS IN THE DECAY TABLE FROM TAPE ITDKY
C          AND SETS UP /DKYTAB/.
C          QUARK-BASED IDENT CODE
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/FORCE/NFORCE,IFORCE(20),MFORCE(5,20)
C     LOOK MUST BE DIMENSIONED TO THE MAXIMUM VALUE OF INDEX
      COMMON/DKYTAB/LOOK(400),CBR(600),MODE(5,600)
      LOGICAL NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOKA0
      COMMON/NODCAY/NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOKA0
      DIMENSION IMODE(6)
      LOGICAL LPRINT
      DATA IQUIT/' END'/,IBLANK/'     '/
      IF(LPRINT) WRITE(ITLIS,10)
10    FORMAT(1H1,30(1H*)/2H *,28X,1H*/
     1 2H *,5X,17HCOLLI DECAY TABLE,5X,1H*/
     2 2H *,28X,1H*/1X,30(1H*)//
     3 6X,4HPART,18X,10HDECAY MODE,19X,6HCUM BR,15X,5HIDENT,17X,
     4 11HDECAY IDENT/)
      LOOP=0
      IOLD=0
      DO 100 I=1,400
      LOOK(I)=0
100   CONTINUE
      IF(NODCAY) RETURN
      REWIND ITDKY
200   LOOP=LOOP+1
      IF(LOOP.GT.600) GO TO 9999
220   DO 210 I=1,5
      IMODE(I)=0
      LMODE(I)=IBLANK
210   CONTINUE
      READ(ITDKY,*,END=300) IRES,ITYPE,BR,(IMODE(I),I=1,5)
      IF(IRES.EQ.0) GO TO 300
C     IF(NOPI0.AND.IRES.EQ.110) GO TO 220
C     IF(NOETA.AND.IRES.EQ.220) GO TO 220
      IF(IRES.EQ.IOLD) GO TO 230
      CALL FLAVOR(IRES,IFL1,IFL2,IFL3,JSPIN,INDEX)
      LOOK(INDEX)=LOOP
230   IOLD=IRES
      CBR(LOOP)=BR
      DO 240 I=1,5
      MODE(I,LOOP)=IMODE(I)
      IF(IMODE(I).NE.0) CALL LABEL(LMODE(I),IMODE(I))
240   CONTINUE
      CALL LABEL(LRES,IRES)
      IF(LPRINT) WRITE(ITLIS,20) LRES,(LMODE(K),K=1,5),
     1BR,IRES,(IMODE(K),K=1,5)
20    FORMAT(6X,A5,6X,5(A5,2X),3X,F8.5,15X,I5,4X,5(I5,2X))
      GO TO 200
C          SET FORCED DECAY MODES
300   IF(NFORCE.EQ.0) GO TO 400
      DO 310 I=1,NFORCE
      LOOP=LOOP+1
      IF(LOOP.GT.600) GO TO 9999
      CALL FLAVOR(IFORCE(I),IFL1,IFL2,IFL3,JSPIN,INDEX)
      LOOK(INDEX)=LOOP
      DO 320 K=1,5
320   MODE(K,LOOP)=MFORCE(K,I)
      CBR(LOOP)=1.
310   CONTINUE
C          READ AND PRINT NOTES FROM DECAYTABLE FILE
400   IF(.NOT.LPRINT) RETURN
410   READ(ITDKY,1002,END=9998) LREAD
1002  FORMAT(10A8)
      IF(LREAD(1).EQ.IQUIT) RETURN
      WRITE(ITLIS,1003) LREAD
1003  FORMAT(1X,10A8)
      GO TO 410
9998  RETURN
9999  WRITE(ITLIS,30)
30    FORMAT(//44H ***** ERROR IN SETDKY ... LOOP.GT.600 *****)
      RETURN
      END
      
c **********************************************************************
    
      SUBROUTINE inigam0 (fgam, egamma)

c ======================================================================
c
c    Main routine to extract ds/do for channel 1-22:
c
c    Written by K. K. Gudima, Fall 2003?
c    Modified by AJS, May, 2004.
c    Modified by KKG, Nov., 2004
c ======================================================================

      implicit real*8 (a-h, o-z), integer (i-n)

      character*60 fgam


c ======================================================================
 
      common /ixsgpn/ thetai(181),cthetai(181),si(22,182),ri(22,181)
      common /xsecgn0/ xsectd(22,50,0:18),ecm(22,50),elg(22,50)

      dimension theta(19), ctheta(19), s(22,19), st(22), r(22,19)

      data  dtheta, dtheti /10.d0, 1.0d0/
      data zro, two /0.d0, 2.d0/, pi/3.141592d0/, emn/0.939d0/
      data ifirst/1/ 
c ======================================================================
      twpi =2.0d0*pi
	degrad = pi/180.0d0 
      if (ifirst.ne.1) go  to  9
	ifirst = 0
c ***      Read differential cross section data file
!      open (17, file=fgam, type='old')  ! REP 04 Aug 2009
      open (17, file=fgam, status='old')
      do j = 1,19
        theta(j) = dble(j-1)*dtheta
        ctheta(j) = cos(theta(j)*degrad)
      enddo
      do jch = 1,22
        read (17, 40)
        do inw = 1,50
          read (17, 30) ecm(jch,inw), (xsectd(jch,inw,inth), inth=0,4)
          read (17, 20) (xsectd(jch,inw,inth), inth=5,18)
          elg(jch,inw) = (ecm(jch,inw)**2 - emn**2)/(two*emn)
        end do
      end do
      close (17)
    9 continue
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
          write (*,*) '  stop in INIGAM: ieg1, ieg2 = ', ieg1, ieg2
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
c ======================================================================

c   7 format (20x,5e10.3)
   20 format (7e10.3)
   30 format (7x,f5.2,8x,5e10.3)
   40 format (1x)
c 100 format(/12x,' differential- and r-sections for Egamma = ',f7.4/
c    &              12x,'g+p=>pi+ + n',11x,'g+p=>pi0 + p',11x,
c    &                  'g+n=>pi- + p',7x,'g+n=>pi0 + n'/
c    &   1x'theta',' costhet',4(1x,'ds/do(mkb/sr)',3x,'r',4x))           
c 101 format(1x,f5.1,1x,f8.5,8(1pe11.4))
c 102 format(1x,'sint = ',9x,4(1pe11.4,11x))
      return

c ======================================================================
      END                                           
