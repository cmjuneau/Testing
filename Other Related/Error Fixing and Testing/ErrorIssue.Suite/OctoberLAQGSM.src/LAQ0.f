      SUBROUTINE CASCAW(NEL,RN,DELTA,MV,IRET,kstart)
c
c     Modified: 19-DEC-2003 BY NVM
c     Edited by KKG, May 2007 to include GDR region fot gamma + A 
c     interactoin below 20 MeV
c     Modified: 09-MAY-2007 BY NVM
!     Modified: 08/24/16 by CMJ to allow CEM use

      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      REAL *8 RN,DELTA
! CMJ 08/16
      integer*4 casc, indx
      
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE,
     *VEP(3),VET(3),GEP,GET
      COMMON/ACTIV/MPA(300),MYP(5999),MYT(5999),MYY(5999)
      COMMON/CENTER/XC(2,300),YC(2,300),ZC(2,300),IZ(2,300)
      COMMON/MEMORY/PMEMO(9,5999),IMEMO(5,5999)
     */KYPT/ KPT,KYP,KYT,KPT1,KYP1,KYT1,KYY,KYY1
      COMMON/EXCIT/TEX1,TEX2,HEX1,HEX2,PEX1,PEX2
      COMMON/JGSTAR/ JG
      COMMON/NRAPI/NRAPI
      COMMON/NUCOLL/ NUCOLL,MVCOLL(5999)
      COMMON/ACTIM/TINT
      COMMON/INTCC/INTCC
      COMMON/INTCEN/IPO,INT1(300),INT2(5999),INT3(5999),
     *              INT4(100000),IJPA,IST,JST,NRST
      COMMON/IPAUL/IP
      COMMON/TIMGO/TIMGO,MINUTE,ITMGO
      COMMON/IPRIMC/IPRIMC
      COMMON/IDPME/IDPME(5999)
      COMMON/MVUP/MVU
      COMMON/BIMP/BIMP,BIMX,BIMY
      COMMON/NCOLLT/ ARES1,ARES2,COLLT(4)
      COMMON /STIN/ STIN,AMIN
      COMMON/GEOCRS/ SIG1,SIG2
c  kkg 10/14/03
      COMMON/INTTYP/ ITYP
      COMMON /NCASCA/ NCAS,NCPRI
      COMMON/BGLAB/BLAB,GLAB,KOBR
      common /gbrems/ tgmin, tgmax, teqv, sxabs, ibrems
c CMJ 08/2016
		  common /TNEL/    TNEL
      common /casc/    casc
      common /blok77/  spt(5,150)
      common /sptxyz/  sptr(3,150)
      common /zapp/    parz(6,150)
      common /counter/ icntr, icnt
      common /gener/   ngen(150), ing, igs(150), meso
      common /pi/      pi, twpi
     
      DIMENSION P1(9),P2(9),IP1(5),IP2(5),V12(3),P3(9),IP3(5), PTEMP(3)
      dimension partin(9), ipatin(5)

      data zro, hlf, one, two, ten /0.d0, 0.5d0, 1.d0, 2.d0, 10.d0/
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NEL=0
      IRET=0

 9    MVU=0
      NUCOLL=0
      AN1=ANUCL1	! Projectile Size
      ZN1=ZNUCL1
      AN2=ANUCL2 ! Target Size
      ZN2=ZNUCL2
      ENEXT1=ZRO ! Exit Energies
      ENEXT2=ZRO
      PNUCL1(1)=ZRO ! X-Momenta (Projectile)
      PNUCL1(2)=ZRO ! Y-Momenta
      PNUCL1(3)=ZRO ! Z-Momenta
      PNUCL2(1)=ZRO ! X-Momenta (Target)
      PNUCL2(2)=ZRO ! Y-Momenta
      PNUCL2(3)=ZRO ! Z-Momenta
      AMNUC1(1)=ZRO ! X-Angular Momenta (Projectile)
      AMNUC1(2)=ZRO ! Y-Angular Momenta
      AMNUC1(3)=ZRO ! Z-Angular Momenta
      AMNUC2(1)=ZRO ! X-Angular Momenta (Target)
      AMNUC2(2)=ZRO ! Y-Angular Momenta
      AMNUC2(3)=ZRO ! Z-Angular Momenta
      TINT=ZRO ! Interaction Time
      MV=0
      KPT=0
      KYP=0
      KYT=0
      KYY=0
      KPT1=0
      KYP1=0
      KYT1=0
      KYY1=0
      TEX1=0 ! Exciton Information
      TEX2=ZRO
      PEX1=ZRO
      PEX2=ZRO
      HEX1=ZRO
      HEX2=ZRO
      DO L=1,300
	      MPA(L)=1
      END DO
      ing = 0						! A generation index, assigned to PARZ(5,k)
      kstart = 1

! Determining entry point of particle in nucleus
      CALL PINPNQ(RN,R0X,R0Y,R0Z,MV,DELTA,NEL)
      OBR1=ZRO
      IF(RM1.GT.0.1D0) OBR1=0.00144D0*ZN1/RM1
      OBR2=ZRO
      IF(RM2.GT.0.1D0) OBR2=0.00144D0*ZN2/RM2
! Initial velocities
      CALL  VINIT(VPR,VTA,ANUCL1,ANUCL2,T0)
      RADP(1)=R0X
      RADP(2)=R0Y
      RADP(3)=R0Z
      RADT(1)=ZRO
      RADT(2)=ZRO
      RADT(3)=ZRO

c  extension  to Egamma < 20 MeV
      if(AMIN.lt.0.0001D0.and.T0.lt.0.020D0)  then
        t0mev = T0*1000.D0
        sgabs = gabs(t0mev, ANUCL2)
        pint = sgabs/SIG1
        if(RNDM(-1.).gt.pint)  then
          go  to  21
        else  
          AN2=ANUCL2
          ZN2=ZNUCL2
          ENEXT2 = T0
          PNUCL2(1)=ZRO
          PNUCL2(2)=ZRO
          PNUCL2(3)= T0 
          AMNUC2(1)=1.0D0
          AMNUC2(2)=ZRO
          AMNUC2(3)=ZRO
          HEX2 = 1.0D0
          TEX2 = 1.0D0
          PEX2 = ZRO
          if(RNDM(-1.).le.(ZNUCL2/ANUCL2)) PEX2 = 1.0D0
          if(ibrems.eq.1)  then
            teqv = teqv + T0
            sxabs = sxabs + sgabs
          endif
          return
        endif
      endif
c 
      IP=-1
      IJPA=0
 11   NA1=AN1+0.1D0
      NA2=AN2+0.1D0
c  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(NA1.ge.1.and.NA2.ge.1)  then    !
         do  l=1,NA1                     !
           MPA(l)=1                      !
         enddo                           !
      endif 
c   kkg 11/05/03                             !
      if(MV.ge.1.and.NUCOLL.ne.0)   then ! 05.09.97
         do  l=1,MV                      !
           if(NA1.ge.1)  MYP(l)=1        !
           if(NA2.ge.1)  MYT(l)=1        !
           if(IABS(IDPME(l)).le.20) then !
                         MYP(l)=0        !
                         MYT(l)=0        !
                         MYY(l)=0        !
           endif                         !
         enddo                           !
      endif                              !
c  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IPRIMC=0
CNVM      CALL  STOPGU(TIMGO,ISTOP)
CNVM      IF(ISTOP.NE.0)   IRET=1
CNVM      IF(ISTOP.NE.0)   RETURN
      CALL  UPACOW(MV,MVU,0)
      CALL  DISSIP(IP,T0) ! Changing momenta due to coulomb interaction
c
c	if(NCAS.gt.NCPRI) then
c       write(16,*) ' NCAS,AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2=',
c    &                NCAS,AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2 
c       write( *,*) ' NCAS,AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2=',
c    &                NCAS,AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2 
c       endif 
c
c
c Determining interaction hadron pair and collision type
      CALL POINTN(M,P1,P2,IP1,IP2,V12,U12,T12,SIG,SAB,TINT,NR,N1,N2,
     *NA1,NA2,MV,DELTA,R0X,R0Y,R0Z)
      NTRY=0
      IF(M)   203,203,12
  203 IF(MV.EQ.0)   GO  TO  205
      DO  204  L=1,MV
      MYP(L)=0
      MYT(L)=0
      MYY(L)=0
  204 CONTINUE
  205 CONTINUE
c  kkg 17.06.02
c     IF(MV.eq.0.and.AN2.eq.ANUCL2.and.ZN2.eq.ZNUCL2.and.
c    &   ENEXT2.lt.0.001)  GO  TO  21
      IF(ENEXT2.LT.1.D-7)  GO  TO  21
      CALL  UPACOW(MV,MVU,1)
      MV=MVU
      GO  TO  20
   12 GO  TO  (13,14,13,112,15),NR
   13 NU=2
      GO TO 15
   14 NU=1
   15 CONTINUE
      NTRY=NTRY+1
c  kkg  09/05/08
c Determines interaction type
      CALL TYPNEW(P1,IP1,P2,IP2,V12,U12,T12,SAB,MV,NP,NABS,
     *P3,IP3,N3,NU,N2,NA1,NA2,IRET)
      if(IRET.eq.1)   go  to  11
c
      NRAPI=NR
C
C
      IP=0
c  kkg 10/14/03
      IF(ITYP.eq.0)  go  to  11
c 
      IF(NP)11,11,16
   16 GO TO  (17,18,19,112,113),NR
   17 CONTINUE
C
      CALL PAULI1(P1,P2,IP1,IP2,N1,N2,V12,NP,MV,TINT,IP)
		  !write(*,*) "-> Pauli1"
C
      KPT1=KPT1+IP
      KPT=KPT+1
      GO  TO  11
   18 CONTINUE
C
      CALL PAULI2(P1,P2,P3,IP1,IP2,IP3,N1,N2,N3,V12,NP,MV,TINT,
     *IP,OBR1)
		  !write(*,*) "-> Pauli2"
C
      KYP1=KYP1+IP
      KYP=KYP+1
      GO  TO  11
   19 CONTINUE
      IF(ANUCL1.LT.1.1D0.AND.NUCOLL.LE.1) IPRIMC=1
C
      CALL PAULI3(P1,P2,P3,IP1,IP2,IP3,N1,N2,N3,V12,NP,MV,TINT,
     *IP,OBR2)
		  !write(*,*) "-> Pauli3"
C
      KYT1=KYT1+IP
      KYT=KYT+1
      GO  TO  11
  112 CALL  DECREN(N1,NUR,MV,NP)
! Exciton absorbed, storing data
		  !write(*,*) "-> Decren"
      DO 213 L=1,9
      P1(L)=PMEMO(L,N1)
      IF(L.LE.5) IP1(L)=IMEMO(L,N1)
      if(N1.le.150) then
        partin(L) = pmemo(L,N1)
        if(L.le.5) ipatin(L) = imemo(L,N1)
      endif
  213 CONTINUE
      CALL  PAULID(MV,N1,NUR,1,NP,IP)
		  !write(*,*) "-> PauliD"
      GO  TO  11
  113 CONTINUE
C
      CALL PAULI4(P1,P2,IP1,IP2,N1,N2,V12,NP,MV,TINT,IP)
		 !write(*,*) "-> Pauli4"
C
      KYY1=KYY1+IP
      KYY=KYY+1
C
      GO  TO  11
c  kkg 17.06.02
   20 continue
c  kkg  18/04/03
      if(ANUCL1.ge.0.9D0)  then 
        Exmax2=ANUCL1*(T0+EPS2) ! Exmax2 is the overal energy, GeV
      else 
        Exmax2=T0+AMIN
      endif
	    Exmin=0.0001D0 ! 100 MeV
      if((AN1.LT.ZN1.OR.ZN1.LT.ZRO.or.ENEXT1.lt.ZRO).or.  
     &   (AN2.LT.ZN2.OR.ZN2.LT.ZRO.or.ENEXT2.lt.Exmin).or.
     &   (ENEXT2.gt.Exmax2)) then  
c         write(*,*) 'AN1,ZN1,AN2,ZN2,ENEXT1,ENEXT2=',
c     &               AN1,ZN1,AN2,ZN2,ENEXT1,ENEXT2
         IRET=1
	       RETURN
	     else 
	       go  to  22
	     endif               
c  kkg  17.06.02
   21 NEL=NEL+1
      GO TO 9
   22 DO 23 L=1,3
!   Convert angular momentum to units of h-bar.  Constant =
!   1000./(h-bar*c = 197.32858)
      AMNUC1(L)=AMNUC1(L)*5.06D0
   23 AMNUC2(L)=AMNUC2(L)*5.06D0
c   kkg 12/10/04
      if(AMIN.lt.0.0001D0.and.ibrems.eq.1)  then   ! bremss. gamma
         teqv = teqv + T0
      endif   
      IF(ANUCL1-2.1D0)  24,26,26
   26 AM1=.94D0*AN1 ! new projectile rest mass
! assigning prjectil momenta
      CALL CINEMA(PNUCL1,VPR,PTEMP,CT,ST,CF,SF,TTEMP,AM1)
      PNUCL1(1)=PTEMP(1)
      PNUCL1(2)=PTEMP(2)
      PNUCL1(3)=PTEMP(3)
c  kkg 07/13/04
      if(KOBR.eq.1)  then
        PNUCL1(3)=-GLAB*(PNUCL1(3)-BLAB*(TTEMP+AM1))
      endif
   24 CONTINUE
! Assigning target momenta
      IF(ANUCL2.LT.2.1D0)   GO  TO  27
      AM2=.94D0*AN2
      CALL CINEMA(PNUCL2,VTA,PTEMP,CT,ST,CF,SF,TTEMP,AM2)
      PNUCL2(1)=PTEMP(1)
      PNUCL2(2)=PTEMP(2)
      PNUCL2(3)=PTEMP(3)
c  kkg 07/13/04
      if(KOBR.eq.1)  then
        PNUCL2(3)=-GLAB*(PNUCL2(3)-BLAB*(TTEMP+AM2))
      endif
   27 CONTINUE
      COLLT(1)=COLLT(1)+KPT1
      COLLT(2)=COLLT(2)+KYP1
      COLLT(3)=COLLT(3)+KYT1
      COLLT(4)=COLLT(4)+KYY1
      ARES1=ARES1+AN1
      ARES2=ARES2+AN2
      
! CMJ 08/16 - created for CEM use
      if(casc.eq.1) then
! partin(9) array
        do indx = 1, 9
          partin(indx) = P2(indx)
        enddo
! ipatin(5) array
        do indx = 1, 5
          ipatin(indx) = IP2(indx)
        enddo
! spt(5,150) array
        if (kstart.gt.150) then
          write(*,*) 'Array out of bounds in CASCAW'
! spt is out of bounds if this is the case
!          write (16, 100) icntr
!          return
        else
          temp = dble(ipatin(1)) ! Charge of particle
          spt(1,kstart) = partin(4)
          spt(2,kstart) = partin(5)
          spt(3,kstart) = partin(8)
          spt(4,kstart) = temp
          spt(5,kstart) = partin(9)
        endif
! parz(6,150) array
        if (kstart.gt.150) then
! parz is out of bounds if this is the case
!          write (16, 100) icntr
!          return
        else
          parz(6,kstart) = temp
          parz(5,kstart) = dble(ing) ! Particle Index
          if (partin(9).le.hlf) then
!   Pion:
            if (temp.lt.zro) then
!   pi-
              parz(1,kstart) = 7.d0
            elseif (temp.eq.zro) then
!   pi0
              parz(1,kstart) = 8.d0
            else
!   pi+
              parz(1,kstart) = 9.d0
            endif
          else
!    Nucleon:
            if (temp.le.zro) then
!   Neutron:
              parz(1,kstart) = one
            else
!   Proton:
              parz(1,kstart) = two
            endif
            if (indi.eq.3) parz(5,kstart) = -parz(5,kstart)
          endif
          parz(3,kstart) = atan2(partin(4), partin(5))
          fi1 = atan2(partin(6), partin(7))
          if (fi1.lt.zro) fi1 = twpi + fi1
          parz(4,kstart) = fi1
          parz(2,kstart) = partin(8)
        endif
! sptr array and misc.
        if (kstart.gt.150) then
! sptr is out of bounds if this is the case
!          write (16, 100) icntr
!          return
        else
          igs(kstart) = ing
          sptr(1,kstart) = RM2*partin(1)
          sptr(2,kstart) = RM2*partin(2)
          sptr(3,kstart) = RM2*partin(3)
        endif
        kstart = kstart + 1
      endif

      if(NCAS.gt.NCPRI) then
	 write(*,*) ' CASCAW: NCAS,AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2=',
     &                NCAS,AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2 
      RETURN
      endif
			
      RETURN
      
! ======================================================================

! 100 format  (/5x,'More than 150 particles in spt after CASCAD in ', &
!       'trial number ',i7,'.')

! ======================================================================
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  UPACOW(MV,MVU,IU)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c      REAL*4 UP1
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
c  kkg  11/05/03
      if(IM(2,M).ne.0.and.NUCOLL.eq.0)               go to 20
c  kkg  11/05/03
   10 CONTINUE
      I=MVU+1
      IF((11*I).LT.66000)       GO  TO  201
      IPER=1
      write(16,'(15X,''ARRAY  UP1 IS EXCEEDED IN UPACOW''/)')
      write( *,'(15X,''ARRAY  UP1 IS EXCEEDED IN UPACOW''/)')
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
*      UP1(11*I   )=1.           !!!  weight=1.
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
c
C  *************************************************************
C
      SUBROUTINE  PRSPE1(LIM1,LIM2,NC,NC1,INE,IW11)
      implicit real*8 (a-h, o-z), integer (i-n)
      CHARACTER*4  AAAA(6),AAAA1(8), SYST(3)
      CHARACTER*4 TEX1(3),TEX2(6),TEX3(6),TEXT(22)
      logical prtz(151),printit
c  kkg 04/07/04
      INTEGER*8 INE 
c  kkg 12/14/05
      common /fitaf2/fitpt(2),fitpt1(2)
c
      COMMON/INSP/INSP
      COMMON/SPE/SPE(11,11,205)
      COMMON/MULT/AMUL(11),MULT(11)
      COMMON/DELEN/WF,Nfiss /AZT/AT(2),ZT(2)
      COMMON/VGSYS/T0A,VLA,GLA,VEV,GEV,VCM,GCM,ISYS
      COMMON/COARAD/P0D,P0T,P0A
      COMMON/FISSION/W1,W2,NF1,NF2
      COMMON/INDZAP/INDZAP
      COMMON/TGRLE/ T0le,tgr(4),dtgr(4),jgr(4)
      COMMON/GEOCRS/ SIG1,SIG2
c
	common /pmultip/ pmult(5,9)
      common /targ0/ atar0,ztar0,ener0
      common /proj0/ apro0,zpro0
      common /nopreco/ sigpre,noprec
      common /angwin/ ewin(2,11),angw(25,11),dthw
      common /gbrems/ tgmin, tgmax, teqv, sxabs, ibrems
c
c SGM, 03/16/05
      common /indyeld/ iyeld 
c KKG     12/14/05
      common /indgsi/ IGSI
      common /nevtyp/  nevtype
      dimension ar(351), z(151)
      dimension dasum(151)
      dimension iar(351),izr(151),apm(5)
c kkg

      DIMENSION C(11),CC(11),NBA(11),DOM(10)
      DATA NBA/2*0,2*1,2,3,3,4,0,0,1/
      REAL MA(11)/2*0.140,2*0.940,1.88,2*2.82,3.76,2*0.492,0.940/
      DATA  SYST/'LAB.','E.V.','CMS.'/
      DATA  AAAA/'CASC','AD  ','+PRE','COMP','+FER','MI  '/
      DATA  AAAA1/'CASC','AD  ','+COA','LES ','+PRE','COMP',
     *'+FER','MI  '/
c
C SGM, 04/01/13
c     *TEXT/'PI-M','INUS',' PI-','PLUS','NEUT','RONS',
c     *     ' PRO','TONS','DEUT','RONS',' TRI','TONS ','    ','HE-3',
c     *     '    ','HE-4','    ','K+  ','    ','K-  ','ANTI','PROT'/,
c SGM, 04/01/13
c
      data zro/0.0d0/, one/1.0d0/
c********************************************************************
c     INSP   X          U      FUNCTION         WEIGHT
c
c      1     P         TETA    E*dS/d3P     SQRT(P**2+M**2)/P**2
c      2     T=E-M     TETA    dS/dT/dO     1.
c      3     P         TETA    dN/dP/dO     1.
c
c***************************************************************
c     IGSI
c      1     produced fragments distribution includes only fragments from
c            projectile fragmentation. in this case the kinetic energy, 
c            emission angle, z-velocity are calculated in the projectile 
c            rest frame 
c      2     produced fragments distribution includes only fragments from
c            target     fragmentation. in this case the kinetic energy, 
c            emission angle, z-velocity are calculated in the target 
c            rest frame (lab. system)
c      3     produced fragments distribution includes all fragments from
c            both projectile and target fragmentation. the kinetic energy, 
c            emission angle, z-velocity are calculated in the target 
c            rest frame (lab. system)
c************************************************************************  
c
			LIM1=1
      LIM2=1
      P0D=0.090d0
      P0T=0.108d0
      P0A=0.115d0
      CALL  GITAB
c
      T0le=5.0
      tgr(1)=0.020d0
      tgr(2)=0.100d0
      tgr(3)=1.000d0
      tgr(4)=5.150d0
      dtgr(1)=0.002d0
      dtgr(2)=0.005d0
      dtgr(3)=0.010d0
      dtgr(4)=0.050d0
      jgr(1)=1
      jgr(2)=11
      jgr(3)=27
      jgr(4)=117
c
      ewin(1,1)=0.
      ewin(2,1)=1.0 
      ewin(1,2)=0.
      ewin(2,2)=1.0 
      ewin(1,3)=0.0028
      ewin(2,3)=0.200
      ewin(1,4)=0.0028
      ewin(2,4)=0.200
      ewin(1,5)=0.0036
      ewin(2,5)=0.2 
      ewin(1,6)=0.0041
      ewin(2,6)=0.2 
      ewin(1,7)=0.015
      ewin(2,7)=0.2 
      ewin(1,8)=0.015
      ewin(2,8)=0.2 
      ewin(1,9)=0.026
      ewin(2,9)=0.2 
      ewin(1,10)=0.026
      ewin(2,10)=0.2 
      ewin(1,11)=0.026
      ewin(2,11)=0.2
      dthw=10. 
c
!      read(15,*)   JJJ,ISYS,INSP,iyeld,nevtype,IGSI
c   kkg  12/12/05
      if(ztar0.ge.67.0.and.ztar0.le.88.0) then
	    call  FITAFPAQ(atar0,ztar0,ener0,fitpt(2),fitpt1(2))
      elseif(ztar0.gt.88.0)  then
	    call  FITAFACQ(atar0,ztar0,ener0,fitpt(2),fitpt1(2))
      else
            fitpt(2)  = 1.0d0
            fitpt1(2) = 1.0d0       
      endif
c
      if(zpro0.ge.67.0.and.zpro0.le.88.0) then
	    call  FITAFPAQ(apro0,zpro0,ener0,fitpt(1),fitpt1(1))
      elseif(zpro0.gt.88.0)  then
	    call  FITAFACQ(apro0,zpro0,ener0,fitpt(1),fitpt1(1))
      else
            fitpt(1)  = 1.0d0
            fitpt1(1) = 1.0d0       
      endif
c  fitaf = a_f(CEM)/a_f(RAL)
c  fitaf1 = C(Z)[CEM]/C(Z)[RAL]
c
      if(noprec.eq.1) then
        write(16,*) ' Precompound emission is excluded'
      else
        write(16,*) ' Precompound emission is included: sigpre=',
     &                sigpre
      endif
      MSYS=IABS(ISYS)
      W1=0.
      W2=0.
c      CALL  PRIJC1
c
      NC=0
      NC1=0
      INE=0
!      read(15,*) UF,DU
!      write(*,*) UF,DU
!      read(15,*) DXS,XMIN
!      write(*,*) DXS,XMIN
      DO  K=1,10
!        T1=(UF(K)  -DU)/180.*3.141592
!        if(T1.le.0.) T1=0.
!        T2=(UF(K)  +DU)/180.*3.141592
!        DOM(K) = 6.283185*(cos(T1)-cos(T2))
        DO  L=1,11
          MULT(L)=0
          DO  M=1,205
            SPE(L,11,M)=0.
            SPE(L,K,M)=0.
          ENDDO
	      ENDDO
      ENDDO
      do  k=1,5
        do  l=1,9
           pmult(k,l)=0.
        enddo
      enddo
      do  k=1,11
        do  l=1,25
          angw(l,k)=0.
        enddo
      enddo  
c  mbi
c      call  distr
c  mbi	        
c
      RETURN
      end subroutine PRSPE1
      
c ======================================================================
      
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