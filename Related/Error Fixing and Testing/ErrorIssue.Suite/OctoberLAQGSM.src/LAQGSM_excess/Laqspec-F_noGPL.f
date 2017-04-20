c ======================================================================
c
c                              NOTICE
c
c    This software and ancillary information (herein called "software")
c    named LAQGSM03.03 is made available under the terms described here.
c    The software has been approved for controlled limited release with
c    associated LA-CC number LA-CC-09-087.
c
c    Copyright (2009). Los Alamos National Security, LLC (LANS).
c    This material was produced under U. S. Government contract
c    DE-AC52-06NA25396 at Los Alamos National Laboratory, which is
c    operated by LANS for the National Nuclear Security Agency of the
c    U. S. Department of Energy.
c    The U. S. Government has rights to use, reproduce, and distribute
c    this software. NEITHER THE GOVERNMENT NOR LANS MAKES ANY WARRANTY,
c    EXPRESS OR IMPLIED OR, ASSUMES ANY LIABILITY FOR THE USE OF THIS
c    SOFTWARE.
c
c ======================================================================
c
c    The primary authors of LAQGSM03.03 are:  K. K. Gudima (LANL
c    consultant),  S. G. Mashnik (LANL), and A. J. Sierk (LANL);
c    with important contributions from R. E. Prael (LANL).
c
c ======================================================================
c
C*****************************************************************C 
*      PROGRAM LAQGSM  04/27/2007 ===> MV=5999 <====* 
c  Updated for MARS by K.K. Gudima, March-April 2007
C**INCLUDE (INTCC=1) INTERACTIONS OF CASCADE PARTICLES (M+B,B+B)
C**INCLUDE (INTCC=2) INTERACTIONS OF CASCADE PARTICLES (M+B,B+B,M+M)
C** NEW RNDM(-1.)  instead of RNDMD(-1)
c    GEM+PRECO for deexitation
c last correction in  CASCAW   10/14/03
c include gamma as projectile,  kkg 29.10.03
c KOBR=1 if ANUCL1>ANUCL2,      kkg 13.07.04
c extended to Egamma up to 10 GeV, KKG,  Dec. 2004 
c extended for bremsstrahlung gamma, KKG,Dec. 2004    
c  Dimensions of arrays XC,YC,ZC,IZ,MPA,INT1 are 300,kkg 23.11.04    
c edited by SGM, 06/08/06
c tested by KKG, 06/09/06
c edited by KKG, 06/26/06
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      REAL*8 HOUR,MINU,SECO,HUSE
      CHARACTER*4 TDATE(3),STAT
      CHARACTER*60 inptfl,FINP,FOUT,Laq1,Laq2,FTAB,FMAS,FLEV,FSHE,FGAM
c                                      ! changed from 80 !
c  kkg 04/07/04
      INTEGER*8 INTEL
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      COMMON /STIN/ STIN,AMIN
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/IFSPE/IFSPE,NCAS1
	COMMON/RNDELT/ RN,DELTA
      COMMON/JGSTAR/ JG
      COMMON/TIMGO/TIMGO,MIN,ITMGO
      common /gbmul/   gbm(4)
      common /gbrems/ tgmin, tgmax, teqv, sxabs, ibrems 
      common /filegam/ FGAM 
      common /targ0/ atar0,ztar0,ener0
      common /proj0/ apro0,zpro0
c sgm 12/04/05
      common /limita/jamin(135), jamax(135)
c  CMJ 07/2016 - adding common statements here due to moving when read
      COMMON/DEXIND/JJJ
      COMMON/VGSYS/T0A,VLA,GLA,VEV,GEV,VCM,GCM,ISYS
      COMMON/INSP/INSP
      common /indyeld/ iyeld 
      common /indgsi/ IGSI
      common /nevtyp/  nevtype
      COMMON/DOMEGA/DOM(10),UF(10),DXS(10),DU,XMIN
      DATA  TDATE/' 09/','17/','2009'/
C------------------------------------
c	WRITE( *,'('' TYPE NAME OF INPUT FILE:''/)')
C------------------------------------
				 write(*,*) 'INPUT FILE NAME: '
         read(*,*) FINP
         write(*,*)'FINP =',FINP
ccc, sgm, 02/23/06
       OPEN(15,FILE=FINP, status='old') !SGM, 06/11/08
c        OPEN(15,FILE=FINP, status='old')
         read(15,'(A60)') FOUT
        OPEN(16,FILE=FOUT,status='unknown')
c        write(16,'(A60)') FINP
c
        write (16, 600)
        write (16, 700)
c
        write(16,'(A60)') FOUT
         read(15,'(A60)') Laq1
        write(16,'(A60)') Laq1
         read(15,'(A60)') Laq2
        write(16,'(A60)') Laq2
         read(15,'(A60)') FTAB
        write(16,'(A60)') FTAB
C------------------------------------
	READ(15,'(A60)')  FMAS
	OPEN(19,FILE=FMAS,STATUS='OLD')
        WRITE(16,'(A60)') FMAS
	READ(15,'(A60)')  FLEV
	OPEN(20,FILE=FLEV,STATUS='OLD')
        WRITE(16,'(A60)') FLEV
	READ(15,'(A60)')  FSHE
	OPEN(21,FILE=FSHE,STATUS='OLD')
        WRITE(16,'(A60)') FSHE
c  kkg 03/15/04
c 	CMJ 07/2016 - removed need for FGAM in input, didn't remove from input for ease of changing
	READ(15,'(A60)')  FGAM
	FGAM = 'channel1.tab'
	OPEN(17,FILE=FGAM,STATUS='OLD')
        WRITE(16,'(A60)') FGAM
C------------------------------------
C
      READ(15,*)  JG,STAT
	open(40,FILE=Laq1,STATUS=STAT,FORM='UNFORMATTED')
	OPEN(41,FILE=Laq2,STATUS=STAT)
	OPEN(18,FILE=FTAB,STATUS='OLD')
C
C
C
      write(16,*)  'JG=',JG,' STAT=',STAT
      write( *,*)  'JG=',JG,' STAT=',STAT
      write(16,2001) TDATE
      write( *,2001) TDATE
 2001 FORMAT(2X,'PROGRAM VERSION FROM ',3A4)
C
   10 CONTINUE
      CALL  STIMER
      READ(15,*)  ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,STIN,AMIN,T0,LIMC
c   modified for bremss. gamma   12/10/04
      if(AMIN.lt.0.0001)  then
        read(15,*) ibrems
        if(ibrems.eq.1)  then
         tgmin = 0.030
c          tgmin = 0.010  ! SGM, 6/24/2011
          tgmax = T0
          teqv  = 0.0
          sxabs = 0.0
        endif  
      endif

      read(15,*) JJJ,ISYS,INSP,iyeld,nevtype,IGSI
      read(15,*) UF,DU
      read(15,*) DXS,XMIN
c
      CALL  INITGU
c  KKG 03/15/04, 12/12/05
      apro0=ANUCL1
      zpro0=ZNUCL1      
      if(ANUCL1.lt.0.1.and.AMIN.lt.0.001) then
         call inigam0(FGAM,T0)
         call inigamn(T0)
      else
	  		close(17)
      endif
c 
      if(IFSPE.NE.0.and.JG.eq.-2)  then
        CALL  PRSPE1(LIM1,LIM2,NCAS,NCAS1,INTEL,0)
        CALL  PRSPE2(LIM1,LIM2,NCAS,NCAS1,INTEL,1)
        stop
      endif
      call  iniprec
C
      CALL  RTAP10(NCAS,INTEL,NZON)

C
      IF(IFSPE.NE.0)  CALL  PRSPE1(LIM1,LIM2,NCAS,NCAS1,INTEL,JG)
c
c   kkg 01/11/05
c   Determine semiempirical multipliers gbm for gamma_beta
      if(ANUCL1.eq.1.0d0.and.ZNUCL1.eq.1.0d0)  then
        call gambetp(gbm)
      elseif(ANUCL1.eq.1.0d0.and.ZNUCL1.eq.0.0d0)  then
        call gambetn(gbm)
      else
        do  k=1,4
          gbm(k) = 1.0d0
        enddo
      endif
c
      IF(NCAS.GE.LIMC)   GO  TO  13
      MNC=NCAS+1
C
   11 continue
      IPER=0
      if(NCAS.ge.NCPRI) write( *,*) 'start Ncas=',NCAS
c
c      write(*,*) ' curent NCAS,INTEL=',NCAS,INTEL
c
      CALL CASCAW(NEL,RN,DELTA,MV,IRET)
      IF(IRET.EQ.1)  GO  TO  11
      IF(IPER.EQ.1)  GO  TO  11
c sgm, 12/04/05
      jz1 = nint(ZN1) 
      jz2 = nint(ZN2)
      ja1 = nint(AN1) 
      ja2 = nint(AN2)

c      write(*,*) ' NCAS= ',NCAS,
c     *' Z1=',jz1,' A1=',ja1,'[',jamin(jz1),'-',jamax(jz1),
c     *'] Z2=',jz2,' A2=',ja2,'[',jamin(jz2),'-',jamax(jz2),']'
c      write(16,*) ' NCAS= ',NCAS,
c     *' Z1=',jz1,' A1=',ja1,'[',jamin(jz1),'-',jamax(jz1),
c     *'] Z2=',jz2,' A2=',ja2,'[',jamin(jz2),'-',jamax(jz2),']'

      if((ja1.ge.12.and.jz1.eq.0).or.(ja2.ge.12.and.jz2.eq.0)) go to 11
      if(jz1.eq.0.or.jz2.eq.0) go to 1111
      if((ja1.ge.12.and.ja1.lt.jamin(jz1)).or.
     *   (ja1.ge.12.and.ja1.gt.jamax(jz1))) go to 11
      if((ja2.ge.12.and.ja2.lt.jamin(jz2)).or.
     *   (ja2.ge.12.and.ja2.gt.jamax(jz2))) go to 11
1111  continue
c sgm, 12/04/05

      INTEL=INTEL+NEL
      IF(ENEXT2.GT.1.E-6.AND.IRET.EQ.0)    NCAS=NCAS+1
      IF(IFSPE.NE.0)  CALL DIFSPE(NCAS,INTEL,NCAS1,MV)
      CALL  COLEC(NCAS,INTEL,MV,NZON,LUP,0)
      IF(IFSPE.NE.0)  CALL PRSPE2(LIM1,LIM2,NCAS,NCAS1,INTEL,0)
      IF(NCAS.GE.LIMC)   then
        GO  TO  12
      ELSE
        GO  TO  11
      ENDIF
   12 CONTINUE
      CALL  COLEC(NCAS,INTEL,MV,NZON,LUP,1)
C
      IF(IFSPE.NE.0)  CALL  PRSPE2(LIM1,LIM2,NCAS,NCAS1,INTEL,1)
c
   13 CONTINUE
      NEVENT=NCAS-MNC+1
      CALL  DTIMER(HOUR,MINU,SECO,HUSE,1)
      TOTT=3600.*HOUR+60.*MINU+SECO+HUSE/100.
      ETIM=TOTT
      IF(NEVENT.NE.0)  THEN
      	 ETIM=TOTT/DBLE(NEVENT)
        write( *,701)   HOUR,MINU,SECO,HUSE,ETIM
        write(16,701)   HOUR,MINU,SECO,HUSE,ETIM
 701    FORMAT(//2X,'EXEC TIME : ',4(F3.0,1X)//
     *  2X,'TIME PER EVENT=',F10.3,' SEC')
      ENDIF
C
c ======================================================================
  600 format (/30x,'NOTICE'//3x,'This software and ancillary ',
     &       'information (herein called "software") named'/3x,'LAQGSM',
     &       '03.03 is made available under the terms described here. ',
     &       'The software has'/3x,'been approved for release with ',
     &       'associated LA-CC number LA-CC-09-087.'//3x,
     &       'Copyright (2009). Los Alamos National Security, LLC '
     &       ,'(LANS).'/3x,'This material was produced under U. S. ',
     &       'Government contract'/3x,'DE-AC52-06NA25396 at Los Alamos '
     &       ,'National Laboratory, which is operated'/3x,'by LANS for',
     &       ' the National Nuclear Security Agency of the U. S. ',
     &       'Department of Energy.'/3x,'The U. S. Government has ',
     &       'rights to use, reproduce, and distribute'/3x,'this ',
     &       'software. NEITHER THE GOVERNMENT NOR LANS MAKES ANY'/3x,
     &       'WARRANTY, EXPRESS OR IMPLIED OR ASSUMES ANY LIABILITY ',
     &       'FOR THE USE'/3x,'OF THIS SOFTWARE.')
  700 format (/3x,'Primary authors of LAQGSM03.03 are: K. K. Gudima ',
     &        '(LANL consultant),'/3x,'S. G. Mashnik (LANL), and ',
     &        'A. J. Sierk (LANL); with ',
     &        'important contributions'/3x,'from R. E. Prael (LANL),'
     &        /3x,'---------------------------------------------------',
     &        '----------------------------'/)

      CLOSE(40)
      CLOSE(41)
      CLOSE(16)
C
      STOP
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
