c  This is code from the Laq03-F_noGPL.f file used in LAQGSM.
c  This code was deemed not necessary to be used within CEM, 
c  but is not being deleted for the time being to help aid
c  in merging CEM and LAQGSM.
c
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
    SUBROUTINE  COLEC(NCAS,INTEL,MV,NZON,LUP,ICOL)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c   
      REAL*8 MMES,MBAR
c  kkg 04/07/04
      INTEGER*8 INTEL 
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,
     *A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,EN1,EN2,PN1(3),PN2(3),AM1(3),AM2(3)
      COMMON/BIMP/BIMP,BIMX,BIMY/KYPT/KPT,KYP,KYT,KPT1,KYP1,KYT1,KPZ,KTZ
      COMMON/ACTIM/TINT
      COMMON/PIDABS/PIDABS
      COMMON/RCOR/ RCOR
      COMMON/MMATUR/MMES,MBAR
      COMMON /STIN/ STIN,AMIN
      COMMON/INDINT/ INDINT
      COMMON/INTCC/ INTCC
      COMMON /IACT/ IACT/INDDEC/INDDEC/ISOB3/ISOB3/ISOB2/ISOB2
      COMMON/KSYST/ KSYST
      COMMON/EXCIT/TEX1,TEX2,HEX1,HEX2,PEX1,PEX2
      COMMON/UPAC/UP1(66000),IPER/UPAC2/UP2(70000),LU2
      COMMON /NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE,
     *VEP(3),VET(3),GEP,GET
      COMMON/INDZAP/INDZAP
      COMMON/BGLAB/BLAB,GLAB,KOBR
      INDZAP=0
      IF(ICOL.EQ.1)   GO  TO  604
C
      IF((LU2+50+MV*11+1).GE.70000) GO  TO  604
C
  600 UP2(LU2)=DBLE(MV)
      UP2(LU2+1)=AN1
      UP2(LU2+2)=AN2
      UP2(LU2+3)=ZN1
      UP2(LU2+4)=ZN2
      UP2(LU2+5)=EN1
      UP2(LU2+6)=EN2
      DO  601  KL=1,3
      UP2(LU2+6+KL)=PN1(KL)
      UP2(LU2+9+KL)=PN2(KL)
      UP2(LU2+12+KL)=AM1(KL)
  601 UP2(LU2+15+KL)=AM2(KL)
      IF(KOBR.NE.1)  GO  TO  9
      IF(AN1.LT.0.1)  GO  TO  8
      EA1=SQRT(PN1(1)**2+PN1(2)**2+PN1(3)**2+(0.94*AN1)**2)
      UP2(LU2+9)=-GLAB*(PN1(3)-BLAB*EA1)
    8 IF(AN2.LT.0.1)  GO  TO  9
      EA2=SQRT(PN2(1)**2+PN2(2)**2+PN2(3)**2+(0.940*AN2)**2)
      UP2(LU2+12)=-GLAB*(PN2(3)-BLAB*EA2)
    9 CONTINUE
      UP2(LU2+19)=TEX1
      UP2(LU2+20)=TEX2
      UP2(LU2+21)=HEX1
      UP2(LU2+22)=HEX2
      UP2(LU2+23)=PEX1
      UP2(LU2+24)=PEX2
      UP2(LU2+25)=FLOAT(KPT)
      UP2(LU2+26)=FLOAT(KYP)
      UP2(LU2+27)=FLOAT(KYT)
      UP2(LU2+28)=FLOAT(KPT1)
      UP2(LU2+29)=FLOAT(KYP1)
      UP2(LU2+30)=FLOAT(KYT1)
      UP2(LU2+31)=FLOAT(KPZ)
      UP2(LU2+32)=FLOAT(KTZ)
      UP2(LU2+33)=BIMP
      UP2(LU2+34)=BIMX
      UP2(LU2+35)=BIMY
      UP2(LU2+36)=TINT
      DO  10  LB=1,3
      UP2(LU2+36+LB)=RADP(LB)
   10 UP2(LU2+39+LB)=RADT(LB)
      UP2(LU2+43)=PIDABS
      UP2(LU2+44)=RCOR
      UP2(LU2+45)=MMES
      UP2(LU2+46)=MBAR
      UP2(LU2+47)=STIN
      UP2(LU2+48)=AMIN
      UP2(LU2+49)=DBLE(10*INDINT+INTCC)
      UP2(LU2+50)=DBLE(KSYST*10000+IACT*1000+INDDEC*100+ISOB2*10+
     *ISOB3)
      MVK=11*MV
      IF(MV.EQ.0)   GO  TO  603
      DO  602  KL=1,MVK
      UP2(LU2+50+KL)=UP1(KL)
  602 CONTINUE
  603 LU2=LU2+50+MVK+1
      GO  TO  606
  604 IF(LU2.EQ.1)   GO  TO  606
      INDZAP=1
      CALL  WTAP10(NCAS,INTEL,NZON)
      IF(ICOL.EQ.0)  GO  TO  600
  606 CONTINUE
      LUP=LU2
      RETURN
    END
cc   **************************************************************
c
      SUBROUTINE RTAP10(NCAS,INTEL,NZON)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c  kkg 04/07/04
      INTEGER*8 INTEL,IS7 
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/UPAC/UP1(66000),IPER
      COMMON/UPAC2/UP2(70000),LU2
      COMMON/JGSTAR/ JG
      COMMON /IW10/ IW10
      COMMON/IFSPE/ IFSPE,NCAS1
      common /gbrems/ tgmin, tgmax, teqv, sxabs, ibrems 
      NUNIT=2000
      NZON=0
      LU2=1
      REWIND 40
      IF(JG.eq.0)   JG=500000
      IF(JG.GT.0) then
      DO  KG=1,JG
c      write( *,*) ' KG=', KG
        READ(40,END=98) S1,S2,S3,S4,S5,IS6,IS7,IS8,IS91,IS92,IS93,IS10
C  !!!
        if(IW10.EQ.1)  then
C  !!!
          NUN=IS10/NUNIT
          IF(NUN.ge.1) THEN
            do  LUN=1,NUN
              K1=(LUN-1)*NUNIT+1
              K2=LUN*NUNIT
              READ(40) (UP2(KL),KL=K1,K2)
            enddo
          ENDIF
          K3=NUN*NUNIT+1
c          write( *,*) 'K3,IS10,NUN=', K3,IS10,NUN
          IF(K3.LE.IS10) READ(40) (UP2(KL),KL=K3,IS10)
        endif
        NZON=NZON+1
      ENDDO
      GO TO 99
   98 CONTINUE
c      BACKSPACE 10
 99   A10=S1
      A20=S2
      Z10=S3
      Z20=S4
      T00=S5
      NCAS=IS6
      INTEL=IS7
      NZON=IS8
      IS100=IS10
      write(16,'(1X,''FROM10: '',4(F5.0,1X),F7.3,1X,4I9)')
     & A10,A20,Z10,Z20,T00,NCAS,INTEL,NZON,IS100
      write( *,'(1X,''FROM10: '',4(F5.0,1X),F7.3,1X,4I9)')
     & A10,A20,Z10,Z20,T00,NCAS,INTEL,NZON,IS100
      JG=NZON        
        IF(ABS(ANUCL1-A10).le.0.1.and.
     &  ABS(ANUCL2-A20).le.0.1.and.
     &  ABS(ZNUCL1-Z10).le.0.1.and.
     &  ABS(ZNUCL2-Z20).le.0.1.and.
     &  ABS(T0-T00).le. 0.0001.or.ibrems.eq.1)     then
          IBRAN=IS91
c         write(*,*) ' IBRAN=',IBRAN
          CALL  RDMIN(IS91,IS92,IS93)
        else
          NCAS=0
          INTEL=0
          JG=-1
        endif
      ELSE
        NCAS=0
        INTEL=0
        JG=-1
      ENDIF
      NZONS=NZON+1
      write(16,'(/10X,''START NZON='',I6,'' NCAS='',I9,'' INTEL='',I9)')
     & NZONS,NCAS,INTEL
      write( *,'(/10X,''START NZON='',I6,'' NCAS='',I9,'' INTEL='',I9)')
     & NZONS,NCAS,INTEL
      RETURN
C---------------------------------------------
      ENTRY WTAP10(NCAS,INTEL,NZON)
c
      LU2=LU2-1
      NZON=NZON+1
      CALL  RDMOUT(IBR1,IBR2,IBR3)
      S1=ANUCL1
      S2=ANUCL2
      S3=ZNUCL1
      S4=ZNUCL2
      S5=T0
      WRITE(40) S1,S2,S3,S4,S5,NCAS,INTEL,NZON,IBR1,IBR2,IBR3,LU2
C  !!!
      IF(IW10.EQ.1) THEN
C  !!!
        NUN=LU2/NUNIT
        IF(NUN.ge.1)  then
          DO  LUN=1,NUN
            K1=(LUN-1)*NUNIT+1
            K2=LUN*NUNIT
            WRITE(40) (UP2(KL),KL=K1,K2)
          ENDDO
        endif
        K3=NUN*NUNIT+1
        IF(K3.LE.LU2) WRITE(40) (UP2(KL),KL=K3,LU2)
      ENDIF
c     write(16,'(1X,4(F5.0,1X),F8.3,1X,3(I10,1X),I7)')
c    & ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,NCAS,INTEL,NZON,LU2
      write( *,'(1X,4(F5.0,1X),F8.3,1X,3(I10,1X),I7)')
     & ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,NCAS,INTEL,NZON,LU2
c     write( *,*) ' WRITE on 10: IBR1=', IBR1
      LU2=1
      RETURN
      END
C
C  *************************************************************
C   * * * * * * * * * * * * * * * * * * * * * * *

    subroutine  PRSPE2(LIM1,LIM2,NC,NC1,INE,IW11)
!      ENTRY  PRSPE2(LIM1,LIM2,NC,NC1,INE,IW11)
			write(*,*) "Into PRSPE2"
			!stop
c
      if(INDZAP.eq.1)  then
        if(ibrems.eq.0)  then
          WRITE(41,*) NC,NC1,INE,MULT,NF1,NF2,SPE,pmult,dadz,angw,
     &    W1,W2,dadz4,rdis,opan,disnm
        else
          WRITE(41,*) NC,NC1,INE,MULT,NF1,NF2,SPE,pmult,dadz,angw,teqv,
     &    W1,W2,dadz4,rdis,opan,disnm
        endif
        REWIND 41
      endif
      if(IW11.eq.0)  RETURN
c
      SIN=SIG1*FLOAT(NC1)/FLOAT(NC+INE)
      if(ibrems.eq.1)  then
        eqqv = teqv/tgmax/float(NC1)
        sineq=SIN/eqqv
        SIN = sineq
        write(16,'(/5X,''Ncas='',I8,1X,''Selected Ncas='',I8,1X,
     &''INTel='',I11/
     &  5X,''Number of eq. quanta per inelastic event='',E13.6,'' mb''/
     &  5X,''INEL. CROSS SECTION per eqqv='',E13.6,'' mb''/
     &''  All present results are normalized on eqqv'')') 
     &   NC,NC1,INE,eqqv,SIN
       else
        write(16,'(/5X,''Ncas='',I8,1X,''Selected Ncas='',I8,1X,
     &  ''INTel='',I11/5X,''INEL. CROSS SECTION ='',E13.6,'' mb'')')
     &   NC,NC1,INE,SIN
       endif 
c
      WNF=FLOAT(NC1)
      SIGF1=SIn*FLOAT(NF1)/WNF
      SIGF2=SIN*FLOAT(NF2)/WNF
      dsfi1=SIN*SQRT(FLOAT(NF1))/WNF
      dsfi2=SIN*SQRT(FLOAT(NF2))/WNF
      SIGF1w=SIN*W1/WNF
      SIGF2w=SIN*W2/WNF
	write(16,'(/10X,''FISSION  Probabilities '',      
     &''Projectile='',1PE13.6,1X,''Target='',1PE13.6)')
     & W1/WNF,W2/WNF 
      write(16,'(/1X,''FISSION  CROSS  SECTIONS(direct  ) '', 
     &''SIGFISP='',F13.6,''+-'',F13.6,1X,
     &''SIGFIST='',F13.6,''+-'',F13.6)')
     & SIGF1,dsfi1,SIGF2,dsfi2
      write(16,'(/1X,''FISSION  CROSS  SECTIONS(weighted) '', 
     &''SIGFISP='',F13.6,1X,''SIGFIST='',F13.6)')
     & SIGF1w,SIGF2w
ccc      
      write(16,3400)
 3400 format(5x,'Selected particle multiplicity'/
     &       5x,'   preco   ','   evapor  ','   Fermi   ',
     &          'casc+coales','   summ   ')
      do  ip=1,9
        apm(5)=0.0
        do  io=1,4
          apm(io)=pmult(io,ip)/FLOAT(NC1)
          apm(5) =apm(5)+apm(io)
        enddo
        if(ip.eq.1) write(16,3401) apm
 3401   format(1x,' n  ',5(1PE11.4))
        if(ip.eq.2) write(16,3402) apm
 3402   format(1x,' p  ',5(1PE11.4))
        if(ip.eq.3) write(16,3403) apm
 3403   format(1x,' d  ',5(1PE11.4))
        if(ip.eq.4) write(16,3404) apm
 3404   format(1x,' t  ',5(1PE11.4))
        if(ip.eq.5) write(16,3405) apm
 3405   format(1x,'He-3',5(1PE11.4))
        if(ip.eq.6) write(16,3406) apm
 3406   format(1x,'He-4',5(1PE11.4))
        if(ip.eq.7) write(16,3407) apm
 3407   format(1x,'pi- ',5(1PE11.4))
        if(ip.eq.8) write(16,3408) apm
 3408   format(1x,'pi+ ',5(1PE11.4))
        if(ip.eq.9) write(16,3409) apm
 3409   format(1x,'pi0 ',5(1PE11.4))
      enddo
c
c     LP=1
      LP=2
  100 CONTINUE
      DO    JP=1,11
        AMUL(JP)=FLOAT(MULT(JP))/FLOAT(NC1)
      ENDDO
      IF(LP.EQ.1) then
        WRITE(16,101) TEXT
  101   FORMAT(10x,'Number of produced particles'/1x,11(3x,2A4))
        WRITE(16,102) MULT
  102   FORMAT(1x,11(2x,I9))
      ELSE
        WRITE(16,103) TEXT
  103   FORMAT(/10x,'Multiplicity of produced particles'/
     &  1x,11(3x,2A4))
        WRITE(16,104) AMUL
  104   FORMAT(1x,11(1PE11.4))
      ENDIF
      if(INSP.eq.0)  go  to   2002  
      its=1+(INSP-1)*2
      DO  IU=1,10
        IF(LP.EQ.1) then
          WRITE(16,105) UF(IU),DU,TEX3(its),TEX3(its+1),TEXT
  105     FORMAT(/10x,'Number of produced particles at theta=',
     &    F6.2,'+/-',F6.2/1x,2A4,11(3x,2A4))
        ELSE
          WRITE(16,106) TEX2(its),TEX2(its+1),UF(IU),DU,
     &                  TEX3(its),TEX3(its+1),TEXT
  106     FORMAT(/10x,'Spectra ( ',2A4,' )',
     &    ' of produced particles at theta=',F6.2,'+/-',F6.2/
     &    1x,2A4,11(3x,2A4))
        ENDIF
        DO  k=1,205
          X= (k-1)*DXS(IU)+DXS(IU)/2.+XMIN
          WX=1.
          cps=0.0
          DO  JP=1,11
            if(INSP.eq.1)  then
              xa = X
c             if(JP.ge.5)  xa = xa*NBA(JP)   ! 03/28/05
              WX=SQRT(xa**2+MA(JP)**2)/xa**2
              DX1=DXS(IU)
c             if(JP.ge.5)  DX1=DX1*NBA(JP)   ! 03/28/05          
              IF(LP.eq.1)  THEN
               C(JP)=SPE(JP,IU,k)
              ELSE
               C(JP)=SPE(JP,IU,k)*SIN/FLOAT(NC1)/DOM(IU)/DX1*WX
              ENDIF
            elseif(INSP.eq.2)  then
              WX=1.
              if(T0A.lt.T0le)  then
                if(k.ge.jgr(1).and.k.lt.jgr(2))            then
                  dx2=dtgr(1)
                  X=(k-jgr(1))*dx2+dx2/2.
                elseif(k.ge.jgr(2).and.k.lt.jgr(3))        then
                  dx2=dtgr(2)
                  X=(k-jgr(2))*dx2+dx2/2.+tgr(1)
                elseif(k.ge.jgr(3).and.k.lt.jgr(4))        then
                  dx2=dtgr(3)
                  X=(k-jgr(3))*dx2+dx2/2.+tgr(2)
                elseif(k.ge.jgr(4))                        then
                  dx2=dtgr(4)
                  X=(k-jgr(4))*dx2+dx2/2.+tgr(3)
                else
                endif
              else
                dx2=DXS(IU)
              endif
              IF(LP.eq.1)  THEN
               C(JP)=SPE(JP,IU,k)
              ELSE
               C(JP)=SPE(JP,IU,k)*SIN/FLOAT(NC1)/DOM(IU)/dx2*WX
              ENDIF
            else
              IF(LP.eq.1)  THEN
               C(JP)=SPE(JP,IU,k)
              ELSE
               WX=1.
               C(JP)=SPE(JP,IU,k)*1.0/FLOAT(NC1)/DOM(IU)/DXS(IU)*WX
              ENDIF
            endif
            if(k.gt.200) C(JP)=SPE(JP,IU,k)
            if(k.le.200) cps=cps+C(JP)
          ENDDO
          IF(LP.eq.1)  THEN
            if(k.le.200.and.cps.gt.0.0) WRITE(16,107) X,C
  107       FORMAT(F9.4,11(F11.0))
c           if(k.gt.200) WRITE(16,108)   C
  108       FORMAT(9x,11(F11.0))
          ELSE
            if(k.le.200.and.cps.gt.0.0) WRITE(16,109) X,C
  109       FORMAT(F9.4,11(1PE11.4))
c           if(k.gt.200) WRITE(16,110)   C
  110       FORMAT(9x,11(F11.0))
          ENDIF
        ENDDO
      ENDDO
C
      IF(INSP.EQ.2)  THEN
        IF(LP.EQ.1)  then
         WRITE(16,111) TEX3(its),TEX3(its+1),TEXT
  111    FORMAT(/10X,'Angle inegrated number of produced particles'/
     &   1X,2A4,11(3X,2A4))
        ELSE
         WRITE(16,112) TEXT
  112    FORMAT(/10X,'Angle integrated spectra of produced particles'/
     &   9X,11(3X,2A4))
        ENDIF
        DO  k=1,205
          X=(k-1)*DXS(1)+DXS(1)/2.+XMIN
          WX=1.
          cps=0.0
          DO  JP=1,11
            if(T0A.lt.T0le)  then
              if(k.ge.jgr(1).and.k.lt.jgr(2))            then
                dx2=dtgr(1)
                X=(k-jgr(1))*dx2+dx2/2.
              elseif(k.ge.jgr(2).and.k.lt.jgr(3))        then
                dx2=dtgr(2)
                X=(k-jgr(2))*dx2+dx2/2.+tgr(1)
              elseif(k.ge.jgr(3).and.k.lt.jgr(4))        then
                dx2=dtgr(3)
                X=(k-jgr(3))*dx2+dx2/2.+tgr(2)
              elseif(k.ge.jgr(4))                        then
                dx2=dtgr(4)
                X=(k-jgr(4))*dx2+dx2/2.+tgr(3)
              else
              endif
            else
              dx2=DXS(1)           ! 17.04.2003
            endif
            IF(LP.EQ.1)  THEN
              C(JP)=SPE(JP,11,k)
            ELSE
              C(JP)=SPE(JP,11,k)*SIN/FLOAT(NC1)/dx2*WX
              if(k.gt.200)  C(JP)=SPE(JP,11,k)
            ENDIF
            if(k.le.200) cps=cps+C(JP)
          ENDDO
          IF(LP.EQ.1)  THEN
            if(k.le.200.and.cps.gt.0.0) WRITE(16,107) X,C
c           if(k.gt.200) WRITE(16,108)   C
          ELSE
            if(k.le.200.and.cps.gt.0.0) WRITE(16,109) X,C
c           if(k.gt.200) WRITE(16,110)   C
          ENDIF
        ENDDO
      ENDIF
C
c
      IF(LP.EQ.1)  then
       WRITE(16,113) TEXT
  113  FORMAT(/10X,'Energy inegrated number of produced particles'/
     & 1X,'theta ',11(3X,2A4))
      ELSE
       WRITE(16,114) TEXT
  114  FORMAT(/10X,'Energy integ. ang. distrib. of prod. particles'/
     & 1X,'theta ',11(3X,2A4))
      ENDIF
      do  k=1,18
        x=(k-1)*dthw+dthw/2.
        x1=(x-dthw/2.)*3.141592/180.  
        x2=(x+dthw/2.)*3.141592/180.  
        dome=6.283185*(COS(x1)-COS(x2))
        do  l=1,11
          if(LP.eq.1)  then
            CC(l)=angw(k,l)
          else
            CC(l)=angw(k,l)*SIN/FLOAT(NC1)/dome
          endif
        enddo
        if(LP.eq.1)  then
          write(16,115) x,CC
  115     FORMAT(F9.4,11(F11.0))
        else
          write(16,116) x,CC
  116     FORMAT(F9.4,11(1PE11.4))
        endif 
      enddo       
c        
      LP=LP+1
      IF(LP.LE.2)  GO  TO  100
 2002 CONTINUE
      if(iyeld.eq.0)  go  to  2003
c
c   Print out isotope formation cross sections:
c
c  kkg
        abeg=atar0 + 5.
	  zbeg=ztar0 + 5.
        izmax=int(zbeg+0.1)
        if(izmax.gt.149)  izmax=149
        iamax=int(abeg+0.1)
        if(iamax.gt.349)  iamax=349
        do 590 iz = 1,izmax
          z(iz) = max (zro, zbeg - real(iz-1))
          dasum(iz) = zro
          do 580 ia = 1,iamax
            ar(ia) = max( abeg - real(ia-1), zro)
            dasum(iz) = dasum(iz) + dadz(ia,iz)
  580     continue
        prtz(iz) = z(iz).ge.zro .and. dasum(iz).gt.zro
  590   continue
        write (16, 3500)
        lprm=(izmax+3)/3
        if(lprm.lt.1) lprm=1
        do 610 lpr=1,lprm
          iz1=3*lpr-2
          iz2=3*lpr-1
          iz3=3*lpr
          if (prtz(iz3).or.prtz(iz2).or.prtz(iz1)) then
            write (16, 3600) z(iz1), z(iz2), z(iz3)
c         elseif (prtz(iz2)) then
c           write (16, 3610) z(iz1), z(iz2)
c         elseif (prtz(iz1)) then
c           write (16, 3620) z(iz1)
          endif
          fac = SIN/FLOAT(NC1)
          ppt=zro
          do 600 k = 1,iamax
            iaa = nint(ar(k))
            y1 = dadz(k,iz1)*fac
            dy1 = sqrt(abs(dadz(k,iz1)))*fac
            y2 = dadz(k,iz2)*fac
            dy2 = sqrt(abs(dadz(k,iz2)))*fac
            y3 = dadz(k,iz3)*fac
            dy3 = sqrt(abs(dadz(k,iz3)))*fac
            printit = y1.gt.zro .or. y2.gt.zro .or. y3.gt.zro
            if (printit) then
              if (ar(k).gt.zro .and. ar(k).gt.z(iz3)) then
                if (prtz(iz3).or.prtz(iz2).or.prtz(iz1)) then
                   ppt=ppt+1.
                   write (16, 3700) iaa, y1, dy1, y2, dy2, y3, dy3
c               elseif (prtz(iz2)) then
c                  write (16, 3710) iaa, y1, dy1, y2, dy2
c               elseif (prtz(iz1)) then
c                  write (16, 3720) iaa, y1, dy1
                endif
              endif
            endif
  600     continue
          if(ppt.gt.zro)  then
            ys1 = dasum(iz1)*fac
            dys1= sqrt(abs(dasum(iz1)))*fac
            ys2 = dasum(iz2)*fac
            dys2= sqrt(abs(dasum(iz2)))*fac
            ys3 = dasum(iz3)*fac
            dys3= sqrt(abs(dasum(iz3)))*fac
            ipt=nint(ppt)
            write(16,3701) ipt, ys1,dys1, ys2,dys2, ys3,dys3
          endif 
  610   continue
        write(16,3702)
 3500 format (//1x,'*************** Nuclide yields [mb]  (zero values ',
     &       'suppressed) *****************'/1x,
     &       '(normalized to the Total Reaction Cross Section.)')
 3600 format (/17x,'Z = ',f4.0,16x,'Z = ',f4.0,16x,'Z = ',f4.0)
 3610 format (/17x,'Z = ',f4.0,16x,'Z = ',f4.0)
 3620 format (/17x,'Z = ',f4.0)
 3700 format (1x,'A =',i4,3(1x,1pe9.3,' +/- ',e8.2))
 3701 format (1x,'S =',i4,3(1x,1pe9.3,' +/- ',e8.2))
 3702 format (1x,'End of nuclide yields.') 
 3710 format (1x,'A =',i4,2(1x,1pe9.3,' +/- ',e8.2))
 3720 format (1x,'A =',i4,1(1x,1pe9.3,' +/- ',e8.2))
c
        ip = nint(atar0) + 5
        ipz = nint(ztar0) + 5
        do 980 i=1,ip
980     iar(i) = ip - i + 1
        do 915 i=1,ipz
915     izr(i) = ipz - i + 1

        ss = fac
        write(16,951)
951     format(/'Mass-Yield [mb] and Kinetic Energy [MeV]',
     *  ' Distributions of Residual Nuclei')
        rese=0.
	  yys=0.
	  reses=0.
	  npt=0
        do 960 i=1,ip
          yy=dadz(i,150)*ss
          dyy=sqrt(abs(dadz(i,150)))*ss
	    yys=yys+dadz(i,150)
	    reses=reses+dadz(i,151)
          if(dadz(i,150).gt.0.) then
             rese = dadz(i,151)/dadz(i,150)
	    else
	       rese = 0.
	    endif
          printit = yy.gt.zro
	    if (printit) then
	       npt=npt+1
             write(16,952) iar(i),yy,dyy,rese
          endif
960     continue
        dyys=sqrt(abs(yys))*ss
	  if(yys.gt.0) reses=reses/yys
	  yys=yys*ss
        write(16,955) npt,yys,dyys,reses
955     format(1x,'S =',I4,1X,1PE9.3,'+/-',1PE8.2,5X,1PE9.3)
952     FORMAT(1X,'A =',I4,1X,1PE9.3,'+/-',1PE8.2,5X,1PE9.3)
        write(16,953)
953     format(/'Charge-Yield [mb] and Kinetic Energy [MeV]',
     *  ' Distributions of Residual Nuclei')
	  yys=0.
	  reses=0.
	  npt=0
        do 961 i=1,ipz
          yy=dadz(351,i)*ss
          dyy=sqrt(abs(dadz(351,i)))*ss
	    yys=yys+dadz(351,i)
	    reses=reses+dadz(350,i)
          if(dadz(351,i).gt.0.) then
            rese = dadz(350,i)/dadz(351,i)
	    else
	      rese = 0.
	    endif
          printit = yy.gt.zro
	    if (printit) then
	      npt=npt+1 
            write(16,954) izr(i),yy,dyy,rese
	    endif
961     continue
        dyys=sqrt(abs(yys))*ss
	  if(yys.gt.0) reses=reses/yys
	  yys=yys*ss
        write(16,956) npt,yys,dyys,reses
954     FORMAT(1X,'Z =',I3,1X,1PE9.3,'+/-',1PE8.2,5x,1PE9.3)
956     FORMAT(1X,'S =',I3,1X,1PE9.3,'+/-',1PE8.2,5x,1PE9.3)
c  kkg
 2003 continue
      fn = dble(NC1)
      if(iyeld.ge.2) call prdadzq(atar0,ztar0,0,SIN,fn)
      if(iyeld.eq.3) call   prrdisq(atar0,ztar0,0,fn)
      if(iyeld.eq.3.and.NF2.ne.0)  call propanq(fn)
      if(iyeld.eq.3) call  pdisnmq(fn) 
c
      WRITE(16,*) ' to 11:  NC=',NC,' NC1=',NC1,' INE=',INE
      WRITE( *,*) ' to 11:  NC=',NC,' NC1=',NC1,' INE=',INE
      RETURN
    END
    
C  *************************************************************
    
      subroutine prdadzq(anucl,znucl,mb0,sigin,fn)

c ======================================================================
c
c   This subroutine prints out average kinetic  energies  of produced 
c   fragments and fragment yields emitted if forward/backward direction
c
c   kkg 11/13/04
      implicit real*8 (a-h, o-z), integer (i-n)
c
      logical printit, prtz(151), prtz1, prtz2, prtz3
      common /yeldaz/ dadz(351,151)
      common /yeldaz1/ dadz4(4,351,151)
      dimension ar(351), z(151)
      dimension dasum(151)
      dimension iar(351), izr(151)
      data zro, one, two, mil /0.d0, 1.d0, 2.d0, 1.d3/

      abeg = anucl + 5.d0
      zbeg = znucl + 5.d0
      izmax = nint(zbeg)
      izmax = min(izmax, 148)
      iamax = nint(abeg) + mb0
      iamax = min(iamax, 348)
      lprm = (izmax + 3)/3
      lprm = max(lprm, 1)
      faca = sigin/fn
      ip = nint(abeg)  
      ipz = nint(zbeg) 
      do i = 1,ip
        iar(i) = ip - i + 1
      end do
      do i= 1,ipz
        izr(i) = ipz - i + 1
      end do
c
c   Print out forward isotope formation cross sections:
c
      do iz = 1,izmax + 1
        z(iz) = max (-one, zbeg - dble(iz-1))
        dasum(iz) = zro
        do ia = 1,iamax
          ar(ia) = max(abeg - dble(ia-1), zro)
          dasum(iz) = dasum(iz) + dadz4(3,ia,iz)
        end do
        prtz(iz) = z(iz).ge.zro .and. dasum(iz).gt.zro
      end do
      write (16, 3701)
      do lpr = 1,lprm
        iz1 = 3*lpr - 2
        iz2 = 3*lpr - 1
        iz3 = 3*lpr
        prtz1 = prtz(iz3)
        prtz2 = .not.prtz(iz3) .and. prtz(iz2)
        prtz3 = (.not.prtz(iz2) .and. .not.prtz(iz3)) .and. prtz(iz1)
        if (prtz1) then
          write (16, 3800) z(iz1), z(iz2), z(iz3)
        elseif (prtz2) then
          write (16, 3900) z(iz1), z(iz2)
        elseif (prtz3) then
          write (16, 4000) z(iz1)
        endif
        ppt = zro
        do k = 1,iamax
          y1 = dadz4(3,k,iz1)*faca
          dy1 = sqrt(abs(dadz4(3,k,iz1)))*faca
          y2 = dadz4(3,k,iz2)*faca
          dy2 = sqrt(abs(dadz4(3,k,iz2)))*faca
          y3 = dadz4(3,k,iz3)*faca
          dy3 = sqrt(abs(dadz4(3,k,iz3)))*faca
          printit = y1.gt.zro .or. y2.gt.zro .or. y3.gt.zro
          if (printit) then
c  kkg 11/15/04
            if (ar(k).gt.zro .and. ar(k).ge.z(iz3)) then
              iaa = nint(ar(k))
              prtz1 = prtz(iz3)
              prtz2 = .not.prtz(iz3) .and. prtz(iz2)
              prtz3 = (.not.prtz(iz2) .and. .not.prtz(iz3)) .and. 
     &            		  prtz(iz1)
              if (prtz1) then
                ppt = ppt + one
                write (16, 4100) iaa, y1, dy1, y2, dy2, y3, dy3
              elseif (prtz2) then
                ppt = ppt + one
                write (16, 4200) iaa, y1, dy1, y2, dy2
              elseif (prtz3) then
                ppt = ppt + one
                write (16, 4300) iaa, y1, dy1
              endif
            endif
          endif
        end do
        if (ppt.gt.zro) then
          ipt = nint(ppt)
          ys1  = dasum(iz1)*faca
          dys1 = sqrt(abs(dasum(iz1)))*faca
          ys2  = dasum(iz2)*faca
          dys2 = sqrt(abs(dasum(iz2)))*faca
          ys3  = dasum(iz3)*faca
          dys3 = sqrt(abs(dasum(iz3)))*faca
          if (ys3.gt.zro) then
            write (16, 4400) ipt, ys1, dys1, ys2, dys2, ys3, dys3
          elseif (ys2.gt.zro) then
            write (16, 4500) ipt, ys1, dys1, ys2, dys2
          elseif (ys1.gt.zro) then
            write (16, 4600) ipt, ys1, dys1
          endif 
        endif 
      end do
      write (16, 4701)
      ss = faca
      write (16, 1101)
	atks = zro
	datks= zro
      yys = zro
      npt = 0
      do i = 1,ip
        yy = dadz4(3,i,150)*ss
        dyy = sqrt(abs(dadz4(3,i,150)))*ss
        yys = yys + dadz4(3,i,150)
        if (dadz4(3,i,150).gt.zro) then
          atk = dadz4(3,i,151)/dadz4(3,i,150)
	    atks= atks + dadz4(3,i,151)
	    datk = sqrt(abs(dadz4(3,i,149)/dadz4(3,i,150)-atk**2))
		datks= datks + dadz4(3,i,149)  
        else
          atk = zro
	    datk= zro
        endif
        printit = yy.gt.zro 
        if (printit) then
          npt = npt + 1
          write (16, 1200) iar(i), yy, dyy, atk, datk
        endif
      end do
      dyys = sqrt(abs(yys))*ss
      if (yys.gt.zro) then
	  atks = atks/yys
	  datks= sqrt(abs(datks/yys - atks**2))
      endif
      yys = yys*ss
      write (16, 1500) npt, yys, dyys, atks, datks
      write (16, 1301)
      yys = zro
	atks = zro
	datks= zro
      npt = 0
      do i = 1,ipz+1
        yy = dadz4(3,351,i)*ss
        dyy = sqrt(abs(dadz4(3,351,i)))*ss
        yys = yys + dadz4(3,351,i)
        if (dadz4(3,351,i).gt.zro) then
          atk = dadz4(3,350,i)/dadz4(3,351,i)
	    atks= atks + dadz4(3,350,i)
	    datk = sqrt(abs(dadz4(3,349,i)/dadz4(3,350,i)-atk**2))
		datks= datks + dadz4(3,349,i)  
        else
          atk = zro
	    datk= zro
        endif
        printit = yy.gt.zro 
        if (printit) then
          npt = npt + 1 
          write (16, 1400) izr(i), yy, dyy, atk, datk
        endif
      end do
      dyys = sqrt(abs(yys))*ss
      if (yys.gt.zro) then
	  atks = atks/yys
	  datks= sqrt(abs(datks/yys - atks**2))
      endif
      yys = yys*ss
      write (16, 1600) npt, yys, dyys, atks, datks
c
c   Print out backward isotope formation cross sections:
c
      do iz = 1,izmax + 1
        z(iz) = max (-one, zbeg - dble(iz-1))
        dasum(iz) = zro
        do ia = 1,iamax
          ar(ia) = max(abeg - dble(ia-1), zro)
          dasum(iz) = dasum(iz) + dadz4(4,ia,iz)
        end do
        prtz(iz) = z(iz).ge.zro .and. dasum(iz).gt.zro
      end do
      write (16, 3702)
      do lpr = 1,lprm
        iz1 = 3*lpr - 2
        iz2 = 3*lpr - 1
        iz3 = 3*lpr
        prtz1 = prtz(iz3)
        prtz2 = .not.prtz(iz3) .and. prtz(iz2)
        prtz3 = (.not.prtz(iz2) .and. .not.prtz(iz3)) .and. prtz(iz1)
        if (prtz1) then
          write (16, 3800) z(iz1), z(iz2), z(iz3)
        elseif (prtz2) then
          write (16, 3900) z(iz1), z(iz2)
        elseif (prtz3) then
          write (16, 4000) z(iz1)
        endif
        ppt = zro
        do k = 1,iamax
          y1 = dadz4(4,k,iz1)*faca
          dy1 = sqrt(abs(dadz4(4,k,iz1)))*faca
          y2 = dadz4(4,k,iz2)*faca
          dy2 = sqrt(abs(dadz4(4,k,iz2)))*faca
          y3 = dadz4(4,k,iz3)*faca
          dy3 = sqrt(abs(dadz4(4,k,iz3)))*faca
          printit = y1.gt.zro .or. y2.gt.zro .or. y3.gt.zro
          if (printit) then
c  kkg 11/15/04
            if (ar(k).gt.zro .and. ar(k).ge.z(iz3)) then
              iaa = nint(ar(k))
              prtz1 = prtz(iz3)
              prtz2 = .not.prtz(iz3) .and. prtz(iz2)
              prtz3 = (.not.prtz(iz2) .and. .not.prtz(iz3)) .and. 
     &            		  prtz(iz1)
              if (prtz1) then
                ppt = ppt + one
                write (16, 4100) iaa, y1, dy1, y2, dy2, y3, dy3
              elseif (prtz2) then
                ppt = ppt + one
                write (16, 4200) iaa, y1, dy1, y2, dy2
              elseif (prtz3) then
                ppt = ppt + one
                write (16, 4300) iaa, y1, dy1
              endif
            endif
          endif
        end do
        if (ppt.gt.zro) then
          ipt = nint(ppt)
          ys1  = dasum(iz1)*faca
          dys1 = sqrt(abs(dasum(iz1)))*faca
          ys2  = dasum(iz2)*faca
          dys2 = sqrt(abs(dasum(iz2)))*faca
          ys3  = dasum(iz3)*faca
          dys3 = sqrt(abs(dasum(iz3)))*faca
          if (ys3.gt.zro) then
            write (16, 4400) ipt, ys1, dys1, ys2, dys2, ys3, dys3
          elseif (ys2.gt.zro) then
            write (16, 4500) ipt, ys1, dys1, ys2, dys2
          elseif (ys1.gt.zro) then
            write (16, 4600) ipt, ys1, dys1
          endif 
        endif 
      end do
      write (16, 4702)
      ss = faca
      write (16, 1102)
	atks = zro
	datks= zro
      yys = zro
      npt = 0
      do i = 1,ip
        yy = dadz4(4,i,150)*ss
        dyy = sqrt(abs(dadz4(4,i,150)))*ss
        yys = yys + dadz4(4,i,150)
        if (dadz4(4,i,150).gt.zro) then
          atk = dadz4(4,i,151)/dadz4(4,i,150)
	    atks= atks + dadz4(4,i,151)
	    datk = sqrt(abs(dadz4(4,i,149)/dadz4(4,i,150)-atk**2))
		datks= datks + dadz4(4,i,149)  
        else
          atk = zro
	    datk= zro
        endif
        printit = yy.gt.zro 
        if (printit) then
          npt = npt + 1
          write (16, 1200) iar(i), yy, dyy, atk, datk
        endif
      end do
      dyys = sqrt(abs(yys))*ss
      if (yys.gt.zro) then
	  atks = atks/yys
	  datks= sqrt(abs(datks/yys - atks**2))
      endif
      yys = yys*ss
      write (16, 1500) npt, yys, dyys, atks, datks
      write (16, 1302)
      yys = zro
	atks = zro
	datks= zro
      npt = 0
      do i = 1,ipz+1
        yy = dadz4(4,351,i)*ss
        dyy = sqrt(abs(dadz4(4,351,i)))*ss
        yys = yys + dadz4(4,351,i)
        if (dadz4(4,351,i).gt.zro) then
          atk = dadz4(4,350,i)/dadz4(4,351,i)
	    atks= atks + dadz4(4,350,i)
	    datk = sqrt(abs(dadz4(4,349,i)/dadz4(4,351,i)-atk**2))
		datks= datks + dadz4(4,349,i)  
        else
          atk = zro
	    datk= zro
        endif
        printit = yy.gt.zro 
        if (printit) then
          npt = npt + 1 
          write (16, 1400) izr(i), yy, dyy, atk, datk
        endif
      end do
      dyys = sqrt(abs(yys))*ss
      if (yys.gt.zro) then
	  atks = atks/yys
	  datks= sqrt(abs(datks/yys - atks**2))
      endif
      yys = yys*ss
      write (16, 1600) npt, yys, dyys, atks, datks
c
c   Print out isotope average kinetic energies:
c
      do iz = 1,izmax + 1
        z(iz) = max (-one, zbeg - dble(iz-1))
        dasum(iz) = zro
        do ia = 1,iamax
          ar(ia) = max(abeg - dble(ia-1), zro)
          dasum(iz) = dasum(iz) + dadz(ia,iz)
        end do
        prtz(iz) = z(iz).ge.zro .and. dasum(iz).gt.zro
      end do
      write (16, 3703)
      do lpr = 1,lprm
        iz1 = 3*lpr - 2
        iz2 = 3*lpr - 1
        iz3 = 3*lpr
        prtz1 = prtz(iz3)
        prtz2 = .not.prtz(iz3) .and. prtz(iz2)
        prtz3 = (.not.prtz(iz2) .and. .not.prtz(iz3)) .and. prtz(iz1)
        if (prtz1) then
          write (16, 3800) z(iz1), z(iz2), z(iz3)
        elseif (prtz2) then
          write (16, 3900) z(iz1), z(iz2)
        elseif (prtz3) then
          write (16, 4000) z(iz1)
        endif
        ppt = zro
	  ys1 = zro
	  dys1=zro
	  ys2 = zro
	  dys2=zro
	  ys3 = zro
	  dys3=zro
        do k = 1,iamax
	    if(dadz(k,iz1).gt.zro)  then
		  y1 = dadz4(1,k,iz1)/dadz(k,iz1)
		  dy1 = sqrt(abs(dadz4(2,k,iz1)/dadz(k,iz1)-y1**2))
	      ys1 = ys1 + dadz4(1,k,iz1)
	      dys1= dys1+ dadz4(2,k,iz1)
		else   
	      y1 = zro
	      dy1 = zro
          endif
	    if(dadz(k,iz2).gt.zro)  then
		  y2 = dadz4(1,k,iz2)/dadz(k,iz2)
		  dy2 = sqrt(abs(dadz4(2,k,iz2)/dadz(k,iz2)-y2**2))
	      ys2 = ys2 + dadz4(1,k,iz2)
	      dys2= dys2+ dadz4(2,k,iz2)
		else   
	      y2 = zro
	      dy2 = zro
          endif
	    if(dadz(k,iz3).gt.zro)  then
		  y3 = dadz4(1,k,iz3)/dadz(k,iz3)
		  dy3 = sqrt(abs(dadz4(2,k,iz3)/dadz(k,iz3)-y3**2))
	      ys3 = ys3 + dadz4(1,k,iz3)
	      dys3= dys3+ dadz4(2,k,iz3)
		else   
	      y3 = zro
	      dy3 = zro
          endif
          printit = y1.gt.zro .or. y2.gt.zro .or. y3.gt.zro
          if (printit) then
            if (ar(k).gt.zro .and. ar(k).ge.z(iz3)) then
              iaa = nint(ar(k))
              prtz1 = prtz(iz3)
              prtz2 = .not.prtz(iz3) .and. prtz(iz2)
              prtz3 = (.not.prtz(iz2) .and. .not.prtz(iz3)) .and. 
     &            		  prtz(iz1)
              if (prtz1) then
                ppt = ppt + one
                write (16, 4100) iaa, y1, dy1, y2, dy2, y3, dy3
              elseif (prtz2) then
                ppt = ppt + one
                write (16, 4200) iaa, y1, dy1, y2, dy2
              elseif (prtz3) then
                ppt = ppt + one
                write (16, 4300) iaa, y1, dy1
              endif
            endif
          endif
        end do
        if (ppt.gt.zro) then
          ipt = nint(ppt)
          if(dasum(iz1).gt.zro)  then
	      ys1 = ys1/dasum(iz1)
	      dys1 = sqrt(abs(dys1/dasum(iz1)-ys1**2))
          else
	      ys1 = zro
	      dys1= zro
          endif
          if(dasum(iz2).gt.zro)  then
	      ys2 = ys2/dasum(iz2)
	      dys2 = sqrt(abs(dys2/dasum(iz2)-ys2**2))
          else
	      ys2 = zro
	      dys2= zro
          endif
          if(dasum(iz3).gt.zro)  then
	      ys3 = ys3/dasum(iz3)
	      dys3 = sqrt(abs(dys3/dasum(iz3)-ys3**2))
          else
	      ys3 = zro
	      dys3= zro
          endif
          if (ys3.gt.zro) then
            write (16, 4400) ipt, ys1, dys1, ys2, dys2, ys3, dys3
          elseif (ys2.gt.zro) then
            write (16, 4500) ipt, ys1, dys1, ys2, dys2
          elseif (ys1.gt.zro) then
            write (16, 4600) ipt, ys1, dys1
          endif 
        endif 
      end do
      write (16, 4703)
      ss = faca
      write (16, 1103)
      yys = zro
	tets= zro
	dtets=zro
	vzs  =zro
	dvzs =zro 
	yyfs =zro
	yybs =zro
      npt = 0
      do i = 1,ip
        yy = dadz(i,150)*ss
        dyy = sqrt(abs(dadz(i,150)))*ss
        yys = yys + dadz(i,150)
	  yyfs= yyfs+ dadz4(3,i,150)
	  yybs= yybs+ dadz4(4,i,150) 
 	  if(dadz(i,150).gt.zro)  then
	    atet = dadz4(1,i,151)/dadz(i,150)
	    datet= sqrt(abs(dadz4(2,i,151)/dadz(i,150)-atet**2))
	    tets = tets + dadz4(1,i,151)
	    dtets= dtets+ dadz4(2,i,151) 
	    avz  = dadz4(1,i,150)/dadz(i,150)
	    davz = sqrt(abs(dadz4(2,i,150)/dadz(i,150)- avz**2))
	    vzs = vzs + dadz4(1,i,150)
	    dvzs= dvzs+ dadz4(2,i,150)
	    yyf = dadz4(3,i,150)
	    yyb = dadz4(4,i,150)
          dyyf= sqrt(abs(dadz4(3,i,150))) 
          dyyb= sqrt(abs(dadz4(4,i,150)))
c   calculation of forward/backward ratio
		if(yyb.gt.zro)  then
		  fb = yyf/yyb
		  dfb=(dyyf*yyb + dyyb*yyf)/yyb**2
		else
		  fb = one
		  dfb= zro
		endif   
c		      
        else
	    atet = zro
	    datet= zro
	    avz  = zro
	    davz = zro 	 
		fb = one
		dfb= zro
        endif
        printit = yy.gt.zro 
        if (printit) then
          npt = npt + 1
          write (16, 1201) iar(i), yy, dyy, atet, datet, avz, davz,
     &                     fb,dfb 
        endif
      end do
      dyys = sqrt(abs(yys))*ss
      if (yys.gt.zro) then
	  tets = tets/yys
	  dtets= sqrt(abs(dtets/yys - tets**2))
	  vzs  = vzs/yys
	  dvzs = sqrt(abs(dvzs/yys - vzs**2)) 
	endif    
	if(yybs.gt.zro)  then
	  fbs=yyfs/yybs
	  dyyfs=sqrt(abs(yyfs))
	  dyybs=sqrt(abs(yybs))
	  dfbs=(dyyfs*yybs+dyybs*yyfs)/yybs**2
	else
	   fbs = one
	   dfbs= zro
      endif
      yys = yys*ss
      write (16, 1501) npt, yys, dyys, tets, dtets, vzs, dvzs,
     &                 fbs,dfbs
c
      write (16, 1303)
      yys = zro
	tets= zro
	dtets=zro
	vzs  =zro
	dvzs =zro 
	yyfs =zro
	yybs =zro
      npt = 0
      do i = 1,ipz+1
        yy = dadz(351,i)*ss
        dyy = sqrt(abs(dadz(351,i)))*ss
        yys = yys + dadz(351,i)
	  yyfs= yyfs+ dadz4(3,351,i)
	  yybs= yybs+ dadz4(4,351,i) 
	  if(dadz(351,i).gt.zro)  then
	    atet = dadz4(1,351,i)/dadz(351,i)
	    datet= sqrt(abs(dadz4(2,351,i)/dadz(351,i)-atet**2))
	    tets = tets + dadz4(1,351,i)
	    dtets= dtets+ dadz4(2,351,i) 
	    avz  = dadz4(1,350,i)/dadz(351,i)
	    davz = sqrt(abs(dadz4(2,350,i)/dadz(351,i)- avz**2))
	    vzs = vzs + dadz4(1,350,i)
	    dvzs= dvzs+ dadz4(2,350,i)
	    yyf = dadz4(3,350,i)
	    yyb = dadz4(4,350,i)
          dyyf= sqrt(abs(dadz4(3,350,i))) 
          dyyb= sqrt(abs(dadz4(4,350,i)))
c   calculation of forward/backward ratio
		if(yyb.gt.zro)  then
		  fb = yyf/yyb
		  dfb=(dyyf*yyb + dyyb*yyf)/yyb**2
		else
		  fb = one
		  dfb= zro
		endif   
        else
	    atet = zro
	    datet= zro
	    avz  = zro
	    davz = zro 	 
		fb = one
		dfb= zro
        endif
        printit = yy.gt.zro 
        if (printit) then
          npt = npt + 1 
          write (16, 1401) izr(i), yy, dyy, atet, datet, avz, davz,
     &                     fb,dfb 
        endif
      end do
      dyys = sqrt(abs(yys))*ss
      if (yys.gt.zro) then
	  tets = tets/yys
	  dtets= sqrt(abs(dtets/yys - tets**2))
	  vzs  = vzs/yys
	  dvzs = sqrt(abs(dvzs/yys - vzs**2)) 
	endif    
	if(yybs.gt.zro)  then
	  fbs=yyfs/yybs
	  dyyfs=sqrt(abs(yyfs))
	  dyybs=sqrt(abs(yybs))
	  dfbs=(dyyfs*yybs+dyybs*yyfs)/yybs**2
	else
	   fbs = one
	   dfbs= zro
      endif
      yys = yys*ss
      write (16, 1601) npt, yys, dyys, tets, dtets, vzs, dvzs,
     &                 fbs, dfbs	                
c
      return  

c ======================================================================
 3701 format (/1x,'*************** Nuclide yields [mb] in Forward ',
     &'direction (theta_lab < 90) (zero values suppressed) *******') 
 3702 format (/1x,'*************** Nuclide yields [mb] in Backward ',
     &'direction (theta_lab > 90) (zero values suppressed) *******') 
 3703 format (/1x,'*************** Nuclide average kinetic energies ',
     &'(MeV (zero yield suppressed) *******') 
 3800 format (/17x,'Z = ',f4.0,15x,'Z = ',f4.0,15x,'Z = ',f4.0)
 3900 format (/17x,'Z = ',f4.0,15x,'Z = ',f4.0)
 4000 format (/17x,'Z = ',f4.0)
 4100 format (1x,'A =',i4,3(1x,1pe9.3,' +/- ',e8.2))
 4200 format (1x,'A =',i4,2(1x,1pe9.3,' +/- ',e8.2))
 4300 format (1x,'A =',i4,1(1x,1pe9.3,' +/- ',e8.2))
 4400 format (1x,'S =',i4,3(1x,1pe9.3,' +/- ',e8.2))
 4500 format (1x,'S =',i4,2(1x,1pe9.3,' +/- ',e8.2))
 4600 format (1x,'S =',i4,1(1x,1pe9.3,' +/- ',e8.2))
 4701 format (/1x,'End of nuclide yields (Forward direction).') 
 4702 format (/1x,'End of nuclide yields (Backward direction).') 
 4703 format (/1x,'End of nuclide average kinetic energies.') 
 1101 format (/'Mass Yield[mb] and Mean Kin. Energy[MeV]',
     %         ' of Resid. Nucl. produced in Forward direction:')
 1102 format (/'Mass Yield[mb] and Mean Kin. Energy[MeV]',
     %         ' of Resid. Nucl. produced in Backward direction:')
 1103 format (/'Mass Yield[mb], Mean emiss. angle[deg.]',
     %         ' and Mean z-velocity[1/c] of Resid. Nucl.',
     &           ' and Forward/Backward ratio :')
 1200 format (1x,'A =',i4,1x,1pe9.3,'+/-',e8.2,5x,e9.3,'+/-',e8.2)
 1201 format (1x,'A =',i4,1x,1pe9.3,'+/-',e8.2,4(2x,e10.3,'+/-',e9.2))
 1500 format (1x,'S =',i4,1x,1pe9.3,'+/-',e8.2,5x,e9.3,'+/-',e8.2)
 1501 format (1x,'S =',i4,1x,1pe9.3,'+/-',e8.2,4(2x,e10.3,'+/-',e9.2))
 1301 format (/'Charge Yield[mb] and Mean Kin. Energy[MeV]',
     &         ' of Resid. Nucl. produced in Forward direction:')
 1302 format (/'Charge Yield[mb] and Mean Kin. Energy[MeV]',
     &         ' of Resid. Nucl. produced in Backward direction:')
 1303 format (/'Charge Yield[mb], Mean emiss. angle[deg.]',
     %         ' and Mean z-velocity[1/c] of Resid. Nucl.',
     &           ' and Forward/Backward ratio :')
 1400 format (1x,'Z =',i3,1x,1pe9.3,'+/-',e8.2,5x,e9.3,'+/-',e8.2)
 1401 format (1x,'Z =',i3,1x,1pe9.3,'+/-',e8.2,4(2x,e10.3,'+/-',e9.2))
 1600 format (1x,'S =',i3,1x,1pe9.3,'+/-',e8.2,5x,e9.3,'+/-',e8.2)
 1601 format (1x,'S =',i3,1x,1pe9.3,'+/-',e8.2,4(2x,e10.3,'+/-',e9.2))

c ======================================================================
      end
    
C  *************************************************************
      

    subroutine prrdisq(anucl,znucl,mb0,fn)

c ======================================================================
c
c   This subroutine prints out the information of residual nuclei after
c   cascade, preequilibrium, evaporation and fission stages of reaction 
c
c   kkg 12/02/04
      implicit real*8 (a-h, o-z), integer (i-n)
	common /rezdis/  rdis(5,5,250), dex, dpm
	dimension pro(5), iar(250), izr(150)
	data zro /0.0d0/
c
      abeg = anucl + 5.d0
      zbeg = znucl + 5.d0
      izmax = nint(zbeg)
      izmax = min(izmax, 247)
      iamax = nint(abeg) + mb0
      iamax = min(iamax, 247)
      ip = nint(abeg)  
      ipz = nint(zbeg) 
      do i = 1,ip
        iar(i) = ip - i + 1
      end do
      do i= 1,ipz
        izr(i) = ipz - i + 1
      end do
c   Print out the A-distribution of residual nuclei
      write(16,100)  
  100 format(/24x,' Mass distributions of residual nuclei'/
     & 12x,'After cascade',' After preco ','Dec. by evap ',
     & 'Dec. by fiss.',' Before fission'/) 
      do  i=1,ip
	   s = zro  	     
	   if(i.le.247)  then
	     do  j=1,5
	       pro(j) = rdis(j,1,i)/fn
		   s = s + pro(j) 
           enddo
		 if(s.gt.1.0d-10) write(16,101) iar(i),pro
  101      format(1x,'  A =',i4,5(4x,1pe9.3))		  
         endif
      enddo
	do j=1,5
	  pro(j) = zro 
	  if(rdis(j,1,250).gt.zro) pro(j) = rdis(j,1,248)/rdis(j,1,250)
      enddo
	write(16,102) pro
  102 format(4x,'<A> =',1x,5(4x,1pe9.3))
	do j=1,5
	  pro(j) = zro 
	  if(rdis(j,1,250).gt.zro) then 
	    av  = rdis(j,1,248)/rdis(j,1,250)
	    pro(j) = sqrt(abs(rdis(j,1,249)/rdis(j,1,250) - av**2))
        endif
      enddo
	write(16,103) pro
  103 format(4x,'DisA=',1x,5(4x,1pe9.3))
	do j=1,5
	  pro(j) = zro 
	  pro(j) = rdis(j,1,250)/fn
      enddo
	write(16,104) pro
  104 format(4x,'norm=',1x,5(4x,1pe9.3))
c
c   Print out the Z-distribution of residual nuclei
      write(16,200)  
  200 format(/24x,' Charge distributions of residual nuclei'/
     & 12x,'After cascade',' After preco ','Dec. by evap ',
     & 'Dec. by fiss.',' Before fission'/) 
      do  i=1,ipz
	   s = zro  	     
	   if(i.le.247)  then
	     do  j=1,5
	       pro(j) = rdis(j,2,i)/fn
		   s = s + pro(j) 
           enddo
		 if(s.gt.1.0d-10) write(16,201) izr(i),pro
  201      format(1x,'  Z =',i4,5(4x,1pe9.3))		  
         endif
      enddo
	do j=1,5
	  pro(j) = zro 
	  if(rdis(j,2,250).gt.zro) pro(j) = rdis(j,2,248)/rdis(j,2,250)
      enddo
	write(16,202) pro
  202 format(4x,'<Z> =',1x,5(4x,1pe9.3))
	do j=1,5
	  pro(j) = zro 
	  if(rdis(j,2,250).gt.zro) then 
	    av  = rdis(j,2,248)/rdis(j,2,250)
	    pro(j) = sqrt(abs(rdis(j,2,249)/rdis(j,2,250) - av**2))
        endif
      enddo
	write(16,203) pro
  203 format(4x,'DisZ=',1x,5(4x,1pe9.3))
	do j=1,5
	  pro(j) = zro 
	  pro(j) = rdis(j,2,250)/fn
      enddo
	write(16,204) pro
  204 format(4x,'norm=',1x,5(4x,1pe9.3))
c
c   Print out the E*-distribution of residual nuclei
      write(16,300)  
  300 format(/20x,' Excitation energy distributions [1/MeV]',
     & ' of residual nuclei'/
     & 6x,'E*(MeV)',3x,'After cascade',' After preco ','Dec. by evap ',
     & 'Dec. by fiss.',' Before fission'/) 
      do  i=1,247
	   ex1 = (i-1)*dex 
	   ex2 = ex1 + dex 
	   s = zro  	     
	   do  j=1,5
	     pro(j) = rdis(j,3,i)/fn/dex
		 s = s + pro(j) 
         enddo
	   if(s.gt.1.0d-10) write(16,301) ex1,ex2,pro
  301      format(1x,f6.0,'-',f6.0,5(4x,1pe9.3))		  
      enddo
	do j=1,5
	  pro(j) = zro 
	  if(rdis(j,3,250).gt.zro) pro(j) = rdis(j,3,248)/rdis(j,3,250)
      enddo
	write(16,302) pro
  302 format(7x,'<E*> =',1x,5(4x,1pe9.3))
	do j=1,5
	  pro(j) = zro 
	  if(rdis(j,3,250).gt.zro) then 
	    av  = rdis(j,3,248)/rdis(j,3,250)
	    pro(j) = sqrt(abs(rdis(j,3,249)/rdis(j,3,250) - av**2))
        endif
      enddo
	write(16,303) pro
  303 format(7x,'DisE*=',1x,5(4x,1pe9.3))
	do j=1,5
	  pro(j) = zro 
	  pro(j) = rdis(j,3,250)/fn
      enddo
	write(16,304) pro
  304 format(7x,' norm=',1x,5(4x,1pe9.3))
c
c   Print out the P-distribution of residual nuclei
c  KKG  05/23/07
      write(16,400)  
  400 format(/20x,' Linear momentum per nucleon distributions [1/MeV/c]'
     & ,' of residual nuclei'/
     &3x,'P/A(MeV/c)',3x,'After cascade',' After preco ','Dec. by evap '
     & ,'Dec. by fiss.',' Before fission'/) 
      do  i=1,247
	   pm1 = (i-1)*dpm 
	   pm2 = pm1 + dpm 
	   s = zro  	     
	   do  j=1,5
	     pro(j) = rdis(j,4,i)/fn/dpm
		 s = s + pro(j) 
         enddo
	   if(s.gt.1.0d-10) write(16,401) pm1,pm2,pro
  401      format(1x,f6.0,'-',f6.0,5(4x,1pe9.3))		  
      enddo
	do j=1,5
	  pro(j) = zro 
	  if(rdis(j,4,250).gt.zro) pro(j) = rdis(j,4,248)/rdis(j,4,250)
      enddo
	write(16,402) pro
  402 format(6x,'<P/A> =',1x,5(4x,1pe9.3))
	do j=1,5
	  pro(j) = zro 
	  if(rdis(j,4,250).gt.zro) then 
	    av  = rdis(j,4,248)/rdis(j,4,250)
	    pro(j) = sqrt(abs(rdis(j,4,249)/rdis(j,4,250) - av**2))
        endif
      enddo
	write(16,403) pro
  403 format(4x,'Dis(P/A)=',1x,5(4x,1pe9.3))
	do j=1,5
	  pro(j) = zro 
	  pro(j) = rdis(j,4,250)/fn
      enddo
	write(16,404) pro
  404 format(7x,' norm=',1x,5(4x,1pe9.3))
c
c   Print out the L-distribution of residual nuclei
      write(16,500)  
  500 format(/20x,' Angular momentum distributions [1/hc]',
     & ' of residual nuclei'/
     & 5x,'  L(hc)',4x,'After cascade',' After preco ','Dec. by evap ',
     & 'Dec. by fiss.',' Before fission'/) 
	dam = 1.0d0 
      do  i=1,247
	   am1 = (i-1)*dam 
	   am2 = am1 + dam 
	   s = zro  	     
	   do  j=1,5
	     pro(j) = rdis(j,5,i)/fn/dam
		 s = s + pro(j) 
         enddo
	   if(s.gt.1.0d-10) write(16,501) am1,am2,pro
  501      format(1x,f6.0,'-',f6.0,5(4x,1pe9.3))		  
      enddo
	do j=1,5
	  pro(j) = zro 
	  if(rdis(j,5,250).gt.zro) pro(j) = rdis(j,5,248)/rdis(j,5,250)
      enddo
	write(16,502) pro
  502 format(7x,' <L> =',1x,5(4x,1pe9.3))
	do j=1,5
	  pro(j) = zro 
	  if(rdis(j,5,250).gt.zro) then 
	    av  = rdis(j,5,248)/rdis(j,5,250)
	    pro(j) = sqrt(abs(rdis(j,5,249)/rdis(j,5,250) - av**2))
        endif
      enddo
	write(16,503) pro
  503 format(7x,' DisL=',1x,5(4x,1pe9.3))
	do j=1,5
	  pro(j) = zro 
	  pro(j) = rdis(j,5,250)/fn
      enddo
	write(16,504) pro
  504 format(7x,' norm=',1x,5(4x,1pe9.3))
c
      return
c ======================================================================
    end
    
c ======================================================================

      subroutine propanq(fn)

c ======================================================================
c
c   This subroutine prints out the distribution of fission fragments
c   opening angle (lab.sys.) in different bin of the produced neutron
c   multiplicity
c
c   kkg 12/02/04
      implicit real*8 (a-h, o-z), integer (i-n)
	common /fisopa/  opan(7,185),dth12
	dimension pro(7)
	data zro /0.0d0/
c 
      write(16,100)
  100 format(/6x,'Distribution of angle between fission',
     &' fragments[1/deg.]',' (lab.sys.) in different'/
     &30x,' bin of neutron',' multiplicity:'/1x,'theta(deg.)',
     &1x,' All events','  n=0 - 5  ','  n=6 - 8  ','  n=9 - 12 ',
     &' n=13 - 15 ',' n=16 - 19 ','   n > 20')	              
      do  i=1,180
	  tet1 = (i-1)*dth12
	  tet2 = tet1 + dth12
	  s = zro  
        do j=1,7
      	pro(j) = opan(j,i)/fn/dth12
		s = s + pro(j)
	  enddo
	  if(s.gt.zro) write(16,101) tet1,tet2,pro 
  101   format(1x,f4.0,' - ',f4.0,7(2x,1pe9.3))
      enddo
	do  j=1,7
	  pro(j) = zro 
	  if(opan(j,185).gt.zro) pro(j) = opan(j,183)/opan(j,185)
      enddo
	write(16,102) pro
  102 format(5x,'<thet>=',7(2x,1pe9.3))
	do  j=1,7
	  pro(j) = zro 
	  if(opan(j,185).gt.zro) then
	    av = opan(j,183)/opan(j,185)
	    pro(j) = sqrt(abs(opan(j,184)/opan(j,185) - av**2))
        endif
      enddo
	write(16,103) pro
  103 format(5x,'Dthet =',7(2x,1pe9.3))
      do  j=1,7
	  pro(j) = opan(j,185)/fn
      enddo
	write(16,104) pro
  104 format(5x,'norm. =',7(2x,1pe9.3)/)
      return
c ======================================================================
      end

c ======================================================================
      
      subroutine pdisnmq(fn)

c ======================================================================
c
c   This subroutine prints out the multiplicity distribution of neutrons 
c   produced :
c              1) after all stages, 
c              2) in cascade stage, 
c              3) in preeq. stage, 
c              4) in evaporation without fission,
c              5) in prefission stage, 
c              6) in postfission stage  
c
c   kkg 12/02/04
      implicit real*8 (a-h, o-z), integer (i-n)
      common /disnmu/ disnm(6,155)
	dimension pro(6)
	data zro /0.0d0/
c 
      write(16,100)
  100 format(/26x,'Multiplicity distribution of neutrons produced:'/
     &2x,'Nn',
     &' After all stages',' At cascade stage',' At preeq. stage ', 
     &'At evap. w/o fis.',' At prefiss stage',' At postfiss. stage')
      do  i=1,152
	  nn = i-1
	  s = zro  
        do j=1,6
      	pro(j) = disnm(j,i)/fn
		s = s + pro(j)
	  enddo
	  if(s.gt.zro) write(16,101) nn,pro 
  101   format(1x,i3,2x,6(4x,1pe9.3,4x))
      enddo
	do  j=1,6
	  pro(j) = zro 
	  if(disnm(j,155).gt.zro) pro(j) = disnm(j,153)/disnm(j,155)
      enddo
	write(16,102) pro
  102 format(1x,'<n>=',1x,6(4x,1pe9.3,4x))
	do  j=1,6
	  pro(j) = zro 
	  if(disnm(j,155).gt.zro) then
	    av = disnm(j,153)/disnm(j,155)
	    pro(j) = sqrt(abs(disnm(j,154)/disnm(j,155) - av**2))
        endif
      enddo
	write(16,103) pro
  103 format(1x,' Dn=',1x,6(4x,1pe9.3,4x))
      do  j=1,6
	  pro(j) = disnm(j,155)/fn
      enddo
	write(16,104) pro
  104 format(1x,'norm=',6(4x,1pe9.3,4x)/)
      return
c ======================================================================
      end