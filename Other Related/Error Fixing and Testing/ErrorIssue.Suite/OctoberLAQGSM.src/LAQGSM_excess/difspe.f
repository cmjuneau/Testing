c  **********************************************************************
c
c   Updated for MARS by K.K. Gudima, March 2007
c   FERMID ==> fermid, rezdist==>rezdisq, opandis==>opandiq.
c   disnmul==> disnmuq,prtdadz==>prdadzq
c   prrdis==> prrdisq, propan==>propanq, pdisnm==>pdisnmq,
c   HELP4 ==> HELPQ, CINEM8==> CINEMQ'
c   FITAFPA==>FITAFPAQ, FITAFAC==>FITAFACQ
      SUBROUTINE DIFSPE(NCAS,INTEL,NCAS1,MV)
      implicit real*8 (a-h, o-z), integer (i-n)
c     06/13/02 (LA version with GEM2 and new Fermi decay)
c   pi,p,n,d,t,He-3,He-4,K+,K- and pbar spectra
c                DIFSPE is called by laqgsmg
C     last corrections 11/27/02,10/23/03, 01/06/05
c last editing by SGM, 03/16/05
ccc, modified by SGM, 06/08/06
C
c*******************************************************************
c
c      JJJ = 1 - only cascade hadrons
c            2 - cascade + evaporation
c            3 - cascade + evaporation + Fermi decay
c      JJJ =11 - cascade + coalescense
c           12 - cascade + coalescense + evaporation
c           13 - cascade + coalescense + evaporation + Fermi decay
c
c  ---- observer's frame is that where real calculations have been done;
c  ---- observables may be calculated in three different frames:
c  ------------- Laboratory System (LAB.S)
c  ------------- Equal-Velicity System (EVS)
c  ------------- Center-of-Mass System (CMS)
c    ISYS =1  - LAB.S is the observer's frame (no boost)
c    ISYS=+2  - to boost into LAB.S if EVS is observer's frame
c    ISYS=-2  - to boost into EVS if LAB.S is observer's frame
c    ISYS=+3  - to boost into LAB.S if CMS is observer's frame
c    ISYS=-3  - to boost into CMS if LAB.S is observer's frame
c
c     ID=0 for coalescence and particles emitted from hot nuclei
c                          IORI(1,M)  IORI(2,M)   IORI(3,M)
c     residual proj.nuclei  0          0             1
c     residual targ.nuclei  0          0             2
c     pre-equilibrium       A-parent   Z-parent      1
c     equilibrium           A-parent   Z-parent      2
c     Fermi-decay           A-parent   Z-parent      3
c     coalescence           0          0             4
c     multifragmentation    A-parent   Z-parent      5
c     fission fragment      A-parent   Z-parent      6
c
c***************************************************************
C
c  kkg 04/07/04
      INTEGER*8 INTEL
      COMMON/RESULT/ANU1,ANU2,ZNU1,ZNU2,ENU1,ENU2,PNU1(3),PNU2(3),
     *AMNU1(3),AMNU2(3)
      COMMON/RESUL4/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      COMMON /EXCIT/ TEXC1,TEXC2,HEXC1,HEXC2,PEXC1,PEXC2
      COMMON/BL1003/U,A,Z
      COMMON/DELEN/WF,Nfiss
      COMMON/FISSION/W1,W2,NF1,NF2
      COMMON/AZT/AT(2),ZT(2)
      COMMON /EXCIT4/ TEX1,TEX2,HEX1,HEX2,PEX1,PEX2
      COMMON/COARAD/P0D,P0T,P0A
      COMMON/DEXIND/JJJ
      common /adbf/   amf, r0m, ijsp, nhump
      common /targ0/ atar0,ztar0,ener0
      COMMON/NCASCA/NCASC,NCPRI
      common /fitaf/fitaf,fitaf1
      common /fitaf2/fitpt(2),fitpt1(2)
      common /gbrems/ tgmin, tgmax, teqv, sxabs, ibrems 
      COMMON /MEMORY/ PME(9,5999),IME(5,5999)
      COMMON/IDPME/ IDPME(5999)
      COMMON/PORIG/ IORI(3,5999)
      dimension UT(2)
c CMJ 08/2016 - used to retrieve kstart in CEM
			COMMON /KST1/ KST1
c
      AN1=ANU1
      AN2=ANU2
      ZN1=ZNU1
      ZN2=ZNU2
      ENEXT1=ENU1
      ENEXT2=ENU2
      do  k=1,3
        PNUCL1(k)=PNU1(k)
        PNUCL2(k)=PNU2(k)
        AMNUC1(k)=AMNU1(k)
        AMNUC2(k)=AMNU2(k)
      enddo
      TEX1=TEXC1
      TEX2=TEXC2
      HEX1=HEXC1
      HEX2=HEXC2
      PEX1=PEXC1
      PEX2=PEXC2
c
c      if(AN2.gt.atar0)   write(*,100) AN2,ZN2,ENEXT2
 100  format(1x,'AN2,ZN2,ENEXT2=',2F8.1,E11.4)
c
      AT(1)=0.
      ZT(1)=0.
      UT(1)=0.
      AT(2)=0.
      ZT(2)=0.
      UT(2)=0.
      CALL  UPPME(MV)
      KST1=MV+1
      IF(JJJ.EQ.1)     GO  TO  14
      if(JJJ.ge.11)    then
        if(NCASC.ge.NCPRI) write( *,*) ' to coales',MV
        CALL  COALES(MV,P0D,P0T,P0A)
        if(NCASC.ge.NCPRI) write( *,*) ' from coales',MV
        KST1=MV+1
        IF(JJJ.EQ.11)  GO  TO  14
      endif
      IF(AN1.LT.0.9)   GO  TO  13
ccc, SGM, 06/08/06      IF(AN1.le.11.0.and.(JJJ.eq.3.or.JJJ.eq.13)) then
      IF(AN1.le.12.0.and.(JJJ.eq.3.or.JJJ.eq.13)) then
        KST1=KST1-1
        if(NCASC.ge.NCPRI) write( *,*) ' to fermid1',AN1,ZN1,ENEXT1
        CALL  fermid(AN1,ZN1,ENEXT1,PNUCL1,KST1)
        if(NCASC.ge.NCPRI) write( *,*) ' from fermid1'
        KST1=KST1+1
        GO  TO  13
      else
        GO  TO  113
      endif
  113 CONTINUE
      WF=0.
	Nfiss=0
      wam1=1.
      NEX1=INT(TEX1)
      NHEX1=INT(HEX1)
      NZEX1=INT(PEX1)
      NPEX1=NEX1-NHEX1
      elx1=AMNUC1(1)
      ely1=AMNUC1(2)
      elz1=AMNUC1(3)
      ENEX1M = ENEXT1*1000.0d0
      if(NCASC.ge.NCPRI) write( *,*) ' to precof',AN1,ZN1,ENEXT1
c    
      fitaf = fitpt(1)
      fitaf1= fitpt1(1) 
c
      CALL  precof(ENEX1M,AN1,ZN1,PNUCL1(1),PNUCL1(2),PNUCL1(3),
     *              elx1,ely1,elz1,r0m,KST1,
     *              NEX1,NPEX1,NHEX1,NZEX1,wam1)
      if(NCASC.ge.NCPRI) write( *,*) ' from precof'
      AT(1)=A
      ZT(1)=Z
c   kkg  12/14/05
c     include in MEMORY projectile residual nucleus
      if(AT(1).gt.0.9)  then !SGM, 12/30/09, restored back per KKG 03.01.10
c      if(AT(1).gt.0.9.and.NINT(AT(1)).ne.IME(4,KST1-1))  then !SGM, 12/30/09, per Marcus
       PME(1,KST1) = 0.0        ! x    
       PME(2,KST1) = 0.0        ! y    
       PME(3,KST1) = 0.0        ! z    
       PME(4,KST1) = PNUCL1(1)   ! px    
       PME(5,KST1) = PNUCL1(2)   ! py    
       PME(6,KST1) = PNUCL1(3)   ! pz
c   kkg 06/27/06
       PME(7,KST1) = UT(1)     ! excit. energy of proj. residue (MeV)
       amass = 0.940*AT(1)   
       PME(8,KST1) = 
     & SQRT(PNUCL1(1)**2+PNUCL1(2)**2+PNUCL1(3)**2+amass**2)-amass
       PME(9,KST1) = amass
       IME(1,KST1) = NINT(ZT(1))
       IME(2,KST1) = 0
       IME(3,KST1) = 0
       IME(4,KST1) = NINT(AT(1))
       IME(5,KST1) = 0
       IDPME(KST1) = 0
       IORI(1,KST1)= 0
       IORI(2,KST1)= 0
       IORI(3,KST1)= 1
       KST1 = KST1 + 1
      endif 
      W1=W1+WF
	NF1=NF1+Nfiss
   13 CONTINUE
      MV1=KST1-1
      IF(AN2.LT.0.9)   GO  TO  14
ccc, SGM, 06/08/06      IF(AN2.le.11.0.and.(JJJ.eq.3.or.JJJ.eq.13)) then
      IF(AN2.le.12.0.and.(JJJ.eq.3.or.JJJ.eq.13)) then
        KST1=KST1-1
        if(NCASC.ge.NCPRI) write( *,*) ' to fermid2',AN2,ZN2,ENEXT2
        CALL  fermid(AN2,ZN2,ENEXT2,PNUCL2,KST1)
        if(NCASC.ge.NCPRI) write( *,*) ' from fermid2'

        KST1=KST1+1
        GO  TO  14
      else
        GO  TO  114
      endif
  114 CONTINUE
      WF=0.
	Nfiss=0
      wam2=1.
      NEX2=INT(TEX2)
      NHEX2=INT(HEX2)
      NZEX2=INT(PEX2)
      NPEX2=NEX2-NHEX2
      elx2=AMNUC2(1)
      ely2=AMNUC2(2)
      elz2=AMNUC2(3)
      ENEX2M = ENEXT2*1000.0d0
      if(NCASC.ge.NCPRI) write( *,*) ' to precof2',AN2,ZN2,ENEXT2
      if(ibrems.eq.1)  then
c   Determine semiempirical af/an fits for fission cross sections:
c   kkg 12/12/05
        if(ztar0.ge.67.0.and.ztar0.le.88.d0)  then 
          call FITAFPAQ (atar0,ztar0,ener0,fitaf, fitaf1)
        elseif(ztar0.gt.88.0)  then
          call FITAFACQ (atar0,ztar0,ener0,fitaf, fitaf1)
        else
          fitaf = 1.0d0
          fitaf1= 1.0d0
	      endif  
      else 
c    
      fitaf = fitpt(2)
      fitaf1= fitpt1(2) 
c
c                         fitaf = a_f(CEM)/a_f(RAL)
c                         fitaf1 = C(Z)[CEM]/C(Z)[RAL]
c 
      endif
      CALL  precof(ENEX2M,AN2,ZN2,PNUCL2(1),PNUCL2(2),PNUCL2(3),
     *              elx2,ely2,elz2,r0m,KST1,
     *              NEX2,NPEX2,NHEX2,NZEX2,wam2)
      if(NCASC.ge.NCPRI) write( *,*) ' from precof2'
      AT(2)=A
      ZT(2)=Z
c   kkg  12/14/05
c     include in the MEMORY target residual nucleus
      if(AT(2).gt.0.9)  then !SGM, 12/31/09, restored back per KKG 03.01.10
c      if(AT(2).gt.0.9.and.NINT(AT(2)).ne.IME(4,KST1-1)) then !SGM, 12/31/09, per Marcus
       PME(1,KST1) = 0.0        ! x    
       PME(2,KST1) = 0.0        ! y    
       PME(3,KST1) = 0.0        ! z    
       PME(4,KST1) = PNUCL2(1)   ! px    
       PME(5,KST1) = PNUCL2(2)   ! py    
       PME(6,KST1) = PNUCL2(3)   ! pz
c   kkg 06/27/06
       PME(7,KST1) = UT(2)      ! excit. energy of targ. residue (MeV)
       amass = 0.940*AT(2)   
       PME(8,KST1) = 
     & SQRT(PNUCL2(1)**2+PNUCL2(2)**2+PNUCL2(3)**2+amass**2)-amass
       PME(9,KST1) = amass
       IME(1,KST1) = NINT(ZT(2))
       IME(2,KST1) = 0
       IME(3,KST1) = 0
       IME(4,KST1) = NINT(AT(2))
       IME(5,KST1) = 0
       IDPME(KST1) = 0
       IORI(1,KST1)= 0
       IORI(2,KST1)= 0
       IORI(3,KST1)= 2
       KST1 = KST1 + 1
      endif 
      W2=W2+WF
      NF2=NF2+Nfiss
   14 continue
c  kkg  08/02/07
      call  corrdec(MV,KST1)
      MV2=KST1-1       
c     
      CALL  PMEUP(MV2)
      CALL  SPECTR(MV,MV1,MV2,NCAS1)
c  kkg  06/27/06
      MV=MV2
c  after this point all products (MV) are in PME and IME !!!
c  for products with mass(baryon) number IB=IME(4,M)>4, exitation energy
c  (MeV) is in PME(7,M). 
      RETURN
      END
      
c************************************************
C
      SUBROUTINE  corrdec(MV,KST1)
c  decays instable framents, KKG 08/02/07 
c  editted by SGM, 08/18/07; 06/11/08   
      implicit real*8 (a-h, o-z), integer (i-n)
      COMMON /MEMORY/ PME(9,5999),IME(5,5999)
      COMMON /IDPME/ IDPME(5999)
      COMMON/PORIG/ IORI(3,5999)
      dimension pn(3)
c  list of 'stable' nuclei with A<13
c Z=0               n
c Z=1               p     d      t     
c Z=2             He3   He4    He6    He8
c Z=3             Li6   Li7    Li8    Li9
c Z=3             Li6   Li7    Li8    Li9 Li11             !SGM, 06/11/08
c Z=4             Be7   Be9    Be10   Be11  Be12 ...   
c Z=5             B8    B10    B11    B12 ...    
cc Z=6             C10   C11    C12    
c Z=6             C9 C10   C11    C12 ...                  !SGM, 06/11/08   
c Z=7             N12    
c                   
c
      m1=MV+1
      m2=KST1-1
      ks=m2
      if(m2.le.m1)    return
      do  m=m1,m2
        ia = IME(4,m)
        iz = IME(1,m)
        if(ia.ge.2.and.ia.le.12)     then
 	   if((iz.eq.0.and.ia.gt.1).or.
     &      (iz.eq.1.and.ia.gt.3).or. 		 
     &      (iz.eq.2.and.(ia.lt.3.or.ia.eq.5.or.ia.eq.7.or.ia.gt.8)).or. 		 
c     &      (iz.eq.3.and.(ia.lt.6.or.ia.gt.9)).or. 
     &      (iz.eq.3.and.(ia.lt.6.or.ia.eq.10.or.ia.gt.11)).or. !SGM, 06/11/08 		 
     &      (iz.eq.4.and.(ia.lt.7.or.ia.eq.8.or.ia.gt.12)).or. 		 
     &      (iz.eq.5.and.(ia.lt.8.or.ia.eq.9.or.ia.gt.12)).or. 
c     &      (iz.eq.6.and.(ia.lt.10.or.ia.gt.12)).or.
     &      (iz.eq.6.and.(ia.lt.9.or.ia.gt.12)).or.            !GSM, 06/11/08
     &      (iz.eq.7.and.(ia.lt.12.or.ia.gt.12)))    then		 
             a = dble(ia)
             z = dble(iz)
             exg = PME(7,m)/1000.0d0
             if(exg.lt.0.001d0)  exg=0.001d0
             pn(1) = PME(4,m)
             pn(2) = PME(5,m)
             pn(3) = PME(6,m)
             call  fermid(a,z,exg,pn,ks,1)
c            
ccc, SGM            write(*,*) ' decay of a,z,m,m2,ks=',a,z,m,m2,ks
c
             do  i=1,9
               PME(i,m)=PME(i,ks)
               if(i.le.5)  IME(i,m)=IME(i,ks)
             enddo
             IDPME(m)=IDPME(ks)
             IORI(1,m)=IORI(1,ks)    
             IORI(2,m)=IORI(2,ks)    
             IORI(3,m)=IORI(3,ks)
             ks=ks-1    
          endif
        endif          
      enddo
      KST1=ks+1   
      return
      end           
C
C  *************************************************************
C
      SUBROUTINE  SPECTR(MV,MV1,MV2,NCAS1)
      implicit real*8 (a-h, o-z), integer (i-n)
c      REAL*4 UP1,UP2
      logical fisevent
c kkg
      common /yeldaz/ dadz(351,151)
      common /yeldaz1/ dadz4(4,351,151)
	common /pmultip/ pmult(5,9)
      common /targ0/ atar0,ztar0,ener0
c  kkg  10/27/05
c     common /tarres/ teres
      common /trec/ teres
c kkg
      COMMON/INSP/INSP
      COMMON /UPAC/UP1(66000),IPER
      COMMON /UPAC2/UP2(70000),LU2
      COMMON/SPE/SPE(11,11,205)
      COMMON/DOMEGA/DOM(10),UF(10),DXS(10),DU,XMIN
      COMMON/MULT/AMUL(11),MULT(11)
      COMMON/DELEN/WF,Nfiss /AZT/AT(2),ZT(2)
      COMMON/VGSYS/T0A,VLA,GLA,VEV,GEV,VCM,GCM,ISYS
      COMMON/IDPME/ IDPME(5999)
      COMMON/PORIG/ IORI(3,5999)
      COMMON/TGRLE/ T0le,tgr(4),dtgr(4),jgr(4)
      common /angwin/ ewin(2,11),angw(25,11),dthw
      COMMON/RESUL4/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
c  kkg  12/14/05
      common /indgsi/ IGSI
c
      IF(MV2.EQ.0)      RETURN
         nn   = 0
         nnc  = 0
         nnp  = 0
         nne  = 0
         nnprf= 0
         nnpof= 0 
c
         abeg=atar0+5.
         zbeg=ztar0+5.
c
      MSYS=IABS(ISYS)
      GO  TO  (1,2,3),MSYS
    1 VSYS=0.
      GSYS=1.
      GO  TO  4
    2 VSYS=ISIGN(1,ISYS)*VEV
      GSYS=GEV
      GO  TO  4
    3 VSYS=ISIGN(1,ISYS)*VCM
      GSYS=GCM
    4 CONTINUE
      NCAS1=NCAS1+1
C
      DO  12  M=1,MV2
      LU=11*M
      CM=UP1(LU- 7)
      IB=INT(UP1(LU- 8)*1.0001)
      IQ=INT(UP1(LU-10)*1.0001)
      IS=INT(UP1(LU-9)*1.0001)
      PX=UP1(LU-3)
      PY=UP1(LU-2)
      PZS=UP1(LU- 1)
      PS=SQRT(PX**2+PY**2+PZS**2)
      ES=SQRT(PS**2+CM**2)
      PZ=GSYS*(PZS+ES*VSYS)
      P=SQRT(PX**2+PY**2+PZ**2)
c kkg 11/14/05 -> SGM, 11/28/05
            if(P.lt.1.0d-6) P=1.0d-6
      E=SQRT(P**2+CM**2)
c  kkg  12/14/05
      if(IS.eq.0.and.IB.ne.0.and.M.gt.MV)  then
c                                 ! cascade particle  are excluded
c   for IGSI=1 only fragments from projetile fragmentation are included
c   for IGSI=2 only fragments from target    fragmentation are included
c   for IGSI=3 all  fragments are included
        if((IGSI.eq.1.and.M.gt.MV1).or.(IGSI.eq.2.and.M.le.MV1))
     &                          go  to  5
        if(IGSI.eq.1)  then
c    tkfr,tet,vz are calculated in the projectile rest frame
          PZRF = GLA*(PZ-E*VLA)
          PRF  =SQRT(PX**2+PY**2+PZRF**2)
          if(PRF.lt.1.0d-6)  PRF=1.0d-6
          ERF   = SQRT(PRF**2+CM**2)
          tkfr  = (ERF - CM)*1000.
          ctrf =  PZRF/PRF
          tet  =  ACOS(ctrf)*180./3.141592
          vz   =  PZRF/ERF
        else                     
c    tkfr,tet,vz are calculated in the target rest frame(lab.system)
          tkfr=(E-CM)*1000.
          ctfr=PZ/P
          tet = ACOS(ctfr)*180./3.141592
          vz  = PZ/E
        endif             
            afr=FLOAT(IB)
            zfr=FLOAT(IQ)
            da = abeg - afr
            dz = zbeg - zfr
            ida = int(da + 0.1) + 1
            ida = min(int(abeg + 0.1), ida)
            ida = min(ida,349)
            idz = int(dz + 0.1) + 1
            idz = min(int(zbeg + 1.1), idz)
            idz = min(idz,149)
c
            if(ida.lt.1.or.idz.lt.1)  then
	         write(16,*) ' afr,zfr=',afr,zfr
               write(*,*) ' afr,zfr=',afr,zfr
            endif
c
            if (afr.ge.zfr)
     &      dadz(ida,idz) = dadz(ida,idz) + 1.0
            dadz(351,idz)=dadz(351,idz) + 1.0
            dadz(350,idz)=dadz(350,idz) + tkfr
            dadz(ida,150)=dadz(ida,150) + 1.0
            dadz(ida,151)=dadz(ida,151) + tkfr
            dadz4(1,ida,idz) = dadz4(1,ida,idz) + tkfr
            dadz4(1,351,idz) = dadz4(1,351,idz) + tet
            dadz4(1,350,idz) = dadz4(1,350,idz) + vz
            dadz4(1,ida,151) = dadz4(1,ida,151) + tet 
            dadz4(1,ida,150) = dadz4(1,ida,150) + vz
            dadz4(2,ida,idz) = dadz4(2,ida,idz) + tkfr**2
            dadz4(2,351,idz) = dadz4(2,351,idz) + tet**2
            dadz4(2,350,idz) = dadz4(2,350,idz) + vz**2
            dadz4(2,ida,151) = dadz4(2,ida,151) + tet**2 
            dadz4(2,ida,150) = dadz4(2,ida,150) + vz**2
            if(tet.le.90.0d0) then
              dadz4(3,ida,idz) = dadz4(3,ida,idz) + 1.0
              dadz4(3,351,idz) = dadz4(3,351,idz) + 1.0
              dadz4(3,350,idz) = dadz4(3,350,idz) + tkfr
              dadz4(3,349,idz) = dadz4(3,349,idz) + tkfr**2
              dadz4(3,ida,150) = dadz4(3,ida,150) + 1.0
              dadz4(3,ida,151) = dadz4(3,ida,151) + tkfr
              dadz4(3,ida,149) = dadz4(3,ida,149) + tkfr**2
	    else   
              dadz4(4,ida,idz) = dadz4(4,ida,idz) + 1.0
              dadz4(4,351,idz) = dadz4(4,351,idz) + 1.0
              dadz4(4,350,idz) = dadz4(4,350,idz) + tkfr
              dadz4(4,349,idz) = dadz4(4,349,idz) + tkfr**2
              dadz4(4,ida,150) = dadz4(4,ida,150) + 1.0
              dadz4(4,ida,151) = dadz4(4,ida,151) + tkfr
              dadz4(4,ida,149) = dadz4(4,ida,149) + tkfr**2
            endif
    5     continue 
      else
c  kkg  08/10/04
        IF(IB.EQ.0.AND.IS .EQ. 0.AND.IQ.EQ.-1.AND.ABS(CM).lt.0.2)
     &  pmult(4,7)=pmult(4,7)+1. ! pi-
        IF(IB.EQ.0.AND.IS .EQ. 0.AND.IQ.EQ. 1.AND.ABS(CM).lt.0.2)
     &  pmult(4,8)=pmult(4,8)+1. ! pi+
        IF(IB.EQ.0.AND.IS .EQ. 0.AND.IQ.EQ. 0.AND.ABS(CM).lt.0.2)
     &  pmult(4,9)=pmult(4,9)+1. ! pi0
      endif
c
c  kkg
      JS=0
c  kkg  08/10/04
C SGM, 04/01/13
      IF(IB.EQ.0.AND.IS .EQ. 0.AND.ABS(CM).GT.9.)  GO  TO  12
      IF(IB.EQ.6.AND.IS .EQ. 0.AND.IQ.EQ. 2)   JS=1    ! He6
      IF(IB.EQ.6.AND.IS .EQ. 0.AND.IQ.EQ. 3)   JS=2    ! Li6
      IF(IB.EQ.7.AND.IS .EQ. 0.AND.IQ.EQ. 3)   JS=3    ! Li7
      IF(IB.EQ.7.AND.IS .EQ. 0.AND.IQ.EQ. 4)   JS=4    ! Be7
      IF(IB.EQ.9.AND.IS .EQ. 0.AND.IQ.EQ. 4)   JS=5    ! Be9
      IF(IB.EQ.10.AND.IS .EQ. 0.AND.IQ.EQ. 4)   JS=6    ! Be10
      IF(IB .EQ.8.AND.IS .EQ. 0.AND.IQ.EQ. 5)   JS=7    ! B8  
      IF(IB.EQ.10.AND.IS .EQ. 0.AND.IQ.EQ. 5)   JS=8    ! B10
c      IF(IB.EQ.6.AND.IS .EQ. 0.AND.IQ.EQ. 3)   JS=9    ! Li6 !SGM, 11/15/10 
      IF(IB.EQ.11.AND.IS .EQ. 0.AND.IQ.EQ. 5)   JS=9    ! B11
      IF(IB.EQ.10.AND.IS .EQ. 0.AND.IQ.EQ. 6)   JS=10   ! C10
      IF(IB.EQ.11.AND.IS .EQ. 0.AND.IQ.EQ. 6)   JS=11   ! C11        
C     IF(JS.EQ.0)   write(16,*) 'GA3,JS=0: Q,B,S,M=',IQ,IB,IS,CM
      IF(JS.EQ.0)   GO  TO  12
C  ---------------------------------------------------------
      MULT(JS)=MULT(JS)+1
c
      if(JS.eq.3)  then
        nn = nn + 1
        if(M.le.MV)                       nnc = nnc + 1
        if(M.gt.MV.and.IORI(3,M).eq.1)    nnp = nnp + 1
        if(.not.fisevent)  then
          if(M.gt.MV.and.IORI(3,M).eq.2)  nne = nne + 1
        else
          if(M.gt.MV.and.IORI(3,M).eq.2)  nnprf = nnprf + 1
          if(M.gt.MV.and.IORI(3,M).eq.6)  nnpof = nnpof + 1
        endif
      endif
      if(JS.ge.3.and.JS.le.8)  then
        if(M.le.MV)  then
                              pmult(4,JS-2)=pmult(4,JS-2)+1.
        else
          if(IORI(3,M).eq.3)  pmult(3,JS-2)=pmult(3,JS-2)+1.
        endif
      endif
      IF(P.LT.1.E-05)  GO  TO  10
      CT=PZ/P
      ST2=1.-CT**2
      IF(ST2.LE.0.)   GO  TO  10
      ST=SQRT(ST2)
      CF=PX/P/ST
      IF(ABS(CF).GT.1.)   GO  TO  13
      SF=PY/P/ST
      GO  TO  11
   10 ST=1.E-17
!      CT=SIGN(1.,PZ)
			if(PZ.LT.0) then
				CT = -1.0
			else
				CT = 1.0
			endif
   13 continue !CF=SIGN(1.,PX)
			if(PX.LT.0) then
				CF = -1.0
			else
				CF = 1.0
			endif
      SF=0.
   11 E=SQRT(P**2+CM**2)
      TETA=ACOS(CT)*180./3.141592
      FI=ACOS(CF)*180./3.141592
      IF(SF.LT.0.)   FI=360.-FI
C__________________________________________________________
      IF(INSP.eq.3)  then
        X=P
        U=TETA
        WI=1.
C__________________________________________________________
      ELSEIF(INSP.eq.2)  then
        X=E-CM
        U=TETA
        WI=1.
C__________________________________________________________
      ELSEIF(INSP.eq.1)  then
        X=P
        U=TETA
        WI=1.
C__________________________________________________________
      ELSE
        go  to  12
C__________________________________________________________
      
      ENDIF
      DO  117  K=1,10
      KU=K
      IF(U.GT.(UF(K)-DU).AND.U.LT.(UF(K)+DU))  GO  TO  116
      GO  TO  117
  116 CONTINUE
C
      if(INSP.eq.2.and.T0A.le.T0le)  then
        if(x.le.tgr(1))      then
          kx=INT(x/dtgr(1))+jgr(1)
        elseif(x.lt.tgr(2))  then
          kx=INT((x-tgr(1))/dtgr(2))+jgr(2)
        elseif(x.lt.tgr(3))  then
          kx=INT((x-tgr(2))/dtgr(3))+jgr(3)
        elseif(x.lt.tgr(4))  then
          kx=INT((x-tgr(3))/dtgr(4))+jgr(4)
        else
c  kkg 12/13/04
          kx=200
        endif
        SPE(JS,KU,kx)=SPE(JS,KU,kx)+WI
      else
        DX=DXS(KU)
c       
c        if(IB.gt.1)  DX=DX*FLOAT(IB)      ! 03/28/05
        XMA=200.*DX+XMIN
        CALL  HIST4(X,XMIN,XMA,DX,SPE,11,11,205,WI,JS,KU)
      endif
  117 CONTINUE
C__________________________________________________________
      if(INSP.eq.2.and.T0A.le.T0le)  then
        if(x.le.tgr(1))      then
          kx=INT(x/dtgr(1))+jgr(1)
        elseif(x.lt.tgr(2))  then
          kx=INT((x-tgr(1))/dtgr(2))+jgr(2)
        elseif(x.lt.tgr(3))  then
          kx=INT((x-tgr(2))/dtgr(3))+jgr(3)
        elseif(x.lt.tgr(4))  then
          kx=INT((x-tgr(3))/dtgr(4))+jgr(4)
        else
c  kkg 11/23/04
          kx=200
        endif
        SPE(JS,11,kx)=SPE(JS,11,kx)+WI
      else
        DX11=DXS(1)
c
c        if(IB.gt.1)  DX11=DX11*FLOAT(IB)    ! 03/28/05       
        XMA=200.*DX11+XMIN
        w = 1.0d0
        CALL  HIST4(X,XMIN,XMA,DX11,SPE,11,11,205,w,JS,11)
      endif
c
      Tkin=E-CM
c     if(Tkin.ge.ewin(1,JS).and.Tkin.le.ewin(2,JS))  then
        x=TETA
        xmin=0.
        xmax=180.
        dx=dthw 
        w = 1.0d0
        call  HIST3(x,xmin,xmax,dx,angw,25,11,w,JS)
c     endif
c
   12 CONTINUE
c   kkg  01/18/05   Accumulate residual nuclei information
        fisevent = Nfiss.eq.1
        call  rezdisq(zbeg,abeg,fisevent,nn)  
c      Accumulate distribution of fission fragments oppening angle
        call  opandiq(fisevent,nn) 
c      Accumulate distribution of neutron multiplicity
        call  disnmuq(fisevent,nn,nnc,nnp,nne,nnprf,nnpof)
      RETURN
      END
      
C  *************************************************************
C
	subroutine rezdisq(zbeg,abeg,fisevent,nn)
      implicit real*8 (a-h, o-z), integer (i-n)
      logical  fisevent
c   kkg 12/02/04
      common /residf/  atf, ztf, exf, pmf  
      common /resid0/  at0(2),zt0(2),ex0(2),pm0(2),am0(2)
	common /rezdis/  rdis(5,5,250), dex, dpm
	common /fisopa/  opan(7,185),dth12
      common /ifiss/   af12(2),zf12(2),tf12(2),ex12(2),bf12(2,3), ifiss
      common /degrad/  degrad

      data zro, one/0.d0, 1.d0/

c  kkg  12/02/04
c         Accumulate residual nuclei information
c         nuclei after cascade
        da = abeg - at0(1)
        dz = zbeg - zt0(1)
        ia = nint(da) + 1
        ia = min(nint(abeg), ia)
        ia = min(ia, 247)
        iz = nint(dz) + 1
        iz = min(nint(zbeg) + 1, iz)
        iz = min(iz, 247)
	  iex= int(ex0(1)/dex) + 1
c   KKG 05/23/07
        pm0a=pm0(1)
        if(at0(1).gt.one) pm0a=pm0(1)/at0(1)
        ipm= int(pm0a/dpm) + 1
	  iam= int(am0(1))     + 1
	  if(ia.le.247)  then
	    rdis(1,1,ia) = rdis(1,1,ia) + one
		rdis(1,1,248)= rdis(1,1,248)+ at0(1)
		rdis(1,1,249)= rdis(1,1,249)+ at0(1)**2
		rdis(1,1,250)= rdis(1,1,250)+ one
        endif
	  if(iz.le.247)  then
	    rdis(1,2,iz) = rdis(1,2,iz) + one
		rdis(1,2,248)= rdis(1,2,248)+ zt0(1)
		rdis(1,2,249)= rdis(1,2,249)+ zt0(1)**2
		rdis(1,2,250)= rdis(1,2,250)+ one
        endif
	  if(iex.gt.0.and.iex.le.247)  then
	    rdis(1,3,iex) = rdis(1,3,iex) + one
		rdis(1,3,248) = rdis(1,3,248)+ ex0(1)
		rdis(1,3,249) = rdis(1,3,249)+ ex0(1)**2
		rdis(1,3,250) = rdis(1,3,250)+ one
        endif
	  if(ipm.gt.0.and.ipm.le.247)  then
	    rdis(1,4,ipm) = rdis(1,4,ipm) + one
c   KKG 05/23/07
		rdis(1,4,248) = rdis(1,4,248)+ pm0a
		rdis(1,4,249) = rdis(1,4,249)+ pm0a**2
		rdis(1,4,250) = rdis(1,4,250)+ one
        endif
	  if(iam.gt.0.and.iam.le.247)  then
	    rdis(1,5,iam) = rdis(1,5,iam) + one
		rdis(1,5,248) = rdis(1,5,248)+ am0(1)
		rdis(1,5,249) = rdis(1,5,249)+ am0(1)**2
		rdis(1,5,250) = rdis(1,5,250)+ one
        endif
c         nuclei after preequilibrium
        da = abeg - at0(2)
        dz = zbeg - zt0(2)
        ia = nint(da) + 1
        ia = min(nint(abeg), ia)
        ia = min(ia, 247)
        iz = nint(dz) + 1
        iz = min(nint(zbeg) + 1, iz)
        iz = min(iz, 247)
	  iex= int(ex0(2)/dex) + 1
c   KKG 05/23/07
        pm0a=pm0(2)
        if(at0(2).gt.one) pm0a=pm0(2)/at0(2)
        ipm= int(pm0a/dpm) + 1
	  iam= int(am0(2))     + 1
	  if(ia.le.247)  then
	    rdis(2,1,ia) = rdis(2,1,ia) + one
		rdis(2,1,248)= rdis(2,1,248)+ at0(2)
		rdis(2,1,249)= rdis(2,1,249)+ at0(2)**2
		rdis(2,1,250)= rdis(2,1,250)+ one
        endif
	  if(iz.le.247)  then
	    rdis(2,2,iz) = rdis(2,2,iz) + one
		rdis(2,2,248)= rdis(2,2,248)+ zt0(2)
		rdis(2,2,249)= rdis(2,2,249)+ zt0(2)**2
		rdis(2,2,250)= rdis(2,2,250)+ one
        endif
	  if(iex.gt.0.and.iex.le.247)  then
	    rdis(2,3,iex) = rdis(2,3,iex) + one
		rdis(2,3,248) = rdis(2,3,248)+ ex0(2)
		rdis(2,3,249) = rdis(2,3,249)+ ex0(2)**2
		rdis(2,3,250) = rdis(2,3,250)+ one
        endif
	  if(ipm.gt.0.and.ipm.le.247)  then
	    rdis(2,4,ipm) = rdis(2,4,ipm) + one
		rdis(2,4,248) = rdis(2,4,248)+ pm0a
		rdis(2,4,249) = rdis(2,4,249)+ pm0a**2
		rdis(2,4,250) = rdis(2,4,250)+ one
        endif
	  if(iam.gt.0.and.iam.le.247)  then
	    rdis(2,5,iam) = rdis(2,5,iam) + one
		rdis(2,5,248) = rdis(2,5,248)+ am0(2)
		rdis(2,5,249) = rdis(2,5,249)+ am0(2)**2
		rdis(2,5,250) = rdis(2,5,250)+ one
        endif
        if(.not. fisevent)  then
c         nuclei decaying by evaporation only
	   if(ia.le.247)  then
	    rdis(3,1,ia) = rdis(3,1,ia) + one
		rdis(3,1,248)= rdis(3,1,248)+ at0(2)
		rdis(3,1,249)= rdis(3,1,249)+ at0(2)**2
		rdis(3,1,250)= rdis(3,1,250)+ one
         endif
	   if(iz.le.247)  then
	    rdis(3,2,iz) = rdis(3,2,iz) + one
		rdis(3,2,248)= rdis(3,2,248)+ zt0(2)
		rdis(3,2,249)= rdis(3,2,249)+ zt0(2)**2
		rdis(3,2,250)= rdis(3,2,250)+ one
         endif
	   if(iex.gt.0.and.iex.le.247)  then
	    rdis(3,3,iex) = rdis(3,3,iex) + one
		rdis(3,3,248) = rdis(3,3,248)+ ex0(2)
		rdis(3,3,249) = rdis(3,3,249)+ ex0(2)**2
		rdis(3,3,250) = rdis(3,3,250)+ one
         endif
	   if(ipm.gt.0.and.ipm.le.247)  then
	    rdis(3,4,ipm) = rdis(3,4,ipm) + one
c   KKG  05/23/07
		rdis(3,4,248) = rdis(3,4,248)+ pm0a
		rdis(3,4,249) = rdis(3,4,249)+ pm0a**2
		rdis(3,4,250) = rdis(3,4,250)+ one
         endif
	   if(iam.gt.0.and.iam.le.247)  then
	    rdis(3,5,iam) = rdis(3,5,iam) + one
		rdis(3,5,248) = rdis(3,5,248)+ am0(2)
		rdis(3,5,249) = rdis(3,5,249)+ am0(2)**2
		rdis(3,5,250) = rdis(3,5,250)+ one
         endif
        else 
c         nuclei decaying by fission
	   if(ia.le.247)  then
	    rdis(4,1,ia) = rdis(4,1,ia) + one
		rdis(4,1,248)= rdis(4,1,248)+ at0(2)
		rdis(4,1,249)= rdis(4,1,249)+ at0(2)**2
		rdis(4,1,250)= rdis(4,1,250)+ one
         endif
	   if(iz.le.247)  then
	    rdis(4,2,iz) = rdis(4,2,iz) + one
		rdis(4,2,248)= rdis(4,2,248)+ zt0(2)
		rdis(4,2,249)= rdis(4,2,249)+ zt0(2)**2
		rdis(4,2,250)= rdis(4,2,250)+ one
         endif
	   if(iex.gt.0.and.iex.le.247)  then
	    rdis(4,3,iex) = rdis(4,3,iex) + one
		rdis(4,3,248) = rdis(4,3,248)+ ex0(2)
		rdis(4,3,249) = rdis(4,3,249)+ ex0(2)**2
		rdis(4,3,250) = rdis(4,3,250)+ one
         endif
	   if(ipm.gt.0.and.ipm.le.247)  then
	    rdis(4,4,ipm) = rdis(4,4,ipm) + one
c   KKG  05/23/07
		rdis(4,4,248) = rdis(4,4,248)+ pm0a
		rdis(4,4,249) = rdis(4,4,249)+ pm0a**2
		rdis(4,4,250) = rdis(4,4,250)+ one
         endif
	   if(iam.gt.0.and.iam.le.247)  then
	    rdis(4,5,iam) = rdis(4,5,iam) + one
		rdis(4,5,248) = rdis(4,5,248)+ am0(2)
		rdis(4,5,249) = rdis(4,5,249)+ am0(2)**2
		rdis(4,5,250) = rdis(4,5,250)+ one
         endif
c  nuclei just before fission
c  in GEM2 angular momentum is not changed
         da = abeg - atf
         dz = zbeg - ztf
         ia = nint(da) + 1
         ia = min(nint(abeg), ia)
         ia = min(ia, 247)
         iz = nint(dz) + 1
         iz = min(nint(zbeg) + 1, iz)
         iz = min(iz, 247)
	   iex= int(exf/dex) + 1
c   KKG  05/23/07
         pmfa = pmf
	   if(atf.gt.one)  pmfa = pmf/atf 
	   ipm= int(pmfa/dpm) + 1
	   iam= int(am0(2))     + 1
	   if(ia.le.247)  then
	    rdis(5,1,ia) = rdis(5,1,ia) + one
		rdis(5,1,248)= rdis(5,1,248)+ atf
		rdis(5,1,249)= rdis(5,1,249)+ atf**2
		rdis(5,1,250)= rdis(5,1,250)+ one
         endif
	   if(iz.le.247)  then
	    rdis(5,2,iz) = rdis(5,2,iz) + one
		rdis(5,2,248)= rdis(5,2,248)+ ztf
		rdis(5,2,249)= rdis(5,2,249)+ ztf
		rdis(5,2,250)= rdis(5,2,250)+ one
         endif
	   if(iex.gt.0.and.iex.le.247)  then
	    rdis(5,3,iex) = rdis(5,3,iex) + one
		rdis(5,3,248) = rdis(5,3,248)+ exf
		rdis(5,3,249) = rdis(5,3,249)+ exf**2
		rdis(5,3,250) = rdis(5,3,250)+ one
         endif
	   if(ipm.gt.0.and.ipm.le.247)  then
	    rdis(5,4,ipm) = rdis(5,4,ipm) + one
c   KKG  05/23/07
		rdis(5,4,248) = rdis(5,4,248)+ pmfa
		rdis(5,4,249) = rdis(5,4,249)+ pmfa**2
		rdis(5,4,250) = rdis(5,4,250)+ one
         endif
	   if(iam.gt.0.and.iam.le.247)  then
	    rdis(5,5,iam) = rdis(5,5,iam) + one
		rdis(5,5,248) = rdis(5,5,248)+ am0(2)
		rdis(5,5,249) = rdis(5,5,249)+ am0(2)**2
		rdis(5,5,250) = rdis(5,5,250)+ one
         endif
        endif
	return   	
c     
	entry opandiq(fisevent,nn) 	  
c      Accumulate distribution of fission fragments oppening angle
c      nn is total number of produced neutrons
      if(.not.fisevent)  return
        b1 = sqrt(bf12(1,1)**2 + bf12(1,2)**2 + bf12(1,3)**2) 
        b2 = sqrt(bf12(2,1)**2 + bf12(2,2)**2 + bf12(2,3)**2) 
        ct12 = (bf12(1,1)*bf12(2,1)+bf12(1,2)*bf12(2,2)+
     +          bf12(1,3)*bf12(2,3))/b1/b2
	  st12 = sqrt(abs(one-ct12**2))
	  t12 = atan2(st12,ct12)*degrad
	  it12= int(t12/dth12) + 1
	  opan(1,it12) = opan(1,it12) + one
	  opan(1,183)  = opan(1,183)  + t12
	  opan(1,184)  = opan(1,184)  + t12*t12
	  opan(1,185)  = opan(1,185)  + one
	  if(nn.le.5)  then
	    opan(2,it12) = opan(2,it12) + one
	    opan(2,183)  = opan(2,183)  + t12
	    opan(2,184)  = opan(2,184)  + t12*t12
	    opan(2,185)  = opan(2,185)  + one
        elseif(nn.le.8)  then
	    opan(3,it12) = opan(3,it12) + one
	    opan(3,183)  = opan(3,183)  + t12
	    opan(3,184)  = opan(3,184)  + t12*t12
	    opan(3,185)  = opan(3,185)  + one
        elseif(nn.le.12)  then
	    opan(4,it12) = opan(4,it12) + one
	    opan(4,183)  = opan(4,183)  + t12
	    opan(4,184)  = opan(4,184)  + t12*t12
	    opan(4,185)  = opan(4,185)  + one
        elseif(nn.le.15)  then
	    opan(5,it12) = opan(5,it12) + one
	    opan(5,183)  = opan(5,183)  + t12
	    opan(5,184)  = opan(5,184)  + t12*t12
	    opan(5,185)  = opan(5,185)  + one
        elseif(nn.le.19)  then
	    opan(6,it12) = opan(6,it12) + one
	    opan(6,183)  = opan(6,183)  + t12
	    opan(6,184)  = opan(6,184)  + t12*t12
	    opan(6,185)  = opan(6,185)  + one
        elseif(nn.ge.20)  then
	    opan(7,it12) = opan(7,it12) + one
	    opan(7,183)  = opan(7,183)  + t12
	    opan(7,184)  = opan(7,184)  + t12*t12
	    opan(7,185)  = opan(7,185)  + one
        endif
	return 
c ======================================================================
      end      
C
C
C  *************************************************************
C
      subroutine  disnmuq(fisevent,nn,nnc,nnp,nne,nnprf,nnpof)
c   This subroutine accumulate the distribution of the multiplicity of
c   produced neutrons
c   nn   -  total number of neutrons produced in event
c   nnc  -  number of neutrons produced in cascade stage
c   nnp  -  number of neutrons produced in preequilibrium stage
c   nne  -  number of neutrons produced in evaporation(no fission) stage
c   nnprf-  number of prefission neutrons produced in evaporation stage
c   nnpof-  number of postfission neutrons produced in evaporation of
c           fission fragments
c
      implicit real*8 (a-h, o-z), integer (i-n)
      logical  fisevent
      common /disnmu/ disnm(6,155)
	dimension n(6)
      data zro, one/0.d0, 1.d0/
c
      do  k=1,6
	  if(k.eq.1) n(k) = nn
	  if(k.eq.2) n(k) = nnc
	  if(k.eq.3) n(k) = nnp
	  if(k.eq.4) n(k) = nne
	  if(k.eq.5) n(k) = nnprf
	  if(k.eq.6) n(k) = nnpof
	  if((k.le.3).or.(k.eq.4.and. .not. fisevent).or.
     &     (k.ge.5.and. fisevent))  then
	    if(n(k).ge.0.and.n(k).le.151)  then
	      disnm(k,n(k)+1) = disnm(k,n(k)+1) + one
	      disnm(k,153)  = disnm(k,153)  + dble(n(k))  
	      disnm(k,154)  = disnm(k,154)  + dble(n(k)*n(k))  
	      disnm(k,155)  = disnm(k,155)  + one  
          endif
        endif
	enddo     
	return 
c ======================================================================
      end  
CC  *************************************************************
C
      SUBROUTINE HIST1(X,A,B,H,RX,N,W)
C     BLOCK FOR BILDING OF HISTOGRAMS.
      implicit real*8 (a-h, o-z), integer (i-n)
      DIMENSION RX(N)
      RX(N) = RX(N)+X*X*W
      IF(X.lt.A) then
        RX(N-4) = RX(N-4)+W
        RETURN
      ENDIF
      IF (X.ge.B) then
        RX(N-3) = RX(N-3)+W
      ELSE
        L=(X-A)/H
        NL=N-5
        IF(L.GT.NL) then
          write(16,6)  X,A,B,H,W
    6 FORMAT(25X,'MISTAKE IN DIMENSION OF HIST '/15X,5(F10.3,2X))
          RETURN
	ENDIF
      RX(L+1)=RX(L+1)+W
      RX(N-1)=RX(N-1)+X*W
      RX(N-2)=RX(N-2)+W
      ENDIF
      RETURN
      END
C**********************************************************************
      SUBROUTINE HIST2(X,A,B,H,RX,N,W,J)
C       BLOCK FOR BUILDING 2-DIMENSIONAL HISTOGRAMS
      implicit real*8 (a-h, o-z), integer (i-n)
      DIMENSION  RX(N,2)
      RX(N,J) = RX(N,J)+X*X*W
      IF(X.lt.A) then
        RX(N-4,J) = RX(N-4,J)+W
        RETURN
      ENDIF
      IF(X.ge.B) then
        RX(N-3,J) = RX(N-3,J)+W
      ELSE
        L = (X-A)/H
        NL=N-5
        IF(L.GT.NL) then
          write(16,6)  X,A,B,H,W
    6 FORMAT(25X,'MISTAKE IN DIMENSION OF HIST2 '/15X,5(F10.3,2X))
          RETURN
	ENDIF
        RX(L+1,J) = RX(L+1,J)+W
        RX(N-1,J) = RX(N-1,J)+X
        RX(N-2,J) = RX(N-2,J)+W
      ENDIF
      RETURN
      END
C**********************************************************************
      SUBROUTINE HIST3(X,A,B,H,RX,N,M,W,J)
C       BLOCK FOR BUILDING M-DIMENSIONAL HISTOGRAMS
      implicit real*8 (a-h, o-z), integer (i-n)
      DIMENSION  RX(N,M)
C
      IF(X.lt.A) then
        RX(N-4,J) = RX(N-4,J)+W
        RETURN
      ENDIF
      IF(X.ge.B) then
      RX(N-3,J) = RX(N-3,J)+W
      ELSE
        L=(X-A)/H
        NL=N-5
        IF(L.GT.NL) then
          write(16,6)  X,A,B,H,W,J
    6 FORMAT(25X,'MISTAKE IN DIMENSION OF HIST3'/1X,5(F10.3,3X),I3)
          RETURN
	ENDIF
        RX(L+1,J)=RX(L+1,J)+W
        RX(N,J)=RX(N,J)+X*X*W
        RX(N-1,J)=RX(N-1,J)+X*W
        RX(N-2,J)=RX(N-2,J)+W
      ENDIF
      RETURN
      END
C**********************************************************************
      SUBROUTINE HIST4(X,A,B,H,RX,N,M,L,W,I,J)
      implicit real*8 (a-h, o-z), integer (i-n)
      DIMENSION RX(N,M,L)
      RX(I,J,L)=RX(I,J,L)+X*X*W
      IF(X.lt.A) then
        RX(I,J,L-4)=RX(I,J,L-4)+W
        RETURN
      ENDIF
      IF(X.ge.B) then
        RX(I,J,L-3)=RX(I,J,L-3)+W
      ELSE
        L1=(X-A)/H
        NL=L-5
        IF(L1.GT.NL) then
          write(16,6)  X,A,B,H,W,I,J
    6 FORMAT(25X,'MISTAKE IN DIMENSION OF HIST4'/1X,5(F10.3,3X),2I3)
          RETURN
	ENDIF
        RX(I,J,L1+1) = RX(I,J,L1+1)+W
        RX(I,J,L-1)=RX(I,J,L-1)+X*W
        RX(I,J,L-2)=RX(I,J,L-2)+W
      ENDIF
      RETURN
      END
C**********************************************************************
      SUBROUTINE PARAM(L,LIM1,LPCON)
      implicit real*8 (a-h, o-z), integer (i-n)
c      REAL*4 UP1,UP2
      REAL*8 MMES,MBAR
      COMMON /UPAC/UP1(66000),IPER
      COMMON /UPAC2/UP2(70000),LU2
      IF(L.GT.LIM1.OR.LPCON.NE.1) RETURN
      PIDABS=UP2(44)
      RCOR  =UP2(45)
      MMES  =UP2(46)
      MBAR  =UP2(47)
      STIN  =UP2(48)
      AMIN  =UP2(49)
      !DSG=SIGN(0.01,UP2(50))
			if(UP2(50).LT.0) then
				DSG = -0.01
			else
				DSG = 0.01
			endif
      I50   =INT(UP2(50)+DSG)
      !DSG=SIGN(0.01,UP2(51))
			if(UP2(51).LT.0) then
				DSG = -0.01
			else
				DSG = 0.01
			endif
      I51   =INT(UP2(51)+DSG)
      write(16,100) PIDABS,RCOR,MMES,MBAR,STIN,AMIN,I50,I51
  100 FORMAT(2X,'PARAMETERS: PIDABS=',F5.3,2X,'RCOR=',F5.3,2X,
     *'MMES=',F7.3,2X,'MBAR=',F7.3/
     *14X,'STIN=',F5.0,2X,'AMIN=',F5.3/
     *14X,'INDINT*10+INTCC=',I3/
     *14X,'KSYST*10000+IACT*1000+INDEC*100+ISOB2*10+ISOB3=',I6//)
      RETURN
      END
      

      SUBROUTINE  UPPME(MV)
      implicit real*8 (a-h, o-z), integer (i-n)
c      real*4 UP1,UP2
      COMMON /UPAC/UP1(66000),IPER
      COMMON /UPAC2/UP2(70000),LU2
      COMMON /MEMORY/ PM(9,5999),IM(5,5999)
      COMMON/RESUL4/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      COMMON/AZT/AT(2),ZT(2)
      COMMON/IDPME/ IDPME(5999)
      COMMON/PORIG/ IORI(3,5999)
      CHARACTER*8 PNAR
      IF(MV.GT.5999)  write(16,'(2X,''UPPME : MV>5999'',I7)') MV
      IF(MV.EQ.0)   RETURN
      DO  M=1,MV
        LU=11*M
	IDPME(M)=NINT(UP1(LU-10))           ! particle identificator
	IORI(1,M)=NINT(UP1(LU-9))           ! 1-st parent ID
        IORI(3,M)=NINT(ABS(UP1(LU- 8)))/10**4  ! number of rescatter.
	icode=(NINT(ABS(UP1(LU-8)))-IORI(3,M)*10**4)
	IORI(2,M)=ISIGN(icode,INT(UP1(LU-8)))  ! 2-nd parent ID
*                              for resonance decay IORI(2,M)=0
        AM=UP1(LU- 7)
        IF(IDPME(M).EQ.0)       THEN
C  cascade Deuteron
	  IM(1,M)=1
	  IM(2,M)=0
	  IM(3,M)=0
	  IM(4,M)=2
	  IM(5,M)=0
	ELSEIF(IDPME(M).EQ.10)  THEN
C  photon
	  IM(1,M)=0
	  IM(2,M)=0
	  IM(3,M)=0
	  IM(4,M)=0
	  IM(5,M)=0
c	  AM=0.
	ELSE
          IM(1,M)=INT(1.001*CHARGE(IDPME(M)))    ! electric charge
          IM(2,M)=0                           ! (for leptonic ID)
          IM(3,M)=IS(IDPME(M))                ! strangeness
          IM(4,M)=IB(IDPME(M))                ! baryonic charge
          IM(5,M)=0                       ! free (flag of resonance ?)
	ENDIF
        IF(ABS(AM).GT.9..AND.IM(4,M).EQ.0) THEN
          IM(2,M)=NINT(AM)                       ! for w=1. leptons
c         AM=0.                                  ! for w=1. leptons
	ENDIF
        PM(1,M)=UP1(LU- 6)                  ! X   (fm)
        PM(2,M)=UP1(LU- 5)                  ! Y
        PM(3,M)=UP1(LU- 4)                  ! Z
        PM(4,M)=UP1(LU- 3)                  ! p_x  (GeV/c)
        PM(5,M)=UP1(LU- 2)                  ! p_y
        PM(6,M)=UP1(LU- 1)                  ! p_z
        PM(7,M)=UP1(LU   )                  ! weight
        PM(8,M)=SQRT(PM(4,M)**2+PM(5,M)**2+PM(6,M)**2+AM**2)-AM  ! Tkin
        PM(9,M)=AM 
c
c        write(16,100) M,(PM(j,M),j=1,9),(IM(j,M),j=1,5)
  100   format(1x,I5,9(1PE11.4),5I5)	  	                           ! mass (in GeV)
      ENDDO
c      write(16,*) ' UPPME: MV=',MV
      RETURN
C   * * * * * * * * * * * * * * * * * * * * * * **
      ENTRY  PMEUP(MV)
      if(MV.gt.0) then
        IF((11*MV).GT.66000)
     &  write(16,'(2X,''PMEUP: 11*MV>66000'')') MV
        DO  M=1,MV
          LU=11*M
c      ---------------  for output arrays ---------------
          UP1(LU-10)=FLOAT(IM(1,M))     ! electric charge
          UP1(LU- 9)=FLOAT(IM(3,M))     ! strangeness
          UP1(LU- 8)=FLOAT(IM(4,M))     ! baryonic charge
          UP1(LU- 7)=PM(9,M)            ! mass (GeV)
          UP1(LU- 6)=PM(1,M)            ! X (fm)
          UP1(LU- 5)=PM(2,M)            ! Y
          UP1(LU- 4)=PM(3,M)            ! Z
          UP1(LU- 3)=PM(4,M)            ! p_x (GeV/c)
          UP1(LU- 2)=PM(5,M)            ! p_y
          UP1(LU- 1)=PM(6,M)            ! p_z
          UP1(LU   )=PM(7,M)            ! weight
c        write(16,100) M,(PM(j,M),j=1,9),(IM(j,M),j=1,5)
        ENDDO
c      write(16,*) ' PMEUP: MV=',MV
      endif
      RETURN
      END

C   * * * * * * * * * * * * * * * * * * * * * * **