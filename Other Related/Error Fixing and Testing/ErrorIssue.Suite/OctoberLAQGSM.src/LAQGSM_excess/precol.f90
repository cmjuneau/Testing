
      subroutine precol (enext, atwght, charge, pnx, pny, pnz, elx, ely, &
                      &  elz, rm, kstart, n0, np0, nh0, npz0, wam)

! ======================================================================
!
!   Calculates pre-equlibrium and equlibrium particle emission.
!   Various parts removed to separate subroutines to streamline
!   the structure.  10/20/03
!
!   Called by: LAQGSM
!
!   Calls: AUXL DELTAM EQDECY FAM MOLNIX PEQEMT PREQAUX STSTCS TKINM3
!          VHELP
!
!   CEM95 written by S. G. Mashnik
!   Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!   Edited by A. J. Sierk,  LANL  T-2  November, 1997.
!   Modified by AJS, February-March, 1999.
!   Modified by SGM to include reduced masses, 1999
!   Modified by SGM at 06/16/2000
!   Last modification by SGM of 04/27/01, to get the module for CEM2k
!   "Last" change: 13-AUG-2003 by NVMokhov
!   Modified by A. J. Sierk, LANL T-16, October, 2003.  
!   Modified by K. K. Gudima, December, 2004.  
!   Edited by AJS, January, 2005.
!   Modified     : 02-Aug-2005 by NVM
!   Updated      : 07-Jun-2006 by NVM according to 05/11/06, kkg, rep, ajs, sgm
!   Edited by SGM 07/09/06 to account for the KKG 06/23/06
!   changes to use Fermi break-up model in Preco and Evap when A < 13
!   Updated : Mar-2007 by KKG for MARS
! ======================================================================

      implicit real*8 (a-h, o-z), integer (i-n)

      real*8 molnix

      logical first, preqgo

! ======================================================================

      common /ato3rd/  a3rd(500)
      common /azi/     iz, in
      common /bl1003/  u, a, z
      common /bl1008/  dl, dlmn(6)
      common /bl1018/  afjthr(7), athrd
      common /blac/    ac
      common /blexn/   exn
      common /blr0/    r0
! NVM      common /counter/ icntr
      common /dele/    sfu, wf, fusion, sigfw
      common /fiss5/   z5, a5, u5, ifis5
      common /gemspt/  ares, zres, ures, eres, sptg(6,300), numf 
      common /kktot/   ktot
      common /mnuc/    emnucm, emnucg, emnucb, emnuct, emneut, emprot
      common /nopreco/ sigpre, noprec
      common /probf/   probf 
      common /resid/   angmom(3), v(3), remn

      common /stopr/   istp
      common /trec/    trec
!   KKG  12/01/04
      common /resid0/  at0(2), zt0(2), ex0(2), pm0(2), am0(2)
      COMMON/DELEN/    Wfiss,Nfiss
      dimension pnuc(3)

      data zro, one, two, thr, for, thsn /0.d0, 1.d0, 2.d0, 3.d0, 4.d0, 1.d3/
      data first /.true./

      save first
      
! ======================================================================

      if (first) then
         first = .false.
         r0 = rm
      endif
      fusion = -one
      istart = 1
      u = enext
      a = atwght
      ia = nint(atwght)
      athrd = a3rd(ia)
      z = charge
      un = a - z
      in = nint(un)
      iz = nint(z)
      in = max(1,in)
      iz = max(1,iz)
      emx = molnix (iz, in, 2)
      if (iz.gt.7 .and. in.gt.7) then
        remn = a*emnucb + emx/thsn
      else
        remn = a*emnuct + emx/thsn
      endif
!  Total recoil energy and kineti! energy of entire nucleus
      e = sqrt(pnx**2 + pny**2 + pnz**2 + remn**2)
      trec = (e - remn)*thsn
      if (e < 1.0d-12 .and. e > -1.0d-12) then
        e = 1.0d-12
        print *, 'Divide by zero error in precof.f90 line 106,107,108'
      end if
      v(1) = pnx/e
      v(2) = pny/e
      v(3) = pnz/e
      angmom(1) = elx
      angmom(2) = ely
      angmom(3) = elz
      ktot = kstart - 1
      n = n0
      h = dble(nh0)
      p = dble(np0)
      pz = dble(npz0)
!      u = u - tre!         KKG   05/22/07    
      call auxl (angmom, bf0, ln, erotev, delu)
      u = u + delu
      pevap = molnix (iz, in, 3)
      ue = u - pevap - erotev
      if (ue.lt.thr) then
!  Keep track of nuclei returned with less than 3 MeV of E*:
!       call ststcs (a, z, ue, ln, bf0, 2)
!   Also add them to the residual nuclei total:
!       call ststcs (a, z, trec, ln, bf0, 4)
!        
        call restorl (a, z, ue, ktot)
        kstart = ktot + 1
        return
      endif
!   KKG 12/01/04; for distribution of A, Z, E*, P, L after the cascade:
       at0(1) = a
       zt0(1) = z
       ex0(1) = ue
       pm0(1) = sqrt(pnx**2 + pny**2 + pnz**2)*thsn
       am0(1) = dble(ln)
!   Beginning of main decay loop:
!   SGM, 07/23/07       do 100 k = kstart,150         
!      do 100 k = kstart,300  !SGM 07/23/07
      do 100 k = kstart,5999  !SGM 06/11/08
!   kkg  06/26/06
! SGM 01/12/09          if (a.le.12.d0) then
! SGM 01/12/09
          unn = a - z
		if (a.le.12.d0.or.z.lt.one.or.unn.lt.one) then
! SGM 01/12/09
            uf = ue/thsn
            pnuc(1) = pnx
            pnuc(2) = pny
            pnuc(3) = pnz
            ksf = ktot
            call  FERMIQ(a,z,uf,pnuc,ksf,1)
            ktot = ksf 
            kstart = ktot + 1
            a = zro
            z = zro
            return
          endif
        dl = molnix (iz, in, 2)
        call vhelp
        pevap = molnix (iz, in, 3)
!   ue = "Thermal" energy of nucleus at rotating ground state:
        ue = u - pevap - erotev
        if (ue.le.0.1d0) then
          call restorl (a, z, ue, ktot)
          kstart = ktot + 1
          return
        endif
!   Set up auxiliary quantities for finding prequilibrium
!   emission widths:
        am = fam (a, z, ue, 0)
        ac = 0.595d0*am
        call preqaux (a, z, am, erotev, rr)
!  nsp is the number of excitons at which a compound nucleus is
!  assumed to be formed.
        nsp = int(sqrt(1.19d0*am*a*u + 0.5d0))
        preqgo = noprec.eq.0 .and. rr.gt.zro
        lm = 0
   10   continue
!  Restore old logic of sharp transition when n.gt.nsp AJS 07/08/05.
        rnp = dble(n)/dble(nsp)
        ppre = one - exp(-0.5d0*((rnp - one)/sigpre)**2)
!  Dick's proposal, KKG, 05/23/07 
        dick = rndm(-1.)
        if (n.lt.nsp .and. dick.lt.ppre .and. preqgo) then

!  **************** Pre-equilibrium emission *************************
          ipflg = 0
          exn = dble(n)
          call peqemtl (n, p, h, pz, lm, ipflg, rm, pnx, &
                    &  pny, pnz, bf0, ln, erotev)
!  ipflg = 1 if number of excitons has changed;
          if (ipflg.eq.1) go to 10
!  Otherwise, prequilibrium particle of type lm was emitted.
        else
!  ****************** equilibrium emission ***************************

!   Calculate statistics on average E*, Z, A, L when nucleus enters
!   statistical decay phase.
 
          if (istart.eq.1) then
!   KKG 12/01/04; for distribution of A, Z, E*, P, L after preeq. decay:
            at0(2) = a
            zt0(2) = z
            ex0(2) = ue
            pm0(2) = sqrt(pnx**2 + pny**2 + pnz**2)*thsn
            am0(2) = dble(ln)
!           call ststcs (a, z, ue, ln, bf0, 3)
            istart = 0
          endif
!  kkg 01/29/04
            a0 = a
            z0 = z
!            u0 = ue          !  KKG 05/23/07 
          gasum = zro
          gzsum = zro
          call eqdecyl (a, z, ue, pnx, pny, pnz, gasum, gzsum)
          wf = probf
	  Wfiss=wf
	  Nfiss=ifis5
          if (ifis5.eq.1) then
!   Note that since GEM2 ignores angular momentum, the actual L of
!   the fissioning nucleus will not be identical to ln.  Since we
!   are going to replace GEM2, ignore this discrepancy. AJS 10/28/03.
            fusion = one
!    kkg 12/14/05
            u = zro
            a = zro
            z = zro

!           call ststcs (a5, z5, u5, ln, bf0, 5)
          else
            fusion = -one
            u = ue
            a = ares
            z = zres
            gasum = gasum + a
            gzsum = gzsum + z
            trec = eres
!           call ststcs (a, z, trec, ln, bf0, 4) 
          endif
!  kkg 01/29/04
          if(abs(gasum-a0).gt.0.1.or.abs(gzsum-z0).gt.0.1) then
            write(*,*) 'eq.dec: a0,gasum,z0,gzsum,ares,zres=', &
          & a0,gasum,z0,gzsum,ares,zres 
          endif                 
  
          kstart = ktot + 1
          return
        endif
  100   continue
      wam = -13.d0
      write (16, 1200)  u, a, z
      write ( *, 1200)  u, a, z
      kstart = ktot + 1
      return

! ======================================================================

! 1000 format (1x,'Icntr, u, trec, pnx, pny, pnz, remn are:'/i9,6f14.8)
! 1010 format (1x,'Icntr, ue, pevap, erotev are:'/i9,3f14.8)
 1200 format (/5x,'In PRECOL, number of emitted particles exceeds ', &
     &       '5999 after evaporation;'     &            !GSM, 06/11/08
     &        /5x,'U = ',f10.5,', A = ',f5.1,',  Z = ',f4.1,', & ', &
     &       'icntr = ',i7,'.')

! ======================================================================
      end

