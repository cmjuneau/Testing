!  ========================================================

  subroutine make_Inputs
    implicit none
    
!   This routine is meant to create equivalent outputs to be used
!   by the CEM, LAQGSM, and CEM-LAQ codes.  This routine will allow
!   for ease in trying to create files for each code.
    
    integer :: i, j, k, l, m, n, jg, num_bin, ZAID
    integer :: tot_files, cem, claq, laq, ang_width, ang_start, temp
    real    :: tstep
    
! Variables in CEM/CEMLAQ
    integer :: aproj, zproj, STIN, aproj_l, zproj_l
    integer :: atarg, ztarg
    real    :: T0, dt0, t0max, t0mev, t0gev, dt0_c, t0mev_c, t0max_c, cm0
    integer :: limc
    real    :: dteta
    integer :: mspec, mpyld, mchy, misy, mdubl, mang
    integer :: ipar1, ipar2
    integer :: ang1(10), ang2(10)
    real    :: tmin(4), tmax(4)
    real    :: dt(4)
    integer :: nevtype, npreqtype, ipisa, ihist, ityp, nh
    real    :: stp
    
! Variables in LAQGSM
    integer :: ibrems, JJJ, ISYS, INSP, iyeld, igsi
    real    :: theta(10), dtheta, dp(10), pmin
    integer :: ang(10), dang, test
    
! Character variables used in CEM/LAQGSM
    character(LEN = 100) :: pHeading(2)
    character(LEN = 17) :: laux1_s, laux2_s, name_s
    character(LEN = 16) :: clinp_s, clout_s, claux_s
    character(LEN = 15) :: cinp_s, cout_s, linp_s, lout_s, caux_s
    character(LEN = 12) :: channel1
    character(LEN = 10) :: input, output
    character(LEN =  9) :: lvl,shell
    character(LEN =  8) :: atab, mass, test_inp, temp_inp
    character(LEN =  7) :: la1_s, la2_s
    character(LEN =  6) :: cli_s, cl_s, cla_s, ZAID_s
    character(LEN =  5) :: ci_s, c_s, li_s, l_s, ca_s, stat
    character(LEN =  4) :: pname
    character(LEN =  1) :: inc_Old, custom_ang
    
    logical :: include_CEM, include_CEMLAQ, include_LAQ, get_Info
    logical :: get_constants, print_ibrems, choice, spaced, write_test
    logical :: first

    data pHeading / '     ZAID, prot - proton, neut - neutron, pipl - pi+, pimi - pi-, pize - pi0,', &
         & '     gamm - gamma with fixed energy, gamb - bremss. gamma, stop - no more calc.' /
    
    ! Default values
    include_CEM    = .false.
    include_CEMLAQ = .true.
    include_LAQ    = .true.
    get_Info       = .true.
    get_constants  = .true.
    print_ibrems   = .false.
    choice         = .false.
    spaced         = .true.
    write_test     = .true.
    first          = .true.
    
    ! File Reference No.
    cem  = 100
    claq = 200
    laq  = 300
    test = 400
    ! Creating input extensions
    if(get_Info) then
      ! Input extensions
      ci_s  = 'c.inp'
      cli_s = 'cl.inp'
      li_s  = 'l.inp'
      
      ! Output extensions
      c_s  = 'c.out'
      cl_s = 'cl.out'
      l_s  = 'l.out'
      
      ! Aux. extensions
      ca_s  = 'c.aux'
      cla_s = 'cl.aux'
      la1_s = 'l10.aux'
      la2_s = 'l11.aux'
      
      ! Files for test file
      test_inp = "store.prm"
      temp_inp = "temp.prm"
    endif
    
    if(get_constants) then
      
! CEM/CEMLAQ Data
      STIN      = 0.0d0
      dt0       = 0.0d0
      dt0_c     = -50.d0
      dteta     = 5.0d0
      mspec     = 0 ! Energy spectra
      mpyld     = 0 ! 2ndary particle yields, multiplicities, and KE (MeV)
      mchy      = 0 !
      misy      = 0 ! yields (<=3)
      mdubl     = 0 ! double differential spectra of ejectiles
      mang      = 0 ! Energy integrated angular spectra of ejectiles
      ipar1     = 1
      ipar2     = 9 ! used in conjunction with mdubl (picking which spectra we want)
      num_bin   = 13
      tstep     = 2.d0/1000.d0
      nevtype   = 66
      npreqtype = 66
      ipisa     = 1
      if(ipisa.eq.1) then
        ! Including ang(10), dang - defaults
        dang    = 5.0d0
        ang(1)  = 12.0d0
        ang(2)  = 15.0d0
        ang(3)  = 20.0d0
        ang(4)  = 35.0d0
        ang(5)  = 50.0d0
        ang(6)  = 65.0d0
        ang(7)  = 80.0d0
        ang(8)  = 90.0d0
        ang(9)  = 100.0d0
        ang(10) = 120.0d0
      endif
      ihist     = 0
      ityp      = 1
      nh        = 3
      stp       = -2.0
      
! LAQGSM Data
      stat     = "'new'"
      atab     = 'atab.dat'
      mass     = 'mass.tbl'
      lvl      = 'level.tbl'
      shell    = 'shell.tbl'
      channel1 = 'channel1.tab'
      jg       = -1
      ibrems   = 1
      JJJ      = 13
      ISYS     = 1
      INSP     = 2
      iyeld    = 3
      IGSI     = 3
    endif
      
    
! Getting preliminary information

    print *, ''
    print *, 'How many sets of inputs to make? '
    read(*,*) tot_files
    
! Creating files
    do i = 1, tot_files
      include_CEM    = .true.
      include_CEMLAQ = .true.
      include_LAQ    = .true.
      cm0 = 0.0d0
      
      ! Getting base information for inputs
      if(get_Info) then
        print *, ''
        print *, '================================'
        print *, "Name of input file (<11 char.)? "
        read(*,*) input
        
        ! Concatenating files for input name
        if(get_Info) then
          ! For input files
          cinp_s  = trim(input) // ci_s
          clinp_s = trim(input) // cli_s
          linp_s  = trim(input) // li_s
        endif
        
        ! Output file names
        print *, ''
        if(choice) then
          print *, 'Name of output file (< 11 char.)? '
          read(*,*) output
        else
          print *, 'Using input file name for ouput files...'
          output = input
        endif
        
        ! Concatenating files for output/aux name
        if(get_Info) then
          ! For output files
          cout_s  = trim(output) // c_s
          clout_s = trim(output) // cl_s
          lout_s  = trim(output) // l_s
          
          ! For aux files
          caux_s  = trim(output) // ca_s
          claux_s = trim(output) // cla_s
          laux1_s = trim(output) // la1_s
          laux2_s = trim(output) // la2_s
        endif
        
        ! Projectile
10      continue
        pname = 'unkn'
        print *, ''
        print *, "A, Z, STIN (Projectile) = ? "
        read(*,*) aproj, zproj, stin
        ! print , "Enter ZAID or PNAME below..."
        ! read(*,*) ZAID

        aproj_l = aproj
        zproj_l = zproj
        if(zproj.gt.aproj .AND. aproj.gt.0) then
          print *, "Projectile Mass Error, please re-try..."
          go to 10
        endif
        if(aproj.gt.1) then
          print *, "ERROR: Classic CEM cannot handle this projectile,"
          print *, "error making classic CEM input; continuing on..."
          include_CEM = .false.
        elseif(aproj.eq.-1.d0) then
          ! Pions
          print *, ''
          print *, 'ERROR: LAQGSM cannot handle this projectile,'
          print *, 'error making LAQGSM input; continuing on...'
          include_LAQ = .false.
        elseif(aproj.eq.0 .AND. zproj.eq.1) then
          ! Gamma - brems.
          pname = 'gamb';
          zproj_l = 0
          ibrems = 1.d0
        elseif(aproj.eq.0 .AND. zproj.eq.0) then
          ! Gamma - mono-energetic
          pname = 'gamm';
          ibrems = 0.d0
        endif
        if(aproj.le.1) then
          if(aproj.eq.1) then
            ! proton or neutron
            if(zproj.eq.1) then
              pname = 'prot';
            elseif(zproj.eq.0) then
              pname = 'neut';
            else
              print *, 'ERROR: invalid target selection!'
              go to 10
            endif
          else
            ! Other particle
            print *, 'Projectile: ', trim(pname)
          endif
        endif
        cm0 = zproj*0.9382723d0+(aproj-zproj)*0.9395656d0
        if(aproj.eq.0 .AND. cm0.lt.0.00001d0) then
          print_ibrems = .true.
        endif

        if ( pname.eq.'unkn' ) then

           ! CEM not applicable; not nucleon projectile
           ZAID = (1000*zproj) + (aproj)
           write(ZAID_s,"(i6)") ZAID

        else

           ZAID_s = trim(pname)

        endif

        print *, "ZAID is ", trim(ZAID_s)

        ! Target
20      continue
        print *, ''
        print *, 'A, Z (Target) = ? '
        read(*,*) atarg, ztarg
        if(ztarg.gt.atarg) then
          print *, "Target Mass Error, please re-try..."
          go to 20
        endif
        
        ! Energy
        print *, ''
        print *, 'Projectile Energy (GeV/A) = ? '
        print *, "(Enter a Negative No. for MeV)"
        read(*,*) T0
        if(T0.lt.0) then
          ! Units are MeV
          t0mev_c = -1*T0
          T0 = abs(T0/(aproj*1000.d0))
        else
          t0mev_c = T0*1000.d0
        endif
        t0mev = -1*T0
        t0max = t0mev*1.01
        t0max_c = 1.01*t0mev_c
        
        ! No. Simulations
        print *, ''
        print *, 'limc = ? '
        print *, "(Default min value is 500; enter negative value to over-ride)"
        read(*,*) limc
        if(limc.le.500 .AND. limc.ge.0) then
          limc = 500
        elseif(limc.lt.0) then
          limc = -1*limc
        endif
        

        ! Energy bins for angle integrated & double differential cross sections
        if(t0max_c.lt.13) then
          print *, ""
          print *, "Energy is small (relatively), please"
          print *, "manually enter data."
          print *, ""
          print *, "Enter energy step size (MeV): "
          print *, "(Negative value for GeV)"
          read(*,*) tstep
          if(tstep.le.0) then
            tstep = -1*tstep
          else
            tstep = tstep/1000.d0
          endif
          
          print *, "Enter energy bins:"
          print *, "(tmin, tmax, dt)"
          do j = 1,4
            write(*,50) j
            read(*,*) tmin(j), tmax(j), dt(j)
          enddo
          
        else
          k = t0max_c/num_bin
          k = max(k,1)
          tmin(1) = 0
          tmin(2) = (num_bin*k/13)
          tmin(3) = (3*num_bin*k/13)
          tmin(4) = (7*num_bin*k/13)
        
          do j = 1, 10
            dp(j) = k
          enddo
  
          if(spaced) then
            dt(1) = k/8
            dt(2) = k/4
            dt(3) = k/2
            dt(4) = k
          else
            do j = 1, 4
              dt(j) = k*1000.d0 ! in MeV
            enddo
          endif
        endif
        pmin    = 0
        tmax(1) = tmin(2)
        tmax(2) = tmin(3)
        tmax(3) = tmin(4)
        tmax(4) = t0max_c
        
        ! Angle information
        print *, "Custom angles (y/n)?"
        read(*,*) custom_ang
        if(custom_ang.eq.'Y' .OR. custom_ang.eq.'y') then
          print *, ""
          print *, "Please enter values that are either positive"
          print *, "or negative for that angle set (issue with angle"
          print *, "storage)"
          print *, ""
          print *, "Enter 10 angle bins (start)"
          do j = 1, 10
            write(*,50) j
            read(*,*) ang1(j), ang2(j)
          enddo
          
          if(include_LAQ) then
            ang_width = 0
            do j = 1, 10
              theta(j) = (ang1(j) + ang2(j))/2
              ang_width = ang_width + abs(theta(j)-ang2(j))
            enddo
            dtheta = ang_width/j
          endif
          
        else
          ang_start = 16
          ang_width = ang_start/2
          
          do j = 1, 10
            ang1(j) = j*ang_start - (ang_width)
            ang2(j) = j*ang_start + (ang_width)
          enddo
          
          if(include_LAQ) then
            do j = 1, 10
              theta(j) = j*ang_start
              ang(j)   = theta(j)
            enddo
            dtheta = ang_width
          endif
        endif
        
!        print *, ''
!        print *, 'The last angle used is', ang_start*10
          
        
      endif
      
      ! Creating input for CEM (old format)
      if(include_CEM) then
        cem = cem + i
        print *, ''
        print *, "Creating CEM input"
        
        ! Writing to CEM input file
        open (cem, file = cinp_s, status = 'unknown')
        write(cem,100) trim(caux_s)
        write(cem,100) trim(cout_s)
        write(cem,100) pname
        write(cem,100) trim(pHeading(1))
        write(cem,100) trim(pHeading(2))
        write(cem,200) t0mev_c
        write(cem,300) atarg
        write(cem,300) ztarg
        write(cem,400) limc
        write(cem,200) dt0_c
        write(cem,200) t0max_c
        write(cem,1300) dteta
        write(cem,300) mspec
        write(cem,300) mpyld
        write(cem,300) mchy
        write(cem,300) misy
        write(cem,300) mdubl
        write(cem,300) mang
        write(cem,600) ipar1, ipar2
        write(cem,700) ang1(1), ang2(1), ang1(2), ang2(2), ang1(3), ang2(3), &
                   & ang1(4), ang2(4), ang1(5), ang2(5), ang1(6), ang2(6), &
                   & ang1(7), ang2(7), ang1(8), ang2(8), ang1(9), ang2(9), &
                   & ang1(10), ang2(10)
        write(cem,800) tmin(1), tmax(1), dt(1), tmin(2), tmax(2), dt(2), &
                   & tmin(3), tmax(3), dt(3), tmin(4), tmax(4), dt(4)
        write(cem,300) nevtype
        write(cem,300) npreqtype
        write(cem,300) ipisa
        if(ipisa.eq.1.d0) then
          write(cem,950) theta, dtheta
        endif
        write(cem,300) ihist
        write(cem,300) ityp
        write(cem,300) nh
        write(cem,1600) aproj, zproj, pname
        write(cem,1610) atarg, ztarg
        write(cem,1620) T0, limc
        write(cem,100) 'stop'
        close(cem)
      endif
      
      ! Creating input for CEMLAQ (new format)
      if(include_CEMLAQ) then
        claq = claq + 1
        print *, ''
        print *, "Creating CEMLAQ input"
        
        ! Writing to CEM-LAQ input file
        open (claq, file = clinp_s, status = 'unknown')
        write(claq,100) trim(claux_s)
        write(claq,100) trim(clout_s)
        write(claq,100) trim(ZAID_s)
        write(claq,100) trim(pHeading(1))
        write(claq,100) trim(pHeading(2))
        ! write(claq,500) aproj, zproj, STIN
        write(claq,200) t0mev_c
        write(claq,300) atarg
        write(claq,300) ztarg
        write(claq,400) limc
        write(claq,200) dt0_c
        write(claq,200) t0max_c
        ! write(claq,600) atarg, ztarg
        ! write(claq,200) t0mev
        ! write(claq,200) dt0
        ! write(claq,200) t0max
        ! write(claq,400) limc
        write(claq,1300) dteta
        write(claq,300) mspec
        write(claq,300) mpyld
        write(claq,300) mchy
        write(claq,300) misy
        write(claq,300) mdubl
        write(claq,300) mang
        write(claq,600) ipar1, ipar2
        write(claq,700) ang1(1), ang2(1), ang1(2), ang2(2), ang1(3), ang2(3), &
                   & ang1(4), ang2(4), ang1(5), ang2(5), ang1(6), ang2(6), &
                   & ang1(7), ang2(7), ang1(8), ang2(8), ang1(9), ang2(9), &
                   & ang1(10), ang2(10)
        write(claq,800) tmin(1), tmax(1), dt(1), tmin(2), tmax(2), dt(2), &
                   & tmin(3), tmax(3), dt(3), tmin(4), tmax(4), dt(4)
        write(claq,300) nevtype
        write(claq,300) npreqtype
        write(claq,300) ipisa
        if(ipisa.eq.1.d0) then
          write(claq,950) theta, dtheta
        endif
        write(claq,300) ihist
        write(claq,300) ityp
        write(claq,300) nh
        write(claq,1600) aproj, zproj, pname
        write(claq,1610) atarg, ztarg
        write(claq,1620) T0, limc
        write(claq,100) "stop"
        ! write(claq,500) stp, stp, stp
        close(claq)
      endif
      
      ! Creating input for LAQ (standard format)
      if(include_LAQ) then
        laq = laq + 1
        print *, ''
        print *, "Creating LAQ input"
        
        ! Writing to CEM input file
        open (laq, file = linp_s, status = 'unknown')
        write(laq,100) trim(lout_s)
        write(laq,100) trim(laux1_s)
        write(laq,100) trim(laux2_s)
        write(laq,100) trim(atab)
        write(laq,100) trim(mass)
        write(laq,100) trim(lvl)
        write(laq,100) trim(shell)
        write(laq,100) trim(channel1)
        write(laq,1000) jg, trim(stat)
        write(laq,1100) aproj_l, atarg, zproj_l, ztarg, STIN, cm0, T0, limc
        if(print_ibrems) then
          write(laq,300) ibrems
        endif
        write(laq,1200) JJJ, ISYS, INSP, iyeld, nevtype, IGSI
        write(laq,950) theta, dtheta
        write(laq,950) tstep, tstep, tstep, tstep, tstep, tstep,  &
                   & tstep, tstep, tstep, tstep, pmin
        close(laq)
      endif
      
      ! Writing file names to a future input file, test.inp
      if(write_test) then
        print *, ''
        print *, "Adding infomration to ", trim(test_inp)
        
        ! Writing beginning of test file
        if(first) then
          first = .false.
          open (test, file = test_inp, status = 'unknown')
          do k = 1, 4
            write(test,1500) k
          enddo
          ! write(test,  *) ""
          ! write(test,300) tot_files
        endif
        
        ! Writing file namess into output
        write(test,100) ""
        write(test,300) 3
        if(include_CEM) then
          write(test,100) trim(cout_s)
        else
          write(test,100) trim(temp_inp)
        endif
        if(include_CEMLAQ) then
          write(test,100) trim(clout_s)
        else
          write(test,100) trim(temp_inp)
        endif
        if(include_LAQ) then
          write(test,100) trim(lout_s)
        else
          write(test,100) trim(temp_inp)
        endif
        name_s = 't' // output
        write(test,100) name_s
        
      endif
    
    enddo
    close(test)
    
    print *, ''
    return
!  ====================================
50  format ("Bin No. ", i2)
100 format (a)  
200 format (f14.5)
300 format (i3)
350 format (i6)
400 format (i9.0)
500 format (i2, ', ', i2, ', ', i2)
600 format (i3, ', ', i3)
700 format (20i5.2)
800 format (12f13.3)
900 format (11i6.1)
950 format (11f9.4)
1000 format (i3.1, ', ', a)
1100 format (5i3.1, f7.3, f10.3, i12)
1200 format (6i3.1)
1300 format (f5.2)
1400 format (f5.2, ', ', f5.2)
1500 format ("Instruction No. ", i2)
1600 format ("Projectile A, Z: ", i3, ', 'i3, "   (", A, ")")
1610 format ("Target A, Z: ", i3, ', ',i3)
1620 format ("Energy (GeV/A), limc: ", f8.4, ', ' i9.0)

  end subroutine make_Inputs

!  ========================================================
