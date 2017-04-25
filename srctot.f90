  program compare
    implicit none

!  ========================================================
!
!  Program created by CMJ, 08/16 to compare output files line by line
!  Modified by CMJ, 08/16 to allow for output reading
!
!
!  This program can do 2 things.  First, when the variable
!  'storing' is set to false, this program will compare 2
!  files and store differences into an output while displaying
!  results on the terminal window.  This option can only compare
!  2 files at this time (additional files would not be difficult
!  to add).
!  The 2nd option, when 'storing' is set to true, allows us to
!  quickly take data from CEM output files and display the results
!  side to side.  This is mostly meant to 'prep' the data, and later
!  MATLAB or ROOT can be used to extract the data and produce graphs.
!
!  These tools were created due to the frequency that this kind of
!  storing will take place.  This quickly and effectively compares
!  files, and in addition stores the output data in a way that MATLAB
!  and ROOT can easily read the data (in addition, we can too).
!
!  The value of 'read_input' is always meant to be true; it was
!  created for ease in programming.  'getNum' acts the same way.
!
!  Output(1,XXXXX) = First test file
!  Output(2,XXXXX) = Second test file
!  Output(3,XXXXX) = Third test file
!  The max number of characters in an output line is 1,000
!  The max number of lines Output can hold is 10000
!
!  ========================================================

!		----- Setting up variables -----
    integer :: MaxLines, i, j, temp, numFiles, maxFiles, numSets, num, &
             & headmult, nl
             
    real :: ddif, yieldmult
    
    character(LEN = 300) :: output
    character(LEN =  30) :: input, outFile, dummy, comp(2), instructions(4)
    character(LEN =  10) :: type
    character(LEN =   1) :: testcomp
    
    logical :: storing, store_mult, store_dds, make_own, valid, tooBig, comparing
		
!		----- Setting up common block for output -----
    common  /doublediff/  ddif(10,21,5001,361)
    common  /maxline/     MaxLines
    common  /outp/        output(10,30000), nl(3)
    common  /yieldmult/   yieldmult(10,9,9,5)
		
!		----- Setting up data blocks -----
    data valid      /.true. / ! For "while" loop; goes until EOF for input file.
    data make_own   /.false. / ! Enter own input name - requested by user
    data store_mult /.false./ ! to store multiplicities and yields
    data store_dds  /.true. / ! to store double differntial spectra
    data MaxLines   / 30000 / ! Maximum lines allowed; arbitrary really
		data numSets    /     0 / ! The total number of sets used in a given run.
    
!  ========================================================

    ! Determining program use
    ! Input file name is set by default depending on use.
    comparing = .false.
    storing = .false.
    write(*,1000)
10  continue
    read(*,*) testcomp
    if (testcomp.eq."S" .OR. testcomp.eq."s") then
      storing = .true.
      type = "store"
      input = 'store.prm'
      maxFiles = 6
    elseif(testcomp.eq."C" .OR. testcomp.eq."c") then
      comparing = .true.
      type = "compare"
      input = 'comp.prm'
      maxFiles = 2
    elseif(testcomp.eq.'M' .OR. testcomp.eq.'m') then
      write(*,1010) "make"
      call make_Inputs
      write(*,1020)
      go to 100
    else
      write(*,3000)
      go to 10
    endif
    write(*,1010) trim(type)
    
    
    ! Getting input name if not default
    if(make_own) then
      write(*,1030)
      read(*,*) input
    endif
    
    
    ! -----------------------------------------------------
    ! Main loop; reads input file andd performs other tasks
    ! -----------------------------------------------------
    
    ! Opening file
    open(1, file = input, status = 'old')
    write(*,1200) "Reading", "", trim(input)
    write(*,2000)
    write(*,   *) ""
    
    ! Input file has 4 lines of instructions; must skip over these lines
    ! Dummy lines essentially - variables not needed.
    do i = 1, 4
      read(1,*) instructions(i)
    enddo
    
    ! Previously, read the number of comparisons.  Instead, will read file
    ! as we go and exit once the "end of file" signal is obtained.
    ! 
    ! The input's "EOF" marker SHOULD be reached when reading a number or
    ! the ouput file name; otherwise there is an error with the input file
    ! and it needs fixed.
    ! 
    ! This will repeat over and over until the EOF marker is reached.
    !
    
    
    num = 100
    call setheader ! Used in Double Differential Cross Sections mostly
    do while (valid)

      write(*,*) "" ! Blank line for each set of data
      
20    continue

      ! Getting the "numFiles" variable and error checking it.
      tooBig = .false.
      if(storing) then
        ! Getting number of files to read
        read(1, *, end = 80) numFiles
        
        ! Ensuring within limits
        if(numFiles.lt.0) then
          write(*,3010) trim(type)
          
          ! Will use absolute value of the read number
          numFiles = abs(numFiles)
          write(*,3015) numFiles
          
        elseif(numFiles.eq.0) then
          
          numFiles = 3 ! A minimum default when user entered 0.
          write(*,3220) numFiles
          
        else
          
          ! Number is either too large, or in operating limits
          
        endif
          
        if(numFiles.gt.maxFiles) then
          
          write(*,3020) trim(type), maxFiles ! Error message
          temp = numFiles
          numFiles = maxFiles ! Ensures program will not crash when reading input
          tooBig = .true.
          
        else
          
          ! In operating limits
          
        endif
        
        write(*,1405) numFiles
        write(*,2000)
      else
        
        numFiles = 2
        
      endif
      
      
      ! Resetting arrays; ensures no incorrect data prints for a set
      if(storing) then
        yieldmult = 0.d0
        ddif = 0.d0
      endif
      output = "- - - - NO OUTPUT INFORMATION - - - -"
      nl = 0
      
      ! Reading input file and operating on it; iterates for each input file.
      do i = 1, numFiles
        
        ! Acquiring input file name
30      continue
        read(1, '(30A)', end = 90) input ! Name restricted to 30 characters
        if(len_trim(input).eq.0) go to 30 ! Ensures no blank lines are read.
        
        
        ! Reading Input file (the output file created by CEM or LAQGSM)
        write(*,1210) i, "Reading", "", trim(input)
        call readOutput(input, i, num)
        call incr(num)
        
        ! -----  Getting information from input file  -----
        
        if(storing) then
          ! Multiplicities/yields
          if(store_mult) then
            
            call store_yieldmult(i)
            write(*,1400)
            
          endif
          
          ! Double Differential Cross Sections
          if(store_dds) then
            
            ! Storing Double Differential Information from input file
            write(*,1410) trim(input)
            call store_doubledif(i)
            
          endif
        
        else
          
          ! Storing input names temporarily
          comp(i) = input
          
        endif
        
      enddo
      
      ! Skipping over extra file names when input formatted incorrectly.
      if(tooBig) then
        do j = numFiles + 1, temp
          read(1, '(30A)', end = 90) dummy
          write(*,1420) trim(dummy)
        enddo
      endif
      
      
      ! Output file information
40    continue

      ! Creating output file
      read(1, '(30A)', end = 100) outFile
      write(*,1200), "->Creating output", "(s) from", trim(outFile)
      
      if(storing) then
        
        ! Creating multiplicity data
        if(store_mult) then
          call make_multyields(outFile, num)
          call incr(num)
        endif
        
        ! Creating differential data
        if(store_dds) then
          call make_doubledif(outFile, num)
          call incr(num)
        endif
        
      else
        
        call compare_multiple_files(comp(1), comp(2), outFile, num)
        call incr(num)
        
      endif
      
      ! Counts the number of file sets used
      numSets = numSets + 1
      write(*,2005) ! Informing we are done with this set.
      write(*,2000) ! Displaying straight line
      
    enddo
    
    ! -----  Program is done (effectively)  -----
    ! End Signal for user
    
    go to 100
    
    
 80 continue
    ! EOF reached when reading a number; most likely due to extra blank lines.
    go to 100
    
 90 continue
 
    ! Incorrectly formatted input file, unless we are comparing
    if ( comparing ) then
      ! End of File Reached
      go to 100
    endif
    
    ! Displaying Error Message
    write(*,3200)
    
    if ( storing ) then
      ! Continue On
    else
      write(*,2000)
      stop
    endif
    
    write(*,3210) trim(input)
    
    outFile = input ! The last read item was the output name; NOT an input.
    valid = .false. ! Ensures program stops; EOF reached
    go to 40
    

100 continue

    write(*,2000)
    write(*,2010)
    write(*,2020) numSets
    write(*,2000)
    write(*,*) "" ! Blank line
    
    
! ============================================================================
! -----  Format for user dialogue (1000-1199)  -----
1000 format ("Store, Compare, or Make files (S/C/M)?")
1010 format ("Proceeding to ", A, " files")
1020 format ('Input files created')
1030 format ("Enter the input name: ")

! -----  Format for input/output file messages (1200-1399)  -----
1200 format (A, ' file', A, ' "', A, '".')
1210 format (i3, '. ', A, ' file', A, ' "', A, '".')

! -----  Format for acquiring input data (1400-1599)  -----
1400 format ('-Acquiring multiplicities from "', A, '".')
1405 format ('Storing data from the next ', i2, ' files.')
1410 format ('-Acquiring double differential cross section from "', A, '".')
1420 format ('Skipping file "', A, '"...')

! -----  Formats for program end (2000-2199)  -----
2000 format ("=====================================================")
2005 format ("Done with this set of files.")
2010 format ("Program done, see resulting output files for details.")
2020 format ("Process repeated ", i2, " times.")

! -----  Format for error messages  -----
3000 format ("ERROR: Please enter 'S' to store, 'C' to compare, or 'M' to make input files: ")
3010 format ("ERROR: A positive value for the number of files to ", A, " should be used in the input.")
3015 format ("We will try using the absolute value of the number (", i2, ") as the number of input files to be used.")
3020 format ("ERROR: Cannot ", A, " more than ", i2, " file(s).")
3200 format ("ERROR: Incorrectly formatted input file; missing input or output file name at end of file.")
3210 format ("Using last read file, '", A, "', to store resulting output.")
3220 format ("ERROR: Must use a number of files greater than 0; changing to default value (", i2, ").")

    
! ============================================================================

  end program compare
  
! ============================================================================

  subroutine incr(i)
    implicit none
    
    integer :: i
    
    i = i + 1
    
  end subroutine incr
  
! ============================================================================
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
!  ========================================================

  subroutine compare_multiple_files(input1, input2, outFile, num)
    implicit none

!		This program compares files that are input by the user.
!		One issue that has not been fixed (or been bothered to be fixed)
!		includes having to manually enter the number of lines in the file
!		as the very first line in the file.  Otherwise, this little
!		diddy does the job!


!		----- Setting up variables -----
    integer :: MaxLines, nl, nlmax, num, i, j, totmismatch, start
		integer :: mismatch(30000), test1, test2, test3
    
    character(LEN = 300) :: output, line1(30000), line2(30000)
    character(LEN =  34) :: outname
		character(LEN =  30) :: input1, input2, outFile
    character(LEN =   4) :: ext
    
		logical :: makeOutput
    
    common  /maxline/  MaxLines
    common  /outp/     output(10,30000), nl(3)
    
    
		line1 = "- - - - - NO LINE INFORMATION - - - - -"
		line2 = "- - - - - NO LINE INFORMATION - - - - -"

		
		print *, '========================================'
		print *, '---   Comparing files for ', trim(outFile), '   ---'
		print *, '========================================'
		print *, "---   Files being compared   ---"
		print *, "   1) ", trim(input1)
		print *, "   2) ", trim(input2)
		print *, '========================================'
		print *, " "	

    ! Transferring variables from output variable to line1 variable
		do j = 1, nl(1)
      line1(j) = output(1,j)
    enddo
    
    ! Transferring variables from output variable to line2 variable
    do j = 1, nl(2)
      line2(j) = output(2,j)
    enddo
    
	
!  comparing lines in files
		totmismatch = 0
		mismatch = 0
		start = 31
    nlmax= max(nl(1), nl(2))
	  do j = start, nlmax
			if( line1(j).ne.line2(j) .AND. (j.ne.39 .AND. j.ne.51 .AND. j.ne.88)) then
				print *, 'Mismatch at line', j
				mismatch(j) = 1
				totmismatch = totmismatch + 1
			endif
		enddo
		
		print *, ''
		print *, '--------------------'
		if( totmismatch.eq.0) then
			print *, '--->No Mismatches!'
		else
			print *, 'Total mismatches: ', totmismatch
		endif
		print *, '--------------------'
		print *, ''
    call sleep(1)
			
!		----- Formatting output file -----
    makeOutput = .true.
    if (makeOutput) then
      ext = ".txt"
      outname = trim(outFile) // ext
			open (num, file = outname, status = 'new')
      
			write(num,*) 'Comparing the following files: '
			write(num,*) "     1) ", trim(input1), "		('New' Version)"
			write(num,*) "     2) ", trim(input2), "		('Old' Version)"
			write(num,*) ''		
			write(num,*)'Lines in ', trim(input1), ': ', nl(1)
			write(num,*)'Lines in ', trim(input2), ': ', nl(2)
			write(num,*) ''
			if ( totmismatch.eq.0 ) then
				write(num,*) '--------------------'
				write(num,*) '---> No mismatches'
				write(num,*) '--------------------'
				write(num,*) ''
				write(num,*) "----------"
				write(num,*) '| Proof: |'
				write(num,*) "----------"
				write(num,*) ' '
				test1 = 252
				test2 = 550
				test3 = 1105
				write(num,*) 'Line ', test1
				write(num,*) line1(test1)
				write(num,*) line2(test1)
				write(num,*) ' '

				write(num,*) 'Line ', test2
				write(num,*) line1(test2)
				write(num,*) line2(test2)
				write(num,*) ' '

				write(num,*) 'Line ', test3
				write(num,*) line1(test3)
				write(num,*) line2(test3)
				write(num,*) ' '
			else
				write(num,*) ''
				write(num,*) '--------------------'
				write(num,*) 'Total mismatches: ', totmismatch
				write(num,*) '--------------------'
				write(num,*) ''
				write(num,*) 'Differences occurred at the following lines:'
				write(num,*) ''
				do j = start, nlmax
					if ( mismatch(j).eq.1 ) then
						write(num,*) 'Mismatch at line', j
						write(num,*) '   -> From ', trim(input1), ':		', trim(line1(j))
						write(num,*) '   -> From ', trim(input2), ':		', trim(line2(j))
						write(num,*) ''
					endif
				enddo
			endif
			close(num)	
			print *, "--  Output file ", trim(outname), " created with results  --"
			print *, " "
		endif
    
		
		return
	end subroutine compare_multiple_files
	 
!  ========================================================!  ========================================================

  subroutine store_doubledif(i)
    implicit none
    
!  This routine is meant to grab the double differential data 
!  from the ouput file.  This information is stored in
!  ouput(i,LINES).  Data will be stored such that:
!
!  ddif(file no., No. Ejectiles, T (MeV), Theta)
!
!  file = file no. referenced in main program "do loop"
!       = 1, 2, or 3
!
!  No. Ejectiles (1-21)
!      1 = n
!      2 = p
!      3 = d
!      4 = t
!      5 = He-3
!      6 = He-4
!      etc...
!
!  T (MeV) - Max is ~ 2500 MeV
!    -when printing, don't print 0 values (0 for all files)
!    -will have MATLAB graphing exlcude 0 values w/ x(x == 0) NaN; plot(x) ...
!
!  Theta (Degrees) - seems max is ~ 120 degrees - allocated 360 degrees
!    -for Angle Integrated Spectra, use theta = 361
!
!		----- Setting up variables -----
    integer :: i, ej, T, theta, headmult, MaxLines, nl, j, k, l, eject, indx, ddang(10)
    integer :: tindx, Ti, Tf, corr, ejtype, ejindx, ejecta(11), energi, thetai, stat
    integer :: offset, N, dataSet, partNum(10), tempNum
    real :: ddif, ang(360), tempang, angint, Tir, Tfr, temp, thetar, lang(11),energr
    real :: dubdiff, dubdifftemp, sum, dtemp, d_ang_int, pi, correction
    real :: z, y, x

    character(LEN =  300) :: header, mischeader, output
    character(LEN =  300) :: temp0, temp1, temp2, temp3, temp_1
    character(LEN =   10) :: Tis, Tfs, angs(360), angints, temps, temps2, thetas, langs(11), energs
    character(LEN =    5) :: ejectileL
    character(LEN =    4) :: ejectile
    
    logical :: storedatac, storedatal, storeangl, store_enrg, change_3, storeangc, store_nucleons
    
!		----- Setting up common block for output -----
    common  /header/     header(100), mischeader(50), headmult
    common  /outp/       output(10,30000), nl(3)
    common  /maxline/    MaxLines
    common /doublediff/  ddif(10,21,5001,361)
    
    data pi /3.141592653584/
    
!  Setting initial values to 0
    do ej = 1, 21
      do T = 1, 5001
        do theta = 1,361
          ddif(i,ej,T,theta) = 0.d0
        enddo
      enddo
    enddo
    ejectile = " NaN "
    storedatac = .false.
    storedatal = .false.
    storeangl  = .false.
    store_enrg = .false.
    storeangc  = .false.
    store_nucleons = .false.
    change_3   = .false. ! Last modified 10/2016; if true, yields NaN values for 3rd data aset.
    corr = 0
    dataSet = 0;
    
!  Main loop
!    print *, " "
    do j = 1, nl(i)
      temp3 = output(i,j-3)
      temp2 = output(i,j-2)
      temp1 = output(i,j-1)
      temp0 = output(i,j  )
      temp_1= output(i,j+1)
      
! Checking if at CEM section start
      if(temp3(1:83).eq.header(30)) then
        storedatac = .true.
        ejectile = trim(temp3(84:88))
!        print *, "- - Ejectile: ", trim(ejectile), " - -"
        Ti = 0
        Tf = 0
        
! Storing ejectile type
        eject = 0
        if (ejectile.eq."   n") then
          eject = 1
        elseif(ejectile.eq."   p") then
          eject = 2
        elseif(ejectile.eq."   d") then
          eject = 3
        elseif(ejectile.eq."   t") then
          eject = 4
        elseif(ejectile.eq." He3") then
          eject = 5
        elseif(ejectile.eq." He4") then
          eject = 6
        elseif(ejectile.eq." He6") then
          eject = 7
        elseif(ejectile.eq." Li6") then
          eject = 8
        elseif(ejectile.eq." Li7") then
          eject = 9
        elseif(ejectile.eq." Li8") then
          eject = 10
        elseif(ejectile.eq." Li9") then
          eject = 11
        elseif(ejectile.eq." Be7") then
          eject = 12
        elseif(ejectile.eq." Be9") then
          eject = 13
        elseif(ejectile.eq."Be10") then
          eject = 14
        elseif(ejectile.eq." B10") then
          eject = 15
        elseif(ejectile.eq." B11") then
          eject = 16
        elseif(ejectile.eq." B12") then
          eject = 17
        elseif(ejectile.eq." C11") then
          eject = 18
        elseif(ejectile.eq." C12") then
          eject = 19
        elseif(ejectile.eq." C13") then
          eject = 20
        elseif(ejectile.eq." C14") then
          eject = 21
        elseif(ejectile.eq." NaN")  then
          !  Intermiediate reading; should never reach here
          print *, "Mis-reading in store_doubledif, stopping"
          print *, "Issue on line:", j
          stop
        else
!          print *, "Ejectile not applicable for comparison"
!          print *,
          storedatac = .false.
        endif
        
        ! Getting angles for data
        do indx = 1, 10
          k = 11*(indx-1)+17
          temps = temp1(k:k+4)
          read(temps,*) tempang
          ddang(indx) = nint(tempang)
        enddo
      else if(temp3(1:90).eq.header(35)) then
        ! print *, 'Acquiring angle integrated spectra'
        storeangc = .true.
        dataSet = dataSet + 1
        
        ! Gets ejectile numbers based on table formatting, as of 2/2017 version of CEM03.03.
        if ( dataSet.eq.1 ) then
          do k = 1, 10
            partNum(k) = k
          enddo
        elseif ( dataSet.eq.2 ) then
          do k = 11, 20
            l = k
            if(l.eq.15) l = 0 ! B9 not included in program
            if(l.gt.15) l = l - 1
            partNum(k-10) = l
          enddo
        elseif ( dataSet.eq.3 ) then
          do k = 21, 30
            l = k - 1 ! the "-1" corrects for lack of B9.
            if(l.gt.21) l = 0
            partNum(k-20) = l
          enddo
        else
          print *, 'Unknown Error - Stopping.'
          stop
        endif
      endif

! Checking if at CEM section end
      if(temp0.eq. " " .AND. temp1(1:13).eq.mischeader(42)) then
        if(storedatac) then
!          print *, "Done with ejectile: ", trim(ejectile)
!          print *, " "
        endif
        storedatac = .false.
        ejectile = "NaN"
        ddang = 0
        eject = 0
        corr  = 0
      elseif(temp0(1:12).eq.header(36) .OR. temp0(1:11).eq.header(37)) then
        ! print *, '-Done with Angle Integrated Spectra'
        storeangc = .false.
      endif

! Getting values from the CEM string
      if(storedatac) then
        
! Getting data as a string
        Tis = temp0(1:5 )
        Tfs = temp0(7:12)
        if(Tis.eq." ener     ") then
          corr = 1
          Ti = 1
          Tf = 5001
          go to 10 ! Not storing this information - this is the energy
          ! integrated spectra given at the end of each table.
        else
          corr = 0
          read (Tis    ,*) Tir
          read (Tfs    ,*) Tfr
          Ti = nint(Tir)
          Tf = nint(Tfr)
          if(Tf.gt.5000) then
            print *, "ERROR: Final Kinetic Energy found to be greater"
            print *, "than 5 GeV for a ", trim(ejectile), "ejectile at line", j
            print *, "Please change the limit of Tf in store_doubledif and "
            print *, "make_doubledif subroutines."
!            stop
            go to 10
          endif
        endif
        
        ! Gets spectra values for specific angles
        do indx = 1, 10
          tindx = 11*(indx-1) + 13 + corr
          angs(ddang(indx)) = temp0(tindx:tindx+10)
        enddo
        angints = temp0(123+corr:133+corr)
        
! Getting data as a real
        do indx = 1, 10
          tindx = ddang(indx) ! this is the theta for the data
          if(tindx.gt.360 .OR. tindx.lt.0) then
            print *, "ERROR: angle out of bounds!"
            print *, "Line No.", j
            stop
          endif
          temps = angs(tindx)
          read (temps,*) ang(tindx)
          ddif(i,eject,Tf,tindx) = ang(tindx) + ddif(i,eject,Tf,tindx)
        enddo
        read (angints,*) ddif(i,eject,Tf,361) ! Angle integrated data
        
      elseif(storeangc) then
        ! Reading energy integrated values.  This data is presented in 3 sets - perhaps we can count them.
        ! The 1st data set is differnt from the other 2.  The first 2 also have headers describing the particle,
        ! while the 3rd set mostly describes particles by the Z.
        ! print *, 'Line No.', j
        
        ! Need to now read in values from file; we are at the correct location.
        ! Getting average angle
        temps  = temp0(1:4)
        temps2 = temp0(6:10)
        
        ! Converting to number - storing data in right edge of angle bin
        read(temps2,*) tempang
        
        ! Now we have the angle at which we have data.  The energy value is 5001.
        Tf = 5001
        
        ! Reading in data from file.  Need location of data first, then must 
        do k = 1,10
          
          ! Storing data as number in ang(k)
          temps = temp0(k*11:k*11+9)
          read(temps,*) ang(k)
          
          ! Need to store data in ddif array now.  If particle number is 0, particle not allowed.
          tempNum = partNum(k)
          if ( tempNum.gt.0 ) then
            ddif(i,tempNum,Tf,tempang) = ang(k)
            ! print *, k, ddif(i,tempNum,Tf,tempang)
          else
            ! Particle not in list
          endif
        enddo
        
      endif
      
! Checking if at LAQGSM section start
      if(temp2(1:62).eq.header(31) .OR. temp2.eq.header(32) .OR. temp2.eq.header(33)) then
        if(temp1(1:108).eq.mischeader(43)) then
          if(temp0.ne." ") then
            storedatal = .true.
            offset = 15
          else
!            print *, "- - No values for the Angle at output line - -", j
!            print *, " "
            storedatal = .false.
          endif
        elseif(temp1(1:108).eq.mischeader(44)) then
          if(temp0.ne." ") then
            storeangl = .true.
            offset = 15
!            print *, "- - Reading Angle Integrated Spectra - -"
          else
            storeangl = .false.
!            print *, "- - No values for Angle Integrated spectra - -"
          endif
        elseif(temp1(1:106).eq.mischeader(45)) then
          if(temp0.ne." ") then
            store_enrg = .true.
            offset = 13
!            print *, "- - Reading Energy Integrated Spectra - -"
          else
            store_enrg = .false.
!            print *, "- - No values for energy integrated spectra - -"
          endif
        elseif(temp1(1:97).eq.mischeader(46)) then
          store_nucleons = .true.
          offset = 13
!          print *, "- - Storing Values for Nucleons of LAQGSM - -"
        endif
        
        if(storedatal .OR. store_nucleons) then
          thetas = trim(temp2(63:68))
          read(thetas,*) thetar
          thetai = nint(thetar)
!          print *, "- - Angle (degrees): ", trim(thetas), " - -"
        endif
        
! Getting particles for data storage
        if(storedatal .OR. storeangl .OR. store_enrg .OR. store_nucleons) then
          do ejtype = 1, 11
            
            ! Making a variable offset for nucleons
            if ( store_nucleons ) then
              if ( ejtype.eq.2 ) then
                offset = 14
              elseif ( ejtype.eq.3 ) then
                offset = 13
              elseif ( ejtype.eq.4 ) then
                offset = 14
              elseif ( ejtype.eq.5 ) then
                offset = 13
              elseif ( ejtype.eq.6 ) then
                offset = 14
              elseif ( ejtype.eq.7 ) then
                offset = 17
              elseif ( ejtype.eq.8 ) then
                offset = 17
              elseif ( ejtype.eq.9 ) then
                offset = 17
              elseif ( ejtype.eq.10 ) then
                offset = 17
              elseif ( ejtype.eq.11 ) then
                offset = 15
              endif
            endif
            
            if(ejtype.eq.11) then
              ejindx = 108 + offset
            else
              ejindx = offset + (ejtype-1)*11
            endif
            
            ejectileL = temp1(ejindx:ejindx+5)
        
! Storing ejectile type with a do loop
            eject = 0
            if (ejectileL.eq."NEUTR") then
              eject = 1
            elseif(ejectileL.eq."PROTO") then
              eject = 2
            elseif(ejectileL.eq."DEUTR") then
              eject = 3
            elseif(ejectileL.eq."TRITO") then
              eject = 4
            elseif(ejectileL.eq."HE-3 ") then
              eject = 5
            elseif(ejectileL.eq."HE-4 ") then
              eject = 6
            elseif(ejectileL.eq."He6  ") then
              eject = 7
            elseif(ejectileL.eq."Li6  ") then
              eject = 8
            elseif(ejectileL.eq."Li7  ") then
              eject = 9
            elseif(ejectileL.eq."Li8  ") then
              eject = 10
            elseif(ejectileL.eq."Li9  ") then
              eject = 11
            elseif(ejectileL.eq."Be7  ") then
              eject = 12
            elseif(ejectileL.eq."Be9  ") then
              eject = 13
            elseif(ejectileL.eq."Be10 ") then
              eject = 14
            elseif(ejectileL.eq."  B10") then
              eject = 15
            elseif(ejectileL.eq."  B11") then
              eject = 16
            elseif(ejectileL.eq."  B12") then
              eject = 17
            elseif(ejectileL.eq."C11  ") then
              eject = 18
            elseif(ejectileL.eq."  C12") then
              eject = 19
            elseif(ejectileL.eq."  C13") then
              eject = 20
            elseif(ejectileL.eq."  C14") then
              eject = 21
            elseif(ejectileL.eq." NaN ") then
              !  Intermiediate reading; should never reach here
              print *, "Mis-reading in store_doubledif, stopping"
              print *, "Issue on line:", j
              stop
            else
!              print *, "Ejectile ", trim(ejectileL), " not applicable for comparison"
            endif
            ejecta(ejtype) = eject
          enddo
          
        endif
        
      endif
      
! Checking if at LAQGSM section end
      if(temp0.eq." " .OR. temp0.eq."") then
        
!        write(*,*) "Checking for end of file at line ", j
        if(storedatal .AND. (temp_1(1:62).eq.header(31) .OR. temp_1.eq.header(32))) then
!          print *, "- - Done with Angle: ", trim(thetas), " - -"
!          print *, " "
          storedatal = .false.
        endif
        if(storeangl  .AND. temp_1.eq.header(33)) then
!          print *, "- - Done with Angle Integrated Spectra - -"
!          print *, " "
          storeangl  = .false.
        endif
        if(store_enrg .AND. temp_1.eq." ") then
!          print *, "- - Done with Energy Integrated Spectra - -"
!          print *, " "
          store_enrg = .false.
        endif
        if(store_nucleons) then
!          print *, "- - Done with Angle: ", trim(thetas), " - -"
!          print *, " "
          store_nucleons = .false.
        endif
        ejectileL = " NaN "
        ddang = 0
        eject = 0
        corr  = 0
      endif
      
! Getting values from the LAQGSM string
      if(storedatal .OR. storeangl .OR. store_enrg .OR. store_nucleons) then
!        print *, j, trim(temp0)
        ! Storing data as string
        if(storedatal .OR. storeangl .OR. store_nucleons) then
          energs = temp0(1:10)
          
!          write(*,*) "Energy (GeV): ", trim(energs)
          read(energs,*,iostat=stat) energr  ! in GeV
          if(stat.ne.0) then
            print *, "Failure converting energy (string) to int on output line", j
            print *, "Line Info: ", trim(temp0)
            print *, "Value read is:", energs
            stop
          endif
          energr = energr*1000.d0 ! in MeV
          energi = nint(energr) ! Converting to an integer for array storage
! Ensuring energy < 500 MeV
          if(energi.gt.5001) then
            print *, "ERROR: Final Kinetic Energy found to be greater"
            print *, "than 5 GeV for a ", trim(ejectile), "ejectile at line", j
            print *, "Please change the limit of Tf in store_doubledif and "
            print *, "make_doubledif subroutines."
!            stop
            go to 10
          endif
          if(storeangl) then
            thetai = 361
          endif
        elseif(store_enrg) then
          energi = 5001
          thetas = temp0(1:10)
          read(thetas,*,iostat=stat) thetar
          if(stat.ne.0) then
            print *, "Failure converting angle (string) to int on output line", j
            stop
          endif
          thetai = nint(thetar)
        else
          write(*,*) "Incorrect Specifier in 'doublediff.f90', line 488."
          stop
        endif
        do indx = 1, 11
          tindx = 11 * indx
          langs(indx) = temp0(tindx:tindx + 10)
        enddo
        
! Converting data into real
        do indx = 1, 11
          tindx = 11 * indx
          temps = langs(indx)
          read(temps,*,iostat=stat) lang(indx)
! Converting from (mb/GeV/sr) to (mb/MeV/sr) the match CEM
          lang(indx) = lang(indx)/1000.d0
          if(isnan(lang(indx))) then
            print *, "Failure reading cross section at line", j, ",", indx, " cross sections in..."
            print *, "Value read is: ", trim(temp0(tindx:tindx + 10))
            print *, "Value obtained is:", lang(indx)
            stop
          endif
          if(stat.ne.0) then
            print *, "Failure converting dS/(dTdO) value (string) to real on output line", j
            stop
          endif
          
        enddo
        
! Storing data into ddif
        do indx = 1, 11
          eject = ejecta(indx)
          dubdiff = lang(indx)
          
!          write(*,*) "Eject No: ", eject, ", Angle: ", thetai, ", Energy: ", energi, ", Value: ", dubdiff
          
          if(eject.ne.0) then
            dubdifftemp = ddif(i,eject,energi,thetai)
            ddif(i,eject,energi,thetai) = dubdiff + dubdifftemp
          endif
        enddo
      endif
      
   10 continue
    enddo
          
    return
    
 100  format("NaN value created for ejectile no. ", i2,", energy of ", i4," [MeV], and angle ", i3)
  end subroutine store_doubledif
  
!  ========================================================

  subroutine make_doubledif(outname, num)
    implicit none
    
!
! This program will output the double differential data
! into an approporiately named output file.  This routine
! is one of the most important to use for data analysis.
! The program will store this data based on ejectile type
! and angle of emission.
!
		
!		----- Setting up variables -----
    integer :: i, ej, Tf, theta, value
    integer :: num
    real :: ddif, temp(6)
    
    character(LEN = 37) :: output(3)
    character(LEN = 30) :: outname
    character(LEN =  6) :: ext(3)
    
    logical setup, doubledif, angint, energint
    
		
!		----- Setting up common block for output -----
    common /doublediff/  ddif(10,21,5001,361)

    setup = .true.
    
! Initial stuff
    if(setup) then

! Sections
      doubledif = .true. ! All double differential data
      angint    = .false. ! Angle integrated spectra
      energint  = .false. ! Energy integrated spectra
      
! Output file extensions
      num = num + 999
      ext(1) = "dd.txt"
      ext(2) = "ai.txt" ! Angle Integrated spectra
      ext(3) = "ei.txt" ! Energy intgrated spectra
      do i = 1, 3
        output(i) = trim(outname) // ext(i)
      enddo
    endif

! All double differential data
    if(doubledif) then
      num = num + 1
      print *, "Creating output file ", trim(output(1)), "."
      open(num, file = trim(output(1)), status = 'new')

! Writing to file
!      write(num,  *) "Double Differential Data"
!      write(num,100)
!      write(num,  *) " "
! writing double differential spectra
      do ej = 1, 21
        do theta = 1, 361
          do Tf = 1, 5001
        ! do Tf = 1, 5001
         !  do theta = 1, 361
            value = 0;
            do i = 1, 6
              temp(i) = ddif(i,ej,Tf,theta)
              if(temp(i).gt.0) then
                value = 1;
              endif
            enddo
            if(value.eq.1) then
!              write(*,150) ej, Tf, theta, temp1, temp2, temp3
              write(num,150) ej, Tf, theta, temp(1), temp(2), temp(3), temp(4), temp(5), temp(6)
            endif
          enddo
        enddo
      enddo
      close(num)
    endif
    
! Angle Integrated Spectra
    if(angint) then
      num = num + 1
      print *, "Creating output file ", trim(output(2)), "."
      open(num, file = trim(output(2)), status = 'new')

! Writing to file
!      write(num,  *) "Angel Integrated Data for Double Differential Cross Sections"
!      write(num,200)
!      write(num,  *) " "
! writing angle integrated spectra spectra
      theta = 361
      do ej = 1, 21
        do Tf = 1, 5000
          value = 0;
          do i = 1, 6
            temp(i) = ddif(i,ej,Tf,theta)
              if(temp(i).gt.0) then
              value = 1
            endif
          enddo
          if(value.gt.0) then
            write(num,250) ej, Tf, theta, temp(1), temp(2), temp(3), temp(4), temp(5), temp(6)
          endif
        enddo
      enddo
      close(num)
    endif
    
! Energy Integrated Spectra
    if(energint) then
      num = num + 1
      print *, "Creating output file ", trim(output(3)), "."
      open(num, file = trim(output(3)), status = 'new')

! Writing to file
!      write(num,  *) "Energy Integrated Data for Double Differential Cross Sections"
!      write(num,300)
!      write(num,  *) " "
! writing angle integrated spectra spectra
      Tf = 5001
      do ej = 1, 21
        do theta = 1, 360
          value = 0
          do i = 1, 6
            temp(i) = ddif(i,ej,Tf,theta)
            if(temp(i).gt.0) then
              value = 1;
            endif
          enddo
          if(value.gt.0) then
            write(num,350) ej, Tf, theta, temp(1), temp(2), temp(3), temp(4), temp(5), temp(6)
          endif
        enddo
      enddo
      close(num)
    endif
    
    return

!  ========================================================
100 format("Ejectile   E (MeV)   Theta (Deg)   :         Files - up to 6")
150 format(i7.1, " ", i8.1, " ", i10.1, "        :  ", 6es17.7)
200 format("Ejectile   E (MeV)    :         Files - up to 6")
250 format(i7.1, " ", i9.1, "     :  ", 6es17.7)
300 format("Ejectile    Theta (Deg)   :         Files - up to 6")
350 format(i7.1, " ", " ", i9.1, "        :  ", 6es17.7)

  end subroutine make_doubledif

!  ========================================================!  ========================================================

  subroutine store_yieldmult(i)
    implicit none
!
!  This routine is meant to store multiplicity data from
!  the output.  The goal of this routine is to quickly
!  and efficiently store all the data into one array.
!  The data is stored in this manner:
!
!   yieldmult(1-3 , 1-9     , 1-9 ,  1-5)
!   yieldmult(file, particle, type, [-multiplicity -mult. error -yield -yield error -tke)
!
!  file = file no. referenced in main program "do loop"
!       = 1, 2, or 3
!  Type =
!         1 for T  (all production mechanisms)
!         2 for C  (cascade)
!         3 for P  (pre-equilibrium)
!         4 for Sp (spallation)
!         5 for Pf (pre-fision)
!         6 for F  (fission)
!         7 for E  (total evaporation)
!         8 for Co (coalescence)
!         9 for pions
!  Particle = 
!         1 for n
!         2 for p
!         3 for d
!         4 for  t
!         5 for He-3
!         6 for He-4
!         7 for pi-
!         8 for pi0
!         9 for pi+
    
!		----- Setting up variables -----
    integer :: i, j, k, l
    integer :: nl, headmult, MaxLines, tpe, part
    real ::    yieldmult
    
    character(LEN = 300) :: output, header, mischeader, temp
    character(LEN = 7) :: mults, multerrs, yields, yielderrs, tkes
    character(LEN = 5) :: type, prtcl, uncert1, uncert2
    
    logical :: storedata
    
!		----- Setting up common block for output -----
    common  /outp/       output(10,30000), nl(3)
    common  /header/     header(100), mischeader(50), headmult
    common  /maxline/    MaxLines
    common  /yieldmult/  yieldmult(10,9,9,5)
    
! setting all values to 0
    do j = 1, 9
      do k = 1, 9
        do l = 1, 5
          yieldmult(i,j,k,l) = 0.d0
        enddo
      enddo
    enddo
    storedata = .false.
        
    do j = 1, nl(i)
      
! Checking if at section start
      if(output(i,j-2).eq.header(1)) then
        if(output(i,j-1).eq.mischeader(1)) then
				   storedata = .true.
        endif
      endif
      
! Checking if at section end
      if(output(i,j-1).eq.mischeader(1)) then
        if(output(i,j).eq." ") then
          go to 10
        endif
      endif
            
! Getting values from the string
      if(storedata) then
        part = 0
        tpe = 0
        if(output(i,j).eq.mischeader(1) .OR. output(i,j).eq.' ') then
! This indicates a new particle section, all *'s
          go to 5
        endif
        temp = output(i,j)
!  Storing data into temp strings
        type = temp(2:3)
        prtcl = temp(4:7)
        mults = temp(10:18)
        multerrs = temp(21:29)
        yields = temp(32:41)
        yielderrs = temp(45:53)
        tkes = temp(54:62)
        
! 7 different types
        if(type.eq."T ") then
          tpe = 1.d0
        elseif(type.eq."C ") then
          tpe = 2.d0
        elseif(type.eq."P ") then
          tpe = 3.d0
        elseif(type.eq."Sp") then
          tpe = 4.d0
        elseif(type.eq."Pf") then
          tpe = 5.d0
        elseif(type.eq."F ") then
          tpe = 6.d0
        elseif(type.eq."E ") then
          tpe = 7.d0
        elseif(type.eq."Co") then
          tpe = 8.d0
        elseif(type.eq."pi") then
          tpe = 9.d0
          prtcl = temp(2:5)
        else
          print *, "Error Reading Type of Info in multiplicity/yield reading"
          stop
        endif
        
! 6 different particles
        if(prtcl.eq." n ") then
          part = 1.d0
        elseif(prtcl.eq." p ") then
          part = 2.d0
        elseif(prtcl.eq." d ") then
          part = 3.d0
        elseif(prtcl.eq." t ") then
          part = 4.d0
        elseif(prtcl.eq."He3") then
          part = 5.d0
        elseif(prtcl.eq."He4") then
          part = 6.d0
        elseif(prtcl.eq."pi-") then
          part = 7.d0
        elseif(prtcl.eq."pi0") then
          part = 8.d0
        elseif(prtcl.eq."pi+") then
          part = 9.d0   
        else
          print *, "Error reading particle name in multiplicity/yield reading"
        endif
        
! Storing data from line into yieldmult array
        read(mults,*)      yieldmult(i,part,tpe,1) ! Multiplicity
        read(multerrs,*)   yieldmult(i,part,tpe,2)
        read(yields,*)     yieldmult(i,part,tpe,3) ! Yield
        read(yielderrs,*)  yieldmult(i,part,tpe,4)
        if(tkes.eq.mischeader(41)) then
          yieldmult(i,part,tpe,5) = 999999999
        else
          read(tkes,*)       yieldmult(i,part,tpe,5) ! Average KE
        endif
      endif
5     continue
    enddo
10  continue
    
    return
  end subroutine store_yieldmult

!  ========================================================

  subroutine make_multyields(outname, num)
    implicit none
    
!
! This program will output our multiplicity and yields
! in a new file so MATLAB or ROOT can create histograms
! using this data.  Output name is passed into code and
! given an extension "mlt.txt", "yld.txt", and "TKE.txt"
!
! Format of the output is meant to be as:
!
!  particle, type  :  value1 value2 value3 -1 error1 error2 error3
!
! This simply outputs the ENTIRE yieldmult array into a txt file.
! This is the easiest way to output all of the information and 
! allow MATLAB/ROOT to read all the information before creating
! a graph of the data.
!
		
!		----- Setting up variables -----
    integer :: i, j, k, l
    integer :: num
    real :: yieldmult
    
    character(LEN = 37) :: output(3)
    character(LEN = 30) :: outname
    character(LEN =  7) :: ext(3)
    
    logical setup, multiplicity, yields, tkes
    
		
!		----- Setting up common block for output -----
    common  /yieldmult/  yieldmult(10,9,9,5)

    setup = .true.
    
! Initial stuff
    if(setup) then

! Sections
      multiplicity = .true.
      yields = .true.
      tkes = .true.
      
! Output file extensions
      num = num + 999
      ext(1) = "mlt.txt"
      ext(2) = "yld.txt"
      ext(3) = "TKE.txt"
      do i = 1, 3
        output(i) = trim(outname) // ext(i)
      enddo
    endif

! Multiplicities
    if(multiplicity) then
      num = num + 1
      print *, "Creating output file ", trim(output(1)), "."
      open(num, file = trim(output(1)), status = 'new')

! Writing to file
!      write(num,  *) "Multiplicity Data:"
!      write(num,100)
!      write(num,  *) " "
! writing multiplicities
      do j = 1, 9
        do k = 1, 9
          write(num,150) j, k, yieldmult(1,j,k,1), yieldmult(2,j,k,1), yieldmult(3,j,k,1), &
                          &  yieldmult(1,j,k,2), yieldmult(2,j,k,2), yieldmult(3,j,k,2)   
        enddo
      enddo
      close(num)
    endif
    
! Yields
    if(yields) then
      num = num + 1
      print *, "Creating output file ", trim(output(2)), "."
      open(num, file = trim(output(2)), status = 'new')

! Writing to file
!      write(num,  *) "Yield Data:"
!      write(num,100)
!      write(num,  *) " "
! writing yields
      do j = 1, 9
        do k = 1, 9
          write(num,150) j, k, yieldmult(1,j,k,3), yieldmult(2,j,k,3), yieldmult(3,j,k,3), &
                          &  yieldmult(1,j,k,4), yieldmult(2,j,k,4), yieldmult(3,j,k,4)   
        enddo
      enddo
      close(num)
    endif
    
! TKE
    if(tkes) then
      num = num + 1
      print *, "Creating output file ", trim(output(3)), "."
      open(num, file = trim(output(3)), status = 'new')

! Writing to file
!      write(num,  *) "Average Kinetic Energy Data:"
!      write(num,200)
!      write(num,  *) " "
! writing tke
      do j = 1, 9
        do k = 1, 9
          write(num,250) j, k, yieldmult(1,j,k,5), yieldmult(2,j,k,5), yieldmult(3,j,k,5)
        enddo
      enddo
      close(num)
    endif
    
    print *, " "
    return

!  ========================================================
100 format("Part./  Type.  :    Test 1    Test 2    Test 3   ::  &
                            & Error 1    Error 2    Error 3")
150 format(i6.1, i6.1, "   :  ", f9.4, "  ", f9.4, "  ", f9.4, "   :: ", &
                              &  f9.4, "  ", f9.4, "  ", f9.4)
200 format("Part./  Type.  :            Test 1           Test 2           Test 3")
250 format(i6.1, i6.1, "   :  ", f15.3, "  ", f15.3, "  ", f15.3)

  end subroutine make_multyields

!  ========================================================!  ========================================================

  subroutine readOutput(files, i, numi)
    implicit none
    
!    
!   This program will read the output and store
!   the information to be used later
!

!		----- Setting up variables -----
    integer :: numi, i, j, nli, nl, MaxLines, headmult, filenum
    
    character(LEN = 300) :: output, header, mischeader
    character(LEN =  30) ::  files
		
!		----- Setting up common block for output -----
    common  /outp/     output(10,30000), nl(3)
    common  /header/   header(100), mischeader(50), headmult
    common  /maxline/  MaxLines

! Various Initializations
      nli = 1 ! Number of lines
      
!  reading from first file listed
			open (numi, file = trim(files), status ='old')
	    do j = 1, MaxLines
				read(numi, '(300A)', end = 5) output(i,j)
				if (output(i,j).eq."" .OR. output(i,j).eq.mischeader(40)) then 
          output(i,j) = " "
        endif
        nli = nli + 1
		  enddo
      
	5   continue
      close(numi)
      nl(i) = nli-1

    return
  end subroutine readOutput

!  ========================================================!  ========================================================

  subroutine setheader
    implicit none
    
!
!  This subroutine will create a list of section headers
!  used in CEM.  These headers will indicate to the code
!  when to read output values and store appropriately.

    integer :: headmult
    
    character(LEN = 300) :: header, mischeader
    
    logical :: makeheader, makemisc
    
    
    common  /header/ header(100), mischeader(50), headmult
! Initialize several variables
		header = "- - - - - NO HEADER INFORMATION - - - - -"
    mischeader = "- - - - - NO MISC. INFORMATION - - - - -"
    makeheader = .true.
    makemisc = .true.
    
! Header information
    if(makeheader) then
      header(1)  = " Part.     Multiplicities           Yields [mb]     <TKE> [MeV]"
      header(2)  = " Yields of different channels (with > 1 mb):"
      header(3)  = " Yields of several summed channels:"
      header(4)  = " *************** Nuclide yields [mb]  (zero values suppressed) *****************"
      header(5)  = "Mass yield [mb] and the mean and variance of the kinetic energy [MeV]"
      header(6)  = "Charge yield [mb] and the mean and variance of the  kinetic energy [MeV]"
      header(7)  = " --------------------------- Energy Spectrum [mb/MeV] --------------------------"
      header(8)  = " ----------------Normalized Energy Probability Spectrum [1/MeV] ----------------"
      header(9)  = " ------------------------ Angular Distributions [mb/sr] ------------------------"
      header(10) = " Double differential cross sections [mb/MeV/sr];" ! Lab angle proceeds
      header(11) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of    n"
      header(12) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of    p"
      header(13) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of    d"
      header(14) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of    t"
      header(15) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  He3"
      header(16) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  He4"
      header(17) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  He6"
      header(18) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  Li6"
      header(19) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  Li7"
      header(20) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  Li8"
      header(21) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  Li9"
      header(22) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  Be7"
      header(23) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  Be9"
      header(24) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of Be10"
      header(25) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  B10"
      header(26) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  B11"
      header(27) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  C11"
      header(28) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  C12"
      header(29) = "          Angular distribution of produced fragments dS/dOm [mb/sr] for energy range(MeV)"
      header(30) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of"
      header(31) = "          Spectra ( dS/dT/dO ) of produced particles at theta="
      header(32) = "          Angle integrated spectra of produced particles"
      header(33) = "          Energy integ. ang. distrib. of prod. particles"
      header(34) = "          Energy integ. ang. distrib. of prod. particles"
      header(35) = "          Angular distribution of produced fragments dS/dOm [mb/sr] for energy range(MeV)"
      header(36) = " Int. x sec "
      header(37) = " Int. xsec "
    endif
    
! Misc. / Complimentary Header
    if(makemisc) then
      mischeader(1)  = " *****************************************************************"
      mischeader(2)  = " End of nuclide yields."
      mischeader(3)  = " of residual nuclei:" ! used with Mass and Charge yields
! For spectra
      mischeader(4)  = "+/-" ! for uncertainties
      mischeader(5)  = "  [deg.]"
      mischeader(6)  = " Lab. angle =" ! for double differential cross sections
      mischeader(7)  = " **********************************  neutrons  *********************************"
      mischeader(8)  = " **********************************  protons   *********************************"
      mischeader(9)  = " ********************************** deuterons  *********************************"
      mischeader(10) = " **********************************  tritons   *********************************"
      mischeader(11) = " **********************************  Helium-3  *********************************"
      mischeader(12) = " **********************************   alphas   *********************************"
      mischeader(13) = " ********************************** neg. pions *********************************"
      mischeader(14) = " ********************************** neut pions *********************************"
      mischeader(15) = " ********************************** pos. pions *********************************"
      mischeader(16) = "     Energy spectrum from   0.0 to  396.0 MeV (zero values suppressed)."
      mischeader(17) = "    Tn  [MeV]            Total                Cascade               Precompound        Total Evaporation"
      mischeader(18) = "    Tp  [MeV]            Total                Casc"
      mischeader(19) = "    Tt  [MeV]           Total               Coalescence            Precompound         Total Evaporation"
      mischeader(20) = "    THe3[MeV]           Total               Coalescence            Precompound         Total Evaporation"
      mischeader(21) = "    THe4[MeV]           Total               Coalescence            Precompound         Total Evaporation"
      mischeader(22) = "    Tpi-[MeV]      Total = Cascade"
      mischeader(23) = "    Tpi0[MeV]      Total = Cascade"
      mischeader(24) = "    Tpi+[MeV]      Total = Cascade"
      mischeader(25) = "  Ang.n           Total                Cascade               Precompound         Total Evaporation"
      mischeader(26) = "  Ang.p           Total                Cascade               Precompound         Total Evaporation"
      mischeader(27) = "  Ang.t           Total              Coalescence             Precompound         Total Evaporation"
      mischeader(28) = "  Ang.He3         Total              Coalescence             Precompound         Total Evaporation"
      mischeader(29) = "  Ang.He4         Total              Coalescence             Precompound         Total Evaporation"
      mischeader(30) = "  Ang.pi-[deg.]     Total = Cascade"
      mischeader(31) = "  Ang.pi0[deg.]     Total = Cascade"
      mischeader(32) = "  Ang.pi+[deg.]     Total = Cascade"
!      mischeader(33) = "  T(MeV)/angle:  12.        15.        20.        35.        50.        65.        80.        90.       100.       120.    dS/dT(mb/MeV)"
      mischeader(34) = "            prod. xsec for all energies="
!      mischeader(35) = "Tmin-Tmax   0.- 160.   2.- 160.   3.- 215.   3.- 250.   2.- 580.   2.- 650.   2.-2500.   4.-2500.   4.-2500.   4.-2500."
!      mischeader(36) = "ang1-ang2      n          p          d          t        He3        He4        He6        Li6        Li7        Li8"
!      mischeader(37) = "ang1-ang2    Li9        Be7        Be9       Be10         B9        B10        B11        B12        C11        C12"
!      mischeader(38) = "ang1-ang2    C13        C14        Z=7        Z=8        Z=9       Z=10       Z=11       Z=12       Z=13       Z=14"
      mischeader(39) = " The program called Fermi breakup"
      mischeader(40) = "- - - - NO INFORMATION - - - -"
      mischeader(41) = "*******"
      mischeader(42) = " energ. int."
      mischeader(43) = "  T(GeV)      He6        Li6        Li7        Be7        Be9        Be10         B8         B10        B11"
      mischeader(44) = "              He6        Li6        Li7        Be7        Be9        Be10         B8         B10        B11"
      mischeader(45) = " theta      He6        Li6        Li7        Be7        Be9        Be10         B8         B10        B11"
      mischeader(46) = "  T(GeV)    PI-MINUS    PI-PLUS   NEUTRONS    PROTONS   DEUTRONS    TRITONS       HE-3       HE-4"
    endif
    
    return
  end subroutine setheader

!  ========================================================