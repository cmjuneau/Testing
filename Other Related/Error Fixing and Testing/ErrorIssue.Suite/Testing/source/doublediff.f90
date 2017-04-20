!  ========================================================

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

!  ========================================================