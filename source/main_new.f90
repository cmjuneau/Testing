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
