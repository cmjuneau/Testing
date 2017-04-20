  program compare
    implicit none

!  ========================================================
!
!  Program created by CMJ, 08/16 to compare output files line by line
!  Modified by CMJ, 08/16 to allow for output reading
!
!
!  This program can do 2 things.  First, when the variable
!  'testing' is set to false, this program will compare 2
!  files and store differences into an output while displaying
!  results on the terminal window.  This option can only compare
!  2 files at this time (additional files would not be difficult
!  to add).
!  The 2nd option, when 'testing' is set to true, allows us to
!  quickly take data from CEM output files and display the results
!  side to side.  This is mostly meant to 'prep' the data, and later
!  MATLAB or ROOT can be used to extract the data and produce graphs.
!
!  These tools were created due to the frequency that this kind of
!  testing will take place.  This quickly and effectively compares
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

!		----- Setting up variables for Comparisons -----
    integer :: MaxLines, MaxChar, nl1, nl2, nlmin, i, j, totmismatch, start, temp
		integer :: mismatch(10000), NumComp, choice, num1, num2, num3
		integer :: totaltime, incrementlength, numincrement, numfiles, NumTests(100)
    
		character(LEN = 200) :: input, instructions(4)
		character(LEN =  30) :: files(10000), outname
    character(LEN = 1) ::   testcomp
		
!		----- Setting up variables for Output Numerical Analysis -----
    integer :: num4, numf(4), numi, headmult, nl, num
    real :: yieldmult
    
    character(LEN = 300) :: header, mischeader, output
    
    logical :: read_input, testing, getNum, store_mult, store_dds
    logical :: read_output
		
!		----- Setting up common block for output -----
    common  /header/     header(100), mischeader(50), headmult
    common  /outp/       output(3,30000), nl(3)
    common  /maxline/    MaxLines
    common  /yieldmult/  yieldmult(10,9,9,5)
		
!  ========================================================

    read_input = .true.   ! For ease in reading code; keep value as ".true."
    getNum = .true.       ! For ease in reading code; keep value as ".true."
    read_output = .true.  ! Reads output, always meant to be on.
    store_mult = .false.  ! to store multiplicities and yields
    store_dds = .true.    ! to store double differntial spectra
    output = "- - - - NO INFORMATION - - - -"
    files  = "- - - - NO INFORMATION - - - -"
    MaxLines = 30000
  
! Deciding if we are comparing or testing
10  continue
    print *, "Test, compare, or make files (T/C/M)?"
    read(*,*) testcomp
!    testcomp = 'm'
    if (testcomp.eq."T" .OR. testcomp.eq."t") then
      testing = .true.
      print *, "Proceeding to test files"
      input = 'test.prm'
    elseif(testcomp.eq."C" .OR. testcomp.eq."c") then
      testing = .false.
      print *, "Proceeding to compare files"
      input = 'comp.prm'
    elseif(testcomp.eq.'M' .OR. testcomp.eq.'m') then
      print *, 'Creating input files for CEM and LAQGSM...'
      print *, ' '
      call make_Inputs
      print *, 'Input files created'
      print *, " "
      go to 20
      stop
    else
      print *, "Please enter 'T' to test, 'C' to compare, or 'M' to make input files..."
      print *,
      go to 10
    endif

!		----- Getting Files from user -----
    print *, '========================================'
!    choice = 0  
! Getting input name
!    if (choice.eq.0) then
!  		input = 'input.inp'
!    else
!      print *, 'Name of input file:'
!      read(*,*) input
!    endif
  
    if(read_input) then
       open(999, file = input, status = 'old')
       print *, 'Reading input file'
    
!  reading instructions in input
	    do i = 1, 4
		    read(999,*) instructions(i)	! Dummy variables
	    enddo
      
!  getting number of comparisons
		  read(999, *) NumComp
	    print *, 'Comparisons to be made: ', NumComp
    
!  error message, allowed to compare up to 100 files
	    if (NumComp.gt.100) then
		    print *, "Too many files to compare, allowed max of 100."
	      stop
	    endif
    
!  reading file names from input
	    do numfiles = 1,  NumComp
        if(getNum) then
! Max of 3 files to read, 1 output included.  file no. are x1,x2,x3,x4, > 20
		      num1 = 1 + (numfiles + 1) * 10
		      num2 = 2 + (numfiles + 1) * 10
		      num3 = 3 + (numfiles + 1) * 10
		      num4 = 4 + (numfiles + 1) * 10
          numf(1) = num1
          numf(2) = num2
          numf(3) = num3
          numf(4) = num4
        endif
        
! Indicates the number of numerical comparisons requested
		    read(999, *) NumTests(numfiles)
        if(NumTests(numfiles).gt.10) go to 15
        do i = 1, NumTests(numfiles)+1 ! the "+1" term represents the output name
		      read(999, *) files(numf(i))
        enddo
	    enddo
    endif
	  print *, "Done reading input..."
	  print *, " "

    if (testing) then
      
! Sets various headings for sections
      call setheader
      
! Start of main testing loop
		  do numfiles = 1, NumComp
        yieldmult = 0.d0
        
        print *, "==============================="
			  print *, "Test No. ", numfiles
        
! Setting file numbers in numf
        if(getNum) then
! Max of 3 files to compare/compute numerical analysis, 1 output
		      num1 = 1 + (numfiles + 1) * 10
		      num2 = 2 + (numfiles + 1) * 10
		      num3 = 3 + (numfiles + 1) * 10
		      num4 = 4 + (numfiles + 1) * 10
          numf(1) = num1
          numf(2) = num2
          numf(3) = num3
          numf(4) = num4
        endif
        
! Reads files and stores data
        do i = 1, NumTests(numfiles)
          headmult = 0 ! for testing
          numi = numf(i)
        
! Reading output
          if(read_output) then
            print *, " "
            print *, 'Reading output file ', trim(files(numi))
            call readOutput(files, numf, i)
          endif
        
! Storing multiplicity data
          if(store_mult) then
            print *, "Storing multiplicities..."
            call store_yieldmult(i)
          endif
        
! Storing double differential spectra
          if(store_dds) then
            print *, "Storing double differential spectra..."
            call store_doubledif(i)
          endif
        enddo
        

! Determining name of output file
        if(files(numf(4)).ne."- - - - NO INFORMATION - - - -") then
          num = numf(4)
          outname = files(num)
        else
          num = numf(3)
          outname = files(num)
        endif
        
! multiplicities and yields in output stored
        print *, " "
        if(store_mult) then
          call make_multyields(outname, num)
        endif
! Double differential spectra for all particles
        if(store_dds) then
          call make_doubledif (outname, num)
        endif
      enddo
    else
      
! Start of main comparing loop
		  do numfiles = 1,  NumComp
			  print *, "Comparison No. ", numfiles
        if(getNum) then
! Max of 3 files to compare/compute numerical analysis, 1 output
		      num1 = 1 + (numfiles + 1) * 10
		      num2 = 2 + (numfiles + 1) * 10
		      num3 = 3 + (numfiles + 1) * 10
		      num4 = 4 + (numfiles + 1) * 10
        endif
!  comparing sets of files
			  call compare_multiple_files(files, num1, num2, num3)
		  enddo
    endif

    go to 20
    
 15 continue
    print *, '***   WARNING: Input Error   ***'
    print *, '***   End of File reached before done getting'
    print *, 'file names.  Issue may be with number before'
    print *, 'file names or lack of file names, please'
    print *, 'correct   ***'
    stop

 20 continue
    print *, "====================================================="
    print *, "Program done, see resulting output files for details."
    print *, "====================================================="

  end program compare

!  ========================================================
