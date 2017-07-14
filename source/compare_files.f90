!  ========================================================

  subroutine compare_multiple_files(input1, input2, outFile, num)
    implicit none

! This program compares files that are input by the user.
! This program reads in 2 file names, and opens and reads
! the files.  The data is stored here, and each string is
! compared character for character.  Empty lines in the
! files are filled with "No Line Information" so it is
! obvious that no information was in the line.


! ----- Setting up variables -----
    integer :: MaxLines, nl, nlmax, num, i, j, k, totmismatch, start
    integer :: temp, temp2
    integer :: mismatch(30000), test1, test2, test3

    character(LEN = 300) :: output, line1(30000), line2(30000)
    character(LEN =  34) :: outname
    character(LEN =  30) :: input1, input2, outFile
    character(LEN =   4) :: ext

    logical :: makeOutput, first

    common  /maxline/  MaxLines
    common  /outp/     output(10,30000), nl(3)
    data first / .true. /

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

    ! comparing lines in files
    totmismatch = 0
    mismatch = 0
    start = 31
    nlmax= max(nl(1), nl(2))
    do j = start, nlmax
       if( line1(j).ne.line2(j) .AND. (j.ne.39 .AND. j.ne.51 .AND. j.ne.88)) then
          ! print *, 'Mismatch at line', j
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
			
    ! ----- Formatting output file -----
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
          temp = 0
          do j = start, nlmax
             if ( mismatch(j).eq.1 ) then
                temp = temp + 1
                temp2 = totmismatch - temp
                if ( temp.gt.50 .AND. temp2.gt.10) then
                   if ( first ) then
                      first = .false.

                      do k = 1, 5
                         write(num,*) " ."
                      enddo

                      write(num,*) 'ETC.'

                      do k = 1, 5
                         write(num,*) " ."
                      enddo

                   endif
                else
                   write(num,*) 'Mismatch at line', j
                   write(num,*) '   -> From ', trim(input1), ':		', trim(line1(j))
                   write(num,*) '   -> From ', trim(input2), ':		', trim(line2(j))
                   write(num,*) ''
                endif
             endif
          enddo
100       continue
       endif
       close(num)	
       print *, "--  Output file ", trim(outname), " created with results  --"
       print *, " "
    endif

    return
  end subroutine compare_multiple_files
	 
!  ========================================================
