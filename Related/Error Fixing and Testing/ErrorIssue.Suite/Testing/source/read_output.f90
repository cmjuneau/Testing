!  ========================================================

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

!  ========================================================