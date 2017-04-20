!  ========================================================

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

!  ========================================================