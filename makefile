#		Compiller and various options
F77   = gfortran
FFLGS = #-O1 -ffpe-trap=invalid,zero,overflow -fno-automatic
# only for debugging

SRC = source


#		Directories
CEM = cem_test.EXEs
LAQ = laqgsm_test.EXEs
#
SRCTOT = $(SRC)/main_new.f90 $(SRC)/make_input.f90 $(SRC)/compare_files.f90 \
         $(SRC)/doublediff.f90 $(SRC)/yields.f90 $(SRC)/read_output.f90 \
				 $(SRC)/headers.f90
#		Objects
OBJECTS = srctot.o

#		Flags
OFLGS = $(FFLGS)
F77FLGSC = $(OFLGS) -c
F77FLGS = $(OFLGS) 
CHKFLGS = $(FFLGS)
#


#
test.EXE: $(OBJECTS)
	$(F77) -o $@ $(F77FLGS) $(OBJECTS) 
	@ clear
#
compilation: test.EXE               #Compile and make executable only
test:
	@ make clean
	@ make
	@ ./test.EXE
#
clean:
	@ clear
	@ echo '---> clean:   remove   *.aux  *.r10  *.r11  *.txt *.png *.dat'
	@ rm -f  *.aux *.inf  *.r10 *.r11 *.txt *.png *.dat *~ *#
#
veryclean:
	@ clear
	@ make clean
	@ echo '--->		test.EXE *.inp *.out'
	@ rm -f test.EXE *.inp *.out
#

srctot.f90: $(SRCTOT)
	cat $(SRCTOT) > $@
#
srctot.o : srctot.f90
	$(F77) $(F77FLGSC) srctot.f90
#
