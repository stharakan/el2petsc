Eventually will be an interface between Elemental and PETSC C++ libraries
To use, make sure that PETSC_DIR and EL_DIR are defined (these are referenced
in the makefile)


Current functions

	- El2Petsc_vec: Elemental vector (VC,STAR only) --> PETSC vector

	- Petsc2El_vec: PETSC vector --> Elemental vector (VC,STAR only)

See el_petsc_utils.hpp for function signatures and el_petsc_utils.cpp for
source code

Testing

	Run the following on maverick
	
	make test 
	cd runs/
	./test_el2petsc
	
	from this directory to create test_funcs.exe, an 
	executable that will test functions above by creating a PETSC vector 
	x = {0,1,2 ..}, and transforming it to Elemental and back. The next 
	two lines will submit a couple jobs to test performance 


