make:
	gfortran -ffree-form -c src/mod_modifiedbessel.f90 src/mod_maths.f90 src/Gmatrix.f90
	gfortran -o 00_exe Gmatrix.o mod_modifiedbessel.o mod_maths.o
	./00_exe

bessel:
	gfortran src/mod_modifiedbessel.f90