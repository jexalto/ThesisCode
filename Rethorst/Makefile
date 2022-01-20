make:
	gfortran -ffree-form -c src/mod_modifiedbessel.f90 src/mod_integrals.f90 src/Gmatrix.f90 src/mod_math.f90
	gfortran -O -Wall -fcheck=all -g -fbacktrace Gmatrix.o mod_modifiedbessel.o mod_integrals.o mod_math.o
	valgrind --leak-check=full --track-origins=yes --show-leak-kinds=all ./a.out

bessel:
	gfortran src/mod_modifiedbessel.f90