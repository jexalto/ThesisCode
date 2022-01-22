make:
	gfortran -c src/main.f90 src/mod_gvalues_even.f90 src/mod_gvalues_odd.f90 src/mod_integrals.f90 src/mod_modifiedbessel.f90 src/mod_math.f90
	# gfortran -v -g -fbounds-check src/mod_modifiedbessel.f90 src/mod_integrals.f90 src/Gmatrix.f90 src/mod_math.f90
	gfortran -O3 -Wall -fcheck=all -g -fbacktrace main.o mod_gvalues_even.o mod_gvalues_odd.o mod_integrals.o mod_modifiedbessel.o  mod_math.o
	# valgrind --leak-check=full --track-origins=yes --show-leak-kinds=all  -ffast-math -Ofast 
	./a.out