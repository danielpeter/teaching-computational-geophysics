# Makefile

# GNU compiler
F90 = gfortran
FLAGS = -O2 -g -I. -Wall

# directories
O = obj
E = bin
S = src

all: wave

wave: constants.h $O/wave.o $O/gll_library.o $O/lagrange_poly.o $O/define_derivative_matrix.o
	${F90} $(FLAGS) -o $E/xwave $O/wave.o $O/gll_library.o $O/lagrange_poly.o $O/define_derivative_matrix.o

clean:
	rm -f $O/*.o *.o $E/xwave 

# rule
$O/%.o: $S/%.f90 constants.h
	${F90} $(FLAGS) -c -o $@ $<

