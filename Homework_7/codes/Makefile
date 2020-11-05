# Makefile

# GNU compiler
F90 = gfortran
FLAGS = -O2 -g -I.

# directories
O = obj
E = bin
S = src

all: diffusion

diffusion: constants.h $O/diffusion.o $O/gll_library.o $O/lagrange_poly.o $O/define_derivative_matrix.o
	${F90} $(FLAGS) -o $E/xdiffusion $O/diffusion.o $O/gll_library.o $O/lagrange_poly.o $O/define_derivative_matrix.o

clean:
	rm -f $O/*.o $E/xdiffusion

# rule
$O/%.o: $S/%.f90 constants.h
	${F90} $(FLAGS) -c -o $@ $<

