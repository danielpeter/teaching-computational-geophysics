# Makefile

# compilers
CC=gcc
F90=gfortran
MPIF90=mpif90
NVCC=nvcc

# targets
default: xhello xhello_thread xhello_openmp xhello_mpi xhello_cuda

# initial, single threaded
xhello: hello.f90
	$(F90) -o xhello hello.f90

# parallelized versions
xhello_thread: hello_thread.c
	$(CC) -o xhello_thread hello_thread.c	-lpthread

xhello_openmp: hello_openmp.f90
	$(F90) -fopenmp -o xhello_openmp hello_openmp.f90

xhello_mpi: hello_mpi.f90
	$(MPIF90) -o xhello_mpi hello_mpi.f90

xhello_cuda: hello_cuda.cu
	$(NVCC) -arch=sm_20 -o xhello_cuda hello_cuda.cu

clean:
	rm -f *.o x* a.out

