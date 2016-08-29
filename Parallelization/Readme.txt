------------------------------------------------

Parallelizations

------------------------------------------------


how to parallelize a code, or make use of multi- and many-core architectures:

- multi-core architecture

  currently available CPUs, e.g. Intel chips
  contain multiple cores
  
- many-core architecture

  GPUs like Nvidia Tesla cards, 
  contain Streaming Multiprocessor (SM) <-> cores
  

parallelization models:

1. threading

  shared memory, multi-core architectures
  
  multiple threads, mainly C/C++ programming

  software:
  - most operating systems (unix-like) support, set by header file 
     include <pthread.h>, may need library linking -lpthread


2. OpenMP

  shared memory, multi-core architectures
  
  multiple threads, both Fortran and C/C++
  uses either #pragma or directives !$OMP

  software: 
  - most compilers support, set by compilation flag, 
     e.g. gfortran (-fopenmp), intel (-openmp), portland (-mp) 
 

3. MPI

  distributed memory, multi-core architectures 
  
  standard on linux clusters and multi-core architectures
  message-passing schemes, both Fortran and C/C++

  software:
  - most cluster support, or install package e.g. http://www.open-mpi.org/  


4. GPUs
  
  SIMD (single instruction multiple data) memory, many-core architectures
  
  CUDA as multi-threading language, indexed array accesses in C programming
  alternatives: acceleration directives !$ACC or compiler based PGI function calls

  software:
  - Nvidia CUDA http://developer.nvidia.com/cuda-downloads
  
  
5. different language extensions

  Cilk as C extension for multi threading:
    http://software.intel.com/en-us/articles/intel-cilk-plus/
  
  Co-array as Fortran extension mainly on Cray machines:
    http://caf.rice.edu/

