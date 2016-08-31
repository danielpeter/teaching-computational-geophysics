# Parallelizations


How to parallelize a code, or make use of multi- and many-core architectures? Here are some options:

- *multi-core* architecture

  currently available CPUs, e.g. Intel chips
  contain multiple cores

- *many-core* architecture

  GPUs like NVIDIA (Tesla) cards,
  contain Streaming Multiprocessor (SM) <-> cores


## Parallelization models

- **Threading**

  *shared* memory, multi-core architectures

  multiple threads, mainly C/C++ programming

  *Software requirements:*  
  - most operating systems (unix-like) support, set by header file
     include <pthread.h>, may need library linking -lpthread


- **OpenMP**

  *shared* memory, multi-core architectures

  multiple threads, both Fortran and C/C++
  uses either #pragma or directives !$OMP

  *Software requirements:*  
  - most compilers already support, set by compilation flag,
     e.g. gfortran (-fopenmp), intel (-openmp), portland (-mp)


- **MPI**

  *distributed* memory, multi-core architectures

  standard on linux clusters and multi-core architectures
  message-passing schemes, both Fortran and C/C++

  *Software requirements:*  
  - most cluster support, or install package e.g. http://www.open-mpi.org/  


- **GPUs**

  SIMD (single instruction multiple data) memory, many-core architectures

  **CUDA** as multi-threading language, indexed array accesses in C programming
  alternatives: acceleration directives !$ACC or compiler based PGI function calls

  *Software requirements:*  
  - NVIDIA CUDA http://developer.nvidia.com/cuda-downloads

  **OpenCL** as a standard for parallel programming of heterogeneous systems:  
    https://www.khronos.org/opencl/


- different language extensions

  **Cilk** as C extension for multi threading:  
    http://software.intel.com/en-us/articles/intel-cilk-plus/

  **Co-array** as Fortran extension mainly on Cray machines:  
    http://caf.rice.edu/

  **Chapel** as new parallel programming language (resembles C) developed at Cray:  
    http://chapel.cray.com
