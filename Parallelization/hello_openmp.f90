!
! OpenMP version for fortran90
!
!
! optional, instead of routine call omp_set_num_threads(),
! set environment variable OMP_NUM_THREAD:
!     bash: 
!         > export OMP_NUM_THREAD=2
!     csh:  
!         > setenv OMP_NUM_THREAD 2
!
! compile with:
!
! > gfortran -fopenmp hello_openmp.f90
!
! run:
! > time ./a.out
!

program hello
  
  use omp_lib
  
  implicit none

  integer,parameter :: NUM_THREADS = 2

  integer,parameter :: N = 10
  real,dimension(N) :: A,B,C
  
  real :: sum
  integer :: i,isum  
  integer:: id, nthreads  

  ! sets number of parallel OpenMP threads
  call omp_set_num_threads(NUM_THREADS)
    
  !$OMP parallel private(id)  
  id = omp_get_thread_num()
  print*, 'Hello World from thread', id
  
  !$OMP barrier
  if ( id == 0 ) then
    nthreads = omp_get_num_threads()
    print*
    print*, 'There are', nthreads, 'threads'
    print*
  end if
  !$OMP end parallel
  
  ! initialization
  A(:) = 1.0
  B(:) = 2.0
  C(:) = 0.0
  sum = 0.0
  isum = 0

  ! i is private by default    
  !$OMP parallel shared(A,B,C,sum)  private(isum)   
  
  !$OMP do 
  do i=1,N
    C(i) = A(i) + B(i)
    sum = sum + C(i)
    isum = isum + 1

    print*,'isum = ',isum
  enddo
  !$OMP end do

  !$OMP end parallel
  
  print*
  print*,'result:'
  print*,'  isum = ',isum
  print*,'  sum  = ',sum

end program
