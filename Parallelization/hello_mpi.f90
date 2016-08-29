!
! MPI version for fortran90
!
!
! compile with:
!
! > mpif90 hello_mpi.f90
!
! run with:
! > mpirun -np 2 ./a.out
!

program hello
  
  implicit none
  
  include 'mpif.h'
  
  integer,parameter :: N = 10
  
  real,dimension(N) :: A,B,C
  
  real:: sum,total_sum
  integer :: i,istart,iend,isum

  ! mpi parameters
  integer :: myrank, size, ier, tag
  integer :: status(MPI_STATUS_SIZE)
  
  ! initializes MPI
  call MPI_INIT(ier)
  
  ! total number of mpi processes
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ier)

  ! accessing "active/local" process id 
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ier)
  
  print*, 'node', myrank, ': Hello world'

  ! initialization
  A(:) = 1.0
  B(:) = 2.0
  C(:) = 0.0
  sum = 0.0
  isum = 0
  
  ! divides loop 
  istart = myrank * N/size + 1
  iend = myrank * N/size + N/size
  
  do i = istart,iend
    C(i) = A(i) + B(i)
    sum = sum + C(i)
    isum = isum + 1

    print*,'isum = ',isum
  enddo

  ! sums all results from different processes
  call MPI_REDUCE(sum,total_sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ier)
  
  ! user output
  if( myrank == 0 ) then
    print*
    print*,'result:'
    print*,'  sum = ',sum
    print*,'  total sum = ',total_sum
  endif
  
  ! exiting MPI processes
  call MPI_FINALIZE(ier)
   
end program