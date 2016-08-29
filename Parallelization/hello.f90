!
! single-process version for fortran90
!
!
! compile with:
!
! > gfortran hello_single.f90
!
! run with:
! > ./a.out
!

program hello
  
  implicit none
  
  integer,parameter :: N = 10
  
  real,dimension(N) :: A,B,C
  
  real:: sum,total_sum
  integer :: i,istart,iend,isum

  print *,'Hello world'

  ! initialization
  A(:) = 1.0
  B(:) = 2.0
  C(:) = 0.0
  sum = 0.0
  isum = 0
  
  ! divides loop 
  istart = 1
  iend = N
  
  do i = istart,iend
    C(i) = A(i) + B(i)
    sum = sum + C(i)
    isum = isum + 1

    print *,'isum = ',isum
  enddo
  total_sum = sum

  ! user output
  print *
  print *,'result:'
  print *,'  sum = ',sum
  print *,'  total sum = ',total_sum

end program