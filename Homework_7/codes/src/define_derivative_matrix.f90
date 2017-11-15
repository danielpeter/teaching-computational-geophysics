  subroutine define_derivative_matrix(xigll,wgll,hprime)

  implicit none

  include "constants.h"
 
  ! Gauss-Lobatto-Legendre points of integration
  double precision, dimension(NGLL) :: xigll

  ! weights
  double precision, dimension(NGLL) :: wgll

  ! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLL,NGLL) :: hprime

  ! function for calculating derivatives of Lagrange polynomials
  double precision, external :: lagrange_deriv_GLL

  integer i1,i2

  ! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wgll,NGLL,GAUSSALPHA,GAUSSBETA)

  ! if number of points is odd, the middle abscissa is exactly zero
  if(mod(NGLL,2) /= 0) xigll((NGLL-1)/2+1) = 0.d0

  ! calculate derivatives of the Lagrange polynomials
  ! and precalculate some products in double precision
  ! hprime(i,j) = h'_i(xigll_j) by definition of the derivative matrix
  do i1=1,NGLL
    do i2=1,NGLL
      hprime(i1,i2) = lagrange_deriv_GLL(i1-1,i2-1,xigll,NGLL)
    enddo
  enddo

  end subroutine define_derivative_matrix

