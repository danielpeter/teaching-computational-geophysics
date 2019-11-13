! number of spectral elements
  integer, parameter :: NSPEC = 20

! number of GLL points (polynomial degree plus one)
  integer, parameter :: NGLL = 7

! number of global points
  integer, parameter :: NGLOB = (NGLL-1) * NSPEC + 1

! for the Gauss-Lobatto-Legendre points and weights
  double precision, parameter :: GAUSSALPHA = 0.d0,GAUSSBETA = 0.d0
