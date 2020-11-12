  program wave

  implicit none

  !---------------------------------------------------------------------------------
  ! User parameters

  include "constants.h"

  ! number of timesteps
  integer, parameter :: NSTEP = 400

  ! time step in seconds
  double precision, parameter :: DT = 0.25

  ! fixed boundary conditions
  logical, parameter :: FIXED_BC = .true.

  ! model parameters  (SI)
  double precision, parameter :: LENGTH = 100
  double precision, parameter :: DENSITY = 1
  double precision, parameter :: SHEARMODULUS = 1

  !---------------------------------------------------------------------------------

  integer :: ispec,i,j,k,iglob,itime
  integer :: iglobj,iglob1,iglob2

  ! Gauss-Lobatto-Legendre points of integration
  double precision, dimension(NGLL) :: xigll

  ! weights
  double precision, dimension(NGLL) :: wgll

  ! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLL,NGLL) :: hprime

  ! anchors
  double precision, dimension(NSPEC) :: x1,x2

  ! global grid points
  double precision, dimension(NGLOB) :: x

  ! material properties
  double precision, dimension(NGLL,NSPEC) :: rho,shear

  ! Jacobian `matrix' and Jacobian
  double precision, dimension(NGLL,NSPEC) :: dxidx,jacobian

  ! local mass matrix
  double precision :: mass_local

  ! global mass matrix
  double precision, dimension(NGLOB) :: mass_global

  ! temperature and temperature time derivative
  double precision, dimension(NGLOB) :: displ,veloc,accel

  ! local to global numbering
  integer, dimension(NGLL,NSPEC) :: ibool

  ! time marching
  double precision :: deltat,deltatover2,deltatsquareover2

  ! end gradients
  double precision :: grad_1,grad_NGLOB

  ! end temperatures
  double precision :: displ_1,displ_NGLOB

  double precision :: rhs_local,stiffness_local

  ! for double-loops
  double precision :: tmp1l,tmp2l,fac1
  double precision :: tmp1(NGLL)

  ! movie
  character(len=50) :: moviefile

  !++++++++++++++++++++++++++++++++++++++++++++++++++

  ! define polynomial derivatives & weights
  call define_derivative_matrix(xigll,wgll,hprime)

  ! evenly spaced anchors between 0 and 1
  do ispec = 1,NSPEC
    x1(ispec) = LENGTH*dble(ispec-1)/dble(NSPEC)
    x2(ispec) = LENGTH*dble(ispec)/dble(NSPEC)
  enddo

  ! set up the mesh properties
  do ispec = 1,NSPEC
    do i = 1,NGLL
      rho(i,ispec) = DENSITY
      shear(i,ispec) = SHEARMODULUS
    enddo
  enddo

  ! Jacobian
  do ispec = 1,NSPEC
    do i = 1,NGLL
      dxidx(i,ispec) = 2.0 / (x2(ispec)-x1(ispec))
      jacobian(i,ispec) = (x2(ispec)-x1(ispec)) / 2.0
    enddo
  enddo

  ! set up local to global numbering
  iglob = 1
  do ispec = 1,NSPEC
    do i = 1,NGLL
      if(i > 1) iglob = iglob + 1
      ibool(i,ispec) = iglob
    enddo
  enddo

  ! get the global grid points
  do ispec = 1,NSPEC
    do i = 1,NGLL
      iglob = ibool(i,ispec)
      x(iglob) = 0.5*(1.-xigll(i))*x1(ispec)+0.5*(1.+xigll(i))*x2(ispec)
    enddo
  enddo

  ! calculate the global mass matrix 'mass_global'
  !>TODO: put your code here
  mass_global(:) = 0.0
  ..

  ! estimate the time step DT, DT/2, DT^2/2
  deltat = DT
  deltatover2 = DT / 2.0
  deltatsquareover2 = DT*DT / 2.0
  print *,'time step estimate: ',deltat,' seconds'

  ! set up the boundary conditions
  !>TODO: put your code here
  displ_1 = ..
  displ_NGLOB = ..
  grad_1 = ..
  grad_NGLOB = ..
  
  ! initialize
  displ(:) = 0.0
  veloc(:) = 0.0
  accel(:) = 0.0

  ! add source
  ! initial displacement
  !>TODO: put your code here
  ..

  ! time loop
  do itime = 1,NSTEP

    ! "predict" displacement, velocity, initialize acceleration
    ! Newmark time scheme: predictor terms
    !>TODO: put your code here
    displ(:) = ..
    veloc(:) = ..
    accel(:) = ..

    ! boundary conditions: Dirichlet boundary condition
    ! (zero displacement at boundary points)
    ..

    ! element loop
    do ispec = 1,NSPEC

      ! contribution from local element & assembly on global points
      !>TODO: put your code here
      ..

      ! boundary conditions
      !>TODO: put your code here
      ..

    enddo
        
    ! updates acceleration
    !>TODO: put your code here
    accel(:) = ..

    ! "correct" acceleration, velocity, displacement
    ! Newmark time scheme: corrector term
    !>TODO: put your code here
    veloc(:) = ..

    ! write out snapshots
    if(mod(itime-1,25) == 0) then
      ! file output
      write(moviefile,'("figures/snapshot",i5.5,".dat")') itime
      open(unit=10,file=moviefile,status='unknown')
      do iglob = 1,NGLOB
        write(10,*) sngl(x(iglob)),sngl(displ(iglob))
      enddo
      close(10)
      print *, 'time step = ', itime,' - file written: ',trim(moviefile)

      ! plotting by gnuplot
      if (itime == 1) then
        call system('echo "#!/sw/bin/gnuplot" > tmp.gnuplot')
        call system('echo "set yrange [-1:1]" >> tmp.gnuplot')
      endif
      call system('echo "plot \"'//trim(moviefile)//'\" w lp" >> tmp.gnuplot')
      call system('echo "pause -1 \"Hit return to continue\"" >> tmp.gnuplot')
    endif

  enddo ! end time loop

  ! plotting w/ gnuplot
  print *, ""
  print *, "for visualization, use for example command: "
  print *, "  > gnuplot tmp.gnuplot"
  print *, ""
  !call system('gnuplot tmp.gnuplot')

  end program wave
