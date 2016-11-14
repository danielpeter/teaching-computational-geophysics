program diffusion

  implicit none

  !---------------------------------------------------------------------------------
  ! User parameters

  include "constants.h"

  ! number of timesteps
  integer, parameter :: NSTEP = 40000

  ! time step in seconds
  double precision, parameter :: DT = 100000000. ! s

  ! fixed boundary conditions
  logical, parameter :: FIXED_BC = .true.

  ! even/uneven element spacing
  logical, parameter :: UNIFORM_ELEMENTS = .true.

  ! constant/bi-material properties
  logical, parameter :: HOMOGENEOUS_MATERIAL = .true.

  ! model parameters  (SI)
  double precision, parameter :: LENGTH = 3.0d+03 ! m
  double precision, parameter :: DENSITY = 2.5d+03 ! kg/m^3
  double precision, parameter :: THERMALCONDUCTIVITY_1 = 10.0d-01 ! cal/m/s/K
  double precision, parameter :: THERMALCONDUCTIVITY_2 = 2.0d-01 ! cal/m/s/K
  double precision, parameter :: HEATCAPACITY = 0.3d+03 ! cal/kg/K

  !---------------------------------------------------------------------------------

  integer :: ispec,i,j,iglob,itime,iglob1,iglob2,iglobj,nmid
  integer :: nleft,nright
  double precision :: left_length,right_length

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
  double precision, dimension(NGLL,NSPEC) :: rho,heat_capacity,thermal_conductivity

  ! Jacobian `matrix' and Jacobian
  double precision, dimension(NGLL,NSPEC) :: dxidx,jacobian

  ! local mass matrix
  double precision :: mass_local

  ! global mass matrix
  double precision, dimension(NGLOB) :: mass_global

  ! rhs matrix
  double precision, dimension(NGLOB) :: rhs_global

  ! temperature and temperature time derivative
  double precision, dimension(NGLOB) :: temperature,dtemperature_dt

  ! local to global numbering
  integer, dimension(NGLL,NSPEC) :: ibool

  ! time marching
  double precision :: deltat,deltatover2
  double precision :: dh,diffusivity,time_step

  ! end fluxes
  double precision :: flux_1,flux_NGLOB

  ! end temperatures
  double precision :: temperature_1,temperature_NGLOB,stiffness_local

  ! derivatives
  double precision :: dtdx,flux,rhs_local,temp(NGLL)

  ! movie
  character(len=50) :: moviefile

  !++++++++++++++++++++++++++++++++++++++++++++++++++

  ! define derivatives
  call define_derivative_matrix(xigll,wgll,hprime)

  ! uniform grid
  print *,'uniform elements: ',UNIFORM_ELEMENTS
  if (UNIFORM_ELEMENTS) then
    ! evenly spaced achors between 0 and 1
    print *,'  even element spacing:'
    do ispec = 1,NSPEC
      x1(ispec) = LENGTH*dble(ispec-1)/dble(NSPEC)
      x2(ispec) = LENGTH*dble(ispec)/dble(NSPEC)
      print *, ispec,sngl(x1(ispec)),sngl(x2(ispec)),'size = ',sngl(x2(ispec)-x1(ispec))
    enddo
  else
    ! unevenly spaced achors between 0 and 1
    print *,'  uneven element spacing:'
    left_length = LENGTH/3
    right_length = LENGTH-left_length
    nleft = floor(NSPEC/2.)
    nright = NSPEC-nleft
    do ispec = 1,NSPEC
      if (ispec <= nleft) then
        x1(ispec) = left_length*dble(ispec-1)/dble(nleft)
        x2(ispec) = left_length*dble(ispec)/dble(nleft)
      else
        x1(ispec) = left_length + right_length*dble(ispec-nleft-1)/dble(nright)
        x2(ispec) = left_length + right_length*dble(ispec-nleft)/dble(nright)
      endif
      print *, ispec,sngl(x1(ispec)),sngl(x2(ispec)),'size = ',sngl(x2(ispec)-x1(ispec))
    enddo
  endif

  ! set up the mesh properties
  print *,'homogeneous material: ',HOMOGENEOUS_MATERIAL
  if (HOMOGENEOUS_MATERIAL) then
    nmid = NSPEC
  else
    nmid = floor(dble(NSPEC/2))
  endif
  do ispec = 1,NSPEC
    do i = 1,NGLL
      if (ispec <= nmid) then
        thermal_conductivity(i,ispec) = THERMALCONDUCTIVITY_1
      else
        thermal_conductivity(i,ispec) = THERMALCONDUCTIVITY_2
      endif
      rho(i,ispec) = DENSITY
      heat_capacity(i,ispec) = HEATCAPACITY
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
  mass_global(:) = 0.
  ..

  ! estimate the time step 'time_step'
  deltat = DT
  deltatover2 = DT / 2
  print *,'time step estimate: ',deltat,' seconds'


  ! set up the boundary conditions
  !>TODO: put your code here
  temperature_1  = ..
  temperature_NGLOB = ..

  ! initialize
  temperature(:) = 0.0
  dtemperature_dt(:) = 0.0

  ! time loop
  do itime = 1,NSTEP

    ! fixed temperature
    temperature(1) = temperature_1
    temperature(NGLOB) = temperature_NGLOB

    ! update temperature
    !>TODO: put your code here
    temperature(:) = ..

    rhs_global(:) = 0.0

    ! element loop
    do ispec = 1,NSPEC

      ! GLL point loop
      do i = 1,NGLL
        iglob = ibool(i,ispec)

        ! local contribution
        !>TODO: put your code here
        rhs_local = 0.0
        ..

        ! boundary conditions
        !>TODO: put your code here
        if (FIXED_BC) then
          if (ispec == 1 .and. i == 1) then
            ! left side
            iglob1 = ibool(1,ispec)
            iglob2 = ibool(NGLL,ispec)
            rhs_local = ..

          else if (ispec == NSPEC .and. i == NGLL) then
            ! right side
            iglob1 = ibool(1,ispec)
            iglob2 = ibool(NGLL,ispec)
            rhs_local = ..

          endif
        endif

        ! assembly
        !>TODO: put your code here
        rhs_global(iglob) = ..
           
      enddo
    enddo

    ! temperature increment
    !>TODO: put your code here
    dtemperature_dt(:) = ..

    ! time scheme: corrector term
    !>TODO: put your code here
    temperature(:) = ..

    !if (maxval(abs(temperature)) > 1000) then
    !  print *, ' at time step itime = ', itime
    !  stop 'Error : temperature larger than 1000'
    !endif

    ! write out snapshots
    if(mod(itime-1,1000) == 0) then
      ! file output
      write(moviefile,'("figures/snapshot",i5.5,".dat")') itime
      open(unit=10,file=moviefile,status='unknown')
      do iglob = 1,NGLOB
        write(10,*) sngl(x(iglob)),sngl(temperature(iglob))
      enddo
      close(10)
      print *, 'time step = ', itime,' - file written: ',trim(moviefile)
    endif

  enddo ! end time loop

end program diffusion
