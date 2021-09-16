! Homework2
  program wave

  integer :: nx,nt
  double precision :: L,dx,dt,Time,mindt,FACTOR
  double precision :: c(100000),kappa(100000),rho(100000)
  double precision :: disp(100000),veloc(100000),stress(100000)
  double precision :: disp1(100000),disp2(100000),veloc1(100000),veloc0(100000),stress1(100000)
  integer :: it,ix,nout,count,ier
  character(len=64) :: filename

  ! meshing parameters
  dx = 0.1d0
  L = 100.d0
  nx = int(100.d0/dx) + 1

  ! model parameters
  do ix = 1,nx
    kappa(ix) = 1.d0
    rho(ix) = 1.d0
    c(ix) = sqrt(kappa(ix)/rho(ix))
  enddo

  ! user output
  print*,'finite-difference scheme: '
  print *,'  space dx = ',sngl(dx)

  ! timing parameters
  FACTOR = 0.25

  ! determines time step size
  mindt = 999999.
  do ix = 1,nx
    dt = FACTOR * dx**2 / c(ix)**2
    if(mindt > dt) mindt = dt
  enddo
  dt = mindt

  print *,'  time step dt = ',sngl(dt)

  Time = 200.d0
  nt = int(Time/dt) + 1
  print *,'  number of time steps nt = ',nt

  ! output file
  nout = 4000                 ! output time step
  count = 0

  ! initial conditions 
  ! 2nd order
  disp(:) = 0.d0
  do ix = 1,nx
    disp(ix) = exp(-dx**2*(ix-501)**2)
  enddo
  disp1(:) = disp(:)

  ! 1st order
  veloc(:) = 0.d0
  stress(:) = 0.d0
  do ix = 1,nx
    stress(ix) = - 2*kappa(ix)*dx**2*(ix-501)*exp(-dx**2*(ix-501)**2)
  enddo

  ! time marching
  do it=1,nt

    !------------------- 2nd order: (u) ---------------------
    !>TODO: implement your FD scheme
    ..
    !<TODO

    ! boundary conditions
    !>TODO: implement your boundary condition
    ..
    !>TODO

    ! for comparison with 1st order velocity
    veloc0(:) = (disp(:) - disp1(:))/dt

    !------------------- 1st order: (v,T) ---------------------
    !>TODO: implement your FD scheme
    ..
    !<TODO

    ! boundary conditions
    !>TODO: implement your boundary condition
    ..
    !<TODO

    !----- Output files written, to be plotted with plot_wave.m
    if (mod(it-1,nout) == 0 .or. it == 1) then
      count = count + 1
      if(count < 7) then
        write(filename,"('S',i6.6,'.dat')") it
        filename = 'figures/'//trim(filename)
        print*,'  plotting: ',trim(filename)

        open(unit=11,file=trim(filename),status='unknown',iostat=ier)
        if (ier /= 0) then
          print *,'could not open figure ',trim(filename)
          print *,'Please check if directory figures/ exists...'
          stop 'Error opening file'
        endif

        do i = 1,nx
          ! format: #position #velocity(2nd-order equations) #velocity(1st-order equations)
          write(11,*) (i-1)*dx, veloc0(i), veloc(i)
        enddo

        close(11)
      endif
    endif

  ! end of time marching
  enddo

  end
