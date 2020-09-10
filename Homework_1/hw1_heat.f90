! Homework1
  program heat

  implicit none
  integer :: nx,nt
  double precision :: L,dx,D,dt,Time,FACTOR
  double precision,dimension(:),allocatable :: T,Told
  integer :: it,ix
  character(len=64) :: arg(1)

  ! meshing parameters
  dx = 1.0d0
  L = 100.d0
  nx = int(100.d0/dx)+1
  allocate(T(nx),Told(nx))
  
  ! model parameters (diffusivity)
  D = 1.d0

  ! timing parameters
  call getarg(1,arg(1))
  if( trim(arg(1)) == "" ) then
    FACTOR = 0.3d0
  else
    read(arg(1),*) FACTOR
  endif

  ! user output
  print*,'finite-difference scheme: '
  print*,'  factor = ',sngl(FACTOR)

  ! scales time step size
  dt = FACTOR * dx**2 / D
  
  Time = 25.d0
  nt = int(Time/dt) + 1

  ! initial conditions (explicit loop)
  T(:) = 0.d0
  T(int(nx/2)) = 1.d0

  print*,'  time  dt = ',sngl(dt)
  print*,'  space dx = ',sngl(dx)

  ! time marching
  do it = 1,nt

    ! update of the temperature field
    Told(:) = T(:)     
    
    ! finite-difference scheme
    !>TODO: implement your FD scheme
    do ix = 2,nx-1
      T(ix) = Told(ix) + ...
    enddo
    !<TODO

    ! boundary conditions
    T(1) = 0.d0          
    T(nx) = 0.d0

    ! Output files written, to be plotted with plot_heat.m  
    call plot_snapshot(Told,nx,it,dx,dt,nt)      
    
  enddo                    

  ! last snapshot
  call plot_snapshot(Told,nx,it,dx,dt,nt)

  end

!
!-------------------------------------------------------------------------------------------------
!

  subroutine plot_snapshot(Told,nx,it,dx,dt,nt)
  
  implicit none
  integer :: nx,nt
  double precision:: Told(nx),dx,dt
  integer :: it
  
  ! local parameters
  character(len=64) filename
  integer :: i,ier,nout
  double precision :: Tanalytical(nx),err(nx)
  double precision :: t,D
  double precision,parameter:: pi = 3.1415926535897931
  integer,save :: count = 0
  
  nout = nt/10
  
  if (mod(it-1,nout) == 0 .or. it == 1 .or. it == nt) then
    
    ! analytical solution
    ! heat equation: impulse-response function
    t = (it-1)*dt
    if( t > 0.0 ) then
      D = 1.0d0
      do i = 1,nx
        Tanalytical(i) = dx * 1.0d0 / sqrt(4.d0*pi*D*t) * exp( -((i-1)*dx-(int(nx/2)-1)*dx)**2/(4.d0*D*t))
      enddo
    else
      Tanalytical(:) = 0.0d0
      Tanalytical(int(nx/2)) = 1.0d0
    endif
    err(:) = Told(:) - Tanalytical(:)
    err(:) = err(:) / maxval(Tanalytical(:))
    
    
    ! Output files written, to be plotted with plot_heat.m
    write(filename,"('HW1_',i6.6,'.dat')") count
    filename = 'figures/'//trim(filename)
    print*,'  plotted: ',trim(filename)
    
    open(unit=11,file=trim(filename),status='unknown',iostat=ier)
    if (ier /= 0) then
      print *,'could not open figure ',trim(filename)
      print *,'Please check if directory figures/ exists...'
      stop 'Error opening file'
    endif
    do i = 1,nx
      ! format: #position #temperature #analytical-solution #error (difference)
      write(11,*) (i-1)*dx, Told(i),Tanalytical(i),err(i)
    enddo
    close(11)

    count = count + 1
  
  endif
  
  end subroutine
