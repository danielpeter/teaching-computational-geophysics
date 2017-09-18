/* 
  Homework 3

  Program to simulate tsunami waves on a 2D Cartesian grid.  
  The program uses a 2th-order finite difference solution of the equation
  Ptt = div * GH grad P
  where P is the height of the Tsunami wave above sea level, G is the
  acceleration due to gravity (a constant), and H is the ocean depth.
  The speed of the wave is v = sqrt(GH).

  original version from R. Clayton, Caltech, Jan 2005

  Compile as:  
    cc class_tsunami.c
  or:
    gcc -std=c99 class_tsunami.c
  or:
    cc -o class_tsunami class_tsunami.c -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>

// order of finite-difference scheme (2 = 2nd-order, 4 = 4th-order .. not implemented yet!)
int	ord = 2;

// 2D mesh dimensions
int	nx = 1000;
int	ny = 800;

// finite differences
float	h = 10.0;
float	dt = 10.0;

// number of time steps
int	nt = 4000;

// bathymetry
char vmodel[] = "bathy.out.linux";
//char vmodel[] = "bathy.out";

// snapshot output
char output[] = "figures/slices.out";

// time steps between print messages
int	itprint	= 100;
// time steps between slice outputs
int	itslice	= 100;

// station location
float	latref = -40.0;
float	lonref = 35;

// source location
float	slat = 3.30;
float	slon = 95.87;

// minimum depth to consider still ocean (m)
float	mindepth = 10.0;

// mapping from 2-d indexing to 1-d array
#define V(ix,iy)	v[(ix) + nx*(iy)]
#define P1(ix,iy)	p1[(ix) + nx*(iy)]
#define P2(ix,iy)	p2[(ix) + nx*(iy)]

// acc of gravity in km/sec**2
#define G 0.00985

float *v, *p1, *p2;
float *f1, *f2;

int ixref	= 0;
int iyref	= 0;

/* -------------------------------------------------------------- */

// helper functions

/* -------------------------------------------------------------- */

/* convert (lat,lon) to grid coords (ixs,iys) */
#define DEG2KM	111.195
#define DEG2R	  0.0174532

void xcoord_convert(double slat,double slon,int *ixs,int *iys){
  double cos();
  double x, y;

  y = (slat-latref)*DEG2KM;
	x = (slon-lonref)*DEG2KM*cos(slat*DEG2R);
  *ixs = x/h + ixref;
  *iys = y/h + iyref;
}

/* compute the norm of a vector or plane */
double norm(float *x, int n){
  double sum, val, sqrt();

  sum = 0.0;
	for(int i=0; i<n; i++) sum += x[i]*x[i];
  val = sqrt( sum/(double)(n) );
  return(val);
}

/* zero out a field */
void zap(float *x, int n){
  for(int i=0; i<n; i++) x[i] = 0.0;
}

/* find the absolute max of a plane */
double getmax(float *x, int n){
  double fabs(), max;

  max = fabs(x[0]);
  for(int i=1; i<n; i++) if(fabs(x[i]) > max) max = fabs(x[i]);
  return(max);
}

/* output a slice of the field.
 * Note: the field is output at every istep samples
 *       the field is reversed in y
 */
float line[1000];
int outfd	= -1;
int	istep	= 2;

void output_slice(float *x,int nx,int ny,double t){
  int ix,iy,i;

  // file
  if(outfd < 0) outfd = creat(output,0664);
  if(outfd < 0){
    fprintf(stderr,"cannot create plot file = %s\n",output);
    exit(-1);
  }
  // output
  for(iy=ny-1; iy >= 0; iy -= istep){
    for(ix=0, i=0; ix <nx; ix += istep, i++) line[i] = x[ix+iy*nx];
    write(outfd,line,4*i);
  }
}

/* -------------------------------------------------------------- */

// main program

/* -------------------------------------------------------------- */

int main(int ac, char **av){
  int it, ix, iy, fd;
  int ixs, iys;
  float *tmp;
  float vel, velmax, f, val, bathymin, bathymax;
  double norm(), sqrt();

  fprintf(stdout,"Tsunami simulation:\n");
  fprintf(stdout,"finite-difference order = %d\n\n",ord);

  // bathymetry
  v = (float *)(malloc(4*nx*ny));

  // wave-fields
  f1 = (float *)(malloc(4*nx*ny));
  f2 = (float *)(malloc(4*nx*ny));

  if(v == NULL || f1 == NULL || f2 == NULL){
    fprintf(stderr,"cannot alloc memory\n");
    exit(-1);
  }

  // reads in bathymetry file
  fprintf(stdout,"reading bathymetry file = %s\n",vmodel);
	if( (fd = open(vmodel,0)) < 0){
    fprintf(stderr,"cannot open bathymetry file = %s\n",vmodel);
    exit(-1);
  }
	if( read(fd,v,4*nx*ny) != 4*nx*ny ){
    fprintf(stderr,"read error in bathymetry file = %s\n",vmodel);
    exit(-1);
  }
  close(fd);
  output_slice(v,nx,ny,-1.0);

  // bathymetry min/max
  bathymin =  99999.9;
  bathymax = -99999.9;
	for(iy=0; iy<ny; iy++){
    for(ix=0; ix<nx; ix++){
      // elevation (m)
      val = V(ix,iy);
      if (val < bathymin) bathymin = val;
      if (val > bathymax) bathymax = val;
    }
  }
  fprintf(stdout,"bathymetry min/max = %8.4f / %8.4f (m)\n\n",bathymin,bathymax);

	/* convert depth to velocity v = sqrt(g*depth).
	 * set values for land (pos. depths) to negative as flag
	 */
	velmax = 0.0;
	for(iy=0; iy<ny; iy++){
    for(ix=0; ix<nx; ix++){
      // make depth positive
      val = -V(ix,iy);

      // note 0.001 to convert depth (m) to km
      if(val > mindepth){
        vel = sqrt(G*val*0.001);
      }else{
        vel = -0.001;
      }
      if(vel > velmax) velmax = vel;
      if(vel > 0.0){
        V(ix,iy) = vel*vel*dt*dt/(h*h);
      }else{
        V(ix,iy) = -0.001;
      }
    }
  }
  fprintf(stdout,"maximum velocity = %8.4f (km/s)\n",velmax);
  fprintf(stdout,"nx = %d ny = %d nt = %d h = %8.4f dt = %8.4f\n\n",nx,ny,nt,h,dt);

  /* point the memory planes to real memory and zero it */
  p1 = f1;
  p2 = f2;
  zap(p1,nx*ny);
  zap(p2,nx*ny);

  // add source (as initial condition)
  xcoord_convert(slat,slon,&ixs,&iys);
  fprintf(stdout,"source %8.3f %9.3f %4d %4d\n\n",slat,slon,ixs,iys);

  /*  source is placed on a grid:
	 *    1/4  1/2  1/4
	 *    1/2   1   1/2
	 *    1/4  1/2  1/4
	 */
  P2(ixs,iys) = 1.0;
  P2(ixs +1,iys  ) = P2(ixs -1,iys  ) = P2(ixs   ,iys+1) = P2(ixs   ,iys-1) = 0.5;
  P2(ixs+1,iys+1) = P2(ixs-1,iys+1) = P2(ixs+1,iys-1) = P2(ixs-1,iys-1) = 0.25;

  // loop over time steps
  for(it= 0; it<nt; it++){

		// loop over x-y plane
		for(iy=1; iy < ny-1; iy++){
      for(ix=1; ix < nx-1; ix++){

        // ignore points on land
        if(V(ix,iy) < 0.0){
          P1(ix,iy) = 0.0;
          continue;
        }

        // FD scheme
        if(ord == 2 || ix == 1 || ix == nx-2 || iy ==1 || iy == ny-2){
          // 2nd-order
          P1(ix,iy) = 2.0*(1.0 + -2.0*V(ix,iy))*P2(ix,iy) - P1(ix,iy)
                      + 1.0*V(ix,iy)*(P2(ix+1,iy) + P2(ix-1,iy) + P2(ix,iy+1) + P2(ix,iy-1));

        }else{
          //>TODO: implement your 4th-order scheme here
          // 4th-order interior
          fprintf(stderr,"Add 4th order code here\n");
          exit(-1);
          //<TODO
        }
      }
    }
	
		// Dirichlet boundary conditions
		for(ix=0,    iy=0;    ix<nx; ix++) P1(ix,iy) = 0.0;
		for(ix=0,    iy=ny-1; ix<nx; ix++) P1(ix,iy) = 0.0;
		for(ix=0,    iy=0;    iy<ny; iy++) P1(ix,iy) = 0.0;
		for(ix=nx-1, iy=0;    iy<ny; iy++) P1(ix,iy) = 0.0;

    // user output
		if(it%itprint == 0) fprintf(stdout,"done it = %3d norm = %14.3e \n",it,norm(&P1(0,0),nx*ny));

    // time snapshot of wavefield
		if(it%itslice == 0) output_slice(p1,nx,ny,(double)(it*dt));

		// rotate the memory pointers
		tmp = p1; p1 = p2; p2 = tmp;
  }

  return 0;
}

