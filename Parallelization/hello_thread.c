/*
 
 threading in C
 

 compile with:
 
 > cc hello_thread.c
 
 run with:
 
 > ./a.out
 
*/

#include <pthread.h>

#include <stdlib.h>
#include <stdio.h>

#define NUM_THREADS     2

#define N 10


// routine execution within threads

void *do_things(void *threadid)
{
  long tid;

  // thread id
  tid = (long)threadid;  
  printf("  it's me, thread %ld \n", tid);

  float *A, *B, *C;
  float sum;
  int i;
  
  // array allocation
  A = (float *) malloc(N*sizeof(float));
  B = (float *) malloc(N*sizeof(float));
  C = (float *) malloc(N*sizeof(float));
  
  // initializes
  for(i=0; i<N; i++){
    A[i] = 1.0;
    B[i] = 2.0;
    C[i] = 0.0;
  }
  sum = 0.0f;
  
  // only master process calculates, others do whatever
  if( tid == 0 ){
    // some calculations
    for(i=0; i<N; i++){
      C[i] = A[i] + B[i];
      sum += C[i];  
    }  
  }

  // frees memory
  free(A);
  free(B);
  free(C);
  
  printf("  thread %ld result: %f \n",tid,sum);
  pthread_exit(NULL);
}


// main program routine

int main (int argc, char *argv[])
{

  pthread_t threads[NUM_THREADS];
  int rc;
  long t;
  
  printf("Hello world\n");
  
  // starts multiple threads
  for(t=0; t<NUM_THREADS; t++){
  
    // creates threads    
    rc = pthread_create(&threads[t], NULL, do_things, (void *)t);
    if (rc){
      printf("ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }
  
  // Last thing that main() should do
  pthread_exit(NULL);
}

