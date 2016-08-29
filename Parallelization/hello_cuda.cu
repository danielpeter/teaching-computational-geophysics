/* 
 
 helloWorld example for CUDA
 
 compile with:
 
 > nvcc -arch=sm_20 hello_cuda.cu
 
 run with:
 
 > ./a.out
 
*/

#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>

#define N 10

// cuda kernel (runs on GPU)

__global__ void sum_kernel(float* A,float* B, float* C, float* sum, int nmax)
{

  // thread id      
  int id = blockIdx.x*blockDim.x + threadIdx.x;
  
  // sums values
  if( id < nmax ){
    C[id] = A[id] + B[id];
    
    // no atomic summation -> wrong results...
    //*sum += C[id];

    // atomic operation to avoid race-conditions
    atomicAdd(sum,C[id]);    
    
  }
}

// main program (runs on CPU)

int main(void)
{

  float *A, *B, *C;
  float *A_d, *B_d, *C_d;

  float sum;
  float *sum_d;

  printf("hello CUDA: \n");
  
  // array on CPU (host)
  A = (float *) malloc(N*sizeof(float));
  B = (float *) malloc(N*sizeof(float));
  C = (float *) malloc(N*sizeof(float));
  
  // initializes
  for(int i=0; i<N; i++){
    A[i] = 1.0;
    B[i] = 2.0;
    C[i] = 0.0;
  }
  
  sum = 0.0f;
  
  // GPU preparation:
  
  // array on GPU (device)
  cudaMalloc((void **) &A_d, N*sizeof(float));
  cudaMalloc((void **) &B_d, N*sizeof(float));
  cudaMalloc((void **) &C_d, N*sizeof(float));
  cudaMalloc((void **) &sum_d, sizeof(float));
  
  // initializes on GPU with zero
  cudaMemset(sum_d,0,sizeof(float));
  
  // copies arrays from CPU to GPU
  cudaMemcpy(A_d,A,N*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(B_d,B,N*sizeof(float),cudaMemcpyHostToDevice);
  
  // cuda kernel dimensions ( 3 blocks x 4 threads )
  int blocksize = 4;
  int nblock = N/blocksize+(N%blocksize==0?0:1);

  // launches cuda kernel  
  sum_kernel<<<nblock,blocksize>>>(A_d,B_d,C_d,sum_d,N);
  
  // copies back from GPU to CPU
  cudaMemcpy(C,C_d,N*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(&sum,sum_d,sizeof(float),cudaMemcpyDeviceToHost);
  
  // user output
  printf("result: \n");  
  for(int i=0;i<N;i++){
    printf("  %f \n",C[i]);
  }  
  printf("\n");
  printf("  sum = %f\n\n",sum);
  
  // releases memory on CPU
  free(A);
  free(B);
  free(C);
  
  // releases memory on GPU
  cudaFree(A_d);   
  cudaFree(B_d);   
  cudaFree(C_d);   
  cudaFree(sum_d);
}
