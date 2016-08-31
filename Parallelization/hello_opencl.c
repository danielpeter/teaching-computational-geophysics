/* 
 
 helloWorld example for OpenCL
 
 compile with:
 
 MAC OsX
 > cc hello_opencl.c -framework OpenCL
 
 Linux: (needs OpenCL libraries installed)
 > cc hello_opencl.c -lOpenCL

 run with:
 
 > ./a.out
 
*/

#include <stdio.h>
#include <stdlib.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif


#define N 10

// OpenCL kernel (runs on GPU)
const char * sum_kernel_program = "\
#pragma OPENCL EXTENSION cl_khr_fp64: enable\n\
inline void atomicAdd(volatile __global float *source, const float val) {\n\
  union {\n\
    unsigned int iVal;\n\
    float fVal;\n\
  } res, orig;\n\
  do {\n\
    orig.fVal = *source;\n\
    res.fVal = orig.fVal + val;\n\
  } while (atomic_cmpxchg((volatile __global unsigned int *)source, orig.iVal, res.iVal) != orig.iVal);\n\
}\n\
\n\
__kernel void sum_kernel(__global float* A, __global float* B, __global float* C, __global float* sum, const int nmax){\n\
  int id = get_global_id(0) + (get_group_id(1)) * (get_global_size(0));\n\
  if (id < nmax) {\n\
    C[id] = A[id] + B[id];\n\
    // no atomic summation -> wrong results...\n\
    //*sum += C[id];\n\
    \n\
    // atomic operation to avoid race-conditions\n\
    atomicAdd(sum,C[id]);\n\
  }\n\
}\n\
";


// main program (runs on CPU)

int main(void)
{

  float *A, *B, *C;
  float sum;
  int err;

  printf("hello OpenCL: \n");
  
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
  
  // OpenCL GPU preparation:
  cl_uint num_platforms,num_devices;
  cl_device_id device_id;
  char info[512];

  cl_context context;
  cl_command_queue commands;
  cl_program program;
  cl_kernel kernel;

  cl_mem A_d,B_d,C_d,sum_d;

  // first OpenCL call
  // only gets number of platforms
  err = clGetPlatformIDs(0, NULL, &num_platforms);
  if (err != CL_SUCCESS){ printf("Error: Failed to get platform ids!\n"); return EXIT_FAILURE;}
  printf("  number of OpenCL platforms = %d\n",num_platforms);

  // checks if OpenCL platforms available
  if (num_platforms == 0) {
    fprintf(stderr,"OpenCL error: No OpenCL platform available!\n");
    exit(1);
  }else{
    // just shows what devices would be available
    cl_platform_id *platform_ids;
    cl_device_id *devices;

    platform_ids = (cl_platform_id *) malloc(num_platforms * sizeof(cl_platform_id));

    // gets platform infos
    err = clGetPlatformIDs(num_platforms, platform_ids, NULL);
    if (err != CL_SUCCESS){ printf("Error: Failed to get platform ids!\n"); return EXIT_FAILURE;}

    // loops over available devices
    for (int i = 0; i < num_platforms; i++) {
      printf("  platform:\n");
      err = clGetPlatformInfo(platform_ids[i], CL_PLATFORM_NAME, sizeof(info), info, NULL);
      if (err != CL_SUCCESS){ printf("Error: Failed to get platform ids!\n"); return EXIT_FAILURE;}
      printf("    platform name = %s\n",info);

      // gets number of devices for this platform
      err = clGetDeviceIDs(platform_ids[i], CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices);
      if (err != CL_SUCCESS){ printf("Error: Failed to get platform devices!\n"); return EXIT_FAILURE;}

      printf("    platform number of GPU devices = %d\n",num_devices);
      if (num_devices > 0) {
        // fills device infos
        devices = (cl_device_id *) malloc(num_devices * sizeof(cl_device_id));

        // gets device infos
        err = clGetDeviceIDs(platform_ids[i], CL_DEVICE_TYPE_GPU, num_devices, devices, NULL);
        if (err != CL_SUCCESS){ printf("Error: Failed to get devices!\n"); return EXIT_FAILURE;}

        // loops over all devices
        for (int j = 0; j < num_devices; j++) {
          printf("    device %i:\n",j);
          err = clGetDeviceInfo(devices[j], CL_DEVICE_NAME, sizeof(info), info, NULL);
          printf("      device Name = %s\n", info);
          err = clGetDeviceInfo(devices[j], CL_DEVICE_VENDOR, sizeof(info), info, NULL);
          printf("      device Vendor = %s\n", info);
        }
      }
    }
  }

  // initializes first available OpenCL device
  //err = clGetDeviceIDs(NULL, 1? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, 1, &device_id, NULL);
  err = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, 1, &device_id, NULL);
  if (err != CL_SUCCESS){ printf("Error: Failed to get device ids!\n"); return EXIT_FAILURE;}

  err = clGetDeviceInfo(device_id, CL_DEVICE_NAME, sizeof(info), info, NULL);
  if (err != CL_SUCCESS){ printf("Error: Failed to get device info!\n"); return EXIT_FAILURE;}
  printf("running on device: name = %s\n",info);

  // creates an OpenCL context
  context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
  if (!context){ printf("Error: Failed to create context!\n"); return EXIT_FAILURE; }

  // command kernel queues
  // clCreateCommandQueue feature in OpenCL 1.2, will be deprecated in OpenCL 2.0
#ifdef CL_VERSION_2_0
  // version 2.0
  commands = clCreateCommandQueueWithProperties(context, device_id, 0, &err);
#else
  // version 1.2, CL_VERSION_1_2
  commands = clCreateCommandQueue(context, device_id, 0, &err);
#endif
  if (!commands){ printf("Error: Failed to create commands!\n"); return EXIT_FAILURE; }
 
  // creates the compute program
  program = clCreateProgramWithSource(context, 1, (const char **) &sum_kernel_program, NULL, &err);
  if (!program){ printf("Error: Failed to create compute program!\n"); return EXIT_FAILURE; }
 
  // builds the program executable
  err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
  if (err != CL_SUCCESS){
    printf("Error: Failed to build program!\n");
    // returns error message
    size_t len;
    char buffer[2048];
    clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
    printf("%s\n", buffer);
    exit(1);
  }

  // creates the compute kernel
  kernel = clCreateKernel(program, "sum_kernel", &err);
  if (!kernel || err != CL_SUCCESS){ printf("Error: Failed to create compute kernel!\n"); exit(1); }
 
  // array on GPU (device)
  A_d = clCreateBuffer(context,  CL_MEM_READ_ONLY,  N*sizeof(float), NULL, NULL);
  B_d = clCreateBuffer(context,  CL_MEM_READ_ONLY,  N*sizeof(float), NULL, NULL);
  C_d = clCreateBuffer(context, CL_MEM_READ_WRITE, N*sizeof(float), NULL, NULL);
  sum_d = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float), NULL, NULL);
  if (!A_d || !B_d || !C_d || !sum_d){ printf("Error: Failed to allocate device memory!\n"); exit(1); }

  // copies arrays from CPU to GPU
  err = clEnqueueWriteBuffer(commands, A_d, CL_TRUE, 0, N*sizeof(float), A, 0, NULL, NULL);
  if (err != CL_SUCCESS){ printf("Error: Failed to write to source array A!\n"); exit(1); }
  err = clEnqueueWriteBuffer(commands, B_d, CL_TRUE, 0, N*sizeof(float), B, 0, NULL, NULL);
  if (err != CL_SUCCESS){ printf("Error: Failed to write to source array B!\n"); exit(1); }
  // initializes on GPU with zero
  err = clEnqueueWriteBuffer(commands, C_d, CL_TRUE, 0, N*sizeof(float), C, 0, NULL, NULL);
  if (err != CL_SUCCESS){ printf("Error: Failed to write to source array C!\n"); exit(1); }
  err = clEnqueueWriteBuffer(commands, sum_d, CL_TRUE, 0, sizeof(float), &sum, 0, NULL, NULL);
  if (err != CL_SUCCESS){ printf("Error: Failed to write to source sum!\n"); exit(1); }

  // OpenCL kernel dimensions ( 3 blocks x 4 threads )
  size_t global_work_size;
  size_t local_work_size;

  int nmax = N;
  int blocksize = 4;
  int nblock = N/blocksize+(N%blocksize==0?0:1);

  local_work_size = blocksize;
  global_work_size = nblock * blocksize;

  // sets OpenCL kernel arguments
  err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &A_d);
  err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &B_d);
  err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &C_d);
  err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &sum_d);
  err |= clSetKernelArg(kernel, 4, sizeof(int), &nmax);
  if (err != CL_SUCCESS){ printf("Error: Failed to set kernel arguments! %d\n", err); exit(1); }

  // launches OpenCL kernel
  err = clEnqueueNDRangeKernel(commands, kernel, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL);
  if (err != CL_SUCCESS){ printf("Error: Failed to execute kernel!\n"); return EXIT_FAILURE; }

  // wait for finish
  clFinish(commands);

  // copies back from GPU to CPU
  err = clEnqueueReadBuffer(commands, C_d, CL_TRUE, 0, N*sizeof(float), C, 0, NULL, NULL);
  if (err != CL_SUCCESS){ printf("Error: Failed to read output! %d\n",err); return EXIT_FAILURE; }
  err = clEnqueueReadBuffer(commands, sum_d, CL_TRUE, 0, sizeof(float), &sum, 0, NULL, NULL);
  if (err != CL_SUCCESS){ printf("Error: Failed to read output! %d\n",err); return EXIT_FAILURE; }

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
  clReleaseMemObject(A_d);
  clReleaseMemObject(B_d);
  clReleaseMemObject(C_d);
  clReleaseMemObject(sum_d);

  clReleaseProgram(program);
  clReleaseKernel(kernel);
  clReleaseCommandQueue(commands);
  clReleaseContext(context);

  return 0;
}
