
#include <cstdio>



// device kernel def
__global__ void Action_noImage_GPU(double *D_,double *maskCenter,double *SolventMols_,double maxD, int Nmols , int NAtoms);

////////////////////////





void Action_NoImage_Center(double *SolventMols_,double *D_, double maskCenter[3],double maxD,int  NMols, int NAtoms, float &time_gpu)
{


  cudaEvent_t start_event, stop_event;
  float elapsed_time_gpu;

  double *devI2Ptr;
  double *devI1Ptr;
  double *devO1Ptr;
  int t4;
  int t2;
  double Dist;
  int solventMol;
  int solventAtom;



  cudaMalloc(((void **)(&devO1Ptr)),NMols * sizeof(double ));
  cudaMalloc(((void **)(&devI1Ptr)),3 * sizeof(double ));
  cudaMemcpy(devI1Ptr,maskCenter,3 * sizeof(double ),cudaMemcpyHostToDevice);
  cudaMalloc(((void **)(&devI2Ptr)),NMols * NAtoms * 3 * sizeof(double ));
  cudaMemcpy(devI2Ptr,SolventMols_,NMols * NAtoms * 3 * sizeof(double ),cudaMemcpyHostToDevice);


  dim3 dimGrid0 = dim3(NMols,1);
  dim3 dimBlock0 = dim3(NAtoms,1);


  printf("NMols =  %d, NAtoms = %d\n", NMols, NAtoms); 
  printf("About to launch kernel.\n");


  cudaEventCreate(&start_event);
  cudaEventCreate(&stop_event);
  cudaEventRecord(start_event, 0);

  Action_noImage_GPU<<<dimGrid0,dimBlock0>>>(devO1Ptr,devI1Ptr, devI2Ptr, maxD, NMols, NAtoms);
  
  cudaThreadSynchronize();
  cudaEventRecord(stop_event, 0);
  cudaEventSynchronize(stop_event);
  cudaEventElapsedTime(&elapsed_time_gpu,start_event, stop_event );


  printf("Done with kernel CUDA Kernel Time: %.2f\n", elapsed_time_gpu);

  time_gpu  = elapsed_time_gpu;
  
  cudaMemcpy(D_,devO1Ptr,NMols * sizeof(double ),cudaMemcpyDeviceToHost);
  cudaFree(devO1Ptr);
  cudaFree(devI1Ptr);
  cudaFree(devI2Ptr);
}
