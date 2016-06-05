
#include <cstdio>
#include <cmath>

#define BLOCKDIM 1024


// device kernel def
__global__ void Action_noImage_center_GPU(double *D_,double *maskCenter,double *SolventMols_,double maxD, int Nmols , int NAtoms, int active_size);
__global__ void Action_noImage_no_center_GPU(double *D_,double *SolventMols_,double *Solute_atoms ,double maxD, int Nmols , int NAtoms,int NSAtoms , int active_size);


//for imaging with ortho
__global__ void Action_ImageOrtho_center_GPU(double *D_,double *maskCenter,double *SolventMols_,double maxD, double *box, int Nmols , int NAtoms, int active_size);
__global__ void Action_ImageOrtho_no_center_GPU(double *D_,double *SolventMols_,double *Solute_atoms ,double maxD, double *box, int Nmols , int NAtoms,int NSAtoms , int active_size);

//for imaging with NONortho
//TODO

////////////////////////





void Action_NoImage_Center(double *SolventMols_,double *D_, double maskCenter[3],double maxD,int  NMols, int NAtoms, float &time_gpu,int type, double box[3])
{


  cudaEvent_t start_event, stop_event;
  float elapsed_time_gpu;

  double *devI2Ptr;
  double *devI1Ptr;
  double *devO1Ptr;
  double *boxDev;
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

  cudaMalloc(((void**)(&boxDev)), 3 * sizeof(double));
  cudaMemcpy(boxDev,box, 3 * sizeof(double), cudaMemcpyHostToDevice);



  //figue out the decomposition here
  //we need to pad as well

 // due to lack to  using center, each thread is going  rocess the solvent mol 
 //instead of atoms (make it alot easier)   (speacially for the imaging case)




  //figure out how many active thread in a block
  int active_size  =  BLOCKDIM/NAtoms * NAtoms;
  //int NBlocks =  ceil(NMols * NAtoms / float(active_size));  //having unroll factor
  int NBlocks = ceil(float(NMols)/ BLOCKDIM);

  // printf("Nmols = %d; Natoms = %d\n", NMols, NAtoms);
  // printf("active_size =  %d\n", active_size);
  //printf("NBlocks =  %d\n", NBlocks);
  //printf("sezeof(double) = %d\n", sizeof(double));
  //exit(0);



  dim3 dimGrid0 = dim3(NBlocks,1);
  dim3 dimBlock0 = dim3(BLOCKDIM,1);


  printf("NMols =  %d, NAtoms = %d\n", NMols, NAtoms); 
  printf("About to launch kernel.\n");


  cudaEventCreate(&start_event);
  cudaEventCreate(&stop_event);
  cudaEventRecord(start_event, 0);

  if(type == 0)
    Action_noImage_center_GPU<<<dimGrid0,dimBlock0>>>(devO1Ptr,devI1Ptr, devI2Ptr, maxD, NMols, NAtoms,active_size);
  else if (type == 1)
    Action_ImageOrtho_center_GPU<<<dimGrid0,dimBlock0>>>(devO1Ptr,devI1Ptr, devI2Ptr, maxD,boxDev, NMols, NAtoms,active_size);
  else
    printf("kernel_wrapper: error in type\n");

  
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
  cudaFree(boxDev);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Action_NoImage_no_Center(double *SolventMols_,double *D_, double *Solute_atoms,double maxD,int  NMols, int NAtoms,int NSAtoms, float &time_gpu, int type,double box[3])
{


  cudaEvent_t start_event, stop_event;
  float elapsed_time_gpu;

  double *devI3Ptr;
  double *devI2Ptr;
  double *devI1Ptr;
  double *devO1Ptr;
  double *boxDev;
  int t4;
  int t2;
  double Dist;
  int solventMol;
  int solventAtom;



  cudaMalloc(((void **)(&devO1Ptr)),NMols * sizeof(double ));

  //cudaMalloc(((void **)(&devI1Ptr)),3 * sizeof(double ));
  //cudaMemcpy(devI1Ptr,maskCenter,3 * sizeof(double ),cudaMemcpyHostToDevice);
  cudaMalloc(((void **)(&devI2Ptr)),NMols * NAtoms * 3 * sizeof(double ));
  cudaMemcpy(devI2Ptr,SolventMols_,NMols * NAtoms * 3 * sizeof(double ),cudaMemcpyHostToDevice);
  
  cudaMalloc(((void **)(&devI3Ptr)), NSAtoms * 3 * sizeof(double ));
  cudaMemcpy(devI3Ptr,Solute_atoms,NSAtoms * 3 * sizeof(double ),cudaMemcpyHostToDevice);

  cudaMalloc(((void**)(&boxDev)), 3 * sizeof(double));
  cudaMemcpy(boxDev,box, 3 * sizeof(double), cudaMemcpyHostToDevice);



  //figue out the decomposition here
  //we need to pad as well

  //figure out how many active thread in a block
  int active_size  =  BLOCKDIM/NAtoms * NAtoms;
  int NBlocks =  ceil(NMols * NAtoms / float(active_size));
  // printf("Nmols = %d; Natoms = %d\n", NMols, NAtoms);
  // printf("active_size =  %d\n", active_size);
  // printf("NBlocks =  %d\n", NBlocks);
  //printf("sezeof(double) = %d\n", sizeof(double));
  //exit(0);



  dim3 dimGrid0 = dim3(NBlocks,1);
  dim3 dimBlock0 = dim3(BLOCKDIM,1);


  printf("NMols =  %d, NAtoms = %d\n", NMols, NAtoms); 
  printf("About to launch kernel.\n");


  cudaEventCreate(&start_event);
  cudaEventCreate(&stop_event);
  cudaEventRecord(start_event, 0);

  if(type == 0)
    Action_noImage_no_center_GPU<<<dimGrid0,dimBlock0>>>(devO1Ptr, devI2Ptr,devI3Ptr, maxD, NMols, NAtoms,NSAtoms,active_size);
  else if(type == 1)
    Action_ImageOrtho_no_center_GPU<<<dimGrid0,dimBlock0>>>(devO1Ptr, devI2Ptr,devI3Ptr, maxD, boxDev,  NMols, NAtoms,NSAtoms,active_size);
  else
    printf("kernel_wrapper: error in type no center version\n");
  
  cudaThreadSynchronize();
  cudaEventRecord(stop_event, 0);
  cudaEventSynchronize(stop_event);
  cudaEventElapsedTime(&elapsed_time_gpu,start_event, stop_event );


  printf("Done with kernel CUDA Kernel Time: %.2f\n", elapsed_time_gpu);

  time_gpu  = elapsed_time_gpu;
  
  cudaMemcpy(D_,devO1Ptr,NMols * sizeof(double ),cudaMemcpyDeviceToHost);
  cudaFree(devO1Ptr);
  //cudaFree(devI1Ptr);
  cudaFree(devI2Ptr);
  cudaFree(devI3Ptr);
  cudaFree(boxDev);
}