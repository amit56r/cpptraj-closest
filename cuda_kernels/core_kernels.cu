#include <stdio.h>




__global__ void Action_noImage_GPU(double *D_,double *maskCenter,double *SolventMols_,double maxD, int Nmols , int NAtoms)
{
  int bx;
  bx = blockIdx.x;
  int tx;
  tx = threadIdx.x;
  double Dist;
  int t2;
  int t4;

  if(tx == 0 && bx == 0)
  	D_[bx] = maxD;
  __syncthreads();


  int sIndex =  bx*NAtoms*3 + tx*3;

  double x =  maskCenter[0] - SolventMols_[sIndex + 0];
  double y = maskCenter[1] - SolventMols_[sIndex + 1];
  double z =  maskCenter[2] - SolventMols_[sIndex + 2];
  Dist = x*x + y*y + z*z;
  if (Dist  < D_[bx]) 
    D_[bx] = Dist;

  if(tx == 0 && bx == 0 )
	printf("end of kernel");
}
