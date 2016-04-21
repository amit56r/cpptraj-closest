#include <stdio.h>




__global__ void Action_noImage_GPU(double *D_,double *maskCenter,double *SolventMols_,double maxD, int Nmols , int NAtoms)
{
	int bx;
	bx = blockIdx.x;
	int tx;
	tx = threadIdx.x;
	double Dist;

	if(tx == 0 )
		D_[bx] = maxD;
	__syncthreads();


	int sIndex =  bx*NAtoms*3 + tx*3;

	double x =  maskCenter[0] - SolventMols_[sIndex + 0];
	double y = maskCenter[1] - SolventMols_[sIndex + 1];
	double z =  maskCenter[2] - SolventMols_[sIndex + 2];
	//Dist = x*x + y*y + z*z;
	SolventMols_[sIndex] = x*x + y*y + z*z;
	//printf(" dist  =  %f\n", Dist);

	__syncthreads();

	//first thread
	//naive approach to a reduction algorithm
	//this works if NAtoms is small other wise you need split
	//and do some of log(n) parallel reduction 
	int i;
	if( tx ==0 )
	{
		for(i  = 0 ; i < NAtoms ; i++ ){
			sIndex = bx*NAtoms*3 + i*3;
			if (SolventMols_[sIndex]  < D_[bx]) 
				D_[bx] = SolventMols_[sIndex];
		}
	}

	//if(tx == 0 && bx == 0 )
	//	printf("end of kernel");
}
