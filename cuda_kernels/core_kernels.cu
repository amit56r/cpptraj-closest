#include <stdio.h>




__global__ void Action_noImage_center_GPU(double *D_,double *maskCenter,double *SolventMols_,double maxD, int Nmols , int NAtoms, int active_size)
{

	double Dist;

	int mol  =  (blockIdx.x * active_size + threadIdx.x)/NAtoms; 
	int atom  = (blockIdx.x * active_size + threadIdx.x) - (mol * NAtoms);

	double a0 = maskCenter[0];
	double a1 = maskCenter[1];
	double a2 = maskCenter[2];



	if ( threadIdx.x < active_size && mol*NAtoms + atom < Nmols*NAtoms )
	{

		if(atom == 0 )
			D_[mol] = maxD;
		//__syncthreads();


		int sIndex =  mol*NAtoms*3 + atom*3;

		double x =  a0 - SolventMols_[sIndex + 0];
		double y = a1 - SolventMols_[sIndex + 1];
		double z =  a2 - SolventMols_[sIndex + 2];
	//Dist = x*x + y*y + z*z;
		SolventMols_[sIndex] = x*x + y*y + z*z;
	//printf(" dist  =  %f\n", Dist);

		__syncthreads();

	//first thread
	//naive approach to a reduction algorithm
	//this works if NAtoms is small other wise you need split
	//and do some of log(n) parallel reduction 
		int i;
		if( atom ==0 )
		{
			for(i  = 0 ; i < NAtoms ; i++ ){
				sIndex = mol*NAtoms*3 + i*3;
				if (SolventMols_[sIndex]  < D_[mol]) 
					D_[mol] = SolventMols_[sIndex];
			}
		}

	//if(tx == 0 && bx == 0 )
	//	printf("end of kernel");
	}
}
