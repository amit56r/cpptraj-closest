#include <stdio.h>

#define BLOCKDIM 1024
#define RSIZE 1024
#define C_FACTOR 4

//------------------------------------------------------------------------------------------------------------------------------------------------
//try thread coarsening 

__global__ void Action_noImage_center_GPU(double *D_,double *maskCenter,double *SolventMols_,double maxD, int Nmols , int NAtoms, int active_size)
{

	__shared__ double dist_array[BLOCKDIM];



	int mol  =  (blockIdx.x * active_size + threadIdx.x)/NAtoms; 
	int atom  = (blockIdx.x * active_size + threadIdx.x) - (mol * NAtoms);
	//int mol_in_block = threadIdx.x/NAtoms;

	//advantage of register
	double a0 = maskCenter[0];
	double a1 = maskCenter[1];
	double a2 = maskCenter[2];



	if ( threadIdx.x < active_size && mol*NAtoms + atom < Nmols*NAtoms )
	{

		// if(atom == 0 )
		// 	D_[mol] = maxD;
		//__syncthreads();


		int sIndex =  mol*NAtoms*3 + atom*3;

		double x =  a0 - SolventMols_[sIndex + 0];
		double y = a1 - SolventMols_[sIndex + 1];
		double z =  a2 - SolventMols_[sIndex + 2];
	//Dist = x*x + y*y + z*z;
		dist_array[threadIdx.x] = x*x + y*y + z*z;
	//printf(" dist  =  %f\n", Dist);

		__syncthreads();

	//first thread
	//naive approach to a reduction algorithm
	//this works if NAtoms is small other wise you need split
	//and do some of log(n) parallel reduction 
		int i;



		double min_val  = maxD;
		if( atom ==0 )
		{
			for(i  = 0 ; i < NAtoms ; i++ ){
				//sIndex = mol*NAtoms*3 + i*3;
				//if (dist_array[threadIdx.x + i]  < min_val) 
				//	min_val = dist_array[threadIdx.x + i] ;
				min_val =  min(min_val, dist_array[threadIdx.x + i]);
			}
			D_[mol] = min_val;
		}

	//if(tx == 0 && bx == 0 )
	//	printf("end of kernel");
	}
}


	// int i;
	// 	double min_val  = maxD;
	// 	if( atom ==0 )
	// 	{
	// 		for(i  = 0 ; i < NAtoms ; i++ ){
	// 			//sIndex = mol*NAtoms*3 + i*3;
	// 			if (dist_array[threadIdx.x + i]  < min_val) 
	// 				min_val = dist_array[threadIdx.x + i] ;
	// 		}
	// 		D_[mol] = min_val;
	// 	}


	// double min_val  = maxD;
	// if( threadIdx.x < active_size/NAtoms )
	// {

	// 	for(i  = threadIdx.x*NAtoms ; i <threadIdx.x*NAtoms + NAtoms ; i++ ){
	// 		//sIndex = mol*NAtoms*3 + i*3;
	// 		if (dist_array[i]  < min_val) 
	// 			min_val = dist_array[i] ;
	// 	}
	// 	D_[blockIdx.x * active_size/NAtoms + threadIdx.x  ] = min_val;
	// }

//------------------------------------------------------------------------------------------------------------------------------------------------

__global__ void Action_noImage_no_center_GPU(double *D_,double *SolventMols_,double *Solute_atoms ,double maxD, int Nmols , int NAtoms,int NSAtoms , int active_size)
{

	__shared__ double dist_array[BLOCKDIM];
	__shared__ double sAtom_shared[RSIZE];




	int mol  =  (blockIdx.x * active_size + threadIdx.x)/NAtoms; 
	int atom  = (blockIdx.x * active_size + threadIdx.x) - (mol * NAtoms);
	//int mol_in_block = threadIdx.x/NAtoms;

	

	//handling the chunks for  solute_atoms
	int chunksize,start,end, NChunks,i,j;

	if(NSAtoms*3 > RSIZE)
	{
		chunksize = (RSIZE/3)*3;
		NChunks = ceil(double(NSAtoms*3)/chunksize);
		start = 0;
		end = chunksize;
	}
	else
	{
		chunksize = NSAtoms*3;
		NChunks = 1;
		start = 0;
		end = NSAtoms*3;
	}

	// if(threadIdx.x == 0 && blockIdx.x == 0 )
	// 	printf("chunkszize = %d ; Nchunk =  %d; start = %d; end = %d\n ",
	// 		chunksize,NChunks,start,end);



	if ( threadIdx.x < active_size && mol*NAtoms + atom < Nmols*NAtoms )
	{

		// if(atom == 0 )
		// 	D_[mol] = maxD;
		//__syncthreads(); 
		double min_val  = maxD;
		double dist;
		int sIndex =  mol*NAtoms*3 + atom*3;
		double a0 = SolventMols_[sIndex + 0];
		double a1 = SolventMols_[sIndex + 1];
		double a2 = SolventMols_[sIndex + 2];


		for(i  = 0 ; i  < NChunks ; i++)
		{
			//copying to shared
			//if (threadIdx.x < (end - start))
			//	sAtom_shared[threadIdx.x] = Solute_atoms[start + threadIdx.x];

			//__syncthreads();

			//TODO - add skew per thread 
			for (j = start ; j < end; j+=3 )
			{
				//int offset = start + (j + threadIdx.x)%(end - start);
				double x = Solute_atoms[j + 0]  - a0;
				double y = Solute_atoms[j + 1]  - a1;
				double z = Solute_atoms[j + 2]  - a2;
				dist =  x*x + y*y + z*z;
				//if (mol ==  11)
				//	printf("min  = %f\n",min_val);
				min_val = min(min_val,dist);


			}

			start = end;
			end = min(end + chunksize, NSAtoms*3);


		}

		dist_array[threadIdx.x] = min_val;
		//if (threadIdx.x == 0)
		//	printf("min_val  = %f\n",min_val);
	//printf(" dist  =  %f\n", Dist);

		__syncthreads();

	//first thread
	//naive approach to a reduction algorithm
	//this works if NAtoms is small other wise you need split
	//and do some of log(n) parallel reduction 
		//min_val  = maxD;
		if( atom ==0 )
		{
			for(i  = 0 ; i < NAtoms ; i++ ){
				//sIndex = mol*NAtoms*3 + i*3;
				//if (dist_array[threadIdx.x + i]  < min_val) 
				//	min_val = dist_array[threadIdx.x + i] ;
				min_val =  min(min_val, dist_array[threadIdx.x + i]);
			}
			D_[mol] = min_val;
		}

	//if(tx == 0 && bx == 0 )
	//	printf("end of kernel");
	}
}





//------------------------------------------------------------------------------------------------------------------------------------------------
__global__ void Action_ImageOrtho_center_GPU(double *D_,double *maskCenter,double *SolventMols_,double maxD, double *box, int Nmols , int NAtoms, int active_size)
{
	__shared__ double dist_array[BLOCKDIM];



	int mol  =  (blockIdx.x * active_size + threadIdx.x)/NAtoms; 
	int atom  = (blockIdx.x * active_size + threadIdx.x) - (mol * NAtoms);
	//int mol_in_block = threadIdx.x/NAtoms;

	//advantage of register
	double a0 = maskCenter[0];
	double a1 = maskCenter[1];
	double a2 = maskCenter[2];



	if ( threadIdx.x < active_size && mol*NAtoms + atom < Nmols*NAtoms )
	{

		// if(atom == 0 )
		// 	D_[mol] = maxD;
		//__syncthreads();


		int sIndex =  mol*NAtoms*3 + atom*3;

		double x =  a0 - SolventMols_[sIndex + 0];
		double y =  a1 - SolventMols_[sIndex + 1];
		double z =  a2 - SolventMols_[sIndex + 2];

		// Get rid of sign info
		if (x<0) x=-x;
		if (y<0) y=-y;
		if (z<0) z=-z;
		  // Get rid of multiples of box lengths 
		//TODO  WIERD that should be a way to simplify it
		while (x > box[0]) x = x - box[0];
		while (y > box[1]) y = y - box[1];
		while (z > box[2]) z = z - box[2];
		  // Find shortest distance in periodic reference
		double D = box[0] - x;
		if (D < x) x = D;
		D = box[1] - y;
		if (D < y) y = D;  
		D = box[2] - z;
		if (D < z) z = D;




	//Dist = x*x + y*y + z*z;
		dist_array[threadIdx.x] = x*x + y*y + z*z;
		if (box[0]==0.0 || box[1]==0.0 || box[2]==0.0)
			dist_array[threadIdx.x] = -1.0;
	//printf(" dist  =  %f\n", Dist);

		__syncthreads();

	//first thread
	//naive approach to a reduction algorithm
	//this works if NAtoms is small other wise you need split
	//and do some of log(n) parallel reduction 
		int i;



		double min_val  = maxD;
		if( atom ==0 )
		{
			for(i  = 0 ; i < NAtoms ; i++ ){
				//sIndex = mol*NAtoms*3 + i*3;
				//if (dist_array[threadIdx.x + i]  < min_val) 
				//	min_val = dist_array[threadIdx.x + i] ;
				min_val =  min(min_val, dist_array[threadIdx.x + i]);
			}
			D_[mol] = min_val;
		}

	//if(tx == 0 && bx == 0 )
	//	printf("end of kernel");
	}
}

//------------------------------------------------------------------------------------------------------------------------------------------------
__global__ void Action_ImageOrtho_no_center_GPU(double *D_,double *SolventMols_,double *Solute_atoms ,double maxD, double *box, int Nmols , int NAtoms,int NSAtoms , int active_size)
{

	__shared__ double dist_array[BLOCKDIM];
	__shared__ double sAtom_shared[RSIZE];




	int mol  =  (blockIdx.x * active_size + threadIdx.x)/NAtoms; 
	int atom  = (blockIdx.x * active_size + threadIdx.x) - (mol * NAtoms);
	//int mol_in_block = threadIdx.x/NAtoms;

	

	//handling the chunks for  solute_atoms
	int chunksize,start,end, NChunks,i,j;

	if(NSAtoms*3 > RSIZE)
	{
		chunksize = (RSIZE/3)*3;
		NChunks = ceil(double(NSAtoms*3)/chunksize);
		start = 0;
		end = chunksize;
	}
	else
	{
		chunksize = NSAtoms*3;
		NChunks = 1;
		start = 0;
		end = NSAtoms*3;
	}

	// if(threadIdx.x == 0 && blockIdx.x == 0 )
	// 	printf("chunkszize = %d ; Nchunk =  %d; start = %d; end = %d\n ",
	// 		chunksize,NChunks,start,end);



	if ( threadIdx.x < active_size && mol*NAtoms + atom < Nmols*NAtoms )
	{

		// if(atom == 0 )
		// 	D_[mol] = maxD;
		//__syncthreads(); 
		double min_val  = maxD;
		double dist;
		int sIndex =  mol*NAtoms*3 + atom*3;
		double a0 = SolventMols_[sIndex + 0];
		double a1 = SolventMols_[sIndex + 1];
		double a2 = SolventMols_[sIndex + 2];


		for(i  = 0 ; i  < NChunks ; i++)
		{
			//copying to shared
			//if (threadIdx.x < (end - start))
			//	sAtom_shared[threadIdx.x] = Solute_atoms[start + threadIdx.x];

			//__syncthreads();

			//TODO - add skew per thread 
			for (j = start ; j < end; j+=3 )
			{
				//int offset = start + (j + threadIdx.x)%(end - start);
				double x = Solute_atoms[j + 0]  - a0;
				double y = Solute_atoms[j + 1]  - a1;
				double z = Solute_atoms[j + 2]  - a2;


				// Get rid of sign info
				if (x<0) x=-x;
				if (y<0) y=-y;
				if (z<0) z=-z;
		  		// Get rid of multiples of box lengths 
				//TODO  WIERD that should be a way to simplify it
				while (x > box[0]) x = x - box[0];
				while (y > box[1]) y = y - box[1];
				while (z > box[2]) z = z - box[2];

				//below is actually slower! 
				//x = x - box[0]*((int)x/box[0]);
				//y = y - box[0]*((int)y/box[1]);
				//z = z - box[0]*((int)z/box[2]);
		  	// Find shortest distance in periodic reference
				double D = box[0] - x;
				if (D < x) x = D;
				D = box[1] - y;
				if (D < y) y = D;  
				D = box[2] - z;
				if (D < z) z = D;


				//Dist = x*x + y*y + z*z;
				dist =  x*x + y*y + z*z;
				if (box[0]==0.0 || box[1]==0.0 || box[2]==0.0)
					dist = -1.0;


		


				//if (mol ==  11)
				//	printf("min  = %f\n",min_val);
				min_val = min(min_val,dist);


			}

			start = end;
			end = min(end + chunksize, NSAtoms*3);


		}

		dist_array[threadIdx.x] = min_val;
		//if (threadIdx.x == 0)
		//	printf("min_val  = %f\n",min_val);
	//printf(" dist  =  %f\n", Dist);

		__syncthreads();

	//first thread
	//naive approach to a reduction algorithm
	//this works if NAtoms is small other wise you need split
	//and do some of log(n) parallel reduction 
		//min_val  = maxD;
		if( atom ==0 )
		{
			for(i  = 0 ; i < NAtoms ; i++ ){
				//sIndex = mol*NAtoms*3 + i*3;
				//if (dist_array[threadIdx.x + i]  < min_val) 
				//	min_val = dist_array[threadIdx.x + i] ;
				min_val =  min(min_val, dist_array[threadIdx.x + i]);
			}
			D_[mol] = min_val;
		}

	//if(tx == 0 && bx == 0 )
	//	printf("end of kernel");
	}
}





