#define __rose_lt(x,y) ((x)<(y)?(x):(y))
#define __rose_gt(x,y) ((x)>(y)?(x):(y))
#define D__(solventMol) D_[solventMol]
__global__ void Action_noImage_GPU(double *D_,double *maskCenter,double (*SolventMols_)[1024][3]);
//this is only  used for cuda-chill
//heavy simplification
#define NsolventMolecules_ 1024
#define NsolventAtoms_ 1024
// struct MolDist {
//       int mol;        ///< Original solvent molecule number (starts from 1).
//       double D;       ///< Closest distance of solvent molecule to atoms in distanceMask.
//       //AtomMask mask;  ///< Original topology solvent molecule atom mask.
//       double solventAtoms[NsolventAtoms_][3]; ///< Actual solvent atom #s to loop over.
//   };
//using dist for no image 
// and kernel for when we use solute molecule center
//extracting pulling out arrays out from struct 

void Action_NoImage_Center(double SolventMols_[1024][1024][3],double D_[1024],double maskCenter[3],double maxD)
{
  double *devI2Ptr;
  double *devI1Ptr;
  double *devO1Ptr;
  int t4;
  int t2;
  double Dist;
  int solventMol;
  int solventAtom;
  cudaMalloc(((void **)(&devO1Ptr)),3072 * sizeof(double ));
  cudaMalloc(((void **)(&devI1Ptr)),3 * sizeof(double ));
  cudaMemcpy(devI1Ptr,maskCenter,3 * sizeof(double ),cudaMemcpyHostToDevice);
  cudaMalloc(((void **)(&devI2Ptr)),3145728 * sizeof(double ));
  cudaMemcpy(devI2Ptr,SolventMols_,3145728 * sizeof(double ),cudaMemcpyHostToDevice);
  dim3 dimGrid0 = dim3(1024,1);
  dim3 dimBlock0 = dim3(1024,1);
  Action_noImage_GPU<<<dimGrid0,dimBlock0>>>(devO1Ptr,devI1Ptr,((double (*)[1024][3])devI2Ptr));
  cudaMemcpy(D_,devO1Ptr,3072 * sizeof(double ),cudaMemcpyDeviceToHost);
  cudaFree(devO1Ptr);
  cudaFree(devI1Ptr);
  cudaFree(devI2Ptr);
}

__global__ void Action_noImage_GPU(double *D_,double *maskCenter,double (*SolventMols_)[1024][3])
{
  int bx;
  bx = blockIdx.x;
  int tx;
  tx = threadIdx.x;
  double maxD;
  double Dist;
  int t2;
  int t4;
  D_[bx] = maxD;
//main dist2_noImage code
//double *a1 = maskCenter.Dptr(); //center of solute molecule
//double *a2 = frmIn.XYZ(*solvent_atom);
//double *a1 = maskCenter; //center of solute molecule
//double *a2 = SolventMols_[solventMol][solventAtom];  
//double x = a1[0] - a2[0];
//double y = a1[1] - a2[1];
//double z = a1[2] - a2[2];
//Dist = (x*x + y*y + z*z);
  Dist = maskCenter[0] * SolventMols_[bx][tx][0] + maskCenter[1] * SolventMols_[bx][tx][1] + maskCenter[2] * SolventMols_[bx][tx][2];
  if (Dist + 1 <= D__(bx)) 
    D_[bx] = Dist;
}
