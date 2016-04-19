




__global__ void Action_noImage_GPU(double *D_,double *maskCenter,double *SolventMols_,double maxD, int Nmols , int NAtoms)
{
  int bx;
  bx = blockIdx.x;
  int tx;
  tx = threadIdx.x;
  double Dist;
  int t2;
  int t4;
  D_[bx] = maxD;
  int sIndex =  bx*NAtoms*3 + tx*3;

  Dist = maskCenter[0] * SolventMols_[index + 0] + maskCenter[1] * SolventMols_[index + 1] + maskCenter[2] * SolventMols_[index + 2];
  if (Dist  < D_[bx]) 
    D_[bx] = Dist;
}