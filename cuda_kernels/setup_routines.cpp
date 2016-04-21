//cpp file to hold the setup routine 
//

#include "Action_Closest.h"
#include <cfloat> // DBL_MAX
#include <cstdio>
#include <cmath>
#define TOL 1e-7

//kernel_wrapper defs
void Action_NoImage_Center(double *SolventMols_,double *D_, double maskCenter[3],double maxD,int  NMols,int NAtoms, float &time_gpu);

///////////////////////////

//rsposibility of this function
// - pull out corrdinate by frame 
// - rewrite data 
// - gpu is only responsible for finding distance and doing the MIN reduction operator
bool Action_Closest::cuda_action_center(Frame& frmIn, double maxD, Matrix_3x3 ucell, Matrix_3x3 recip,int type, float &time_gpu)
{
	Vec3 maskCenter_holder =  frmIn.VGeometricCenter( distanceMask_ );
	double* maskCenter = maskCenter_holder.Dptr();

	//allocate space and rewrite 
	int NMols = SolventMols_.size();
	int NAtoms = SolventMols_[0].solventAtoms.size();  //guaranteed to same size  -  due to setup  
	int Nlinear_Solvent = 3 * NMols *  NAtoms;
	double *linear_Solvent = new double[Nlinear_Solvent];
	double *D_ = new double[NMols];

	//write in 
	for(int sMol =0; sMol  < NMols; sMol++){
		for(int sAtom = 0; sAtom < NAtoms; sAtom++){
			int index =  (sAtom * 3 ) + (sMol * 3 * NAtoms);
			const double *a = frmIn.XYZ(SolventMols_[sMol].solventAtoms[sAtom]);

			linear_Solvent[index + 0] = a[0];
			linear_Solvent[index + 1] = a[1];
			linear_Solvent[index + 2] = a[2];

		}
	}


	//call the correct function 
	//TODO
	//need to handle cases as well 
	//TODO

	Action_NoImage_Center(linear_Solvent, D_, maskCenter, maxD, NMols, NAtoms,time_gpu);


	//



	//copying back the D__ into the right place
	bool flag = true;
	for(int sMol = 0; sMol < NMols; sMol++){
		if(fabs(SolventMols_[sMol].D - D_[sMol]) > TOL)
		{
			flag = false;
			printf("ACTUAL = %f ; GPU = %f; Where = %d/%d \n",SolventMols_[sMol].D,D_[sMol], sMol+1, NMols);
			break;
		}
		//printf("lhs = %f ; rhs = %f ",SolventMols_[sMol].D,D_[sMol]);
	}

	delete linear_Solvent;
	delete D_;

	return flag;
	//done



}


//TODO
//need to fix this 
bool Action_Closest::cuda_action_no_center(Frame& frmIn, double maxD, Matrix_3x3 ucell, Matrix_3x3 recip,int type, float &time_gpu)
{
	Vec3 maskCenter_holder =  frmIn.VGeometricCenter( distanceMask_ );
	double* maskCenter = maskCenter_holder.Dptr();

	//allocate space and rewrite 
	int NMols = SolventMols_.size();
	int NAtoms = SolventMols_[0].solventAtoms.size();  //guaranteed to same size  -  due to setup  
	int Nlinear_Solvent = 3 * NMols *  NAtoms;
	double *linear_Solvent = new double[Nlinear_Solvent];
	double *D_ = new double[NMols];

	//write in 
	for(int sMol =0; sMol  < NMols; sMol++){
		for(int sAtom = 0; sAtom < NAtoms; sAtom++){
			int index =  (sAtom * 3 ) + (sMol * 3 * NAtoms);
			const double *a = frmIn.XYZ(SolventMols_[sMol].solventAtoms[sAtom]);

			linear_Solvent[index + 0] = a[0];
			linear_Solvent[index + 1] = a[1];
			linear_Solvent[index + 2] = a[2];

		}
	}


	//call the correct function 
	//TODO
	//need to handle cases as well 
	//TODO

	Action_NoImage_Center(linear_Solvent, D_, maskCenter, maxD, NMols, NAtoms,time_gpu);


	//



	//copying back the D__ into the right place
	bool flag = true;
	for(int sMol = 0; sMol < NMols; sMol++){
		if(fabs(SolventMols_[sMol].D - D_[sMol]) > TOL)
		{
			flag = false;
			printf("ACTUAL = %f ; GPU = %f; Where = %d/%d \n",SolventMols_[sMol].D,D_[sMol], sMol+1, NMols);
			break;
		}
		//printf("lhs = %f ; rhs = %f ",SolventMols_[sMol].D,D_[sMol]);
	}

	delete linear_Solvent;
	delete D_;

	return flag;
	//done



}