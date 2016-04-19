//cpp file to hold the setup routine 
//

#include "Action_Closest.h"
#include <cfloat> // DBL_MAX



//kernel_wrapper defs
void Action_NoImage_Center(double *SolventMols_,double *D_, double maskCenter[3],double maxD,int  NMols,int NAtoms);

///////////////////////////

//rsposibility of this function
// - pull out corrdinate by frame 
// - rewrite data 
// - gpu is only responsible for finding distance and doing the MIN reduction operator
void Action_Closest::cuda_action(Frame& frmIn, double maxD, Matrix_3x3 ucell, Matrix_3x3 recip,int type, bool imaginEnabed)
{
	Vec3 maskCenter_holder =  frmIn.VGeometricCenter( distanceMask_ );
	double* maskCenter = maskCenter_holder.Dptr();

	//allocate space and rewrite 
	int NMols = SolventMols_.size();
	int NAtoms = SolventMols_[0].solventAtoms.size();  //guaranteed to same size  -  due to setup  
	int Nlinear_Solvent = 3 * NMols *  NAtoms
	double *linear_Solvent = new double[Nlinear_Solvent]
	double *D_ = new double[NMols]

	//write in 
	for(int sMol =0; sMol  < NMols; sMol++){
		for(int sAtom = 0, sAtom < NAtoms; sAtom++){
			index =  (sAtom * 3 ) + (sMol * 3 * NAtoms);
			double *a = frmIn.XYZ(SolventMols[sMol].solventAtoms[sAtom]);

			linear_Solvent[index + 0] = a[0];
			linear_Solvent[index + 1] = a[1];
			linear_Solvent[index + 2] = a[2];

		}
	}


	//call the correct function 
	//TODO
	//need to handle cases as well 
	//TODO

	Action_NoImage_Center(linear_Solvent, D_, maskCenter, maxD, NMols, NAtoms);


	//



	//copying back the D__ into the right place
	for(int sMol = 0; sMol < NMols; sMol++){
		SolventMols[sMol].D = D_[sMol];
	}


	//done



}