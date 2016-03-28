#include <cmath>
#include <algorithm> // sort
#include <cfloat> // DBL_MAX
#include "Action_Closest.h"
#include "CpptrajStdio.h"

//**************************************************
// issues that need to be addressed before cuda-chill can be used
// 1. .begin() and .end() need to be replaced by dummy variables (i.e see next point)
// 2. need to remove the vector and convert it into a static array ( will fix above issue)
// 3. SolventMols_ is already an array of struct (no issue)
// 4. 



    // struct MolDist {
    //   int mol;        ///< Original solvent molecule number (starts from 1).
    //   double D;       ///< Closest distance of solvent molecule to atoms in distanceMask.
    //   AtomMask mask;  ///< Original topology solvent molecule atom mask.
    //   Iarray solventAtoms; ///< Actual solvent atom #s to loop over.
    // };


    // std::vector<MolDist> SolventMols_;



//code kernel staged for cuda chill


//using dist for no image 
// and kernel for when we use solute molecule center
void Action_Closest::Action_NoImage_Center(Frame&, double maxD)
{
	double Dist;
	int solventMol;

	AtomMask::const_iterator solute_atom;
	Iarray::const_iterator solvent_atom;

	Vec3 maskCenter = frmIn.VGeometricCenter( distanceMask_ );
	for (solventMol=0; solventMol < NsolventMolecules_; solventMol++) {  //standard loop 
		SolventMols_[solventMol].D = maxD;
		for (solvent_atom = SolventMols_[solventMol].solventAtoms.begin();
			solvent_atom != SolventMols_[solventMol].solventAtoms.end(); ++solvent_atom)
		{

			//main dist2_noImage code
			double *a1 = maskCenter.Dptr(); //center of solute molecule
			double *a2 = frmIn.XYZ(*solvent_atom);

			double x = a1[0] - a2[0];
			double y = a1[1] - a2[1];
			double z = a1[2] - a2[2];

			Dist = (x*x + y*y + z*z);

			if (Dist < SolventMols_[solventMol].D) 
				SolventMols_[solventMol].D = Dist;
		}
	}

}



//using dist for no image 
// and kernel for when we find distance with respect to each solute atom
void Action_Closest::Action_NoImage_NoCenter(Frame&, double maxD)
{
	double Dist;
	int solventMol;

	AtomMask::const_iterator solute_atom;
	Iarray::const_iterator solvent_atom;
	
	for (solventMol=0; solventMol < NsolventMolecules_; solventMol++) {
      //if (debug_ > 1)
        //mprintf("DEBUG: Calculating distance for molecule %i\n", solventMol);
      // Set the initial minimum distance for this solvent mol to be the
      // max possible distance.
		SolventMols_[solventMol].D = maxD;
      // Calculate distance between each atom in distanceMask and atoms in solvent Mask
		for (solvent_atom = SolventMols_[solventMol].solventAtoms.begin();
			solvent_atom != SolventMols_[solventMol].solventAtoms.end(); ++solvent_atom)
		{
			for (solute_atom = distanceMask_.begin(); 
				solute_atom != distanceMask_.end(); ++solute_atom)
			{

			//main dist2_noImage code
				double *a1 = frmIn.XYZ(*solute_atom);  // each atom of solute molecule 
				double *a2 = frmIn.XYZ(*solvent_atom);

				double x = a1[0] - a2[0];
				double y = a1[1] - a2[1];
				double z = a1[2] - a2[2];

				Dist = (x*x + y*y + z*z);



				if (Dist < SolventMols_[solventMol].D) 
					SolventMols_[solventMol].D = Dist;

          // if (debug_ > 2)
          //   mprintf("DEBUG: SolvMol %i, soluteAtom %i, solventAtom %i, D= %f, minD= %f\n",
          //           solventMol, *solute_atom, *solvent_atom, Dist,
          //           sqrt(SolventMols_[solventMol].D));
			}
		}
      //f (debug_ > 1) mprintf("DEBUG:\tMol %8i minD= %lf\n",solventMol, SolventMols_[solventMol].D);
    } // END fo

}