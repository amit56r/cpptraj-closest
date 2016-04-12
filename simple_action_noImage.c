
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
  void Action_NoImage_Center(double SolventMols_[NsolventMolecules_][NsolventAtoms_][3],
  				 double D_[NsolventMolecules_], double maskCenter[3] ,double maxD)
  {
  	double Dist;
  	int solventMol, solventAtom;

  	//Vec3 maskCenter = frmIn.VGeometricCenter( distanceMask_ );
	for (solventMol=0; solventMol < NsolventMolecules_; solventMol++) {  //standard loop 
		D_[solventMol] = maxD;
		for (solventAtom = 0; solventAtom < NsolventAtoms_; solventAtom++)
		{
			//main dist2_noImage code
			//double *a1 = maskCenter.Dptr(); //center of solute molecule
			//double *a2 = frmIn.XYZ(*solvent_atom);

			double *a1 = maskCenter; //center of solute molecule
			//double *a2 = SolventMols_[solventMol][solventAtom];  

			//double x = a1[0] - a2[0];
			//double y = a1[1] - a2[1];
			//double z = a1[2] - a2[2];

			//Dist = (x*x + y*y + z*z);
			Dist  = (SolventMols_[solventMol][solventAtom][0] * SolventMols_[solventMol][solventAtom][0]) +
				(SolventMols_[solventMol][solventAtom][1] * SolventMols_[solventMol][solventAtom][1]) +
   				(SolventMols_[solventMol][solventAtom][2] *SolventMols_[solventMol][solventAtom][2]);

			D_[solventMol] = Dist < D_[solventMol] ?  Dist : D_[solventMol];
			//if (Dist < D_[solventMol]) 
			//	D_[solventMol] = Dist;

		}
	}

}



