

//this is only  used for cuda-chill
//heavy simplification

extern int  N;
extern int NsolventMolecules_
extern int NsolventAtoms_


struct MolDist {
      int mol;        ///< Original solvent molecule number (starts from 1).
      double D;       ///< Closest distance of solvent molecule to atoms in distanceMask.
      AtomMask mask;  ///< Original topology solvent molecule atom mask.
      double solventAtoms[NsolventAtoms_][3]; ///< Actual solvent atom #s to loop over.
  };

//using dist for no image 
// and kernel for when we use solute molecule center
  void Action_NoImage_Center(MolDist SolventMols_[NsolventMolecules_],double maskCenter[3] ,double maxD)
  {
  	double Dist;
  	int solventMol, solvent_atom;

  	//Vec3 maskCenter = frmIn.VGeometricCenter( distanceMask_ );
	for (solventMol=0; solventMol < NsolventMolecules_; solventMol++) {  //standard loop 
		SolventMols_[solventMol].D = maxD;
		for (solvent_atom = 0; solvent_atom < NsolventAtoms_; solvent_atom++)
		{
			//main dist2_noImage code
			//double *a1 = maskCenter.Dptr(); //center of solute molecule
			//double *a2 = frmIn.XYZ(*solvent_atom);

			double *a1 = maskCenter; //center of solute molecule
			double *a2 = SolventMols_[solventMol].solventAtoms[solvent_atom]  

			double x = a1[0] - a2[0];
			double y = a1[1] - a2[1];
			double z = a1[2] - a2[2];

			Dist = (x*x + y*y + z*z);

			if (Dist < SolventMols_[solventMol].D) 
				SolventMols_[solventMol].D = Dist;
		}
	}

}