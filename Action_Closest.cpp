#include <cmath>
#include <algorithm> // sort
#include <cfloat> // DBL_MAX
#include "Action_Closest.h"
#include "CpptrajStdio.h"
#include <cstdio>
#include <cuda_runtime_api.h>
#include <cuda.h>

// CONSTRUCTOR
Action_Closest::Action_Closest() :
  closestWaters_(0),
  firstAtom_(false),
  useMaskCenter_(false),
  NsolventMolecules_(0),
  debug_(0)
{} 

void Action_Closest::Help() const {
  mprintf("\t<# to keep> <mask> [noimage] [first | oxygen] [center]\n"
//          "\t[closestout <filename> [name <setname>]] [outprefix <parmprefix>]\n"
//          "\t[parmout <file>]\n"
          "  Keep only the closest <# to keep> solvent molecules to atoms in <mask>.\n"
          "  Molecules can be marked as solvent with the 'solvent' command.\n"
          "  If 'center' specified use geometric center of atoms in <mask>.\n");
}

// Action_Closest::Init()
Action_Closest::RetType Action_Closest::Init(ArgList& actionArgs, int debugIn)
{
  debug_ = debugIn;
  // Get Keywords
  closestWaters_ = actionArgs.getNextInteger(-1);
  if (closestWaters_ < 0) {
    mprinterr("Error: Invalid # solvent molecules to keep (%i).\n",
              closestWaters_);
    return Action_Closest::ERR;
  }
  if ( actionArgs.hasKey("oxygen") || actionArgs.hasKey("first") )
    firstAtom_=true;
  useMaskCenter_ = actionArgs.hasKey("center");
  image_.InitImaging( !(actionArgs.hasKey("noimage")) );

  // Get Masks
  std::string mask1 = actionArgs.GetMaskNext();
  if (mask1.empty()) {
    mprinterr("Error: No mask specified.\n");
    return Action_Closest::ERR;
  }
  distanceMask_.SetMaskString(mask1);

  mprintf("    CLOSEST: Finding closest %i solvent molecules to atoms in mask %s\n",
          closestWaters_, distanceMask_.MaskString());
  if (useMaskCenter_)
    mprintf("\tGeometric center of atoms in mask will be used.\n");
  if (!image_.UseImage()) 
    mprintf("\tImaging will be turned off.\n");
  if (firstAtom_)
    mprintf("\tOnly first atom of solvent molecule used for distance calc.\n");
  return Action_Closest::OK;
}

// Action_Closest::Setup()
/** Like the strip action, closest will modify the current parm keeping info
  * for atoms in mask plus the closestWaters solvent molecules. Set up the
  * vector of MolDist objects, one for every solvent molecule in the original
  * parm file. Atom masks for each solvent molecule will be set up.
  */
Action_Closest::RetType Action_Closest::Setup(Topology const& topIn, CoordinateInfo const& cInfoIn)
{
  // If there are no solvent molecules this action is not valid.
  if (topIn.Nsolvent()==0) {
    mprintf("Warning: Parm %s does not contain solvent.\n",topIn.c_str());
    return Action_Closest::SKIP;
  }
  // If # solvent to keep >= solvent in this parm the action is not valid.
  if (closestWaters_ >= topIn.Nsolvent()) {
    mprintf("Warning: # solvent to keep (%i) >= # solvent molecules in '%s' (%i)\n",
            closestWaters_, topIn.c_str(), topIn.Nsolvent());
    return Action_Closest::SKIP;
  }
  image_.SetupImaging( cInfoIn.TrajBox().Type() );
  if (image_.ImagingEnabled())
    mprintf("\tDistances will be imaged.\n");
  else
    mprintf("\tImaging off.\n"); 
  // LOOP OVER MOLECULES
  // 1: Check that all solvent molecules contain same # atoms. Solvent 
  //    molecules must be identical for the command to work properly; 
  //    the prmtop strip occurs only once so the solvent params become fixed.
  // 2: Set up a mask for all solvent molecules.
  SolventMols_.clear();
  // NOTE: May not be necessary to init 'solvent'
  MolDist solvent;
  solvent.D = 0.0;
  solvent.mol = 0;
  SolventMols_.resize(topIn.Nsolvent(), solvent);
  std::vector<MolDist>::iterator mdist = SolventMols_.begin();
  // 3: Set up the soluteMask for all non-solvent molecules.
  int molnum = 1;
  int nclosest = 0;
  int NsolventAtoms = -1;
  for (Topology::mol_iterator Mol = topIn.MolStart();
                              Mol != topIn.MolEnd(); ++Mol)
  {
    if ( Mol->IsSolvent() ) {
      // Solvent, check for same # of atoms.
      if (NsolventAtoms == -1)
        NsolventAtoms = Mol->NumAtoms();
      else if ( NsolventAtoms != Mol->NumAtoms() ) {
        mprinterr("Error: Solvent molecules in '%s' are not of uniform size.\n"
                  "Error:   First solvent mol = %i atoms, solvent mol %i = %i atoms.\n",
                  topIn.c_str(), NsolventAtoms, molnum, (*Mol).NumAtoms());
        return Action_Closest::ERR;
      }
      // mol here is the output molecule number which is why it starts from 1.
      mdist->mol = molnum;
      // Solvent molecule mask
      mdist->mask.AddAtomRange( Mol->BeginAtom(), Mol->EndAtom() );
      // Atoms in the solvent molecule to actually calculate distances to.
      if (firstAtom_) {
        mdist->solventAtoms.assign(1, Mol->BeginAtom() );
      } else {
        mdist->solventAtoms.clear();
        mdist->solventAtoms.reserve( Mol->NumAtoms() );
        for (int svatom = Mol->BeginAtom(); svatom < Mol->EndAtom(); svatom++)
          mdist->solventAtoms.push_back( svatom );
      }
      if (debug_ > 0) {
        mprintf("DEBUG:\tSet up mol %i:", mdist->mol); // DEBUG
        mdist->mask.PrintMaskAtoms("solvent"); // DEBUG
        mprintf("\n"); // DEBUG
      }
      ++mdist;
    }
    ++molnum;
  }

  // Setup distance atom mask
  // NOTE: Should ensure that no solvent atoms are selected!
  if ( topIn.SetupIntegerMask(distanceMask_) ) return Action_Closest::ERR;
  if (distanceMask_.None()) {
    mprintf("Warning: Distance mask '%s' contains no atoms.\n",
            distanceMask_.MaskString());
    return Action_Closest::SKIP;
  }
  distanceMask_.MaskInfo();

  // Check the total number of solvent atoms to be kept.
  NsolventAtoms *= closestWaters_;
  mprintf("\tKeeping %i solvent atoms.\n",NsolventAtoms);
  if (NsolventAtoms < 1) {
    mprintf("Warning: # of solvent atoms to be kept is < 1.\n");
    return Action_Closest::SKIP;
  }
  NsolventMolecules_ = (int)SolventMols_.size();

  return Action_Closest::OK;
}

// Action_Closest::DoAction()
/** Find the minimum distance between atoms in distanceMask and each 
  * solvent Mask.
  */
Action_Closest::RetType Action_Closest::DoAction(int frameNum, Frame& frmIn) { 
  double maxD;
  Matrix_3x3 ucell, recip;
  AtomMask::const_iterator solute_atom;
  Iarray::const_iterator solvent_atom;

  if (image_.ImagingEnabled()) {
    frmIn.BoxCrd().ToRecip(ucell, recip);
    // Calculate max possible imaged distance
    maxD = frmIn.BoxCrd().BoxX() + frmIn.BoxCrd().BoxY() + 
           frmIn.BoxCrd().BoxZ();
    maxD *= maxD;
  } else {
    // If not imaging, set max distance to an arbitrarily large number
    maxD = DBL_MAX;
  }

  //subroutines to find the distance
  // if (image_.ImageType() == NOIMAGE)
  //   Action_NoImage(frmIn,maxD);
  // else if (image_.ImageType() == ORTHO)
  //   Action_ImageOrtho(frmIn,maxD);
  // else
  //   Action_ImageNonOrtho(frmIn,maxD, ucell,recip);

//remove this ..TODOi

useMaskCenter_ = false;
cudaEvent_t start_event, stop_event;
float elapsed_time_seq;

bool v[2] = { true, false };
int type = 1;   //keep it no imaging (as of now)
char* dict[3] = {"NONE", "ORTHO", "NON-ORTHO"};

for(int k =0 ; k < 2 ; k++)
{
  printf("Solute Center : %s\n", v[k] ? "YES" : "NO");
  printf("Imaging : %s\n", dict[type]); 
  useMaskCenter_ = v[k];

cudaEventCreate(&start_event);
cudaEventCreate(&stop_event);
cudaEventRecord(start_event, 0);

//serial section of the code
if (type == 0 )
  Action_NoImage(frmIn,maxD);
else if (type == 1)
  Action_ImageOrtho(frmIn,maxD);
else{
  printf("Error, invalid imaging type\n");
  exit(1);
}



cudaThreadSynchronize();
cudaEventRecord(stop_event, 0);
cudaEventSynchronize(stop_event);
cudaEventElapsedTime(&elapsed_time_seq,start_event, stop_event );
printf("Done with kernel SEQ Kernel Time: %.2f\n", elapsed_time_seq);



bool result = true;
float elapsed_time_gpu;
if (useMaskCenter_)
  result = cuda_action_center(frmIn,maxD,ucell ,recip ,type ,elapsed_time_gpu); 
else
  result = cuda_action_no_center(frmIn,maxD,ucell ,recip ,type ,elapsed_time_gpu);//handling all the data formatting and copying etc
// we will only care about kernel time
//fixing the overhead will be later

  if(result){
    printf("CUDA PASS\n");
  }
  else{
    printf("CUDA FAIL!\n");
    exit(0);
  }

  printf("Seq Time:  = %0.2f\n", elapsed_time_seq);
  printf("CUDA Time: = %0.2f\n", elapsed_time_gpu);
  printf("Speedup =  %0.2f\n", elapsed_time_seq/elapsed_time_gpu);
}

  // Sort distances
  std::sort( SolventMols_.begin(), SolventMols_.end(), moldist_cmp() );
  // Add first closestWaters solvent atoms to stripMask
  std::vector<MolDist>::iterator solventend = SolventMols_.begin() + closestWaters_;
  for ( std::vector<MolDist>::const_iterator solvent = SolventMols_.begin();
                                             solvent != solventend;
                                           ++solvent ) 
  {
    solvent_atom = solvent->mask.begin();
    mprintf("\tMol= %8i  Atom= %8i  Dist= %10.4f\n", solvent->mol,
            *solvent_atom + 1, sqrt( solvent->D ));
  }

  return Action_Closest::OK;
}




void Action_Closest::Action_ImageNonOrtho(Frame& frmIn, double maxD, Matrix_3x3 ucell, Matrix_3x3 recip)
{
  double Dist;
  int solventMol;
  AtomMask::const_iterator solute_atom;
  Iarray::const_iterator solvent_atom;
    // Loop over all solvent molecules in original frame
  if (useMaskCenter_) {
    Vec3 maskCenter = frmIn.VGeometricCenter( distanceMask_ ); //can be calculated outside
    for (solventMol=0; solventMol < NsolventMolecules_; solventMol++) {
      SolventMols_[solventMol].D = maxD;
      for (solvent_atom = SolventMols_[solventMol].solventAtoms.begin();
           solvent_atom != SolventMols_[solventMol].solventAtoms.end(); ++solvent_atom)
      {

        Dist = DIST2_ImageNonOrtho( maskCenter.Dptr(),
                      frmIn.XYZ(*solvent_atom),ucell, recip);  //frame translation can be done inside gpu
        if (Dist < SolventMols_[solventMol].D) 
          SolventMols_[solventMol].D = Dist;
      }
    }
  } else {
    for (solventMol=0; solventMol < NsolventMolecules_; solventMol++) {
      if (debug_ > 1)
        mprintf("DEBUG: Calculating distance for molecule %i\n", solventMol);
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
          Dist = DIST2_ImageNonOrtho(frmIn.XYZ(*solute_atom),
                       frmIn.XYZ(*solvent_atom), ucell, recip);
          if (Dist < SolventMols_[solventMol].D) 
            SolventMols_[solventMol].D = Dist;
          if (debug_ > 2)
            mprintf("DEBUG: SolvMol %i, soluteAtom %i, solventAtom %i, D= %f, minD= %f\n",
                    solventMol, *solute_atom, *solvent_atom, Dist,
                    sqrt(SolventMols_[solventMol].D));
        }
      }
      if (debug_ > 1) mprintf("DEBUG:\tMol %8i minD= %lf\n",solventMol, SolventMols_[solventMol].D);
    } // END for loop over solventMol
  }

}

void Action_Closest::Action_ImageOrtho(Frame& frmIn, double maxD)
{
  double Dist;
  int solventMol;
  AtomMask::const_iterator solute_atom;
  Iarray::const_iterator solvent_atom;
    // Loop over all solvent molecules in original frame
  if (useMaskCenter_) {
    Vec3 maskCenter = frmIn.VGeometricCenter( distanceMask_ );
    for (solventMol=0; solventMol < NsolventMolecules_; solventMol++) {
      SolventMols_[solventMol].D = maxD;
      for (solvent_atom = SolventMols_[solventMol].solventAtoms.begin();
           solvent_atom != SolventMols_[solventMol].solventAtoms.end(); ++solvent_atom)
      {

        Dist = DIST2_ImageOrtho( maskCenter.Dptr(),
                      frmIn.XYZ(*solvent_atom),frmIn.BoxCrd());
        if (Dist < SolventMols_[solventMol].D) 
          SolventMols_[solventMol].D = Dist;
      }
    }
  } else {
    for (solventMol=0; solventMol < NsolventMolecules_; solventMol++) {
      if (debug_ > 1)
        mprintf("DEBUG: Calculating distance for molecule %i\n", solventMol);
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
          Dist = DIST2_ImageOrtho(frmIn.XYZ(*solute_atom),
                       frmIn.XYZ(*solvent_atom), frmIn.BoxCrd());
          if (Dist < SolventMols_[solventMol].D) 
            SolventMols_[solventMol].D = Dist;
          if (debug_ > 2)
            mprintf("DEBUG: SolvMol %i, soluteAtom %i, solventAtom %i, D= %f, minD= %f\n",
                    solventMol, *solute_atom, *solvent_atom, Dist,
                    sqrt(SolventMols_[solventMol].D));
        }
      }
      if (debug_ > 1) mprintf("DEBUG:\tMol %8i minD= %lf\n",solventMol, SolventMols_[solventMol].D);
    } // END for loop over solventMol
  }

}


//pulling out the dist control statement
void Action_Closest::Action_NoImage(Frame& frmIn,double maxD)
{
  double Dist;
  int solventMol;
  AtomMask::const_iterator solute_atom;
  Iarray::const_iterator solvent_atom;
    // Loop over all solvent molecules in original frame
  if (useMaskCenter_) {
    Vec3 maskCenter = frmIn.VGeometricCenter( distanceMask_ );
    for (solventMol=0; solventMol < NsolventMolecules_; solventMol++) {
      SolventMols_[solventMol].D = maxD;
      for (solvent_atom = SolventMols_[solventMol].solventAtoms.begin();
           solvent_atom != SolventMols_[solventMol].solventAtoms.end(); ++solvent_atom)
      {

        Dist = DIST2_NoImage( maskCenter.Dptr(),
                      frmIn.XYZ(*solvent_atom));
        //printf("DIST  = %f\n", Dist);
	if (Dist < SolventMols_[solventMol].D) 
          SolventMols_[solventMol].D = Dist;
      }
    }
  } else {
    for (solventMol=0; solventMol < NsolventMolecules_; solventMol++) {
      if (debug_ > 1)
        mprintf("DEBUG: Calculating distance for molecule %i\n", solventMol);
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
          Dist = DIST2_NoImage(frmIn.XYZ(*solute_atom),
                       frmIn.XYZ(*solvent_atom));
          //printf("no center DIST  = %f\n", Dist);
	if (Dist < SolventMols_[solventMol].D) 
            SolventMols_[solventMol].D = Dist;
          if (debug_ > 2)
            mprintf("DEBUG: SolvMol %i, soluteAtom %i, solventAtom %i, D= %f, minD= %f\n",
                    solventMol, *solute_atom, *solvent_atom, Dist,
                    sqrt(SolventMols_[solventMol].D));
        }
      }
      if (debug_ > 1) mprintf("DEBUG:\tMol %8i minD= %lf\n",solventMol, SolventMols_[solventMol].D);
    } // END for loop over solventMol
  }

}


