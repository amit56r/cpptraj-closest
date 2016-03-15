#ifndef INC_ACTION_CLOSEST_H
#define INC_ACTION_CLOSEST_H
#include "ImagedAction.h"
#include "Topology.h"
#include "ArgList.h"
// NOTE: Simplified version of Action_Closest for testing with GPU code.
/// Modify the state so that only the closest solvent molecules are kept.
class Action_Closest {
  public:
    Action_Closest();
    void Help() const;
    enum RetType { OK, ERR, USE_ORIGINAL_FRAME, SUPPRESS_COORD_OUTPUT,
                   SKIP, MODIFY_TOPOLOGY, MODIFY_COORDS };
    RetType Init(ArgList&, int);
    RetType Setup(Topology const&, CoordinateInfo const&);
    RetType DoAction(int, Frame&);
  private:

    ImagedAction image_;    ///< Imaging routines.
    int closestWaters_;     ///< Closest # of molecules to keep.
    bool firstAtom_;        ///< If true just calc based on molecule first atom.
    bool useMaskCenter_;    ///< If true use geometric center of mask.
    AtomMask distanceMask_; ///< Mask of atoms to calculate distance from solvent to.
    int NsolventMolecules_; ///< # of solvent molecules in SolventMols.
    int debug_;
    typedef std::vector<int> Iarray;
    /** The moldist structure is used in order to preserve the original
      * solvent molecule numbers after sorting. */
    struct MolDist {
      int mol;        ///< Original solvent molecule number (starts from 1).
      double D;       ///< Closest distance of solvent molecule to atoms in distanceMask.
      AtomMask mask;  ///< Original topology solvent molecule atom mask.
      Iarray solventAtoms; ///< Actual solvent atom #s to loop over.
    };
    /// Return true if the first molecule is closer than the second
    struct moldist_cmp {
      inline bool operator()(MolDist const& first, MolDist const& second) const {
        return (first.D < second.D);
      }
    };
    std::vector<MolDist> SolventMols_;
};
#endif  
