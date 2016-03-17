/*! \file main.cpp
    \brief A reduced version of the closest command from Cpptraj.
    \author Daniel R. Roe
    \date 2016-01-22
 */
#include <cstdio>
#include "Action_Closest.h"
#include "Parm_Amber.h"
#include "Traj_AmberNetcdf.h"

static void Help() {
  printf("Options: [-p <Amber topology file>] [-y <Amber NetCDF trajectory file>]\n"
         "         [-n <# closest waters>] [-m <solute mask expression>]\n");
}

int main(int argc, char** argv) {
  fprintf(stderr,"\n     CLOSEST COMMAND from CPPTRAJ\n");
  // Default options
  int debug = 0;
  std::string topName("tz2.truncoct.parm7");
  std::string trajName("tz2.truncoct.nc");
  std::string Nclosest("10");
  std::string maskExp(":1-13");

  // Get options
  for (int i = 1; i < argc; i++) {
    std::string Arg(argv[i]);
    if (Arg == "-p" && i+1 != argc)
      topName.assign( argv[++i] );
    else if (Arg == "-y" && i+1 != argc)
      trajName.assign( argv[++i] );
    else if (Arg == "-n" && i+1 != argc)
      Nclosest.assign( argv[++i] );
    else if (Arg == "-m" && i+1 != argc)
      maskExp.assign( argv[++i] );
    else if (Arg == "-h" || Arg == "--help") {
      Help();
      return 0;
    } else {
      fprintf(stderr,"Error: Unrecognized option: %s\n", argv[i]);
      Help();
      return 1;
    }
  }
  printf("  Topology name: %s\n", topName.c_str());
  printf("  Trajectory name: %s\n", trajName.c_str());
  printf("  # closest: %s\n", Nclosest.c_str());
  printf("  Solute mask: '%s'\n", maskExp.c_str());

  // Load topology from file
  Parm_Amber parmFile;
  CpptrajFile fileIn;
  fileIn.SetupRead( topName, debug );
  if (!parmFile.ID_ParmFormat(fileIn)) {
    fprintf(stderr,"Error: Topology file '%s' is not Amber topology.\n", topName.c_str());
    return 1;
  }
  Topology top;
  if (parmFile.ReadParm(fileIn.Filename(), top)) return 1;
  top.CommonSetup();
  top.Summary();

  // Set up trajectory file for reading
  Traj_AmberNetcdf trajFile;
  fileIn.SetupRead( trajName, debug );
  if (!trajFile.ID_TrajFormat( fileIn )) {
    fprintf(stderr,"Error: Trajectory file '%s' is not Amber NetCDF.\n", trajName.c_str());
    return 1;
  }
  int Nframes = trajFile.setupTrajin( fileIn.Filename(), &top );
  if (Nframes == Traj_AmberNetcdf::TRAJIN_ERR) return 1;
  printf("\tTrajectory: '%s'", trajFile.Filename().full());
  trajFile.Info();
  printf(", %i frames.\n", Nframes);

  // Allocate frame for reading in coordinates
  Frame frm;
  frm.SetupFrameV( top.Atoms(), trajFile.CoordInfo() );
  frm.Info("Input frame");

  // Initialize and set up action. Keep closest 10 to residues 1-13, calc
  // distances between all solvent and solute atoms.
  Action_Closest closest;
  ArgList actionArgs(Nclosest + " " + maskExp);
  if (closest.Init(actionArgs, debug) != Action_Closest::OK) return 1;
  if (closest.Setup(top, trajFile.CoordInfo()) != Action_Closest::OK) return 1;

  // Loop over all frames in the trajectory.
  if (trajFile.openTrajin()) return 1;
  for (int set = 0; set != Nframes; set++) {
    printf("Frame %i\n", set+1);
    if (trajFile.readFrame(set, frm)) break;
    closest.DoAction(set, frm);  //new frame each time  -  maybe preloading then for gpu
  }
  trajFile.closeTraj();
  fprintf(stderr,"\nComplete.\n\n");
  return 0;
}
