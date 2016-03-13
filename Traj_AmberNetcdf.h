#ifndef INC_TRAJ_AMBERNETCDF_H
#define INC_TRAJ_AMBERNETCDF_H
#ifdef BINTRAJ
#include "NetcdfFile.h"
#include "CpptrajFile.h"
#include "Topology.h"
#include "ArgList.h"
/// Reads and writes Amber Netcdf format trajectories. 
class Traj_AmberNetcdf : private NetcdfFile {
  public:
    static const int TRAJIN_ERR;
    Traj_AmberNetcdf();
    ~Traj_AmberNetcdf();
    static void ReadHelp();
    static void WriteHelp();
    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(FileName const&, Topology*);
    int setupTrajout(FileName const&, Topology*, CoordinateInfo const&,int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int readVelocity(int, Frame&);
    int readForce(int, Frame&);
    int writeFrame(int,Frame const&);
    void Info();
    int processWriteArgs(ArgList&);
    int processReadArgs(ArgList&);
    FileName const& Filename() const { return filename_; }
    CoordinateInfo const& CoordInfo() const { return cInfo_; }
    std::string const& Title() { return title_; } 
  private:
    CoordinateInfo cInfo_;
    CoordinateInfo& SetCoordInfo(CoordinateInfo const& c) { cInfo_ = c; }
   
    std::string title_;
    void SetTitle( std::string const& t ) { title_ = t; }

    int debug_;
    float *Coord_;
    FileName filename_;
    int eptotVID_;
    int binsVID_;
    bool useVelAsCoords_;
    bool readAccess_;
    bool outputTemp_;
    bool outputVel_;
    bool outputFrc_;
};
#endif
#endif
