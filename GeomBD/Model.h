#ifndef Model_h_
#define Model_h_

#include "Main.h"
#include "Body.h"
#include "Grid.h"
#include "Grid_Type.h"
#include "Grid_EX.h"
#include "Grid_ES.h"
#include "Grid_D.h"
#include "Session.h"


 
class Model {
  friend class Body;
  friend class Receptor;

  public:
    // {con|de}structors
    Model();
    Model(string inputfn, string outputfn, string logfn);
    ~Model();

  public:
    // input file
    string ifn;
    void parseInputFile();
    void parseReceptorPQR(string rfn);
    void parseLigandPQR(string lfn);

    // trajectory output
    string ofn;
    void writeCoordinatesPQR();

    // log output
    string lfn;
    fstream lout;
    void printRateConstant();

  private:
    // random number generator
    VSLStreamStatePtr *rngCPU;
    vertex *rand;
    void initializeRNG(int seed);
    void generateNormal();

  public:
    // simulation parameters
    int threads;
    int k_trj;
    int k_log;
    double T;
    double viscosity;
    double convergence;
    int max_simulations;
    int fd_order;

  public:
    vector< Session* > sessions;
    vector< Body* > ligands;
    void populateLigands();

  public:
    // grids
    vector< Grid_ES* > esmaps;
    vector< Grid_D* > dmaps;
    vector< Grid_Type* > typemaps;
    vector< Grid_EX* > exmaps;
    Grid_EX *debug_map;

  public:
    // system/receptor definitions
    vertex center;
    double system_r;
    double receptor_radius;
    vertex bounds_min, bounds_max;

  public:
    // time stepping
    int step;
    bool done;
    double dt_fine, dt_coarse;
    double dt_scale_start, dt_scale_end;
    void run();
    void integrate();

};


#endif
