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
    string ifn;
    void parseInputFile();
    void parseReceptorPQR(string rfn);
    void parseLigandPQR(string lfn);

  private:
    string ofn;
    void writeCoordinatesPQR();

  private:
    void printRateConstant();

  public:
    string lfn;
    fstream lout;

  private:
    VSLStreamStatePtr *rngCPU;
    vertex *rand;

  public:
    void initializeRNG(int seed);
    void generateNormal();

  public:
    double T;
    double viscosity;
    double receptorRhyd;

    vector< Session* > sessions;

  public:
    vector< Body* > ligands;
    void populateLigands();

  public:
    vector< Grid_ES* > esmaps;
    vector< Grid_D* > dmaps;
    vector< Grid_Type* > typemaps;
    vector< Grid_EX* > exmaps;

    Grid_EX *debug_map;

    vertex center;
    double system_r;
    vertex bounds_min, bounds_max;

  public:
    int Nthreads;
    int Vtraj;
    int Vprint;
    int active;
    int step;
    double dt_fine, dt_coarse;
    double dt_scale_start, dt_scale_end;
    void integrate();

  public:
    double convergence;
    int max_simulations;

  public:
    Model();
    Model(string inputfn, string outputfn, string logfn);
    ~Model();

    bool done;
    void run();

};


#endif
