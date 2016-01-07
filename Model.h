#ifndef Model_h_
#define Model_h_

#include "Main.h"
#include "Body.h"
#include "BinaryPotentialMap.h"
#include "TypePotentialMap.h"
#include "ESPotentialMap.h"
#include "APBSPotentialMap.h"
#include "Session.h"
#include "ExclusionMap.h"


 
class Model {
  friend class Body;
  friend class Receptor;

  public:
    string ifn;
    void parseInputFile();
    void parseReceptorPDBQT(string rfn);
    void parseLigandPDBQT(string lfn);

  private:
    string ofn;
    void writeCoordinatesPQR();

  private:
    VSLStreamStatePtr *rngCPU;
    vertex *rand;

  public:
    void initializeRNG(int seed);
    void generateNormal();

  public:
    double T;
    double viscosity;
    double receptorRoG;

    vector< Session* > sessions;

  public:
    vector< Body* > ligands;
    void populateLigands();

  public:
    vector< ESPotentialMap* > esmaps;
    vector< APBSPotentialMap* > apbsmaps;
    vector< TypePotentialMap* > typemaps;
    vector< ExclusionMap* > xmaps;

    vertex center;

  public:
    int Nthreads;
    int Vtraj;
    int active;
    int step;
    double dt_fine, dt_coarse;
    void integrate();

  public:
    Model();
    Model(string inputfn, string outputfn);
    ~Model();

  public:
    bool done;
    void run();

    void printRateConstant();
};


#endif
