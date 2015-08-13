#ifndef Model_h_
#define Model_h_

#include "Main.h"
#include "Body.h"
#include "BinaryPotentialMap.h"
#include "TypePotentialMap.h"
#include "ESPotentialMap.h"

enum LigandPosition {
  LIGAND_POSITION_RADIAL,
  LIGAND_POSITION_ABSOLUTE,
  LIGAND_POSITION_RANDOM
};

struct BindingSite {
  double r2;
  double x, y, z;
};

 
class Model {
  friend class Body;
  friend class Receptor;

  public:
    string fldfn, dlgfn;
    void parseFLD();
    void parseDLG();

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
    double t_limit;
    vertex bounds;
    double r2_escape;
    double receptorRoG;

    LigandPosition ligandPosition;
    double r_ligand;
    vertex R_ligand;

  public:
    vector< ESPotentialMap* > esmaps;
    vector< TypePotentialMap* > typemaps;
    vector< Body* > conformations;
    vector< Body* > ligands;
    //vector< BindingSite > bindingSites;
    BindingSite bindingSite;

  public:
    int Nthreads;
    int Nreplicates;
    int active;
    int step;
    double dt_fine, dt_coarse;
    void integrate();

  public:
    Model();
    Model(string fldfn, string dlgfn, string outputfn, double rRoG);
    ~Model();

  public:
    bool done;
    void run();

    void printRateConstant();
};


#endif
