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
    // constructors
    Model();
    Model(string inputfn, string outputfn, string logfn);
    // destructor
    ~Model();

  public:
    // input file
    string ifn;
    void parseInputFile();
    void parseReceptorPQR(string rfn);
    void parseLigandPQR(string lfn);

    // trajectory output
    string ofn;   //trajectory filename
    fstream ofd;  //trajectory file handle
    int Natoms;   //neccesary storage
    int Nframes;  //neccesary storage
    void writeCoordinatesPQR();
    void openTrajectoryDCD();
    void writeCoordinatesDCD();
    void closeTrajectoryDCD();

    // bound conformation output
    bool writeBinders;
    string bfn;   //bound pqr filename

    // log output
    string lfn;   //log filename
    fstream lout; //log file stream for writing
    void printRateConstant();

  private:
    // random number generator
    VSLStreamStatePtr *rngCPU;
    vertex *rand;
    void initializeRNG(int seed);
    void generateNormal();

  public:
    // simulation parameters
    double T;                 //temperature of system
    double viscosity;         //viscosity of solvent
    int threads;              //number of threads to parallelize operations
    int fd_order;             //order of finite difference approximation (default = 2)
    int rate_trj;             //frquency of trajectory writes (steps)
    int rate_beta;            //frequency of beta calculations, and thuslog writes (units are number of steps)
    double convergence;       //terminate run once a beta/k value has converged by a certain order of magnitude (default = 1e-4)
    int convergence_window;   //number of recorded beta values required to check for convergence (default = 100)
    int max_simulations;      //terminate run after a certain number of total simulations have been completed (default = none)
    void checkConvergence();

  public:
    vector< Session* > sessions;    //BD sessions
    vector< Body* > ligands;        //all ligands across all sessions
    void populateLigands();

  public:
    // grids
    vector< Grid_ES* > esmaps;      //electrostatic maps
    vector< Grid_D* > dmaps;        //desolvation maps
    vector< Grid_Type* > typemaps;  //atom typed maps (currently just LJ has been implemented)
    vector< Grid_EX* > exmaps;      //exclusion maps
    Grid_EX *debug_map;             //debug maps for recording force magnitudes

  public:
    // system/receptor definitions
    vertex center;                  //center of receptor
    double system_extent;           //receptor maximum extension from center
    double receptor_radius;         //receptor radius of gyration
    vertex bounds_min, bounds_max;  //bounding coordinates

  public:
    // time stepping
    int step;                            //current step of the simulation progression
    bool done;                           //boolean to kill simulation when a termination condition is met
    double dt_fine, dt_coarse;           //timestep scaling values(min->max = fine->coarse)
    double dt_scale_start, dt_scale_end; //timestep scaning distance
    void run();                          //main loop
    void integrate();                    //time step function

};


#endif
