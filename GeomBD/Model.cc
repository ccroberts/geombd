
#include "Model.h"
#include "Timer.h"
#include "BindingCriteria.h"



Model::Model() {
  viscosity = 1.0/*cP*/ * 0.001/*cP->Pa.s->J.s/m^3*/ * 2.39e-4 /*J->kcal*/ * 1e12 /*s->ps*/ * 1e-30 /*1/m^3->1/A^3*/ * Na /*kcal->kcal/mol*/;
  T = 298.;
  receptor_radius = 1e20;
  k_trj = 10000;
  k_log = 5000;
  convergence = 0.;
  max_simulations = -1;
  fd_order = 2;

  memset(&center, 0, sizeof(vertex));
  memset(&bounds_min, 0, sizeof(vertex));
  memset(&bounds_max, 0, sizeof(vertex));
  step = 0;
  done = false;
  dt_fine = 0.010;
  dt_coarse = 10.000;
  dt_scale_start = 100.;
  dt_scale_end = 500.;

  //debug
  debug_map = NULL;
}



Model::Model(string inputfn, string outputfn, string logfn) : Model() {
  ifn = inputfn;
  if(!file_exists(ifn)) {
    cout << "! Input file \"" << ifn << "\" does not exist. Exiting." << endl;
    exit(EXIT_FAILURE);
  }

  ofn = outputfn;
  lfn = logfn;
  lout.open(lfn, ios_base::out);

  parseInputFile();
  threads = __cilkrts_get_nworkers();

  populateLigands();
  initializeRNG(time(NULL));
}



Model::~Model() {
  for(int k=0; k < threads; k++) 
    vslDeleteStream(&rngCPU[k]);
}



void Model::initializeRNG(int seed) {
  lout << "* Initializing random number generator (N=" << threads << ")" << endl;

  srand(seed);

  rngCPU = (VSLStreamStatePtr*)calloc(threads, sizeof(VSLStreamStatePtr));
  for(int k=0; k < threads; k++) {
    vslNewStream(&rngCPU[k], VSL_BRNG_MT19937, seed+k);
  }

  rand = (vertex*)calloc(2 * ligands.size(), sizeof(vertex));
}



void Model::generateNormal() {
  int blockSize = ((2*ligands.size()*3) / threads);
  cilk_for(int k=0; k < threads; k++) {
    vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, rngCPU[k], blockSize, (double*)rand+(k*blockSize), 0., 1.);
  }
}



void Model::run() {
  Timer t;

  t.start();
  if(step == 0) writeCoordinatesPQR();

  while(!done) {
    integrate();
    step++;

    // Write trajectory and statistics every 'k_trj' steps
    if(step % k_trj == 0) {
      writeCoordinatesPQR();
    }
    if(step % k_log == 0) {
      t.stop();
      lout << "* Step " << step << " (" << (t.duration/k_trj) << " s/step)" << endl;
      printRateConstant();
      t.start();
    }
  }

  // Final information
  printRateConstant();
}


void Model::printRateConstant() {
  for(int i=0; i < sessions.size(); i++) {
    sessions[i]->printRateConstant();
  }
}


void Model::populateLigands() {
  for(int i=0; i < sessions.size(); i++) {
    sessions[i]->populateLigands();
  }
}
