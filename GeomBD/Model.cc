
#include "Model.h"
#include "Timer.h"
#include "BindingCriteria.h"



Model::Model() {
  writeBinders = false;

  rngCPU = NULL;
  rand = NULL;

  T = 298.;
  viscosity = 1.0/*cP*/ * 0.001/*cP->Pa.s->J.s/m^3*/ * 2.39e-4 /*J->kcal*/ * 1e12 /*s->ps*/ * 1e-30 /*1/m^3->1/A^3*/ * Na /*kcal->kcal/mol*/;
  threads = 1;
  fd_order = 2;
  rate_trj = 10000;
  rate_beta = 5000;
  convergence = 1e-4;
  convergence_window = 100;
  max_simulations = -1;

  debug_map = NULL;

  center.x = center.y = center.z = 0.;
  system_extent = 0.;
  receptor_radius = 1e20;
  bounds_min.x = bounds_min.y = bounds_min.z = 0.;
  bounds_max.x = bounds_max.y = bounds_max.z = 0.;

  step = 0;
  done = false;
  dt_fine = 0.010;
  dt_coarse = 1.000;
  dt_scale_start = 100.;
  dt_scale_end = 500.;

  Natoms = 0;
  Nframes = 0;
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

  openTrajectoryDCD();

  t.start();

  while(!done) {
    integrate();
    step++;

    if(step % rate_trj == 0) {
      writeCoordinatesDCD();
    }
    if(step % rate_beta == 0) {
      t.stop();
      lout << "* Step " << step << " (" << (t.duration/rate_trj) << " s/step)" << endl;
      printRateConstant();
      t.start();
    }
  }

  // Final information
  printRateConstant();
  closeTrajectoryDCD();
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

void Model::checkConvergence() {
  for(int i=0; i < sessions.size(); i++) {
    sessions[i]->checkConvergence();
  }
  // check to see if all sessions are done
  bool _done = true;
  for(int i=0; i < sessions.size(); i++) {
    if(! sessions[i]->done) { _done = false; break; }
  }
  done = _done;
}


