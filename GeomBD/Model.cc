
#include "Model.h"
#include "Timer.h"
#include "BindingCriteria.h"



Model::Model() {
  viscosity = 1.0/*cP*/ * 0.001/*cP->Pa.s->J.s/m^3*/ * 2.39e-4 /*J->kcal*/ * 1e12 /*s->ps*/ * 1e-30 /*1/m^3->1/A^3*/ * Na /*kcal->kcal/mol*/;
  T = 298.;
  receptorRhyd = 1e20;
  Vtraj = 10000;
  Vprint = 5000;
  step = 0;
  memset(&center, 0, sizeof(vertex));
  memset(&bounds_min, 0, sizeof(vertex));
  memset(&bounds_max, 0, sizeof(vertex));
  done = false;
  dt_fine = 0.010;
  dt_coarse = 10.000;
  dt_scale_start = 100.;
  dt_scale_end = 500.;

  //ending conditions
  convergence = 0.;
  max_simulations = -1;
}



Model::Model(string inputfn, string outputfn, string logfn) : Model() {
  ifn = inputfn;
  ofn = outputfn;
  lfn = logfn;
  lout.open(lfn, ios_base::out);

  parseInputFile();
  Nthreads = __cilkrts_get_nworkers();

  populateLigands();
  initializeRNG(time(NULL));

}



Model::~Model() {
  for(int k=0; k < Nthreads; k++) 
    vslDeleteStream(&rngCPU[k]);
}



void Model::initializeRNG(int seed) {
  lout << "* Initializing random number generator (N=" << Nthreads << ")" << endl;

  srand(seed);

  rngCPU = (VSLStreamStatePtr*)calloc(Nthreads, sizeof(VSLStreamStatePtr));
  for(int k=0; k < Nthreads; k++) {
    vslNewStream(&rngCPU[k], VSL_BRNG_MT19937, seed+k);
  }

  rand = (vertex*)calloc(2 * ligands.size(), sizeof(vertex));
}



void Model::generateNormal() {
  int blockSize = ((2*ligands.size()*3) / Nthreads);
  cilk_for(int k=0; k < Nthreads; k++) {
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

    // Write trajectory and statistics every 'Vtraj' steps
    if(step % Vtraj == 0) {
      writeCoordinatesPQR();
    }
    if(step % Vprint == 0) {
      t.stop();
      lout << "* Step " << step << " (" << (t.duration/Vtraj) << " s/step)" << endl;
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
