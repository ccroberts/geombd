
#include "Model.h"
#include "Timer.h"



Model::Model() {
  viscosity = 1.0/*cP*/ * 0.001/*cP->Pa.s->J.s/m^3*/ * 2.39e-4 /*J->kcal*/ * 1e12 /*s->ps*/ * 1e-30 /*1/m^3->1/A^3*/ * Na /*kcal->kcal/mol*/;
  T = 298.;
  receptorRoG = 1e20;
  Nthreads = __cilkrts_get_nworkers();
  Vtraj = 10000;
  step = 0;
  center.x = 0.;
  center.y = 0.;
  center.z = 0.;
  done = false;
  dt_fine = 0.010;
  dt_coarse = 10.000;
}



Model::Model(string inputfn, string outputfn) : Model() {
  ifn = inputfn;
  ofn = outputfn;

  parseInputFile();
  populateLigands();
  initializeRNG(time(NULL));
}



Model::~Model() {
  for(int k=0; k < Nthreads; k++) 
    vslDeleteStream(&rngCPU[k]);
}



void Model::initializeRNG(int seed) {
  cout << "* Initializing random number generator (N=" << Nthreads << ")" << ligands.size() << endl;

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

  if(step == 0) writeCoordinatesPQR();

  t.start();

  while(!done) {
    integrate();
    step++;

    if(step % Vtraj == 0) {
      t.stop(); // end step timer

      writeCoordinatesPQR();

      done = true;
      for(int i=0; i < ligands.size(); i++) {
        if(ligands[i]->done == false) {
          done = false;
          break;
        }
      }

      cout << "* Step " << step << " (" << (t.duration/Vtraj) << " s/step)" << endl;
      printRateConstant();

      t.start(); // start step timer
    }
  }

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
