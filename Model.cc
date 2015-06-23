
#include "Model.h"
#include "Timer.h"



Model::Model() {
  viscosity = 1.0/*cP*/ * 0.001/*cP->Pa.s->J.s/m^3*/ * 2.39e-4 /*J->kcal*/ * 1e12 /*s->ps*/ * 1e-30 /*1/m^3->1/A^3*/ * Na /*kcal->kcal/mol*/;
  T = 298.;
  t_limit = 1e9;
  bounds.x = 1000.;
  bounds.y = 1000.;
  bounds.z = 1000.;
  r2_escape = 0.;
  receptorRoG = 1e20;
  active = 0;
  Nthreads = 1;
  Nreplicates = 10000;
  step = 0;
  r_ligand = 0.;
  R_ligand = { 0., 0., 0. };
  done = false;
  dt_fine = 0.2;
  dt_coarse = 10.0;
}



Model::Model(string fieldfn, string docklogfn, string outputfn, double rRoG) : Model() {
  fldfn = fieldfn;
  dlgfn = docklogfn;
  ofn = outputfn;
  receptorRoG = rRoG;
}



Model::~Model() {
  for(int k=0; k < Nthreads; k++) 
    vslDeleteStream(&rngCPU[k]);
}



void Model::initializeRNG(int seed) {
  srand(seed);

  rngCPU = (VSLStreamStatePtr*)calloc(Nthreads, sizeof(VSLStreamStatePtr));
  for(int k=0; k < Nthreads; k++) {
    vslNewStream(&rngCPU[k], VSL_BRNG_MT19937, seed+k);
  }

  rand = (vertex*)calloc(2 * ligands.size(), sizeof(vertex));

  active = ligands.size();

  generateNormal();
}



void Model::generateNormal() {
  int blockSize = ((2*ligands.size()*3) / Nthreads) + 1;
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

    if(step % 50000 == 0) {
      t.stop();
      writeCoordinatesPQR();
      printRateConstant();

      done = true;
      active = 0;
      for(int i=0; i < ligands.size(); i++) {
        if(ligands[i]->done == false) {
          done = false;
          active++;
        }
      }
      cout << "Step #" << step << ". Rate: " << (t.duration/50000.) << " s/step. " << active << " simulations still active." << endl;
      t.start();
    }
  }

  printRateConstant();
}


void Model::printRateConstant() {
  double B = 0., Bb = 0., t = 0., tb = 0., k = 0., V = 0., C = 0.;

  // calculate volume and concentration
  if(ligandPosition == LIGAND_POSITION_RADIAL) {
    V = (4. / 3.) * M_PI * pow(r_ligand, 3);
  } else {
    V = bounds.x * bounds.y * bounds.z;
  }
  V *= LperA3;
  C = (1. / Na) / V;

  // calculate bound fraction and average binding time
  for(int i=0; i < ligands.size(); i++) {
    if(ligandPosition == LIGAND_POSITION_RANDOM) {
      if(ligands[i]->bound) {
        Bb += 1.;
        tb += ligands[i]->t;
      }
    }
    if(ligandPosition == LIGAND_POSITION_ABSOLUTE) {
      if(ligands[i]->bound) {
        if(!ligands[i]->bulk) {
          B += 1.;
          t += ligands[i]->t;
        } else {
          Bb += 1.;
          tb += ligands[i]->t;
        }
      }
    }
    if(ligandPosition == LIGAND_POSITION_RADIAL) {
      if(ligands[i]->bound) {
        B += 1.;
        t += ligands[i]->t;
      }
      if(ligands[i]->bulk) {
        Bb += 1;
      }
    }
  }

  //Association constants
  if(ligandPosition == LIGAND_POSITION_ABSOLUTE) {
    if(B > 0.) {
      t *= 1e-12 / B; //average binding time, converted to seconds
      B /= ligands.size();
      k = B * (1. / (t));
      printf("k_iet:   %.5e s⁻¹ (β = %.4f, t_avg = %.4e)\n", k, B, t);
    }
    if(Bb > 0.) {
      tb *= 1e-12 / Bb;
      Bb /= ligands.size();
      k = Bb * (1. / (tb * C * C));
      printf("k_bulk:     %.5e M⁻²s⁻¹ (β = %.4f, t_avg = %.4e s, rate = %.5e s⁻¹)\n", k, Bb, tb, Bb*(1. / tb));
    }
  }
  if(ligandPosition == LIGAND_POSITION_RANDOM and Bb > 0.) {
    tb *= 1e-12 / Bb;
    Bb /= ligands.size();
    k = Bb * (1. / (tb * C * C));
    printf("k_bulk:     %.5e M⁻²s⁻¹ (β = %.4f, t_avg = %.4e s, rate = %.5e s⁻¹)\n", k, Bb, tb, Bb*(1. / tb));
  }
  if(ligandPosition == LIGAND_POSITION_RADIAL and B > 0.) {
    double kb = 4. * M_PI * r_ligand * ligands[0]->D;
    double kd = 4. * M_PI * 5.*r_ligand * ligands[0]->D;
    B /= ligands.size();
    k = (kd * B) / (1 - ((1 - B)*kb/kd));
    k *= Na * 1e12 * 1e-27;
    printf("k_on = %.5e M⁻¹s⁻¹ (β = %.4f)\n", k, B);
  }
}


