
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
  dt_fine = 0.005;
  dt_coarse = 0.500;
  cout << "> Timesteps: Fine="<<dt_fine<<" Coarse=" <<dt_coarse<<endl;
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

    if(step % 20000 == 0) {
      t.stop();
      writeCoordinatesPQR();
      printRateConstant();

      done = true;
      active = 0;
      double within500 = 0;
      double within200 = 0;
      for(int i=0; i < ligands.size(); i++) {
        if(ligands[i]->done == false) {
          done = false;
          active++;
        }
        if(ligands[i]->bound == false) {
          vertex dr;
          dr.x = bindingSite.x - ligands[i]->R.x;
          dr.y = bindingSite.y - ligands[i]->R.y;
          dr.z = bindingSite.z - ligands[i]->R.z;
          double l = sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
          if(l < 500.) within500++;
          if(l < 200.) within200++;
        }
      }
      within500 /= ligands.size();
      within200 /= ligands.size();
      cout << "> Step #" << step << ". Rate: " << (t.duration/20000.) << " s/step. " << active << " simulations still active. ";
      printf("L<500=%f \tL<200=%f", within500, within200);
      cout << endl;
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
    if(ligands[i]->bound) {
      B += 1.;
      t += ligands[i]->t;
    }
  }

  //Association constants
  if(ligandPosition == LIGAND_POSITION_ABSOLUTE and B > 0.) {
    double kb = 4. * M_PI * r_ligand * ligands[0]->D;
    double kd = 4. * M_PI * 5.*r_ligand * ligands[0]->D;

    t *= 1e-12 / B; //average binding time, converted to seconds
    B /= ligands.size();
    k = B * (1. / (t));
    printf("k_iet:   %.5e s⁻¹ (β = %.4f, t_avg = %.4e)\n", k, B, t);

    k = (kd * B) / (1 - ((1 - B)*kb/kd));
    k *= Na * 1e12 * 1e-27;

    printf("k_on = %.5e M⁻¹s⁻¹ (β = %.4f, k_b = %.5e M⁻¹s⁻¹, k_d = %.5e M⁻¹s⁻¹)\n", k, B, kb * Na * 1e12 * 1e-27, kd * Na * 1e12 * 1e-27);
  }
  if(ligandPosition == LIGAND_POSITION_RADIAL and B > 0.) {
    double kb = 4. * M_PI * r_ligand * ligands[0]->D;
    double kd = 4. * M_PI * 5.*r_ligand * ligands[0]->D;
    B /= ligands.size();
    k = (kd * B) / (1 - ((1 - B)*kb/kd));
    k *= Na * 1e12 * 1e-27;
    printf("k_on = %.5e M⁻¹s⁻¹ (β = %.4f, k_b = %.5e M⁻¹s⁻¹, k_d = %.5e M⁻¹s⁻¹)\n", k, B, kb * Na * 1e12 * 1e-27, kd * Na * 1e12 * 1e-27);
  }
}


