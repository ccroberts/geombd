#include "Session.h"
#include "Model.h"
#include "Body.h"
#include "BindingCriteria.h"

 
Session::Session(Model *m, SimulationConfig s) {
  id = 0;
  model = m;
  type = s;
  Nreplicates = 0;
  Nbind = 0;
  Nexit = 0;
  Ntlim = 0;
  Davg = 0.;
}

Session::~Session() {
}

void Session::populateLigands() {
}

void Session::printRateConstant() {
}





SessionRadial::SessionRadial(Model *m) : Session(m, CONFIGURATION_RADIAL) {
  b = 0.;
  q = 0.;
  q2 = 0.;
}

void SessionRadial::populateLigands() {
  for(int i=0; i < Nreplicates; i++) {
    int rndConf = floor(random(0.0, (double)conformations.size()));
    Body *rc = conformations[rndConf];
    Body *bi = new Body(model, this);
    for(int j=0; j < rc->beads.size(); j++) {
      Bead *nb = new Bead();
      Bead *rb = rc->beads[j];
      nb->R.x = rb->R.x;
      nb->R.y = rb->R.y;
      nb->R.z = rb->R.z;
      nb->q = rb->q;
      nb->r = rb->r;
      nb->m = rb->m;
      nb->type = rc->beads[j]->type;
      bi->beads.push_back(nb);
    }
    bi->define();
   
    vertex Q = { random(-1., 1.), random(-1., 1.), random(-1., 1.) };
    double l = sqrt(Q.x*Q.x + Q.y*Q.y + Q.z*Q.z);
    Q.x *= b / l;
    Q.y *= b / l;
    Q.z *= b / l;
   
    bi->center();
    bi->translate(model->center.x + Q.x, model->center.y + Q.y, model->center.z + Q.z);
    bi->rotate(random(0., M_PI), random(0., M_PI), random(0., M_PI));
    model->ligands.push_back(bi);
    ligands.push_back(bi);
  }

  Davg = 0.;
  for(int i=0; i < conformations.size(); i++) {
    Davg += conformations[i]->D;
  }
  Davg /= conformations.size();
}

void SessionRadial::printRateConstant() {
  int Ndone = Nbind + Nexit + Ntlim;
  if(Ndone == 0) return;

  double B = ((double)Nbind) / ((double)Ndone);
  double kb = 4. * M_PI * b * Davg;
  double k = (kb * B) / (1 - ((1 - B)*b/q));
  k *= Na * 1e12 * 1e-27;
  printf("   (session %d)   k_on = %.5e M⁻¹s⁻¹ (Nbind=%d Ndone=%d β=%.4f b=%.1f kd(b)=%.5e q=%.1f Davg=%.5e)", id, k, Nbind, Ndone, B, b, kb * Na * 1e12 * 1e-27, q, Davg);
  cout << endl;
}



SessionAbsolutePeriodic::SessionAbsolutePeriodic(Model *m) : Session(m, CONFIGURATION_ABSOLUTE_PERIODIC) {
  b = 0.;
  t_max = 0.;
}


void SessionAbsolutePeriodic::populateLigands() {
  for(int i=0; i < Nreplicates; i++) {
    int rndConf = floor(random(0.0, (double)conformations.size()));
    Body *rc = conformations[rndConf];
    Body *bi = new Body(model, this);
    for(int j=0; j < rc->beads.size(); j++) {
      Bead *nb = new Bead();
      Bead *rb = rc->beads[j];
      nb->R.x = rb->R.x;
      nb->R.y = rb->R.y;
      nb->R.z = rb->R.z;
      nb->q = rb->q;
      nb->r = rb->r;
      nb->m = rb->m;
      nb->type = rc->beads[j]->type;
      bi->beads.push_back(nb);
    }
    bi->define();
    bi->center();
    bi->translate(start.x, start.y, start.z);
    model->ligands.push_back(bi);
    ligands.push_back(bi);
  }

  Davg = 0.;
  for(int i=0; i < conformations.size(); i++) {
    Davg += conformations[i]->D;
  }
  Davg /= conformations.size();
}

void SessionAbsolutePeriodic::printRateConstant() {
  int Ndone = Nbind + Nexit + Ntlim;
  if(Ndone == 0) return;

  double B = ((double)Nbind) / ((double)Ndone);
  double V = bounds.x * bounds.y * bounds.z * LperA3;
  double C = (1. / Na) / V;
  double tavg = 0.;

  for(int i=0; i < ligands.size(); i++)
    if(ligands[i]->bound)
      tavg += ligands[i]->t;
  tavg /= Nbind;

  double rate = B * C / (tavg * 1e-12);

  printf("   (session %d)   rate = %.5e Ms⁻¹ (Nbind=%d Ndone=%d β=%.4f C=%.1f k=%.5e s⁻¹ tavg=%.5e Davg=%.5e)", id, rate, Nbind, Ndone, B, C, B / (tavg * 1e-12), tavg, Davg);
  cout << endl;
}



SessionAbsoluteRadial::SessionAbsoluteRadial(Model *m) : Session(m, CONFIGURATION_ABSOLUTE_RADIAL) {
}


void SessionAbsoluteRadial::populateLigands() {
  for(int i=0; i < Nreplicates; i++) {
    int rndConf = floor(random(0.0, (double)conformations.size()));
    Body *rc = conformations[rndConf];
    Body *bi = new Body(model, this);
    for(int j=0; j < rc->beads.size(); j++) {
      Bead *nb = new Bead();
      Bead *rb = rc->beads[j];
      nb->R.x = rb->R.x;
      nb->R.y = rb->R.y;
      nb->R.z = rb->R.z;
      nb->q = rb->q;
      nb->r = rb->r;
      nb->m = rb->m;
      nb->type = rc->beads[j]->type;
      bi->beads.push_back(nb);
    }
    bi->define();
    bi->center();
    bi->translate(start.x, start.y, start.z);
    model->ligands.push_back(bi);
    ligands.push_back(bi);
  }

  Davg = 0.;
  for(int i=0; i < conformations.size(); i++) {
    Davg += conformations[i]->D;
  }
  Davg /= conformations.size();
}

void SessionAbsoluteRadial::printRateConstant() {
  int Ndone = Nbind + Nexit + Ntlim;
  if(Ndone == 0) return;

  double B = ((double)Nbind) / ((double)Ndone);
  double tavg = 0.;

  for(int i=0; i < ligands.size(); i++)
    if(ligands[i]->bound)
      tavg += ligands[i]->t;
  tavg /= Nbind;

  double k = B / (tavg * 1e-12);

  printf("   (session %d)   k_direct = %.5e s⁻¹ (Nbind=%d Ndone=%d βdirect=%.4f tavg=%.5e Davg=%.5e)", id, k, Nbind, Ndone, B, tavg, Davg);
  cout << endl;
}
