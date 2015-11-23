#include "Session.h"
#include "Model.h"
#include "Body.h"
#include "BindingCriteria.h"

 
Session::Session(Model *m, SimulationConfig s) {
  id = 0;
  model = m;
  type = s;
  Nreplicates = 0;
  Nbind.set_value(0);
  Nexit.set_value(0);
  Ntlim.set_value(0);
  Davg = 0.;
}

Session::~Session() {
}

void Session::populateLigands() {
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
    positionLigand(bi);
    model->ligands.push_back(bi);
    ligands.push_back(bi);
  }

  Davg = 0.;
  for(int i=0; i < conformations.size(); i++) {
    Davg += conformations[i]->D;
  }
  Davg /= conformations.size();
}

void Session::positionLigand(Body *body) {
}

void Session::printRateConstant() {
}





SessionRadial::SessionRadial(Model *m) : Session(m, CONFIGURATION_RADIAL) {
  b = 0.;
  q = 0.;
  q2 = 0.;
}


void SessionRadial::positionLigand(Body *bi) {
   
  vertex Q = { random(-1., 1.), random(-1., 1.), random(-1., 1.) };
  double l = sqrt(Q.x*Q.x + Q.y*Q.y + Q.z*Q.z);
  Q.x *= b / l;
  Q.y *= b / l;
  Q.z *= b / l;
 
  bi->center();
  bi->translate(model->center.x + Q.x, model->center.y + Q.y, model->center.z + Q.z);
  bi->rotate(random(0., M_PI), random(0., M_PI), random(0., M_PI));

}

void SessionRadial::printRateConstant() {
  int Ndone = Nbind.get_value() + Nexit.get_value();
  if(Ndone == 0) return;

  for(int bsi=0; bsi < bindingCriteria.size(); bsi++) {
    double B = ((double)bindingCriteria[bsi]->Nbind.get_value()) / ((double)Ndone);
    double kb = 4. * M_PI * b * Davg;
    double k = (kb * B) / (1 - ((1 - B)*b/q));
    k *= Na * 1e12 * 1e-27;
    printf("   (session %d bs %d)   k_on = %.5e M⁻¹s⁻¹ (Nbind=%d Ndone=%d β=%.4f b=%.1f kd(b)=%.5e q=%.1f Davg=%.5e)", id, bsi, k, bindingCriteria[bsi]->Nbind.get_value(), Ndone, B, b, kb * Na * 1e12 * 1e-27, q, Davg);
    cout << endl;
  }

  double B = ((double)Nbind.get_value()) / ((double)Ndone);
  double kb = 4. * M_PI * b * Davg;
  double k = (kb * B) / (1 - ((1 - B)*b/q));
  k *= Na * 1e12 * 1e-27;
  printf("   (session %d)   k_on = %.5e M⁻¹s⁻¹ (Nbind=%d Ndone=%d β=%.4f b=%.1f kd(b)=%.5e q=%.1f Davg=%.5e)", id, k, Nbind.get_value(), Ndone, B, b, kb * Na * 1e12 * 1e-27, q, Davg);
  cout << endl;
}



SessionAbsolutePeriodic::SessionAbsolutePeriodic(Model *m) : Session(m, CONFIGURATION_ABSOLUTE_PERIODIC) {
  b = 0.;
  t_max = 0.;
  t_avgt.set_value(0.);
}


void SessionAbsolutePeriodic::positionLigand(Body *bi) {
  bi->center();
  bi->translate(start.x, start.y, start.z);
}

void SessionAbsolutePeriodic::printRateConstant() {
  int Ndone = Nbind.get_value() + Ntlim.get_value();
  if(Ndone == 0) return;

  for(int bsi=0; bsi < bindingCriteria.size(); bsi++) {
    double B = ((double)bindingCriteria[bsi]->Nbind.get_value()) / ((double)Ndone);
    double V = bounds.x * bounds.y * bounds.z * LperA3;
    double C = (1. / Na) / V;
    double tavg = bindingCriteria[bsi]->t_avgt.get_value() / bindingCriteria[bsi]->Nbind.get_value();
    double rate = B * C / (tavg * 1e-12);

    printf("   (session %d bs %d)   rate = %.5e Ms⁻¹ (Nbind=%d Ndone=%d β=%.4f C=%.1f k=%.5e s⁻¹ tavg=%.5e Davg=%.5e)", id, bsi, rate, bindingCriteria[bsi]->Nbind.get_value(), Ndone, B, C, B / (tavg * 1e-12), tavg, Davg);
    cout << endl;
  }

  double B = ((double)Nbind.get_value()) / ((double)Ndone);
  double V = bounds.x * bounds.y * bounds.z * LperA3;
  double C = (1. / Na) / V;
  double tavg = t_avgt.get_value() / Nbind.get_value();

  double rate = B * C / (tavg * 1e-12);

  printf("   (session %d)   rate = %.5e Ms⁻¹ (Nbind=%d Ndone=%d β=%.4f C=%.1f k=%.5e s⁻¹ tavg=%.5e Davg=%.5e)", id, rate, Nbind.get_value(), Ndone, B, C, B / (tavg * 1e-12), tavg, Davg);
  cout << endl;
}



SessionAbsoluteRadial::SessionAbsoluteRadial(Model *m) : Session(m, CONFIGURATION_ABSOLUTE_RADIAL) {
  t_avgt.set_value(0.);
}


void SessionAbsoluteRadial::positionLigand(Body *bi) {
  bi->center();
  bi->translate(start.x, start.y, start.z);
}


void SessionAbsoluteRadial::printRateConstant() {
  int Ndone = Nbind.get_value() + Nexit.get_value();
  if(Ndone == 0) return;

  for(int bsi=0; bsi < bindingCriteria.size(); bsi++) {
    double B = ((double)bindingCriteria[bsi]->Nbind.get_value()) / ((double)Ndone);
    double tavg = bindingCriteria[bsi]->t_avgt.get_value() / bindingCriteria[bsi]->Nbind.get_value();
    double k = B / (tavg * 1e-12);

    printf("   (session %d bs %d)   k_direct = %.5e s⁻¹ (Nbind=%d Ndone=%d βdirect=%.4f tavg=%.5e Davg=%.5e)", id, bsi, k, bindingCriteria[bsi]->Nbind.get_value(), Ndone, B, tavg, Davg);
    cout << endl;
  }

  double B = ((double)Nbind.get_value()) / ((double)Ndone);
  double tavg = t_avgt.get_value() / Nbind.get_value();
  double k = B / (tavg * 1e-12);

  printf("   (session %d)   k_direct = %.5e s⁻¹ (Nbind=%d Ndone=%d βdirect=%.4f tavg=%.5e Davg=%.5e)", id, k, Nbind.get_value(), Ndone, B, tavg, Davg);
  cout << endl;
}


