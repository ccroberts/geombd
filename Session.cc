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
    printf("   (session %d bs %d)   k_on = %.5e M⁻¹s⁻¹ ± %.1f%% (Nbind=%d Ndone=%d β=%.4f b=%.1f kd(b)=%.5e q=%.1f Davg=%.5e)", id, bsi, k, 100.*pow(Nbind.get_value(), -0.5), bindingCriteria[bsi]->Nbind.get_value(), Ndone, B, b, kb * Na * 1e12 * 1e-27, q, Davg);
    cout << endl;
  }

  double B = ((double)Nbind.get_value()) / ((double)Ndone);
  double kb = 4. * M_PI * b * Davg;
  double k = (kb * B) / (1 - ((1 - B)*b/q));
  k *= Na * 1e12 * 1e-27;
  printf("   (session %d)   k_on = %.5e M⁻¹s⁻¹ ± %.1f%% (Nbind=%d Ndone=%d β=%.4f b=%.1f kd(b)=%.5e q=%.1f Davg=%.5e)", id, k, 100.*pow(Nbind.get_value(), -0.5), Nbind.get_value(), Ndone, B, b, kb * Na * 1e12 * 1e-27, q, Davg);
  cout << endl;
}


void SessionRadial::checkLigand(Body *bi) {
  double dr[3], l2;
  dr[0] = bi->R.x - model->center.x;
  dr[1] = bi->R.y - model->center.y;
  dr[2] = bi->R.z - model->center.z;
  l2 = (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
  if(l2 >= q2) {
    //bi->done = true;
    *Nexit += 1;
    bi->session->positionLigand(bi);
    cout << "#" << id << "\t Escape event at t=" << bi->t << " ps  (t_dwell=" << bi->t_dwell << "ps, max=" << bi->t_dwell_max << "ps, total=" << bi->t_dwell_total << "ps)" << endl;
    bi->t = 0.;
    bi->t_dwell = 0.;
    bi->t_dwell_max = 0.;
    bi->t_dwell_total = 0.;
  }
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


void SessionAbsolutePeriodic::checkLigand(Body *bi) {
  if(bi->R.x > 0.5*bounds.x)  { bi->translate(-bounds.x, 0., 0.); }
  if(bi->R.y > 0.5*bounds.y)  { bi->translate(0., -bounds.y, 0.); }
  if(bi->R.z > 0.5*bounds.z)  { bi->translate(0., 0., -bounds.z); }
  if(bi->R.x < -0.5*bounds.x) { bi->translate(bounds.x, 0., 0.);  }
  if(bi->R.y < -0.5*bounds.y) { bi->translate(0., bounds.y, 0.);  }
  if(bi->R.z < -0.5*bounds.z) { bi->translate(0., 0., bounds.z);  }

  if(bi->t >= t_max) {
    //bi->done = true;
    *Ntlim += 1;
    positionLigand(bi);
    cout << "#" << id << "\t Time-out event  (t_dwell=" << bi->t_dwell << "ps, max=" << bi->t_dwell_max << "ps, total=" << bi->t_dwell_total << "ps)" << endl;
    bi->t = 0.;
    bi->t_dwell = 0.;
    bi->t_dwell_max = 0.;
    bi->t_dwell_total = 0.;
  }
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

    printf("   (session %d bs %d)   k_direct = %.5e s⁻¹ ± %.1f%% (Nbind=%d Ndone=%d βdirect=%.4f tavg=%.5e Davg=%.5e)", id, bsi, k, 100.*pow(Nbind.get_value(), -0.5), bindingCriteria[bsi]->Nbind.get_value(), Ndone, B, tavg, Davg);
    cout << endl;
  }

  double B = ((double)Nbind.get_value()) / ((double)Ndone);
  double tavg = t_avgt.get_value() / Nbind.get_value();
  double k = B / (tavg * 1e-12);

  printf("   (session %d)   k_direct = %.5e s⁻¹ ± %.1f%% (Nbind=%d Ndone=%d βdirect=%.4f tavg=%.5e Davg=%.5e)", id, k, 100.*pow(Nbind.get_value(), -0.5), Nbind.get_value(), Ndone, B, tavg, Davg);
  cout << endl;
}


void SessionAbsoluteRadial::checkLigand(Body *bi) {
  double dr[3], l2;
  dr[0] = bi->R.x - model->center.x;
  dr[1] = bi->R.y - model->center.y;
  dr[2] = bi->R.z - model->center.z;
  l2 = (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);

  if(l2 >= q2) {
    //bi->done = true;
    *Nexit += 1;
    bi->session->positionLigand(bi);
    cout << "#" << id << "\t Escape event at t=" << bi->t << " ps  l=" << sqrt(l2) << " (t_dwell=" << bi->t_dwell << "ps, max=" << bi->t_dwell_max << "ps, total=" << bi->t_dwell_total << "ps)" << endl;
    bi->t = 0.;
    bi->t_dwell = 0.;
    bi->t_dwell_max = 0.;
    bi->t_dwell_total = 0.;
  }
}


SessionMilestone::SessionMilestone(Model *m) : Session(m, CONFIGURATION_MILESTONE) {
  state = NULL;
  spacing = 100.;
  reaction = 15.;
}

void SessionMilestone::positionLigand(Body *bi) {
  vertex Q = { random(-1., 1.), random(-1., 1.), random(-1., 1.) };
  double l = sqrt(Q.x*Q.x + Q.y*Q.y + Q.z*Q.z);
  Q.x *= state->r0 / l;
  Q.y *= state->r0 / l;
  Q.z *= state->r0 / l;
 
  bi->center();
  bi->translate(model->center.x + Q.x, model->center.y + Q.y, model->center.z + Q.z);
  bi->rotate(random(0., M_PI), random(0., M_PI), random(0., M_PI));
}


void SessionMilestone::printRateConstant() {
  const std::list<double> &tfwds = state->tfwd.get_value();
  int Nfwd = tfwds.size();
  const std::list<double> &tbcks = state->tbck.get_value();
  int Nbck = tbcks.size();

  printf("   (session %d)   Sampled %d progressions (%d fwd, %d bck)", id, Nfwd+Nbck, Nfwd, Nbck);
  cout << endl;
}


void SessionMilestone::checkLigand(Body *bi) {
  double dr[3], l2, l;
  dr[0] = bi->R.x - model->center.x;
  dr[1] = bi->R.y - model->center.y;
  dr[2] = bi->R.z - model->center.z;
  l2 = (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
  l = sqrt(l2);

  if(l >= state->r0 + 0.5*spacing) {
    state->tbck->push_back(bi->t);    
    positionLigand(bi);
  } else if(l <= state->r0 - 0.5*spacing) {
    state->tfwd->push_back(bi->t);    
    positionLigand(bi);
  }
}


