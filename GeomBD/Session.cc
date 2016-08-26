#include "Session.h"
#include "Model.h"
#include "Body.h"
#include "BindingCriteria.h"

///Session definition 
Session::Session(Model *m, SimulationConfig s) {
  id = 0;
  model = m;
  type = s;
  Nreplicates = 0;
  Nbind.set_value(0);
  Nexit.set_value(0);
  Ntlim.set_value(0);
  Davg = 0.;
  done = false;
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


void Session::recordBeta(double beta) {
  beta_history.push_back(beta);
  while(beta_history.size() > 1000.) {
    beta_history.pop_front();
  }
}


void Session::checkConvergence() {
  double m = 0., s = 0.;
  for(int i=0; i < beta_history.size(); i++) {
    m += beta_history[i];
  }
  m /= beta_history.size();
  if(m == 0) return;

  for(int i=0; i < beta_history.size(); i++) {
    s += pow(beta_history[i]-m, 2);
  }
  s = sqrt(s / beta_history.size());
  if(s / m <= model->convergence) {
    done = true;
    for(int i=0; i < ligands.size(); i++) ligands[i]->done = true;
  }
}




//SessionRadial definition
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

  if(Ndone == 0) {
    double kb = 4. * M_PI * b * Davg   *Na*1e12*1e-27;
    model->lout << "   (session " << id << ") ";
    model->lout << "b=" << b << " ";
    model->lout << "kd(b)=" << kb << " ";
    model->lout << "q=" << q << " ";
    model->lout << "Davg=" << Davg << " ";
    model->lout << endl;
    return;
  }

  for(int bsi=0; bsi < bindingCriteria.size(); bsi++) {
    double B = ((double)bindingCriteria[bsi]->Nbind.get_value()) / ((double)Ndone);
    double kb = 4. * M_PI * b * Davg;
    double k = (kb * B) / (1 - ((1 - B)*b/q));
    k *= Na * 1e12 * 1e-27;
    model->lout << "   (session " << id << " bs " << bsi << ") ";
    model->lout << "k_on = " << k << " M⁻¹s⁻¹ ";
    model->lout << "Nbind=" << bindingCriteria[bsi]->Nbind.get_value() << " ";
    model->lout << "Ndone=" << Ndone << " ";
    model->lout << "β=" << B << " ";
    model->lout << "b=" << b << " ";
    model->lout << "kd(b)=" << kb * Na * 1e12 * 1e-27 << " ";
    model->lout << "q=" << q << " ";
    model->lout << "Davg=" << Davg << " ";
    model->lout << endl;
  }

  if(bindingCriteria.size() > 1) {
    double B = ((double)Nbind.get_value()) / ((double)Ndone);
    recordBeta(B);
    double kb = 4. * M_PI * b * Davg;
    double k = (kb * B) / (1 - ((1 - B)*b/q));
    k *= Na * 1e12 * 1e-27;
    model->lout << "   (session " << id << ") ";
    model->lout << "k_on = " << k << " M⁻¹s⁻¹ ";
    model->lout << "Nbind=" << Nbind.get_value() << " ";
    model->lout << "Ndone=" << Ndone << " ";
    model->lout << "β=" << B << " ";
    model->lout << "b=" << b << " ";
    model->lout << "kd(b)=" << kb * Na * 1e12 * 1e-27 << " ";
    model->lout << "q=" << q << " ";
    model->lout << "Davg=" << Davg << " ";
    model->lout << endl;
  }

  if(model->max_simulations > 0 and Ndone >= model->max_simulations) {
    cout << "> Maximum simulations reached. Exiting." << endl;
    model->done = true;
  }
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
    model->lout << "#" << id << "\t Escape event at t=" << bi->t << " ps  (t_dwell=" << bi->t_dwell << "ps, max=" << bi->t_dwell_max << "ps, total=" << bi->t_dwell_total << "ps)" << endl;
    bi->t = 0.;
    bi->t_dwell = 0.;
    bi->t_dwell_max = 0.;
    bi->t_dwell_total = 0.;
  }
}





//SessionAbsolutePeriodic definition
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
    model->lout << endl;
  }

  if(bindingCriteria.size() > 1) {
    double B = ((double)Nbind.get_value()) / ((double)Ndone);
    recordBeta(B);
    double V = bounds.x * bounds.y * bounds.z * LperA3;
    double C = (1. / Na) / V;
    double tavg = t_avgt.get_value() / Nbind.get_value();

    double rate = B * C / (tavg * 1e-12);

    printf("   (session %d)   rate = %.5e Ms⁻¹ (Nbind=%d Ndone=%d β=%.4f C=%.1f k=%.5e s⁻¹ tavg=%.5e Davg=%.5e)", id, rate, Nbind.get_value(), Ndone, B, C, B / (tavg * 1e-12), tavg, Davg);
    model->lout << endl;
  }

  if(model->max_simulations > 0 and Ndone >= model->max_simulations) {
    cout << "> Maximum simulations reached. Exiting." << endl;
    model->done = true;
  }
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
    model->lout << "#" << id << "\t Time-out event  (t_dwell=" << bi->t_dwell << "ps, max=" << bi->t_dwell_max << "ps, total=" << bi->t_dwell_total << "ps)" << endl;
    bi->t = 0.;
    bi->t_dwell = 0.;
    bi->t_dwell_max = 0.;
    bi->t_dwell_total = 0.;
  }
}




//SessionAbsoluteRadial definition
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

    model->lout << "   (session " << id << " bs " << bsi << ") ";
    model->lout << "k_direct = " << k << " s⁻¹ ";
    model->lout << "Nbind=" << bindingCriteria[bsi]->Nbind.get_value() << " ";
    model->lout << "Ndone=" << Ndone << " ";
    model->lout << "t_avg=" << tavg << " ";
    model->lout << "β=" << B << " ";
    model->lout << "q=" << q << " ";
    model->lout << "Davg=" << Davg << " ";
    model->lout << endl;
  }

  if(bindingCriteria.size() > 1) {
    double B = ((double)Nbind.get_value()) / ((double)Ndone);
    recordBeta(B);
    double tavg = t_avgt.get_value() / Nbind.get_value();
    double k = B / (tavg * 1e-12);

    model->lout << "   (session " << id << ") ";
    model->lout << "k_direct = " << k << " s⁻¹ ";
    model->lout << "Nbind=" << Nbind.get_value() << " ";
    model->lout << "Ndone=" << Ndone << " ";
    model->lout << "t_avg=" << tavg << " ";
    model->lout << "β=" << B << " ";
    model->lout << "q=" << q << " ";
    model->lout << "Davg=" << Davg << " ";
    model->lout << endl;
  }

  if(model->max_simulations > 0 and Ndone >= model->max_simulations) {
    cout << "> Maximum simulations reached. Exiting." << endl;
    model->done = true;
  }
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
    model->lout << "#" << id << "\t Escape event at t=" << bi->t << " ps  l=" << sqrt(l2) << " (t_dwell=" << bi->t_dwell << "ps, max=" << bi->t_dwell_max << "ps, total=" << bi->t_dwell_total << "ps)" << endl;
    bi->t = 0.;
    bi->t_dwell = 0.;
    bi->t_dwell_max = 0.;
    bi->t_dwell_total = 0.;
  }
}



