
#include "Body.h"
#include "Strings.h"
#include "Model.h"
#include "Session.h"

 
Body::Body() {
  model = NULL;
  session = NULL;

  t = 0.;
  dt = 0.100;

  t_dwell = 0.;
  t_dwell_max = 0.;
  t_dwell_total = 0.;

  r = 0.;
  r_max = 0.;
}


 
Body::Body(Model *m, Session *s) : Body() {
  model = m;
  session = s;
}


 
Body::~Body() {
  beads.clear();
}



void Body::define() {
  m = 0;
  I = 0;
  D = 0;
  Da = 0;
  R.x = 0;
  R.y = 0;
  R.z = 0;
  r = 0;
  r_max = 0;

  N = beads.size();

  // Define properties
  if(N == 1) {
    Bead *bi = beads[0];
    m = bi->m;
    R.x = bi->R.x;
    R.y = bi->R.y;
    R.z = bi->R.z;
    I = 0;
    r = bi->r;
    r_max = bi->r;
  } else {
    //mass and center of mass
    for(int i=0; i < N; i++) {
      Bead *bi = beads[i];
      m += bi->m;
      R.x += bi->R.x * bi->m;
      R.y += bi->R.y * bi->m;
      R.z += bi->R.z * bi->m;
    }

    R.x /= m;
    R.y /= m;
    R.z /= m;

    vertex dr;
    //moment of intertia
    for(int i=0; i < N; i++) {
      Bead *bi = beads[i];

      dr.x = bi->R.x - R.x;
      dr.y = bi->R.y - R.y;
      dr.z = bi->R.z - R.z;
      double distSqr = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

      I += distSqr * bi->m;
      r += distSqr;

      if(sqrt(distSqr) + bi->r > r_max) r_max = sqrt(distSqr) + bi->r;
    }

    r = sqrt(r / N);
    /*
    double invSumRij = 0.;
    //hydrodynamic radius
    for(int i=0; i < N; i++) {
      Bead *bi = beads[i];
      for(int j=i+1; j < N; j++) {
        Bead *bj = beads[j];
        dr.x = bi->R.x - bj->R.x;
        dr.y = bi->R.y - bj->R.y;
        dr.z = bi->R.z - bj->R.z;
        invSumRij += 1. / sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
      }
    }

    r = 1.0 / (invSumRij / (N * N));
    */
  }

  if(model != NULL) {
    // Define Dco
    double Pl = (kB * model->T) / (6. * M_PI * model->viscosity);
    double Pa = (kB * model->T) / (8. * M_PI * model->viscosity);
    D  = (Pl / r) + (Pl / model->receptorRhyd);
    Da = (Pa / pow(r, 3)) + (Pa / pow(model->receptorRhyd, 3));
  } else {
    model->lout << "! Warning: You're using a Body object outside the context of a Model. Diffusion coefficients not calculated." << endl;
  }

  //lout << ">> m=" << m << " I=" << I << " r=" << r << " r_max=" << r_max << " D=" << D << " Da=" << Da << endl;
}


void Body::save() {
  _R.x = R.x;
  _R.y = R.y;
  _R.z = R.z;
  _Ra.x = Ra.x;
  _Ra.y = Ra.y;
  _Ra.z = Ra.z;

  for(int i=0; i < beads.size(); i++) {
    beads[i]->saveR();
  }
}


void Body::restore() {
  R.x = _R.x;
  R.y = _R.y;
  R.z = _R.z;
  Ra.x = _Ra.x;
  Ra.y = _Ra.y;
  Ra.z = _Ra.z;

  for(int i=0; i < beads.size(); i++) {
    beads[i]->restoreR();
  }
}


 
bool Body::translate(double dx, double dy, double dz, bool suppressWarning) {
  R.x += dx;
  R.y += dy;
  R.z += dz;

  for(int i=0; i < N; i++) {
    Bead *bi = beads[i];
    bi->translate(dx, dy, dz);
  }

  double jump = dx*dx + dy*dy + dz*dz;
  if(jump > 25.0 and not suppressWarning) {
    model->lout << "! Warning: Body translation greater than 5.0A in a single step." << endl;
    return false;
  }

  return true;
}


 
void Body::center() {
  translate(-R.x, -R.y, -R.z, true);
}


 
void Body::rotate(double dax, double day, double daz) {
  Ra.x += dax;
  Ra.y += day;
  Ra.z += daz;

  for(int i=0; i < N; i++) {
    Bead *bi = beads[i];
    bi->rotateAbout(R.x, R.y, R.z, dax, day, daz);
  }
}



void writePDBBead(Bead *bi, unt index, char chain, fstream &outf) {
  outf << "ATOM  ";
  outf.width(5);
  outf << right << index;
  outf << "      ";
  outf.width(3);
  outf << " P " << " ";
  outf << chain;
  outf << "        ";
  outf.width(8);
  outf.precision(3);
  outf << fixed << bi->R.x/1.;
  outf.width(8);
  outf.precision(3);
  outf << fixed << bi->R.y/1.;
  outf.width(8);
  outf.precision(3);
  outf << fixed << bi->R.z/1.;
  outf << " ";
  outf.setf(ios::right, ios::adjustfield);
  outf.precision(4);
  outf.width(7);
  outf << bi->q << " " << (bi->r/1.);
  outf << endl;
}

void Body::writePDB(fstream &outf, char chain) {
  outf << "REMARK ";
  outf << endl;

  for(int i=0; i < beads.size(); i++) {
    Bead *bi = beads[i];
    writePDBBead(bi, i+1, chain, outf);
  }
}


