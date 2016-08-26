
#include "Body.h"
#include "Strings.h"
#include "Model.h"
#include "Session.h"

 
Body::Body() {
  model = NULL;
  session = NULL;

  N = 0;

  m = 0;
  I = 0;
  D = 0;
  Da = 0;
  R.x = Ra.y = Ra.z = 0.;
  Ra.x = Ra.y = Ra.z = 0.;
  F.x = Fa.y = Fa.z = 0.;
  Fa.x = Fa.y = Fa.z = 0.;
  mF = 0.;
  r = 0;
  r_max = 0;

  t = 0.;
  dt = 0.;

  t_dwell = 0.;
  t_dwell_max = 0.;
  t_dwell_total = 0.;

  r = 0.;
  r_max = 0.;

  mF = 0.;

  done = false;
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
  R.x = Ra.y = Ra.z = 0.;
  Ra.x = Ra.y = Ra.z = 0.;
  F.x = Fa.y = Fa.z = 0.;
  Fa.x = Fa.y = Fa.z = 0.;
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
    D  = (Pl / r) + (Pl / model->receptor_radius);
    Da = (Pa / pow(r, 3)) + (Pa / pow(model->receptor_radius, 3));
  } else {
    model->lout << "! Warning: You're using a Body object outside the context of a Model. Diffusion coefficients not calculated." << endl;
  }

  //lout << ">> R=(" << R.x << ", " << R.y << ", " << R.z << ") m=" << m << " I=" << I << " r=" << r << " r_max=" << r_max << " D=" << D << " Da=" << Da << endl;
}


void Body::save() {
  _t = t;
  _R.x = R.x;
  _R.y = R.y;
  _R.z = R.z;
  _Ra.x = Ra.x;
  _Ra.y = Ra.y;
  _Ra.z = Ra.z;
  _F.x = F.x;
  _F.y = F.y;
  _F.z = F.z;
  _Fa.x = Fa.x;
  _Fa.y = Fa.y;
  _Fa.z = Fa.z;

  for(int i=0; i < beads.size(); i++) {
    beads[i]->save();
  }
}


void Body::restore() {
  t = _t;
  R.x = _R.x;
  R.y = _R.y;
  R.z = _R.z;
  Ra.x = _Ra.x;
  Ra.y = _Ra.y;
  Ra.z = _Ra.z;
  F.x = _F.x;
  F.y = _F.y;
  F.z = _F.z;
  Fa.x = _Fa.x;
  Fa.y = _Fa.y;
  Fa.z = _Fa.z;

  for(int i=0; i < beads.size(); i++) {
    beads[i]->restore();
  }
}


 
void Body::translate(double dx, double dy, double dz) {
  R.x += dx;
  R.y += dy;
  R.z += dz;

  for(int i=0; i < N; i++) {
    Bead *bi = beads[i];
    bi->translate(dx, dy, dz);
  }
}


 
void Body::center() {
  translate(-R.x, -R.y, -R.z);
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
  outf << fixed << bi->R.x;
  outf.width(8);
  outf.precision(3);
  outf << fixed << bi->R.y;
  outf.width(8);
  outf.precision(3);
  outf << fixed << bi->R.z;
  outf << " ";
  outf.setf(ios::right, ios::adjustfield);
  outf.precision(4);
  outf.width(7);
  outf << bi->q << " " << (bi->r);
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


