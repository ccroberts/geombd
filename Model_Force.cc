
#include "Model.h"
#include "Body.h"
#include "MartiniNB.h"



void Model::integrate() {
  generateNormal();

  cilk_for(int il=0; il < ligands.size(); il++) {
    Body *Bi = ligands[il];

    if(! Bi->done) {
      bool onGrid = false;
      double E;

      // Calculate forces
      for(int i=0; i < Bi->beads.size(); i++) {
        Bead *bi = Bi->beads[i];
        bi->F.x = 0.;
        bi->F.y = 0.;
        bi->F.z = 0.;
        
        if(bi->q != 0.) {
          for(int es=0; es < esmaps.size(); es++) {
            if(esmaps[es]->force(&bi->R, &bi->F, bi->q, &E)) {
              onGrid = true;
            }
          }
        }

        for(int tmap=0; tmap < typemaps.size(); tmap++) {
          if(typemaps[tmap]->type == bi->type) {
            if(typemaps[tmap]->force(&bi->R, &bi->F, &E)) {
              onGrid = true;
            }
          }
        }
      }


      Bi->F.x = 0.;
      Bi->F.y = 0.;
      Bi->F.z = 0.;
      Bi->Fa.x = 0.;
      Bi->Fa.y = 0.;
      Bi->Fa.z = 0.;

      // Propogate bead forces to body
      for(int k=0; k < Bi->beads.size(); k++) {
        Bead *bk = Bi->beads[k];

        if(bk->F.x == 0. and bk->F.y == 0. and bk->F.z == 0.) continue;

        Bi->F.x += bk->F.x;
        Bi->F.y += bk->F.y;
        Bi->F.z += bk->F.z;

        vertex A = { bk->R.x - Bi->R.x, bk->R.y - Bi->R.y, bk->R.z - Bi->R.z }, B;
        B.x = A.y * bk->F.z - A.z * bk->F.y;
        B.y = A.z * bk->F.x - A.x * bk->F.z;
        B.z = A.x * bk->F.y - A.y * bk->F.x;

        Bi->Fa.x += B.x;
        Bi->Fa.y += B.y;
        Bi->Fa.z += B.z;
      }

      double dt = (onGrid) ? dt_fine : dt_coarse;
      double dr[3], l2;
      dr[0] = Bi->R.x - bindingSite.x;
      dr[1] = Bi->R.y - bindingSite.y;
      dr[2] = Bi->R.z - bindingSite.z;
      l2 = (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
      if(l2 < bindingSite.r2) {
        cout << "bind: " << Bi->t << endl;
        Bi->bound = true;
        Bi->done = true;
      }
      if(ligandPosition == LIGAND_POSITION_RADIAL and l2 >= r2_escape) {
        Bi->bulk = true;
        Bi->done = true;
      }

      double dtOVERkBT = dt / (kB * T);

      // Integrate
      vertex *Si = &rand[il];
      vertex *Sai = &rand[ligands.size()+il];

      vertex dR, dRa;
      double A = sqrt(2. * Bi->D * dt);
      double B = Bi->D * dtOVERkBT;

      dR.x = (A * Si->x) + (B * Bi->F.x);
      dR.y = (A * Si->y) + (B * Bi->F.y);
      dR.z = (A * Si->z) + (B * Bi->F.z);

      Bi->translate(dR.x, dR.y, dR.z);

      double C = sqrt(2 * Bi->Da * dt);
      double D = Bi->Da * dtOVERkBT;

      dRa.x = (C * Sai->x) + (D * Bi->Fa.x);
      dRa.y = (C * Sai->y) + (D * Bi->Fa.y);
      dRa.z = (C * Sai->z) + (D * Bi->Fa.z);

      Bi->rotate(dRa.x, dRa.y, dRa.z);

      Bi->t += dt;
      if(Bi->t > t_limit) Bi->done = true;

      if(ligandPosition != LIGAND_POSITION_RADIAL) {
        if(Bi->R.x > 0.5*bounds.x)  { Bi->translate(-bounds.x, 0., 0.); Bi->bulk = true; }
        if(Bi->R.y > 0.5*bounds.y)  { Bi->translate(0., -bounds.y, 0.); Bi->bulk = true; }
        if(Bi->R.z > 0.5*bounds.z)  { Bi->translate(0., 0., -bounds.z); Bi->bulk = true; }
        if(Bi->R.x < -0.5*bounds.x) { Bi->translate(bounds.x, 0., 0.);  Bi->bulk = true; }
        if(Bi->R.y < -0.5*bounds.y) { Bi->translate(0., bounds.y, 0.);  Bi->bulk = true; }
        if(Bi->R.z < -0.5*bounds.z) { Bi->translate(0., 0., bounds.z);  Bi->bulk = true; }
      }
    }
  }

  cilk_sync;
}


