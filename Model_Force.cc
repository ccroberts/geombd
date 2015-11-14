
#include "Model.h"
#include "Body.h"
#include "BindingCriteria.h"



void Model::integrate() {

  generateNormal();

  cilk_for(int il=0; il < ligands.size(); il++) {
    Body *Bi = ligands[il];

    if(! Bi->done) {
      bool onGrid = false, associated = false;
      double E;

      // Calculate forces
      for(int i=0; i < Bi->beads.size(); i++) {
        Bead *bi = Bi->beads[i];
        bi->F.x = 0.;
        bi->F.y = 0.;
        bi->F.z = 0.;
        vertex dF;
        dF.x = 0.;
        dF.y = 0.;
        dF.z = 0.;
        
        if(bi->q != 0.) {
          for(int es=0; es < esmaps.size(); es++) {
            if(esmaps[es]->force(&bi->R, &dF, bi->q, &E)) {
              onGrid = true;
              bi->F.x += dF.x;
              bi->F.y += dF.y;
              bi->F.z += dF.z;
            }
          }
          for(int apbs=0; apbs < apbsmaps.size(); apbs++) {
            if(apbsmaps[apbs]->force(&bi->R, &dF, bi->q, &E)) {
              onGrid = true;
              bi->F.x += dF.x;
              bi->F.y += dF.y;
              bi->F.z += dF.z;
            }
          }
        }

        for(int tmap=0; tmap < typemaps.size(); tmap++) {
          if(typemaps[tmap]->type == bi->type) {
            if(typemaps[tmap]->force(&bi->R, &dF, &E)) {
              onGrid = true;
              if(fabs(E) > 0) {
                associated = true;
              }
              bi->F.x += dF.x;
              bi->F.y += dF.y;
              bi->F.z += dF.z;
            }
          }
        }
      }

      // determine timestep
      double dt = (onGrid) ? dt_fine : dt_coarse;
      double dtOVERkBT = dt / (kB * T);


      // Propogate bead forces to body
      Bi->F.x = 0.;
      Bi->F.y = 0.;
      Bi->F.z = 0.;
      Bi->Fa.x = 0.;
      Bi->Fa.y = 0.;
      Bi->Fa.z = 0.;

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

      bool valid = true;
      for(int i=0; i < Bi->beads.size(); i++) {
        Bead *bi = Bi->beads[i];
        vertex tR = { bi->R.x + dR.x, bi->R.y + dR.y, bi->R.z + dR.z };

        for(int ex=0; ex < xmaps.size(); ex++) {
          ExclusionMap *em = xmaps[ex];
          if(em->onGrid(&tR)) {
            short v = xmaps[ex]->value(&tR);
            if(v == 1)  {
              Bi->translate(-dR.x, -dR.y, -dR.z);
              Bi->rotate(-dRa.x, -dRa.y, -dRa.z);
              valid = false;
              break;
            }
          }
        }
        if(!valid) break;
      }
      if(valid) {
        // increment time
        Bi->t += dt;
        // record dwell-time
        if(onGrid and associated) {
          Bi->t_dwell += dt;
          Bi->t_dwell_total += dt;
          if(Bi->t_dwell > Bi->t_dwell_max) Bi->t_dwell_max = Bi->t_dwell;
        } else {
          if(Bi->t_dwell > 0.) {
            if(Bi->t_dwell > Bi->t_dwell_max) Bi->t_dwell_max = Bi->t_dwell;
            Bi->t_dwell = 0.;
          }
        }
      }


      if(Bi->session->type == CONFIGURATION_ABSOLUTE_PERIODIC) {
        // check for periodic wrapping
        SessionAbsolutePeriodic *sap = (SessionAbsolutePeriodic*)Bi->session;
        if(Bi->R.x > 0.5*sap->bounds.x)  { Bi->translate(-sap->bounds.x, 0., 0.); }
        if(Bi->R.y > 0.5*sap->bounds.y)  { Bi->translate(0., -sap->bounds.y, 0.); }
        if(Bi->R.z > 0.5*sap->bounds.z)  { Bi->translate(0., 0., -sap->bounds.z); }
        if(Bi->R.x < -0.5*sap->bounds.x) { Bi->translate(sap->bounds.x, 0., 0.);  }
        if(Bi->R.y < -0.5*sap->bounds.y) { Bi->translate(0., sap->bounds.y, 0.);  }
        if(Bi->R.z < -0.5*sap->bounds.z) { Bi->translate(0., 0., sap->bounds.z);  }

        // check for timed-out ligands
        if(Bi->t >= sap->t_max) {
          cout << "#" << Bi->session->id << "\t Time-out event" << endl;
          Bi->done = true;
          *Bi->session->Ntlim += 1;
        }
      }

      // check for binding
      for(int bs=0; bs < Bi->session->bindingCriteria.size(); bs++) {
        BindingCriteria* bc = Bi->session->bindingCriteria[bs];
        if(bc->checkBinding(Bi)) {
          cout << "#" << Bi->session->id << "\t Binding event at t=" << Bi->t << " ps" << endl;
          Bi->bound = true;
          Bi->done = true;
          *Bi->session->Nbind += 1;
          *bc->Nbind += 1;
        }
      }

      // check for escape
      double dr[3], l2;
      if(Bi->session->type == CONFIGURATION_RADIAL or Bi->session->type == CONFIGURATION_ABSOLUTE_RADIAL) {
        dr[0] = Bi->R.x - center.x;
        dr[1] = Bi->R.y - center.y;
        dr[2] = Bi->R.z - center.z;
        l2 = (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        if(Bi->session->type == CONFIGURATION_RADIAL) {
          SessionRadial *sr = (SessionRadial*)Bi->session;
          if(l2 >= sr->q2) {
            Bi->done = true;
            *sr->Nexit += 1;
            cout << "#" << Bi->session->id << "\t Escape event at t=" << Bi->t << " ps" << endl;
          }
        }
        if(Bi->session->type == CONFIGURATION_ABSOLUTE_RADIAL) {
          SessionAbsoluteRadial *sar = (SessionAbsoluteRadial*)Bi->session;
          if(l2 >= sar->q2) {
            Bi->done = true;
            *sar->Nexit += 1;
            cout << "#" << Bi->session->id << "\t Escape event at t=" << Bi->t << " ps" << endl;
          }
        }
      }
    }
  }

  cilk_sync;
}


