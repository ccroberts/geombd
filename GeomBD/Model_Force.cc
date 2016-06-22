
#include "Model.h"
#include "Body.h"
#include "BindingCriteria.h"



void Model::integrate() {

  generateNormal();

  cilk_for(int il=0; il < ligands.size(); il++) {
    Body *Bi = ligands[il];

    //if(! Bi->done) {
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
            if(esmaps[es]->approximate_force(&bi->R, &dF, &E, 2, bi->q)) {
              onGrid = true;
              bi->F.x += dF.x;
              bi->F.y += dF.y;
              bi->F.z += dF.z;
            }
          }
          for(int ds=0; ds < dmaps.size(); ds++) {
            if(dmaps[ds]->approximate_force(&bi->R, &dF, &E, 2, fabs(bi->q))) {
              onGrid = true;
              bi->F.x += dF.x;
              bi->F.y += dF.y;
              bi->F.z += dF.z;
              if(fabs(E) > 0) {
                associated = true;
              }
            }
          }
        }

        for(int tmap=0; tmap < typemaps.size(); tmap++) {
          if(typemaps[tmap]->type == bi->type) {
            if(typemaps[tmap]->approximate_force(&bi->R, &dF, &E, 2, 1.0)) {
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

      // Determine timestep
      double dr[3], radius2, radius;
      dr[0] = Bi->R.x - center.x;
      dr[1] = Bi->R.y - center.y;
      dr[2] = Bi->R.z - center.z;
      radius2 = (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
      if(onGrid) {
        double r_max = 0.5;
        /*
        double mF = vertex_magnitude(Bi->F);
        Bi->dt = min(dt_fine, (r_max*kB*T)/(Bi->D * mF));
        if(Bi->dt == 0. or Bi->dt > dt_fine) Bi->dt = dt_fine;//for first time adaptation
        */
        double mF = vertex_magnitude(Bi->F);
        if(fabs(mF - Bi->mF) > 2. and Bi->dt > 0.0001) {
          Bi->restore();
          Bi->dt /= 2.0;
          if(Bi->dt < 0.0001) Bi->dt = 0.0001;
        } else {
          Bi->mF = mF;
          if(Bi->dt < dt_fine) {
            Bi->dt = min(dt_fine, Bi->dt * 2.0);
          } else {
            Bi->dt = dt_fine;
          }
        }
      } else {
        if(radius2 > dt_scale_start and radius2 < dt_scale_end) {
          double s = (radius2 - dt_scale_start) / (dt_scale_end - dt_scale_start);
          Bi->dt = dt_fine + s * (dt_coarse - dt_fine);
        } else {
          if(radius2 <= dt_scale_start) Bi->dt = dt_fine;
          if(radius2 >= dt_scale_end) Bi->dt = dt_coarse;
        }
      }

      // Backup coordinates in case of interpenetration
      Bi->save();

      // Integrate
      vertex *Si = &rand[il];
      vertex *Sai = &rand[ligands.size()+il];

      vertex dR, dRa;
      double A = sqrt(2. * Bi->D * Bi->dt);
      double dtOVERkBT = Bi->dt / (kB * T);
      double B = Bi->D * dtOVERkBT;
      dR.x = (A * Si->x) + (B * Bi->F.x);
      dR.y = (A * Si->y) + (B * Bi->F.y);
      dR.z = (A * Si->z) + (B * Bi->F.z);
      if(! Bi->translate(dR.x, dR.y, dR.z)) {
        cout << "this shouldn't happen (5A step" << endl;
      }

      double C = sqrt(2 * Bi->Da * Bi->dt);
      double D = Bi->Da * dtOVERkBT;
      dRa.x = (C * Sai->x) + (D * Bi->Fa.x);
      dRa.y = (C * Sai->y) + (D * Bi->Fa.y);
      dRa.z = (C * Sai->z) + (D * Bi->Fa.z);
      Bi->rotate(dRa.x, dRa.y, dRa.z);
      
      // Exclusion check
      bool penetrating = false;
      for(int i=0; i < Bi->beads.size(); i++) {
        Bead *bi = Bi->beads[i];
        for(int ex=0; ex < exmaps.size(); ex++) {
          if(exmaps[ex]->value(&bi->R) > 0) {
            penetrating = true;
            break;
          }
        }
        if(penetrating) {
        }
      }
      // Excluded step, do not increment time
      if(penetrating) {
        Bi->restore();
        continue;
      }

      // Increment time, record dwell-time
      Bi->t += Bi->dt;

      if(onGrid and associated) {
        Bi->t_dwell += Bi->dt;
        Bi->t_dwell_total += Bi->dt;
        if(Bi->t_dwell > Bi->t_dwell_max) Bi->t_dwell_max = Bi->t_dwell;
      } else {
        if(Bi->t_dwell > 0.) {
          if(Bi->t_dwell > Bi->t_dwell_max) Bi->t_dwell_max = Bi->t_dwell;
          Bi->t_dwell = 0.;
        }
      }

      // Session-specific checks
      Bi->session->checkLigand(Bi);

      // Binding criteria check
      for(int bs=0; bs < Bi->session->bindingCriteria.size(); bs++) {
        BindingCriteria* bc = Bi->session->bindingCriteria[bs];
        if(bc->checkBinding(Bi)) {
          *Bi->session->Nbind += 1;
          *Bi->session->t_avgt += Bi->t;
          *bc->Nbind += 1;
          *bc->t_avgt += Bi->t;
          Bi->session->positionLigand(Bi);
          lout << "#" << Bi->session->id << "\t Binding event at t=" << Bi->t << " ps  (t_dwell=" << Bi->t_dwell << "ps, max=" << Bi->t_dwell_max << "ps, total=" << Bi->t_dwell_total << "ps)" << endl;
          Bi->t = 0.;
          Bi->t_dwell = 0.;
          Bi->t_dwell_max = 0.;
          Bi->t_dwell_total = 0.;
        }
      }

  }

  cilk_sync;
}

