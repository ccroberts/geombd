
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
            if(esmaps[es]->approximate_force(&bi->R, &dF, &E, fd_order, bi->q)) {
              onGrid = true;
              bi->F.x += dF.x;
              bi->F.y += dF.y;
              bi->F.z += dF.z;
            }
          }
          for(int ds=0; ds < dmaps.size(); ds++) {
            if(dmaps[ds]->approximate_force(&bi->R, &dF, &E, fd_order, fabs(bi->q))) {
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
            if(typemaps[tmap]->approximate_force(&bi->R, &dF, &E, fd_order, 1.0)) {
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

      // Calculate force magnitude
      vertex Rcom = { Bi->R.x - center.x, Bi->R.y - center.y, Bi->R.z - center.z };
      double Rcom_mag = sqrt(Rcom.x*Rcom.x + Rcom.y*Rcom.y + Rcom.z*Rcom.z);
      if(Rcom_mag != 0.) {
        Rcom.x /= Rcom_mag;
        Rcom.y /= Rcom_mag;
        Rcom.z /= Rcom_mag;
      }
      vertex Fi = Bi->F;
      double Fi_mag = sqrt(Fi.x*Fi.x + Fi.y*Fi.y + Fi.z*Fi.z);
      if(Fi_mag != 0.) {
        Fi.x /= Fi_mag;
        Fi.y /= Fi_mag;
        Fi.z /= Fi_mag;
      }
      double Rcom_dot_Fi = Rcom.x*Fi.x + Rcom.y*Fi.y + Rcom.z*Fi.z;
      double New_mF = Rcom_dot_Fi * Fi_mag;

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

      if(radius2 > dt_scale_start and radius2 < dt_scale_end) {
        double s = (radius2 - dt_scale_start) / (dt_scale_end - dt_scale_start);
        Bi->dt = dt_fine + s * (dt_coarse - dt_fine);
      } else {
        if(radius2 <= dt_scale_start) {
          Bi->dt = dt_fine;
          if(New_mF < -1.0) {
            Bi->dt /= (New_mF / -1.0);
          }
        }
        if(radius2 >= dt_scale_end) Bi->dt = dt_coarse;
      }

      // Backup coordinates in case of interpenetration
      if(exmaps.size() > 0)
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
      Bi->translate(dR.x, dR.y, dR.z);

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
          }
        }
        if(penetrating) break;
      }
      if(penetrating) {
        Bi->restore();
        continue;
      }
      
      // Debug check
      if(debug_map) {
        bool debug_penetrate = false;
        for(int i=0; i < Bi->beads.size(); i++) {
          Bead *bi = Bi->beads[i];
          if(debug_map->value(&bi->R) > 0) {
            debug_penetrate = true;
          }
          if(debug_penetrate) break;
        }
        if(debug_penetrate) {
          Bi->mF = New_mF;
        } else {
          Bi->mF = 0.;
        }
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
          // Should we write the bound conformation?
          if(writeBinders) {
            fstream boutf(bfn.c_str(), ios::out | ios::app);
            Bi->writePDB(boutf, 'A');
            boutf << "TER" << endl;
            boutf.close();
          }
          if(logBinders)
            lout << "#" << Bi->session->id << "\t Binding event at t=" << Bi->t << " ps  (t_dwell=" << Bi->t_dwell << "ps, max=" << Bi->t_dwell_max << "ps, total=" << Bi->t_dwell_total << "ps)" << endl;
          Bi->bound = true;
          *Bi->session->t_avgt += Bi->t;
          *bc->Nbind += 1;
          *bc->t_avgt += Bi->t;
          Bi->session->positionLigand(Bi);
          Bi->t = 0.;
          Bi->t_dwell = 0.;
          Bi->t_dwell_max = 0.;
          Bi->t_dwell_total = 0.;
        }
      }

    }
  }

  cilk_sync;

  for(int il=0; il < ligands.size(); il++) {
    Body *Bi = ligands[il];
    if(Bi->bound == true) {
      *Bi->session->Nbind += 1;
      Bi->session->recordBeta();
      Bi->bound = false;
    }
    if(Bi->exited == true) {
      *Bi->session->Nexit += 1;
      Bi->session->recordBeta();
      Bi->exited = false;
    }
    if(Bi->timedout == true) {
      *Bi->session->Nexit += 1;
      Bi->session->recordBeta();
      Bi->timedout = false;
    }
  }
}


