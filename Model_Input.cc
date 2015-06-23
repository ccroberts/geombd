
#include "Model.h"
#include "Body.h"
#include "Strings.h"
#include "PotentialMap.h"
#include "ESPotentialMap.h"
#include "DPotentialMap.h"
#include "Timer.h"



void Model::parseFLD() {
  string bfn, line, token;
  ifstream cfd;

  cfd.open(fldfn.c_str(), ifstream::in);
  while(getline(cfd, line)) {
    if(line[0] == '#') continue;
    parseNextValue(&line, &token);

    if(token.length() > 0) {
      if(token == "variable") {
        parseNextValue(&line, &token);
        while(parseNextValue(&line, &token)) {
          if(token.substr(0,5) == "file=") {
            char type = token[token.length()-5];
            if(type == 'e') {
              cout << "* Loading electrostatic potential map...";
              cout.flush();
              Timer *t = new Timer();
              t->start();
              esmaps.push_back(new ESPotentialMap(token.substr(5, token.length()-5)));
              t->stop();
              cout << " done. ";
              t->print(&cout);
            } else {
              if(type == 'd') {
                cout << "* Loading desolvation potential map...";
                cout.flush();
                Timer *t = new Timer();
                t->start();
                esmaps.push_back(new DPotentialMap(token.substr(5, token.length()-5)));
                t->stop();
                cout << " done. ";
                t->print(&cout);
              } else {
                cout << "* Loading potential map for atom type " << type << "...";
                cout.flush();
                Timer *t = new Timer();
                t->start();
                typemaps.push_back(new TypePotentialMap(token.substr(5, token.length()-5), type));
                t->stop();
                cout << " done. ";
                t->print(&cout);
              }
            }
          }
        }
      }
    }
  }

  cfd.close();
}




void Model::parseDLG() {
  string line, token;
  Body *bi;
  Bead *bj;

  ifstream fd(dlgfn);

  while(getline(fd, line)) {
    if(parseNextValue(&line, &token)) {
      if(token == "DOCKED:") {
        parseNextValue(&line, &token);
        if(token == "MODEL") {
          bi = new Body();
        }
        if(token == "ATOM") {
          bj = new Bead();
          parseNextValue(&line, &token);
          parseNextValue(&line, &token);
          if(token == "C") {
            bj->m = 12.;
            bj->r = 1.7; 
          }
          if(token == "N") {
            bj->m = 14.;
            bj->r = 1.55; 
          }
          if(token == "O") {
            bj->m = 16.;
            bj->r = 1.52; 
          }
          if(token == "H") {
            bj->m = 1.;
            bj->r = 1.2; 
          }
          if(bj->m == 0.) {
            cout << "!!!! FATAL: Unassigned ligand bead type: " << token << endl;
            exit(-1);
          }
          parseNextValue(&line, &token);
          parseNextValue(&line, &token);
          parseNextValue(&line, &token);
          parseNextValue(&line, &token);
          bj->R.x = stringToDouble(token);
          parseNextValue(&line, &token);
          bj->R.y = stringToDouble(token);
          parseNextValue(&line, &token);
          bj->R.z = stringToDouble(token);
          parseNextValue(&line, &token);
          parseNextValue(&line, &token);
          parseNextValue(&line, &token);
          bj->q = stringToDouble(token);
          parseNextValue(&line, &token);
          bj->type = token[0];
          bi->beads.push_back(bj);
        }
        if(token == "ENDMDL") {
          cout << "* Loaded a ligand conformation (" << conformations.size() << ")" << endl;
          conformations.push_back(bi);
        }
      }
    }
  }

  // Create our ligand replicates
  cout << "* Replicating ligand x " << Nreplicates << endl;
  for(int i=0; i < Nreplicates; i++) {
    int rndConf = floor(random(0.0, (double)conformations.size()));
    Body *rc = conformations[rndConf];
    Body *bi = new Body(this);
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
    ligands.push_back(bi);
  }

  // Position our ligands
  if(ligandPosition == LIGAND_POSITION_RADIAL) {
    for(int i=0; i < ligands.size(); i++) {
      Body *lig = ligands[i];
      lig->center();
      lig->translate(random(-0.5*100, 0.5*100), random(-0.5*100, 0.5*100), random(-0.5*100, 0.5*100));

      vertex Q = lig->R;
      double l = sqrt(Q.x*Q.x + Q.y*Q.y + Q.z*Q.z);
      Q.x *= r_ligand / l;
      Q.y *= r_ligand / l;
      Q.z *= r_ligand / l;

      lig->center();
      lig->translate(bindingSite.x + Q.x, bindingSite.y + Q.y, bindingSite.z + Q.z);
      lig->rotate(random(0., M_PI), random(0., M_PI), random(0., M_PI));
    }
    r2_escape = pow(r_ligand * 5., 2.);
    cout << "* Radial mode set. Bradius: " << r_ligand <<"A -- Qradius: " << sqrt(r2_escape) << "A" << endl;
  }

  if(ligandPosition == LIGAND_POSITION_ABSOLUTE) {
    if(!(R_ligand.x == 0 and R_ligand.y == 0 and R_ligand.z == 0)) {

      for(int i=0; i < ligands.size(); i++) {
        Body *lig = ligands[i];
        lig->center();
        lig->translate(R_ligand.x, R_ligand.y, R_ligand.z);
      }

      vertex dr;
      dr.x = bindingSite.x - R_ligand.x;
      dr.y = bindingSite.y - R_ligand.y;
      dr.z = bindingSite.z - R_ligand.z;
      r_ligand = sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
    } else {
      cout << "Absolute positioning, but with no initial translation. Relying on RBBD coordinates" << endl;
    }
  }

  if(ligandPosition == LIGAND_POSITION_RANDOM) {
    for(int i=0; i < ligands.size(); i++) {
      Body *lig = ligands[i];
      lig->center();
      lig->translate(random(-0.5*bounds.x, 0.5*bounds.x), random(-0.5*bounds.y, 0.5*bounds.y), random(-0.5*bounds.z, 0.5*bounds.z));
    }
  }
}

