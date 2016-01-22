
#include "Model.h"
#include "Body.h"
#include "Strings.h"
#include "BinaryPotentialMap.h"
#include "ESPotentialMap.h"
#include "APBSPotentialMap.h"
#include "DPotentialMap.h"
#include "Timer.h"
#include "Session.h"
#include "BindingCriteria.h"



void Model::parseInputFile() {
  string lfn, line, token, type;
  ifstream cfd;

  cfd.open(ifn.c_str(), ifstream::in);
  while(getline(cfd, line)) {
    if(line[0] == '#') continue;
    parseNextValue(&line, &token);

    if(token.length() > 0) {
      if(token == "bpm") {
        parseNextValue(&line, &type);
        parseNextValue(&line, &token);
            if(type == "e") {
              cout << "* Loading electrostatic potential map...";
              cout.flush();
              Timer *t = new Timer();
              t->start();
              esmaps.push_back(new ESPotentialMap(token, type));
              t->stop();
              cout << " done. ";
              t->print(&cout);
            } else {
              if(type == "d") {
                cout << "* Loading desolvation potential map...";
                cout.flush();
                Timer *t = new Timer();
                t->start();
                esmaps.push_back(new DPotentialMap(token, type));
                t->stop();
                cout << " done. ";
                t->print(&cout);
              } else {
                cout << "* Loading potential map for atom type " << type << "...";
                cout.flush();
                Timer *t = new Timer();
                t->start();
                typemaps.push_back(new TypePotentialMap(token, type));
                t->stop();
                cout << " done. ";
                t->print(&cout);
              }
            }
      }
      if(token == "apbs") {
        parseNextValue(&line, &token);
        cout << "* Loading APBS electrostatic potential map...";
        cout.flush();
        Timer *t = new Timer();
        t->start();
        apbsmaps.push_back(new APBSPotentialMap(token, kB * T));
        t->stop();
        cout << " done. ";
        t->print(&cout);
      }
      /*
      if(token == "ex") {
        parseNextValue(&line, &token);
        cout << "* Loading exclusion map \"" << token << "\"...";
        cout.flush();
        Timer *t = new Timer();
        t->start();
        xmaps.push_back(new ExclusionMap(token, "x"));
        t->stop();
        cout << " done. ";
        t->print(&cout);
      }
      */
      if(token == "temperature") {
        parseNextValue(&line, &token);
        T = stringToDouble(token);
        cout << "* System temperature: " << T << endl;
      }
      if(token == "writetraj") {
        parseNextValue(&line, &token);
        Vtraj = stringToInt(token);
        cout << "* Writing trajectory every " << Vtraj << " steps." << endl;
      }
      if(token == "timestep") {
        parseNextValue(&line, &token);
        dt_fine = stringToDouble(token);
        parseNextValue(&line, &token);
        dt_coarse = stringToDouble(token);
        cout << "* Timesteps: fine="<<dt_fine<<" ps, coarse=" <<dt_coarse << " ps" <<endl;
      }
      if(token == "receptor") {
        parseNextValue(&line, &token);
        parseReceptorPDBQE(token);
      }
      if(token == "associate") {
        SessionRadial *sr = new SessionRadial(this);
        sessions.push_back(sr);
        sr->id = sessions.size();
        cout << "* Defining session: Association" << endl;
      }
      if(token == "transfer") {
        parseNextValue(&line, &token); //"periodic" or "radial"
        if(token == "periodic") {
          SessionAbsolutePeriodic *sap = new SessionAbsolutePeriodic(this);
          sessions.push_back(sap);
          sap->id = sessions.size();
          cout << "* Defining session: fixed-volume interenzyme transfer" << endl;
        }
        if(token == "radial") {
          SessionAbsoluteRadial *sar = new SessionAbsoluteRadial(this);
          sessions.push_back(sar);
          sar->id = sessions.size();
          cout << "* Defining session: interenzyme transfer with radial boundary" << endl;
        }
      }

      if(token == "ligand") {
        parseNextValue(&line, &token); //ligand filename
        parseLigandPDBQE(token);
        parseNextValue(&line, &token); //Nreplicates
        sessions[sessions.size()-1]->Nreplicates = stringToInt(token);
        cout<< " + Ligand replicates: " << token << endl;
      }
      if(token == "from") {
        parseNextValue(&line, &token);
        double _fx = stringToDouble(token);
        parseNextValue(&line, &token);
        double _fy = stringToDouble(token);
        parseNextValue(&line, &token);
        double _fz = stringToDouble(token);
        SessionAbsolutePeriodic *sap = dynamic_cast< SessionAbsolutePeriodic* >(sessions[sessions.size()-1]);
        SessionAbsoluteRadial *sar = dynamic_cast< SessionAbsoluteRadial* >(sessions[sessions.size()-1]);
        if(sap) {
          sap->start.x = _fx;
          sap->start.y = _fy;
          sap->start.z = _fz;
          cout << " + Ligand starting position: " << _fx << ", " << _fy << ", " << _fz << endl;
        }
        if(sar) {
          sar->start.x = _fx;
          sar->start.y = _fy;
          sar->start.z = _fz;
          cout << " + Ligand starting position: " << _fx << ", " << _fy << ", " << _fz << endl;
        }
      }
      if(token == "bind") {
        BindingCriteria *bc = new BindingCriteria();
        cout << " + Binding criteria (AND): " << endl;
        parseNextValue(&line, &token);
        double bx = stringToDouble(token);
        parseNextValue(&line, &token);
        double by = stringToDouble(token);
        parseNextValue(&line, &token);
        double bz = stringToDouble(token);
        parseNextValue(&line, &token);
        int laid = stringToInt(token);
        parseNextValue(&line, &token);
        double r = stringToDouble(token);
        bc->addPair(bx, by, bz, laid-1, r);
        cout << "    - LAID: " << laid << " within " << r << "A of (" << bx << " " << by << " " << bz << ")" << endl;
        sessions[sessions.size()-1]->bindingCriteria.push_back(bc);
      }
      if(token == "bindand") {
        BindingCriteria *bc = new BindingCriteria();
        cout << " + Binding criteria (AND): " << endl;
        while(parseNextValue(&line, &token)) {
          double bx = stringToDouble(token);
          parseNextValue(&line, &token);
          double by = stringToDouble(token);
          parseNextValue(&line, &token);
          double bz = stringToDouble(token);
          parseNextValue(&line, &token);
          int laid = stringToInt(token);
          parseNextValue(&line, &token);
          double r = stringToDouble(token);
          bc->addPair(bx, by, bz, laid-1, r);
          cout << "    - LAID: " << laid << " within " << r << "A of (" << bx << " " << by << " " << bz << ")" << endl;
        }
        sessions[sessions.size()-1]->bindingCriteria.push_back(bc);
      }
      if(token == "bindor") {
        BindingCriteria *bc = new BindingCriteria(false);
        cout << " + Binding criteria (OR): " << endl;
        while(parseNextValue(&line, &token)) {
          double bx = stringToDouble(token);
          parseNextValue(&line, &token);
          double by = stringToDouble(token);
          parseNextValue(&line, &token);
          double bz = stringToDouble(token);
          parseNextValue(&line, &token);
          int laid = stringToInt(token);
          parseNextValue(&line, &token);
          double r = stringToDouble(token);
          bc->addPair(bx, by, bz, laid-1, r);
          cout << "    - LAID: " << laid << " within " << r << "A of (" << bx << " " << by << " " << bz << ")" << endl;
        }
        sessions[sessions.size()-1]->bindingCriteria.push_back(bc);
      }
      if(token == "bounds") {
        parseNextValue(&line, &token);
        double _bx = stringToDouble(token);
        parseNextValue(&line, &token);
        double _by = stringToDouble(token);
        parseNextValue(&line, &token);
        double _bz = stringToDouble(token);
        SessionAbsolutePeriodic *sap = dynamic_cast< SessionAbsolutePeriodic* >(sessions[sessions.size()-1]);
        if(sap) {
          sap->bounds.x = _bx;
          sap->bounds.y = _by;
          sap->bounds.z = _bz;
          cout << " + Periodic boundary: " << _bx << ", " << _by << ", " << _bz << endl;
        }
      }
      if(token == "b") {
        parseNextValue(&line, &token);
        SessionRadial *sr = dynamic_cast< SessionRadial* >(sessions[sessions.size()-1]);
        if(sr) {
          sr->b = stringToDouble(token);
          cout << " + Starting radius: " << sr->b << " A" << endl;
        }
      }
      if(token == "q") {
        parseNextValue(&line, &token);
        SessionRadial *sr = dynamic_cast< SessionRadial* >(sessions[sessions.size()-1]);
        if(sr) {
          sr->q = stringToDouble(token);
          sr->q2 = sr->q * sr->q;
          cout << " + Exit radius: " << sr->q << " A" << endl;
        }
        SessionAbsoluteRadial *sar = dynamic_cast< SessionAbsoluteRadial* >(sessions[sessions.size()-1]);
        if(sar) {
          sar->q = stringToDouble(token);
          sar->q2 = sar->q * sar->q;
          cout << " + Exit radius: " << sar->q << " A" << endl;
        }
      }
      if(token == "time") {
        SessionAbsolutePeriodic *sap = dynamic_cast< SessionAbsolutePeriodic* >(sessions[sessions.size()-1]);
        if(sap) {
          parseNextValue(&line, &token);
          sap->t_max = stringToDouble(token);
          cout << " + Time limit set to " << sap->t_max << " ps" << endl;
        }
      }
      // Miletone specific
      if(token == "milestone") {
        SessionMilestone *sm = new SessionMilestone(this);
        sessions.push_back(sm);
        sm->id = sessions.size();
        cout << "* Defining session: Milestone" << endl;
        parseNextValue(&line, &token);
        sm->spacing = stringToDouble(token);
      }
      if(token == "reaction") {
        SessionMilestone *sm = dynamic_cast< SessionMilestone* >(sessions[sessions.size()-1]);
        if(sm) {
          parseNextValue(&line, &token);
          center.x = stringToDouble(token);
          parseNextValue(&line, &token);
          center.y = stringToDouble(token);
          parseNextValue(&line, &token);
          center.z = stringToDouble(token);
          parseNextValue(&line, &token);
          sm->states.push_back(new MilestoneState(stringToDouble(token) + sm->spacing*0.5));
          sm->state = sm->states[0];
        }
      }
    }
  }

  cfd.close();
}


void Model::parseReceptorPDBQE(string rfn) {
  string line, token;
  vector<double> rx;
  vector<double> ry;
  vector<double> rz;
  double cr[3] = { 0., 0., 0. };

  ifstream fd(rfn);

  while(getline(fd, line)) {
    if(starts_with(&line, "ATOM")) {
      double x = stringToDouble(line.substr(30, 8));
      double y = stringToDouble(line.substr(38, 8));
      double z = stringToDouble(line.substr(46, 8));
      rx.push_back(x);
      ry.push_back(y);
      rz.push_back(z);
      cr[0] += x;
      cr[1] += y;
      cr[2] += z;
    }
  }

  cr[0] /= rx.size();
  cr[1] /= rx.size();
  cr[2] /= rx.size();

  double rmsd = 0., dr = 0.;

  for(int i=0; i < rx.size(); i++) {
    dr = rx[i] - cr[0];
    rmsd += dr * dr;
    dr = ry[i] - cr[1];
    rmsd += dr * dr;
    dr = rz[i] - cr[2];
    rmsd += dr * dr;
  }

  rmsd /= rx.size();
  receptorRoG = sqrt(rmsd);

  center.x = cr[0];
  center.y = cr[1];
  center.z = cr[2];

  rx.clear();
  ry.clear();
  rz.clear();
}


void Model::parseLigandPDBQE(string lfn) {
  string line, token;
  Body *bi = new Body(this, sessions[sessions.size()-1]);
  Bead *bj = NULL;

  int Nconfs = 0;

  ifstream fd(lfn);

  while(getline(fd, line)) {
    if(starts_with(&line, "ATOM") or starts_with(&line, "HETATM")) {
      if(bi == NULL) {
        bi = new Body(this, sessions[sessions.size()-1]);
      }
      bj = new Bead();

      char element = line[13];
      double x = stringToDouble(line.substr(30, 8));
      double y = stringToDouble(line.substr(38, 8));
      double z = stringToDouble(line.substr(46, 8));
      double q = stringToDouble(line.substr(70, 7));
      string at = line.substr(78, 2);
      if(at == "A")  at = "C";
      if(at == "HD") at = "H";
      if(at == "HS") at = "H";
      if(at == "NA") at = "N";
      if(at == "NS") at = "N";
      if(at == "OA") at = "O";
      if(at == "OS") at = "O";
      if(at == "SA") at = "S";
      if(at == "CL") at = "Cl";
      if(at == "BR") at = "Br";
      if(at == "MG") at = "Mg";
      if(at == "CA") at = "Ca";
      if(at == "MN") at = "Mn";
      if(at == "FE") at = "Fe";
      if(at == "ZN") at = "Zn";

      if(at == "C " or at == "C") {
        bj->m = 12.;
        bj->r = 1.7; 
      }
      if(at == "N " or at == "N") {
        bj->m = 14.;
        bj->r = 1.55; 
      }
      if(at == "O " or at == "O") {
        bj->m = 16.;
        bj->r = 1.52; 
      }
      if(at == "H " or at == "H") {
        bj->m = 1.;
        bj->r = 1.2; 
      }
      if(at == "S " or at == "S") {
        bj->m = 32.;
        bj->r = 1.8; 
      }
      if(bj->m == 0.) {
        cout << "!!!! FATAL: Unassigned ligand bead type: " << at << endl;
        exit(-1);
      }
      bj->R.x = x;
      bj->R.y = y;
      bj->R.z = z;
      bj->q = q;
      bj->type = at;
      cout << "> Atom type: " << at << endl;
      bi->beads.push_back(bj);
    }
    if(starts_with(&line, "END")) {
      if(bi != NULL) {
        sessions[sessions.size()-1]->conformations.push_back(bi);
        bi->define();
        bi = NULL;
        Nconfs++;
      }
    }
  }
  cout << " + Loaded " << Nconfs << " ligand conformation" << endl;
}


/*void Model::parseLigandPDB(string lfn) {
  string line, token;
  Body *bi;
  Bead *bj;

  ifstream fd(lfn);

  while(getline(fd, line)) {
    if(starts_with(&line, "DOCKED: MODEL")) {
      bi = new Body(this, sessions[sessions.size()-1]);
    }
    if(starts_with(&line, "DOCKED: ATOM")) {
      bj = new Bead();

      char element = line[21];
      double x = stringToDouble(line.substr(38, 8));
      double y = stringToDouble(line.substr(46, 8));
      double z = stringToDouble(line.substr(54, 8));
      double q = stringToDouble(line.substr(78, 6));
      string at = line.substr(85, 2);
      if(at == "A")  at = "C";
      if(at == "HD") at = "H";
      if(at == "HS") at = "H";
      if(at == "NA") at = "N";
      if(at == "NS") at = "N";
      if(at == "OA") at = "O";
      if(at == "OS") at = "O";
      if(at == "SA") at = "S";
      if(at == "CL") at = "Cl";
      if(at == "BR") at = "Br";
      if(at == "MG") at = "Mg";
      if(at == "CA") at = "Ca";
      if(at == "MN") at = "Mn";
      if(at == "FE") at = "Fe";
      if(at == "ZN") at = "Zn";

      if(at == "C") {
        bj->m = 12.;
        bj->r = 1.7; 
      }
      if(at == "N") {
        bj->m = 14.;
        bj->r = 1.55; 
      }
      if(at == "O") {
        bj->m = 16.;
        bj->r = 1.52; 
      }
      if(at == "H") {
        bj->m = 1.;
        bj->r = 1.2; 
      }
      if(at == "S") {
        bj->m = 32.;
        bj->r = 1.8; 
      }
      if(bj->m == 0.) {
        cout << "!!!! FATAL: Unassigned ligand bead type: " << at << endl;
        exit(-1);
      }
      bj->R.x = x;
      bj->R.y = y;
      bj->R.z = z;
      bj->q = q;
      bj->type = at;
      cout << "> Atom type: " << at << endl;
      bi->beads.push_back(bj);
    }
    if(starts_with(&line, "DOCKED: ENDMDL")) {
      sessions[sessions.size()-1]->conformations.push_back(bi);
    }
  }

  for(int i=0; i < sessions[sessions.size()-1]->conformations.size(); i++) {
    sessions[sessions.size()-1]->conformations[i]->define();
  }

  cout << " + Loaded " << sessions[sessions.size()-1]->conformations.size() << " ligand conformations" << endl;
}*/




