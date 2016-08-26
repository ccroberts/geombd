
#include "Model.h"
#include "Body.h"
#include "Strings.h"
#include "Grid.h"
#include "Grid_ES.h"
#include "Grid_EX.h"
#include "Grid_D.h"
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
      if(token == "grid") {
        parseNextValue(&line, &type);
        parseNextValue(&line, &token);
            if(type == "es") {
              lout << "* Loading electrostatic potential grid \"" << token << "\"...";
              lout.flush();
              Timer *t = new Timer();
              t->start();
              esmaps.push_back(new Grid_ES(token, type));
              t->stop();
              lout << " done. ";
              t->print(&lout);
            } else {
              if(type == "ex") {
                lout << "* Loading exclusion grid \"" << token << "\"...";
                lout.flush();
                Timer *t = new Timer();
                t->start();
                exmaps.push_back(new Grid_EX(token, "x"));
                t->stop();
                lout << " done. ";
                t->print(&lout);
              } else {
                if(type == "d") {
                  lout << "* Loading desolvation potential grid \"" << token << "\"...";
                  lout.flush();
                  Timer *t = new Timer();
                  t->start();
                  dmaps.push_back(new Grid_D(token, type));
                  t->stop();
                  lout << " done. ";
                  t->print(&lout);
                } else {
                  lout << "* Loading LJ potential energy grid \"" << token << "\" for atom type \"" << type << "\"...";
                  lout.flush();
                  Timer *t = new Timer();
                  t->start();
                  typemaps.push_back(new Grid_Type(token, type));
                  t->stop();
                  lout << " done. ";
                  t->print(&lout);
                }
              }
            }
      }
      if(token == "debug") {
        parseNextValue(&line, &token);
        lout << "* Loading debug map \"" << token << "\"" << endl;
        debug_map = new Grid_EX(token, "x");
      }
      if(token == "convergence") {
        parseNextValue(&line, &token);
        lout << "* Convergence criteria: " << token << endl;
        convergence = stringToDouble(token);
      }
      if(token == "calcconv") {
        parseNextValue(&line, &token);
        rate_conv = stringToInt(token);
        lout << "* Checking convergence every " << rate_conv<< " steps." << endl;
      }
      if(token == "threads") {
        parseNextValue(&line, &token);
        __cilkrts_set_param("nworkers", token.c_str());
        lout << "* Attempting to set number of threads to " << token << endl;
      }
      if(token == "temperature") {
        parseNextValue(&line, &token);
        T = stringToDouble(token);
        lout << "* System temperature: " << T << endl;
      }
      if(token == "writetraj") {
        parseNextValue(&line, &token);
        rate_trj = stringToInt(token);
        lout << "* Writing trajectory every " << rate_trj << " steps." << endl;
      }
      if(token == "writelog") {
        parseNextValue(&line, &token);
        rate_log = stringToInt(token);
        lout << "* Writing association rate information to logfile every " << rate_log << " steps." << endl;
      }
      if(token == "timestep") {
        parseNextValue(&line, &token);
        dt_fine = stringToDouble(token);
        parseNextValue(&line, &token);
        dt_coarse = stringToDouble(token);
        parseNextValue(&line, &token);
        dt_scale_start = pow(stringToDouble(token), 2.); //SQUARED VALUE!
        parseNextValue(&line, &token);
        dt_scale_end = pow(stringToDouble(token), 2.);
        lout << "* Timesteps: fine="<<dt_fine<<" ps, coarse=" <<dt_coarse << " ps, scaling from a radius of " << sqrt(dt_scale_start) << "A to a radius of " << sqrt(dt_scale_end) << "A." <<endl;
      }
      if(token == "order") {
        parseNextValue(&line, &token);
        fd_order = stringToInt(token);
        lout << "* Finite difference force approximation using order " << fd_order << "." << endl;
      }
      if(token == "receptor") {
        parseNextValue(&line, &token);
        lout<< "* Receptor filename: " << token << endl;
        if(file_exists(token)) {
          parseReceptorPQR(token);
        } else {
          cout <<"! Receptor file does not exist! Exiting." << endl;
          exit(-1);
        }
      }
      if(token == "associate") {
        SessionRadial *sr = new SessionRadial(this);
        sessions.push_back(sr);
        sr->id = sessions.size();
        lout << "* Defining session: Association" << endl;
      }
      if(token == "transfer") {
        parseNextValue(&line, &token); //"periodic" or "radial"
        if(token == "periodic") {
          SessionAbsolutePeriodic *sap = new SessionAbsolutePeriodic(this);
          sessions.push_back(sap);
          sap->id = sessions.size();
          lout << "* Defining session: fixed-volume interenzyme transfer" << endl;
        }
        if(token == "radial") {
          SessionAbsoluteRadial *sar = new SessionAbsoluteRadial(this);
          sessions.push_back(sar);
          sar->id = sessions.size();
          lout << "* Defining session: interenzyme transfer with radial boundary" << endl;
        }
      }

      if(token == "ligand") {
        parseNextValue(&line, &token); //ligand filename
        lout<< "* Ligand filename: " << token << endl;
        if(file_exists(token)) {
          parseLigandPQR(token);
        } else {
          cout <<"! Ligand file does not exist! Exiting." << endl;
          exit(-1);
        }
        parseNextValue(&line, &token); //Nreplicates
        sessions[sessions.size()-1]->Nreplicates = stringToInt(token);
        lout<< "* Ligand replicates: " << token << endl;
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
          lout << "* Ligand starting position: " << _fx << ", " << _fy << ", " << _fz << endl;
        }
        if(sar) {
          sar->start.x = _fx;
          sar->start.y = _fy;
          sar->start.z = _fz;
          lout << "* Ligand starting position: " << _fx << ", " << _fy << ", " << _fz << endl;
        }
      }
      if(token == "bind") {
        BindingCriteria *bc = new BindingCriteria();
        lout << "* Binding criteria (Single Ligand Atom): " << endl;
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
        lout << "    - LAID: " << laid << " within " << r << "A of (" << bx << " " << by << " " << bz << ")" << endl;
        sessions[sessions.size()-1]->bindingCriteria.push_back(bc);
      }
      if(token == "bindand") {
        BindingCriteria *bc = new BindingCriteria();
        lout << "* Binding criteria (AND): " << endl;
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
          lout << "    - LAID: " << laid << " within " << r << "A of (" << bx << " " << by << " " << bz << ")" << endl;
        }
        sessions[sessions.size()-1]->bindingCriteria.push_back(bc);
      }
      if(token == "bindor") {
        BindingCriteria *bc = new BindingCriteria(false);
        lout << " + Binding criteria (OR): " << endl;
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
          lout << "    - LAID: " << laid << " within " << r << "A of (" << bx << " " << by << " " << bz << ")" << endl;
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
          lout << " + Periodic boundary: " << _bx << ", " << _by << ", " << _bz << endl;
        }
      }
      if(token == "b") {
        parseNextValue(&line, &token);
        SessionRadial *sr = dynamic_cast< SessionRadial* >(sessions[sessions.size()-1]);
        if(sr) {
          sr->b = stringToDouble(token);
          lout << " + Starting radius: " << sr->b << " A" << endl;
        }
      }
      if(token == "q") {
        parseNextValue(&line, &token);
        SessionRadial *sr = dynamic_cast< SessionRadial* >(sessions[sessions.size()-1]);
        if(sr) {
          sr->q = stringToDouble(token);
          sr->q2 = sr->q * sr->q;
          lout << " + Exit radius: " << sr->q << " A" << endl;
        }
        SessionAbsoluteRadial *sar = dynamic_cast< SessionAbsoluteRadial* >(sessions[sessions.size()-1]);
        if(sar) {
          sar->q = stringToDouble(token);
          sar->q2 = sar->q * sar->q;
          lout << " + Exit radius: " << sar->q << " A" << endl;
        }
      }
      if(token == "time") {
        SessionAbsolutePeriodic *sap = dynamic_cast< SessionAbsolutePeriodic* >(sessions[sessions.size()-1]);
        if(sap) {
          parseNextValue(&line, &token);
          sap->t_max = stringToDouble(token);
          lout << " + Time limit set to " << sap->t_max << " ps" << endl;
        }
      }
      if(token == "maxsims") {
        parseNextValue(&line, &token);
        max_simulations = stringToInt(token);
        lout << " + Setting maximum number of completed replicate simulations to " << max_simulations << "." << endl;
      }
    }
  }

  cfd.close();

  // Determine system geometry
  for(int i=0; i < esmaps.size(); i++) {
    if(esmaps[i]->origin[0] < bounds_min.x) bounds_min.x = esmaps[i]->origin[0];
    if(esmaps[i]->origin[1] < bounds_min.y) bounds_min.y = esmaps[i]->origin[1];
    if(esmaps[i]->origin[2] < bounds_min.z) bounds_min.z = esmaps[i]->origin[2];
    if(esmaps[i]->origin[0] > bounds_max.x) bounds_max.x = esmaps[i]->origin[0];
    if(esmaps[i]->origin[1] > bounds_max.y) bounds_max.y = esmaps[i]->origin[1];
    if(esmaps[i]->origin[2] > bounds_max.z) bounds_max.z = esmaps[i]->origin[2];
  }
  for(int i=0; i < typemaps.size(); i++) {
    if(typemaps[i]->origin[0] < bounds_min.x) bounds_min.x = typemaps[i]->origin[0];
    if(typemaps[i]->origin[1] < bounds_min.y) bounds_min.y = typemaps[i]->origin[1];
    if(typemaps[i]->origin[2] < bounds_min.z) bounds_min.z = typemaps[i]->origin[2];
    if(typemaps[i]->origin[0] > bounds_max.x) bounds_max.x = typemaps[i]->origin[0];
    if(typemaps[i]->origin[1] > bounds_max.y) bounds_max.y = typemaps[i]->origin[1];
    if(typemaps[i]->origin[2] > bounds_max.z) bounds_max.z = typemaps[i]->origin[2];
  }
  for(int i=0; i < exmaps.size(); i++) {
    if(exmaps[i]->origin[0] < bounds_min.x) bounds_min.x = exmaps[i]->origin[0];
    if(exmaps[i]->origin[1] < bounds_min.y) bounds_min.y = exmaps[i]->origin[1];
    if(exmaps[i]->origin[2] < bounds_min.z) bounds_min.z = exmaps[i]->origin[2];
    if(exmaps[i]->origin[0] > bounds_max.x) bounds_max.x = exmaps[i]->origin[0];
    if(exmaps[i]->origin[1] > bounds_max.y) bounds_max.y = exmaps[i]->origin[1];
    if(exmaps[i]->origin[2] > bounds_max.z) bounds_max.z = exmaps[i]->origin[2];
  }

  double dr[3], rmin, rmax;
  dr[0] = bounds_min.x - center.x;
  dr[1] = bounds_min.y - center.y;
  dr[2] = bounds_min.z - center.z;
  rmin = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);

  dr[0] = bounds_max.x - center.x;
  dr[1] = bounds_max.y - center.y;
  dr[2] = bounds_max.z - center.z;
  rmax = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);

  system_extent = max(rmin, rmax);
  lout << "* System extends to a maximum distance of " << system_extent << "A" << endl;
}


void Model::parseReceptorPQR(string rfn) {
  string line, token;
  vector<double> rx;
  vector<double> ry;
  vector<double> rz;
  vector<double> radii;
  double cr[3] = { 0., 0., 0. };

  ifstream fd(rfn);

  while(getline(fd, line)) {
    if(starts_with(&line, "ATOM")) {
      string at_raw = line.substr(12, 4);
      string at = trim(at_raw);
      double rd = stringToDouble(line.substr(69, 6));
      radii.push_back(rd);
      double x = stringToDouble(line.substr(30, 10));
      double y = stringToDouble(line.substr(40, 10));
      double z = stringToDouble(line.substr(50, 10));
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

  if(rx.size() == 1) {
    receptor_radius = radii[0];
  } else {
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
    receptor_radius = sqrt(rmsd);
    /*
    double invSumRij = 0.;
    //hydrodynamic radius
    for(int i=0; i < rx.size(); i++) {
      for(int j=i+1; j < rx.size(); j++) {
        dr.x = rx[i] - rx[j];
        dr.y = ry[i] - ry[j];
        dr.z = rz[i] - rz[j];
        invSumRij += 1. / sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
      }
    }

    receptor_radius = 1.0 / (invSumRij / (rx.size()*rx.size()));
    */
  }


  center.x = cr[0];
  center.y = cr[1];
  center.z = cr[2];

  rx.clear();
  ry.clear();
  rz.clear();
}


void Model::parseLigandPQR(string lfn) {
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
      double x = stringToDouble(line.substr(30, 10));
      double y = stringToDouble(line.substr(40, 10));
      double z = stringToDouble(line.substr(50, 10));
      double q = stringToDouble(line.substr(70, 7));
      string at_raw = line.substr(12, 4);
      string at = trim(at_raw);
      bj->r = stringToDouble(line.substr(69, 6));

      //TODO: Complete table
      if(at == "C") bj->m = 12.;
      if(at == "N") bj->m = 14.;
      if(at == "O") bj->m = 16.;
      if(at == "H") bj->m = 1.;
      if(at == "S") bj->m = 32.;
      if(at == "Na") bj->m = 22.9898;
      if(at == "Cl") bj->m = 35.45;

      if(bj->m == 0.) {
        lout << "! FATAL: No mass for ligand bead type: " << at << endl;
        exit(-1);
      }
      bj->R.x = x;
      bj->R.y = y;
      bj->R.z = z;
      bj->q = q;
      bj->type = at;
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
  lout << " + Loaded " << Nconfs << " ligand conformation" << endl;
}






