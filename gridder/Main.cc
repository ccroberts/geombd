#include "../Main.h"
#include "Strings.h"
#include "GAFF.h"
#include "Timer.h"

bool getInputWithFlag(int argc, char **argv, char flag, string *value) {
  int  opt;                                                                                                                                                                                                     
  bool bopt = false;
  char gopt[2] = { flag, ':' };

  for(int i=1; i < argc; i++) {
    if(argv[i][0] == '-' && argv[i][1] == flag) {
      if(argv[i+1][0] != '-') {
        *value = argv[i+1];
        bopt = true;
        i++;
        break;
      }
    }
  }

  return bopt;
}


bool coordinateToGrid(double x, double y, double z, int *Gx, int *Gy, int *Gz, vertex *origin, double delta) {
  *Gx = (int)floor(((x-origin->x)/delta) + 0.5);
  *Gy = (int)floor(((y-origin->y)/delta) + 0.5);
  *Gz = (int)floor(((z-origin->z)/delta) + 0.5);
  return true;
}


void usage() {
  printf("Usage: gridder -d [AD4.1.DAT] -n NTHREADS(=max) -r [Receptor.PDBQT] -l [Ligand.PDBQT] (Optional: -p GRID_PADDING(=40A) -s GRID_SPACING(=0.375A) -w OutputFilenamePrefix)\n");
}


int main(int argc, char **argv) {
  string datfn, ligfn, recfn, fldfn, stoken, token, name_prefix;
  double T = 298.;
  double diel_rec = 78.5;
  double grid_resolution = 0.375;
  double padding = 40., padding_sqr;

  // GAFF parameter file
  if(!getInputWithFlag(argc, argv, 'd', &datfn)) { usage(); return -1; }
  // receptor pdbqt
  if(!getInputWithFlag(argc, argv, 'r', &recfn)) { usage(); return -1; }
  // ligand pdbqt
  if(!getInputWithFlag(argc, argv, 'l', &ligfn)) { usage(); return -1; }
  // number of processors/threads
  if(getInputWithFlag(argc, argv, 'n', &stoken)) {
    __cilkrts_set_param("nworkers", stoken.c_str());
  }
  // grid padding
  if(getInputWithFlag(argc, argv, 'p', &stoken)) {
    padding = stringToDouble(stoken);
  }
  padding_sqr = padding * padding;
  // grid spacing
  if(getInputWithFlag(argc, argv, 's', &stoken)) {
    grid_resolution = stringToDouble(stoken);
  }
  // grid padding
  if(getInputWithFlag(argc, argv, 'w', &name_prefix)) {
    cout << "> Output filename prefix: " << name_prefix << endl;
  }

  // Load AD parameters
  GAFFParameters *adp = new GAFFParameters(datfn);
  LigandPDBQT *lig = new LigandPDBQT(ligfn);

  cout << "> Types in ligand:";
  set<string>::iterator it;
  for(it=lig->types_set.begin(); it!=lig->types_set.end(); ++it) {
    bool found = false;
    for(int j=0; j < adp->types.size(); j++) {
      if(adp->types[j] == *it) found = true;
    }
    if(found) cout << ' ' << *it;
    else {
      cout << "! Error: Ligand atom type " << *it << " not found in the GAFF parameters." << endl;
      exit(-1);
    }
  }
  cout << endl;

  // Load receptor file
  cout << "> Loading receptor PDBQT..." << endl;
  ReceptorPDBQT *rec = new ReceptorPDBQT(recfn, adp);
  cout << "> Receptor center: " << rec->center.x << ", " << rec->center.y << ", " << rec->center.z << endl;
  cout << "> Receptor minimum: " << rec->min.x << ", " << rec->min.y << ", " << rec->min.z << endl;
  cout << "> Receptor maximum coordinates: " << rec->max.x << ", " << rec->max.y << ", " << rec->max.z << endl;
  cout << "> Types in receptor:";
  for(it=rec->types_set.begin(); it!=rec->types_set.end(); ++it) {
    cout << ' ' << *it;
  }
  cout << endl;

  // Check to make sure all receptor atom types are in our parameter data set
  for(it=rec->types_set.begin(); it!=rec->types_set.end(); ++it) {
    bool found = false;
    for(int j=0; j < adp->types.size(); j++) {
      if(adp->types[j] == *it) found = true;
    }
    if(!found) {
      cout << "! Error: Ligand atom type " << *it << " not found in the GAFF parameters." << endl;
      return EXIT_FAILURE;
    }
  }

  // Calculate grid geometries
  vertex origin = { rec->min.x - padding, rec->min.y - padding, rec->min.z - padding };
  vertex dimensions = { (rec->max.x + padding) - origin.x, (rec->max.y + padding) - origin.y, (rec->max.z + padding) - origin.z };
  int Npoints[3] = { dimensions.x / grid_resolution, dimensions.y / grid_resolution, dimensions.z / grid_resolution };

  //// Data for all grids
  // Single point calculation storage
  vector<double> u_t(lig->types_set.size());
  // Exclusion grid
  cout << "Allocating data" << endl;
  bool ***ex = (bool***)calloc(Npoints[0], sizeof(bool**));
  for(int i=0; i < Npoints[0]; i++) {
    ex[i] = (bool**)calloc(Npoints[1], sizeof(bool*));
    for(int j=0; j < Npoints[1]; j++) {
      ex[i][j] = (bool*)calloc(Npoints[2], sizeof(bool));
    }
  }

  // Dynamic storage for VDW/HB calcs
  vector<double*> data_t;
  vector<int> type_t;
  vector<ofstream*> bpm_t;
  string typefn;
  for(it=lig->types_set.begin(); it!=lig->types_set.end(); ++it) {
    double *new_data_t = (double*)calloc(Npoints[2], sizeof(double));
    if(!new_data_t) { cout << "! Error: Could not allocate memory for grid calculation." << endl; return EXIT_FAILURE; }
    data_t.push_back(new_data_t);
    type_t.push_back(adp->index_for_type(*it));
    typefn = name_prefix;
    typefn.append(*it);
    typefn.append(".bpm");
    bpm_t.push_back(new ofstream(typefn, ios::out | ios::binary));
    // write header
    bpm_t[bpm_t.size()-1]->write((char*)&origin, sizeof(vertex));
    bpm_t[bpm_t.size()-1]->write((char*)&Npoints, sizeof(int)*3);
    bpm_t[bpm_t.size()-1]->write((char*)&grid_resolution, sizeof(double));
  }

  // Exclusion grid
  cout << "> Creating exclusion grid...";
  cout.flush();
  for(int i=0; i < rec->coordinates.size(); i++) {
    vertex Rrec = rec->coordinates[i];
    double radius = rec->radii[i];
    int Grec[3];
    if(coordinateToGrid(Rrec.x, Rrec.y, Rrec.z, &Grec[0], &Grec[1], &Grec[2], &origin, grid_resolution)) {
      int searchr = ((padding + radius) / grid_resolution) + 1;
      for(int gx=(Grec[0]-searchr); gx <= Grec[0]+searchr; gx++) {
        if(gx < 0 or gx >= Npoints[0]) continue;
        for(int gy=(Grec[1]-searchr); gy <= Grec[1]+searchr; gy++) {
          if(gy < 0 or gy >= Npoints[1]) continue;
          for(int gz=(Grec[2]-searchr); gz <= Grec[2]+searchr; gz++) {
            if(gz < 0 or gz >= Npoints[2]) continue;

            double dx = (origin.x + (gx * grid_resolution)) - Rrec.x;
            double dy = (origin.y + (gy * grid_resolution)) - Rrec.y;
            double dz = (origin.z + (gz * grid_resolution)) - Rrec.z;
            if(sqrt(dx*dx + dy*dy + dz*dz) <= radius + padding) ex[gx][gy][gz] = true;
          }
        }
      }
    }
  }
  cout << " done." << endl;
  cout << "> Starting grid calculation." << endl;

  // Start a timer and start our calculations
  Timer *timer = new Timer();
  timer->start();
  for(int nx=0, nt=0; nx < Npoints[0]; nx++) {
    double X = (nx * grid_resolution) + origin.x;
    for(int ny=0; ny < Npoints[1]; ny++) {
      double Y = (ny * grid_resolution) + origin.y;
      cilk_for(int nz=0; nz < Npoints[2]; nz++) {
        double Z = (nz * grid_resolution) + origin.z;

        if(ex[nx][ny][nz] == false) continue;

        for(int at=0; at < type_t.size(); at++)
          data_t[at][nz] = 0;

        for(int i=0; i < rec->coordinates.size(); i++) {
          vertex Rrec = rec->coordinates[i];
          vertex dr;
          dr.x = X - Rrec.x; 
          dr.y = Y - Rrec.y; 
          dr.z = Z - Rrec.z; 
          double dist_sqr = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
          if(dist_sqr > padding_sqr) continue;
          double dist = sqrt(dist_sqr);
          // Iterate of each ligand atom type and calc VDW/HB
          int rec_type = rec->types[i];
          double dist_6 = dist_sqr * dist_sqr * dist_sqr;
          double dist_12 = dist_6 * dist_6;
          for(int at=0; at < type_t.size(); at++) {
            int lig_type = type_t[at];
            pair_parameter parm;
            if(lig_type > rec_type) parm = adp->lj_map[lig_type][rec_type];
            else parm = adp->lj_map[rec_type][lig_type];
            double du_t = (parm.A / dist_12) - (parm.B / dist_6);
            data_t[at][nz] += du_t;
          }
        }
      }

      //bpm_e.write((char*)data_e, sizeof(double) * Npoints[2]);
      //bpm_d.write((char*)data_d, sizeof(double) * Npoints[2]);
      for(int at=0; at < type_t.size(); at++) {
        bpm_t[at]->write((char*)data_t[at], sizeof(double) * Npoints[2]);
      }
    }

    timer->stop();
    double percent_complete = (((double)nx+1) / (double)Npoints[0]);
    cout << "\r> " << (100.*percent_complete) << " percent complete.\t\t\t" << flush;
  }

  //bpm_e.close();
  //bpm_d.close();
  for(int at=0; at < type_t.size(); at++)
    bpm_t[at]->close();

  cout << endl << "> Done." << endl;
  return EXIT_SUCCESS;
}
