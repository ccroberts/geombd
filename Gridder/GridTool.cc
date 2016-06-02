#include "Main.h"
#include "Strings.h"
#include "GBD2Parameters.h"
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
  printf("Usage: GridTool -t [LJ=0, sLJ=1] -d [FORECFIELD.ff] -n NTHREADS(=max) -r [Receptor.PDBQE] -l [Ligand.PDBQE] (Optional: -p GRID_PADDING(=40A) -s GRID_SPACING(=0.375A) -w OutputFilenamePrefix)\n");
}


int main(int argc, char **argv) {
  string p_grid_type_s,
         p_param_fn,
         p_lig_fn,
         p_rec_fn, 
         p_name_prefix;
  string token;
  int p_grid_type = -1;
  double p_T = 298.;
  double p_diel_rec = 78.5;
  double p_grid_resolution = 0.375;
  double p_padding = 40., p_padding_sqr;
  Timer *timer = new Timer();
  const double pow_two_16 = pow(2., 1./6.);

  // Calculation type
  if(!getInputWithFlag(argc, argv, 't', &p_grid_type_s)) { usage(); return -1; }
  p_grid_type = stringToInt(p_grid_type_s);
  if(p_grid_type < 0 or p_grid_type > 1) {
    cout << "! Error: GridTool flag -t requires a value of 0 (LJ), or 1 (shifted LJ). " << endl;
  }
  // GBD2 parameter file
  if(!getInputWithFlag(argc, argv, 'd', &p_param_fn)) { usage(); return -1; }
  // Receptor pdbqe
  if(!getInputWithFlag(argc, argv, 'r', &p_rec_fn)) { usage(); return -1; }
  // Ligand pdbqe
  if(!getInputWithFlag(argc, argv, 'l', &p_lig_fn)) { usage(); return -1; }
  // Number of processors/threads
  if(getInputWithFlag(argc, argv, 'n', &token)) {
    __cilkrts_set_param("nworkers", token.c_str());
  }
  // Grid padding
  if(getInputWithFlag(argc, argv, 'p', &token)) {
    p_padding = stringToDouble(token);
  }
  p_padding_sqr = pow(p_padding, 2);
  // Grid spacing
  if(getInputWithFlag(argc, argv, 's', &token)) {
    p_grid_resolution = stringToDouble(token);
  }
  // Grid output filename prefix
  if(getInputWithFlag(argc, argv, 'w', &p_name_prefix)) {
    cout << "> Output filename prefix: " << p_name_prefix << endl;
  }

  // Load GeomBD2 parameters
  GBD2Parameters *gbdparm = new GBD2Parameters(p_param_fn);

  // Load ligand file
  cout << "> Loading ligand PDBQE..." << endl;
  LigandPDBQE *lig = new LigandPDBQE(p_lig_fn);
  if(lig->types_set.size() == 0) {
    cout << "! Error: No ligand atom types found." << endl;
    exit(-1);
  }
  cout << "> Types in ligand:";
  set<string>::iterator it;
  for(it=lig->types_set.begin(); it!=lig->types_set.end(); ++it) {
    bool found = false;
    for(int j=0; j < gbdparm->types.size(); j++) {
      if(gbdparm->types[j] == *it) found = true;
    }
    if(found) cout << ' ' << *it;
    else {
      cout << "! Error: Ligand atom type " << *it << " not found in the GBD2 parameters." << endl;
      exit(-1);
    }
  }
  cout << endl;

  // Load receptor file
  cout << "> Loading receptor PDBQE..." << endl;
  ReceptorPDBQE *rec = new ReceptorPDBQE(p_rec_fn, gbdparm);
  cout << "> Types in receptor:";
  for(it=rec->types_set.begin(); it!=rec->types_set.end(); ++it) {
    cout << ' ' << *it;
  }
  cout << endl;

  // Check to make sure all receptor atom types are in our parameter data set
  for(it=rec->types_set.begin(); it!=rec->types_set.end(); ++it) {
    bool found = false;
    for(int j=0; j < gbdparm->types.size(); j++) {
      if(gbdparm->types[j] == *it) found = true;
    }
    if(!found) {
      cout << "! Error: Ligand atom type " << *it << " not found in the GBD2 parameters." << endl;
      return EXIT_FAILURE;
    }
  }

  // Calculate grid geometries
  vertex origin = { rec->min.x - p_padding, rec->min.y - p_padding, rec->min.z - p_padding };
  vertex dimensions = { (rec->max.x + p_padding) - origin.x, (rec->max.y + p_padding) - origin.y, (rec->max.z + p_padding) - origin.z };
  int Npoints[3] = { dimensions.x / p_grid_resolution, dimensions.y / p_grid_resolution, dimensions.z / p_grid_resolution };
  int Ntotal = Npoints[0] * Npoints[1] * Npoints[2];

  // Calculation storage
  vector<int> relevant_atoms; // Dynamic relevant atom list
  vector<double> u_t(lig->types_set.size()); //Parallel calculation
  vector<double*> data_t;
  vector<int> type_t;
  vector<ofstream*> bpm_t;
  string typefn;

  if(p_grid_type == 0 or p_grid_type == 1) {
    for(it=lig->types_set.begin(); it!=lig->types_set.end(); ++it) {
      double *new_data_t = (double*)calloc(Npoints[2], sizeof(double));
      if(!new_data_t) { cout << "! Error: Could not allocate memory for grid calculation." << endl; return EXIT_FAILURE; }
      data_t.push_back(new_data_t);
      type_t.push_back(gbdparm->index_for_type(*it));
      typefn = p_name_prefix;
      typefn.append(*it);
      typefn.append(".bpm");
      bpm_t.push_back(new ofstream(typefn, ios::out | ios::binary));
      // write header
      bpm_t[bpm_t.size()-1]->write((char*)&origin, sizeof(vertex));
      bpm_t[bpm_t.size()-1]->write((char*)&Npoints, sizeof(int)*3);
      bpm_t[bpm_t.size()-1]->write((char*)&p_grid_resolution, sizeof(double));
    }
  }

  // Exclusion grid
  cout << "> Allocating exclusion grid memory..." << endl;
  bool ***ex = (bool***)calloc(Npoints[0], sizeof(bool**));
  for(int i=0; i < Npoints[0]; i++) {
    ex[i] = (bool**)calloc(Npoints[1], sizeof(bool*));
    for(int j=0; j < Npoints[1]; j++) {
      ex[i][j] = (bool*)calloc(Npoints[2], sizeof(bool));
    }
  }

  timer->start();
  cout << "> Creating exclusion grid...";
  cout.flush();
  cilk_for(int i=0; i < rec->coordinates.size(); i++) {
    vertex Rrec = rec->coordinates[i];
    double radius = rec->radii[i];
    int Grec[3];
    if(coordinateToGrid(Rrec.x, Rrec.y, Rrec.z, &Grec[0], &Grec[1], &Grec[2], &origin, p_grid_resolution)) {
      int searchr = ((p_padding + radius) / p_grid_resolution) + 1;
      for(int gx=(Grec[0]-searchr); gx <= Grec[0]+searchr; gx++) {
        if(gx < 0 or gx >= Npoints[0]) continue;
        for(int gy=(Grec[1]-searchr); gy <= Grec[1]+searchr; gy++) {
          if(gy < 0 or gy >= Npoints[1]) continue;
          for(int gz=(Grec[2]-searchr); gz <= Grec[2]+searchr; gz++) {
            if(gz < 0 or gz >= Npoints[2]) continue;

            double dx = (origin.x + (gx * p_grid_resolution)) - Rrec.x;
            double dy = (origin.y + (gy * p_grid_resolution)) - Rrec.y;
            double dz = (origin.z + (gz * p_grid_resolution)) - Rrec.z;
            if(sqrt(dx*dx + dy*dy + dz*dz) <= radius + p_padding) ex[gx][gy][gz] = true;
          }
        }
      }
    }
  }
  timer->stop();
  cout << endl << "> Elapsed time (s): ";
  timer->print(&cout);


  // Start a timer and start our calculations
  cout << "> Starting grid calculation." << p_grid_type << endl;
  timer->start();
  for(int nx=0, nt=0; nx < Npoints[0]; nx++) {
    double X = (nx * p_grid_resolution) + origin.x;

    relevant_atoms.clear();
    for(int i=0; i < rec->coordinates.size(); i++) {
      vertex Rrec = rec->coordinates[i];
      if(Rrec.x <= X + p_padding and Rrec.x >= X - p_padding)
        relevant_atoms.push_back(i);
    }

    for(int ny=0; ny < Npoints[1]; ny++) {
      double Y = (ny * p_grid_resolution) + origin.y;
      cilk_for(int nz=0; nz < Npoints[2]; nz++) {
        double Z = (nz * p_grid_resolution) + origin.z;

        if(ex[nx][ny][nz] == false) continue;

        for(int at=0; at < type_t.size(); at++)
          data_t[at][nz] = 0;

        for(int i=0; i < relevant_atoms.size(); i++) {
          vertex Rrec = rec->coordinates[relevant_atoms[i]];
          vertex dr;
          dr.x = X - Rrec.x; 
          dr.y = Y - Rrec.y; 
          if(dr.y > p_padding) continue;
          dr.z = Z - Rrec.z; 
          if(dr.z > p_padding) continue;
          double dist_sqr = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
          if(dist_sqr > p_padding_sqr) continue;

          if(p_grid_type == 0 /*LJ*/ or p_grid_type == 1 /*sLJ*/) {
            // Iterate of each ligand atom type and calc VDW/HB
            int rec_type = rec->types[relevant_atoms[i]];
            double dist = sqrt(dist_sqr);
            double dist_6 = dist_sqr * dist_sqr * dist_sqr;
            double dist_12 = dist_6 * dist_6;
            for(int at=0; at < type_t.size(); at++) {
              int lig_type = type_t[at];
              pair_parameter parm;
              double du_t;
              if(p_grid_type == 0) {
                if(lig_type > rec_type) parm = gbdparm->lj_map[lig_type][rec_type];
                else parm = gbdparm->lj_map[rec_type][lig_type];
                du_t = (parm.A / dist_12) - (parm.B / dist_6);
              }
              if(p_grid_type == 1) {
                if(lig_type > rec_type) parm = gbdparm->slj_map[lig_type][rec_type];
                else parm = gbdparm->slj_map[rec_type][lig_type];
                double s_over_pow_two_16 = parm.s / pow_two_16;
                du_t = (parm.A / dist_12) - (parm.B / dist_6);
              }
              data_t[at][nz] += du_t;
            }
          }
        }
      }

      //bpm_e.write((char*)data_e, sizeof(double) * Npoints[2]);
      //bpm_d.write((char*)data_d, sizeof(double) * Npoints[2]);
      for(int at=0; at < type_t.size(); at++) {
        bpm_t[at]->write((char*)data_t[at], sizeof(double) * Npoints[2]);
      }
      double percent_complete = ( (((double)nx) * Npoints[1] * Npoints[2]) + (((double)ny) * Npoints[2]) ) / (double)Ntotal;
      cout << "\r> " << (100.*percent_complete) << " percent complete.                                            " << flush;
    }

  }
  timer->stop();
  cout << endl << "> Elapsed time (s): ";
  timer->print(&cout);

  for(int at=0; at < type_t.size(); at++)
    bpm_t[at]->close();

  cout << endl << "> Done." << endl;
  return EXIT_SUCCESS;
}
