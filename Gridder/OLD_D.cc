#include "Main.h"
#include "Strings.h"
#include "Timer.h"
#include "GBD2Parameters.h"

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
  printf("Usage: Gridder-D -d [GBD2Paramter.ff] -o [OUT.bpm] -n NTHREADS(=max) -r [Receptor.PDBQE] (Optional: -p GRID_PADDING(=40A) -s GRID_SPACING(=0.375A))\n");
}


int main(int argc, char **argv) {
  string datfn, ligfn, outfn, recfn, fldfn, stoken, token;
  double T = 298.;
  double grid_resolution = 0.375;
  double padding = 40., padding_sqr;
  double multiplier = 0.1322; //FE_coeff_desolv (AD4.1)

  // GBD2 parameter file
  if(!getInputWithFlag(argc, argv, 'd', &datfn)) { usage(); return -1; }
  // Output file
  if(!getInputWithFlag(argc, argv, 'o', &outfn)) { usage(); return -1; }
  // receptor pdbqt
  if(!getInputWithFlag(argc, argv, 'r', &recfn)) { usage(); return -1; }
  // number of processors/threads
  if(getInputWithFlag(argc, argv, 'n', &stoken)) {
    __cilkrts_set_param("nworkers", stoken.c_str());
  }
  // coefficient
  if(getInputWithFlag(argc, argv, 'm', &stoken)) {
    multiplier = stringToDouble(stoken);
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

  // Load GBD2 parameters
  GBD2Parameters *gbdparm = new GBD2Parameters(datfn);
  // Load receptor file
  cout << "> Loading receptor PDBQE..." << endl;
  ReceptorPDBQE *rec = new ReceptorPDBQE(recfn, gbdparm);
  cout << "> Receptor center: " << rec->center.x << ", " << rec->center.y << ", " << rec->center.z << endl;
  cout << "> Receptor minimum: " << rec->min.x << ", " << rec->min.y << ", " << rec->min.z << endl;
  cout << "> Receptor maximum coordinates: " << rec->max.x << ", " << rec->max.y << ", " << rec->max.z << endl;

  // Calculate grid geometries
  vertex origin = { rec->min.x - padding, rec->min.y - padding, rec->min.z - padding };
  vertex dimensions = { (rec->max.x + padding) - origin.x, (rec->max.y + padding) - origin.y, (rec->max.z + padding) - origin.z };
  int Npoints[3] = { dimensions.x / grid_resolution, dimensions.y / grid_resolution, dimensions.z / grid_resolution };
  int Ntotal = Npoints[0] * Npoints[1] * Npoints[2];

  //// Data for all grids
  // Z vector storage for E calcs
  double *data_e = (double*)calloc(Npoints[2], sizeof(double));
  if(!data_e) { cout << "! Error: Could not allocate memory for grid calculation." << endl; return EXIT_FAILURE; }
  ofstream bpm_e(outfn, ios::out | ios::binary);
  bpm_e.write((char*)&origin, sizeof(vertex));
  bpm_e.write((char*)&Npoints, sizeof(int)*3);
  bpm_e.write((char*)&grid_resolution, sizeof(double));

  // Exclusion grid
  cout << "Allocating data" << endl;
  bool ***ex = (bool***)calloc(Npoints[0], sizeof(bool**));
  for(int i=0; i < Npoints[0]; i++) {
    ex[i] = (bool**)calloc(Npoints[1], sizeof(bool*));
    for(int j=0; j < Npoints[1]; j++) {
      ex[i][j] = (bool*)calloc(Npoints[2], sizeof(bool));
    }
  }
  cout << "> Creating exclusion grid...";
  cout.flush();
  cilk_for(int i=0; i < rec->coordinates.size(); i++) {
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
            double rr = sqrt(dx*dx + dy*dy + dz*dz);
            if(rr <= radius + padding) ex[gx][gy][gz] = true;
          }
        }
      }
    }
  }

  // relevant atom list
  vector<int> atoms;

  // Start a timer and start our calculations
  Timer *timer = new Timer();
  timer->start();
  for(int nx=0, nt=0; nx < Npoints[0]; nx++) {
    double X = (nx * grid_resolution) + origin.x;

    atoms.clear();
    for(int i=0; i < rec->coordinates.size(); i++) {
      vertex Rrec = rec->coordinates[i];
      if(Rrec.x <= X + padding and Rrec.x >= X - padding)
        atoms.push_back(i);
    }

    for(int ny=0; ny < Npoints[1]; ny++) {
      double Y = (ny * grid_resolution) + origin.y;
      cilk_for(int nz=0; nz < Npoints[2]; nz++) {
        double Z = (nz * grid_resolution) + origin.z;

        if(ex[nx][ny][nz] == false) continue;

        data_e[nz] = 0;

        for(int i=0; i < atoms.size(); i++) {
          if(rec->charges[atoms[i]] == 0.) continue;
          vertex Rrec = rec->coordinates[atoms[i]];
          vertex dr;
          dr.x = X - Rrec.x; 
          dr.y = Y - Rrec.y; 
          dr.z = Z - Rrec.z; 
          double dist_sqr = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
          if(dist_sqr > padding_sqr) continue;
          double dist = sqrt(dist_sqr);
          if(dist < 0.5*rec->radii[atoms[i]]) dist = 0.5*rec->radii[atoms[i]];
          // Desolvation
          double du_d = gbdparm->vol[rec->types[atoms[i]]] * exp(-dist_sqr/24.5);/*24.5 = 2*s^2 where s = 3.5*/
          data_e[nz] += du_d;
        }
        data_e[nz] *= -0.01097; //Qasp = 0.010 97 stderr=0.002 63 (AD4.1)
        data_e[nz] *=  multiplier;
      }

      bpm_e.write((char*)data_e, sizeof(double) * Npoints[2]);
    }

    timer->stop();
    double percent_complete = (((double)nx+1) / (double)Npoints[0]);
    cout << "\r> " << (100.*percent_complete) << " percent complete.                      " << flush;
  }

  bpm_e.close();

  cout << endl << "> Done." << endl;
  return EXIT_SUCCESS;
}
