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


void usage() {
  printf("Usage: Gridder-EX -o [OUT.spm] -d [FORCEFIELD.ff] -n NTHREADS(=max) -r [Receptor.PQR] (Optional: -p GRID_PADDING(=40A) -s GRID_SPACING(=0.375A) -m RECEPTOR_RADIUS_SCALING(=1.25) -w OutputFilenamePrefix))\n");
}


bool coordinateToGrid(double x, double y, double z, int *Gx, int *Gy, int *Gz, vertex *origin, double delta) {
  *Gx = (int)floor(((x-origin->x)/delta) + 0.5);
  *Gy = (int)floor(((y-origin->y)/delta) + 0.5);
  *Gz = (int)floor(((z-origin->z)/delta) + 0.5);
  return true;
}


int main(int argc, char **argv) {
  string datfn, outfn, ligfn, recfn, fldfn, stoken, token, name_prefix;
  double grid_resolution = 0.375;
  double padding = 40., padding_sqr;
  double scaling = 1.25;

  // GBD2 parameter file
  if(!getInputWithFlag(argc, argv, 'd', &datfn)) { usage(); return -1; }
  // out file
  if(!getInputWithFlag(argc, argv, 'o', &outfn)) { usage(); return -1; }
  // receptor pqr
  if(!getInputWithFlag(argc, argv, 'r', &recfn)) { usage(); return -1; }
  // number of processors/threads
  if(getInputWithFlag(argc, argv, 'n', &stoken)) {
    __cilkrts_set_param("nworkers", stoken.c_str());
  }
  // grid padding
  if(getInputWithFlag(argc, argv, 'p', &stoken)) {
    padding = stringToDouble(stoken) + 1.4;
  }
  padding_sqr = padding * padding;
  // grid spacing
  if(getInputWithFlag(argc, argv, 's', &stoken)) {
    grid_resolution = stringToDouble(stoken);
  }
  // receptor radius scaling
  if(getInputWithFlag(argc, argv, 'm', &stoken)) {
    scaling = stringToDouble(stoken);
  }
  // filename prefix
  if(getInputWithFlag(argc, argv, 'w', &name_prefix)) {
    cout << "> Output filename prefix: " << name_prefix << endl;
  }

  // Load AD parameters
  GBD2Parameters *adp = new GBD2Parameters(datfn);
  // Load receptor file
  cout << "> Loading receptor PDBQE..." << endl;
  ReceptorPDBQE *rec = new ReceptorPDBQE(recfn, adp);
  cout << "> Receptor center: " << rec->center.x << ", " << rec->center.y << ", " << rec->center.z << endl;
  cout << "> Receptor minimum: " << rec->min.x << ", " << rec->min.y << ", " << rec->min.z << endl;
  cout << "> Receptor maximum coordinates: " << rec->max.x << ", " << rec->max.y << ", " << rec->max.z << endl;

  // Calculate grid geometries
  vertex origin = { rec->min.x - padding, rec->min.y - padding, rec->min.z - padding };
  vertex dimensions = { (rec->max.x + padding) - origin.x, (rec->max.y + padding) - origin.y, (rec->max.z + padding) - origin.z };
  int Npoints[3] = { dimensions.x / grid_resolution, dimensions.y / grid_resolution, dimensions.z / grid_resolution };

  cout << "> Npoints: " << Npoints[0] << " " << Npoints[1] << " " << Npoints[2] << endl;

  //// Data for all grids
  cout << "> Allocating data" << endl;
  short ***data = (short***)calloc(Npoints[0], sizeof(short**));
  for(int i=0; i < Npoints[0]; i++) {
    data[i] = (short**)calloc(Npoints[1], sizeof(short*));
    for(int j=0; j < Npoints[1]; j++) {
      data[i][j] = (short*)calloc(Npoints[2], sizeof(short));
    }
  }

  string ofn = name_prefix;
  ofn.append(outfn);
  ofstream fo(ofn);

  // Start a timer and start our calculations
  Timer *timer = new Timer();
  cout << "> Starting calculation..." << endl;
  timer->start();
  for(int i=0; i < rec->coordinates.size(); i++) {
    vertex Rrec = rec->coordinates[i];
    double radius = rec->radii[i];
    int Grec[3];
    if(coordinateToGrid(Rrec.x, Rrec.y, Rrec.z, &Grec[0], &Grec[1], &Grec[2], &origin, grid_resolution)) {
      int searchr = ((1.4 + radius) / grid_resolution) + 1;
      for(int gx=(Grec[0]-searchr); gx <= Grec[0]+searchr; gx++) {
        if(gx < 0 or gx >= Npoints[0]) continue;
        for(int gy=(Grec[1]-searchr); gy <= Grec[1]+searchr; gy++) {
          if(gy < 0 or gy >= Npoints[1]) continue;
          for(int gz=(Grec[2]-searchr); gz <= Grec[2]+searchr; gz++) {
            if(gz < 0 or gz >= Npoints[2]) continue;

            double dx = (origin.x + (gx * grid_resolution)) - Rrec.x;
            double dy = (origin.y + (gy * grid_resolution)) - Rrec.y;
            double dz = (origin.z + (gz * grid_resolution)) - Rrec.z;
            if(sqrt(dx*dx + dy*dy + dz*dz) <= radius * scaling) data[gx][gy][gz] = 1;
          }
        }
      }
    }
  }

  cout << "Writing spm..." << endl;
  fo.write((char*)&origin, sizeof(vertex));
  fo.write((char*)&Npoints, sizeof(int)*3);
  fo.write((char*)&grid_resolution, sizeof(double));

  for(int i=0; i < Npoints[0]; i++) {
    for(int j=0; j < Npoints[1]; j++) {
      fo.write((char*)data[i][j], sizeof(short) * Npoints[2]);
    }
  }

  timer->stop();
  cout << endl << "> Done." << endl;
  timer->print(&cout);
  return EXIT_SUCCESS;
}
