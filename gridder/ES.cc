#include "../Main.h"
#include "Strings.h"
#include "Timer.h"
#include "AutoDock4.1.h"

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
  printf("Usage: griddes -n NTHREADS(=max) -r [Receptor.PDBQT] (Optional: -q [[Ion+1](=0.001M),[Ion-1](=0.001M),[Ion+2](=0M),[Ion-2](=0M)] -p GRID_PADDING(=40A) -s GRID_SPACING(=0.375A))\n");
}


int main(int argc, char **argv) {
  string datfn, ligfn, recfn, fldfn, stoken, token;
  double T = 298.;
  double diel_rec = 78.5;
  double ions_mono[2] = { 0.001, 0.001 };//Molar
  double ions_di[2] = { 0., 0. };
  double grid_resolution = 0.375;
  double padding = 40., padding_sqr;

  // receptor pdbqt
  if(!getInputWithFlag(argc, argv, 'r', &recfn)) { usage(); return -1; }
  // number of processors/threads
  if(getInputWithFlag(argc, argv, 'n', &stoken)) {
    __cilkrts_set_param("nworkers", stoken.c_str());
  }
  // ion concentrations
  if(getInputWithFlag(argc, argv, 'q', &stoken)) {
    parseNextValue(&stoken, &token);
    ions_mono[0] = stringToDouble(token);
    parseNextValue(&stoken, &token);
    ions_mono[1] = stringToDouble(token);
    parseNextValue(&stoken, &token);
    ions_di[0] = stringToDouble(token);
    parseNextValue(&stoken, &token);
    ions_di[1] = stringToDouble(token);
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

  // Load receptor file
  cout << "> Loading receptor PDBQT..." << endl;
  ReceptorPDBQT *rec = new ReceptorPDBQT(recfn, NULL);
  cout << "> Receptor center: " << rec->center.x << ", " << rec->center.y << ", " << rec->center.z << endl;
  cout << "> Receptor minimum: " << rec->min.x << ", " << rec->min.y << ", " << rec->min.z << endl;
  cout << "> Receptor maximum coordinates: " << rec->max.x << ", " << rec->max.y << ", " << rec->max.z << endl;

  // Calculate grid geometries
  vertex origin = { rec->min.x - padding, rec->min.y - padding, rec->min.z - padding };
  vertex dimensions = { (rec->max.x + padding) - origin.x, (rec->max.y + padding) - origin.y, (rec->max.z + padding) - origin.z };
  int Npoints[3] = { dimensions.x / grid_resolution, dimensions.y / grid_resolution, dimensions.z / grid_resolution };

  //// Data for all grids
  // Z vector storage for E calcs
  double *data_e = (double*)calloc(Npoints[2], sizeof(double));
  if(!data_e) { cout << "! Error: Could not allocate memory for grid calculation." << endl; return EXIT_FAILURE; }
  ofstream bpm_e("e.bpm", ios::out | ios::binary);
  bpm_e.write((char*)&origin, sizeof(vertex));
  bpm_e.write((char*)&Npoints, sizeof(int)*3);
  bpm_e.write((char*)&grid_resolution, sizeof(double));


  // Electrostatic grid parameters
  double ions = (/*1+*/1. * (ions_mono[0] * 6.022e23 / 1e27))
              + (/*1-*/1. * (ions_mono[1] * 6.022e23 / 1e27))
              + (/*2+*/4. * (ions_di[0] * 6.022e23 / 1e27))
              + (/*2-*/4. * (ions_di[1] * 6.022e23 / 1e27));
  double kappa = sqrt((4. * M_PI * kC * ions) / (diel_rec * kB * T));

  // Start a timer and start our calculations
  Timer *timer = new Timer();
  timer->start();
  for(int nx=0, nt=0; nx < Npoints[0]; nx++) {
    double X = (nx * grid_resolution) + origin.x;
    for(int ny=0; ny < Npoints[1]; ny++) {
      double Y = (ny * grid_resolution) + origin.y;
      cilk_for(int nz=0; nz < Npoints[2]; nz++) {
        double Z = (nz * grid_resolution) + origin.z;

        data_e[nz] = 0;

        for(int i=0; i < rec->coordinates.size(); i++) {
          if(rec->charges[i] == 0.) continue;
          vertex Rrec = rec->coordinates[i];
          vertex dr;
          dr.x = X - Rrec.x; 
          dr.y = Y - Rrec.y; 
          dr.z = Z - Rrec.z; 
          double dist_sqr = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
          if(dist_sqr > padding_sqr) continue;
          double dist = sqrt(dist_sqr);
          // Electrostatic
          double du_e = kC * rec->charges[i] * exp(-dist * kappa) / (diel_rec * dist);
          data_e[nz] += du_e;
        }
      }

      bpm_e.write((char*)data_e, sizeof(double) * Npoints[2]);
    }

    timer->stop();
    double percent_complete = (((double)nx+1) / (double)Npoints[0]);
    cout << "\r> " << (100.*percent_complete) << " percent complete.\t\t\t" << flush;
  }

  bpm_e.close();

  cout << endl << "> Done." << endl;
  return EXIT_SUCCESS;
}
