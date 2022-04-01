#include "Gridder.h"


void usage() {
  printf("Usage: Gridder-ES -d [FF.gbdp] -r [Receptor.PQR] -l [Ligand.PQR] -o [Output Prefix] -n [N_THREADS] (-p [Padding=40A] -s [Spacing=0.5A] -q [[Ion+1](=0.001M),[Ion-1](=0.001M),[Ion+2](=0M),[Ion-2](=0M)])\n");
}


int main(int argc, char **argv) {
  string Arg, Arg_Q, Arg_ParamFN, Arg_ReceptorFN, Arg_LigandFN, Arg_OutputPrefix;
  int Arg_Padding = 40, Arg_Padding2 = 40*40;
  double Arg_GridSpacing = 0.5;
  double Arg_Ions[4] = { 0.001, 0.001, 0.0, 0.0 };//Molar

  // Obtain command line input
  if(!getInputWithFlag(argc, argv, 'd', &Arg_ParamFN)) { usage(); return -1; }
  if(!getInputWithFlag(argc, argv, 'r', &Arg_ReceptorFN)) { usage(); return -1; }
  if(!getInputWithFlag(argc, argv, 'l', &Arg_LigandFN)) { usage(); return -1; }
  if(getInputWithFlag(argc, argv, 'n', &Arg)) {
    __cilkrts_set_param("nworkers", Arg.c_str());
  }
  if(getInputWithFlag(argc, argv, 'p', &Arg)) {
    Arg_Padding = stringToDouble(Arg);
    Arg_Padding2 = pow(Arg_Padding, 2);
  }
  if(getInputWithFlag(argc, argv, 's', &Arg)) {
    Arg_GridSpacing = stringToDouble(Arg);
  }
  if(getInputWithFlag(argc, argv, 'o', &Arg_OutputPrefix)) {
    cout << "> Output filename prefix: " << Arg_OutputPrefix << endl;
  }
/* POTENTIAL-SPECIFIC */
  if(getInputWithFlag(argc, argv, 'q', &Arg)) {
    parseNextValue(&Arg, &Arg_Q);
    Arg_Ions[0] = stringToDouble(Arg_Q);
    parseNextValue(&Arg, &Arg_Q);
    Arg_Ions[1] = stringToDouble(Arg_Q);
    parseNextValue(&Arg, &Arg_Q);
    Arg_Ions[2] = stringToDouble(Arg_Q);
    parseNextValue(&Arg, &Arg_Q);
    Arg_Ions[3] = stringToDouble(Arg_Q);
  }
/**********************/

  // Load GeomBD3 parameters
  GBD3Parameters *params = new GBD3Parameters(Arg_ParamFN);

  // Load ligand file
  cout << "> Loading ligand PQR..." << endl;
  Molecule_PQRE *lig = new Molecule_PQRE(Arg_LigandFN, params);
  if(lig->types_set.size() == 0) {
    cout << "! Error: No ligand atom types found." << endl;
    exit(-1);
  } else {
    lig->print_types();
    if(! lig->check_types()) exit(EXIT_FAILURE);
  }

  // Load receptor file
  cout << "> Loading receptor PQR..." << endl;
  Molecule_PQRE *rec = new Molecule_PQRE(Arg_ReceptorFN, params);
  if(rec->types_set.size() == 0) {
    cout << "! Error: No receptor atom types found." << endl;
    exit(-1);
  } else {
    rec->print_types();
    if(! rec->check_types()) exit(EXIT_FAILURE);
  }

  // Create exclusion map
  Map_Exclusion *map_ex = new Map_Exclusion(rec, Arg_GridSpacing, Arg_Padding);
  map_ex->calculate();

/* POTENTIAL-SPECIFIC */
  // Create potential map(s)
  string bpmfn = Arg_OutputPrefix;
  bpmfn.append("es.bpm");
  Map_Potential *map_es = new Map_Potential(bpmfn, 0, rec, Arg_GridSpacing, Arg_Padding);

  // Electrostatic grid parameters
  double diel_rec = 78.5; //TODO Make this an input parameter
  double T = 298.;        //TODO Make this an input parameter
  double ions = (/*1+*/1. * (Arg_Ions[0] * 6.022e23 / 1e27))
              + (/*1-*/1. * (Arg_Ions[1] * 6.022e23 / 1e27))
              + (/*2+*/4. * (Arg_Ions[2] * 6.022e23 / 1e27))
              + (/*2-*/4. * (Arg_Ions[3] * 6.022e23 / 1e27));
  double kappa = sqrt((4. * M_PI * kC * ions) / (diel_rec * kB * T));
/**********************/

  // Calculation storage
  vector<int> relevant_atoms; // Dynamic relevant atom list

  // Start our calculations
  cout << "> Starting grid calculation." << endl;
  for(int nx=0, nt=0; nx < map_ex->Npoints[0]; nx++) {
    double X = (nx * Arg_GridSpacing) + map_ex->origin.x;

    relevant_atoms.clear();
    for(int i=0; i < rec->coordinates.size(); i++) {
      vertex Rrec = rec->coordinates[i];
      if(Rrec.x <= X + Arg_Padding and Rrec.x >= X - Arg_Padding)
        relevant_atoms.push_back(i);
    }

    for(int ny=0; ny < map_ex->Npoints[1]; ny++) {
      double Y = (ny * Arg_GridSpacing) + map_ex->origin.y;
      cilk_for(int nz=0; nz < map_ex->Npoints[2]; nz++) {
        double Z = (nz * Arg_GridSpacing) + map_ex->origin.z;

        if(map_ex->data[nx][ny][nz] == false) continue;

/* POTENTIAL-SPECIFIC */
        map_es->data_t[nz] = 0;
/**********************/

        for(int i=0; i < relevant_atoms.size(); i++) {
          vertex Rrec = rec->coordinates[relevant_atoms[i]];
          vertex dr;
          dr.x = X - Rrec.x; 
          dr.y = Y - Rrec.y; 
          if(dr.y > Arg_Padding) continue;
          dr.z = Z - Rrec.z; 
          if(dr.z > Arg_Padding) continue;
          double dist_sqr = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
          if(dist_sqr > Arg_Padding2) continue;

/* POTENTIAL-SPECIFIC CALCULATION */
          double dist = sqrt(dist_sqr);
          double du_t;
          if(ions == 0) du_t = kC * rec->charges[relevant_atoms[i]] / dist;
                   else du_t = kC * rec->charges[relevant_atoms[i]] * exp(-dist * kappa) / (diel_rec * dist);
          map_es->data_t[nz] += du_t;
/**********************************/
        }

        // Clamp potential value at 5000
/* POTENTIAL-SPECIFIC */
        if(map_es->data_t[nz] > 5000) map_es->data_t[nz] = 5000.;
/**********************/
      }

/* POTENTIAL-SPECIFIC */
      map_es->bpm_t->write((char*)map_es->data_t, sizeof(double) * map_es->Npoints[2]);
/**********************/
      double percent_complete = ( (((double)nx) * map_ex->Npoints[1] * map_ex->Npoints[2]) + (((double)ny) * map_ex->Npoints[2]) ) / (double)map_ex->Ntotal;
      cout << "\r> " << (100.*percent_complete) << " percent complete.                                            " << flush;
    }

  }
  cout << endl << "> 100 percent complete." << endl;

/* POTENTIAL-SPECIFIC */
  delete map_es;
/**********************/
  delete map_ex;
  delete rec;
  delete lig;
  return EXIT_SUCCESS;
}
