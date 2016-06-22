#include "Gridder.h"


void usage() {
  printf("Usage: Gridder-EX -r [Receptor.PQR] -o [Output Prefix] -n [N_THREADS] (-p [Padding=5A] -s [Spacing=1.0A])\n");
}


int main(int argc, char **argv) {
  string Arg, Arg_ReceptorFN, Arg_OutputPrefix;
  double Arg_Padding = 5, Arg_Padding2 = 5*5;
  double Arg_GridSpacing = 1.0;
  double Arg_Ions[4] = { 0.001, 0.001, 0.0, 0.0 };//Molar

  // Obtain command line input
  if(!getInputWithFlag(argc, argv, 'r', &Arg_ReceptorFN)) { usage(); return -1; }
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
    Arg_OutputPrefix.append("ex.bxm");
  }
/* POTENTIAL-SPECIFIC */
/**********************/

  // Load receptor file
  cout << "> Loading receptor PQR..." << endl;
  Molecule_PQRE *rec = new Molecule_PQRE(Arg_ReceptorFN, NULL);
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
  map_ex->write(Arg_OutputPrefix);

  delete map_ex;
  delete rec;
  return EXIT_SUCCESS;
}
