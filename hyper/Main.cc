
#include "../Main.h"
#include "Strings.h"
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
  printf("Usage: hyper -d [AD4.1.DAT] -n NTHREADS(=max) -r [Receptor.PDBQT] -l [Ligand.PDBQT]\n");
}


int main(int argc, char **argv) {
  string datfn, ligfn, recfn, stoken;

  if(!getInputWithFlag(argc, argv, 'd', &datfn)) { usage(); return -1; }
  if(!getInputWithFlag(argc, argv, 'r', &recfn)) { usage(); return -1; }
  if(!getInputWithFlag(argc, argv, 'l', &ligfn)) { usage(); return -1; }
  if(getInputWithFlag(argc, argv, 'n', &stoken)) {
    __cilkrts_set_param("nworkers", stoken.c_str());
  }

  // Load AD parameters
  AutoDockParameters *adp = new AutoDockParameters(datfn);
  set<string> atom_types;

  // Load ligand file, check that atom types are in the parameters
  cout << "> Loading ligand PDBQT..." << endl;
  LigandPDBQT *lig = new LigandPDBQT(ligfn);
  cout << "> Types in ligand:";
  set<string>::iterator it;
  for(it=lig->types_set.begin(); it!=lig->types_set.end(); ++it) {
    cout << ' ' << *it;
    atom_types.insert(*it);
  }
  cout << endl;
  for(it=lig->types_set.begin(); it!=lig->types_set.end(); ++it) {
    bool found = false;
    for(int j=0; j < adp->types.size(); j++) {
      if(adp->types[j] == *it) found = true;
    }
    if(!found) {
      cout << "! Error: Ligand atom type " << *it << " not found in the AutoDock parameters." << endl;
    }
  }

  // Load receptor file
  cout << "> Loading receptor PDBQT..." << endl;
  ReceptorPDBQT *rec = new ReceptorPDBQT(recfn, adp);
  cout << "> Receptor center: " << rec->center.x << ", " << rec->center.y << ", " << rec->center.z << endl;
  cout << "> Receptor minimum: " << rec->min.x << ", " << rec->min.y << ", " << rec->min.z << endl;
  cout << "> Receptor maximum coordinates: " << rec->max.x << ", " << rec->max.y << ", " << rec->max.z << endl;
  cout << "> Types in receptor:";
  for(it=rec->types_set.begin(); it!=rec->types_set.end(); ++it) {
    cout << ' ' << *it;
    atom_types.insert(*it);
  }
  cout << endl;
  for(it=rec->types_set.begin(); it!=rec->types_set.end(); ++it) {
    bool found = false;
    for(int j=0; j < adp->types.size(); j++) {
      if(adp->types[j] == *it) found = true;
    }
    if(!found) {
      cout << "! Error: Ligand atom type " << *it << " not found in the AutoDock parameters." << endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
