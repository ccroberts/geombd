
#include "Main.h"
#include "Timer.h"
#include "Model.h"


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
  printf("Usage: add -f FOO.maps.fld -l BAR.dlg -o Trajectory.pdb -b B_SPHERE_RADIUS -s BINDINGSITERadius,X,Y,Z [-R REPLICATES(=10000) -n NTHREADS(=max) -T TEMP_KELVIN(=298,room temp) -v VISCOSITY_cP(=1.,Water) -t MAXIMUM_TIME_PS(=100000000ps)]\n");
}


int main(int argc, char **argv) {
  string stoken, fldfn, trjfn;
  SimulationConfig sconfig;

  srand(time(NULL));

  if(!getInputWithFlag(argc, argv, 'f', &fldfn)) { usage(); return -1; }
  if(!getInputWithFlag(argc, argv, 'o', &trjfn)) { usage(); return -1; }
  if(getInputWithFlag(argc, argv, 'n', &stoken)) {
    __cilkrts_set_param("nworkers", stoken.c_str());
    cout << "* Attempting to set number of threads to " << stoken << endl;
  }

  // Create model
  Model *model = new Model(fldfn, trjfn);

  cout << "* Starting simulation with " << model->Nthreads << " threads" << endl;
  Timer *timer = new Timer();
  timer->start();
  model->run();
  timer->stop();
  timer->print(&cout);


  delete timer;
  delete model;

  return EXIT_SUCCESS;
}





















  /*
  if(sconfig == CONFIGURATION_RADIAL) {
    model->r_ligand = stringToDouble(bradius);
    model->r2_escape = pow(model->r_ligand * 5., 2.);
    cout << "> Radial mode set. Bradius: " << model->r_ligand <<"A -- Qradius: " << sqrt(model->r2_escape) << "A" << endl;
  }

  if(sconfig == CONFIGURATION_ABSOLUTE_RADIAL) {
    parseNextValue(&babs, &stoken);
    model->R_ligand.x = stringToDouble(stoken);
    parseNextValue(&babs, &stoken);
    model->R_ligand.y = stringToDouble(stoken);
    parseNextValue(&babs, &stoken);
    model->R_ligand.z = stringToDouble(stoken);
    parseNextValue(&babs, &stoken);
    double qradius = stringToDouble(stoken);
    model->r2_escape = pow(qradius, 2.);

    vertex dr;
    dr.x = model->R_ligand.x - model->center.x;
    dr.y = model->R_ligand.y - model->center.y;
    dr.z = model->R_ligand.z - model->center.z;
    model->r_ligand = sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);

    printf("> Ligand starting position set to (%f, %f, %f) with escape radius %fA\n", model->R_ligand.x, model->R_ligand.y, model->R_ligand.z, qradius); 
  }

  if(sconfig == CONFIGURATION_ABSOLUTE_PERIODIC) {
    parseNextValue(&babs, &stoken);
    model->R_ligand.x = stringToDouble(stoken);
    parseNextValue(&babs, &stoken);
    model->R_ligand.y = stringToDouble(stoken);
    parseNextValue(&babs, &stoken);
    model->R_ligand.z = stringToDouble(stoken);

    parseNextValue(&babs, &stoken);
    model->bounds.x = stringToDouble(stoken);
    parseNextValue(&babs, &stoken);
    model->bounds.y = stringToDouble(stoken);
    parseNextValue(&babs, &stoken);
    model->bounds.z = stringToDouble(stoken);


    printf("> Ligand starting position set to (%f, %f, %f) with periodic boundary (%f, %f, %f)\n", model->R_ligand.x, model->R_ligand.y, model->R_ligand.z, model->bounds.x, model->bounds.y, model->bounds.z); 
  }
  */
