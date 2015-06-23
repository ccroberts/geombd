
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
  printf("Usage: ADD -f FOO.maps.fld -l BAR.dlg -o Trajectory.pdb -b B_SPHERE_RADIUS -s BINDINGSITERadius,X,Y,Z [-R REPLICATES(=10000) -n NTHREADS(=max) -T TEMP_KELVIN(=298,room temp) -v VISCOSITY_cP(=1.,Water) -t MAXIMUM_TIME_PS(=100000000ps)]\n");
}


int main(int argc, char **argv) {
  string fldfn, dlgfn, trjfn, bradius, bsite, stoken, breceptorRoG;

  srand(time(NULL));

  if(!getInputWithFlag(argc, argv, 'f', &fldfn)) { usage(); return -1; }
  if(!getInputWithFlag(argc, argv, 'l', &dlgfn)) { usage(); return -1; }
  if(!getInputWithFlag(argc, argv, 'o', &trjfn)) { usage(); return -1; }
  if(!getInputWithFlag(argc, argv, 'b', &bradius)) { usage(); return -1; }
  if(!getInputWithFlag(argc, argv, 's', &bsite)) { usage(); return -1; }
  if(!getInputWithFlag(argc, argv, 'g', &breceptorRoG)) { usage(); return -1; }
  if(getInputWithFlag(argc, argv, 'n', &stoken)) {
    __cilkrts_set_param("nworkers", stoken.c_str());
    cout << "* Attempting to set number of threads to " << stoken << endl;
  }

  // Create model
  Model *model = new Model(fldfn, dlgfn, trjfn, stringToDouble(breceptorRoG));
  model->r_ligand = stringToDouble(bradius);
  model->ligandPosition = LIGAND_POSITION_RADIAL;
  parseNextValue(&bsite, &stoken);
  model->bindingSite.r2 = pow(stringToDouble(stoken), 2);
  parseNextValue(&bsite, &stoken);
  model->bindingSite.x = stringToDouble(stoken);
  parseNextValue(&bsite, &stoken);
  model->bindingSite.y = stringToDouble(stoken);
  parseNextValue(&bsite, &stoken);
  model->bindingSite.z = stringToDouble(stoken);
  cout << "* Binding site set to: (" << model->bindingSite.x << ", " << model->bindingSite.y << ", " << model->bindingSite.z << ") Radius: " << sqrt(model->bindingSite.r2) << endl;

  // Set optional parameters
  if(getInputWithFlag(argc, argv, 'R', &stoken)) { 
    model->Nreplicates = stringToInt(stoken);
    cout << "* Number of replicate ligand simulations set to " << model->Nreplicates << endl;
  }
  if(getInputWithFlag(argc, argv, 'T', &stoken)) { 
    model->T = stringToDouble(stoken);
    cout << "* Temperature set to " << model->T << "K" << endl;
  }
  if(getInputWithFlag(argc, argv, 'v', &stoken)) { 
    model->viscosity = stringToDouble(stoken)/*cP*/ * 0.001/*cP->Pa.s->J.s/m^3*/ * 2.39e-4 /*J->kcal*/ * 1e12 /*s->ps*/ * 1e-30 /*1/m^3->1/A^3*/ * Na /*kcal->kcal/mol*/;
    cout << "* Viscosity set to " << stoken << "cP" << endl;
  }
  if(getInputWithFlag(argc, argv, 't', &stoken)) { 
    model->t_limit = stringToDouble(stoken);
    cout << "* Maximum simulation time set to " << (model->t_limit/1000.) << " ns." << endl;
  }

  model->Nthreads = __cilkrts_get_nworkers();
  cout << "* Parsing DLG..." << endl;
  model->parseDLG();
  cout << "* Parsing FLD..." << endl;
  model->parseFLD();
  cout << "* Initializing RNG..." << endl;
  model->initializeRNG(time(NULL));

  cout << "* Starting simulation." << endl;
  Timer *timer = new Timer();
  timer->start();
  model->run();
  timer->stop();
  timer->print(&cout);


  delete timer;
  delete model;

  return EXIT_SUCCESS;
}
