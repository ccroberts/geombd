
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
  printf("Usage: GeomBD2 -f INPUTFILE -o TRAJECTORY.pqr -l LOGFILE.log\n");
}


static Model *model = NULL;

void term(int signal) {
  if(model) {
    model->done = true;
  }
}

int main(int argc, char **argv) {
  string stoken, fldfn, trjfn, logfn;

  srand(time(NULL));

  if(!getInputWithFlag(argc, argv, 'f', &fldfn)) { usage(); return -1; }
  if(!getInputWithFlag(argc, argv, 'o', &trjfn)) { usage(); return -1; }
  if(!getInputWithFlag(argc, argv, 'l', &logfn)) { usage(); return -1; }

  // Exit gracefully if possible
  struct sigaction action;
  memset(&action, 0, sizeof(struct sigaction));
  action.sa_handler = term;
  sigaction(SIGINT, &action, NULL);
  sigaction(SIGTERM, &action, NULL);
  sigaction(SIGQUIT, &action, NULL);

  // Create model
  model = new Model(fldfn, trjfn, logfn);

  model->lout << "* Starting simulation with " << model->Nthreads << " threads" << endl;
  Timer *timer = new Timer();
  timer->start();
  model->run();
  timer->stop();
  timer->print(&model->lout);


  delete timer;
  delete model;

  return EXIT_SUCCESS;
}



