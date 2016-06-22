#include "Gridder.h"


void usage() {
  printf("Usage: bxm2dx -i [IN.BPM] -o [OUT.DX]\n");
}


int main(int argc, char **argv) {
  string Arg_InFN, Arg_OutFN;

  if(!getInputWithFlag(argc, argv, 'i', &Arg_InFN)) { usage(); return -1; }
  if(!getInputWithFlag(argc, argv, 'o', &Arg_OutFN)) { usage(); return -1; }

  Map_Exclusion *map_ex = new Map_Exclusion(Arg_InFN);
  map_ex->write_dx(Arg_OutFN);

  delete map_ex;
  return 0;
}
