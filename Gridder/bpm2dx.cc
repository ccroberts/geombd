#include "Gridder.h"

void usage() {
  printf("Usage: bpm2dx -i [IN.BPM] -o [OUT.DX]\n");
}

int main(int argc, char **argv) {
  string Arg_InFN, Arg_OutFN;

  // Obtain command line input
  if(!getInputWithFlag(argc, argv, 'i', &Arg_InFN)) { usage(); return -1; }
  if(!getInputWithFlag(argc, argv, 'o', &Arg_OutFN)) { usage(); return -1; }

  ifstream fd(Arg_InFN, ios::in | ios::binary);
  double origin[3];
  int Npoints[3];
  int NpointsT = 0;
  double grid_resolution;

  fd.read((char*)&origin, sizeof(double) * 3);
  fd.read((char*)&Npoints, sizeof(int) * 3);
  fd.read((char*)&grid_resolution, sizeof(double));

  NpointsT = Npoints[0] * Npoints[1] * Npoints[2];

  // allocate memory
  double ***data = (double***)calloc(Npoints[0], sizeof(double**));
  if(data == NULL) printf("Error allocating memory.\n");
  for(int nx=0; nx < Npoints[0]; nx++) {
    data[nx] = (double**)calloc(Npoints[1], sizeof(double*));              
    if(data[nx] == NULL) printf("Error allocating memory.\n");
    for(int ny=0; ny < Npoints[1]; ny++) {
      data[nx][ny] = (double*)calloc(Npoints[2], sizeof(double));              
      if(data[nx][ny] == NULL) printf("Error allocating memory.\n");
    }
  }

  for(int nx=0; nx < Npoints[0]; nx++) {
    for(int ny=0; ny < Npoints[1]; ny++) {
      for(int nz=0; nz < Npoints[2]; nz++) {
        fd.read((char*)&data[nx][ny][nz], sizeof(double));
      }
    }
  }

  printf("> Starting to write OpenDX potential map...\n");
  // output opendx data
  FILE *fdo;
  fdo = fopen(Arg_OutFN.c_str(), "w");
  fprintf(fdo, "object 1 class gridpositions counts %d %d %d\n", Npoints[0], Npoints[1], Npoints[2]);
  fprintf(fdo, "origin %12.6e %12.6e %12.6e\n", origin[0], origin[1], origin[2]);
  fprintf(fdo, "delta %12.6e %12.6e %12.6e\n", grid_resolution, 0., 0.);
  fprintf(fdo, "delta %12.6e %12.6e %12.6e\n", 0., grid_resolution, 0.);
  fprintf(fdo, "delta %12.6e %12.6e %12.6e\n", 0., 0., grid_resolution);
  fprintf(fdo, "object 2 class gridconnections counts %d %d %d\n", Npoints[0], Npoints[1], Npoints[2]);
  fprintf(fdo, "object 3 class array type double rank 0 items %d data follows\n", NpointsT);
  int i[3] = { 0, 0, 0 };
  for(int it=0; it < NpointsT; it++) {
    fprintf(fdo, "%12.6e ", data[i[0]][i[1]][i[2]]);
    if((it+1) % 3 == 0) fprintf(fdo, "\n");
    i[2]++;
    if(i[2] >= Npoints[2]) { i[2] = 0; i[1]++; }
    if(i[1] >= Npoints[1]) { i[1] = 0; i[0]++; }
  }
  if(NpointsT % 3 != 0) fprintf(fdo, "\n");
  fprintf(fdo, "attribute \"dep\" string \"positions\"\n");
  fprintf(fdo, "object \"regular positions regular connections\" class field\n");
  fprintf(fdo, "component \"positions\" value 1\n");
  fprintf(fdo, "component \"connections\" value 2\n");
  fprintf(fdo, "component \"data\" value 3\n");
  fclose(fdo);

  // free memory
  for(int nx=0; nx < Npoints[0]; nx++) {
    for(int ny=0; ny < Npoints[1]; ny++) {
      free(data[nx][ny]);
    }
    free(data[nx]);
  }
  free(data);

  printf("> Done.\n");

  return 0;
}
