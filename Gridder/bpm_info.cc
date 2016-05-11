#include <fstream>
#include <cstdio>
#include <cstdlib>
using namespace std;

int main(int argc, char **argv) {
  ifstream fd(argv[1], ios::in | ios::binary);
  double origin[3];
  int Npoints[3];
  int NpointsT = 0;
  double grid_resolution;

  fd.read((char*)&origin, sizeof(double) * 3);
  fd.read((char*)&Npoints, sizeof(int) * 3);
  fd.read((char*)&grid_resolution, sizeof(double));

  NpointsT = Npoints[0] * Npoints[1] * Npoints[2];

  printf("origin %12.6e %12.6e %12.6e\n", origin[0], origin[1], origin[2]);
  printf("N %d %d %d\n", Npoints[0], Npoints[1], Npoints[2]);
  printf("Npoints %d\n", NpointsT);
  printf("delta %12.6e\n", grid_resolution);

  return 0;
}
